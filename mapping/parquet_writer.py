from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import re
import tempfile

from mapping.mapping_result import RawAlignment
from mapping.reference_id import get_group_id


_CIGAR_CODES = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X"}
_GROUP_ID_PATTERN = re.compile(r"^(?P<chrom>.+)_(?:plus|minus)_G\d+$", re.IGNORECASE)
_COLUMNS = (
    "read_id", "read_length", "hit_index", "group_index", "transcript_id",
    "group_id", "other_group_ids", "reference_start", "reference_end",
    "is_cds", "cds_source", "query_start", "query_end", "strand", "cigar",
    "hit_count", "group_count",
)


def cigar_to_string(cigartuples: list[tuple[int, int]]) -> str:
    return "".join(f"{length}{_CIGAR_CODES.get(operation, '?')}" for operation, length in cigartuples)


@lru_cache(maxsize=65_536)
def get_chromosome_from_group_id(group_id: str) -> str:
    """Extract a chromosome, preserving underscores in ALT/scaffold names."""
    match = _GROUP_ID_PATTERN.fullmatch(group_id)
    if match is None:
        raise ValueError(f"Cannot extract chromosome from group ID: {group_id}")
    return match.group("chrom")


def get_mapping_schema():
    """Return the unchanged mapping-part schema."""
    import pyarrow as pa

    return pa.schema([
        ("read_id", pa.string()), ("read_length", pa.int32()),
        ("hit_index", pa.int16()), ("group_index", pa.int16()),
        ("transcript_id", pa.string()), ("group_id", pa.string()),
        ("other_group_ids", pa.list_(pa.string())),
        ("reference_start", pa.int32()), ("reference_end", pa.int32()),
        ("is_cds", pa.bool_()), ("cds_source", pa.string()),
        ("query_start", pa.int32()), ("query_end", pa.int32()),
        ("strand", pa.int8()), ("cigar", pa.string()),
        ("hit_count", pa.int16()), ("group_count", pa.int16()),
    ])


def _alignment_values(read_id, read_length, alignments):
    """Yield column values in the original stable hit order without row dicts."""
    if not alignments:
        raise ValueError("At least one alignment is required.")
    ordered = sorted(alignments, key=lambda a: (
        get_group_id(a.reference_name), a.reference_name,
        a.reference_start, a.reference_end, a.strand,
    ))
    group_ids = [get_group_id(a.reference_name) for a in ordered]
    candidates = list(dict.fromkeys(group_ids))
    group_index = {gid: index for index, gid in enumerate(candidates)}
    others = {gid: [other for other in candidates if other != gid] for gid in candidates}
    for hit_index, (alignment, gid) in enumerate(zip(ordered, group_ids)):
        yield (
            read_id, int(read_length), hit_index, group_index[gid],
            alignment.reference_name, gid, others[gid],
            int(alignment.reference_start), int(alignment.reference_end),
            alignment.is_cds, alignment.cds_source,
            int(alignment.query_start), int(alignment.query_end), int(alignment.strand),
            cigar_to_string(alignment.cigartuples), len(ordered), len(candidates),
        )


def alignments_to_rows(read_id: str, read_length: int, alignments: list[RawAlignment]) -> list[dict]:
    """Compatibility helper; the streaming writer does not create row dicts."""
    rows = []
    for values in _alignment_values(read_id, read_length, alignments):
        row = dict(zip(_COLUMNS, values))
        row["other_group_ids"] = list(row["other_group_ids"])
        rows.append(row)
    return rows


class MappingBatchPartWriter:
    """Write one part per chromosome with a shared, bounded column buffer.

    row_group_size limits total buffered rows across all chromosomes.
    Flush the largest buffer when the limit is reached. Parquet writers remain
    open until the chunk ends; their overhead scales with chromosome count.
    Temporary files are published only after every writer closes successfully.
    """

    def __init__(self, output_directory: str, batch_index: int,
                 compression: str = "zstd", row_group_size: int = 100000) -> None:
        if batch_index < 0:
            raise ValueError("batch_index must not be negative.")
        if row_group_size < 1:
            raise ValueError("row_group_size must be at least 1.")
        self.output_directory = Path(output_directory)
        self.output_directory.mkdir(parents=True, exist_ok=True)
        self.batch_index = batch_index
        self.compression = compression
        self.row_group_size = row_group_size
        self.schema = get_mapping_schema()
        self._buffers = {}
        self._writers = {}
        self._temporary_paths = {}
        self._paths = {}
        self.buffered_row_count = 0
        self.peak_buffered_rows = 0
        self._closed = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None:
            self.write_and_close()
        else:
            self.abort()
        return False

    def add_read(self, read_id: str, read_length: int, alignments: list[RawAlignment]) -> int:
        """Append a read and return its number of candidate groups."""
        if self._closed:
            raise RuntimeError("Mapping writer is closed.")
        group_count = 0
        for values in _alignment_values(read_id, read_length, alignments):
            chromosome = get_chromosome_from_group_id(values[5])
            columns = self._buffers.get(chromosome)
            if columns is None:
                columns = [[] for _ in _COLUMNS]
                self._buffers[chromosome] = columns
            for column, value in zip(columns, values):
                column.append(value)
            group_count = values[-1]
            self.buffered_row_count += 1
            self.peak_buffered_rows = max(self.peak_buffered_rows, self.buffered_row_count)
            if self.buffered_row_count >= self.row_group_size:
                largest = max(self._buffers, key=lambda chrom: len(self._buffers[chrom][0]))
                self._flush(largest)
        return group_count

    def _flush(self, chromosome):
        import pyarrow as pa
        import pyarrow.parquet as pq

        columns = self._buffers[chromosome]
        if chromosome not in self._writers:
            output_path = self.output_directory / f"{chromosome}.part{self.batch_index:06d}.parquet"
            if output_path.exists():
                raise FileExistsError(f"Mapping part already exists: {output_path}")
            with tempfile.NamedTemporaryFile(dir=self.output_directory,
                    prefix=f".{output_path.name}.", suffix=".tmp", delete=False) as handle:
                temporary = Path(handle.name)
            self._temporary_paths[chromosome] = temporary
            self._paths[chromosome] = output_path
            self._writers[chromosome] = pq.ParquetWriter(
                temporary, self.schema, compression=self.compression,
            )
        table = pa.Table.from_pydict(dict(zip(_COLUMNS, columns)), schema=self.schema)
        self._writers[chromosome].write_table(table, row_group_size=self.row_group_size)
        self.buffered_row_count -= len(columns[0])
        del self._buffers[chromosome]

    def write_and_close(self) -> list[Path]:
        """Flush, close, and publish parts; repeated calls return the same paths."""
        if self._closed:
            return [self._paths[chrom] for chrom in sorted(self._paths)]
        try:
            for chromosome in sorted(self._buffers):
                self._flush(chromosome)
            for writer in self._writers.values():
                writer.close()
            self._writers.clear()
            for path in self._paths.values():
                if path.exists():
                    raise FileExistsError(f"Mapping part already exists: {path}")
            for chromosome, temporary in self._temporary_paths.items():
                temporary.rename(self._paths[chromosome])
            self._temporary_paths.clear()
            self._closed = True
        except BaseException:
            self.abort()
            raise
        return [self._paths[chrom] for chrom in sorted(self._paths)]

    def abort(self):
        """Close handles and discard unpublished data after a mapping failure."""
        try:
            for writer in self._writers.values():
                try:
                    writer.close()
                except Exception:
                    pass
        finally:
            self._writers.clear()
            self._buffers.clear()
            self.buffered_row_count = 0
            self._closed = True
            for temporary in self._temporary_paths.values():
                temporary.unlink(missing_ok=True)
            self._temporary_paths.clear()
            self._paths.clear()
