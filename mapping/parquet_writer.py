from __future__ import annotations

from pathlib import Path
from typing import Iterable

from mapping.mapping_result import ProjectedAlignment


class MappingParquetWriter:
    """Write multi-mapping candidate alignments in bounded Parquet batches."""

    def __init__(self, output_parquet: str, batch_size: int = 100000):
        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError as error:
            raise ImportError(
                "Writing mapping candidates to Parquet requires pyarrow. "
                "Install it with: pip install pyarrow"
            ) from error

        if batch_size < 1:
            raise ValueError("batch_size must be at least 1.")

        self.pa = pa
        self.pq = pq
        self.output_parquet = Path(output_parquet)
        self.output_parquet.parent.mkdir(parents=True, exist_ok=True)
        self.batch_size = batch_size
        self.rows: list[dict] = []
        self.writer = None

    def write_read(
        self,
        read_id: str,
        read_length: int,
        genome_status: str,
        genome_hit_count: int,
        alignments: Iterable[ProjectedAlignment],
    ) -> None:
        alignments = list(alignments)
        hit_count = len(alignments)

        for hit_index, alignment in enumerate(alignments, start=1):
            self.rows.append(
                {
                    "read_id": read_id,
                    "read_length": int(read_length),
                    "genome_status": genome_status,
                    "genome_hit_count": int(genome_hit_count),
                    # Counts are based on distinct genomic placements after
                    # transcript hits have been projected and deduplicated.
                    "candidate_count": int(hit_count),
                    "genomic_hit_count": int(hit_count),
                    "candidate_index": int(hit_index),
                    "genomic_hit_index": int(hit_index),
                    "gene_id": alignment.gene_id,
                    "record_type": alignment.record_type,
                    "source_reference": alignment.source_reference,
                    "source_reference_start": int(alignment.source_reference_start),
                    "source_reference_end": int(alignment.source_reference_end),
                    "chrom": alignment.chrom,
                    "reference_start": int(alignment.reference_start),
                    "genomic_blocks": _blocks_to_string(alignment.genomic_blocks),
                    "is_reverse": bool(alignment.is_reverse),
                    "mapq": int(alignment.mapq),
                    "edit_distance": int(alignment.edit_distance),
                    "score": int(alignment.score),
                    "cigar": _cigar_to_string(alignment.cigartuples),
                    "support_sources": ",".join(sorted(alignment.support_sources)),
                }
            )

        if len(self.rows) >= self.batch_size:
            self.flush()

    def flush(self) -> None:
        if not self.rows:
            return

        table = self.pa.Table.from_pylist(self.rows)

        if self.writer is None:
            self.writer = self.pq.ParquetWriter(
                str(self.output_parquet),
                table.schema,
                compression="zstd",
            )

        self.writer.write_table(table)
        self.rows.clear()

    def close(self) -> None:
        self.flush()
        if self.writer is not None:
            self.writer.close()
            self.writer = None

    def __enter__(self) -> "MappingParquetWriter":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


def _cigar_to_string(cigartuples: list[tuple[int, int]]) -> str:
    operation_codes = {
        0: "M",
        1: "I",
        2: "D",
        3: "N",
        4: "S",
        5: "H",
        6: "P",
        7: "=",
        8: "X",
    }
    return "".join(
        f"{length}{operation_codes.get(operation, '?')}"
        for operation, length in cigartuples
    )


def _blocks_to_string(blocks: tuple[tuple[int, int], ...]) -> str:
    return ",".join(f"{start}-{end}" for start, end in blocks)
