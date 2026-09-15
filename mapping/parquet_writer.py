from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re

from mapping.mapping_result import RawAlignment
from mapping.reference_id import get_group_id


_CIGAR_CODES = {
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

_GROUP_ID_PATTERN = re.compile(
    r"^(?P<chrom>.+)_(?:plus|minus)_G\d+$",
    re.IGNORECASE,
)


def cigar_to_string(
    cigartuples: list[tuple[int, int]],
) -> str:
    """Convert CIGAR tuples to a compact CIGAR string."""

    return "".join(
        f"{length}{_CIGAR_CODES.get(operation, '?')}"
        for operation, length in cigartuples
    )


def get_chromosome_from_group_id(
    group_id: str,
) -> str:
    """Extract the chromosome name from a group ID.

    Examples
    --------
    ``chr1_plus_G00001`` becomes ``chr1``.
    ``chrX_KI270913v1_alt_minus_G00001`` becomes
    ``chrX_KI270913v1_alt``.
    """

    match = _GROUP_ID_PATTERN.fullmatch(
        group_id
    )

    if match is None:
        raise ValueError(
            "Cannot extract chromosome from group ID: "
            f"{group_id}"
        )

    return match.group(
        "chrom"
    )


def get_mapping_schema():
    """Return the schema used by all mapping part files."""

    try:
        import pyarrow as pa
    except ImportError as error:
        raise ImportError(
            "Parquet output requires pyarrow. "
            "Install it with: pip install pyarrow"
        ) from error

    return pa.schema(
        [
            ("read_id", pa.string()),
            ("read_length", pa.int32()),
            ("hit_index", pa.int16()),
            ("group_index", pa.int16()),
            ("transcript_id", pa.string()),
            ("group_id", pa.string()),
            (
                "other_group_ids",
                pa.list_(
                    pa.string()
                ),
            ),
            ("reference_start", pa.int32()),
            ("reference_end", pa.int32()),
            ("is_cds", pa.bool_()),
            ("cds_source", pa.string()),
            ("query_start", pa.int32()),
            ("query_end", pa.int32()),
            ("strand", pa.int8()),
            ("cigar", pa.string()),
            ("hit_count", pa.int16()),
            ("group_count", pa.int16()),
        ]
    )


def alignments_to_rows(
    read_id: str,
    read_length: int,
    alignments: list[RawAlignment],
) -> list[dict]:
    """Convert retained best-tie alignments to one row per hit."""

    if not alignments:
        raise ValueError(
            "At least one alignment is required."
        )

    ordered_alignments = sorted(
        alignments,
        key=lambda alignment: (
            get_group_id(
                alignment.reference_name
            ),
            alignment.reference_name,
            alignment.reference_start,
            alignment.reference_end,
            alignment.strand,
        ),
    )

    group_ids = [
        get_group_id(
            alignment.reference_name
        )
        for alignment in ordered_alignments
    ]

    candidate_group_ids = list(
        dict.fromkeys(
            group_ids
        )
    )

    group_index_by_id = {
        group_id: group_index
        for group_index, group_id in enumerate(
            candidate_group_ids
        )
    }

    hit_count = len(
        ordered_alignments
    )
    group_count = len(
        candidate_group_ids
    )

    rows: list[dict] = []

    for hit_index, (
        alignment,
        group_id,
    ) in enumerate(
        zip(
            ordered_alignments,
            group_ids,
        )
    ):
        rows.append(
            {
                "read_id": read_id,
                "read_length": int(
                    read_length
                ),
                "hit_index": int(
                    hit_index
                ),
                "group_index": int(
                    group_index_by_id[
                        group_id
                    ]
                ),
                "transcript_id": (
                    alignment.reference_name
                ),
                "group_id": group_id,
                "is_cds": alignment.is_cds,
                "cds_source": alignment.cds_source,
                "other_group_ids": [
                    candidate_group_id
                    for candidate_group_id
                    in candidate_group_ids
                    if candidate_group_id
                    != group_id
                ],
                "reference_start": int(
                    alignment.reference_start
                ),
                "reference_end": int(
                    alignment.reference_end
                ),
                "query_start": int(
                    alignment.query_start
                ),
                "query_end": int(
                    alignment.query_end
                ),
                "strand": int(
                    alignment.strand
                ),
                "cigar": cigar_to_string(
                    alignment.cigartuples
                ),
                "hit_count": int(
                    hit_count
                ),
                "group_count": int(
                    group_count
                ),
            }
        )

    return rows


class MappingBatchPartWriter:
    """Aggregate one FASTQ batch and write one part per chromosome.

    Rows are not sorted here. Each chromosome file is opened, written,
    and closed once after all worker results for the FASTQ batch have
    been aggregated.
    """

    def __init__(
        self,
        output_directory: str,
        batch_index: int,
        compression: str = "zstd",
        row_group_size: int = 100000,
    ) -> None:
        if batch_index < 0:
            raise ValueError(
                "batch_index must not be negative."
            )

        if row_group_size < 1:
            raise ValueError(
                "row_group_size must be at least 1."
            )

        self.output_directory = Path(
            output_directory
        )
        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.batch_index = batch_index
        self.compression = compression
        self.row_group_size = row_group_size
        self.schema = get_mapping_schema()

        self.rows_by_chromosome: dict[
            str,
            list[dict],
        ] = defaultdict(
            list
        )

    def add_read(
        self,
        read_id: str,
        read_length: int,
        alignments: list[RawAlignment],
    ) -> None:
        """Add all retained hits of one read to this FASTQ batch."""

        for row in alignments_to_rows(
            read_id=read_id,
            read_length=read_length,
            alignments=alignments,
        ):
            chromosome = (
                get_chromosome_from_group_id(
                    row[
                        "group_id"
                    ]
                )
            )

            self.rows_by_chromosome[
                chromosome
            ].append(
                row
            )

    def write_and_close(
        self,
    ) -> list[Path]:
        """Write one unsorted part file per represented chromosome."""

        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError as error:
            raise ImportError(
                "Parquet output requires pyarrow. "
                "Install it with: pip install pyarrow"
            ) from error

        output_files: list[Path] = []

        for chromosome in sorted(
            self.rows_by_chromosome
        ):
            rows = self.rows_by_chromosome[
                chromosome
            ]

            if not rows:
                continue

            output_path = (
                self.output_directory
                / (
                    f"{chromosome}."
                    f"part{self.batch_index:06d}."
                    "parquet"
                )
            )

            if output_path.exists():
                raise FileExistsError(
                    "Mapping part already exists: "
                    f"{output_path}"
                )

            table = pa.Table.from_pylist(
                rows,
                schema=self.schema,
            )

            pq.write_table(
                table,
                output_path,
                compression=self.compression,
                row_group_size=(
                    self.row_group_size
                ),
            )

            output_files.append(
                output_path
            )

        self.rows_by_chromosome.clear()

        return output_files
