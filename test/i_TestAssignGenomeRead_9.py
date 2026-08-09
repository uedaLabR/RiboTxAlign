from __future__ import annotations

import logging
from pathlib import Path

from bam_output.final_alignment_assignment import (
    create_final_genome_assignments,
)


READS_INFO_PARQUET = Path(
    "/mnt/share/project/riboML/sample_ref/"
    "transcript_reads_info.parquet"
)

COORDINATE_TABLE = Path(
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.coordinates.tsv"
)

OUTPUT_DIR = Path(
    "/mnt/share/project/riboML/sample_ref/bam"
)

FINAL_ASSIGNMENT_PARQUET = (
    OUTPUT_DIR / "final_genome_assignments.parquet"
)

ASSIGNMENT_STATS_TSV = (
    OUTPUT_DIR / "final_genome_assignments.stats.tsv"
)


def create_logger() -> logging.Logger:
    """Create the console logger."""

    logger = logging.getLogger("RiboTxAlignAssignment")
    logger.setLevel(logging.INFO)

    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter(
                "%(asctime)s\t%(levelname)s\t%(message)s"
            )
        )
        logger.addHandler(handler)

    return logger


def main() -> None:
    """Resolve one final genomic alignment for every mapped read."""

    create_final_genome_assignments(
        reads_info_parquet=READS_INFO_PARQUET,
        coordinate_table=COORDINATE_TABLE,
        output_parquet=FINAL_ASSIGNMENT_PARQUET,
        output_stats_tsv=ASSIGNMENT_STATS_TSV,
        batch_size=250_000,

        # Weak local balancing window for equivalent multi-map candidates.
        local_bin_size=30,

        logger=create_logger(),
    )


if __name__ == "__main__":
    main()