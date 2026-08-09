from __future__ import annotations

import logging
from pathlib import Path

from bam_output.fastq_order_partition import (
    partition_assignments_by_fastq_order,
)
from bam_output.bam_writer import (
    write_bam_from_fastq_order_partitions,
)


FASTQ_PATH = Path(
    "/mnt/share/project/riboML/riboseq_raw/trimmed_26_34/deduplicated/SRR26942051.26_34.deduplicated.fastq.gz"
)

FINAL_ASSIGNMENT_PARQUET = Path(
    "/mnt/share/project/riboML/sample_ref/bam/"
    "final_genome_assignments.parquet"
)

PARTITION_DIR = Path(
    "/mnt/share/project/riboML/sample_ref/bam/assignment_parts"
)

PARTITION_PREFIX = "final_genome_assignments"
PARTITION_COUNT = 10


GENOME_FAI = Path(
    "/mnt/ssdnas07/pipeline/rna_modscore_v2/source/genome/hg38.fa.fai"
)

OUTPUT_BAM = Path(
    "/mnt/share/project/riboML/sample_ref/bam/"
    "RiboTxAlign.genome.unsorted.bam"
)


def create_logger() -> logging.Logger:
    """Create a console logger for the BAM workflow."""

    logger = logging.getLogger("RiboTxAlignBam")
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
    """Partition final assignments and create an unsorted BAM."""

    logger = create_logger()

    partition_assignments_by_fastq_order(
        assignment_parquet=FINAL_ASSIGNMENT_PARQUET,
        fastq_path=FASTQ_PATH,
        output_dir=PARTITION_DIR,
        output_prefix=PARTITION_PREFIX,
        partition_count=PARTITION_COUNT,
        batch_size=250_000,
        compression="zstd",
        logger=logger,
    )

    write_bam_from_fastq_order_partitions(
        fastq_path=FASTQ_PATH,
        partition_dir=PARTITION_DIR,
        partition_prefix=PARTITION_PREFIX,
        partition_count=PARTITION_COUNT,
        genome_fai=GENOME_FAI,
        output_bam=OUTPUT_BAM,
        sort_output=False,
        index_output=False,
        logger=logger,
    )


if __name__ == "__main__":
    main()