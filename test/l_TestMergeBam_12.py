from __future__ import annotations

import logging
from pathlib import Path

from mapping.mapping_stats import (
    load_mapping_stats,
    write_mapping_stats,
)


MAIN_BAM = Path(
    "/mnt/share/project/riboML/sample_ref/bam/"
    "RiboTxAlign.genome.unsorted.bam"
)

RESCUE_BAM = Path(
    "/mnt/share/project/riboML/sample_ref/bam/"
    "RiboTxAlign.genome_rescued_unique.unsorted.bam"
)

MERGED_UNSORTED_BAM = Path(
    "/mnt/share/project/riboML/sample_ref/bam/"
    "RiboTxAlign.genome.merged.unsorted.bam"
)

OUTPUT_BAM = Path(
    "/mnt/share/project/riboML/sample_ref/bam/"
    "RiboTxAlign.genome.sorted.bam"
)

STATS_JSON = Path(
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.mapping_stats.json"
)

LOG_FILE = Path(
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.mapping_mt.log"
)

SORT_THREADS = 8
SORT_MEMORY = "2G"


def create_logger() -> logging.Logger:
    """Create console and append-mode file logging for final BAM output."""

    logger = logging.getLogger("RiboTxAlignMergeBam")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()

    formatter = logging.Formatter(
        "%(asctime)s\t%(levelname)s\t%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(
        LOG_FILE,
        mode="a",
        encoding="utf-8",
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def count_bam_records(path: Path, pysam_module) -> int:
    """Count primary BAM records without requiring an index."""

    count = 0
    with pysam_module.AlignmentFile(path, "rb") as bam_handle:
        for record in bam_handle.fetch(until_eof=True):
            if record.is_secondary or record.is_supplementary:
                continue
            count += 1

    return count


def main() -> None:
    """Merge the main and rescue BAMs, then sort, index, and validate."""

    try:
        import pysam
    except ImportError as error:
        raise RuntimeError(
            "BAM merge requires pysam."
        ) from error

    logger = create_logger()
    OUTPUT_BAM.parent.mkdir(parents=True, exist_ok=True)

    main_bam_reads = count_bam_records(MAIN_BAM, pysam)
    rescue_bam_reads = count_bam_records(RESCUE_BAM, pysam)

    logger.info("Main BAM primary reads: %d", main_bam_reads)
    logger.info("Rescue BAM primary reads: %d", rescue_bam_reads)

    logger.info("Merging BAM files")
    pysam.merge(
        "-f",
        str(MERGED_UNSORTED_BAM),
        str(MAIN_BAM),
        str(RESCUE_BAM),
    )

    merged_bam_reads = count_bam_records(
        MERGED_UNSORTED_BAM,
        pysam,
    )
    expected_merged_reads = main_bam_reads + rescue_bam_reads

    if merged_bam_reads != expected_merged_reads:
        raise ValueError(
            "Merged BAM read count is inconsistent: "
            f"main={main_bam_reads}, rescue={rescue_bam_reads}, "
            f"merged={merged_bam_reads}"
        )

    logger.info("Sorting merged BAM")
    pysam.sort(
        "-@",
        str(max(SORT_THREADS, 1)),
        "-m",
        SORT_MEMORY,
        "-o",
        str(OUTPUT_BAM),
        str(MERGED_UNSORTED_BAM),
    )

    logger.info("Indexing sorted BAM")
    pysam.index(
        "-@",
        str(max(SORT_THREADS, 1)),
        str(OUTPUT_BAM),
    )

    final_bam_reads = count_bam_records(OUTPUT_BAM, pysam)
    if final_bam_reads != expected_merged_reads:
        raise ValueError(
            "Final BAM read count is inconsistent: "
            f"expected={expected_merged_reads}, "
            f"final={final_bam_reads}"
        )

    stats = load_mapping_stats(STATS_JSON)
    expected_main_reads = int(stats["transcript_mapped"])
    expected_rescue_reads = int(stats["genome_rescued_unique"])
    expected_final_reads = int(stats["final_mapped"])

    if main_bam_reads != expected_main_reads:
        raise ValueError(
            "Main BAM count does not match mapping statistics: "
            f"BAM={main_bam_reads}, statistics={expected_main_reads}"
        )

    if rescue_bam_reads != expected_rescue_reads:
        raise ValueError(
            "Rescue BAM count does not match mapping statistics: "
            f"BAM={rescue_bam_reads}, statistics={expected_rescue_reads}"
        )

    if final_bam_reads != expected_final_reads:
        raise ValueError(
            "Final BAM count does not match mapping statistics: "
            f"BAM={final_bam_reads}, statistics={expected_final_reads}"
        )

    stats.update(
        {
            "main_bam_primary_reads": main_bam_reads,
            "rescue_bam_primary_reads": rescue_bam_reads,
            "final_bam_primary_reads": final_bam_reads,
            "bam_merge_validated": True,
            "final_bam": str(OUTPUT_BAM),
            "final_bam_index": str(OUTPUT_BAM) + ".bai",
        }
    )
    write_mapping_stats(STATS_JSON, stats)

    MERGED_UNSORTED_BAM.unlink(missing_ok=True)

    logger.info("Merged BAM primary reads: %d", merged_bam_reads)
    logger.info("Final BAM primary reads: %d", final_bam_reads)
    logger.info("BAM merge validation: passed")
    logger.info("Final BAM written: %s", OUTPUT_BAM)
    logger.info("Final BAM index written: %s.bai", OUTPUT_BAM)
    logger.info("Mapping statistics JSON updated: %s", STATS_JSON)


if __name__ == "__main__":
    main()
