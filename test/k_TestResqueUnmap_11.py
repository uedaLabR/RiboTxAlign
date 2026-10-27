from __future__ import annotations

import logging
from pathlib import Path

from mapping.genome_rescue import (
    rescue_unique_genome_alignments,
)
from mapping.mapping_stats import (
    load_mapping_stats,
    percentage,
    write_mapping_stats,
)


UNMAPPED_FASTQ = Path(
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.unmapped.fastq.gz"
)

GENOME_FASTA = Path(
    "/mnt/ssdnas07/pipeline/rna_modscore_v2/"
    "source/genome/hg38.fa"
)

GENOME_FAI = Path(
    "/mnt/ssdnas07/pipeline/rna_modscore_v2/"
    "source/genome/hg38.fa.fai"
)

OUTPUT_BAM = Path(
    "/mnt/share/project/riboML/sample_ref/bam/"
    "RiboTxAlign.genome_rescued_unique.unsorted.bam"
)

STATS_JSON = Path(
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.mapping_stats.json"
)

LOG_FILE = Path(
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.mapping_mt.log"
)


def create_logger() -> logging.Logger:
    """Create console and append-mode file logging for genome rescue."""

    logger = logging.getLogger("RiboTxAlignGenomeRescue")
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


def update_mapping_statistics(
    processed_reads: int,
    contiguous_rescued_reads: int,
    splice_rescued_reads: int,
    logger: logging.Logger,
) -> None:
    """Add two-stage genome rescue values to the statistics JSON."""

    stats = load_mapping_stats(STATS_JSON)
    transcript_unmapped = int(stats["transcript_unmapped"])
    if processed_reads != transcript_unmapped:
        raise ValueError(
            "Genome rescue FASTQ count does not match transcript-unmapped "
            f"count: FASTQ={processed_reads}, "
            f"statistics={transcript_unmapped}"
        )

    input_reads = int(stats["input_reads"])
    transcript_mapped = int(stats["transcript_mapped"])
    group_unique = int(stats["group_unique"])
    group_multi = int(stats["group_multi"])
    decoy_filtered = int(stats["decoy_filtered"])
    analyzable_reads = input_reads - decoy_filtered
    rescued_unique_reads = contiguous_rescued_reads + splice_rescued_reads
    final_unmapped = processed_reads - rescued_unique_reads
    final_mapped = transcript_mapped + rescued_unique_reads
    final_unique = group_unique + rescued_unique_reads
    final_multi = group_multi

    stats.update({
        "genome_rescue_input": processed_reads,
        "genome_contiguous_rescued_unique": contiguous_rescued_reads,
        "genome_splice_rescued_unique": splice_rescued_reads,
        "genome_rescued_unique": rescued_unique_reads,
        "genome_rescue_rate_percent": percentage(
            rescued_unique_reads, processed_reads
        ),
        "final_mapped": final_mapped,
        "final_unique": final_unique,
        "final_multi": final_multi,
        "final_unmapped": final_unmapped,
        "final_mapped_rate_percent": percentage(final_mapped, input_reads),
        "final_unique_rate_percent": percentage(final_unique, input_reads),
        "final_multi_rate_percent": percentage(final_multi, input_reads),
        "final_unmapped_rate_percent": percentage(final_unmapped, input_reads),
        "final_mapped_non_decoy_rate_percent": percentage(
            final_mapped, analyzable_reads
        ),
    })
    write_mapping_stats(STATS_JSON, stats)

    logger.info("Genome rescue summary")
    logger.info("Genome rescue input reads: %d", processed_reads)
    logger.info("Contiguous rescued unique reads: %d", contiguous_rescued_reads)
    logger.info("Splice rescued unique reads: %d", splice_rescued_reads)
    logger.info("Total rescued unique reads: %d", rescued_unique_reads)
    logger.info("Remaining unmapped reads: %d", final_unmapped)
    logger.info("Final mapped reads: %d", final_mapped)
    logger.info("Final unique reads: %d", final_unique)
    logger.info("Final multi reads: %d", final_multi)
    logger.info("Mapping statistics JSON updated: %s", STATS_JSON)


def main() -> None:
    """Rescue unique genome mappings and update mapping statistics."""

    logger = create_logger()
    counters = rescue_unique_genome_alignments(
        unmapped_fastq=UNMAPPED_FASTQ,
        genome_fasta=GENOME_FASTA,
        genome_fai=GENOME_FAI,
        output_bam=OUTPUT_BAM,
        remaining_unmapped_fastq=None,
        contiguous_preset="sr",
        splice_preset="splice",
        kmer_size=11,
        window_size=5,
        minimum_minimizers=2,
        minimum_chain_score=15,
        minimum_dp_score=15,
        minimum_query_coverage=0.85,
        minimum_identity=0.90,
        enable_splice_rescue=True,
        minimum_anchor_length=11,
        minimum_intron_length=20,
        maximum_intron_length=2000,
        canonical_motifs=("GT-AG", "GC-AG", "AT-AC"),
        progress_interval=100_000,
        logger=logger,
    )

    update_mapping_statistics(
        processed_reads=counters.processed_reads,
        contiguous_rescued_reads=(
            counters.contiguous_rescued_unique_reads
        ),
        splice_rescued_reads=counters.splice_rescued_unique_reads,
        logger=logger,
    )


if __name__ == "__main__":
    main()
