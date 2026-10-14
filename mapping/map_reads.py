from __future__ import annotations

import time
from pathlib import Path

from mapping.common import (
    MappingCounters,
    create_logger,
    read_fastq_batches,
)
from mapping.mappy_mapper import MappyMapper
from mapping.alt_tie_policy import prefer_primary_alt_ties
from mapping.parquet_writer import (
    MappingBatchPartWriter,
)
from mapping.reference_id import (
    get_group_id,
)


def _count_candidate_groups(
    alignments,
) -> int:
    """Return the number of distinct candidate groups."""

    return len(
        {
            get_group_id(
                alignment.reference_name
            )
            for alignment
            in alignments
        }
    )


def run_transcript_mapping_to_parquet(
    fastq_path: str,
    transcript_fasta: str,
    mapping_parts_directory: str,
    maximum_reads: int | None = None,
    maximum_candidates: int = 16,
    kmer_size: int = 11,
    window_size: int = 5,
    minimum_minimizers: int = 2,
    minimum_chain_score: int = 15,
    minimum_dp_score: int = 15,
    minimum_query_coverage: float = 0.85,
    minimum_identity: float = 0.90,
    mapping_batch_size: int = 100000,
    parquet_row_group_size: int = 100000,
    progress_interval: int = 100000,
    log_file: str | None = None,
    prefer_primary_on_alt_tie: bool = False,
    reference_metadata_tsv: str | None = None,
    reference_gtf: str | None = None,
) -> MappingCounters:
    """
    Map reads and write chromosome/strand Parquet part files.

    Every FASTQ batch creates new part files. Each file is closed before
    the next chromosome/strand file or mapping batch is written.
    """

    if mapping_batch_size < 1:
        raise ValueError(
            "mapping_batch_size must be at least 1."
        )

    output_directory = Path(
        mapping_parts_directory
    )
    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    if log_file is None:
        log_file = str(
            output_directory
            / "mapping.log"
        )

    logger = create_logger(
        "RiboTxAlignSingle",
        log_file,
    )

    mapper = MappyMapper(
        reference_path=transcript_fasta,
        reference_metadata_tsv=reference_metadata_tsv,
        reference_gtf=reference_gtf,
        maximum_candidates=maximum_candidates,
        index_threads=1,
        kmer_size=kmer_size,
        window_size=window_size,
        minimum_minimizers=minimum_minimizers,
        minimum_chain_score=minimum_chain_score,
        minimum_dp_score=minimum_dp_score,
        minimum_query_coverage=minimum_query_coverage,
        minimum_identity=minimum_identity,
    )

    counters = MappingCounters()
    start_time = time.time()
    next_progress = progress_interval

    logger.info(
        "Single-process transcript mapping started"
    )
    logger.info(
        "FASTQ: %s",
        fastq_path,
    )
    logger.info(
        "Transcript FASTA: %s",
        transcript_fasta,
    )
    logger.info(
        "Mapping-parts directory: %s",
        output_directory,
    )

    for batch_index, batch in enumerate(
        read_fastq_batches(
            fastq_path=fastq_path,
            batch_size=mapping_batch_size,
            maximum_reads=maximum_reads,
        )
    ):
        batch_writer = MappingBatchPartWriter(
            output_directory=str(
                output_directory
            ),
            batch_index=batch_index,
            row_group_size=(
                parquet_row_group_size
            ),
        )

        with batch_writer:
            for read in batch:
                counters.total_reads += 1

                alignments, limit_reached = (
                    mapper.map_read(
                        read_id=read.read_id,
                        sequence=read.sequence,
                    )
                )

                if limit_reached:
                    counters.best_hit_overflow_reads += 1

                if not alignments:
                    counters.unmapped_reads += 1
                    continue

                alignments, removed_alt_hits = prefer_primary_alt_ties(
                    alignments, enabled=prefer_primary_on_alt_tie, limit_reached=limit_reached,
                )
                counters.primary_preferred_alt_tie_reads += int(removed_alt_hits > 0)
                counters.removed_alt_tie_hits += removed_alt_hits
                counters.mapped_reads += 1

                group_count = batch_writer.add_read(
                    read_id=read.read_id,
                    read_length=len(
                        read.sequence
                    ),
                    alignments=alignments,
                )

                if group_count == 1:
                    counters.group_unique_reads += 1
                else:
                    counters.group_multi_reads += 1

                if (
                    progress_interval > 0
                    and counters.total_reads
                    >= next_progress
                ):
                    elapsed = (
                        time.time()
                        - start_time
                    )
                    speed = (
                        counters.total_reads
                        / elapsed
                        if elapsed > 0
                        else 0.0
                    )

                    logger.info(
                        "Processed=%d "
                        "Speed=%.1f reads/s "
                        "Mapped=%d "
                        "GroupUnique=%d "
                        "GroupMulti=%d "
                        "Unmapped=%d "
                        "BestHitLimitReached=%d",
                        counters.total_reads,
                        speed,
                        counters.mapped_reads,
                        counters.group_unique_reads,
                        counters.group_multi_reads,
                        counters.unmapped_reads,
                        counters.best_hit_overflow_reads,
                    )

                    while (
                        next_progress
                        <= counters.total_reads
                    ):
                        next_progress += (
                            progress_interval
                        )

            output_files = (
                batch_writer.write_and_close()
            )

        logger.info(
            "Batch=%d InputReads=%d PartFiles=%d",
            batch_index,
            len(batch),
            len(output_files),
        )

    elapsed = (
        time.time()
        - start_time
    )
    speed = (
        counters.total_reads
        / elapsed
        if elapsed > 0
        else 0.0
    )

    logger.info(
        "Single-process transcript mapping finished"
    )
    logger.info(
        "Reads: %d",
        counters.total_reads,
    )
    logger.info(
        "Mapped: %d",
        counters.mapped_reads,
    )
    logger.info(
        "Group unique: %d",
        counters.group_unique_reads,
    )
    logger.info(
        "Group multi: %d",
        counters.group_multi_reads,
    )
    logger.info(
        "Unmapped: %d",
        counters.unmapped_reads,
    )
    logger.info(
        "Best-hit limit reached: %d",
        counters.best_hit_overflow_reads,
    )
    logger.info(
        "Elapsed seconds: %.1f",
        elapsed,
    )
    logger.info(
        "Average speed: %.1f reads/s",
        speed,
    )

    logger.info("Primary-preferred ALT ties: %d reads, %d ALT hits removed",
                counters.primary_preferred_alt_tie_reads, counters.removed_alt_tie_hits)
    return counters
