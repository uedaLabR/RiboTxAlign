from __future__ import annotations

import time
from fnmatch import fnmatch
from concurrent.futures import (
    ProcessPoolExecutor,
)
from dataclasses import dataclass
from multiprocessing import (
    get_context,
)
from pathlib import Path
import gzip
import shutil

from mapping.common import (
    MappingCounters,
    ReadRecord,
    create_logger,
    read_fastq_batches,
)
from mapping.mappy_mapper import MappyMapper
from mapping.mapping_stats import (
    percentage,
    write_mapping_stats,
)
from mapping.parquet_writer import (
    MappingBatchPartWriter,
)
from mapping.reference_id import (
    get_group_id,
)


_WORKER_MAPPER: MappyMapper | None = None


@dataclass(slots=True)
class WorkerTask:
    """One mapping-and-writing task sent to a worker."""

    reads: list[ReadRecord]
    output_directory: str
    part_index: int
    parquet_row_group_size: int
    exclude_decoy_reads: bool
    decoy_reference_patterns: tuple[str, ...]
    unmapped_fastq_part: str | None


@dataclass(slots=True)
class WorkerResult:
    """Small result returned from a worker."""

    counters: MappingCounters
    mapped_read_count: int
    decoy_filtered_read_count: int
    output_files: list[str]
    unmapped_fastq_part: str | None


def _get_reference_chromosome(
    reference_name: str,
) -> str:
    """Extract the chromosome component from a reference ID."""

    for separator in (
        "_plus_",
        "_minus_",
    ):
        if separator in reference_name:
            return reference_name.split(
                separator,
                maxsplit=1,
            )[0]

    return reference_name


def _parse_decoy_reference_patterns(
    pattern_text: str,
) -> tuple[str, ...]:
    """Parse comma-separated decoy reference patterns."""

    return tuple(
        pattern.strip()
        for pattern in pattern_text.split(
            ","
        )
        if pattern.strip()
    )


def _is_decoy_reference(
    reference_name: str,
    decoy_reference_patterns: tuple[str, ...],
) -> bool:
    """Return whether a reference chromosome matches a decoy pattern."""

    chromosome = _get_reference_chromosome(
        reference_name
    )

    return any(
        fnmatch(
            chromosome,
            pattern,
        )
        for pattern in decoy_reference_patterns
    )


def _contains_decoy_hit(
    alignments,
    decoy_reference_patterns: tuple[str, ...],
) -> bool:
    """Return whether any retained best hit maps to a decoy reference."""

    return any(
        _is_decoy_reference(
            alignment.reference_name,
            decoy_reference_patterns,
        )
        for alignment in alignments
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
            for alignment in alignments
        }
    )


def initialize_worker(
    transcript_fasta: str,
    maximum_candidates: int,
    kmer_size: int,
    window_size: int,
    minimum_minimizers: int,
    minimum_chain_score: int,
    minimum_dp_score: int,
    minimum_query_coverage: float,
    minimum_identity: float,
) -> None:
    """Create one independent mapper in each worker process."""

    global _WORKER_MAPPER

    _WORKER_MAPPER = MappyMapper(
        reference_path=transcript_fasta,
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


def map_and_write_chunk(
    task: WorkerTask,
) -> WorkerResult:
    """
    Map one read chunk and write chromosome-level Parquet parts.

    Mapping records remain inside the worker process. Only counters and
    output filenames are returned to the parent process.
    """

    if _WORKER_MAPPER is None:
        raise RuntimeError(
            "Worker mapper is not initialized."
        )

    counters = MappingCounters()

    writer = MappingBatchPartWriter(
        output_directory=(
            task.output_directory
        ),
        batch_index=(
            task.part_index
        ),
        row_group_size=(
            task.parquet_row_group_size
        ),
    )

    mapped_read_count = 0
    decoy_filtered_read_count = 0
    unmapped_handle = None

    if task.unmapped_fastq_part is not None:
        unmapped_part_path = Path(task.unmapped_fastq_part)
        unmapped_part_path.parent.mkdir(parents=True, exist_ok=True)
        unmapped_handle = gzip.open(
            unmapped_part_path,
            "wt",
            encoding="utf-8",
            newline="",
        )

    try:
        for read in task.reads:
            counters.total_reads += 1

            alignments, limit_reached = (
                _WORKER_MAPPER.map_read(
                    read_id=read.read_id,
                    sequence=read.sequence,
                )
            )

            if limit_reached:
                counters.best_hit_overflow_reads += 1

            if not alignments:
                counters.unmapped_reads += 1

                if unmapped_handle is not None:
                    quality = read.quality
                    if quality is None:
                        quality = "I" * len(read.sequence)
                    unmapped_handle.write(
                        f"@{read.read_id}\n"
                        f"{read.sequence}\n"
                        "+\n"
                        f"{quality}\n"
                    )
                continue

            if (
                task.exclude_decoy_reads
                and _contains_decoy_hit(
                    alignments,
                    task.decoy_reference_patterns,
                )
            ):
                decoy_filtered_read_count += 1
                continue

            counters.mapped_reads += 1
            mapped_read_count += 1

            writer.add_read(
                read_id=read.read_id,
                read_length=len(read.sequence),
                alignments=alignments,
            )

            if _count_candidate_groups(alignments) == 1:
                counters.group_unique_reads += 1
            else:
                counters.group_multi_reads += 1
    finally:
        if unmapped_handle is not None:
            unmapped_handle.close()

    output_paths = writer.write_and_close()

    return WorkerResult(
        counters=counters,
        mapped_read_count=(
            mapped_read_count
        ),
        decoy_filtered_read_count=(
            decoy_filtered_read_count
        ),
        output_files=[
            str(output_path)
            for output_path in output_paths
        ],
        unmapped_fastq_part=task.unmapped_fastq_part,
    )


def split_batch(
    reads: list[ReadRecord],
    worker_count: int,
) -> list[list[ReadRecord]]:
    """Split one FASTQ batch into approximately equal worker chunks."""

    if not reads:
        return []

    actual_workers = min(
        worker_count,
        len(reads),
    )

    chunk_size = (
        len(reads)
        + actual_workers
        - 1
    ) // actual_workers

    return [
        reads[
            start:
            start + chunk_size
        ]
        for start in range(
            0,
            len(reads),
            chunk_size,
        )
    ]


def run_transcript_mapping_to_parquet_mt(
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
    worker_count: int = 6,
    reads_per_worker: int = 50000,
    progress_interval: int = 100000,
    parquet_row_group_size: int = 100000,
    exclude_decoy_reads: bool = True,
    decoy_reference_patterns: str = (
        "chrUn_*,*_random,DECOY_*"
    ),
    log_file: str | None = None,
    unmapped_fastq_path: str | None = None,
    stats_json_path: str | None = None,
    multiprocessing_start_method: str = "spawn",
) -> MappingCounters:
    """
    Map reads and write worker-level chromosome Parquet parts.

    Each worker maps one bounded chunk and writes its own part files.
    Mapping records are therefore not transferred back to the parent.

    Output partitioning must be chromosome-only in
    ``MappingBatchPartWriter``:

        chr1.part000000.parquet
        chr2.part000000.parquet

    Strand remains in group IDs and row data, but plus and minus strands
    are not written to separate files.
    """

    if worker_count < 1:
        raise ValueError(
            "worker_count must be at least 1."
        )

    if reads_per_worker < 1:
        raise ValueError(
            "reads_per_worker must be at least 1."
        )


    parsed_decoy_reference_patterns = (
        _parse_decoy_reference_patterns(
            decoy_reference_patterns
        )
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
            / "mapping_mt.log"
        )

    logger = create_logger(
        "RiboTxAlignMT",
        log_file,
    )

    counters = MappingCounters()
    decoy_filtered_read_count = 0
    start_time = time.time()
    next_progress = progress_interval

    fastq_batch_size = (
        worker_count
        * reads_per_worker
    )

    logger.info(
        "Multiprocess transcript mapping started"
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
    logger.info(
        "Worker count: %d",
        worker_count,
    )
    logger.info(
        "Reads per worker: %d",
        reads_per_worker,
    )
    logger.info(
        "Reads per FASTQ batch: %d",
        fastq_batch_size,
    )
    logger.info(
        "Output partition: chromosome only"
    )
    logger.info(
        "Worker output: enabled"
    )
    logger.info(
        "Exclude decoy reads: %s",
        exclude_decoy_reads,
    )
    logger.info(
        "Decoy reference patterns: %s",
        ",".join(parsed_decoy_reference_patterns),
    )
    logger.info(
        "Unmapped FASTQ: %s",
        unmapped_fastq_path if unmapped_fastq_path is not None else "disabled",
    )

    unmapped_parts_directory = None
    if unmapped_fastq_path is not None:
        unmapped_destination = Path(unmapped_fastq_path)
        unmapped_destination.parent.mkdir(parents=True, exist_ok=True)
        unmapped_parts_directory = (
            unmapped_destination.parent
            / f".{unmapped_destination.name}.parts"
        )
        if unmapped_parts_directory.exists():
            shutil.rmtree(unmapped_parts_directory)
        unmapped_parts_directory.mkdir(parents=True, exist_ok=True)

    context = get_context(
        multiprocessing_start_method
    )

    if unmapped_fastq_path is not None:
        Path(unmapped_fastq_path).unlink(missing_ok=True)

    with ProcessPoolExecutor(
        max_workers=worker_count,
        mp_context=context,
        initializer=initialize_worker,
        initargs=(
            transcript_fasta,
            maximum_candidates,
            kmer_size,
            window_size,
            minimum_minimizers,
            minimum_chain_score,
            minimum_dp_score,
            minimum_query_coverage,
            minimum_identity,
        ),
    ) as executor:
        for fastq_batch_index, batch in enumerate(
            read_fastq_batches(
                fastq_path=fastq_path,
                batch_size=fastq_batch_size,
                maximum_reads=maximum_reads,
            )
        ):
            batch_start = time.time()

            chunks = split_batch(
                batch,
                worker_count,
            )

            tasks = [
                WorkerTask(
                    reads=chunk,
                    output_directory=str(
                        output_directory
                    ),
                    part_index=(
                        fastq_batch_index
                        * worker_count
                        + chunk_index
                    ),
                    parquet_row_group_size=(
                        parquet_row_group_size
                    ),
                    exclude_decoy_reads=(
                        exclude_decoy_reads
                    ),
                    decoy_reference_patterns=(
                        parsed_decoy_reference_patterns
                    ),
                    unmapped_fastq_part=(
                        str(
                            unmapped_parts_directory
                            / (
                                f"part_{fastq_batch_index:06d}_"
                                f"{chunk_index:03d}.fastq.gz"
                            )
                        )
                        if unmapped_parts_directory is not None
                        else None
                    ),
                )
                for chunk_index, chunk
                in enumerate(
                    chunks
                )
            ]

            batch_mapped_read_count = 0
            batch_decoy_filtered_read_count = 0
            batch_output_files: list[str] = []
            batch_unmapped_parts: list[str] = []

            for worker_result in executor.map(
                map_and_write_chunk,
                tasks,
            ):
                counters.add(
                    worker_result.counters
                )

                batch_mapped_read_count += (
                    worker_result.mapped_read_count
                )

                batch_decoy_filtered_read_count += (
                    worker_result.decoy_filtered_read_count
                )

                decoy_filtered_read_count += (
                    worker_result.decoy_filtered_read_count
                )

                batch_output_files.extend(worker_result.output_files)

                if worker_result.unmapped_fastq_part is not None:
                    batch_unmapped_parts.append(
                        worker_result.unmapped_fastq_part
                    )

            if unmapped_fastq_path is not None:
                with Path(unmapped_fastq_path).open("ab") as output_handle:
                    for part_path_text in batch_unmapped_parts:
                        part_path = Path(part_path_text)
                        with part_path.open("rb") as input_handle:
                            shutil.copyfileobj(input_handle, output_handle)
                        part_path.unlink(missing_ok=True)

            logger.info(
                "FASTQBatch=%d "
                "InputReads=%d "
                "MappedReads=%d "
                "DecoyFilteredReads=%d "
                "WorkerChunks=%d "
                "ChromosomePartFiles=%d "
                "Elapsed=%.1fs",
                fastq_batch_index,
                len(batch),
                batch_mapped_read_count,
                batch_decoy_filtered_read_count,
                len(tasks),
                len(batch_output_files),
                time.time() - batch_start,
            )

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
                    "DecoyFiltered=%d "
                    "BestHitLimitReached=%d",
                    counters.total_reads,
                    speed,
                    counters.mapped_reads,
                    counters.group_unique_reads,
                    counters.group_multi_reads,
                    counters.unmapped_reads,
                    decoy_filtered_read_count,
                    counters.best_hit_overflow_reads,
                )

                while (
                    next_progress
                    <= counters.total_reads
                ):
                    next_progress += (
                        progress_interval
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
        "Multiprocess transcript mapping finished"
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
        "Decoy-filtered reads: %d",
        decoy_filtered_read_count,
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

    if stats_json_path is not None:
        analyzable_reads = (
            counters.total_reads
            - decoy_filtered_read_count
        )
        stats = {
            "input_reads": counters.total_reads,
            "transcript_mapped": counters.mapped_reads,
            "group_unique": counters.group_unique_reads,
            "group_multi": counters.group_multi_reads,
            "transcript_unmapped": counters.unmapped_reads,
            "decoy_filtered": decoy_filtered_read_count,
            "best_hit_limit_reached": counters.best_hit_overflow_reads,
            "analyzable_reads": analyzable_reads,
            "transcript_mapped_rate_percent": percentage(
                counters.mapped_reads, counters.total_reads
            ),
            "transcript_unique_rate_percent": percentage(
                counters.group_unique_reads, counters.total_reads
            ),
            "transcript_multi_rate_percent": percentage(
                counters.group_multi_reads, counters.total_reads
            ),
            "transcript_unmapped_rate_percent": percentage(
                counters.unmapped_reads, counters.total_reads
            ),
            "decoy_filtered_rate_percent": percentage(
                decoy_filtered_read_count, counters.total_reads
            ),
            "transcript_mapped_non_decoy_rate_percent": percentage(
                counters.mapped_reads, analyzable_reads
            ),
        }
        write_mapping_stats(stats_json_path, stats)
        logger.info(
            "Mapping statistics JSON: %s",
            stats_json_path,
        )

    return counters