from __future__ import annotations

import logging
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, fields
from multiprocessing import get_context
from pathlib import Path
from typing import Iterator

import mappy as mp

from mapping.bam_writer import GenomeBamWriter
from mapping.cigar_utils import add_soft_clips
from mapping.coordinate_converter import CoordinateConverter
from mapping.mappy_mapper import MappyMapper
from mapping.mapping_result import ProjectedAlignment, RawAlignment
from mapping.parquet_writer import MappingParquetWriter


_WORKER_GENOME_MAPPER: MappyMapper | None = None
_WORKER_TRANSCRIPT_MAPPER: MappyMapper | None = None
_WORKER_CONVERTER: CoordinateConverter | None = None
_WORKER_TRANSCRIPT_CANDIDATE_COUNT = 4


@dataclass
class ReadRecord:
    read_id: str
    sequence: str
    quality: str


@dataclass
class MappedRead:
    read_id: str
    sequence: str
    quality: str
    alignments: list[ProjectedAlignment]
    genome_status: str
    genome_hit_count: int
    transcript_over_candidate_limit: bool = False


@dataclass
class MappingCounters:
    read_count: int = 0
    genome_unique_count: int = 0
    transcript_single_count: int = 0
    transcript_multi_count: int = 0
    transcript_multihit_2_count: int = 0
    transcript_multihit_3_count: int = 0
    transcript_multihit_4_count: int = 0
    transcript_multihit_5_count: int = 0
    transcript_multihit_6_count: int = 0
    transcript_multihit_7_count: int = 0
    transcript_multihit_8_count: int = 0
    transcript_multihit_over8_count: int = 0
    unmapped_count: int = 0
    conversion_errors: int = 0
    genome_multi_transcript_unmapped_count: int = 0
    genome_softclip_transcript_unmapped_count: int = 0
    genome_unmapped_transcript_unmapped_count: int = 0
    genome_no_hit_count: int = 0
    genome_single_hit_count: int = 0
    genome_equal_score_count: int = 0
    genome_unequal_score_count: int = 0
    genome_best_not_full_length_count: int = 0

    def add(self, other: "MappingCounters") -> None:
        for field_info in fields(self):
            field_name = field_info.name
            setattr(
                self,
                field_name,
                getattr(self, field_name)
                + getattr(other, field_name),
            )


@dataclass
class ChunkResult:
    mapped_reads: list[MappedRead]
    counters: MappingCounters
    conversion_error_messages: list[str]


def project_genome_alignment(
    alignment: RawAlignment,
    read_length: int,
) -> ProjectedAlignment:
    is_reverse = alignment.strand == -1

    return ProjectedAlignment(
        read_id=alignment.read_id,
        chrom=alignment.reference_name,
        reference_start=alignment.reference_start,
        cigartuples=add_soft_clips(
            cigartuples=alignment.cigartuples,
            query_start=alignment.query_start,
            query_end=alignment.query_end,
            read_length=read_length,
            is_reverse=is_reverse,
        ),
        is_reverse=is_reverse,
        mapq=alignment.mapq,
        edit_distance=alignment.edit_distance,
        score=alignment.score,
        source_reference=alignment.reference_name,
        source_reference_start=alignment.reference_start,
        source_reference_end=alignment.reference_end,
        support_sources={"genome"},
        assignment_weight=1.0,
    )


def classify_genome_hits(
    hits: list[RawAlignment],
    read_length: int,
) -> tuple[str, RawAlignment | None]:
    if not hits:
        return "unmapped", None

    best_hit = hits[0]

    if not best_hit.is_full_length(read_length):
        return "soft_clipped", None

    if len(hits) == 1:
        return "unique", best_hit

    if best_hit.score == hits[1].score:
        return "multi_equal_score", None

    return "unique", best_hit


def select_equal_score_transcript_hits(
    hits: list[RawAlignment],
    read_length: int,
    maximum_candidates: int | None = None,
) -> list[RawAlignment]:
    full_length_hits = [
        hit
        for hit in hits
        if hit.is_full_length(read_length)
    ]

    if not full_length_hits:
        return []

    best_score = full_length_hits[0].score

    equal_score_hits = [
        hit
        for hit in full_length_hits
        if hit.score == best_score
    ]

    if maximum_candidates is None:
        return equal_score_hits

    return equal_score_hits[:maximum_candidates]


def deduplicate_projected_alignments(
    alignments: list[ProjectedAlignment],
) -> list[ProjectedAlignment]:
    deduplicated: dict[tuple, ProjectedAlignment] = {}

    for alignment in alignments:
        key = alignment.deduplication_key

        if key not in deduplicated:
            deduplicated[key] = alignment
        else:
            deduplicated[key].support_sources.update(
                alignment.support_sources
            )

    return list(deduplicated.values())


def create_logger(log_file: str) -> logging.Logger:
    logger = logging.getLogger("RiboTxAlignMT")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()

    formatter = logging.Formatter(
        "%(asctime)s\t%(levelname)s\t%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(
        log_file,
        mode="w",
        encoding="utf-8",
    )
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


def initialize_mapping_worker(
    genome_fasta: str,
    transcript_fasta: str,
    coordinate_tsv: str,
    transcript_candidate_count: int,
    kmer_size: int,
    window_size: int,
    minimum_minimizers: int,
    minimum_chain_score: int,
    minimum_dp_score: int,
    minimum_query_coverage: float,
    minimum_identity: float,
) -> None:
    """
    Initialize one independent mapper pair in each process.
    """

    global _WORKER_GENOME_MAPPER
    global _WORKER_TRANSCRIPT_MAPPER
    global _WORKER_CONVERTER
    global _WORKER_TRANSCRIPT_CANDIDATE_COUNT

    _WORKER_GENOME_MAPPER = MappyMapper(
        reference_path=genome_fasta,
        source="genome",
        preset="sr",
        kmer_size=kmer_size,
        window_size=window_size,
        minimum_minimizers=minimum_minimizers,
        minimum_chain_score=minimum_chain_score,
        minimum_dp_score=minimum_dp_score,
        maximum_hits=2,
        index_threads=1,
        minimum_query_coverage=minimum_query_coverage,
        minimum_identity=minimum_identity,
    )

    _WORKER_TRANSCRIPT_MAPPER = MappyMapper(
        reference_path=transcript_fasta,
        source="transcript",
        preset="sr",
        kmer_size=kmer_size,
        window_size=window_size,
        minimum_minimizers=minimum_minimizers,
        minimum_chain_score=minimum_chain_score,
        minimum_dp_score=minimum_dp_score,
        maximum_hits=max(32, transcript_candidate_count * 16),
        index_threads=1,
        minimum_query_coverage=minimum_query_coverage,
        minimum_identity=minimum_identity,
    )

    _WORKER_CONVERTER = CoordinateConverter(
        coordinate_tsv
    )

    _WORKER_TRANSCRIPT_CANDIDATE_COUNT = (
        transcript_candidate_count
    )


def _require_worker_objects() -> tuple[
    MappyMapper,
    MappyMapper,
    CoordinateConverter,
]:
    if _WORKER_GENOME_MAPPER is None:
        raise RuntimeError(
            "Genome mapper is not initialized."
        )

    if _WORKER_TRANSCRIPT_MAPPER is None:
        raise RuntimeError(
            "Transcript mapper is not initialized."
        )

    if _WORKER_CONVERTER is None:
        raise RuntimeError(
            "Coordinate converter is not initialized."
        )

    return (
        _WORKER_GENOME_MAPPER,
        _WORKER_TRANSCRIPT_MAPPER,
        _WORKER_CONVERTER,
    )


def map_read_chunk(
    reads: list[ReadRecord],
) -> ChunkResult:
    """
    Map one read chunk inside one worker process.
    """

    (
        genome_mapper,
        transcript_mapper,
        converter,
    ) = _require_worker_objects()

    counters = MappingCounters()
    mapped_reads: list[MappedRead] = []
    conversion_error_messages: list[str] = []

    for read in reads:
        counters.read_count += 1
        read_length = len(read.sequence)

        # Reset the per-read flag before any mapping branch.
        transcript_over_candidate_limit = False

        genome_hits = genome_mapper.map_read(
            read_id=read.read_id,
            sequence=read.sequence,
        )

        (
            genome_status,
            genome_best_hit,
        ) = classify_genome_hits(
            hits=genome_hits,
            read_length=read_length,
        )

        if not genome_hits:
            counters.genome_no_hit_count += 1

        else:
            best_hit = genome_hits[0]

            if not best_hit.is_full_length(
                read_length
            ):
                counters.genome_best_not_full_length_count += 1

            elif len(genome_hits) == 1:
                counters.genome_single_hit_count += 1

            elif (
                genome_hits[0].score
                == genome_hits[1].score
            ):
                counters.genome_equal_score_count += 1

            else:
                counters.genome_unequal_score_count += 1

        if (
            genome_status == "unique"
            and genome_best_hit is not None
        ):
            counters.genome_unique_count += 1

            output_alignments = [
                project_genome_alignment(
                    alignment=genome_best_hit,
                    read_length=read_length,
                )
            ]

        else:
            transcript_hits = transcript_mapper.map_read(
                read_id=read.read_id,
                sequence=read.sequence,
            )

            transcript_candidates = (
                select_equal_score_transcript_hits(
                    hits=transcript_hits,
                    read_length=read_length,
                    maximum_candidates=None,
                )
            )

            output_alignments: list[
                ProjectedAlignment
            ] = []

            for transcript_hit in transcript_candidates:
                try:
                    output_alignments.append(
                        converter.project_alignment(
                            alignment=transcript_hit,
                            read_length=read_length,
                        )
                    )

                except (
                    KeyError,
                    ValueError,
                ) as error:
                    counters.conversion_errors += 1

                    if len(
                        conversion_error_messages
                    ) < 10:
                        conversion_error_messages.append(
                            (
                                "Coordinate conversion skipped: "
                                f"read={read.read_id} "
                                "reference="
                                f"{transcript_hit.reference_name} "
                                f"error={error}"
                            )
                        )

            deduplicated_alignments = (
                deduplicate_projected_alignments(
                    output_alignments
                )
            )
            candidate_count = len(deduplicated_alignments)

            if candidate_count == 0:
                counters.unmapped_count += 1

                if (
                    genome_status
                    == "multi_equal_score"
                ):
                    counters.genome_multi_transcript_unmapped_count += 1

                elif genome_status == "soft_clipped":
                    counters.genome_softclip_transcript_unmapped_count += 1

                elif genome_status == "unmapped":
                    counters.genome_unmapped_transcript_unmapped_count += 1

                output_alignments = []

            elif candidate_count == 1:
                counters.transcript_single_count += 1
                output_alignments = deduplicated_alignments
                output_alignments[
                    0
                ].assignment_weight = 1.0

            elif candidate_count > _WORKER_TRANSCRIPT_CANDIDATE_COUNT:
                # Reads with more than the configured number of genomic
                # candidates are recorded as unresolved multi-mappers.
                # They remain in multiplemap.bam but are not written to the
                # EM candidate Parquet.
                counters.transcript_multi_count += 1
                counters.transcript_multihit_over8_count += 1
                transcript_over_candidate_limit = True
                output_alignments = []

            else:
                counters.transcript_multi_count += 1
                output_alignments = deduplicated_alignments

                if candidate_count == 2:
                    counters.transcript_multihit_2_count += 1
                elif candidate_count == 3:
                    counters.transcript_multihit_3_count += 1
                elif candidate_count == 4:
                    counters.transcript_multihit_4_count += 1
                elif candidate_count == 5:
                    counters.transcript_multihit_5_count += 1
                elif candidate_count == 6:
                    counters.transcript_multihit_6_count += 1
                elif candidate_count == 7:
                    counters.transcript_multihit_7_count += 1
                elif candidate_count == 8:
                    counters.transcript_multihit_8_count += 1

                # Multi-mapping candidates are kept unweighted here.
                # Their final fractional weights will be estimated later.
                equal_weight = 1.0 / candidate_count
                for alignment in output_alignments:
                    alignment.assignment_weight = equal_weight

        mapped_reads.append(
            MappedRead(
                read_id=read.read_id,
                sequence=read.sequence,
                quality=read.quality,
                alignments=output_alignments,
                genome_status=genome_status,
                genome_hit_count=len(genome_hits),
                transcript_over_candidate_limit=(
                    transcript_over_candidate_limit
                ),
            )
        )

    return ChunkResult(
        mapped_reads=mapped_reads,
        counters=counters,
        conversion_error_messages=(
            conversion_error_messages
        ),
    )


def read_fastq_batches(
    fastq_path: str,
    batch_size: int,
    maximum_reads: int | None,
) -> Iterator[list[ReadRecord]]:
    """
    Read FASTQ records in bounded batches.
    """

    batch: list[ReadRecord] = []
    total_read_count = 0

    for (
        read_id,
        sequence,
        quality,
    ) in mp.fastx_read(
        fastq_path
    ):
        batch.append(
            ReadRecord(
                read_id=read_id,
                sequence=sequence,
                quality=quality,
            )
        )

        total_read_count += 1

        if len(batch) >= batch_size:
            yield batch
            batch = []

        if (
            maximum_reads is not None
            and total_read_count >= maximum_reads
        ):
            break

    if batch:
        yield batch


def split_batch(
    reads: list[ReadRecord],
    worker_count: int,
) -> list[list[ReadRecord]]:
    """
    Split one batch into approximately equal chunks.
    """

    if not reads:
        return []

    actual_worker_count = min(
        worker_count,
        len(reads),
    )

    chunk_size = (
        len(reads)
        + actual_worker_count
        - 1
    ) // actual_worker_count

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


def run_dual_mapping_to_bam_mt(
    fastq_path: str,
    genome_fasta: str,
    transcript_fasta: str,
    coordinate_tsv: str,
    unique_bam: str,
    multiplemap_bam: str,
    mapping_parquet: str,
    maximum_reads: int | None = None,
    transcript_candidate_count: int = 8,
    kmer_size: int = 11,
    window_size: int = 5,
    minimum_minimizers: int = 2,
    minimum_chain_score: int = 15,
    minimum_dp_score: int = 15,
    minimum_query_coverage: float = 0.85,
    minimum_identity: float = 0.90,
    worker_count: int = 4,
    reads_per_worker: int = 10000,
    sort_threads: int = 8,
    progress_interval: int = 100000,
    parquet_batch_size: int = 100000,
    log_file: str | None = None,
    multiprocessing_start_method: str = "spawn",
) -> None:
    """Map reads and separate unique and multiple mappings.

    Unique reads are written with their genomic alignment to ``unique_bam``.
    Reads with two or more transcript candidates are written once as unmapped
    records to ``multiplemap_bam``. Their candidate mapping information is
    stored separately in ``mapping_parquet``.
    """

    if worker_count < 1:
        raise ValueError("worker_count must be at least 1.")
    if reads_per_worker < 1:
        raise ValueError("reads_per_worker must be at least 1.")

    if log_file is None:
        log_file = str(Path(unique_bam).with_suffix(".mapping_mt.log"))

    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    logger = create_logger(log_file)
    start_time = time.time()
    total_counters = MappingCounters()
    displayed_read_count = 0
    next_progress_count = progress_interval
    batch_size = worker_count * reads_per_worker

    logger.info("Multiprocess mapping started")
    logger.info("FASTQ: %s", fastq_path)
    logger.info("Genome FASTA: %s", genome_fasta)
    logger.info("Transcript FASTA: %s", transcript_fasta)
    logger.info("Coordinate TSV: %s", coordinate_tsv)
    logger.info("Unique BAM: %s", unique_bam)
    logger.info("Multiple-map BAM: %s", multiplemap_bam)
    logger.info("Mapping Parquet: %s", mapping_parquet)
    logger.info("Worker count: %d", worker_count)
    logger.info("Reads per worker: %d", reads_per_worker)
    logger.info("Batch size: %d", batch_size)
    logger.info("Maximum reads: %s", "all" if maximum_reads is None else maximum_reads)

    mp_context = get_context(multiprocessing_start_method)

    with ProcessPoolExecutor(
        max_workers=worker_count,
        mp_context=mp_context,
        initializer=initialize_mapping_worker,
        initargs=(
            genome_fasta,
            transcript_fasta,
            coordinate_tsv,
            transcript_candidate_count,
            kmer_size,
            window_size,
            minimum_minimizers,
            minimum_chain_score,
            minimum_dp_score,
            minimum_query_coverage,
            minimum_identity,
        ),
    ) as executor:
        with (
            GenomeBamWriter(
                genome_fasta=genome_fasta,
                output_bam=unique_bam,
                sort_threads=sort_threads,
            ) as unique_writer,
            GenomeBamWriter(
                genome_fasta=genome_fasta,
                output_bam=multiplemap_bam,
                sort_threads=sort_threads,
            ) as multiple_writer,
            MappingParquetWriter(
                output_parquet=mapping_parquet,
                batch_size=parquet_batch_size,
            ) as parquet_writer,
        ):
            for batch in read_fastq_batches(
                fastq_path=fastq_path,
                batch_size=batch_size,
                maximum_reads=maximum_reads,
            ):
                chunks = split_batch(reads=batch, worker_count=worker_count)
                chunk_results = executor.map(map_read_chunk, chunks)

                for chunk_result in chunk_results:
                    total_counters.add(chunk_result.counters)

                    for message in chunk_result.conversion_error_messages:
                        if total_counters.conversion_errors <= 10:
                            logger.warning(message)

                    for mapped_read in chunk_result.mapped_reads:
                        hit_count = len(mapped_read.alignments)

                        if hit_count == 1:
                            unique_writer.write_read(
                                read_id=mapped_read.read_id,
                                sequence=mapped_read.sequence,
                                quality=mapped_read.quality,
                                alignments=mapped_read.alignments,
                            )
                        elif hit_count >= 2:
                            # Store the read sequence once, without choosing a
                            # genomic placement before EM allocation.
                            multiple_writer.write_read(
                                read_id=mapped_read.read_id,
                                sequence=mapped_read.sequence,
                                quality=mapped_read.quality,
                                alignments=[],
                            )
                            parquet_writer.write_read(
                                read_id=mapped_read.read_id,
                                read_length=len(mapped_read.sequence),
                                genome_status=mapped_read.genome_status,
                                genome_hit_count=mapped_read.genome_hit_count,
                                alignments=mapped_read.alignments,
                            )
                        elif mapped_read.transcript_over_candidate_limit:
                            # Preserve >N-hit reads as unresolved records, but
                            # do not expose a truncated candidate set to EM.
                            multiple_writer.write_read(
                                read_id=mapped_read.read_id,
                                sequence=mapped_read.sequence,
                                quality=mapped_read.quality,
                                alignments=[],
                            )

                        if displayed_read_count < 20:
                            logger.info(
                                "Read=%s length=%d genome_status=%s genome_hits=%d output_hits=%d",
                                mapped_read.read_id,
                                len(mapped_read.sequence),
                                mapped_read.genome_status,
                                mapped_read.genome_hit_count,
                                hit_count,
                            )
                            displayed_read_count += 1

                if progress_interval > 0 and total_counters.read_count >= next_progress_count:
                    elapsed_seconds = time.time() - start_time
                    reads_per_second = (
                        total_counters.read_count / elapsed_seconds
                        if elapsed_seconds > 0 else 0.0
                    )
                    logger.info(
                        "Processed=%d Speed=%.1f reads/s GenomeUnique=%d "
                        "TranscriptSingle=%d TranscriptMulti=%d "
                        "MultiHit2=%d MultiHit3=%d MultiHit4=%d "
                        "MultiHit5=%d MultiHit6=%d MultiHit7=%d MultiHit8=%d MultiHitOver8=%d "
                        "Unmapped=%d",
                        total_counters.read_count,
                        reads_per_second,
                        total_counters.genome_unique_count,
                        total_counters.transcript_single_count,
                        total_counters.transcript_multi_count,
                        total_counters.transcript_multihit_2_count,
                        total_counters.transcript_multihit_3_count,
                        total_counters.transcript_multihit_4_count,
                        total_counters.transcript_multihit_5_count,
                        total_counters.transcript_multihit_6_count,
                        total_counters.transcript_multihit_7_count,
                        total_counters.transcript_multihit_8_count,
                        total_counters.transcript_multihit_over8_count,
                        total_counters.unmapped_count,
                    )
                    while next_progress_count <= total_counters.read_count:
                        next_progress_count += progress_interval

    elapsed_seconds = time.time() - start_time
    reads_per_second = (
        total_counters.read_count / elapsed_seconds
        if elapsed_seconds > 0 else 0.0
    )
    logger.info("Multiprocess mapping finished")
    logger.info("Reads: %d", total_counters.read_count)
    logger.info("Genome unique: %d", total_counters.genome_unique_count)
    logger.info("Transcript single candidate: %d", total_counters.transcript_single_count)
    logger.info("Transcript multiple candidates: %d", total_counters.transcript_multi_count)
    logger.info("TranscriptMulti multihit=2: %d", total_counters.transcript_multihit_2_count)
    logger.info("TranscriptMulti multihit=3: %d", total_counters.transcript_multihit_3_count)
    logger.info("TranscriptMulti multihit=4: %d", total_counters.transcript_multihit_4_count)
    logger.info("TranscriptMulti multihit=5: %d", total_counters.transcript_multihit_5_count)
    logger.info("TranscriptMulti multihit=6: %d", total_counters.transcript_multihit_6_count)
    logger.info("TranscriptMulti multihit=7: %d", total_counters.transcript_multihit_7_count)
    logger.info("TranscriptMulti multihit=8: %d", total_counters.transcript_multihit_8_count)
    logger.info("TranscriptMulti multihit>8: %d", total_counters.transcript_multihit_over8_count)
    logger.info("Unmapped: %d", total_counters.unmapped_count)
    logger.info("Coordinate conversion errors: %d", total_counters.conversion_errors)
    logger.info("Elapsed seconds: %.1f", elapsed_seconds)
    logger.info("Average speed: %.1f reads/s", reads_per_second)
