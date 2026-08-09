from __future__ import annotations

import csv
import logging
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pyarrow as pa
import pyarrow.parquet as pq

from depth_output.transcript_genome_coordinate import (
    TranscriptGenomeCoordinateConverter,
)


@dataclass(slots=True)
class CandidateAlignment:
    """One read candidate after transcript-to-genome projection."""

    read_id: str
    transcript_id: str
    alignment_transcript_id: str
    group_id: str
    chrom: str
    genome_start: int
    genome_end: int
    strand: str
    genome_cigar: str
    posterior_probability: float
    transcript_start: int
    transcript_end: int
    read_length: int
    hit_index: int
    query_start: int
    query_end: int
    mapping_strand: int
    transcript_cigar: str

    @property
    def genome_signature(self) -> tuple[str, int, str, str]:
        """Return the BAM-relevant alignment identity."""

        return (
            self.chrom,
            self.genome_start,
            self.genome_cigar,
            self.strand,
        )


@dataclass(slots=True)
class ReadCandidateSet:
    """All collapsed genomic candidates for one read."""

    read_id: str
    candidates: list[CandidateAlignment]

    @property
    def candidate_count(self) -> int:
        """Return the number of distinct genomic alignments."""

        return len(self.candidates)


def _create_logger(
    logger: logging.Logger | None,
) -> logging.Logger:
    """Return the provided logger or create a module logger."""

    if logger is not None:
        return logger

    result = logging.getLogger("RiboTxAlignAssignment")
    result.setLevel(logging.INFO)

    if not result.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter(
                "%(asctime)s\t%(levelname)s\t%(message)s"
            )
        )
        result.addHandler(handler)

    return result


def _safe_float(value: object) -> float:
    """Convert a possibly missing numeric value to float."""

    if value is None:
        return 0.0

    try:
        result = float(value)
    except (TypeError, ValueError):
        return 0.0

    if not math.isfinite(result):
        return 0.0

    return max(result, 0.0)


def _safe_int(value: object, default: int = 0) -> int:
    """Convert a possibly missing numeric value to int."""

    if value is None:
        return default

    try:
        return int(value)
    except (TypeError, ValueError):
        return default


_CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
_REFERENCE_CONSUMING_OPERATIONS = {"M", "D", "N", "=", "X"}
_QUERY_CONSUMING_OPERATIONS = {"M", "I", "S", "=", "X"}


def _parse_cigar(cigar: str) -> list[tuple[str, int]]:
    """Parse a compact CIGAR string into operation-length pairs."""

    if not cigar or cigar == "*":
        raise ValueError("Mapped alignment requires a non-empty CIGAR.")

    matches = _CIGAR_PATTERN.findall(cigar)

    if not matches:
        raise ValueError(f"Invalid CIGAR string: {cigar}")

    reconstructed = "".join(
        f"{length}{operation}"
        for length, operation in matches
    )

    if reconstructed != cigar:
        raise ValueError(f"Invalid CIGAR string: {cigar}")

    return [
        (operation, int(length))
        for length, operation in matches
    ]


def _merge_adjacent_cigar_operations(
    operations: list[tuple[str, int]],
) -> list[tuple[str, int]]:
    """Merge adjacent CIGAR operations of the same type."""

    merged: list[tuple[str, int]] = []

    for operation, length in operations:
        if length <= 0:
            continue

        if merged and merged[-1][0] == operation:
            previous_operation, previous_length = merged[-1]
            merged[-1] = (
                previous_operation,
                previous_length + length,
            )
        else:
            merged.append((operation, length))

    return merged


def _format_cigar(
    operations: list[tuple[str, int]],
) -> str:
    """Format operation-length pairs as one compact CIGAR string."""

    return "".join(
        f"{length}{operation}"
        for operation, length in _merge_adjacent_cigar_operations(operations)
    )


def _query_consumed_length(
    operations: list[tuple[str, int]],
) -> int:
    """Return the number of query bases consumed by CIGAR operations."""

    return sum(
        length
        for operation, length in operations
        if operation in _QUERY_CONSUMING_OPERATIONS
    )


def _reference_consumed_length(
    operations: list[tuple[str, int]],
) -> int:
    """Return the number of transcript-reference bases consumed."""

    return sum(
        length
        for operation, length in operations
        if operation in _REFERENCE_CONSUMING_OPERATIONS
    )


def _add_soft_clips(
    operations: list[tuple[str, int]],
    read_length: int,
    query_start: int,
    query_end: int,
    mapping_strand: int,
) -> list[tuple[str, int]]:
    """Add unaligned query ends as soft clips in target-coordinate order."""

    if not (0 <= query_start <= query_end <= read_length):
        raise ValueError(
            "Invalid query interval: "
            f"start={query_start}, end={query_end}, "
            f"read_length={read_length}"
        )

    consumed = _query_consumed_length(operations)

    # Some aligners may already include complete soft clipping.
    if consumed == read_length:
        return operations

    aligned_query_length = query_end - query_start

    if consumed != aligned_query_length:
        raise ValueError(
            "Original CIGAR query length does not match the mapped query "
            f"interval: cigar={_format_cigar(operations)}, "
            f"cigar_query_length={consumed}, "
            f"query_interval_length={aligned_query_length}"
        )

    if mapping_strand >= 0:
        left_soft_clip = query_start
        right_soft_clip = read_length - query_end
    else:
        # The CIGAR is ordered along increasing transcript coordinates.
        left_soft_clip = read_length - query_end
        right_soft_clip = query_start

    result: list[tuple[str, int]] = []

    if left_soft_clip > 0:
        result.append(("S", left_soft_clip))

    result.extend(operations)

    if right_soft_clip > 0:
        result.append(("S", right_soft_clip))

    return result


def _find_segment_index(
    segments,
    reference_position: int,
) -> int:
    """Return the segment containing one transcript reference position."""

    for index, segment in enumerate(segments):
        if (
            segment.reference_start
            <= reference_position
            < segment.reference_end
        ):
            return index

    raise ValueError(
        f"Transcript coordinate is outside the reference: "
        f"{reference_position}"
    )


def _genomic_intron_length(
    left_segment,
    right_segment,
) -> int:
    """Return genomic skipped length between adjacent transcript segments."""

    if left_segment.strand == "+":
        intron_length = (
            right_segment.genome_start
            - left_segment.genome_end
        )
    else:
        intron_length = (
            left_segment.genome_start
            - right_segment.genome_end
        )

    if intron_length < 0:
        raise ValueError(
            "Transcript segments overlap in genomic coordinates: "
            f"{left_segment.reference_id}"
        )

    return intron_length


def build_genome_cigar(
    converter: TranscriptGenomeCoordinateConverter,
    transcript_id: str,
    reference_start: int,
    transcript_cigar: str,
    read_length: int,
    query_start: int,
    query_end: int,
    mapping_strand: int,
) -> tuple[str, int, int, str, str] | None:
    """Project a transcript alignment CIGAR onto genomic exon coordinates.

    Original match, mismatch, insertion, deletion, clipping, and padding
    operations are retained. Genomic N operations are inserted only when a
    reference-consuming operation crosses a transcript exon boundary.
    """

    segments = converter.get_segments(transcript_id)

    if not segments:
        return None

    operations = _add_soft_clips(
        operations=_parse_cigar(transcript_cigar),
        read_length=read_length,
        query_start=query_start,
        query_end=query_end,
        mapping_strand=mapping_strand,
    )

    reference_consumed = _reference_consumed_length(operations)

    if reference_consumed <= 0:
        return None

    reference_end = reference_start + reference_consumed
    transcript_length = segments[-1].reference_end

    if reference_start < 0 or reference_end > transcript_length:
        raise ValueError(
            f"Alignment exceeds transcript reference {transcript_id}: "
            f"{reference_start}-{reference_end} / {transcript_length}"
        )

    projected_operations: list[tuple[str, int]] = []
    reference_position = reference_start
    genomic_blocks: list[tuple[int, int]] = []

    for operation, operation_length in operations:
        if operation not in _REFERENCE_CONSUMING_OPERATIONS:
            projected_operations.append((operation, operation_length))
            continue

        remaining = operation_length

        while remaining > 0:
            segment_index = _find_segment_index(
                segments,
                reference_position,
            )
            segment = segments[segment_index]
            available = segment.reference_end - reference_position
            chunk_length = min(remaining, available)

            projected_operations.append((operation, chunk_length))

            relative_start = reference_position - segment.reference_start
            relative_end = relative_start + chunk_length

            if segment.strand == "+":
                genome_block_start = segment.genome_start + relative_start
                genome_block_end = segment.genome_start + relative_end
            else:
                genome_block_start = segment.genome_end - relative_end
                genome_block_end = segment.genome_end - relative_start

            genomic_blocks.append(
                (genome_block_start, genome_block_end)
            )

            reference_position += chunk_length
            remaining -= chunk_length

            if (
                remaining > 0
                and reference_position == segment.reference_end
            ):
                next_segment_index = segment_index + 1

                if next_segment_index >= len(segments):
                    raise ValueError(
                        f"CIGAR exceeds transcript reference {transcript_id}."
                    )

                next_segment = segments[next_segment_index]
                intron_length = _genomic_intron_length(
                    segment,
                    next_segment,
                )

                if intron_length > 0:
                    projected_operations.append(
                        ("N", intron_length)
                    )

    transcript_strand = segments[0].strand

    # Transcript coordinates increase opposite to genomic coordinates on a
    # minus-strand transcript. BAM CIGAR must always be ordered by increasing
    # genomic reference coordinate.
    if transcript_strand == "-":
        projected_operations = list(
            reversed(projected_operations)
        )

    genome_is_reverse = (
        (transcript_strand == "-")
        != (mapping_strand < 0)
    )
    genome_strand = "-" if genome_is_reverse else "+"

    genome_start = min(start for start, _ in genomic_blocks)
    genome_end = max(end for _, end in genomic_blocks)
    chroms = {segment.chrom for segment in segments}

    if len(chroms) != 1:
        raise ValueError(
            f"Inconsistent chromosomes for transcript {transcript_id}."
        )

    genome_cigar = _format_cigar(projected_operations)

    if _query_consumed_length(projected_operations) != read_length:
        raise ValueError(
            "Projected genome CIGAR does not consume the full query: "
            f"read_length={read_length}, cigar={genome_cigar}, "
            f"query_consumed={_query_consumed_length(projected_operations)}"
        )

    return (
        next(iter(chroms)),
        genome_start,
        genome_end,
        genome_strand,
        genome_cigar,
    )


def iter_reads_info_rows(
    reads_info_parquet: str | Path,
    batch_size: int,
) -> Iterable[dict]:
    """Yield read-candidate rows from one Parquet file or directory."""

    input_path = Path(reads_info_parquet)

    if input_path.is_dir():
        parquet_paths = sorted(input_path.glob("*.parquet"))
    else:
        parquet_paths = [input_path]

    if not parquet_paths:
        raise FileNotFoundError(
            f"No reads_info Parquet files found: {input_path}"
        )

    required_columns = {
        "read_id",
        "group_id",
        "transcript_id",
        "alignment_transcript_id",
        "start",
        "end",
        "prior_probability",
        "transcript_probability",
        "hit_index",
        "query_start",
        "query_end",
        "mapping_strand",
        "cigar",
    }

    for parquet_path in parquet_paths:
        parquet_file = pq.ParquetFile(parquet_path)
        schema_names = set(parquet_file.schema_arrow.names)
        missing = sorted(required_columns - schema_names)

        if missing:
            raise ValueError(
                f"{parquet_path} is missing columns: {', '.join(missing)}"
            )

        columns = [
            "read_id",
            "group_id",
            "transcript_id",
            "alignment_transcript_id",
            "start",
            "end",
            "prior_probability",
            "transcript_probability",
            "hit_index",
            "query_start",
            "query_end",
            "mapping_strand",
            "cigar",
        ]

        if "read_length" in schema_names:
            columns.append("read_length")

        for batch in parquet_file.iter_batches(
            batch_size=batch_size,
            columns=columns,
        ):
            yield from batch.to_pylist()


def load_and_project_candidates(
    reads_info_parquet: str | Path,
    coordinate_table: str | Path,
    batch_size: int,
    logger: logging.Logger,
) -> dict[str, list[CandidateAlignment]]:
    """Load candidate rows and project them to genome coordinates."""

    converter = TranscriptGenomeCoordinateConverter.from_coordinate_table(
        coordinate_table
    )

    by_read: dict[str, list[CandidateAlignment]] = defaultdict(list)
    input_row_count = 0
    projected_count = 0
    missing_coordinate_count = 0

    for row in iter_reads_info_rows(
        reads_info_parquet=reads_info_parquet,
        batch_size=batch_size,
    ):
        input_row_count += 1

        read_id = str(row["read_id"])
        transcript_id = str(row["transcript_id"])
        alignment_transcript_id = str(
            row["alignment_transcript_id"]
        )
        reference_start = _safe_int(row["start"])
        reference_end = _safe_int(row["end"])

        read_length = _safe_int(
            row.get("read_length"),
            default=max(reference_end - reference_start, 0),
        )
        hit_index = _safe_int(row.get("hit_index"), default=-1)
        query_start = _safe_int(row.get("query_start"), default=0)
        query_end = _safe_int(
            row.get("query_end"),
            default=read_length,
        )
        mapping_strand = _safe_int(
            row.get("mapping_strand"),
            default=1,
        )
        transcript_cigar = str(row["cigar"])

        projected = build_genome_cigar(
            converter=converter,
            transcript_id=alignment_transcript_id,
            reference_start=reference_start,
            transcript_cigar=transcript_cigar,
            read_length=read_length,
            query_start=query_start,
            query_end=query_end,
            mapping_strand=mapping_strand,
        )

        if projected is None:
            missing_coordinate_count += 1
            continue

        (
            chrom,
            genome_start,
            genome_end,
            strand,
            genome_cigar,
        ) = projected

        posterior = (
            _safe_float(row["prior_probability"])
            * _safe_float(row["transcript_probability"])
        )

        by_read[read_id].append(
            CandidateAlignment(
                read_id=read_id,
                transcript_id=transcript_id,
                alignment_transcript_id=alignment_transcript_id,
                group_id=str(row["group_id"]),
                chrom=chrom,
                genome_start=genome_start,
                genome_end=genome_end,
                strand=strand,
                genome_cigar=genome_cigar,
                posterior_probability=posterior,
                transcript_start=reference_start,
                transcript_end=reference_end,
                read_length=read_length,
                hit_index=hit_index,
                query_start=query_start,
                query_end=query_end,
                mapping_strand=mapping_strand,
                transcript_cigar=transcript_cigar,
            )
        )
        projected_count += 1

    logger.info("Input candidate rows: %d", input_row_count)
    logger.info("Projected candidate rows: %d", projected_count)
    logger.info(
        "Candidates without coordinate mapping: %d",
        missing_coordinate_count,
    )
    logger.info("Reads with projected candidates: %d", len(by_read))

    return by_read


def collapse_same_genome_alignments(
    candidates_by_read: dict[str, list[CandidateAlignment]],
) -> list[ReadCandidateSet]:
    """Collapse transcript candidates producing the same genomic alignment."""

    result: list[ReadCandidateSet] = []

    for read_id, candidates in candidates_by_read.items():
        collapsed: dict[
            tuple[str, int, str, str],
            CandidateAlignment,
        ] = {}

        transcript_ids_by_signature: dict[
            tuple[str, int, str, str],
            list[str],
        ] = defaultdict(list)

        group_ids_by_signature: dict[
            tuple[str, int, str, str],
            list[str],
        ] = defaultdict(list)

        for candidate in candidates:
            signature = candidate.genome_signature
            transcript_ids_by_signature[signature].append(
                candidate.transcript_id
            )
            group_ids_by_signature[signature].append(candidate.group_id)

            existing = collapsed.get(signature)

            if existing is None:
                collapsed[signature] = candidate
            else:
                existing.posterior_probability += (
                    candidate.posterior_probability
                )

        collapsed_candidates = []

        for signature, candidate in collapsed.items():
            candidate.transcript_id = ",".join(
                sorted(set(transcript_ids_by_signature[signature]))
            )
            candidate.group_id = ",".join(
                sorted(set(group_ids_by_signature[signature]))
            )
            collapsed_candidates.append(candidate)

        total_probability = sum(
            candidate.posterior_probability
            for candidate in collapsed_candidates
        )

        if total_probability <= 0.0:
            equal_probability = 1.0 / len(collapsed_candidates)

            for candidate in collapsed_candidates:
                candidate.posterior_probability = equal_probability
        else:
            for candidate in collapsed_candidates:
                candidate.posterior_probability /= total_probability

        collapsed_candidates.sort(
            key=lambda candidate: (
                -candidate.posterior_probability,
                candidate.chrom,
                candidate.genome_start,
                candidate.genome_cigar,
                candidate.strand,
            )
        )

        result.append(
            ReadCandidateSet(
                read_id=read_id,
                candidates=collapsed_candidates,
            )
        )

    return result


def largest_remainder_quota(
    expected_counts: dict[str, float],
    total_read_count: int,
) -> dict[str, int]:
    """Convert expected counts to integer quotas with a fixed total."""

    if total_read_count <= 0 or not expected_counts:
        return {
            transcript_id: 0
            for transcript_id in expected_counts
        }

    total_expected = sum(expected_counts.values())

    if total_expected <= 0.0:
        transcript_ids = sorted(expected_counts)
        base = total_read_count // len(transcript_ids)
        remainder = total_read_count % len(transcript_ids)

        return {
            transcript_id: base + (index < remainder)
            for index, transcript_id in enumerate(transcript_ids)
        }

    scaled = {
        transcript_id: (
            expected_count
            * total_read_count
            / total_expected
        )
        for transcript_id, expected_count in expected_counts.items()
    }

    quotas = {
        transcript_id: int(math.floor(value))
        for transcript_id, value in scaled.items()
    }

    remaining = total_read_count - sum(quotas.values())

    order = sorted(
        scaled,
        key=lambda transcript_id: (
            scaled[transcript_id] - quotas[transcript_id],
            scaled[transcript_id],
            transcript_id,
        ),
        reverse=True,
    )

    for transcript_id in order[:remaining]:
        quotas[transcript_id] += 1

    return quotas


def candidate_target_ids(
    candidate: CandidateAlignment,
) -> tuple[str, ...]:
    """Return group-level targets represented by one genomic candidate."""

    return tuple(
        group_id
        for group_id in candidate.group_id.split(",")
        if group_id
    )


def calculate_relative_loss_delta(
    assigned_count: float,
    target_count: float,
    epsilon: float = 1e-9,
) -> float:
    """Return the relative squared-error change after one assignment."""

    denominator = max(target_count, epsilon)
    current_loss = (assigned_count - target_count) ** 2 / denominator
    next_loss = (assigned_count + 1.0 - target_count) ** 2 / denominator
    return next_loss - current_loss


def assign_multimap_reads(
    multimap_reads: list[ReadCandidateSet],
    local_bin_size: int = 30,
) -> dict[str, CandidateAlignment]:
    """Assign multi-map reads while preserving group-level TPM proportions.

    Expected counts are accumulated from posterior probabilities at group
    level. Candidate selection minimizes the increase in relative squared
    error rather than the absolute count deficit, preventing highly expressed
    groups from absorbing nearly all ambiguous reads.
    """

    target_counts: dict[str, float] = defaultdict(float)

    for read_info in multimap_reads:
        for candidate in read_info.candidates:
            target_ids = candidate_target_ids(candidate)
            if not target_ids:
                continue

            probability_per_group = (
                candidate.posterior_probability / len(target_ids)
            )
            for group_id in target_ids:
                target_counts[group_id] += probability_per_group

    assigned_counts: dict[str, float] = defaultdict(float)

    ordered_reads = sorted(
        multimap_reads,
        key=lambda read_info: (
            -(
                read_info.candidates[0].posterior_probability
                - read_info.candidates[1].posterior_probability
            ),
            read_info.read_id,
        ),
    )

    selected: dict[str, CandidateAlignment] = {}

    for read_info in ordered_reads:
        best_candidate = None
        best_key = None

        for candidate in read_info.candidates:
            target_ids = candidate_target_ids(candidate)
            if not target_ids:
                continue

            fractional_assignment = 1.0 / len(target_ids)
            loss_delta = 0.0

            for group_id in target_ids:
                current = assigned_counts[group_id]
                target = target_counts.get(group_id, 0.0)
                denominator = max(target, 1e-9)
                current_loss = (current - target) ** 2 / denominator
                next_loss = (
                    current + fractional_assignment - target
                ) ** 2 / denominator
                loss_delta += next_loss - current_loss

            # Relative TPM error is primary. Posterior is a deterministic
            # tie-breaker and is not added again as a strong objective.
            selection_key = (
                loss_delta,
                -candidate.posterior_probability,
                candidate.chrom,
                candidate.genome_start,
                candidate.genome_cigar,
            )

            if best_key is None or selection_key < best_key:
                best_key = selection_key
                best_candidate = candidate

        if best_candidate is None:
            raise RuntimeError(
                f"No candidate could be selected for read {read_info.read_id}."
            )

        selected[read_info.read_id] = best_candidate
        selected_ids = candidate_target_ids(best_candidate)
        fractional_assignment = 1.0 / len(selected_ids)

        for group_id in selected_ids:
            assigned_counts[group_id] += fractional_assignment

    return selected


def assignment_row(
    read_info: ReadCandidateSet,
    candidate: CandidateAlignment,
    assignment_type: str,
) -> dict:
    """Create one final assignment output row."""

    return {
        "read_id": candidate.read_id,
        "chrom": candidate.chrom,
        "genome_start": candidate.genome_start,
        "genome_end": candidate.genome_end,
        "strand": candidate.strand,
        "genome_cigar": candidate.genome_cigar,
        "selected_transcript_id": candidate.transcript_id,
        "alignment_transcript_id": candidate.alignment_transcript_id,
        "group_id": candidate.group_id,
        "posterior_probability": candidate.posterior_probability,
        "original_candidate_count": read_info.candidate_count,
        "assignment_type": assignment_type,
        "transcript_start": candidate.transcript_start,
        "transcript_end": candidate.transcript_end,
        "read_length": candidate.read_length,
        "hit_index": candidate.hit_index,
        "query_start": candidate.query_start,
        "query_end": candidate.query_end,
        "mapping_strand": candidate.mapping_strand,
        "transcript_cigar": candidate.transcript_cigar,
    }


def write_assignment_parquet(
    rows: list[dict],
    output_parquet: str | Path,
) -> None:
    """Write final one-alignment-per-read assignments."""

    output_path = Path(output_parquet)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    schema = pa.schema(
        [
            pa.field("read_id", pa.string(), nullable=False),
            pa.field("chrom", pa.string(), nullable=False),
            pa.field("genome_start", pa.int64(), nullable=False),
            pa.field("genome_end", pa.int64(), nullable=False),
            pa.field("strand", pa.string(), nullable=False),
            pa.field("genome_cigar", pa.string(), nullable=False),
            pa.field(
                "selected_transcript_id",
                pa.string(),
                nullable=False,
            ),
            pa.field(
                "alignment_transcript_id",
                pa.string(),
                nullable=False,
            ),
            pa.field("group_id", pa.string(), nullable=False),
            pa.field(
                "posterior_probability",
                pa.float64(),
                nullable=False,
            ),
            pa.field(
                "original_candidate_count",
                pa.int32(),
                nullable=False,
            ),
            pa.field("assignment_type", pa.string(), nullable=False),
            pa.field("transcript_start", pa.int64(), nullable=False),
            pa.field("transcript_end", pa.int64(), nullable=False),
            pa.field("read_length", pa.int32(), nullable=False),
            pa.field("hit_index", pa.int32(), nullable=False),
            pa.field("query_start", pa.int32(), nullable=False),
            pa.field("query_end", pa.int32(), nullable=False),
            pa.field("mapping_strand", pa.int8(), nullable=False),
            pa.field("transcript_cigar", pa.string(), nullable=False),
        ]
    )

    table = pa.Table.from_pylist(rows, schema=schema)
    pq.write_table(
        table,
        output_path,
        compression="zstd",
    )


def write_assignment_stats(
    rows: list[dict],
    output_tsv: str | Path,
) -> None:
    """Write assignment summary statistics."""

    counts: dict[str, int] = defaultdict(int)

    for row in rows:
        counts[row["assignment_type"]] += 1

    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(["metric", "value"])
        writer.writerow(["total_assigned_reads", len(rows)])

        for assignment_type in sorted(counts):
            writer.writerow(
                [
                    f"assignment_type_{assignment_type}",
                    counts[assignment_type],
                ]
            )


def create_final_genome_assignments(
    reads_info_parquet: str | Path,
    coordinate_table: str | Path,
    output_parquet: str | Path,
    output_stats_tsv: str | Path,
    batch_size: int = 250_000,
    local_bin_size: int = 30,
    logger: logging.Logger | None = None,
) -> Path:
    """Create one final genomic alignment assignment per mapped read."""

    logger = _create_logger(logger)
    logger.info("Final genome assignment workflow started")

    candidates_by_read = load_and_project_candidates(
        reads_info_parquet=reads_info_parquet,
        coordinate_table=coordinate_table,
        batch_size=batch_size,
        logger=logger,
    )

    read_sets = collapse_same_genome_alignments(candidates_by_read)

    unique_reads = [
        read_info
        for read_info in read_sets
        if read_info.candidate_count == 1
    ]
    multimap_reads = [
        read_info
        for read_info in read_sets
        if read_info.candidate_count > 1
    ]

    logger.info("Genome-unique reads: %d", len(unique_reads))
    logger.info("Genome multi-map reads: %d", len(multimap_reads))

    selected_multimap = assign_multimap_reads(
        multimap_reads=multimap_reads,
        local_bin_size=local_bin_size,
    )

    rows: list[dict] = []

    for read_info in unique_reads:
        rows.append(
            assignment_row(
                read_info=read_info,
                candidate=read_info.candidates[0],
                assignment_type="genome_unique",
            )
        )

    for read_info in multimap_reads:
        rows.append(
            assignment_row(
                read_info=read_info,
                candidate=selected_multimap[read_info.read_id],
                assignment_type="multimap_assigned",
            )
        )


    write_assignment_parquet(
        rows=rows,
        output_parquet=output_parquet,
    )
    write_assignment_stats(
        rows=rows,
        output_tsv=output_stats_tsv,
    )

    logger.info("Final assignments written: %d", len(rows))
    logger.info("Assignment Parquet: %s", output_parquet)
    logger.info("Assignment statistics: %s", output_stats_tsv)

    return Path(output_parquet)