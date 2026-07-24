from __future__ import annotations

import logging
import hashlib
import random
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
import pysam

from re_distribute.gene_groups import GeneGroupIndex


@dataclass(frozen=True)
class CandidatePatternResult:
    """Compact multimap candidate patterns used by the EM iterations."""

    pattern_counts: Counter[tuple[int, ...]]
    read_count: int
    candidate_row_count: int
    dropped_candidate_count: int
    dropped_read_count: int


@dataclass(frozen=True)
class EmResult:
    """Final gene-group abundance and EM convergence history."""

    final_tpm_table: pd.DataFrame
    history_table: pd.DataFrame
    pattern_result: CandidatePatternResult


@dataclass(frozen=True)
class AssignmentWriteResult:
    """Summary of sequence-group posterior allocation."""

    read_count: int
    selected_read_count: int
    unassigned_read_count: int
    candidate_row_count: int
    sequence_group_count: int
    resolved_sequence_group_count: int
    unresolved_sequence_group_count: int


def _require_pyarrow():
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as error:
        raise ImportError(
            "EM redistribution requires pyarrow. Install it with: "
            "pip install pyarrow"
        ) from error

    return pa, pq


def _parse_genomic_blocks(
    blocks_text: str | None,
) -> tuple[int, int] | None:
    """Return the outer genomic interval represented by a block string."""

    if not blocks_text:
        return None

    starts: list[int] = []
    ends: list[int] = []

    for block in str(blocks_text).split(","):
        if not block:
            continue
        start_text, end_text = block.split("-", maxsplit=1)
        starts.append(int(start_text))
        ends.append(int(end_text))

    if not starts:
        return None

    return min(starts), max(ends)


def _assign_candidate_group_id(
    row: dict,
    group_index: GeneGroupIndex,
) -> str | None:
    """Assign one genomic candidate to a strand-specific gene group."""

    interval = _parse_genomic_blocks(
        row.get("genomic_blocks")
    )

    if interval is not None:
        genome_start, genome_end = interval
    else:
        reference_start = row.get("reference_start")
        if reference_start is None:
            return None
        genome_start = int(reference_start)
        genome_end = genome_start + 1

    strand = "-" if bool(row.get("is_reverse")) else "+"
    group = group_index.find_containing_group(
        chrom=str(row["chrom"]),
        strand=strand,
        genome_start=genome_start,
        genome_end=genome_end,
    )

    if group is not None:
        return group.group_id

    gene_id = str(row.get("gene_id") or "")
    if gene_id:
        return group_index.group_id_by_gene.get(gene_id)

    return None


def _iter_read_rows(
    parquet_path: str,
    columns: list[str] | None = None,
    batch_size: int = 100_000,
) -> Iterator[tuple[str, list[dict]]]:
    """Stream candidate rows while keeping all rows of one read together."""

    _, pq = _require_pyarrow()
    parquet_file = pq.ParquetFile(parquet_path)

    current_read_id: str | None = None
    current_rows: list[dict] = []

    for batch in parquet_file.iter_batches(
        batch_size=batch_size,
        columns=columns,
    ):
        for row in batch.to_pylist():
            read_id = str(row["read_id"])

            if current_read_id is None:
                current_read_id = read_id

            if read_id != current_read_id:
                yield current_read_id, current_rows
                current_read_id = read_id
                current_rows = []

            current_rows.append(row)

    if current_read_id is not None:
        yield current_read_id, current_rows


def collect_candidate_patterns(
    mapping_parquet: str,
    group_index: GeneGroupIndex,
    group_index_by_id: dict[str, int],
    batch_size: int = 100_000,
    logger: logging.Logger | None = None,
) -> CandidatePatternResult:
    """Collapse multimap reads into counts of unique gene-group patterns."""

    mapping_path = Path(mapping_parquet)
    if not mapping_path.exists():
        raise FileNotFoundError(mapping_parquet)

    pattern_counts: Counter[tuple[int, ...]] = Counter()
    read_count = 0
    candidate_row_count = 0
    dropped_candidate_count = 0
    dropped_read_count = 0

    columns = [
        "read_id",
        "chrom",
        "reference_start",
        "genomic_blocks",
        "is_reverse",
        "gene_id",
    ]

    for read_id, rows in _iter_read_rows(
        parquet_path=str(mapping_path),
        columns=columns,
        batch_size=batch_size,
    ):
        del read_id
        read_count += 1
        candidate_row_count += len(rows)
        group_indices: set[int] = set()

        for row in rows:
            group_id = _assign_candidate_group_id(
                row=row,
                group_index=group_index,
            )
            if group_id is None:
                dropped_candidate_count += 1
                continue
            group_indices.add(group_index_by_id[group_id])

        if not group_indices:
            dropped_read_count += 1
            continue

        pattern_counts[tuple(sorted(group_indices))] += 1

    if logger is not None:
        logger.info(
            "EM patterns: reads=%d patterns=%d candidate_rows=%d "
            "dropped_candidates=%d dropped_reads=%d",
            read_count,
            len(pattern_counts),
            candidate_row_count,
            dropped_candidate_count,
            dropped_read_count,
        )

    return CandidatePatternResult(
        pattern_counts=pattern_counts,
        read_count=read_count,
        candidate_row_count=candidate_row_count,
        dropped_candidate_count=dropped_candidate_count,
        dropped_read_count=dropped_read_count,
    )


def run_em(
    initial_tpm_table: pd.DataFrame,
    mapping_parquet: str,
    group_index: GeneGroupIndex,
    pseudocount: float = 0.01,
    maximum_iterations: int = 100,
    convergence_threshold: float = 1e-7,
    parquet_batch_size: int = 100_000,
    logger: logging.Logger | None = None,
) -> EmResult:
    """
    Run gene-group EM using unique and multiple-hit reads.

    Alignment scores are never used. During every iteration, each multimap
    pattern is distributed in proportion to the current gene-group TPM. The
    resulting expected multimap counts are added to the fixed unique counts,
    corrected by the one-to-six-hit EM effective length, and converted back
    to TPM.
    """

    if maximum_iterations < 1:
        raise ValueError("maximum_iterations must be at least 1")
    if convergence_threshold <= 0:
        raise ValueError("convergence_threshold must be positive")
    if pseudocount < 0:
        raise ValueError("pseudocount must be non-negative")

    required_columns = {
        "gene_group_id",
        "unique_read_count",
        "em_effective_length",
        "initial_tpm",
    }
    missing_columns = required_columns - set(initial_tpm_table.columns)
    if missing_columns:
        raise ValueError(
            "Missing initial TPM columns: "
            + ", ".join(sorted(missing_columns))
        )

    table = initial_tpm_table.copy().reset_index(drop=True)
    group_ids = table["gene_group_id"].astype(str).tolist()
    group_index_by_id = {
        group_id: index
        for index, group_id in enumerate(group_ids)
    }

    pattern_result = collect_candidate_patterns(
        mapping_parquet=mapping_parquet,
        group_index=group_index,
        group_index_by_id=group_index_by_id,
        batch_size=parquet_batch_size,
        logger=logger,
    )

    unique_counts = table["unique_read_count"].to_numpy(dtype=float)
    effective_lengths = table["em_effective_length"].to_numpy(dtype=float)
    effective_lengths = np.maximum(effective_lengths, 1.0)
    current_tpm = table["initial_tpm"].to_numpy(dtype=float)

    history_rows: list[dict] = []
    final_expected_multi = np.zeros_like(current_tpm)

    for iteration in range(1, maximum_iterations + 1):
        expected_multi = np.zeros_like(current_tpm)

        for pattern, pattern_count in pattern_result.pattern_counts.items():
            pattern_array = np.fromiter(
                pattern,
                dtype=np.int64,
            )
            weights = current_tpm[pattern_array]
            weight_sum = float(weights.sum())

            if weight_sum <= 0:
                probabilities = np.full(
                    len(pattern_array),
                    1.0 / len(pattern_array),
                    dtype=float,
                )
            else:
                probabilities = weights / weight_sum

            expected_multi[pattern_array] += (
                probabilities
                * float(pattern_count)
            )

        total_counts = (
            unique_counts
            + expected_multi
            + pseudocount
        )
        rpk = total_counts / (effective_lengths / 1000.0)
        rpk_sum = float(rpk.sum())

        if rpk_sum <= 0:
            raise ValueError("EM RPK sum is zero")

        next_tpm = rpk / rpk_sum * 1_000_000.0
        relative_change = float(
            np.abs(next_tpm - current_tpm).sum()
            / 1_000_000.0
        )

        history_rows.append(
            {
                "iteration": iteration,
                "relative_tpm_change": relative_change,
                "total_unique_count": float(unique_counts.sum()),
                "total_multimap_expected_count": float(
                    expected_multi.sum()
                ),
                "total_count": float(total_counts.sum()),
                "total_tpm": float(next_tpm.sum()),
            }
        )

        if logger is not None:
            logger.info(
                "EM iteration=%d relative_change=%.10g "
                "expected_multi=%.3f",
                iteration,
                relative_change,
                float(expected_multi.sum()),
            )

        current_tpm = next_tpm
        final_expected_multi = expected_multi

        if relative_change < convergence_threshold:
            break

    final_total_count = (
        unique_counts
        + final_expected_multi
        + pseudocount
    )
    final_rpk = final_total_count / (effective_lengths / 1000.0)

    final_table = table.copy()
    final_table["expected_multimap_count"] = final_expected_multi
    final_table["final_total_count"] = final_total_count
    final_table["final_rpk"] = final_rpk
    final_table["final_tpm"] = current_tpm
    final_table["em_iterations"] = len(history_rows)
    final_table["em_converged"] = (
        history_rows[-1]["relative_tpm_change"]
        < convergence_threshold
    )

    return EmResult(
        final_tpm_table=final_table,
        history_table=pd.DataFrame(history_rows),
        pattern_result=pattern_result,
    )


def _candidate_signature(row: dict) -> tuple:
    """Return a stable genomic candidate key independent of row order."""

    return (
        str(row.get("chrom") or ""),
        int(row.get("reference_start") or 0),
        bool(row.get("is_reverse")),
        str(row.get("cigar") or ""),
        str(row.get("genomic_blocks") or ""),
        str(row.get("gene_id") or ""),
    )


def _stable_sequence_seed(
    sequence: str,
    random_seed: int,
) -> int:
    """Return a reproducible integer seed for one sequence group."""

    digest = hashlib.blake2b(
        f"{random_seed}:{sequence}".encode("utf-8"),
        digest_size=8,
    ).digest()
    return int.from_bytes(digest, byteorder="big", signed=False)


def _allocate_integer_counts(
    read_count: int,
    probabilities: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Round expected candidate counts while preserving the read total."""

    probabilities = np.asarray(probabilities, dtype=float)
    probability_sum = float(probabilities.sum())

    if read_count < 1 or probability_sum <= 0:
        return np.zeros(len(probabilities), dtype=int), 1.0

    probabilities = probabilities / probability_sum
    expected_counts = probabilities * float(read_count)
    assigned_counts = np.floor(expected_counts).astype(int)
    remaining_count = read_count - int(assigned_counts.sum())

    if remaining_count > 0:
        fractional_parts = expected_counts - assigned_counts
        remainder_order = np.argsort(
            -fractional_parts,
            kind="stable",
        )
        assigned_counts[remainder_order[:remaining_count]] += 1

    allocation_error = float(
        np.abs(assigned_counts - expected_counts).sum()
        / (2.0 * float(read_count))
    )

    return assigned_counts, allocation_error


def _calculate_candidate_posteriors(
    rows: list[dict],
    final_tpm_by_group: dict[str, float],
    group_index: GeneGroupIndex,
) -> tuple[list[str | None], list[float], list[float]]:
    """Calculate group and genomic-candidate posterior probabilities."""

    group_id_by_row: list[str | None] = []
    row_indices_by_group: dict[str, list[int]] = {}

    for row_index, row in enumerate(rows):
        group_id = _assign_candidate_group_id(
            row=row,
            group_index=group_index,
        )
        group_id_by_row.append(group_id)
        if group_id is not None:
            row_indices_by_group.setdefault(group_id, []).append(row_index)

    group_weights = {
        group_id: max(float(final_tpm_by_group.get(group_id, 0.0)), 0.0)
        for group_id in row_indices_by_group
    }
    total_group_weight = float(sum(group_weights.values()))

    if not row_indices_by_group:
        return (
            group_id_by_row,
            [0.0] * len(rows),
            [0.0] * len(rows),
        )

    if total_group_weight <= 0:
        group_posterior_by_id = {
            group_id: 1.0 / len(row_indices_by_group)
            for group_id in row_indices_by_group
        }
    else:
        group_posterior_by_id = {
            group_id: weight / total_group_weight
            for group_id, weight in group_weights.items()
        }

    group_posteriors: list[float] = []
    candidate_posteriors: list[float] = []

    for group_id in group_id_by_row:
        if group_id is None:
            group_posteriors.append(0.0)
            candidate_posteriors.append(0.0)
            continue

        group_posterior = group_posterior_by_id[group_id]
        group_candidate_count = len(row_indices_by_group[group_id])
        group_posteriors.append(float(group_posterior))
        candidate_posteriors.append(
            float(group_posterior / group_candidate_count)
        )

    return group_id_by_row, group_posteriors, candidate_posteriors


def _load_read_sequences(
    multiplemap_bam: str,
) -> dict[str, str]:
    """Load the stored sequence for every unresolved multimap read."""

    read_sequence: dict[str, str] = {}

    with pysam.AlignmentFile(
        multiplemap_bam,
        "rb",
        check_sq=False,
    ) as bam:
        for record in bam.fetch(until_eof=True):
            read_id = str(record.query_name)
            sequence = str(record.query_sequence or "")
            if read_id in read_sequence:
                raise ValueError(
                    "Duplicate read ID in multiple-map BAM: "
                    f"{read_id}"
                )
            read_sequence[read_id] = sequence

    return read_sequence


def write_final_assignments(
    mapping_parquet: str,
    multiplemap_bam: str,
    final_tpm_table: pd.DataFrame,
    group_index: GeneGroupIndex,
    output_parquet: str,
    maximum_allocation_error: float = 0.10,
    random_seed: int = 13,
    batch_size: int = 100_000,
    logger: logging.Logger | None = None,
) -> AssignmentWriteResult:
    """
    Allocate identical read sequences according to final EM posteriors.

    Reads are grouped by exact sequence and genomic candidate set. Candidate
    expected counts are rounded with the largest-remainder method. A sequence
    group is written to genomic positions only when the normalized rounding
    error does not exceed ``maximum_allocation_error``. Otherwise, every read
    in that sequence group remains unresolved in the output BAM.
    """

    if not 0.0 <= maximum_allocation_error <= 1.0:
        raise ValueError(
            "maximum_allocation_error must be between 0 and 1"
        )

    pa, pq = _require_pyarrow()
    output_path = Path(output_parquet)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    final_tpm_by_group = dict(
        zip(
            final_tpm_table["gene_group_id"].astype(str),
            final_tpm_table["final_tpm"].astype(float),
        )
    )

    read_sequence = _load_read_sequences(multiplemap_bam)
    rows_by_read: dict[str, list[dict]] = {}
    candidate_row_count = 0

    for read_id, rows in _iter_read_rows(
        parquet_path=mapping_parquet,
        columns=None,
        batch_size=batch_size,
    ):
        rows_by_read[read_id] = sorted(rows, key=_candidate_signature)
        candidate_row_count += len(rows)

    missing_sequence_ids = set(rows_by_read) - set(read_sequence)
    if missing_sequence_ids:
        example = sorted(missing_sequence_ids)[:5]
        raise ValueError(
            "Reads in mapping Parquet are absent from multiple-map BAM: "
            + ", ".join(example)
        )

    # Exact sequence alone should normally determine the same candidate set.
    # Candidate signature is included in the key as a safety guard.
    grouped_read_ids: dict[tuple[str, tuple[tuple, ...]], list[str]] = {}
    sequence_candidate_sets: dict[str, set[tuple[tuple, ...]]] = {}

    for read_id, rows in rows_by_read.items():
        sequence = read_sequence[read_id]
        signature = tuple(_candidate_signature(row) for row in rows)
        grouped_read_ids.setdefault((sequence, signature), []).append(read_id)
        sequence_candidate_sets.setdefault(sequence, set()).add(signature)

    inconsistent_sequence_count = sum(
        len(signatures) > 1
        for signatures in sequence_candidate_sets.values()
    )
    if logger is not None and inconsistent_sequence_count:
        logger.warning(
            "Identical sequences with different candidate sets: %d; "
            "these were split into separate allocation groups",
            inconsistent_sequence_count,
        )

    writer = None
    output_rows: list[dict] = []
    selected_read_count = 0
    unassigned_read_count = 0
    resolved_sequence_group_count = 0
    unresolved_sequence_group_count = 0

    def flush() -> None:
        nonlocal writer
        if not output_rows:
            return
        table = pa.Table.from_pylist(output_rows)
        if writer is None:
            writer = pq.ParquetWriter(
                str(output_path),
                table.schema,
                compression="zstd",
            )
        writer.write_table(table)
        output_rows.clear()

    # Process sequence groups in insertion order. The assignment itself
    # remains reproducible because each sequence group uses a stable
    # sequence-derived random seed. Avoiding a global sort substantially
    # reduces runtime and memory use for large read sets.
    for sequence_group_index, ((sequence, signature), read_ids) in enumerate(
        grouped_read_ids.items(),
        start=1,
    ):
        del signature
        representative_rows = rows_by_read[read_ids[0]]
        (
            group_id_by_row,
            group_posteriors,
            candidate_posteriors,
        ) = _calculate_candidate_posteriors(
            rows=representative_rows,
            final_tpm_by_group=final_tpm_by_group,
            group_index=group_index,
        )

        assigned_counts, allocation_error = _allocate_integer_counts(
            read_count=len(read_ids),
            probabilities=np.asarray(candidate_posteriors, dtype=float),
        )
        probability_sum = float(sum(candidate_posteriors))
        resolved = (
            probability_sum > 0
            and allocation_error <= maximum_allocation_error
        )

        selected_candidate_by_read: dict[str, int] = {}

        if resolved:
            resolved_sequence_group_count += 1
            # Preserve the deterministic input order and shuffle a copy.
            # No string sorting is required.
            shuffled_read_ids = list(read_ids)
            rng = random.Random(
                _stable_sequence_seed(sequence, random_seed)
            )
            rng.shuffle(shuffled_read_ids)

            offset = 0
            for candidate_index, assigned_count in enumerate(assigned_counts):
                for read_id in shuffled_read_ids[
                    offset:offset + int(assigned_count)
                ]:
                    selected_candidate_by_read[read_id] = candidate_index
                offset += int(assigned_count)

            if offset != len(shuffled_read_ids):
                raise RuntimeError(
                    "Integer candidate allocation did not consume all reads"
                )
            selected_read_count += len(read_ids)
        else:
            unresolved_sequence_group_count += 1
            unassigned_read_count += len(read_ids)

        for read_id in read_ids:
            rows = rows_by_read[read_id]
            selected_candidate_index = selected_candidate_by_read.get(read_id)

            for row_index, row in enumerate(rows):
                enriched_row = dict(row)
                group_id = group_id_by_row[row_index]
                enriched_row["gene_group_id"] = group_id or ""
                enriched_row["group_posterior"] = float(
                    group_posteriors[row_index]
                )
                enriched_row["candidate_posterior"] = float(
                    candidate_posteriors[row_index]
                )
                enriched_row["sequence_group_id"] = int(
                    sequence_group_index
                )
                enriched_row["sequence_group_size"] = int(len(read_ids))
                enriched_row["allocation_error"] = float(allocation_error)
                enriched_row["maximum_allocation_error"] = float(
                    maximum_allocation_error
                )
                enriched_row["resolved"] = bool(resolved)
                enriched_row["selected"] = bool(
                    resolved and row_index == selected_candidate_index
                )
                enriched_row["assignment_method"] = (
                    "SEQUENCE_GROUP_EM"
                    if resolved
                    else "UNRESOLVED_ALLOCATION_ERROR"
                )
                enriched_row["random_seed"] = int(random_seed)
                output_rows.append(enriched_row)

            if len(output_rows) >= batch_size:
                flush()

    flush()
    if writer is not None:
        writer.close()

    read_count = len(rows_by_read)
    sequence_group_count = len(grouped_read_ids)

    if logger is not None:
        logger.info(
            "Final sequence allocation: reads=%d selected=%d unresolved=%d "
            "sequence_groups=%d resolved_groups=%d unresolved_groups=%d "
            "candidate_rows=%d maximum_error=%.4f output=%s",
            read_count,
            selected_read_count,
            unassigned_read_count,
            sequence_group_count,
            resolved_sequence_group_count,
            unresolved_sequence_group_count,
            candidate_row_count,
            maximum_allocation_error,
            output_path,
        )

    return AssignmentWriteResult(
        read_count=read_count,
        selected_read_count=selected_read_count,
        unassigned_read_count=unassigned_read_count,
        candidate_row_count=candidate_row_count,
        sequence_group_count=sequence_group_count,
        resolved_sequence_group_count=resolved_sequence_group_count,
        unresolved_sequence_group_count=unresolved_sequence_group_count,
    )