from __future__ import annotations

import logging
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pysam

from re_distribute.gene_groups import GeneGroup, GeneGroupIndex


@dataclass(frozen=True)
class UniqueReadCountResult:
    """Summary of assigning unique BAM records to gene groups."""

    counts: Counter[str]
    mapped_primary_records: int
    assigned_records: int
    unassigned_records: int
    skipped_records: int


@dataclass(frozen=True)
class InitialTpmResult:
    """Initial abundance table and unique-read assignment statistics."""

    table: pd.DataFrame
    unique_read_counts: UniqueReadCountResult


def count_unique_reads_by_group(
    unique_bam: str,
    group_index: GeneGroupIndex,
    progress_interval: int = 1_000_000,
    logger: logging.Logger | None = None,
) -> UniqueReadCountResult:
    """
    Count primary unique BAM records fully contained in a gene group.

    The BAM strand must match the strand of the gene group. Secondary,
    supplementary, unmapped and duplicate records are ignored.
    """

    bam_path = Path(unique_bam)
    if not bam_path.exists():
        raise FileNotFoundError(unique_bam)

    if logger is None:
        logger = logging.getLogger("RiboTxAlignRedistribute")

    counts: Counter[str] = Counter()
    mapped_primary_records = 0
    assigned_records = 0
    unassigned_records = 0
    skipped_records = 0
    start_time = time.time()

    with pysam.AlignmentFile(
        str(bam_path),
        "rb",
    ) as bam:
        for record_index, record in enumerate(
            bam.fetch(until_eof=True),
            start=1,
        ):
            if (
                record.is_unmapped
                or record.is_secondary
                or record.is_supplementary
                or record.is_duplicate
            ):
                skipped_records += 1
                continue

            mapped_primary_records += 1

            if record.reference_end is None:
                unassigned_records += 1
                continue

            chrom = bam.get_reference_name(record.reference_id)
            strand = "-" if record.is_reverse else "+"
            group = group_index.find_containing_group(
                chrom=chrom,
                strand=strand,
                genome_start=int(record.reference_start),
                genome_end=int(record.reference_end),
            )

            if group is None:
                unassigned_records += 1
            else:
                counts[group.group_id] += 1
                assigned_records += 1

            if (
                progress_interval > 0
                and record_index % progress_interval == 0
            ):
                elapsed_seconds = time.time() - start_time
                speed = (
                    record_index / elapsed_seconds
                    if elapsed_seconds > 0
                    else 0.0
                )
                logger.info(
                    "Unique BAM progress: records=%d speed=%.1f records/s "
                    "assigned=%d unassigned=%d",
                    record_index,
                    speed,
                    assigned_records,
                    unassigned_records,
                )

    return UniqueReadCountResult(
        counts=counts,
        mapped_primary_records=mapped_primary_records,
        assigned_records=assigned_records,
        unassigned_records=unassigned_records,
        skipped_records=skipped_records,
    )


def _interval_overlap_length(
    start1: int,
    end1: int,
    start2: int,
    end2: int,
) -> int:
    """Return the overlap length between two half-open intervals."""

    return max(0, min(end1, end2) - max(start1, start2))


def _calculate_group_effective_lengths(
    group: GeneGroup,
    group_index: GeneGroupIndex,
    gene_unique_mappable_length: dict[str, float],
    gene_em_mappable_length: dict[str, float],
) -> tuple[float, float, float, int]:
    """
    Estimate non-duplicated initial and EM effective lengths for one group.

    Existing gene-level mappability values are converted into per-gene
    mappable fractions using each gene's merged coordinate length. Group
    intervals are split at every gene boundary. Each atomic segment is counted
    once and weighted by the maximum fraction among genes covering it.

    The initial length uses unique-mappable sequence only. The EM length uses
    sequence with one to six genomic hits, because those reads are included in
    the redistribution procedure. The original mappability calculation is not
    rerun or modified.
    """

    unique_fractions: dict[str, float] = {}
    em_fractions: dict[str, float] = {}
    boundaries: set[int] = set()

    for gene_id in group.gene_ids:
        intervals = group_index.gene_intervals.get(gene_id, ())
        gene_length = sum(end - start for start, end in intervals)
        unique_length = max(
            float(gene_unique_mappable_length.get(gene_id, 0.0)),
            0.0,
        )
        em_length = max(
            float(gene_em_mappable_length.get(gene_id, 0.0)),
            0.0,
        )

        if gene_length > 0:
            unique_fractions[gene_id] = min(
                unique_length / gene_length,
                1.0,
            )
            em_fractions[gene_id] = min(
                em_length / gene_length,
                1.0,
            )
        else:
            unique_fractions[gene_id] = 0.0
            em_fractions[gene_id] = 0.0

        for interval_start, interval_end in intervals:
            boundaries.add(interval_start)
            boundaries.add(interval_end)

    sorted_boundaries = sorted(boundaries)
    initial_effective_length = 0.0
    em_effective_length = 0.0
    covered_union_length = 0.0
    atomic_segment_count = 0

    for left, right in zip(
        sorted_boundaries,
        sorted_boundaries[1:],
    ):
        if right <= left:
            continue

        covering_genes: list[str] = []

        for gene_id in group.gene_ids:
            intervals = group_index.gene_intervals.get(gene_id, ())
            if any(
                _interval_overlap_length(
                    left,
                    right,
                    interval_start,
                    interval_end,
                ) == right - left
                for interval_start, interval_end in intervals
            ):
                covering_genes.append(gene_id)

        if not covering_genes:
            continue

        segment_length = right - left
        covered_union_length += segment_length
        initial_effective_length += segment_length * max(
            unique_fractions.get(gene_id, 0.0)
            for gene_id in covering_genes
        )
        em_effective_length += segment_length * max(
            em_fractions.get(gene_id, 0.0)
            for gene_id in covering_genes
        )
        atomic_segment_count += 1

    return (
        initial_effective_length,
        em_effective_length,
        covered_union_length,
        atomic_segment_count,
    )

def calculate_initial_tpm(
    unique_bam: str,
    mappability_parquet: str,
    group_index: GeneGroupIndex,
    pseudocount: float = 0.01,
    progress_interval: int = 1_000_000,
    logger: logging.Logger | None = None,
) -> InitialTpmResult:
    """
    Calculate initial gene-group TPM with non-duplicated group lengths.

    Mappability is read from the existing gene-level Parquet. The expensive
    mappability calculation is not repeated. Length duplication is removed
    only here, during TPM effective-length construction.
    """

    if pseudocount < 0:
        raise ValueError("pseudocount must be non-negative")

    mappability_path = Path(mappability_parquet)
    if not mappability_path.exists():
        raise FileNotFoundError(mappability_parquet)

    unique_result = count_unique_reads_by_group(
        unique_bam=unique_bam,
        group_index=group_index,
        progress_interval=progress_interval,
        logger=logger,
    )

    mappability_df = pd.read_parquet(
        mappability_path,
        columns=[
            "gene_id",
            "unique_mappable_length",
            "hit2_mappable_length",
            "hit3_mappable_length",
            "hit4_mappable_length",
            "hit5_mappable_length",
            "hit6_mappable_length",
            "hit7_mappable_length",
            "hit8_mappable_length",
        ],
    )

    required_columns = {
        "gene_id",
        "unique_mappable_length",
        "hit2_mappable_length",
        "hit3_mappable_length",
        "hit4_mappable_length",
        "hit5_mappable_length",
        "hit6_mappable_length",
        "hit7_mappable_length",
        "hit8_mappable_length",
    }
    missing_columns = required_columns - set(mappability_df.columns)
    if missing_columns:
        raise ValueError(
            "Missing mappability columns: "
            + ", ".join(sorted(missing_columns))
        )

    numeric_columns = [
        "unique_mappable_length",
        "hit2_mappable_length",
        "hit3_mappable_length",
        "hit4_mappable_length",
        "hit5_mappable_length",
        "hit6_mappable_length",
        "hit7_mappable_length",
        "hit8_mappable_length",
    ]
    grouped_mappability = (
        mappability_df.groupby(
            "gene_id",
            sort=False,
        )[numeric_columns]
        .max()
        .astype(float)
    )
    gene_unique_mappable_length = (
        grouped_mappability["unique_mappable_length"].to_dict()
    )
    gene_em_mappable_length = (
        grouped_mappability[numeric_columns].sum(axis=1).to_dict()
    )

    rows: list[dict] = []

    for group in group_index.groups:
        unique_read_count = int(
            unique_result.counts.get(group.group_id, 0)
        )
        raw_gene_unique_mappable_sum = sum(
            max(
                float(gene_unique_mappable_length.get(gene_id, 0.0)),
                0.0,
            )
            for gene_id in group.gene_ids
        )
        raw_gene_em_mappable_sum = sum(
            max(
                float(gene_em_mappable_length.get(gene_id, 0.0)),
                0.0,
            )
            for gene_id in group.gene_ids
        )
        (
            nonduplicated_unique_length,
            nonduplicated_em_length,
            covered_union_length,
            atomic_segment_count,
        ) = _calculate_group_effective_lengths(
            group=group,
            group_index=group_index,
            gene_unique_mappable_length=gene_unique_mappable_length,
            gene_em_mappable_length=gene_em_mappable_length,
        )

        initial_effective_length = max(
            nonduplicated_unique_length,
            1.0,
        )
        em_effective_length = max(
            nonduplicated_em_length,
            1.0,
        )
        initial_count = unique_read_count + pseudocount
        initial_rpk = initial_count / (
            initial_effective_length / 1000.0
        )

        rows.append(
            {
                "gene_group_id": group.group_id,
                "chrom": group.chrom,
                "strand": group.strand,
                "genome_start": group.genome_start,
                "genome_end": group.genome_end,
                "gene_ids": ",".join(group.gene_ids),
                "unique_read_count": unique_read_count,
                "raw_gene_unique_mappable_sum": (
                    raw_gene_unique_mappable_sum
                ),
                "raw_gene_em_mappable_sum": raw_gene_em_mappable_sum,
                "union_gene_length": float(group.union_length),
                "covered_union_length": covered_union_length,
                "nonduplicated_unique_mappable_length": (
                    nonduplicated_unique_length
                ),
                "nonduplicated_em_mappable_length": (
                    nonduplicated_em_length
                ),
                "initial_effective_length": initial_effective_length,
                "em_effective_length": em_effective_length,
                "initial_length_duplicate_removed": (
                    raw_gene_unique_mappable_sum
                    - nonduplicated_unique_length
                ),
                "em_length_duplicate_removed": (
                    raw_gene_em_mappable_sum
                    - nonduplicated_em_length
                ),
                "atomic_segment_count": atomic_segment_count,
                "initial_mappability_missing": (
                    nonduplicated_unique_length <= 0
                ),
                "em_mappability_missing": (
                    nonduplicated_em_length <= 0
                ),
                "pseudocount": pseudocount,
                "initial_count": initial_count,
                "initial_rpk": initial_rpk,
            }
        )

    initial_df = pd.DataFrame(rows)
    total_rpk = float(initial_df["initial_rpk"].sum())

    if total_rpk <= 0:
        raise ValueError("Initial RPK sum is zero")

    initial_df["initial_tpm"] = (
        initial_df["initial_rpk"]
        / total_rpk
        * 1_000_000.0
    )

    return InitialTpmResult(
        table=initial_df,
        unique_read_counts=unique_result,
    )
