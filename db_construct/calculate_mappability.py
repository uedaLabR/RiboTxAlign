from __future__ import annotations

import csv
import logging
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pysam


@dataclass(frozen=True)
class TranscriptRegion:
    """Genomic span assigned to one synthetic transcript reference."""

    reference_id: str
    gene_id: str
    chrom: str
    strand: str
    genome_start: int
    genome_end: int
    group_id: int


def _encode_kmer_2bit(sequence: str) -> int | None:
    """Encode an A/C/G/T k-mer into one integer."""

    code = 0
    for base in sequence:
        code <<= 2
        if base == "A":
            continue
        if base == "C":
            code |= 1
        elif base == "G":
            code |= 2
        elif base == "T":
            code |= 3
        else:
            return None
    return code


def _load_reference_spans(
    coordinate_tsv: str,
) -> dict[str, tuple[str, str, str, int, int]]:
    """Load one genomic span for each synthetic reference."""

    reference_rows: dict[str, list[tuple[str, str, str, int, int]]] = {}

    with open(coordinate_tsv, encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        required_columns = {
            "reference_id",
            "gene_id",
            "chrom",
            "genome_start",
            "genome_end",
            "strand",
        }
        missing_columns = required_columns - set(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(
                "Missing coordinate-table columns: "
                + ", ".join(sorted(missing_columns))
            )

        for row in reader:
            reference_rows.setdefault(row["reference_id"], []).append(
                (
                    row["gene_id"],
                    row["chrom"],
                    row["strand"],
                    int(row["genome_start"]),
                    int(row["genome_end"]),
                )
            )

    reference_spans: dict[str, tuple[str, str, str, int, int]] = {}

    for reference_id, rows in reference_rows.items():
        gene_ids = {row[0] for row in rows}
        chroms = {row[1] for row in rows}
        strands = {row[2] for row in rows}

        if len(gene_ids) != 1:
            raise ValueError(f"Reference maps to multiple genes: {reference_id}")
        if len(chroms) != 1:
            raise ValueError(f"Reference maps to multiple chromosomes: {reference_id}")
        if len(strands) != 1:
            raise ValueError(f"Reference has mixed strands: {reference_id}")

        reference_spans[reference_id] = (
            next(iter(gene_ids)),
            next(iter(chroms)),
            next(iter(strands)),
            min(row[3] for row in rows),
            max(row[4] for row in rows),
        )

    return reference_spans


def _assign_overlapping_reference_groups(
    reference_spans: dict[str, tuple[str, str, str, int, int]],
) -> dict[str, TranscriptRegion]:
    """
    Assign one group ID to overlapping transcript spans.

    Grouping is chromosome- and strand-specific, and overlap is transitive.
    """

    by_chrom_strand: dict[
        tuple[str, str],
        list[tuple[str, str, int, int]],
    ] = {}

    for reference_id, (
        gene_id,
        chrom,
        strand,
        genome_start,
        genome_end,
    ) in reference_spans.items():
        by_chrom_strand.setdefault((chrom, strand), []).append(
            (reference_id, gene_id, genome_start, genome_end)
        )

    transcript_regions: dict[str, TranscriptRegion] = {}
    next_group_id = 0

    for (chrom, strand), rows in by_chrom_strand.items():
        rows.sort(key=lambda item: (item[2], item[3], item[0]))

        current_group_id: int | None = None
        current_group_end: int | None = None

        for reference_id, gene_id, genome_start, genome_end in rows:
            if (
                current_group_id is None
                or current_group_end is None
                or genome_start >= current_group_end
            ):
                current_group_id = next_group_id
                next_group_id += 1
                current_group_end = genome_end
            else:
                current_group_end = max(current_group_end, genome_end)

            transcript_regions[reference_id] = TranscriptRegion(
                reference_id=reference_id,
                gene_id=gene_id,
                chrom=chrom,
                strand=strand,
                genome_start=genome_start,
                genome_end=genome_end,
                group_id=current_group_id,
            )

    return transcript_regions


def calculate_genomic_mappability(
    reference_fasta: str,
    coordinate_tsv: str,
    output_prefix: str,
    read_length: int = 30,
    step: int = 1,
    progress_reference_interval: int = 100,
) -> pd.DataFrame:
    """
    Calculate exact-match mappability from encoded k-mers.

    Each synthetic reference is assigned to an overlapping genomic transcript
    group. For each encoded k-mer, only distinct group IDs are retained.
    Therefore, the same sequence appearing in multiple overlapping transcript
    references is counted as one genomic hit.
    """

    if read_length <= 0:
        raise ValueError("read_length must be positive")
    if read_length > 32:
        raise ValueError("2-bit encoding supports at most 32 bases")
    if step <= 0:
        raise ValueError("step must be positive")
    if progress_reference_interval < 0:
        raise ValueError("progress_reference_interval must be non-negative")

    reference_path = Path(reference_fasta)
    coordinate_path = Path(coordinate_tsv)

    if not reference_path.exists():
        raise FileNotFoundError(reference_fasta)
    if not coordinate_path.exists():
        raise FileNotFoundError(coordinate_tsv)

    Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("RiboTxAlignMappability")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter(
                "%(asctime)s\t%(levelname)s\t%(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    reference_spans = _load_reference_spans(str(coordinate_path))
    transcript_regions = _assign_overlapping_reference_groups(reference_spans)

    # k-mer code -> distinct overlapping-transcript group IDs.
    kmer_to_group_ids: dict[int, set[int]] = {}

    # Source k-mers are retained only until the gene-level summary is built.
    reference_kmers: dict[str, list[int]] = {}
    ambiguous_count_by_gene: Counter[str] = Counter()

    start_time = time.time()

    with pysam.FastaFile(str(reference_path)) as fasta:
        valid_reference_ids = [
            reference_id
            for reference_id in fasta.references
            if reference_id in transcript_regions
        ]
        total_reference_count = len(valid_reference_ids)

        logger.info(
            "Mappability collection started: references=%d read_length=%d step=%d",
            total_reference_count,
            read_length,
            step,
        )

        for processed_reference_count, reference_id in enumerate(
            valid_reference_ids,
            start=1,
        ):
            region = transcript_regions[reference_id]
            sequence = fasta.fetch(reference_id).upper()
            possible_start_count = max(len(sequence) - read_length + 1, 0)
            source_kmers: list[int] = []

            for start in range(0, possible_start_count, step):
                kmer = sequence[start:start + read_length]
                kmer_code = _encode_kmer_2bit(kmer)

                if kmer_code is None:
                    ambiguous_count_by_gene[region.gene_id] += 1
                    continue

                source_kmers.append(kmer_code)
                group_ids = kmer_to_group_ids.get(kmer_code)

                if group_ids is None:
                    kmer_to_group_ids[kmer_code] = {region.group_id}
                else:
                    group_ids.add(region.group_id)

            reference_kmers[reference_id] = source_kmers

            if (
                progress_reference_interval > 0
                and (
                    processed_reference_count % progress_reference_interval == 0
                    or processed_reference_count == total_reference_count
                )
            ):
                elapsed_seconds = time.time() - start_time
                references_per_second = (
                    processed_reference_count / elapsed_seconds
                    if elapsed_seconds > 0
                    else 0.0
                )
                logger.info(
                    "Mappability progress: references=%d/%d (%.1f%%) "
                    "speed=%.2f references/s distinct_kmers=%d",
                    processed_reference_count,
                    total_reference_count,
                    100.0 * processed_reference_count / total_reference_count
                    if total_reference_count > 0
                    else 100.0,
                    references_per_second,
                    len(kmer_to_group_ids),
                )

    # Count each source k-mer only once per gene.
    gene_to_kmers: dict[str, set[int]] = {}
    for reference_id, source_kmers in reference_kmers.items():
        gene_id = transcript_regions[reference_id].gene_id
        gene_to_kmers.setdefault(gene_id, set()).update(source_kmers)

    gene_metadata: dict[str, tuple[str, str]] = {}
    for region in transcript_regions.values():
        gene_metadata.setdefault(region.gene_id, (region.chrom, region.strand))

    summary_rows: list[dict] = []

    for gene_id in sorted(gene_to_kmers):
        hit_counter: Counter[int] = Counter()

        for kmer_code in gene_to_kmers[gene_id]:
            hit_counter[len(kmer_to_group_ids[kmer_code])] += 1

        evaluated_count = sum(hit_counter.values())
        chrom, strand = gene_metadata.get(gene_id, ("", ""))

        summary_rows.append(
            {
                "gene_id": gene_id,
                "chrom": chrom,
                "strand": strand,
                "read_length": read_length,
                "step": step,
                "genomic_start_count": evaluated_count,
                "ambiguous_sequence_count": ambiguous_count_by_gene.get(gene_id, 0),
                "unmapped_mappable_length": hit_counter.get(0, 0),
                "unique_mappable_length": hit_counter.get(1, 0),
                "hit2_mappable_length": hit_counter.get(2, 0),
                "hit3_mappable_length": hit_counter.get(3, 0),
                "hit4_mappable_length": hit_counter.get(4, 0),
                "hit5_mappable_length": hit_counter.get(5, 0),
                "hit6_mappable_length": hit_counter.get(6, 0),
                "hit7_mappable_length": hit_counter.get(7, 0),
                "hit8_mappable_length": hit_counter.get(8, 0),
                "hit_gt8_mappable_length": sum(
                    count
                    for hit_count, count in hit_counter.items()
                    if hit_count >= 9
                ),
                "unique_fraction": (
                    hit_counter.get(1, 0) / evaluated_count
                    if evaluated_count > 0
                    else 0.0
                ),
                "hit_unit": "overlapping_transcript_group",
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_file = f"{output_prefix}.mappability.parquet"
    summary_df.to_parquet(summary_file, index=False)

    total_evaluated_length = int(
        summary_df["genomic_start_count"].sum()
    ) if not summary_df.empty else 0
    total_unique_length = int(
        summary_df["unique_mappable_length"].sum()
    ) if not summary_df.empty else 0

    hit_up_to_8_columns = [
        "unique_mappable_length",
        "hit2_mappable_length",
        "hit3_mappable_length",
        "hit4_mappable_length",
        "hit5_mappable_length",
        "hit6_mappable_length",
        "hit7_mappable_length",
        "hit8_mappable_length",
    ]
    total_hit_up_to_8_length = int(
        summary_df[hit_up_to_8_columns].sum().sum()
    ) if not summary_df.empty else 0
    total_hit_over_8_length = int(
        summary_df["hit_gt8_mappable_length"].sum()
    ) if not summary_df.empty else 0

    unique_percent = (
        100.0 * total_unique_length / total_evaluated_length
        if total_evaluated_length > 0
        else 0.0
    )
    hit_up_to_8_percent = (
        100.0 * total_hit_up_to_8_length / total_evaluated_length
        if total_evaluated_length > 0
        else 0.0
    )
    hit_over_8_percent = (
        100.0 * total_hit_over_8_length / total_evaluated_length
        if total_evaluated_length > 0
        else 0.0
    )

    logger.info(
        "Mappability summary: evaluated_length=%d unique_length=%d "
        "unique_percent=%.2f hit1_to_8_length=%d hit1_to_8_percent=%.2f "
        "hit_over_8_length=%d hit_over_8_percent=%.2f",
        total_evaluated_length,
        total_unique_length,
        unique_percent,
        total_hit_up_to_8_length,
        hit_up_to_8_percent,
        total_hit_over_8_length,
        hit_over_8_percent,
    )
    logger.info(
        "Mappability calculation completed: genes=%d distinct_kmers=%d output=%s",
        len(summary_df),
        len(kmer_to_group_ids),
        summary_file,
    )

    return summary_df


def calculate_reference_mappability(
    reference_fasta: str,
    output_prefix: str,
    coordinate_tsv: str,
    read_length: int = 30,
    step: int = 1,
    progress_gene_interval: int = 100,
    write_position_table: bool = False,
    **_: object,
) -> tuple[pd.DataFrame, None]:
    """Compatibility wrapper for the previous caller."""

    if write_position_table:
        raise ValueError("Position-level mappability output has been removed.")

    summary_df = calculate_genomic_mappability(
        reference_fasta=reference_fasta,
        coordinate_tsv=coordinate_tsv,
        output_prefix=output_prefix,
        read_length=read_length,
        step=step,
        progress_reference_interval=progress_gene_interval,
    )
    return summary_df, None