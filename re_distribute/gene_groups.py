from __future__ import annotations

import csv
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


GenomicInterval = tuple[int, int]


@dataclass(frozen=True)
class ReferenceSpan:
    """Genomic span covered by one synthetic transcript reference."""

    reference_id: str
    gene_id: str
    chrom: str
    strand: str
    genome_start: int
    genome_end: int


@dataclass(frozen=True)
class GeneGroup:
    """One strand-specific genomic group formed from overlapping references."""

    group_id: str
    chrom: str
    strand: str
    genome_start: int
    genome_end: int
    gene_ids: tuple[str, ...]
    reference_ids: tuple[str, ...]
    union_intervals: tuple[GenomicInterval, ...]
    union_length: int


class GeneGroupIndex:
    """Fast lookup and interval metadata for strand-specific gene groups."""

    def __init__(
        self,
        groups: list[GeneGroup],
        gene_intervals: dict[str, tuple[GenomicInterval, ...]],
    ):
        self.groups = groups
        self.group_by_id = {
            group.group_id: group
            for group in groups
        }
        self.gene_intervals = gene_intervals
        self.group_id_by_gene: dict[str, str] = {}
        self._groups_by_chrom_strand: dict[
            tuple[str, str],
            list[GeneGroup],
        ] = {}
        self._starts_by_chrom_strand: dict[
            tuple[str, str],
            list[int],
        ] = {}

        for group in groups:
            key = (group.chrom, group.strand)
            self._groups_by_chrom_strand.setdefault(
                key,
                [],
            ).append(group)

            for gene_id in group.gene_ids:
                previous_group_id = self.group_id_by_gene.get(gene_id)
                if (
                    previous_group_id is not None
                    and previous_group_id != group.group_id
                ):
                    raise ValueError(
                        "One gene belongs to multiple non-overlapping gene "
                        f"groups: {gene_id}"
                    )
                self.group_id_by_gene[gene_id] = group.group_id

        for key, key_groups in self._groups_by_chrom_strand.items():
            key_groups.sort(
                key=lambda item: (
                    item.genome_start,
                    item.genome_end,
                    item.group_id,
                )
            )
            self._starts_by_chrom_strand[key] = [
                group.genome_start
                for group in key_groups
            ]

    def find_containing_group(
        self,
        chrom: str,
        strand: str,
        genome_start: int,
        genome_end: int,
    ) -> GeneGroup | None:
        """Return the group that fully contains the genomic interval."""

        key = (chrom, strand)
        groups = self._groups_by_chrom_strand.get(key)
        starts = self._starts_by_chrom_strand.get(key)

        if not groups or not starts:
            return None

        index = bisect_right(starts, genome_start) - 1
        if index < 0:
            return None

        group = groups[index]
        if (
            genome_start >= group.genome_start
            and genome_end <= group.genome_end
        ):
            return group

        return None


def merge_intervals(
    intervals: list[GenomicInterval],
) -> tuple[GenomicInterval, ...]:
    """Merge overlapping or touching half-open genomic intervals."""

    if not intervals:
        return ()

    sorted_intervals = sorted(intervals)
    merged: list[list[int]] = [
        [sorted_intervals[0][0], sorted_intervals[0][1]]
    ]

    for start, end in sorted_intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])

    return tuple((start, end) for start, end in merged)


def load_coordinate_metadata(
    coordinate_tsv: str,
) -> tuple[
    list[ReferenceSpan],
    dict[str, tuple[GenomicInterval, ...]],
]:
    """Load reference spans and merged exon/block intervals for each gene."""

    coordinate_path = Path(coordinate_tsv)
    if not coordinate_path.exists():
        raise FileNotFoundError(coordinate_tsv)

    rows_by_reference: dict[
        str,
        list[tuple[str, str, str, int, int]],
    ] = {}
    raw_gene_intervals: dict[str, list[GenomicInterval]] = {}
    gene_locus: dict[str, tuple[str, str]] = {}

    with coordinate_path.open(
        encoding="utf-8",
    ) as input_file:
        reader = csv.DictReader(
            input_file,
            delimiter="\t",
        )

        required_columns = {
            "reference_id",
            "gene_id",
            "chrom",
            "genome_start",
            "genome_end",
            "strand",
        }
        missing_columns = (
            required_columns
            - set(reader.fieldnames or [])
        )
        if missing_columns:
            raise ValueError(
                "Missing coordinate-table columns: "
                + ", ".join(sorted(missing_columns))
            )

        for row in reader:
            reference_id = row["reference_id"]
            gene_id = row["gene_id"]
            chrom = row["chrom"]
            strand = row["strand"]
            genome_start = int(row["genome_start"])
            genome_end = int(row["genome_end"])

            if strand not in {"+", "-"}:
                raise ValueError(
                    f"Invalid strand for {reference_id}: {strand}"
                )
            if genome_end <= genome_start:
                raise ValueError(
                    f"Invalid genomic interval for {reference_id}: "
                    f"{genome_start}-{genome_end}"
                )

            locus = (chrom, strand)
            previous_locus = gene_locus.get(gene_id)
            if previous_locus is not None and previous_locus != locus:
                raise ValueError(
                    f"Gene maps to multiple loci: {gene_id}"
                )
            gene_locus[gene_id] = locus

            rows_by_reference.setdefault(
                reference_id,
                [],
            ).append(
                (
                    gene_id,
                    chrom,
                    strand,
                    genome_start,
                    genome_end,
                )
            )
            raw_gene_intervals.setdefault(
                gene_id,
                [],
            ).append(
                (genome_start, genome_end)
            )

    spans: list[ReferenceSpan] = []

    for reference_id, rows in rows_by_reference.items():
        gene_ids = {row[0] for row in rows}
        chroms = {row[1] for row in rows}
        strands = {row[2] for row in rows}

        if len(gene_ids) != 1:
            raise ValueError(
                f"Reference maps to multiple genes: {reference_id}"
            )
        if len(chroms) != 1:
            raise ValueError(
                f"Reference maps to multiple chromosomes: {reference_id}"
            )
        if len(strands) != 1:
            raise ValueError(
                f"Reference has mixed strands: {reference_id}"
            )

        spans.append(
            ReferenceSpan(
                reference_id=reference_id,
                gene_id=next(iter(gene_ids)),
                chrom=next(iter(chroms)),
                strand=next(iter(strands)),
                genome_start=min(row[3] for row in rows),
                genome_end=max(row[4] for row in rows),
            )
        )

    gene_intervals = {
        gene_id: merge_intervals(intervals)
        for gene_id, intervals in raw_gene_intervals.items()
    }

    return spans, gene_intervals


def load_reference_spans(
    coordinate_tsv: str,
) -> list[ReferenceSpan]:
    """Compatibility wrapper returning one span per synthetic reference."""

    spans, _ = load_coordinate_metadata(coordinate_tsv)
    return spans


def build_gene_groups(
    coordinate_tsv: str,
) -> GeneGroupIndex:
    """
    Build transitive overlap groups from synthetic transcript spans.

    Gene-group formation still uses the inexpensive outer reference spans.
    Separately, all coordinate-table exon/block intervals are merged within
    each group. These merged intervals are used only for TPM length correction
    and do not alter the mappability calculation.
    """

    spans, gene_intervals = load_coordinate_metadata(coordinate_tsv)
    spans_by_key: dict[
        tuple[str, str],
        list[ReferenceSpan],
    ] = {}

    for span in spans:
        spans_by_key.setdefault(
            (span.chrom, span.strand),
            [],
        ).append(span)

    grouped_components: list[
        tuple[str, str, int, int, set[str], set[str]]
    ] = []

    for (chrom, strand), key_spans in spans_by_key.items():
        key_spans.sort(
            key=lambda item: (
                item.genome_start,
                item.genome_end,
                item.reference_id,
            )
        )

        current_start: int | None = None
        current_end: int | None = None
        current_gene_ids: set[str] = set()
        current_reference_ids: set[str] = set()

        for span in key_spans:
            if (
                current_start is None
                or current_end is None
                or span.genome_start >= current_end
            ):
                if current_start is not None and current_end is not None:
                    grouped_components.append(
                        (
                            chrom,
                            strand,
                            current_start,
                            current_end,
                            current_gene_ids,
                            current_reference_ids,
                        )
                    )

                current_start = span.genome_start
                current_end = span.genome_end
                current_gene_ids = {span.gene_id}
                current_reference_ids = {span.reference_id}
            else:
                current_end = max(current_end, span.genome_end)
                current_gene_ids.add(span.gene_id)
                current_reference_ids.add(span.reference_id)

        if current_start is not None and current_end is not None:
            grouped_components.append(
                (
                    chrom,
                    strand,
                    current_start,
                    current_end,
                    current_gene_ids,
                    current_reference_ids,
                )
            )

    grouped_components.sort(
        key=lambda item: (
            item[0],
            item[1],
            item[2],
            item[3],
        )
    )

    groups: list[GeneGroup] = []

    for group_number, component in enumerate(
        grouped_components,
        start=1,
    ):
        (
            chrom,
            strand,
            genome_start,
            genome_end,
            gene_ids,
            reference_ids,
        ) = component

        group_intervals = merge_intervals(
            [
                interval
                for gene_id in gene_ids
                for interval in gene_intervals.get(gene_id, ())
            ]
        )
        union_length = sum(
            end - start
            for start, end in group_intervals
        )

        groups.append(
            GeneGroup(
                group_id=f"GG{group_number:08d}",
                chrom=chrom,
                strand=strand,
                genome_start=genome_start,
                genome_end=genome_end,
                gene_ids=tuple(sorted(gene_ids)),
                reference_ids=tuple(sorted(reference_ids)),
                union_intervals=group_intervals,
                union_length=union_length,
            )
        )

    return GeneGroupIndex(
        groups=groups,
        gene_intervals=gene_intervals,
    )


def gene_groups_to_dataframe(
    group_index: GeneGroupIndex,
) -> pd.DataFrame:
    """Convert gene-group definitions to a tabular representation."""

    return pd.DataFrame(
        [
            {
                "gene_group_id": group.group_id,
                "chrom": group.chrom,
                "strand": group.strand,
                "genome_start": group.genome_start,
                "genome_end": group.genome_end,
                "gene_count": len(group.gene_ids),
                "reference_count": len(group.reference_ids),
                "gene_ids": ",".join(group.gene_ids),
                "reference_ids": ",".join(group.reference_ids),
                "union_interval_count": len(group.union_intervals),
                "union_gene_length": group.union_length,
                "union_intervals": ",".join(
                    f"{start}-{end}"
                    for start, end in group.union_intervals
                ),
            }
            for group in group_index.groups
        ]
    )
