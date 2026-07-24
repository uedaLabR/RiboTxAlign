from __future__ import annotations

from collections import defaultdict

from db_construct.data_class import GeneGroup, Interval, Transcript


def strip_version(identifier: str | None) -> str | None:
    if identifier is None:
        return None
    return identifier.split(".", 1)[0]


def group_gencode_transcripts(
    transcripts: dict[str, Transcript],
) -> dict[str, GeneGroup]:
    gene_groups: dict[str, GeneGroup] = {}
    for transcript in transcripts.values():
        if transcript.gene_id is None:
            continue
        gene_id = transcript.gene_id
        if gene_id not in gene_groups:
            gene_groups[gene_id] = GeneGroup(
                gene_id=gene_id,
                chrom=transcript.chrom,
                strand=transcript.strand,
                gene_name=transcript.gene_name,
            )
        gene_groups[gene_id].add_transcript(transcript)
    return gene_groups


def build_gene_id_index(
    gene_groups: dict[str, GeneGroup],
) -> dict[str, str]:
    gene_id_index: dict[str, str] = {}
    for gene_id in gene_groups:
        gene_id_index[gene_id] = gene_id
        base_id = strip_version(gene_id)
        if base_id is not None:
            gene_id_index.setdefault(base_id, gene_id)
    return gene_id_index


def add_stringtie_by_reference_gene(
    gene_groups: dict[str, GeneGroup],
    transcripts: dict[str, Transcript],
) -> list[Transcript]:
    gene_id_index = build_gene_id_index(gene_groups)
    remaining: list[Transcript] = []

    for transcript in transcripts.values():
        reference_gene_id = transcript.reference_gene_id
        canonical_gene_id = None
        if reference_gene_id is not None:
            canonical_gene_id = gene_id_index.get(reference_gene_id)
            if canonical_gene_id is None:
                canonical_gene_id = gene_id_index.get(strip_version(reference_gene_id))

        if canonical_gene_id is None:
            remaining.append(transcript)
            continue

        gene_groups[canonical_gene_id].add_transcript(transcript)

    return remaining


def merge_intervals(intervals: list[Interval]) -> list[Interval]:
    if not intervals:
        return []
    sorted_intervals = sorted(intervals)
    merged: list[list[int]] = [[sorted_intervals[0][0], sorted_intervals[0][1]]]

    for start, end in sorted_intervals[1:]:
        previous = merged[-1]
        if start <= previous[1]:
            previous[1] = max(previous[1], end)
        else:
            merged.append([start, end])

    return [(start, end) for start, end in merged]


def interval_overlap_length(
    first: list[Interval],
    second: list[Interval],
) -> int:
    first = merge_intervals(first)
    second = merge_intervals(second)
    i = 0
    j = 0
    overlap = 0

    while i < len(first) and j < len(second):
        first_start, first_end = first[i]
        second_start, second_end = second[j]
        overlap += max(
            0,
            min(first_end, second_end) - max(first_start, second_start),
        )
        if first_end <= second_end:
            i += 1
        else:
            j += 1

    return overlap


def build_gene_regions(
    gene_groups: dict[str, GeneGroup],
) -> dict[tuple[str, str], list[tuple[int, int, str, list[Interval]]]]:
    regions: dict[
        tuple[str, str],
        list[tuple[int, int, str, list[Interval]]],
    ] = defaultdict(list)

    for gene_id, gene_group in gene_groups.items():
        gene_exons: list[Interval] = []
        for transcript in gene_group.transcripts.values():
            gene_exons.extend(transcript.exons)
        gene_exons = merge_intervals(gene_exons)
        if not gene_exons:
            continue
        start = min(interval[0] for interval in gene_exons)
        end = max(interval[1] for interval in gene_exons)
        regions[(gene_group.chrom, gene_group.strand)].append(
            (start, end, gene_id, gene_exons)
        )

    for key in regions:
        regions[key].sort(key=lambda value: value[0])

    return regions


def assign_transcripts_by_overlap(
    gene_groups: dict[str, GeneGroup],
    transcripts: list[Transcript],
) -> list[Transcript]:
    gene_regions = build_gene_regions(gene_groups)
    transcripts_by_key: dict[tuple[str, str], list[Transcript]] = defaultdict(list)

    for transcript in transcripts:
        transcripts_by_key[(transcript.chrom, transcript.strand)].append(transcript)

    unresolved: list[Transcript] = []

    for key, key_transcripts in transcripts_by_key.items():
        regions = gene_regions.get(key, [])
        if not regions:
            unresolved.extend(key_transcripts)
            continue

        key_transcripts.sort(key=lambda transcript: transcript.start)
        active_regions: list[tuple[int, int, str, list[Interval]]] = []
        region_index = 0

        for transcript in key_transcripts:
            transcript_intervals = transcript.exons or [(transcript.start, transcript.end)]

            while region_index < len(regions) and regions[region_index][0] < transcript.end:
                active_regions.append(regions[region_index])
                region_index += 1

            active_regions = [region for region in active_regions if region[1] > transcript.start]
            scored_candidates: list[tuple[int, str]] = []

            for _, _, gene_id, gene_exons in active_regions:
                overlap = interval_overlap_length(transcript_intervals, gene_exons)
                if overlap > 0:
                    scored_candidates.append((overlap, gene_id))

            if not scored_candidates:
                unresolved.append(transcript)
                continue

            scored_candidates.sort(reverse=True)
            best_overlap, best_gene_id = scored_candidates[0]
            is_ambiguous = (
                len(scored_candidates) > 1
                and scored_candidates[1][0] == best_overlap
            )

            if is_ambiguous:
                unresolved.append(transcript)
                continue

            gene_groups[best_gene_id].add_transcript(transcript)

    return unresolved


def add_novel_stringtie_groups(
    gene_groups: dict[str, GeneGroup],
    transcripts: list[Transcript],
) -> None:
    for transcript in transcripts:
        gene_id = transcript.gene_id or f"NOVEL:{transcript.transcript_id}"
        if gene_id not in gene_groups:
            gene_groups[gene_id] = GeneGroup(
                gene_id=gene_id,
                chrom=transcript.chrom,
                strand=transcript.strand,
                gene_name=transcript.gene_name,
            )
        gene_groups[gene_id].add_transcript(transcript)


def merge_transcripts_by_gene(
    gencode_transcripts: dict[str, Transcript],
    stringtie_transcripts: dict[str, Transcript],
) -> tuple[dict[str, GeneGroup], list[Transcript]]:
    gene_groups = group_gencode_transcripts(gencode_transcripts)
    remaining = add_stringtie_by_reference_gene(
        gene_groups=gene_groups,
        transcripts=stringtie_transcripts,
    )
    novel_transcripts = assign_transcripts_by_overlap(
        gene_groups=gene_groups,
        transcripts=remaining,
    )
    add_novel_stringtie_groups(
        gene_groups=gene_groups,
        transcripts=novel_transcripts,
    )
    return gene_groups, novel_transcripts
