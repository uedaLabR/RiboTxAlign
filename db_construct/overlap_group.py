from __future__ import annotations

import re
from collections import defaultdict

from db_construct.data_class import Transcript
from db_construct.overlap_group_worker import (
    find_groups_in_one_reference_strand,
)
from db_construct.transcript_resolver import (
    resolve_one_overlap_group,
)


def group_transcripts_by_reference_and_strand(
    gencode_transcripts: list[Transcript],
    stringtie_transcripts: list[Transcript],
) -> dict[
    tuple[str, str],
    list[Transcript],
]:
    """Group transcripts by chromosome and strand."""

    grouped_transcripts: dict[
        tuple[str, str],
        list[Transcript],
    ] = defaultdict(list)

    for transcript in (
        gencode_transcripts
        + stringtie_transcripts
    ):
        key = (
            transcript.chrom,
            transcript.strand,
        )

        grouped_transcripts[key].append(
            transcript
        )

    return dict(grouped_transcripts)


def get_strand_label(
    strand: str,
) -> str:
    """Return a readable strand label."""

    if strand == "+":
        return "plus"

    if strand == "-":
        return "minus"

    return "unknown"


def sanitize_id_text(
    value: str,
) -> str:
    """Convert text into a safe ID component."""

    sanitized = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        value,
    )

    sanitized = sanitized.strip(
        "_.-"
    )

    return sanitized or "NA"


def get_overlap_group_gene_label(
    overlap_group: list[Transcript],
    maximum_gene_names: int = 3,
) -> str:
    """Create a readable gene label for one overlap group."""

    gene_names = sorted(
        {
            transcript.gene_name
            for transcript in overlap_group
            if transcript.gene_name
        }
    )

    if not gene_names:
        gene_names = sorted(
            {
                transcript.gene_id
                for transcript in overlap_group
                if transcript.gene_id
            }
        )

    if not gene_names:
        return "NA"

    sanitized_names = [
        sanitize_id_text(name)
        for name in gene_names
    ]

    if (
        len(sanitized_names)
        <= maximum_gene_names
    ):
        return "-".join(
            sanitized_names
        )

    visible_names = sanitized_names[
        :maximum_gene_names
    ]

    hidden_count = (
        len(sanitized_names)
        - maximum_gene_names
    )

    return (
        "-".join(visible_names)
        + f"-plus{hidden_count}"
    )


def create_overlap_group_id(
    chrom: str,
    strand: str,
    group_index: int,
    gene_label: str,
) -> str:
    """Create a readable overlap group ID."""

    strand_label = get_strand_label(
        strand
    )

    return (
        f"{chrom}_"
        f"{strand_label}_"
        f"G{group_index:05d}_"
        f"{gene_label}"
    )


def assign_reference_transcript_ids(
    reference_transcripts: list[Transcript],
    chrom: str,
    strand: str,
    group_index: int,
    gene_label: str,
) -> None:
    """Assign IDs to resolved final reference transcripts."""

    strand_label = get_strand_label(
        strand
    )

    for transcript_index, transcript in enumerate(
        reference_transcripts,
        start=1,
    ):
        previous_reference_id = getattr(
            transcript,
            "reference_id",
            None,
        )

        if previous_reference_id:
            transcript.source_reference_id = str(
                previous_reference_id
            )

            if not getattr(
                transcript,
                "representative_transcript_id",
                None,
            ):
                transcript.representative_transcript_id = str(
                    previous_reference_id
                )
        elif str(
            getattr(
                transcript,
                "source",
                "",
            )
        ).lower().find("gencode") >= 0:
            transcript.representative_transcript_id = str(
                transcript.transcript_id
            )

        transcript.reference_id = (
            f"{chrom}_"
            f"{strand_label}_"
            f"T{transcript_index:03d}_"
            f"G{group_index:05d}_"
            f"{gene_label}"
        )


def build_reference_transcripts(
    gencode_transcripts: list[Transcript],
    stringtie_transcripts: list[Transcript],
) -> list[Transcript]:
    """Build final resolved reference transcripts."""

    grouped_transcripts = (
        group_transcripts_by_reference_and_strand(
            gencode_transcripts=(
                gencode_transcripts
            ),
            stringtie_transcripts=(
                stringtie_transcripts
            ),
        )
    )

    print(
        "Reference-strand groups:",
        len(grouped_transcripts),
    )

    total_transcripts = sum(
        len(transcripts)
        for transcripts
        in grouped_transcripts.values()
    )

    print(
        "Transcripts in groups:",
        total_transcripts,
    )

    final_reference_transcripts: list[
        Transcript
    ] = []

    # Shared across every chromosome/strand and every overlap group, so a
    # GENCODE representative transcript is never added more than once
    # even if its gene's loci end up split across separate overlap
    # groups (e.g. on different strands or in disjoint regions).
    added_gencode_ids: set[str] = set()

    sorted_reference_groups = sorted(
        grouped_transcripts.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
        ),
    )

    for (
        chrom,
        strand,
    ), transcripts in sorted_reference_groups:

        print(
            "Process reference-strand:",
            chrom,
            strand,
            "transcripts:",
            len(transcripts),
        )

        overlap_groups = (
            find_groups_in_one_reference_strand(
                transcripts=transcripts,
                chrom=chrom,
                strand=strand,
            )
        )

        for group_index, overlap_group in enumerate(
            overlap_groups,
            start=1,
        ):
            gene_label = (
                get_overlap_group_gene_label(
                    overlap_group
                )
            )

            overlap_group_id = (
                create_overlap_group_id(
                    chrom=chrom,
                    strand=strand,
                    group_index=group_index,
                    gene_label=gene_label,
                )
            )

            for transcript in overlap_group:
                transcript.overlap_group_id = (
                    overlap_group_id
                )

            resolved_transcripts = (
                resolve_one_overlap_group(
                    overlap_group=overlap_group,
                    overlap_group_id=(
                        overlap_group_id
                    ),
                    added_gencode_ids=(
                        added_gencode_ids
                    ),
                )
            )

            assign_reference_transcript_ids(
                reference_transcripts=(
                    resolved_transcripts
                ),
                chrom=chrom,
                strand=strand,
                group_index=group_index,
                gene_label=gene_label,
            )

            final_reference_transcripts.extend(
                resolved_transcripts
            )

        print(
            "Overlap groups in reference-strand:",
            len(overlap_groups),
        )

    print(
        "Final reference transcripts:",
        len(final_reference_transcripts),
    )

    return final_reference_transcripts
