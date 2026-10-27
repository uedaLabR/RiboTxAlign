from __future__ import annotations

from collections import defaultdict

from db_construct.data_class import Transcript
from db_construct.stringtie_subset_resolver import (
    collapse_stringtie_subsets,
    get_gene_key,
)


def group_by_gene(
    transcripts: list[Transcript],
) -> dict[str, list[Transcript]]:
    """Group transcripts by the same key used by the collapse function."""

    grouped: dict[
        str,
        list[Transcript],
    ] = defaultdict(list)

    for transcript in transcripts:
        grouped[
            get_gene_key(
                transcript
            )
        ].append(
            transcript
        )

    return dict(grouped)


def format_exons(
    transcript: Transcript,
) -> str:
    """Return a compact exon-chain string."""

    return ",".join(
        f"{start}-{end}"
        for start, end in transcript.exons
    )


def print_gene_example(
    gene_key: str,
    before: list[Transcript],
    after: list[Transcript],
) -> None:
    """Print one gene-level before/after example."""

    gene_name = next(
        (
            transcript.gene_name
            for transcript in before
            if transcript.gene_name
        ),
        "NA",
    )

    print()
    print(
        f"Gene={gene_key} "
        f"name={gene_name} "
        f"before={len(before)} "
        f"after={len(after)}"
    )

    print("  Before")

    for transcript in sorted(
        before,
        key=lambda item: (
            item.transcript_id
        ),
    ):
        print(
            "   ",
            transcript.transcript_id,
            format_exons(
                transcript
            ),
        )

    print("  After")

    for transcript in sorted(
        after,
        key=lambda item: (
            item.transcript_id
        ),
    ):
        source_ids = getattr(
            transcript,
            "source_transcript_ids",
            [
                transcript.transcript_id
            ],
        )

        print(
            "   ",
            transcript.transcript_id,
            "members=",
            ",".join(source_ids),
            format_exons(
                transcript
            ),
        )


def print_subset_collapse_examples(
    stringtie_transcripts: list[Transcript],
    maximum_examples: int = 3,
) -> list[Transcript]:
    """
    Collapse transcripts and print a few debug examples.

    Category 1:
        Genes with at least five transcripts before collapse and at most
        two transcripts after collapse.

    Category 2:
        Genes with at least five transcripts before collapse and no
        reduction after collapse.
    """

    collapsed = collapse_stringtie_subsets(
        stringtie_transcripts
    )

    before_by_gene = group_by_gene(
        stringtie_transcripts
    )
    after_by_gene = group_by_gene(
        collapsed
    )

    strongly_collapsed: list[
        tuple[
            str,
            list[Transcript],
            list[Transcript],
        ]
    ] = []

    unchanged: list[
        tuple[
            str,
            list[Transcript],
            list[Transcript],
        ]
    ] = []

    for gene_key, before in (
        before_by_gene.items()
    ):
        if len(before) < 5:
            continue

        after = after_by_gene.get(
            gene_key,
            [],
        )

        if len(after) <= 2:
            strongly_collapsed.append(
                (
                    gene_key,
                    before,
                    after,
                )
            )
        elif len(after) == len(before):
            unchanged.append(
                (
                    gene_key,
                    before,
                    after,
                )
            )

    strongly_collapsed.sort(
        key=lambda item: (
            -len(item[1]),
            item[0],
        )
    )
    unchanged.sort(
        key=lambda item: (
            -len(item[1]),
            item[0],
        )
    )

    print()
    print(
        "=== Collapsed from >=5 to <=2 ==="
    )

    for (
        gene_key,
        before,
        after,
    ) in strongly_collapsed[
        :maximum_examples
    ]:
        print_gene_example(
            gene_key=gene_key,
            before=before,
            after=after,
        )

    print()
    print(
        "=== Unchanged with >=5 transcripts ==="
    )

    for (
        gene_key,
        before,
        after,
    ) in unchanged[
        :maximum_examples
    ]:
        print_gene_example(
            gene_key=gene_key,
            before=before,
            after=after,
        )

    return collapsed
