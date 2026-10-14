from __future__ import annotations

from db_construct.data_class import Transcript
from db_construct.stringtie_subset_resolver import (
    collapse_stringtie_subsets,
)


def exons_overlap(
    first: Transcript,
    second: Transcript,
) -> bool:
    """Return True when any exon pair overlaps."""

    first_exons = sorted(
        first.exons
    )
    second_exons = sorted(
        second.exons
    )

    first_index = 0
    second_index = 0

    while (
        first_index < len(first_exons)
        and second_index < len(second_exons)
    ):
        first_start, first_end = (
            first_exons[first_index]
        )
        second_start, second_end = (
            second_exons[second_index]
        )

        if (
            first_start < second_end
            and second_start < first_end
        ):
            return True

        if first_end <= second_start:
            first_index += 1
        else:
            second_index += 1

    return False


def find_groups_in_one_reference_strand(
    transcripts: list[Transcript],
    chrom: str,
    strand: str,
) -> list[list[Transcript]]:
    """
    Find transitive exon-overlap groups in one chromosome-strand set.

    StringTie models are collapsed conservatively before overlap grouping:
      - exact exon structure;
      - the same inferred CDS structure, using the outermost GENCODE CDS
        window for a uniquely assigned gene;
      - or, when CDS is not uniquely inferable, differences restricted to
        the biological first/last exon outer boundaries.

    The longest compatible StringTie structure is retained and the TPM of
    collapsed StringTie members is summed.
    """

    gencode_transcripts = [
        transcript
        for transcript in transcripts
        if transcript.source == "gencode"
    ]

    stringtie_transcripts = [
        transcript
        for transcript in transcripts
        if transcript.source == "stringtie"
    ]

    collapsed_stringtie = (
        collapse_stringtie_subsets(
            stringtie_transcripts=(
                stringtie_transcripts
            ),
            gencode_transcripts=(
                gencode_transcripts
            ),
        )
    )

    working_transcripts = sorted(
        gencode_transcripts
        + collapsed_stringtie,
        key=lambda transcript: (
            transcript.start
            if transcript.start is not None
            else -1,
            transcript.end
            if transcript.end is not None
            else -1,
            transcript.transcript_id,
        ),
    )

    working_transcripts = [
        transcript for transcript in working_transcripts
        if transcript.start is not None
    ]
    parents = list(range(len(working_transcripts)))
    sizes = [1] * len(working_transcripts)

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(first: int, second: int) -> None:
        first, second = find(first), find(second)
        if first == second:
            return
        if sizes[first] < sizes[second]:
            first, second = second, first
        parents[second] = first
        sizes[first] += sizes[second]

    exon_intervals = sorted(
        (start, end, index)
        for index, transcript in enumerate(working_transcripts)
        for start, end in transcript.exons
        if start < end
    )
    rightmost_end = None
    rightmost_index = None
    for start, end, index in exon_intervals:
        # The rightmost-reaching exon overlaps every next exon that starts
        # before its end. Union transcript IDs, so their separate exons can
        # connect distant interval clusters without merging intervening loci.
        # Exons are half-open: touching endpoints do not count as overlap.
        if rightmost_end is not None and start < rightmost_end:
            union(index, rightmost_index)
            if end > rightmost_end:
                rightmost_end, rightmost_index = end, index
        else:
            rightmost_end, rightmost_index = end, index

    components: dict[int, list[Transcript]] = {}
    for index, transcript in enumerate(working_transcripts):
        components.setdefault(find(index), []).append(transcript)
    # Preserve genomic/transcript-ID ordering for stable group numbering.
    overlap_groups = list(components.values())

    print(
        "StringTie conservative collapse:",
        chrom,
        strand,
        "before=",
        len(stringtie_transcripts),
        "after=",
        len(collapsed_stringtie),
    )

    return overlap_groups
