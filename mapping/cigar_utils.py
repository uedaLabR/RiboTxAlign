from __future__ import annotations

import pysam


REFERENCE_CONSUMING_OPERATIONS = {
    pysam.CMATCH,
    pysam.CDEL,
    pysam.CREF_SKIP,
    pysam.CEQUAL,
    pysam.CDIFF,
}


def merge_adjacent_cigar_operations(
    cigartuples: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """
    Merge adjacent CIGAR operations of the same type.
    """

    merged: list[tuple[int, int]] = []

    for operation, length in cigartuples:
        if length <= 0:
            continue

        if (
            merged
            and merged[-1][0] == operation
        ):
            previous_operation, previous_length = (
                merged[-1]
            )

            merged[-1] = (
                previous_operation,
                previous_length + length,
            )
        else:
            merged.append(
                (
                    operation,
                    length,
                )
            )

    return merged


def calculate_soft_clip_lengths(
    query_start: int,
    query_end: int,
    read_length: int,
    is_reverse: bool,
) -> tuple[int, int]:
    """
    Calculate BAM left and right soft clips.

    mappy q_st/q_en are in the original query orientation.

    For reverse alignments, BAM stores the query in reverse-
    complemented orientation, so the clip sides are exchanged.
    """

    if query_start < 0:
        raise ValueError(
            "query_start must be non-negative"
        )

    if query_end < query_start:
        raise ValueError(
            "query_end must be >= query_start"
        )

    if query_end > read_length:
        raise ValueError(
            "query_end exceeds read length"
        )

    if is_reverse:
        left_clip = (
            read_length
            - query_end
        )

        right_clip = query_start
    else:
        left_clip = query_start

        right_clip = (
            read_length
            - query_end
        )

    return (
        left_clip,
        right_clip,
    )


def add_soft_clips(
    cigartuples: list[tuple[int, int]],
    query_start: int,
    query_end: int,
    read_length: int,
    is_reverse: bool,
) -> list[tuple[int, int]]:
    """
    Add terminal soft clips to a mappy CIGAR.
    """

    left_clip, right_clip = (
        calculate_soft_clip_lengths(
            query_start=query_start,
            query_end=query_end,
            read_length=read_length,
            is_reverse=is_reverse,
        )
    )

    result: list[tuple[int, int]] = []

    if left_clip > 0:
        result.append(
            (
                pysam.CSOFT_CLIP,
                left_clip,
            )
        )

    result.extend(
        cigartuples
    )

    if right_clip > 0:
        result.append(
            (
                pysam.CSOFT_CLIP,
                right_clip,
            )
        )

    return merge_adjacent_cigar_operations(
        result
    )


def cigar_reference_length(
    cigartuples: list[tuple[int, int]],
) -> int:
    """
    Return the number of reference bases consumed.
    """

    return sum(
        length
        for operation, length
        in cigartuples
        if operation
        in REFERENCE_CONSUMING_OPERATIONS
    )
