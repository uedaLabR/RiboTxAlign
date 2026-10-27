from __future__ import annotations

import re
from functools import lru_cache


_TRANSCRIPT_TOKEN = re.compile(
    r"^T\d+$",
    re.IGNORECASE,
)
_GROUP_TOKEN = re.compile(
    r"^G\d+$",
    re.IGNORECASE,
)
_GROUP_ID_PATTERN = re.compile(
    r"^(?P<chr_strand>.+_(?:plus|minus))_G\d+$",
    re.IGNORECASE,
)


@lru_cache(maxsize=65_536)
def get_group_id(
    reference_name: str,
) -> str:
    """Extract the overlap-group ID from a final reference name.

    Example::

        chr1_plus_T001_G00001_GAPDH

    becomes::

        chr1_plus_G00001
    """

    parts = reference_name.split(
        "_"
    )

    transcript_token_index = next(
        (
            index
            for index, part
            in enumerate(parts)
            if _TRANSCRIPT_TOKEN.fullmatch(
                part
            )
        ),
        None,
    )

    if transcript_token_index is None:
        raise ValueError(
            "Cannot extract group ID because no T### token was found in "
            f"reference name: {reference_name}"
        )

    group_token_index = (
        transcript_token_index + 1
    )

    if (
        group_token_index
        >= len(parts)
        or not _GROUP_TOKEN.fullmatch(
            parts[group_token_index]
        )
    ):
        raise ValueError(
            "Cannot extract group ID because no G### token follows the "
            f"T### token in reference name: {reference_name}"
        )

    group_parts = (
        parts[
            :transcript_token_index
        ]
        + [
            parts[
                group_token_index
            ]
        ]
    )

    return "_".join(
        group_parts
    )


def get_chr_strand_from_group_id(
    group_id: str,
) -> str:
    """Extract the chromosome/strand prefix from a group ID."""

    match = _GROUP_ID_PATTERN.fullmatch(
        group_id
    )

    if match is None:
        raise ValueError(
            "Cannot extract chromosome/strand from group ID: "
            f"{group_id}"
        )

    return match.group(
        "chr_strand"
    )
