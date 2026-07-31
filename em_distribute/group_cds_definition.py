from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO


@dataclass(slots=True)
class ReferenceCdsDefinition:
    """CDS intervals associated with one mapping reference."""

    gencode_transcript_ids: set[str] = field(
        default_factory=set
    )
    cds_intervals: list[
        tuple[
            str,
            str,
            int,
            int,
        ]
    ] = field(
        default_factory=list
    )


GroupCdsDefinitions = dict[
    str,
    dict[
        str,
        ReferenceCdsDefinition,
    ],
]


def _open_text(
    path: str | Path,
) -> TextIO:
    """Open a plain-text or gzip-compressed file."""

    file_path = Path(path)

    if file_path.suffix == ".gz":
        return gzip.open(
            file_path,
            mode="rt",
            encoding="utf-8",
        )

    return file_path.open(
        mode="r",
        encoding="utf-8",
    )


def _normalize_transcript_id(
    value: str | None,
) -> str | None:
    """Remove common prefixes and version suffixes."""

    if value is None:
        return None

    transcript_id = value.strip().strip('"')

    if not transcript_id:
        return None

    for prefix in (
        "transcript:",
        "rna:",
    ):
        if transcript_id.startswith(prefix):
            transcript_id = transcript_id[
                len(prefix):
            ]
            break

    return transcript_id.split(
        ".",
        maxsplit=1,
    )[0]


def _split_transcript_ids(
    value: str | None,
) -> set[str]:
    """Split one metadata field containing transcript IDs."""

    if value is None:
        return set()

    text = value.strip()

    if not text:
        return set()

    for separator in (
        ";",
        "|",
    ):
        text = text.replace(
            separator,
            ",",
        )

    transcript_ids: set[str] = set()

    for item in text.split(","):
        transcript_id = _normalize_transcript_id(
            item
        )

        if transcript_id:
            transcript_ids.add(
                transcript_id
            )

    return transcript_ids


def _parse_gtf_attributes(
    attribute_text: str,
) -> dict[str, str]:
    """Parse one GTF attribute field."""

    attributes: dict[str, str] = {}

    for item in attribute_text.split(";"):
        item = item.strip()

        if not item:
            continue

        parts = item.split(
            None,
            maxsplit=1,
        )

        if len(parts) != 2:
            continue

        key, value = parts

        attributes[key] = (
            value.strip().strip('"')
        )

    return attributes


def load_reference_definitions(
    reference_metadata_path: str | Path,
) -> GroupCdsDefinitions:
    """
    Load group and mapping-reference relations from metadata.

    Required columns:

        reference_id
        overlap_group_id

    GENCODE transcript IDs are collected from:

        representative_transcript_id
        gencode_transcript_ids
        source_reference_id
    """

    group_definitions: GroupCdsDefinitions = {}

    with _open_text(
        reference_metadata_path
    ) as input_handle:
        reader = csv.DictReader(
            input_handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise ValueError(
                "Reference metadata has no header."
            )

        required_columns = {
            "reference_id",
            "overlap_group_id",
        }

        missing_columns = (
            required_columns
            - set(reader.fieldnames)
        )

        if missing_columns:
            raise ValueError(
                "Reference metadata is missing columns: "
                + ", ".join(
                    sorted(missing_columns)
                )
            )

        for row in reader:
            reference_id = (
                row.get("reference_id")
                or ""
            ).strip()

            group_id = (
                row.get("overlap_group_id")
                or ""
            ).strip()

            if not reference_id or not group_id:
                continue

            gencode_transcript_ids: set[str] = set()

            gencode_transcript_ids.update(
                _split_transcript_ids(
                    row.get(
                        "representative_transcript_id"
                    )
                )
            )

            gencode_transcript_ids.update(
                _split_transcript_ids(
                    row.get(
                        "gencode_transcript_ids"
                    )
                )
            )

            source_reference_id = (
                _normalize_transcript_id(
                    row.get(
                        "source_reference_id"
                    )
                )
            )

            if (
                source_reference_id is not None
                and source_reference_id.startswith(
                    "ENST"
                )
            ):
                gencode_transcript_ids.add(
                    source_reference_id
                )

            group_definitions.setdefault(
                group_id,
                {},
            )[reference_id] = (
                ReferenceCdsDefinition(
                    gencode_transcript_ids=(
                        gencode_transcript_ids
                    )
                )
            )

    return group_definitions


def collect_needed_gencode_transcripts(
    group_definitions: GroupCdsDefinitions,
) -> set[str]:
    """Collect GENCODE transcript IDs referenced by metadata."""

    transcript_ids: set[str] = set()

    for reference_definitions in (
        group_definitions.values()
    ):
        for definition in (
            reference_definitions.values()
        ):
            transcript_ids.update(
                definition.gencode_transcript_ids
            )

    return transcript_ids


def load_gencode_cds_intervals(
    gencode_gtf_path: str | Path,
    needed_transcript_ids: set[str],
) -> dict[
    str,
    list[
        tuple[
            str,
            str,
            int,
            int,
        ]
    ],
]:
    """
    Load CDS intervals only for transcripts referenced by metadata.

    Coordinates are retained as one-based inclusive genomic intervals.
    """

    cds_by_transcript: dict[
        str,
        list[
            tuple[
                str,
                str,
                int,
                int,
            ]
        ],
    ] = defaultdict(list)

    if not needed_transcript_ids:
        return {}

    with _open_text(
        gencode_gtf_path
    ) as input_handle:
        for line in input_handle:
            if (
                not line
                or line.startswith("#")
            ):
                continue

            fields = line.rstrip("\n").split(
                "\t"
            )

            if (
                len(fields) != 9
                or fields[2] != "CDS"
            ):
                continue

            attributes = _parse_gtf_attributes(
                fields[8]
            )

            transcript_id = (
                _normalize_transcript_id(
                    attributes.get(
                        "transcript_id"
                    )
                )
            )

            if (
                transcript_id is None
                or transcript_id
                not in needed_transcript_ids
            ):
                continue

            cds_by_transcript[
                transcript_id
            ].append(
                (
                    fields[0],
                    fields[6],
                    int(fields[3]),
                    int(fields[4]),
                )
            )

    return dict(cds_by_transcript)


def attach_cds_intervals(
    group_definitions: GroupCdsDefinitions,
    cds_by_transcript: dict[
        str,
        list[
            tuple[
                str,
                str,
                int,
                int,
            ]
        ],
    ],
) -> None:
    """Attach GENCODE CDS intervals to each mapping reference."""

    for reference_definitions in (
        group_definitions.values()
    ):
        for definition in (
            reference_definitions.values()
        ):
            interval_set: set[
                tuple[
                    str,
                    str,
                    int,
                    int,
                ]
            ] = set()

            for transcript_id in (
                definition.gencode_transcript_ids
            ):
                interval_set.update(
                    cds_by_transcript.get(
                        transcript_id,
                        [],
                    )
                )

            definition.cds_intervals = sorted(
                interval_set,
                key=lambda interval: (
                    interval[0],
                    interval[1],
                    interval[2],
                    interval[3],
                ),
            )


def create_group_cds_definitions(
    gencode_gtf_path: str | Path,
    reference_metadata_path: str | Path,
) -> GroupCdsDefinitions:
    """
    Create all group/reference CDS definitions before worker startup.

    This function performs no mapping-Parquet processing.
    """

    group_definitions = (
        load_reference_definitions(
            reference_metadata_path
        )
    )

    needed_transcript_ids = (
        collect_needed_gencode_transcripts(
            group_definitions
        )
    )

    cds_by_transcript = (
        load_gencode_cds_intervals(
            gencode_gtf_path=(
                gencode_gtf_path
            ),
            needed_transcript_ids=(
                needed_transcript_ids
            ),
        )
    )

    attach_cds_intervals(
        group_definitions=group_definitions,
        cds_by_transcript=cds_by_transcript,
    )

    return group_definitions