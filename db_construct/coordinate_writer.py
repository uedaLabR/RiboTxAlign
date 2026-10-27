from __future__ import annotations

import csv
from pathlib import Path

from db_construct.data_class import Transcript
from db_construct.gtf_writer import (
    get_output_gene_id,
    get_output_transcript_id,
    normalize_exons,
)


COORDINATE_COLUMNS = [
    "reference_id",
    "gene_id",
    "record_type",
    "segment_index",
    "reference_start",
    "reference_end",
    "chrom",
    "genome_start",
    "genome_end",
    "strand",
    "transcript_id",
    "block_id",
    "overlap_group_id",
]


def get_exons_in_transcript_order(
    transcript: Transcript,
) -> list[tuple[int, int]]:
    """Return valid exons in biological 5'-to-3' order."""

    exons = normalize_exons(
        transcript
    )

    if transcript.strand == "-":
        exons.reverse()

    return exons


def iter_coordinate_rows(
    transcript: Transcript,
):
    """Yield transcript-to-genome mapping rows for one reference."""

    reference_id = get_output_transcript_id(
        transcript
    )
    gene_id = get_output_gene_id(
        transcript
    )
    overlap_group_id = str(
        getattr(
            transcript,
            "overlap_group_id",
            "",
        )
        or ""
    )

    transcript_position = 0

    for segment_index, (
        genome_start,
        genome_end,
    ) in enumerate(
        get_exons_in_transcript_order(
            transcript
        ),
        start=1,
    ):
        segment_length = (
            genome_end
            - genome_start
        )
        reference_end = (
            transcript_position
            + segment_length
        )

        yield {
            "reference_id": reference_id,
            "gene_id": gene_id,
            "record_type": "TRANSCRIPT",
            "segment_index": segment_index,
            "reference_start": transcript_position,
            "reference_end": reference_end,
            "chrom": transcript.chrom,
            "genome_start": genome_start,
            "genome_end": genome_end,
            "strand": transcript.strand,
            "transcript_id": reference_id,
            "block_id": "",
            "overlap_group_id": overlap_group_id,
        }

        transcript_position = (
            reference_end
        )


def write_coordinate_table(
    reference_transcripts: list[Transcript],
    output_path: str | Path,
) -> None:
    """
    Write transcript-to-genome coordinate intervals.

    All transcript and genomic coordinates are 0-based and half-open.
    Exon segments are written in biological transcript order, including
    reverse genomic order for minus-strand transcripts.
    """

    output_path = Path(
        output_path
    )
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    row_count = 0
    transcript_count = 0
    skipped_without_exons = 0
    seen_reference_ids: set[str] = set()

    sorted_transcripts = sorted(
        reference_transcripts,
        key=lambda transcript: (
            transcript.chrom,
            transcript.strand,
            transcript.start
            if transcript.start is not None
            else -1,
            transcript.end
            if transcript.end is not None
            else -1,
            get_output_transcript_id(
                transcript
            ),
        ),
    )

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=COORDINATE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for transcript in sorted_transcripts:
            reference_id = (
                get_output_transcript_id(
                    transcript
                )
            )

            if reference_id in seen_reference_ids:
                raise ValueError(
                    "Duplicate reference ID in coordinate table: "
                    f"{reference_id}"
                )

            seen_reference_ids.add(
                reference_id
            )

            rows = list(
                iter_coordinate_rows(
                    transcript
                )
            )

            if not rows:
                skipped_without_exons += 1
                continue

            writer.writerows(
                rows
            )
            row_count += len(rows)
            transcript_count += 1

    print(
        "Coordinate table written:",
        output_path,
    )
    print(
        "Coordinate transcripts:",
        transcript_count,
        "segments:",
        row_count,
        "skipped_without_exons:",
        skipped_without_exons,
        "bytes:",
        output_path.stat().st_size,
    )

    if transcript_count == 0:
        raise RuntimeError(
            "The coordinate table contains no transcript segments."
        )
