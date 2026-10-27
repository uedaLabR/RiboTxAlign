from __future__ import annotations

import csv
from pathlib import Path

from db_construct.data_class import Transcript
from db_construct.gtf_writer import (
    get_output_gene_id,
    get_output_transcript_id,
    normalize_exons,
)

METADATA_COLUMNS = [
    "reference_id",
    "overlap_group_id",
    "gene_id",
    "reference_gene_id",
    "gene_name",
    "chromosome",
    "strand",
    "genome_start",
    "genome_end",
    "exon_count",
    "transcript_length",
    "cds_start",
    "cds_end",
    "cds_reference_start",
    "cds_reference_end",
    "cds_length",
    "tpm",
    "resolution_type",
    "original_transcript_id",
    "representative_transcript_id",
    "source_reference_id",
    "gencode_transcript_ids",
    "stringtie_transcript_ids",
    "source_transcript_ids",
    "cds_equivalent_transcript_ids",
]


def join_values(
    values,
) -> str:
    """Join unique non-empty metadata values deterministically."""

    if not values:
        return ""

    return ",".join(
        sorted(
            {
                str(value)
                for value in values
                if value is not None
                and str(value) != ""
            }
        )
    )


def calculate_cds_length(
    transcript: Transcript,
) -> int:
    """Return total valid CDS length."""

    return sum(
        end - start
        for start, end in getattr(
            transcript,
            "cds_regions",
            [],
        )
        if end > start
    )


def project_cds_to_transcript_coordinates(
    transcript: Transcript,
    exons: list[tuple[int, int]],
) -> tuple[int | None, int | None]:
    """Project genomic CDS intervals onto the spliced transcript.

    Input genomic intervals and returned transcript intervals are both
    zero-based and half-open. Transcript coordinates always follow the
    biological 5-prime to 3-prime direction.
    """

    cds_regions = [
        (int(start), int(end))
        for start, end in getattr(
            transcript,
            "cds_regions",
            [],
        )
        if int(end) > int(start)
    ]

    if not exons or not cds_regions:
        return None, None

    transcript_order_exons = list(exons)

    if transcript.strand == "-":
        transcript_order_exons.reverse()

    projected_intervals: list[tuple[int, int]] = []
    transcript_offset = 0

    for exon_start, exon_end in transcript_order_exons:
        exon_length = exon_end - exon_start

        for cds_start, cds_end in cds_regions:
            overlap_start = max(exon_start, cds_start)
            overlap_end = min(exon_end, cds_end)

            if overlap_end <= overlap_start:
                continue

            if transcript.strand == "-":
                reference_start = (
                    transcript_offset
                    + exon_end
                    - overlap_end
                )
                reference_end = (
                    transcript_offset
                    + exon_end
                    - overlap_start
                )
            else:
                reference_start = (
                    transcript_offset
                    + overlap_start
                    - exon_start
                )
                reference_end = (
                    transcript_offset
                    + overlap_end
                    - exon_start
                )

            projected_intervals.append(
                (reference_start, reference_end)
            )

        transcript_offset += exon_length

    if not projected_intervals:
        return None, None

    return (
        min(start for start, _ in projected_intervals),
        max(end for _, end in projected_intervals),
    )


def build_metadata_row(
    transcript: Transcript,
) -> dict[str, object]:
    """Build one human-readable metadata row."""

    exons = normalize_exons(
        transcript
    )
    reference_id = get_output_transcript_id(
        transcript
    )

    genome_start = (
        min(start for start, _ in exons)
        if exons
        else ""
    )
    genome_end = (
        max(end for _, end in exons)
        if exons
        else ""
    )

    cds_regions = sorted(
        getattr(
            transcript,
            "cds_regions",
            [],
        )
    )
    cds_start = (
        min(start for start, _ in cds_regions)
        if cds_regions
        else ""
    )

    cds_end = (
        max(end for _, end in cds_regions)
        if cds_regions
        else ""
    )

    (
        cds_reference_start,
        cds_reference_end,
    ) = project_cds_to_transcript_coordinates(
        transcript=transcript,
        exons=exons,
    )

    return {
        "reference_id": reference_id,
        "cds_equivalent_transcript_ids": join_values(
            getattr(transcript, "cds_equivalent_transcript_ids", [])
        ),
        "overlap_group_id": getattr(
            transcript,
            "overlap_group_id",
            "",
        ) or "",
        "gene_id": get_output_gene_id(
            transcript
        ),
        "reference_gene_id": getattr(
            transcript,
            "reference_gene_id",
            "",
        ) or "",
        "gene_name": getattr(
            transcript,
            "gene_name",
            "",
        ) or "",
        "chromosome": transcript.chrom,
        "strand": transcript.strand,
        "genome_start": genome_start,
        "genome_end": genome_end,
        "exon_count": len(exons),
        "transcript_length": sum(
            end - start
            for start, end in exons
        ),
        "cds_start": cds_start,
        "cds_end": cds_end,
        "cds_reference_start": (
            cds_reference_start
            if cds_reference_start is not None
            else ""
        ),
        "cds_reference_end": (
            cds_reference_end
            if cds_reference_end is not None
            else ""
        ),
        "cds_length": calculate_cds_length(
            transcript
        ),
        "tpm": getattr(
            transcript,
            "tpm",
            0.0,
        ),
        "resolution_type": getattr(
            transcript,
            "resolution_type",
            "",
        ) or "",
        "original_transcript_id": getattr(
            transcript,
            "transcript_id",
            "",
        ) or "",
        "representative_transcript_id": getattr(
            transcript,
            "representative_transcript_id",
            "",
        ) or "",
        "source_reference_id": getattr(
            transcript,
            "source_reference_id",
            "",
        ) or "",
        "gencode_transcript_ids": join_values(
            getattr(
                transcript,
                "gencode_transcript_ids",
                [],
            )
        ),
        "stringtie_transcript_ids": join_values(
            getattr(
                transcript,
                "stringtie_transcript_ids",
                [],
            )
        ),
        "source_transcript_ids": join_values(
            getattr(
                transcript,
                "source_transcript_ids",
                [],
            )
        ),
    }


def write_reference_metadata(
    reference_transcripts: list[Transcript],
    output_path: str | Path,
) -> None:
    """Write human-readable reference construction metadata."""

    output_path = Path(
        output_path
    )
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

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
            fieldnames=METADATA_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for transcript in sorted_transcripts:
            writer.writerow(
                build_metadata_row(
                    transcript
                )
            )

    print(
        "Reference metadata written:",
        output_path,
    )
    print(
        "Metadata transcripts:",
        len(sorted_transcripts),
        "bytes:",
        output_path.stat().st_size,
    )
