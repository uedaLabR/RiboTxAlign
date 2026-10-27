from __future__ import annotations

from pathlib import Path

from db_construct.data_class import Transcript


def escape_gtf_value(
    value: object,
) -> str:
    """Escape a value for use in a GTF attribute."""

    text = str(value)

    return (
        text
        .replace("\\", "\\\\")
        .replace('"', '\\"')
    )


def get_output_transcript_id(
    transcript: Transcript,
) -> str:
    """Return the final transcript ID used in the output GTF."""

    reference_id = getattr(
        transcript,
        "reference_id",
        None,
    )

    if reference_id:
        return str(reference_id)

    transcript_id = getattr(
        transcript,
        "transcript_id",
        None,
    )

    if transcript_id:
        return str(transcript_id)

    raise ValueError(
        "Transcript has neither reference_id nor transcript_id."
    )


def get_output_gene_id(
    transcript: Transcript,
) -> str:
    """Return a stable gene ID for the output GTF."""

    overlap_group_id = getattr(
        transcript,
        "overlap_group_id",
        None,
    )

    if overlap_group_id:
        return str(overlap_group_id)

    reference_gene_id = getattr(
        transcript,
        "reference_gene_id",
        None,
    )

    if reference_gene_id:
        return str(reference_gene_id)

    gene_id = getattr(
        transcript,
        "gene_id",
        None,
    )

    if gene_id:
        return str(gene_id)

    return get_output_transcript_id(
        transcript
    )


def get_gene_name(
    transcript: Transcript,
) -> str | None:
    """Return the gene name when available."""

    gene_name = getattr(
        transcript,
        "gene_name",
        None,
    )

    if gene_name:
        return str(gene_name)

    return None


def normalize_exons(
    transcript: Transcript,
) -> list[tuple[int, int]]:
    """Return valid, unique exons in genomic coordinate order."""

    normalized: list[
        tuple[int, int]
    ] = []

    for exon in getattr(
        transcript,
        "exons",
        [],
    ):
        if len(exon) != 2:
            continue

        start, end = exon

        start = int(start)
        end = int(end)

        if end <= start:
            continue

        normalized.append(
            (
                start,
                end,
            )
        )

    return sorted(
        set(normalized)
    )


def build_gtf_attributes(
    transcript: Transcript,
    transcript_id: str,
    gene_id: str,
    exon_number: int | None = None,
) -> str:
    """Build GTF attributes for a transcript or exon record."""

    attributes = [
        (
            "gene_id",
            gene_id,
        ),
        (
            "transcript_id",
            transcript_id,
        ),
    ]

    gene_name = get_gene_name(
        transcript
    )

    if gene_name:
        attributes.append(
            (
                "gene_name",
                gene_name,
            )
        )

    original_transcript_id = getattr(
        transcript,
        "transcript_id",
        None,
    )

    if (
        original_transcript_id
        and str(original_transcript_id)
        != transcript_id
    ):
        attributes.append(
            (
                "original_transcript_id",
                original_transcript_id,
            )
        )

    source_ids = getattr(
        transcript,
        "source_transcript_ids",
        None,
    )

    if source_ids:
        attributes.append(
            (
                "source_transcript_ids",
                ",".join(
                    str(source_id)
                    for source_id in source_ids
                ),
            )
        )

    cds_aliases = getattr(transcript, "cds_equivalent_transcript_ids", [])
    if cds_aliases:
        attributes.append(("cds_equivalent_transcript_ids", ",".join(cds_aliases)))

    resolution_type = getattr(
        transcript,
        "resolution_type",
        None,
    )

    if resolution_type:
        attributes.append(
            (
                "resolution_type",
                resolution_type,
            )
        )

    if exon_number is not None:
        attributes.append(
            (
                "exon_number",
                exon_number,
            )
        )

    return " ".join(
        f'{key} "{escape_gtf_value(value)}";'
        for key, value in attributes
    )


def write_reference_gtf(
    reference_transcripts: list[Transcript],
    output_path: str | Path,
) -> None:
    """
    Write final reference transcripts as a genome-coordinate GTF.

    Transcript and exon coordinates are stored internally as 0-based,
    half-open intervals and converted here to 1-based, closed GTF
    coordinates.
    """

    output_path = Path(
        output_path
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    transcript_count = 0
    exon_count = 0
    skipped_without_exons = 0

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
    ) as output_file:

        for transcript in sorted_transcripts:
            exons = normalize_exons(
                transcript
            )

            if not exons:
                skipped_without_exons += 1
                continue

            transcript_id = (
                get_output_transcript_id(
                    transcript
                )
            )

            gene_id = get_output_gene_id(
                transcript
            )

            transcript_start = min(
                start
                for start, _ in exons
            )

            transcript_end = max(
                end
                for _, end in exons
            )

            transcript_attributes = (
                build_gtf_attributes(
                    transcript=transcript,
                    transcript_id=transcript_id,
                    gene_id=gene_id,
                )
            )

            output_file.write(
                "\t".join(
                    [
                        str(transcript.chrom),
                        "RiboTxAlign",
                        "transcript",
                        str(
                            transcript_start + 1
                        ),
                        str(transcript_end),
                        ".",
                        str(transcript.strand),
                        ".",
                        transcript_attributes,
                    ]
                )
                + "\n"
            )

            transcript_count += 1

            exon_order = exons

            if transcript.strand == "-":
                exon_order = list(
                    reversed(exons)
                )

            exon_number_by_coordinate = {
                exon: exon_number
                for exon_number, exon in enumerate(
                    exon_order,
                    start=1,
                )
            }

            for exon in exons:
                exon_start, exon_end = exon

                exon_attributes = (
                    build_gtf_attributes(
                        transcript=transcript,
                        transcript_id=(
                            transcript_id
                        ),
                        gene_id=gene_id,
                        exon_number=(
                            exon_number_by_coordinate[
                                exon
                            ]
                        ),
                    )
                )

                output_file.write(
                    "\t".join(
                        [
                            str(
                                transcript.chrom
                            ),
                            "RiboTxAlign",
                            "exon",
                            str(exon_start + 1),
                            str(exon_end),
                            ".",
                            str(
                                transcript.strand
                            ),
                            ".",
                            exon_attributes,
                        ]
                    )
                    + "\n"
                )

                exon_count += 1

    file_size = output_path.stat().st_size

    print(
        "Reference GTF written:",
        output_path,
    )
    print(
        "GTF transcripts:",
        transcript_count,
        "exons:",
        exon_count,
        "skipped_without_exons:",
        skipped_without_exons,
        "bytes:",
        file_size,
    )

    if transcript_count == 0:
        raise RuntimeError(
            "The output GTF contains no transcript records. "
            "Check that reference transcripts contain exon coordinates."
        )
