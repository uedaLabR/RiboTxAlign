import csv
from collections import defaultdict
from pathlib import Path


def load_reference_metadata(
    metadata_tsv,
):
    """Load reference transcript definitions from metadata TSV."""

    reference_definitions = {}

    with Path(metadata_tsv).open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            reference_id = row["reference_id"]

            reference_definitions[reference_id] = {
                "reference_id": reference_id,
                "group_id": row["overlap_group_id"],
                "transcript_length": int(
                    row["transcript_length"]
                ),
                "cds_length": int(
                    row["cds_length"]
                ),
                "representative_transcript_id": row[
                    "representative_transcript_id"
                ],
                "gencode_transcript_ids": [
                    value
                    for value in row[
                        "gencode_transcript_ids"
                    ].split(",")
                    if value
                ],
            }

    return reference_definitions


def determine_transcript_lengths(
    transcript_ids,
    group_records,
    reference_definitions,
    minimum_mapped_reads=10,
):
    """Determine lengths used for transcript-level EM."""

    records_by_transcript = defaultdict(list)

    for record in group_records:
        records_by_transcript[
            record["transcript_id"]
        ].append(record)

    transcript_lengths = {}

    for transcript_id in transcript_ids:
        definition = reference_definitions[
            transcript_id
        ]

        annotated_cds_length = definition[
            "cds_length"
        ]

        if annotated_cds_length > 0:
            transcript_lengths[
                transcript_id
            ] = annotated_cds_length
            continue

        transcript_records = records_by_transcript.get(
            transcript_id,
            [],
        )

        distinct_read_count = len(
            {
                record["read_id"]
                for record in transcript_records
            }
        )

        if (
            distinct_read_count > minimum_mapped_reads
            and transcript_records
        ):
            mapped_start = min(
                int(record["reference_start"])
                for record in transcript_records
            )
            mapped_end = max(
                int(record["reference_end"])
                for record in transcript_records
            )

            transcript_lengths[
                transcript_id
            ] = max(
                mapped_end - mapped_start,
                1,
            )
        else:
            transcript_lengths[
                transcript_id
            ] = max(
                definition["transcript_length"],
                1,
            )

    return transcript_lengths