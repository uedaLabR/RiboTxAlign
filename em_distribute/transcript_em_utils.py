import csv
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def load_group_tpm(
    group_em_tsv: str | Path,
) -> dict[str, float]:
    """Load final TPM values indexed by group ID."""

    results = {}

    with Path(group_em_tsv).open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:

            results[row["group_id"]] = float(
                row["final_tpm"]
            )

    return results


def calculate_group_probability(
    group_id,
    other_group_ids,
    group_tpm,
):
    """Calculate the posterior probability of the current group."""

    candidate_group_ids = [
        group_id,
        *other_group_ids,
    ]

    total_tpm = sum(
        group_tpm.get(
            candidate_group_id,
            0.0,
        )
        for candidate_group_id
        in candidate_group_ids
    )

    if total_tpm <= 0.0:

        return (
            1.0
            / len(candidate_group_ids)
        )

    return (
        group_tpm.get(
            group_id,
            0.0,
        )
        / total_tpm
    )


def process_single_transcript_group(
    group_id,
    transcript_id,
    group_records,
    group_tpm,
    reads_info,
    transcript_tpm,
):
    """Process a group containing only one transcript."""

    # This TPM is used internally.
    # Final output TPM is recalculated after all transcript
    # expected counts have been determined.
    transcript_tpm[transcript_id] = (
        group_tpm[group_id]
    )

    for record in group_records:

        other_group_ids = record.get(
            "other_group_ids",
            [],
        )

        if other_group_ids is None:
            other_group_ids = []

        prior_probability = (
            calculate_group_probability(
                group_id=group_id,
                other_group_ids=other_group_ids,
                group_tpm=group_tpm,
            )
        )

        reads_info.append(
            {
                "group_id": group_id,
                "transcript_id": transcript_id,
                "alignment_transcript_id": record.get(
                    "alignment_transcript_id",
                    record["transcript_id"],
                ),
                "read_id": record["read_id"],
                "start": record["reference_start"],
                "end": record["reference_end"],
                "read_length": int(
                    record["read_length"]
                ),
                "hit_index": int(
                    record["hit_index"]
                ),
                "query_start": int(
                    record["query_start"]
                ),
                "query_end": int(
                    record["query_end"]
                ),
                "mapping_strand": int(
                    record["strand"]
                ),
                "cigar": str(
                    record["cigar"]
                ),
                "prior_probability": (
                    prior_probability
                ),
                "transcript_probability": 1.0,
            }
        )


def write_reads_info(
    reads_info,
    output_parquet,
):
    """Write read assignment information to Parquet."""

    schema = pa.schema(
        [
            ("group_id", pa.string()),
            ("transcript_id", pa.string()),
            (
                "alignment_transcript_id",
                pa.string(),
            ),
            ("read_id", pa.string()),
            ("start", pa.int64()),
            ("end", pa.int64()),
            ("read_length", pa.int32()),
            ("hit_index", pa.int32()),
            ("query_start", pa.int32()),
            ("query_end", pa.int32()),
            ("mapping_strand", pa.int8()),
            ("cigar", pa.string()),
            (
                "prior_probability",
                pa.float64(),
            ),
            (
                "transcript_probability",
                pa.float64(),
            ),
            (
                "other_transcript_ids",
                pa.list_(
                    pa.string()
                ),
            ),
        ]
    )

    normalized_rows = []

    for row in reads_info:

        normalized_row = dict(
            row
        )

        normalized_row.setdefault(
            "other_transcript_ids",
            [],
        )

        normalized_rows.append(
            normalized_row
        )

    table = pa.Table.from_pylist(
        normalized_rows,
        schema=schema,
    )

    pq.write_table(
        table,
        output_parquet,
        compression="zstd",
    )


# ------------------------------------------------------------
# Calculate final transcript TPM
# ------------------------------------------------------------

def calculate_final_transcript_tpm(
    transcript_expected_count,
    transcript_em_cds_length,
):
    """
    Calculate final transcript TPM from expected read counts.

    EM-internal TPM values are not used here.
    Final TPM is calculated globally across all transcripts.
    """

    transcript_rpk = {}

    for (
        transcript_id,
        expected_count,
    ) in transcript_expected_count.items():

        cds_length = (
            transcript_em_cds_length.get(
                transcript_id,
                0,
            )
        )

        if cds_length <= 0:

            transcript_rpk[
                transcript_id
            ] = 0.0

            continue

        transcript_rpk[
            transcript_id
        ] = (
            float(expected_count)
            / float(cds_length)
        )

    total_rpk = sum(
        transcript_rpk.values()
    )

    if total_rpk <= 0.0:

        return {
            transcript_id: 0.0
            for transcript_id
            in transcript_rpk
        }

    transcript_tpm = {}

    for (
        transcript_id,
        rpk,
    ) in transcript_rpk.items():

        transcript_tpm[
            transcript_id
        ] = (
            rpk
            / total_rpk
            * 1_000_000.0
        )

    return transcript_tpm


# ------------------------------------------------------------
# Write detailed RiboTxAlign transcript results
# ------------------------------------------------------------

def write_transcript_tpm(
    transcript_tpm,
    transcript_expected_count,
    transcript_original_cds_length,
    transcript_em_cds_length,
    output_tsv,
):
    """Write transcript-level results using RiboTxAlign IDs."""

    with Path(output_tsv).open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "transcript_id",
                "tpm",
                "expected_count",
                "original_cds_length",
                "em_cds_length",
            ]
        )

        for transcript_id in sorted(
            transcript_tpm
        ):

            writer.writerow(
                [
                    transcript_id,
                    transcript_tpm[
                        transcript_id
                    ],
                    transcript_expected_count.get(
                        transcript_id,
                        0.0,
                    ),
                    transcript_original_cds_length.get(
                        transcript_id,
                        0,
                    ),
                    transcript_em_cds_length.get(
                        transcript_id,
                        0,
                    ),
                ]
            )


# ------------------------------------------------------------
# Determine source of representative transcript ID
# ------------------------------------------------------------

def determine_transcript_source(
    transcript_id,
):
    """Infer transcript source from transcript ID."""

    transcript_id = str(
        transcript_id
    )

    if transcript_id.startswith(
        "ENST"
    ):
        return "GENCODE"

    if (
        transcript_id.startswith(
            "MSTRG"
        )
        or transcript_id.startswith(
            "STRG"
        )
    ):
        return "StringTie"

    return "other"


# ------------------------------------------------------------
# Write representative transcript results
# ------------------------------------------------------------

def write_representative_transcript_tpm(
    transcript_tpm,
    transcript_expected_count,
    reference_metadata_tsv,
    output_tsv,
):
    """
    Write user-facing transcript abundance table using the
    representative GENCODE/StringTie transcript ID.
    """

    metadata_df = pd.read_csv(
        reference_metadata_tsv,
        sep="\t",
    )

    required_columns = {
        "reference_id",
        "representative_transcript_id",
    }

    missing_columns = (
        required_columns
        - set(metadata_df.columns)
    )

    if missing_columns:

        raise ValueError(
            "Missing metadata columns: "
            + ", ".join(
                sorted(
                    missing_columns
                )
            )
        )

    metadata_map = {}

    for _, row in metadata_df.iterrows():

        reference_id = str(
            row["reference_id"]
        )

        representative_id = (
            row[
                "representative_transcript_id"
            ]
        )

        if pd.isna(
            representative_id
        ):
            continue

        metadata_map[
            reference_id
        ] = str(
            representative_id
        )

    with Path(output_tsv).open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "transcript_id",
                "source",
                "tpm",
                "expected_count",
                "reference_id",
            ]
        )

        for reference_id in sorted(
            transcript_tpm
        ):

            representative_id = (
                metadata_map.get(
                    reference_id
                )
            )

            if representative_id is None:
                continue

            source = (
                determine_transcript_source(
                    representative_id
                )
            )

            writer.writerow(
                [
                    representative_id,
                    source,
                    transcript_tpm[
                        reference_id
                    ],
                    transcript_expected_count.get(
                        reference_id,
                        0.0,
                    ),
                    reference_id,
                ]
            )


def append_transcript_merge_log(
    group_id,
    merge_map,
    output_tsv,
):
    """Append transcript merge relationships to a TSV file."""

    output_path = Path(
        output_tsv
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    write_header = (
        not output_path.exists()
        or output_path.stat().st_size == 0
    )

    with output_path.open(
        "a",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        if write_header:

            writer.writerow(
                [
                    "group_id",
                    "source_transcript_id",
                    "target_transcript_id",
                ]
            )

        for (
            source_id,
            target_id,
        ) in merge_map.items():

            writer.writerow(
                [
                    group_id,
                    source_id,
                    target_id,
                ]
            )

        handle.flush()

def write_transcript_tpm_with_source(
    transcript_tpm,
    transcript_expected_count,
    reference_metadata_tsv,
    output_tsv,
):
    """
    Write transcript abundance while preserving RiboTxAlign
    reference IDs and source transcript information.
    """

    metadata = {}

    with open(
        reference_metadata_tsv,
        newline="",
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t",
        )

        for row in reader:

            reference_id = row.get(
                "reference_id",
                ""
            )

            if not reference_id:
                continue

            metadata[
                reference_id
            ] = row

    output_columns = [
        "reference_id",
        "overlap_group_id",
        "gene_id",
        "reference_gene_id",
        "gene_name",
        "tpm",
        "expected_count",
        "resolution_type",
        "original_transcript_id",
        "representative_transcript_id",
        "source_reference_id",
        "gencode_transcript_ids",
        "stringtie_transcript_ids",
        "source_transcript_ids",
    ]

    with open(
        output_tsv,
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=output_columns,
            delimiter="\t",
        )

        writer.writeheader()

        for reference_id in sorted(
            transcript_tpm
        ):

            row = metadata.get(
                reference_id,
                {}
            )

            writer.writerow(
                {
                    "reference_id": reference_id,

                    "overlap_group_id": row.get(
                        "overlap_group_id",
                        "",
                    ),

                    "gene_id": row.get(
                        "gene_id",
                        "",
                    ),

                    "reference_gene_id": row.get(
                        "reference_gene_id",
                        "",
                    ),

                    "gene_name": row.get(
                        "gene_name",
                        "",
                    ),

                    "tpm": transcript_tpm.get(
                        reference_id,
                        0.0,
                    ),

                    "expected_count": (
                        transcript_expected_count.get(
                            reference_id,
                            0.0,
                        )
                    ),

                    "resolution_type": row.get(
                        "resolution_type",
                        "",
                    ),

                    "original_transcript_id": row.get(
                        "original_transcript_id",
                        "",
                    ),

                    "representative_transcript_id": row.get(
                        "representative_transcript_id",
                        "",
                    ),

                    "source_reference_id": row.get(
                        "source_reference_id",
                        "",
                    ),

                    "gencode_transcript_ids": row.get(
                        "gencode_transcript_ids",
                        "",
                    ),

                    "stringtie_transcript_ids": row.get(
                        "stringtie_transcript_ids",
                        "",
                    ),

                    "source_transcript_ids": row.get(
                        "source_transcript_ids",
                        "",
                    ),
                }
            )