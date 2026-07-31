import csv
from pathlib import Path

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
        group_tpm.get(candidate_group_id, 0.0)
        for candidate_group_id in candidate_group_ids
    )

    if total_tpm <= 0.0:
        return 1.0 / len(candidate_group_ids)

    return group_tpm.get(
        group_id,
        0.0,
    ) / total_tpm


def write_reads_info(
    reads_info,
    output_parquet,
):
    """Write read assignment information to Parquet."""

    table = pa.Table.from_pylist(
        reads_info
    )

    pq.write_table(
        table,
        output_parquet,
        compression="zstd",
    )


def write_transcript_tpm(
    transcript_tpm,
    output_tsv,
):
    """Write transcript TPM values to a TSV file."""

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
            ]
        )

        for transcript_id, tpm in transcript_tpm.items():
            writer.writerow(
                [
                    transcript_id,
                    tpm,
                ]
            )