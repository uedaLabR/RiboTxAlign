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


def process_single_transcript_group(
    group_id,
    transcript_id,
    group_records,
    group_tpm,
    reads_info,
    transcript_tpm,
):
    """Process a group containing only one transcript."""

    transcript_tpm[transcript_id] = group_tpm[group_id]

    for record in group_records:
        other_group_ids = record.get(
            "other_group_ids",
            [],
        )

        if other_group_ids is None:
            other_group_ids = []

        prior_probability = calculate_group_probability(
            group_id=group_id,
            other_group_ids=other_group_ids,
            group_tpm=group_tpm,
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
                "hit_index": int(record["hit_index"]),
                "query_start": int(record["query_start"]),
                "query_end": int(record["query_end"]),
                "mapping_strand": int(record["strand"]),
                "cigar": str(record["cigar"]),
                "prior_probability": prior_probability,
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
            ("alignment_transcript_id", pa.string()),
            ("read_id", pa.string()),
            ("start", pa.int64()),
            ("end", pa.int64()),
            ("read_length", pa.int32()),
            ("hit_index", pa.int32()),
            ("query_start", pa.int32()),
            ("query_end", pa.int32()),
            ("mapping_strand", pa.int8()),
            ("cigar", pa.string()),
            ("prior_probability", pa.float64()),
            ("transcript_probability", pa.float64()),
            (
                "other_transcript_ids",
                pa.list_(pa.string()),
            ),
        ]
    )

    normalized_rows = []

    for row in reads_info:
        normalized_row = dict(row)
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

import csv
from pathlib import Path


def append_transcript_merge_log(
    group_id,
    merge_map,
    output_tsv,
):
    """Append transcript merge relationships to a TSV file."""

    output_path = Path(output_tsv)
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

        for source_id, target_id in merge_map.items():
            writer.writerow(
                [
                    group_id,
                    source_id,
                    target_id,
                ]
            )

        handle.flush()
