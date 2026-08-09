import csv
from pathlib import Path
from em_distribute.transcript_em_utils import *



from pathlib import Path

import pyarrow.parquet as pq


def iter_sorted_mapping_groups(
    sorted_mapping_dir: str | Path,
):
    """Yield consecutive records grouped by group_id."""

    current_group_id = None
    current_group_records = []

    parquet_files = sorted(
        Path(sorted_mapping_dir).glob("*.parquet")
    )

    for parquet_path in parquet_files:
        parquet_file = pq.ParquetFile(
            parquet_path
        )

        for batch in parquet_file.iter_batches():
            rows = batch.to_pylist()

            for row in rows:
                group_id = row["group_id"]

                if current_group_id is None:
                    current_group_id = group_id

                if group_id != current_group_id:
                    yield (
                        current_group_id,
                        current_group_records,
                    )

                    current_group_id = group_id
                    current_group_records = []

                current_group_records.append(
                    row
                )

    if current_group_records:
        yield (
            current_group_id,
            current_group_records,
        )

import csv
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq







from em_distribute.transcript_group_multi import (
    process_multi_transcript_group,
)

from em_distribute.reference_utils import *
def em_transcript(
    group_em_tsv: str | Path,
    sorted_mapping_dir: str | Path,
    reads_info_parquet: str | Path,
    transcript_tpm_tsv: str | Path,
    transcript_merge_log_tsv: str | Path,
    reference_metadata_tsv: str | Path,
    reference_gtf: str | Path,

):
    group_tpm = load_group_tpm(
        group_em_tsv
    )
    reference_definitions = load_reference_definitions(
        metadata_tsv=reference_metadata_tsv,
        reference_gtf=reference_gtf,
    )

    reads_info = []
    transcript_tpm = {}

    for group_id, group_records in iter_sorted_mapping_groups(
        sorted_mapping_dir
    ):
        tpm = group_tpm.get(
            group_id,
            0.0,
        )

        if tpm <= 0.0:
            continue

        transcript_ids = {
            record["transcript_id"]
            for record in group_records
        }

        print(
            "group_id,tpm,num_transcripts,num_records",
            group_id,
            tpm,
            len(transcript_ids),
            len(group_records),
        )

        if len(transcript_ids) == 1:
            transcript_id = next(
                iter(transcript_ids)
            )

            process_single_transcript_group(
                group_id=group_id,
                transcript_id=transcript_id,
                group_records=group_records,
                group_tpm=group_tpm,
                reads_info=reads_info,
                transcript_tpm=transcript_tpm,
            )

        elif len(transcript_ids) > 1:
            # Multi-transcript processing will be added next.
            process_multi_transcript_group(
                group_id=group_id,
                transcript_ids=transcript_ids,
                group_records=group_records,
                group_tpm=group_tpm,
                reads_info=reads_info,
                transcript_tpm=transcript_tpm,
                transcript_merge_log_tsv=transcript_merge_log_tsv,
                reference_definitions=reference_definitions
            )
        else:
            print("error not transcript found in group -- should not happen")

    write_reads_info(
        reads_info=reads_info,
        output_parquet=reads_info_parquet,
    )

    write_transcript_tpm(
        transcript_tpm=transcript_tpm,
        output_tsv=transcript_tpm_tsv,
    )

