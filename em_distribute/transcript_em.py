import csv
import logging
from pathlib import Path
import time

import pyarrow.parquet as pq

from em_distribute.transcript_em_utils import *
from em_distribute.transcript_group_multi import (
    process_multi_transcript_group,
)
from em_distribute.reference_utils import *
from em_distribute.cds_candidate_filter import build_group_cds_spans, filter_cds_candidates


def create_transcript_em_logger(log_path: str | Path) -> logging.Logger:
    """Create the transcript-EM file and console logger."""

    path = Path(log_path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"RiboTxAlignTranscriptEM:{path}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s")
    file_handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


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



def add_reference_rna_tpm(
    reference_definitions,
    reference_metadata_tsv: str | Path,
):
    """Add RNA-seq TPM from reference metadata to reference definitions."""

    with open(
        reference_metadata_tsv,
        encoding="utf-8",
    ) as input_file:
        reader = csv.DictReader(
            input_file,
            delimiter="\t",
        )

        for row in reader:
            reference_id = (
                row.get("reference_id")
                or ""
            ).strip()

            if (
                not reference_id
                or reference_id not in reference_definitions
            ):
                continue

            try:
                rna_tpm = float(
                    row.get("tpm")
                    or 0.0
                )
            except (TypeError, ValueError):
                rna_tpm = 0.0

            reference_definitions[
                reference_id
            ]["rna_tpm"] = max(
                rna_tpm,
                0.0,
            )

    return reference_definitions


def em_transcript(
    group_em_tsv: str | Path,
    sorted_mapping_dir: str | Path,
    reads_info_parquet: str | Path,
    transcript_tpm_tsv: str | Path,
    transcript_merge_log_tsv: str | Path,
    reference_metadata_tsv: str | Path,
    reference_gtf: str | Path,
    output_log_path: str | Path | None = None,
):

    workflow_started_at = time.perf_counter()
    if output_log_path is None:
        output_log_path = Path(transcript_tpm_tsv).with_suffix(".transcript_em.log")
    logger = create_transcript_em_logger(output_log_path)
    logger.info("Transcript EM started")

    loading_started_at = time.perf_counter()

    group_tpm = load_group_tpm(
        group_em_tsv
    )

    reference_definitions = (
        load_reference_definitions(
            metadata_tsv=reference_metadata_tsv,
            reference_gtf=reference_gtf,
        )
    )

    reference_definitions = add_reference_rna_tpm(
        reference_definitions=reference_definitions,
        reference_metadata_tsv=reference_metadata_tsv,
    )
    logger.info(
        "Transcript EM input loading elapsed: %.3f s",
        time.perf_counter() - loading_started_at,
    )

    # Start one log per execution; individual groups still append to it.
    # Write a header even when this run produces no transcript merges.
    merge_log_path = Path(transcript_merge_log_tsv)
    merge_log_path.parent.mkdir(parents=True, exist_ok=True)
    with merge_log_path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerow([
            "group_id", "source_transcript_id", "target_transcript_id",
        ])

    group_cds_spans = build_group_cds_spans(reference_definitions)
    cds_region_rows = []
    cds_filter_rows = []

    # This dictionary is used internally by transcript EM.
    transcript_tpm = {}

    transcript_expected_count = {}

    transcript_original_cds_length = {}

    transcript_em_cds_length = {}

    group_processing_started_at = time.perf_counter()
    processed_group_count = 0
    positive_group_count = 0
    with open_reads_info_writer(reads_info_parquet) as reads_info:
        for (
            group_id,
            group_records,
        ) in iter_sorted_mapping_groups(
            sorted_mapping_dir
        ):

            processed_group_count += 1

            tpm = group_tpm.get(
                group_id,
                0.0,
            )

            if tpm <= 0.0:
                continue
            positive_group_count += 1

            original_records = group_records
            group_records, region_rows, lost_reads = filter_cds_candidates(
                original_records, reference_definitions, group_cds_spans
            )
            cds_region_rows.extend(dict(group_id=group_id, **row) for row in region_rows)
            read_priors = {}
            for record in original_records:
                read_id = record["read_id"]
                if read_id not in read_priors:
                    read_priors[read_id] = calculate_group_probability(
                        group_id=group_id,
                        other_group_ids=record.get("other_group_ids") or [],
                        group_tpm=group_tpm,
                    )
            cds_filter_rows.append(dict(
                group_id=group_id, input_reads=len(read_priors),
                retained_reads=len(read_priors) - len(lost_reads),
                unassigned_reads=len(lost_reads),
                unassigned_expected_count=sum(read_priors[r] for r in lost_reads),
                removed_hit_records=len(original_records) - len(group_records),
            ))
            if not group_records:
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

            # ----------------------------------------------------
            # Single-transcript group
            # ----------------------------------------------------

            if len(transcript_ids) == 1:

                transcript_id = next(
                    iter(
                        transcript_ids
                    )
                )

                process_single_transcript_group(
                    group_id=group_id,
                    transcript_id=transcript_id,
                    group_records=group_records,
                    group_tpm=group_tpm,
                    reads_info=reads_info,
                    transcript_tpm=transcript_tpm,
                )

                # ------------------------------------------------
                # Calculate expected read count
                # ------------------------------------------------

                expected_count = 0.0

                seen_reads = set()

                for record in group_records:

                    read_id = record[
                        "read_id"
                    ]

                    if read_id in seen_reads:
                        continue

                    seen_reads.add(
                        read_id
                    )

                    other_group_ids = (
                        record.get(
                            "other_group_ids",
                            [],
                        )
                    )

                    if other_group_ids is None:
                        other_group_ids = []

                    expected_count += (
                        calculate_group_probability(
                            group_id=group_id,
                            other_group_ids=(
                                other_group_ids
                            ),
                            group_tpm=group_tpm,
                        )
                    )

                transcript_expected_count[
                    transcript_id
                ] = expected_count

                # ------------------------------------------------
                # CDS length
                # ------------------------------------------------

                original_cds_length = int(
                    reference_definitions[
                        transcript_id
                    ].get(
                        "cds_length",
                        0,
                    )
                )

                transcript_length = int(
                    reference_definitions[
                        transcript_id
                    ].get(
                        "transcript_length",
                        1,
                    )
                )

                em_cds_length = (
                    reference_definitions[transcript_id]["em_region_length"]
                )

                transcript_original_cds_length[
                    transcript_id
                ] = original_cds_length

                transcript_em_cds_length[
                    transcript_id
                ] = em_cds_length

            # ----------------------------------------------------
            # Multi-transcript group
            # ----------------------------------------------------

            elif len(transcript_ids) > 1:

                result = (
                    process_multi_transcript_group(
                        group_id=group_id,
                        transcript_ids=transcript_ids,
                        group_records=group_records,
                        group_tpm=group_tpm,
                        reads_info=reads_info,
                        transcript_tpm=transcript_tpm,
                        transcript_merge_log_tsv=(
                            transcript_merge_log_tsv
                        ),
                        reference_definitions=(
                            reference_definitions
                        ),
                    )
                )

                if result is None:
                    continue

                (
                    final_expected_counts,
                    original_cds_lengths,
                    em_cds_lengths,
                ) = result

                transcript_expected_count.update(
                    final_expected_counts
                )

                transcript_original_cds_length.update(
                    original_cds_lengths
                )

                transcript_em_cds_length.update(
                    em_cds_lengths
                )

            else:

                print(
                    "ERROR: no transcript found "
                    "in group",
                    group_id,
                )

            if processed_group_count % 1000 == 0:
                logger.info(
                    "Transcript EM progress: groups=%d positive_groups=%d reads_info_rows=%d elapsed=%.3fs",
                    processed_group_count,
                    positive_group_count,
                    reads_info.row_count,
                    time.perf_counter() - group_processing_started_at,
                )

    logger.info(
        "Transcript EM group processing elapsed: %.3f s (groups=%d positive_groups=%d)",
        time.perf_counter() - group_processing_started_at,
        processed_group_count,
        positive_group_count,
    )

    # --------------------------------------------------------
    # Write read-level information
    # --------------------------------------------------------

    output_started_at = time.perf_counter()
    for suffix, rows, fields in (
        ("cds_regions", cds_region_rows,
         ["group_id", "transcript_id", "region_source", "region_intervals", "em_region_length"]),
        ("cds_filter", cds_filter_rows,
         ["group_id", "input_reads", "retained_reads", "unassigned_reads",
          "unassigned_expected_count", "removed_hit_records"]),
    ):
        diagnostic_path = Path(transcript_tpm_tsv).with_suffix(f".{suffix}.tsv")
        with diagnostic_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    logger.info(
        "Transcript EM read output rows: %d (buffer limit=%d)",
        reads_info.row_count, reads_info.batch_size,
    )

    # --------------------------------------------------------
    # IMPORTANT:
    # Recalculate final TPM globally from expected counts.
    #
    # Do not use EM-internal transcript_tpm as final output.
    # --------------------------------------------------------

    final_transcript_tpm = (
        calculate_final_transcript_tpm(
            transcript_expected_count=(
                transcript_expected_count
            ),
            transcript_em_cds_length=(
                transcript_em_cds_length
            ),
        )
    )

    # --------------------------------------------------------
    # Detailed output using RiboTxAlign reference IDs
    # --------------------------------------------------------

    write_transcript_tpm(
        transcript_tpm=final_transcript_tpm,
        transcript_expected_count=(
            transcript_expected_count
        ),
        transcript_original_cds_length=(
            transcript_original_cds_length
        ),
        transcript_em_cds_length=(
            transcript_em_cds_length
        ),
        output_tsv=transcript_tpm_tsv,
    )

    # --------------------------------------------------------
    # User-facing output using representative transcript IDs
    # --------------------------------------------------------

    transcript_tpm_path = Path(
        transcript_tpm_tsv
    )

    transcript_tpm_path = Path(
        transcript_tpm_tsv
    )

    source_tpm_tsv = (
            transcript_tpm_path.parent
            / (
                    transcript_tpm_path.stem
                    + "_with_source.tsv"
            )
    )

    write_transcript_tpm_with_source(
        transcript_tpm=(
            final_transcript_tpm
        ),
        transcript_expected_count=(
            transcript_expected_count
        ),
        reference_metadata_tsv=(
            reference_metadata_tsv
        ),
        output_tsv=(
            source_tpm_tsv
        ),
    )
    logger.info(
        "Transcript EM output writing elapsed: %.3f s",
        time.perf_counter() - output_started_at,
    )
    logger.info(
        "Transcript EM total elapsed: %.3f s",
        time.perf_counter() - workflow_started_at,
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 70)
    print("Transcript abundance output")
    print("=" * 70)

    print(
        "Transcripts:",
        len(final_transcript_tpm),
    )

    print(
        "Expected reads:",
        sum(
            transcript_expected_count.values()
        ),
    )

    print(
        "TPM sum:",
        sum(
            final_transcript_tpm.values()
        ),
    )

    print()

    print(
        "Detailed output:"
    )

    print(
        transcript_tpm_tsv
    )

    print()

    print(
        "Transcript source output:"
    )

    print(
        source_tpm_tsv
    )
