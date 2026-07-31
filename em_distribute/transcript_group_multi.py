from collections import defaultdict


def classify_reads_by_transcript(
    group_records,
):
    """Classify integer read indexes as transcript-unique or multi-hit."""

    read_id_to_index = {}
    index_to_read_id = []
    read_to_transcripts = defaultdict(set)

    for record in group_records:
        read_id = record["read_id"]

        read_index = read_id_to_index.get(
            read_id
        )

        if read_index is None:
            read_index = len(
                index_to_read_id
            )

            read_id_to_index[
                read_id
            ] = read_index

            index_to_read_id.append(
                read_id
            )

        read_to_transcripts[
            read_index
        ].add(
            record["transcript_id"]
        )

    unique_read_indexes = defaultdict(set)
    multi_read_indexes = defaultdict(set)

    for read_index, candidate_transcripts in (
        read_to_transcripts.items()
    ):
        if len(candidate_transcripts) == 1:
            transcript_id = next(
                iter(candidate_transcripts)
            )

            unique_read_indexes[
                transcript_id
            ].add(
                read_index
            )

        else:
            for transcript_id in candidate_transcripts:
                multi_read_indexes[
                    transcript_id
                ].add(
                    read_index
                )

    return (
        unique_read_indexes,
        multi_read_indexes,
        index_to_read_id,
    )


def resolve_merge_map(
    merge_map,
):
    """Resolve chained transcript merges."""

    resolved_map = {}

    for transcript_id in merge_map:
        target_id = merge_map[transcript_id]

        while target_id in merge_map:
            target_id = merge_map[target_id]

        resolved_map[transcript_id] = target_id

    return resolved_map

def find_subset_merges(
    transcript_ids,
    unique_read_indexes,
    multi_read_indexes,
):
    merge_map = {}

    for transcript_id in transcript_ids:
        if unique_read_indexes.get(
            transcript_id,
            set(),
        ):
            continue

        source_reads = multi_read_indexes.get(
            transcript_id,
            set(),
        )

        for other_id in transcript_ids:
            if other_id == transcript_id:
                continue

            other_reads = multi_read_indexes.get(
                other_id,
                set(),
            )

            if source_reads < other_reads:
                merge_map[transcript_id] = other_id
                break

    if len(merge_map) > 1:
        merge_map = resolve_merge_map(merge_map)

    return merge_map


def merge_transcripts(
    group_id,
    transcript_ids,
    group_records,
    merge_map,
):
    """Apply transcript merges and remove duplicate read-transcript hits."""

    before_transcript_count = len(
        transcript_ids
    )

    merged_group_records = apply_transcript_merges(
        group_records=group_records,
        merge_map=merge_map,
    )

    merged_transcript_ids = {
        record["transcript_id"]
        for record in merged_group_records
    }

    after_transcript_count = len(
        merged_transcript_ids
    )

    print(
        "transcript merge",
        group_id,
        "before",
        before_transcript_count,
        "after",
        after_transcript_count,
        "reduced",
        before_transcript_count - after_transcript_count,
    )

    return (
        merged_transcript_ids,
        merged_group_records,
    )


def apply_transcript_merges(
    group_records,
    merge_map,
):
    """Apply transcript merges and remove duplicate read-transcript hits."""

    merged_records = []
    seen_hits = set()

    for record in group_records:
        transcript_id = merge_map.get(
            record["transcript_id"],
            record["transcript_id"],
        )

        hit_key = (
            record["read_id"],
            transcript_id,
        )

        if hit_key in seen_hits:
            continue

        seen_hits.add(hit_key)

        updated_record = record.copy()
        updated_record["transcript_id"] = transcript_id

        merged_records.append(
            updated_record
        )

    return merged_records

from em_distribute.transcript_em_utils import *
from em_distribute.reference_utils import *
from em_distribute.transcript_em_eachgroup import (
    execute_transcript_em_eachgroup,
)

def process_multi_transcript_group(
    group_id,
    transcript_ids,
    group_records,
    group_tpm,
    reads_info,
    transcript_tpm,
    transcript_merge_log_tsv,
    reference_definitions,
):
    """Process a group containing multiple transcripts."""

    (
        unique_read_indexes,
        multi_read_indexes,
        index_to_read_id,
    ) = classify_reads_by_transcript(
        group_records
    )

    merge_map = find_subset_merges(
        transcript_ids=transcript_ids,
        unique_read_indexes=unique_read_indexes,
        multi_read_indexes=multi_read_indexes,
    )

    if merge_map:

        append_transcript_merge_log(
            group_id=group_id,
            merge_map=merge_map,
            output_tsv=transcript_merge_log_tsv,
        )

        (
            transcript_ids,
            group_records,
        ) = merge_transcripts(
            group_id=group_id,
            transcript_ids=transcript_ids,
            group_records=group_records,
            merge_map=merge_map,
        )

        # The group became a single-transcript group after merging.
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
            return

    transcript_cds_lengths = determine_transcript_lengths(
        transcript_ids=transcript_ids,
        group_records=group_records,
        reference_definitions=reference_definitions,
        minimum_mapped_reads=10,
    )

    # Transcript-level EM will start from here.
    execute_transcript_em_eachgroup(
        group_id=group_id,
        transcript_ids=transcript_ids,
        group_records=group_records,
        group_tpm=group_tpm,
        transcript_cds_lengths=transcript_cds_lengths,
        reads_info=reads_info,
        transcript_tpm=transcript_tpm,
    )