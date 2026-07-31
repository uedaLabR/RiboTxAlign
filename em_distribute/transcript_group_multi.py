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

def process_multi_transcript_group(
    group_id,
    transcript_ids,
    group_records,
    group_tpm,
    reads_info,
    transcript_tpm,
):
    """Process a group containing multiple transcripts."""

    (
        unique_read_indexes,
        multi_read_indexes,
        index_to_read_id,
    ) = classify_reads_by_transcript(
        group_records
    )

    for transcript_id in transcript_ids:
        transcript_unique_reads = (
            unique_read_indexes.get(
                transcript_id,
                set(),
            )
        )

        transcript_multi_reads = (
            multi_read_indexes.get(
                transcript_id,
                set(),
            )
        )

        print(
            group_id,
            transcript_id,
            "unique",
            len(transcript_unique_reads),
            "multi",
            len(transcript_multi_reads),
        )

    #find transcript only supoorted by multimap
    merge_map = find_subset_merges(
        transcript_ids,
        unique_read_indexes,
        multi_read_indexes,
    )
    print("merge_map",merge_map)

    return (
        unique_read_indexes,
        multi_read_indexes,
        index_to_read_id,
    )