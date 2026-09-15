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

def has_unique_evidence(
    transcript_id,
    unique_read_indexes,
):
    """Return True when a transcript has at least one observed unique read."""

    return bool(
        unique_read_indexes.get(
            transcript_id,
            set(),
        )
    )


def get_multi_read_evidence(
    transcript_id,
    multi_read_indexes,
):
    """Return the observed multi-read index set for one transcript."""

    return multi_read_indexes.get(
        transcript_id,
        set(),
    )


def have_identical_multi_read_evidence(
    transcript_id_1,
    transcript_id_2,
    multi_read_indexes,
):
    """Return True when two transcripts have the same non-empty multi-read set."""

    reads_1 = get_multi_read_evidence(
        transcript_id_1,
        multi_read_indexes,
    )

    reads_2 = get_multi_read_evidence(
        transcript_id_2,
        multi_read_indexes,
    )

    return bool(
        reads_1
        and reads_1 == reads_2
    )


def get_reference_rna_tpm(
    transcript_id,
    reference_definitions,
):
    """Return non-negative RNA-seq TPM stored in reference metadata."""

    try:
        rna_tpm = float(
            reference_definitions[
                transcript_id
            ].get(
                "rna_tpm",
                0.0,
            )
        )
    except (
        KeyError,
        TypeError,
        ValueError,
    ):
        return 0.0

    return max(
        rna_tpm,
        0.0,
    )


def select_merge_representative(
    transcript_ids,
    reference_definitions,
):
    """Select the highest-RNA-TPM transcript as merge representative."""

    def sort_key(
        transcript_id,
    ):
        definition = reference_definitions.get(
            transcript_id,
            {},
        )

        rna_tpm = get_reference_rna_tpm(
            transcript_id,
            reference_definitions,
        )

        cds_length = int(
            definition.get(
                "cds_length",
                0,
            )
            or 0
        )

        transcript_length = int(
            definition.get(
                "transcript_length",
                0,
            )
            or 0
        )

        effective_length = (
            cds_length
            if cds_length > 0
            else transcript_length
        )

        # Highest RNA TPM first.
        # Length and transcript ID are deterministic tie-breakers only.
        return (
            -rna_tpm,
            effective_length,
            transcript_id,
        )

    return sorted(
        transcript_ids,
        key=sort_key,
    )[0]


def find_indistinguishable_merges(
    transcript_ids,
    unique_read_indexes,
    multi_read_indexes,
    reference_definitions,
):
    """Find transcripts that are observationally indistinguishable.

    A transcript is eligible for merging only when:
      1. it has no observed transcript-unique read;
      2. another transcript also has no observed transcript-unique read;
      3. their observed non-empty multi-read sets are identical.

    Strict-subset relationships are NOT merged. This preserves cases such
    as A-B-C-D versus A-B-C, where reads unique to D can constrain the
    abundance of the longer transcript and shared A-B-C depth can still
    support the shorter transcript during EM.

    Within each indistinguishable set, the transcript with the highest
    RNA-seq TPM is retained as the representative.
    """

    transcript_ids = sorted(
        transcript_ids
    )

    eligible_ids = [
        transcript_id
        for transcript_id in transcript_ids
        if (
            not has_unique_evidence(
                transcript_id,
                unique_read_indexes,
            )
            and get_multi_read_evidence(
                transcript_id,
                multi_read_indexes,
            )
        )
    ]

    evidence_groups = defaultdict(list)

    for transcript_id in eligible_ids:
        evidence_key = frozenset(
            get_multi_read_evidence(
                transcript_id,
                multi_read_indexes,
            )
        )

        evidence_groups[
            evidence_key
        ].append(
            transcript_id
        )

    merge_map = {}

    for indistinguishable_ids in (
        evidence_groups.values()
    ):
        if len(indistinguishable_ids) <= 1:
            continue

        representative_id = (
            select_merge_representative(
                transcript_ids=(
                    indistinguishable_ids
                ),
                reference_definitions=(
                    reference_definitions
                ),
            )
        )

        for transcript_id in indistinguishable_ids:
            if transcript_id == representative_id:
                continue

            merge_map[
                transcript_id
            ] = representative_id

    return merge_map


# Backward-compatible name for callers outside this module.
def find_subset_merges(
    transcript_ids,
    unique_read_indexes,
    multi_read_indexes,
    reference_definitions,
):
    """Compatibility wrapper for the former merge function name."""

    return find_indistinguishable_merges(
        transcript_ids=transcript_ids,
        unique_read_indexes=unique_read_indexes,
        multi_read_indexes=multi_read_indexes,
        reference_definitions=reference_definitions,
    )

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


def _alignment_record_priority(
    record,
    merged_transcript_id,
):
    """Return deterministic priority for a merged alignment record."""

    original_transcript_id = record.get(
        "alignment_transcript_id",
        record["transcript_id"],
    )

    return (
        0 if original_transcript_id == merged_transcript_id else 1,
        int(record.get("hit_index", 0)),
        original_transcript_id,
        int(record.get("reference_start", 0)),
        int(record.get("reference_end", 0)),
    )

def apply_transcript_merges(
    group_records,
    merge_map,
):
    """Apply EM transcript merges while preserving alignment coordinates.

    ``transcript_id`` is the representative ID used by transcript EM.
    ``alignment_transcript_id`` remains the reference against which the
    original mapping coordinates and CIGAR were produced.
    """

    selected_records = {}

    for record in group_records:
        original_transcript_id = record["transcript_id"]
        merged_transcript_id = merge_map.get(
            original_transcript_id,
            original_transcript_id,
        )

        hit_key = (
            record["read_id"],
            merged_transcript_id,
        )

        current_record = selected_records.get(hit_key)

        if current_record is not None:
            current_priority = _alignment_record_priority(
                current_record,
                merged_transcript_id,
            )
            new_priority = _alignment_record_priority(
                record,
                merged_transcript_id,
            )

            if new_priority >= current_priority:
                continue

        updated_record = record.copy()
        updated_record["alignment_transcript_id"] = (
            original_transcript_id
        )
        updated_record["transcript_id"] = merged_transcript_id

        selected_records[hit_key] = updated_record

    return list(selected_records.values())

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

    merge_map = find_indistinguishable_merges(
        transcript_ids=transcript_ids,
        unique_read_indexes=unique_read_indexes,
        multi_read_indexes=multi_read_indexes,
        reference_definitions=reference_definitions,
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

            # Calculate expected read count.
            expected_count = 0.0
            seen_reads = set()

            for record in group_records:

                read_id = record["read_id"]

                if read_id in seen_reads:
                    continue

                seen_reads.add(
                    read_id
                )

                other_group_ids = record.get(
                    "other_group_ids",
                    [],
                )

                if other_group_ids is None:
                    other_group_ids = []

                expected_count += calculate_group_probability(
                    group_id=group_id,
                    other_group_ids=other_group_ids,
                    group_tpm=group_tpm,
                )

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
                reference_definitions[transcript_id].get(
                    "em_region_length",
                    original_cds_length if original_cds_length > 0 else transcript_length,
                )
            )

            return (
                {
                    transcript_id: expected_count,
                },
                {
                    transcript_id: original_cds_length,
                },
                {
                    transcript_id: em_cds_length,
                },
            )

    # --------------------------------------------------------
    # Original CDS lengths
    # --------------------------------------------------------

    original_cds_lengths = {
        transcript_id: int(
            reference_definitions[
                transcript_id
            ].get(
                "cds_length",
                0,
            )
        )
        for transcript_id in transcript_ids
    }

    # --------------------------------------------------------
    # Unique effective lengths
    # --------------------------------------------------------

    transcript_unique_lengths = determine_unique_transcript_lengths(
        transcript_ids=transcript_ids,
        reference_definitions=reference_definitions,
        unique_junction_length=50,
    )

    # --------------------------------------------------------
    # Updated lengths actually used by transcript EM
    # --------------------------------------------------------

    transcript_cds_lengths = determine_transcript_lengths(
        transcript_ids=transcript_ids,
        group_records=group_records,
        reference_definitions=reference_definitions,
        minimum_mapped_reads=10,
        minimum_mapped_span=100,
        minimum_transcript_fraction=1.0 / 3.0,
    )

    # --------------------------------------------------------
    # Transcript-level EM
    # --------------------------------------------------------

    (
        final_tpm,
        final_expected_counts,
    ) = execute_transcript_em_eachgroup(
        group_id=group_id,
        transcript_ids=transcript_ids,
        group_records=group_records,
        group_tpm=group_tpm,
        transcript_unique_lengths=transcript_unique_lengths,
        transcript_cds_lengths=transcript_cds_lengths,
        reads_info=reads_info,
        transcript_tpm=transcript_tpm,
        reference_definitions=reference_definitions,
    )

    return (
        final_expected_counts,
        original_cds_lengths,
        transcript_cds_lengths,
    )