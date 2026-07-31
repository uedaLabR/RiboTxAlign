from collections import defaultdict

from em_distribute.transcript_em_utils import (
    calculate_group_probability,
)


def prepare_read_hits(
    group_id,
    group_records,
    group_tpm,
):
    """
    Group mapping records by read ID.

    One read contains one or more candidate transcript hits.
    Duplicate read-transcript combinations are removed.
    """

    read_data = {}

    for record in group_records:
        read_id = record["read_id"]

        if read_id not in read_data:
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

            read_data[read_id] = {
                "read_id": read_id,
                "prior_probability": prior_probability,
                "candidate_hits": {},
            }

        transcript_id = record["transcript_id"]

        # Keep only one record for each read-transcript pair.
        if transcript_id in read_data[read_id]["candidate_hits"]:
            continue

        read_data[read_id]["candidate_hits"][
            transcript_id
        ] = {
            "transcript_id": transcript_id,
            "start": record["reference_start"],
            "end": record["reference_end"],
        }

    read_hits = []

    for read_info in read_data.values():
        read_info["candidate_hits"] = list(
            read_info["candidate_hits"].values()
        )

        read_hits.append(
            read_info
        )

    return read_hits


def calculate_weighted_unique_counts(
    transcript_ids,
    read_hits,
):
    """
    Count transcript-unique reads weighted by group prior probability.
    """

    unique_counts = {
        transcript_id: 0.0
        for transcript_id in transcript_ids
    }

    for read_info in read_hits:
        candidate_hits = read_info[
            "candidate_hits"
        ]

        if len(candidate_hits) != 1:
            continue

        transcript_id = candidate_hits[0][
            "transcript_id"
        ]

        unique_counts[transcript_id] += read_info[
            "prior_probability"
        ]

    return unique_counts


def initialize_transcript_counts(
    transcript_ids,
    unique_counts,
    zero_start_fraction=1.0 / 3.0,
):
    """
    Initialize expected counts mainly from transcript-unique reads.

    Transcripts without unique reads start from one third of the
    smallest positive unique count.
    """

    positive_unique_counts = [
        count
        for count in unique_counts.values()
        if count > 0.0
    ]

    if positive_unique_counts:
        zero_start_count = (
            min(positive_unique_counts)
            * zero_start_fraction
        )
    else:
        zero_start_count = 1.0

    initial_counts = {}

    for transcript_id in transcript_ids:
        unique_count = unique_counts.get(
            transcript_id,
            0.0,
        )

        if unique_count > 0.0:
            initial_counts[
                transcript_id
            ] = unique_count
        else:
            initial_counts[
                transcript_id
            ] = zero_start_count

    return initial_counts


def calculate_transcript_tpm(
    transcript_ids,
    expected_counts,
    transcript_cds_lengths,
    group_tpm_value,
):
    """
    Calculate transcript TPM values scaled to the group TPM.
    """

    transcript_rates = {}

    for transcript_id in transcript_ids:
        cds_length = max(
            int(
                transcript_cds_lengths[
                    transcript_id
                ]
            ),
            1,
        )

        transcript_rates[
            transcript_id
        ] = (
            expected_counts.get(
                transcript_id,
                0.0,
            )
            / cds_length
        )

    total_rate = sum(
        transcript_rates.values()
    )

    if total_rate <= 0.0:
        equal_tpm = (
            group_tpm_value
            / len(transcript_ids)
        )

        return {
            transcript_id: equal_tpm
            for transcript_id in transcript_ids
        }

    return {
        transcript_id: (
            group_tpm_value
            * transcript_rates[transcript_id]
            / total_rate
        )
        for transcript_id in transcript_ids
    }


def calculate_transcript_probabilities(
    candidate_hits,
    current_tpm,
):
    """
    Calculate transcript probabilities for one read.
    """

    candidate_transcript_ids = [
        hit["transcript_id"]
        for hit in candidate_hits
    ]

    if len(candidate_transcript_ids) == 1:
        return {
            candidate_transcript_ids[0]: 1.0
        }

    total_tpm = sum(
        current_tpm.get(
            transcript_id,
            0.0,
        )
        for transcript_id in candidate_transcript_ids
    )

    if total_tpm <= 0.0:
        equal_probability = (
            1.0
            / len(candidate_transcript_ids)
        )

        return {
            transcript_id: equal_probability
            for transcript_id in candidate_transcript_ids
        }

    return {
        transcript_id: (
            current_tpm.get(
                transcript_id,
                0.0,
            )
            / total_tpm
        )
        for transcript_id in candidate_transcript_ids
    }


def distribute_reads(
    transcript_ids,
    read_hits,
    current_tpm,
):
    """
    Distribute read prior probabilities among candidate transcripts.
    """

    expected_counts = {
        transcript_id: 0.0
        for transcript_id in transcript_ids
    }

    for read_info in read_hits:
        prior_probability = read_info[
            "prior_probability"
        ]

        transcript_probabilities = (
            calculate_transcript_probabilities(
                candidate_hits=read_info[
                    "candidate_hits"
                ],
                current_tpm=current_tpm,
            )
        )

        for transcript_id, probability in (
            transcript_probabilities.items()
        ):
            expected_counts[
                transcript_id
            ] += (
                prior_probability
                * probability
            )

    return expected_counts


def calculate_tpm_difference(
    previous_tpm,
    current_tpm,
):
    """Return the maximum absolute TPM difference."""

    return max(
        abs(
            current_tpm[transcript_id]
            - previous_tpm.get(
                transcript_id,
                0.0,
            )
        )
        for transcript_id in current_tpm
    )


def append_final_read_information(
    group_id,
    read_hits,
    final_tpm,
    reads_info,
):
    """Append final read-transcript probabilities."""

    for read_info in read_hits:
        candidate_hits = read_info[
            "candidate_hits"
        ]

        transcript_probabilities = (
            calculate_transcript_probabilities(
                candidate_hits=candidate_hits,
                current_tpm=final_tpm,
            )
        )

        candidate_transcript_ids = [
            hit["transcript_id"]
            for hit in candidate_hits
        ]

        for hit in candidate_hits:
            transcript_id = hit[
                "transcript_id"
            ]

            other_transcript_ids = [
                candidate_id
                for candidate_id
                in candidate_transcript_ids
                if candidate_id != transcript_id
            ]

            reads_info.append(
                {
                    "group_id": group_id,
                    "read_id": read_info[
                        "read_id"
                    ],
                    "transcript_id": transcript_id,
                    "start": hit["start"],
                    "end": hit["end"],
                    "prior_probability": read_info[
                        "prior_probability"
                    ],
                    "transcript_probability":
                        transcript_probabilities[
                            transcript_id
                        ],
                    "other_transcript_ids":
                        other_transcript_ids,
                }
            )


def execute_transcript_em_eachgroup(
    group_id,
    transcript_ids,
    group_records,
    group_tpm,
    transcript_cds_lengths,
    reads_info,
    transcript_tpm,
    max_iterations=200,
    convergence_threshold=1.0e-6,
    zero_start_fraction=1.0 / 3.0,
):
    """Execute transcript-level EM for one transcript group."""

    transcript_ids = set(
        transcript_ids
    )

    group_tpm_value = group_tpm.get(
        group_id,
        0.0,
    )

    read_hits = prepare_read_hits(
        group_id=group_id,
        group_records=group_records,
        group_tpm=group_tpm,
    )

    unique_counts = calculate_weighted_unique_counts(
        transcript_ids=transcript_ids,
        read_hits=read_hits,
    )

    initial_counts = initialize_transcript_counts(
        transcript_ids=transcript_ids,
        unique_counts=unique_counts,
        zero_start_fraction=zero_start_fraction,
    )

    current_tpm = calculate_transcript_tpm(
        transcript_ids=transcript_ids,
        expected_counts=initial_counts,
        transcript_cds_lengths=transcript_cds_lengths,
        group_tpm_value=group_tpm_value,
    )

    final_expected_counts = initial_counts

    for iteration in range(
        1,
        max_iterations + 1,
    ):
        expected_counts = distribute_reads(
            transcript_ids=transcript_ids,
            read_hits=read_hits,
            current_tpm=current_tpm,
        )

        new_tpm = calculate_transcript_tpm(
            transcript_ids=transcript_ids,
            expected_counts=expected_counts,
            transcript_cds_lengths=transcript_cds_lengths,
            group_tpm_value=group_tpm_value,
        )

        difference = calculate_tpm_difference(
            previous_tpm=current_tpm,
            current_tpm=new_tpm,
        )

        current_tpm = new_tpm
        final_expected_counts = expected_counts

        if difference <= convergence_threshold:
            break

    for transcript_id in transcript_ids:
        transcript_tpm[
            transcript_id
        ] = current_tpm[
            transcript_id
        ]

    append_final_read_information(
        group_id=group_id,
        read_hits=read_hits,
        final_tpm=current_tpm,
        reads_info=reads_info,
    )

    print(
        "transcript EM completed",
        group_id,
        "iterations",
        iteration,
        "transcripts",
        len(transcript_ids),
        "reads",
        len(read_hits),
        "difference",
        difference,
    )

    return (
        current_tpm,
        final_expected_counts,
    )