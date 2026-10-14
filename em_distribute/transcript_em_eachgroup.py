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
    group_probability_cache = {}

    for record in group_records:
        read_id = record["read_id"]

        if read_id not in read_data:
            other_group_ids = record.get(
                "other_group_ids",
                [],
            )

            if other_group_ids is None:
                other_group_ids = []

            # Group TPM is fixed throughout transcript EM. Reuse the same
            # posterior for an identical ordered candidate-group pattern.
            group_pattern = tuple(other_group_ids)
            prior_probability = group_probability_cache.get(group_pattern)
            if prior_probability is None:
                prior_probability = calculate_group_probability(
                    group_id=group_id,
                    other_group_ids=other_group_ids,
                    group_tpm=group_tpm,
                )
                if len(group_probability_cache) >= 4096:
                    group_probability_cache.clear()
                group_probability_cache[group_pattern] = prior_probability

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
            "alignment_transcript_id": record.get(
                "alignment_transcript_id",
                transcript_id,
            ),
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


def calculate_rna_tpm_weights(
    transcript_ids,
    reference_definitions,
    minimum_weight=1.0e-6,
):
    """Normalize RNA-seq TPM within one transcript group.

    A very small floor prevents a transcript with zero RNA TPM from becoming
    an absorbing zero state during transcript EM.
    """

    transcript_ids = sorted(
        transcript_ids
    )

    rna_tpm = {
        transcript_id: get_reference_rna_tpm(
            transcript_id,
            reference_definitions,
        )
        for transcript_id in transcript_ids
    }

    total_tpm = sum(
        rna_tpm.values()
    )

    if total_tpm <= 0.0:
        equal_weight = (
            1.0
            / len(transcript_ids)
        )

        return {
            transcript_id: equal_weight
            for transcript_id in transcript_ids
        }

    weights = {
        transcript_id: (
            rna_tpm[transcript_id]
            / total_tpm
        )
        for transcript_id in transcript_ids
    }

    if minimum_weight > 0.0:
        weights = {
            transcript_id: max(
                weight,
                minimum_weight,
            )
            for transcript_id, weight in (
                weights.items()
            )
        }

        total_weight = sum(
            weights.values()
        )

        weights = {
            transcript_id: (
                weight
                / total_weight
            )
            for transcript_id, weight in (
                weights.items()
            )
        }

    return weights


def calculate_unique_fractions(
    transcript_ids,
    transcript_unique_lengths,
    transcript_cds_lengths,
):
    """Calculate the distinguishable fraction of each effective transcript."""

    unique_fractions = {}

    for transcript_id in transcript_ids:
        unique_length = max(
            float(
                transcript_unique_lengths.get(
                    transcript_id,
                    0,
                )
                or 0
            ),
            0.0,
        )

        transcript_length = max(
            float(
                transcript_cds_lengths.get(
                    transcript_id,
                    0,
                )
                or 0
            ),
            1.0,
        )

        unique_fractions[
            transcript_id
        ] = min(
            unique_length
            / transcript_length,
            1.0,
        )

    return unique_fractions


def calculate_effective_group_read_count(
    read_hits,
):
    """Return group-level expected read count from group prior probabilities."""

    effective_count = sum(
        max(
            float(
                read_info.get(
                    "prior_probability",
                    0.0,
                )
                or 0.0
            ),
            0.0,
        )
        for read_info in read_hits
    )

    if effective_count > 0.0:
        return effective_count

    return float(
        len(read_hits)
    )


def calculate_rna_count_weights(
    transcript_ids,
    rna_tpm_weights,
    transcript_cds_lengths,
):
    """Convert RNA TPM weights to expected Ribo read-count fractions.

    RNA TPM is a length-normalized abundance. Under the initialization model,
    expected total footprint count is proportional to RNA TPM times the
    effective translated length.
    """

    raw_weights = {}

    for transcript_id in transcript_ids:
        transcript_length = max(
            float(
                transcript_cds_lengths.get(
                    transcript_id,
                    0,
                )
                or 0
            ),
            1.0,
        )

        raw_weights[
            transcript_id
        ] = (
            rna_tpm_weights[
                transcript_id
            ]
            * transcript_length
        )

    total_weight = sum(
        raw_weights.values()
    )

    if total_weight <= 0.0:
        equal_weight = (
            1.0
            / len(transcript_ids)
        )

        return {
            transcript_id: equal_weight
            for transcript_id in transcript_ids
        }

    return {
        transcript_id: (
            raw_weights[transcript_id]
            / total_weight
        )
        for transcript_id in transcript_ids
    }


def calculate_poisson_initial_abundance(
    transcript_ids,
    unique_counts,
    unique_fractions,
    rna_count_weights,
    effective_group_read_count,
    rna_prior_strength=5.0,
):
    """Estimate initial transcript read counts with a Gamma-Poisson model.

    Let n_i be the total expected footprint count from transcript i and f_i
    its distinguishable fraction. The observed unique count is modeled as

        U_i ~ Poisson(f_i * n_i)

    RNA-seq supplies a Gamma prior for n_i. ``rna_prior_strength`` is the
    total prior pseudo-count across the transcript group.

    The returned value is the posterior mean of n_i.
    """

    transcript_ids = sorted(
        transcript_ids
    )

    if rna_prior_strength <= 0.0:
        raise ValueError(
            "rna_prior_strength must be positive."
        )

    group_count = max(
        float(
            effective_group_read_count
        ),
        1.0e-12,
    )

    prior_rate = (
        rna_prior_strength
        / group_count
    )

    initial_abundance = {}

    for transcript_id in transcript_ids:
        unique_fraction = max(
            float(
                unique_fractions.get(
                    transcript_id,
                    0.0,
                )
                or 0.0
            ),
            0.0,
        )

        unique_count = max(
            float(
                unique_counts.get(
                    transcript_id,
                    0.0,
                )
                or 0.0
            ),
            0.0,
        )

        # A positive unique count with zero modeled unique fraction is
        # internally inconsistent. Ignore that count for initialization
        # rather than creating an infinite density.
        if unique_fraction <= 0.0:
            unique_count = 0.0

        prior_shape = (
            rna_prior_strength
            * rna_count_weights[
                transcript_id
            ]
        )

        posterior_shape = (
            unique_count
            + prior_shape
        )

        posterior_rate = (
            unique_fraction
            + prior_rate
        )

        initial_abundance[
            transcript_id
        ] = (
            posterior_shape
            / posterior_rate
        )

    return initial_abundance


def convert_abundance_to_density_rates(
    transcript_ids,
    transcript_abundance,
    transcript_cds_lengths,
):
    """Convert expected transcript read counts to length-normalized rates."""

    return {
        transcript_id: (
            transcript_abundance[
                transcript_id
            ]
            / max(
                float(
                    transcript_cds_lengths.get(
                        transcript_id,
                        0,
                    )
                    or 0
                ),
                1.0,
            )
        )
        for transcript_id in transcript_ids
    }


def scale_initial_rates_to_group_tpm(
    transcript_ids,
    initial_rates,
    group_tpm_value,
):
    """Scale positive initial transcript rates so their TPM sum equals group TPM."""

    total_rate = sum(
        max(
            float(
                initial_rates.get(
                    transcript_id,
                    0.0,
                )
                or 0.0
            ),
            0.0,
        )
        for transcript_id in transcript_ids
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
            * max(
                float(
                    initial_rates.get(
                        transcript_id,
                        0.0,
                    )
                    or 0.0
                ),
                0.0,
            )
            / total_rate
        )
        for transcript_id in transcript_ids
    }


def initialize_transcript_tpm(
    transcript_ids,
    unique_counts,
    transcript_unique_lengths,
    transcript_cds_lengths,
    reference_definitions,
    read_hits,
    group_tpm_value,
    rna_prior_strength=5.0,
):
    """Initialize transcript TPM from unique Ribo evidence and RNA-seq TPM.

    The function deliberately delegates each conceptual step to a small helper
    so the initialization model can later be changed without modifying EM.
    """

    transcript_ids = sorted(
        transcript_ids
    )

    rna_tpm_weights = calculate_rna_tpm_weights(
        transcript_ids=transcript_ids,
        reference_definitions=reference_definitions,
    )

    unique_fractions = calculate_unique_fractions(
        transcript_ids=transcript_ids,
        transcript_unique_lengths=(
            transcript_unique_lengths
        ),
        transcript_cds_lengths=(
            transcript_cds_lengths
        ),
    )

    effective_group_read_count = (
        calculate_effective_group_read_count(
            read_hits=read_hits,
        )
    )

    rna_count_weights = (
        calculate_rna_count_weights(
            transcript_ids=transcript_ids,
            rna_tpm_weights=rna_tpm_weights,
            transcript_cds_lengths=(
                transcript_cds_lengths
            ),
        )
    )

    transcript_abundance = (
        calculate_poisson_initial_abundance(
            transcript_ids=transcript_ids,
            unique_counts=unique_counts,
            unique_fractions=unique_fractions,
            rna_count_weights=(
                rna_count_weights
            ),
            effective_group_read_count=(
                effective_group_read_count
            ),
            rna_prior_strength=(
                rna_prior_strength
            ),
        )
    )

    initial_rates = (
        convert_abundance_to_density_rates(
            transcript_ids=transcript_ids,
            transcript_abundance=(
                transcript_abundance
            ),
            transcript_cds_lengths=(
                transcript_cds_lengths
            ),
        )
    )

    return scale_initial_rates_to_group_tpm(
        transcript_ids=transcript_ids,
        initial_rates=initial_rates,
        group_tpm_value=group_tpm_value,
    )

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


def aggregate_em_read_patterns(read_hits):
    """Combine weights with identical candidate sets for EM only.

    prepare_read_hits has already removed duplicate read-transcript pairs.
    The E-step uses only candidate IDs and prior_probability; alignment
    details remain in the original read_hits for final per-read output.
    Summation order differs from the read-wise E-step, so floating-point
    roundoff (and a stopping iteration at the tolerance boundary) can differ.
    """
    pattern_weights = {}
    for read_info in read_hits:
        candidates = tuple(sorted(
            hit["transcript_id"] for hit in read_info["candidate_hits"]
        ))
        pattern_weights[candidates] = (
            pattern_weights.get(candidates, 0.0)
            + read_info["prior_probability"]
        )
    return [
        {
            "prior_probability": weight,
            "candidate_hits": [
                {"transcript_id": transcript_id}
                for transcript_id in candidates
            ],
        }
        for candidates, weight in pattern_weights.items()
    ]


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

    pattern_cache = {}
    for read_info in read_hits:
        candidate_hits = read_info[
            "candidate_hits"
        ]

        # Keep candidate order in the key so posterior sums retain the exact
        # same floating-point order. Alignment coordinates remain read-specific.
        candidate_transcript_ids = tuple(
            hit["transcript_id"]
            for hit in candidate_hits
        )
        cached = pattern_cache.get(candidate_transcript_ids)
        if cached is None:
            transcript_probabilities = calculate_transcript_probabilities(
                candidate_hits=candidate_hits,
                current_tpm=final_tpm,
            )
            others_by_transcript = {
                tid: [other for other in candidate_transcript_ids if other != tid]
                for tid in candidate_transcript_ids
            }
            cached = (transcript_probabilities, others_by_transcript)
            if len(pattern_cache) >= 4096:
                pattern_cache.clear()
            pattern_cache[candidate_transcript_ids] = cached
        transcript_probabilities, others_by_transcript = cached

        for hit in candidate_hits:
            transcript_id = hit[
                "transcript_id"
            ]

            # Do not expose shared mutable lists between output records.
            other_transcript_ids = others_by_transcript[transcript_id].copy()

            reads_info.append(
                {
                    "group_id": group_id,
                    "read_id": read_info[
                        "read_id"
                    ],
                    "transcript_id": transcript_id,
                    "alignment_transcript_id": hit[
                        "alignment_transcript_id"
                    ],
                    "start": hit["start"],
                    "end": hit["end"],
                    "read_length": hit[
                        "read_length"
                    ],
                    "hit_index": hit["hit_index"],
                    "query_start": hit["query_start"],
                    "query_end": hit["query_end"],
                    "mapping_strand": hit["mapping_strand"],
                    "cigar": hit["cigar"],
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
    transcript_unique_lengths,
    transcript_cds_lengths,
    reads_info,
    transcript_tpm,
    reference_definitions,
    max_iterations=200,
    convergence_threshold=1.0e-6,
    rna_prior_strength=5.0,
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

    current_tpm = initialize_transcript_tpm(
        transcript_ids=transcript_ids,
        unique_counts=unique_counts,
        transcript_unique_lengths=(
            transcript_unique_lengths
        ),
        transcript_cds_lengths=(
            transcript_cds_lengths
        ),
        reference_definitions=(
            reference_definitions
        ),
        read_hits=read_hits,
        group_tpm_value=group_tpm_value,
        rna_prior_strength=(
            rna_prior_strength
        ),
    )

    # Build once, outside the iteration loop. Keep read_hits intact for
    # initialization and the final alignment-specific probability output.
    em_read_patterns = aggregate_em_read_patterns(read_hits)
    converged = False

    for iteration in range(
        1,
        max_iterations + 1,
    ):
        expected_counts = distribute_reads(
            transcript_ids=transcript_ids,
            read_hits=em_read_patterns,
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

        if difference <= convergence_threshold:
            converged = True
            break

    # The last M-step changed current_tpm. Recompute counts from that same
    # final state used below for read probabilities, even at the iteration cap.
    final_expected_counts = distribute_reads(
        transcript_ids=transcript_ids,
        read_hits=em_read_patterns,
        current_tpm=current_tpm,
    )

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
        "patterns",
        len(em_read_patterns),
        "difference",
        difference,
        "converged",
        converged,
    )

    return (
        current_tpm,
        final_expected_counts,
    )
