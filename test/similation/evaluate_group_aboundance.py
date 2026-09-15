import pandas as pd
import numpy as np

from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error


# ------------------------------------------------------------
# Build truth read counts per RiboTxAlign group
# ------------------------------------------------------------

def get_truth_group_counts(
    profile_path,
    metadata_path,
):

    # --------------------------------------------------------
    # Read simulation truth
    # --------------------------------------------------------

    records = []

    with open(profile_path) as f:
        tid = None

        for line in f:
            line = line.strip()

            if line.startswith("tid:"):
                tid = (
                    line.split(":", 1)[1]
                    .strip()
                )

            elif line.startswith("rprofile:"):

                profile = [
                    int(x)
                    for x in line.split(":", 1)[1].split()
                ]

                records.append({
                    "transcript_id": tid,
                    "truth_count": sum(profile),
                })

    truth_df = pd.DataFrame(records)

    truth_df["transcript_id"] = (
        truth_df["transcript_id"]
        .str.split(".")
        .str[0]
    )

    # --------------------------------------------------------
    # Read RiboTxAlign metadata
    # --------------------------------------------------------

    meta_df = pd.read_csv(
        metadata_path,
        sep="\t",
    )

    meta_df["transcript_id"] = (
        meta_df["representative_transcript_id"]
        .str.split(".")
        .str[0]
    )

    # Convert:
    # chr10_minus_G00012_LARP4B
    # ->
    # chr10_minus_G00012
    meta_df["group_id"] = (
        meta_df["overlap_group_id"]
        .str.extract(
            r"^(chr[^_]+_(?:plus|minus)_G\d+)",
            expand=False,
        )
    )

    meta_df = meta_df[
        [
            "transcript_id",
            "group_id",
        ]
    ].drop_duplicates()

    # --------------------------------------------------------
    # Number of transcripts in each group
    # --------------------------------------------------------

    group_sizes = (
        meta_df.groupby(
            "group_id"
        )
        .size()
        .reset_index(
            name="num_transcripts"
        )
    )

    # --------------------------------------------------------
    # Map simulation transcripts to RiboTxAlign groups
    # --------------------------------------------------------

    mapped_df = truth_df.merge(
        meta_df,
        on="transcript_id",
        how="inner",
    )

    # --------------------------------------------------------
    # Sum simulated reads within each group
    # --------------------------------------------------------

    group_truth = (
        mapped_df.groupby(
            "group_id",
            as_index=False,
        )["truth_count"]
        .sum()
    )

    group_truth = group_truth.merge(
        group_sizes,
        on="group_id",
        how="left",
    )

    return group_truth


# ------------------------------------------------------------
# Calculate accuracy
# ------------------------------------------------------------

def calculate_accuracy(
    df,
    label,
):

    if len(df) < 2:
        print()
        print(label)
        print("Not enough groups.")
        return

    truth = (
        df["truth_count"]
        .to_numpy()
    )

    pred = (
        df["final_expected_count"]
        .to_numpy()
    )

    pearson = pearsonr(
        truth,
        pred,
    ).statistic

    spearman = spearmanr(
        truth,
        pred,
    ).statistic

    mae = mean_absolute_error(
        truth,
        pred,
    )

    rmse = np.sqrt(
        mean_squared_error(
            truth,
            pred,
        )
    )

    print()
    print("=" * 70)
    print(label)
    print("=" * 70)

    print(
        "Groups:",
        len(df),
    )

    print(
        "Truth reads:",
        truth.sum(),
    )

    print(
        "Estimated reads:",
        pred.sum(),
    )

    print()

    print(
        f"Pearson : {pearson:.6f}"
    )

    print(
        f"Spearman: {spearman:.6f}"
    )

    print(
        f"MAE     : {mae:.6f}"
    )

    print(
        f"RMSE    : {rmse:.6f}"
    )


# ------------------------------------------------------------
# Compare truth with group EM
# ------------------------------------------------------------

def evaluate_group_em(
    profile_path,
    metadata_path,
    group_em_path,
):

    truth_df = get_truth_group_counts(
        profile_path,
        metadata_path,
    )

    em_df = pd.read_csv(
        group_em_path,
        sep="\t",
    )

    # --------------------------------------------------------
    # Check tables before merge
    # --------------------------------------------------------

    print()
    print("=== Truth group table ===")
    print(
        truth_df.head(10)
    )

    print()
    print("=== Group EM table ===")
    print(
        em_df[
            [
                "group_id",
                "unique_read_count",
                "multi_hit_record_count",
                "final_expected_count",
                "initial_tpm",
                "final_tpm",
            ]
        ].head(10)
    )

    # --------------------------------------------------------
    # Merge truth and group EM
    # --------------------------------------------------------

    df = truth_df.merge(
        em_df[
            [
                "group_id",
                "unique_read_count",
                "multi_hit_record_count",
                "final_expected_count",
                "initial_tpm",
                "final_tpm",
            ]
        ],
        on="group_id",
        how="left",
    )

    # --------------------------------------------------------
    # Check matched groups
    # --------------------------------------------------------

    matched = (
        df["final_expected_count"]
        .notna()
        .sum()
    )

    print()
    print(
        "Matched groups:",
        matched,
        "/",
        len(df),
    )

    df["final_expected_count"] = (
        df["final_expected_count"]
        .fillna(0)
    )

    # --------------------------------------------------------
    # Group composition
    # --------------------------------------------------------

    single_groups = df[
        df["num_transcripts"] == 1
    ].copy()

    multi_groups = df[
        df["num_transcripts"] > 1
    ].copy()

    print()
    print("=" * 70)
    print("Group transcript composition")
    print("=" * 70)

    print(
        "Total groups:",
        len(df),
    )

    print(
        "Single-transcript groups:",
        len(single_groups),
        f"({len(single_groups) / len(df) * 100:.2f}%)",
    )

    print(
        "Multi-transcript groups:",
        len(multi_groups),
        f"({len(multi_groups) / len(df) * 100:.2f}%)",
    )

    print()

    print(
        "Truth reads in single-transcript groups:",
        single_groups["truth_count"].sum(),
    )

    print(
        "Truth reads in multi-transcript groups:",
        multi_groups["truth_count"].sum(),
    )

    print()

    print(
        "Fraction of truth reads in single-transcript groups:",
        f"{single_groups['truth_count'].sum() / df['truth_count'].sum() * 100:.2f}%",
    )

    # --------------------------------------------------------
    # Accuracy: all groups
    # --------------------------------------------------------

    calculate_accuracy(
        df,
        "All groups",
    )

    # --------------------------------------------------------
    # Accuracy: single-transcript groups
    # --------------------------------------------------------

    calculate_accuracy(
        single_groups,
        "Single-transcript groups only",
    )

    # --------------------------------------------------------
    # Accuracy: multi-transcript groups
    # --------------------------------------------------------

    calculate_accuracy(
        multi_groups,
        "Multi-transcript groups only",
    )

    # --------------------------------------------------------
    # Groups with no unique reads but multi-hit reads
    # --------------------------------------------------------

    no_unique = df[
        (
            df["unique_read_count"] == 0
        )
        &
        (
            df["multi_hit_record_count"] > 0
        )
    ].copy()

    print()
    print("=" * 70)
    print(
        "Groups with no unique reads but multi-hit reads"
    )
    print("=" * 70)

    print(
        "Number of groups:",
        len(no_unique),
    )

    print(
        "Truth reads:",
        no_unique["truth_count"].sum(),
    )

    print(
        "Multi-hit records:",
        no_unique["multi_hit_record_count"].sum(),
    )

    print(
        "Final expected count:",
        no_unique["final_expected_count"].sum(),
    )

    print()

    print(
        no_unique[
            [
                "group_id",
                "num_transcripts",
                "truth_count",
                "unique_read_count",
                "multi_hit_record_count",
                "final_expected_count",
            ]
        ].head(20)
    )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def compare():

    base = (
        "/mnt/share/project/riboML/"
        "validation_expression/similation/"
        "HEK293T"
    )

    profile_path = (
        f"{base}/"
        "HEK293T.simulation.profile"
    )

    metadata_path = (
        f"{base}/"
        "RiboTxAlign/"
        "HEK293T_simulation/"
        "reference/"
        "HEK293T_simulation.metadata.tsv"
    )

    group_em_path = (
        f"{base}/"
        "RiboTxAlign/"
        "HEK293T_simulation/"
        "em/"
        "HEK293T_simulation.group_tpm.tsv"
    )

    evaluate_group_em(
        profile_path,
        metadata_path,
        group_em_path,
    )


compare()