import pandas as pd
import numpy as np

from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error


# ------------------------------------------------------------
# Read simulation truth count
# ------------------------------------------------------------

def read_simulation_truth(
    profile_path,
):
    """Read transcript-level simulated footprint counts."""

    records = []

    with open(profile_path) as f:

        transcript_id = None

        for line in f:

            line = line.strip()

            if line.startswith("tid:"):

                transcript_id = (
                    line.split(":", 1)[1]
                    .strip()
                )

            elif line.startswith("rprofile:"):

                profile = [
                    int(x)
                    for x in line.split(":", 1)[1].split()
                ]

                records.append(
                    {
                        "transcript_id": transcript_id,
                        "truth_count": sum(profile),
                    }
                )

    df = pd.DataFrame(
        records
    )

    # Remove transcript version
    df["transcript_id"] = (
        df["transcript_id"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    # Sum just in case the same transcript appears more than once
    df = (
        df.groupby(
            "transcript_id",
            as_index=False,
        )["truth_count"]
        .sum()
    )

    return df


# ------------------------------------------------------------
# Read simulator input abundance
# ------------------------------------------------------------

def read_simulation_abundance(
    abundance_path,
):
    """Read TPM supplied to the footprint simulator."""

    df = pd.read_csv(
        abundance_path,
        sep="\t",
    )

    required_columns = {
        "Name",
        "Length",
        "EffectiveLength",
        "TPM",
    }

    missing_columns = (
        required_columns
        - set(df.columns)
    )

    if missing_columns:
        raise ValueError(
            "Missing columns in simulation abundance file: "
            + ", ".join(
                sorted(missing_columns)
            )
        )

    df = df.rename(
        columns={
            "Name": "transcript_id",
            "Length": "simulation_length",
            "EffectiveLength": "simulation_effective_length",
            "TPM": "truth_tpm_file",
        }
    )

    df["transcript_id"] = (
        df["transcript_id"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    return df[
        [
            "transcript_id",
            "simulation_length",
            "simulation_effective_length",
            "truth_tpm_file",
        ]
    ].copy()


# ------------------------------------------------------------
# Read RiboTxAlign metadata
# ------------------------------------------------------------

def read_metadata(
    metadata_path,
):
    """Build RiboTxAlign reference -> representative ENST mapping."""

    df = pd.read_csv(
        metadata_path,
        sep="\t",
    )

    # RiboTxAlign custom transcript ID
    df["reference_id"] = (
        df["reference_id"]
        .astype(str)
    )

    # Representative original transcript ID
    df["transcript_id"] = (
        df["representative_transcript_id"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    # Convert:
    #
    # chr10_minus_G00012_LARP4B
    # ->
    # chr10_minus_G00012
    df["group_id"] = (
        df["overlap_group_id"]
        .astype(str)
        .str.extract(
            r"^(chr[^_]+_(?:plus|minus)_G\d+)",
            expand=False,
        )
    )

    df = df[
        [
            "reference_id",
            "transcript_id",
            "group_id",
        ]
    ].drop_duplicates()

    # --------------------------------------------------------
    # Count number of RiboTxAlign transcripts in each group
    # --------------------------------------------------------

    group_sizes = (
        df.groupby(
            "group_id"
        )["reference_id"]
        .nunique()
        .reset_index(
            name="num_transcripts"
        )
    )

    df = df.merge(
        group_sizes,
        on="group_id",
        how="left",
    )

    return df


# ------------------------------------------------------------
# Read RiboTxAlign transcript results
# ------------------------------------------------------------

def read_ribotxalign_results(
    transcript_tpm_path,
    metadata_path,
):
    """Read RiboTxAlign transcript count and TPM."""

    result_df = pd.read_csv(
        transcript_tpm_path,
        sep="\t",
    )

    metadata_df = read_metadata(
        metadata_path
    )

    # --------------------------------------------------------
    # Join custom RiboTxAlign ID to representative ENST
    # --------------------------------------------------------

    df = result_df.merge(
        metadata_df,
        left_on="transcript_id",
        right_on="reference_id",
        how="left",
        suffixes=(
            "_ribotx",
            "_gencode",
        ),
    )

    # --------------------------------------------------------
    # Check metadata mapping
    # --------------------------------------------------------

    missing = df[
        df["transcript_id_gencode"].isna()
    ]

    print()
    print("=" * 70)
    print("Metadata mapping")
    print("=" * 70)

    print(
        "RiboTxAlign transcripts:",
        len(result_df),
    )

    print(
        "Mapped to metadata:",
        df["transcript_id_gencode"]
        .notna()
        .sum(),
    )

    print(
        "Missing metadata mapping:",
        len(missing),
    )

    if len(missing) > 0:

        print()
        print(
            missing[
                [
                    "transcript_id_ribotx",
                    "expected_count",
                    "tpm",
                ]
            ].head(20)
        )

    df = df[
        df["transcript_id_gencode"].notna()
    ].copy()

    # --------------------------------------------------------
    # Rename representative ENST
    # --------------------------------------------------------

    df = df.rename(
        columns={
            "transcript_id_gencode": "transcript_id",
        }
    )

    # --------------------------------------------------------
    # Check duplicate ENST mappings
    # --------------------------------------------------------

    duplicate_mapping = (
        df.groupby(
            "transcript_id"
        )["reference_id"]
        .nunique()
    )

    duplicate_mapping = duplicate_mapping[
        duplicate_mapping > 1
    ]

    print()
    print(
        "ENST IDs represented by multiple "
        "RiboTxAlign references:",
        len(duplicate_mapping),
    )

    if len(duplicate_mapping) > 0:

        print()
        print(
            duplicate_mapping
            .sort_values(
                ascending=False
            )
            .head(20)
        )

    # --------------------------------------------------------
    # Aggregate if multiple RiboTxAlign references map
    # to the same representative ENST
    # --------------------------------------------------------

    transcript_results = (
        df.groupby(
            "transcript_id",
            as_index=False,
        )
        .agg(
            ribotx_count=(
                "expected_count",
                "sum",
            ),
            ribotx_tpm=(
                "tpm",
                "sum",
            ),
            original_cds_length=(
                "original_cds_length",
                "max",
            ),
            em_cds_length=(
                "em_cds_length",
                "max",
            ),
            num_transcripts=(
                "num_transcripts",
                "max",
            ),
        )
    )

    return transcript_results


# ------------------------------------------------------------
# Calculate count accuracy
# ------------------------------------------------------------

def calculate_count_accuracy(
    df,
    label,
):

    valid_df = df[
        df["ribotx_count"].notna()
    ].copy()

    if len(valid_df) < 2:

        print()
        print(label)
        print("Not enough transcripts.")
        return

    truth = (
        valid_df["truth_count"]
        .to_numpy(
            dtype=float
        )
    )

    pred = (
        valid_df["ribotx_count"]
        .to_numpy(
            dtype=float
        )
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
        "Transcripts:",
        len(valid_df),
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
# Calculate TPM correlation
# ------------------------------------------------------------

def calculate_tpm_correlation(
    df,
    truth_column,
    label,
):

    valid_df = df[
        df[truth_column].notna()
        &
        df["ribotx_tpm"].notna()
    ].copy()

    if len(valid_df) < 2:

        print()
        print(label)
        print("Not enough transcripts.")
        return

    truth = (
        valid_df[truth_column]
        .to_numpy(
            dtype=float
        )
    )

    pred = (
        valid_df["ribotx_tpm"]
        .to_numpy(
            dtype=float
        )
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
        "Transcripts:",
        len(valid_df),
    )

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
# Evaluate transcript-level abundance
# ------------------------------------------------------------

def evaluate_transcript_abundance(
    profile_path,
    abundance_path,
    metadata_path,
    transcript_tpm_path,
    output_path,
):

    # --------------------------------------------------------
    # Read truth count
    # --------------------------------------------------------

    truth_df = read_simulation_truth(
        profile_path
    )

    # --------------------------------------------------------
    # Read simulator input TPM
    # --------------------------------------------------------

    abundance_df = read_simulation_abundance(
        abundance_path
    )

    truth_df = truth_df.merge(
        abundance_df,
        on="transcript_id",
        how="left",
    )

    # --------------------------------------------------------
    # Read RiboTxAlign results
    # --------------------------------------------------------

    estimated_df = read_ribotxalign_results(
        transcript_tpm_path=transcript_tpm_path,
        metadata_path=metadata_path,
    )

    print()
    print("=" * 70)
    print("Input summary")
    print("=" * 70)

    print(
        "Simulation truth transcripts:",
        len(truth_df),
    )

    print(
        "RiboTxAlign estimated transcripts:",
        len(estimated_df),
    )

    print(
        "Simulation total reads:",
        truth_df["truth_count"].sum(),
    )

    print(
        "RiboTxAlign expected reads:",
        estimated_df["ribotx_count"].sum(),
    )

    # --------------------------------------------------------
    # Merge truth and RiboTxAlign by ENST
    # --------------------------------------------------------

    df = truth_df.merge(
        estimated_df,
        on="transcript_id",
        how="left",
    )

    matched = (
        df["ribotx_count"]
        .notna()
    )

    print()
    print("=" * 70)
    print("Truth / RiboTxAlign matching")
    print("=" * 70)

    print(
        "Matched truth transcripts:",
        matched.sum(),
        "/",
        len(df),
        f"({matched.mean() * 100:.2f}%)",
    )

    print(
        "Truth reads in matched transcripts:",
        df.loc[
            matched,
            "truth_count",
        ].sum(),
    )

    print(
        "Truth reads in unmatched transcripts:",
        df.loc[
            ~matched,
            "truth_count",
        ].sum(),
    )

    # --------------------------------------------------------
    # Calculate truth TPM from Ribo-seq truth counts
    #
    # Use the same original CDS length used by RiboTxAlign.
    # --------------------------------------------------------

    df["truth_rpk_cds"] = np.nan

    valid_cds = (
        df["original_cds_length"].notna()
        &
        (
            df["original_cds_length"] > 0
        )
    )

    df.loc[
        valid_cds,
        "truth_rpk_cds",
    ] = (
        df.loc[
            valid_cds,
            "truth_count",
        ]
        /
        df.loc[
            valid_cds,
            "original_cds_length",
        ]
    )

    truth_rpk_sum = (
        df.loc[
            valid_cds,
            "truth_rpk_cds",
        ]
        .sum()
    )

    df["truth_tpm_cds"] = np.nan

    if truth_rpk_sum > 0:

        df.loc[
            valid_cds,
            "truth_tpm_cds",
        ] = (
            df.loc[
                valid_cds,
                "truth_rpk_cds",
            ]
            /
            truth_rpk_sum
            *
            1_000_000
        )

    # --------------------------------------------------------
    # Count accuracy
    # --------------------------------------------------------

    calculate_count_accuracy(
        df,
        "Transcript count accuracy: matched transcripts",
    )

    # --------------------------------------------------------
    # TPM correlation:
    # simulator input TPM
    # --------------------------------------------------------

    calculate_tpm_correlation(
        df,
        truth_column="truth_tpm_file",
        label=(
            "TPM correlation: "
            "simulator input TPM vs RiboTxAlign TPM"
        ),
    )

    # --------------------------------------------------------
    # TPM correlation:
    # TPM recalculated from Ribo-seq truth count + CDS length
    # --------------------------------------------------------

    calculate_tpm_correlation(
        df,
        truth_column="truth_tpm_cds",
        label=(
            "TPM correlation: "
            "truth count/CDS TPM vs RiboTxAlign TPM"
        ),
    )

    # --------------------------------------------------------
    # Single-transcript groups
    # --------------------------------------------------------

    single_df = df[
        matched
        &
        (
            df["num_transcripts"] == 1
        )
    ].copy()

    calculate_count_accuracy(
        single_df,
        "Single-transcript groups: count accuracy",
    )

    calculate_tpm_correlation(
        single_df,
        truth_column="truth_tpm_cds",
        label=(
            "Single-transcript groups: "
            "truth count/CDS TPM vs RiboTxAlign TPM"
        ),
    )

    # --------------------------------------------------------
    # Multi-transcript groups
    # --------------------------------------------------------

    multi_df = df[
        matched
        &
        (
            df["num_transcripts"] > 1
        )
    ].copy()

    calculate_count_accuracy(
        multi_df,
        "Multi-transcript groups: count accuracy",
    )

    calculate_tpm_correlation(
        multi_df,
        truth_column="truth_tpm_cds",
        label=(
            "Multi-transcript groups: "
            "truth count/CDS TPM vs RiboTxAlign TPM"
        ),
    )

    # --------------------------------------------------------
    # CDS length check
    # --------------------------------------------------------

    matched_df = df[
        matched
    ].copy()

    length_changed = matched_df[
        matched_df["original_cds_length"]
        !=
        matched_df["em_cds_length"]
    ]

    print()
    print("=" * 70)
    print("CDS length adjustment")
    print("=" * 70)

    print(
        "Transcripts with changed CDS length:",
        len(length_changed),
        "/",
        len(matched_df),
    )

    if len(matched_df) > 0:

        print(
            "Fraction:",
            f"{len(length_changed) / len(matched_df) * 100:.2f}%",
        )

    # --------------------------------------------------------
    # Errors
    # --------------------------------------------------------

    df["count_error"] = (
        df["ribotx_count"]
        -
        df["truth_count"]
    )

    df["absolute_count_error"] = (
        df["count_error"]
        .abs()
    )

    df["tpm_error_vs_cds_truth"] = (
        df["ribotx_tpm"]
        -
        df["truth_tpm_cds"]
    )

    # --------------------------------------------------------
    # Output table
    # --------------------------------------------------------

    output_columns = [
        "transcript_id",

        "truth_count",
        "ribotx_count",

        "original_cds_length",
        "em_cds_length",

        "simulation_length",
        "simulation_effective_length",

        "ribotx_tpm",
        "truth_tpm_file",
        "truth_tpm_cds",

        "num_transcripts",

        "count_error",
        "absolute_count_error",
        "tpm_error_vs_cds_truth",
    ]

    df[
        output_columns
    ].to_csv(
        output_path,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Show table
    # --------------------------------------------------------

    print()
    print("=" * 70)
    print("Comparison table preview")
    print("=" * 70)

    print(
        df[
            output_columns
        ].head(30).to_string(
            index=False
        )
    )

    # --------------------------------------------------------
    # Largest count errors
    # --------------------------------------------------------

    largest_errors = (
        df.sort_values(
            "absolute_count_error",
            ascending=False,
        )
        .head(30)
    )

    print()
    print("=" * 70)
    print("Largest transcript-level count errors")
    print("=" * 70)

    print(
        largest_errors[
            [
                "transcript_id",
                "truth_count",
                "ribotx_count",
                "original_cds_length",
                "ribotx_tpm",
                "truth_tpm_file",
                "truth_tpm_cds",
                "num_transcripts",
                "count_error",
            ]
        ].to_string(
            index=False
        )
    )

    print()
    print(
        "Comparison table written to:"
    )

    print(
        output_path
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

    abundance_path = (
        f"{base}/"
        "HEK293T.simulation_abundance.sf"
    )

    metadata_path = (
        f"{base}/"
        "RiboTxAlign/"
        "HEK293T_simulation/"
        "reference/"
        "HEK293T_simulation.metadata.tsv"
    )

    transcript_tpm_path = (
        f"{base}/"
        "RiboTxAlign/"
        "HEK293T_simulation/"
        "em/"
        "HEK293T_simulation.transcript_tpm.tsv"
    )

    output_path = (
        f"{base}/"
        "RiboTxAlign/"
        "HEK293T_simulation/"
        "em/"
        "transcript_abundance_comparison.tsv"
    )

    evaluate_transcript_abundance(
        profile_path=profile_path,
        abundance_path=abundance_path,
        metadata_path=metadata_path,
        transcript_tpm_path=transcript_tpm_path,
        output_path=output_path,
    )


compare()