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

    df = pd.DataFrame(records)

    df["transcript_id"] = (
        df["transcript_id"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    df = (
        df.groupby(
            "transcript_id",
            as_index=False,
        )["truth_count"]
        .sum()
    )

    return df


# ------------------------------------------------------------
# Read simulator input TPM
# ------------------------------------------------------------

def read_simulation_abundance(
    abundance_path,
):
    df = pd.read_csv(
        abundance_path,
        sep="\t",
    )

    df = df.rename(
        columns={
            "Name": "transcript_id",
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
            "truth_tpm_file",
        ]
    ].copy()


# ------------------------------------------------------------
# Read RiboTxAlign representative output
# ------------------------------------------------------------

def read_ribotxalign_representative(
    representative_tpm_path,
):
    df = pd.read_csv(
        representative_tpm_path,
        sep="\t",
    )

    # Remove version for ENST IDs
    df["transcript_id"] = (
        df["transcript_id"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    df = df.rename(
        columns={
            "tpm": "ribotx_tpm",
            "expected_count": "ribotx_count",
        }
    )

    return df


# ------------------------------------------------------------
# Read CDS length from detailed RiboTxAlign output
# ------------------------------------------------------------

def read_ribotxalign_cds_length(
    transcript_tpm_path,
    metadata_path,
):
    result_df = pd.read_csv(
        transcript_tpm_path,
        sep="\t",
    )

    metadata_df = pd.read_csv(
        metadata_path,
        sep="\t",
    )

    metadata_df["transcript_id"] = (
        metadata_df["representative_transcript_id"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    df = result_df.merge(
        metadata_df[
            [
                "reference_id",
                "transcript_id",
            ]
        ],
        left_on="transcript_id",
        right_on="reference_id",
        how="left",
    )

    df = (
        df.groupby(
            "transcript_id_y",
            as_index=False,
        )
        .agg(
            cds_length=(
                "em_cds_length",
                "max",
            )
        )
        .rename(
            columns={
                "transcript_id_y": "transcript_id",
            }
        )
    )

    return df


# ------------------------------------------------------------
# Calculate correlation
# ------------------------------------------------------------

def report_correlation(
    df,
    truth_column,
    estimate_column,
    label,
):
    valid = df[
        df[truth_column].notna()
        &
        df[estimate_column].notna()
    ].copy()

    truth = valid[
        truth_column
    ].to_numpy(
        dtype=float
    )

    estimate = valid[
        estimate_column
    ].to_numpy(
        dtype=float
    )

    pearson = pearsonr(
        truth,
        estimate,
    ).statistic

    spearman = spearmanr(
        truth,
        estimate,
    ).statistic

    mae = mean_absolute_error(
        truth,
        estimate,
    )

    rmse = np.sqrt(
        mean_squared_error(
            truth,
            estimate,
        )
    )

    print()
    print("=" * 70)
    print(label)
    print("=" * 70)

    print(
        "Transcripts:",
        len(valid),
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
# Main evaluation
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

    transcript_tpm_path = (
        f"{base}/"
        "RiboTxAlign/"
        "HEK293T_simulation/"
        "em/"
        "HEK293T_simulation.transcript_tpm.tsv"
    )

    representative_tpm_path = (
        f"{base}/"
        "RiboTxAlign/"
        "HEK293T_simulation/"
        "em/"
        "transcript_tpm_representative.tsv"
    )

    metadata_path = (
        f"{base}/"
        "RiboTxAlign/"
        "HEK293T_simulation/"
        "reference/"
        "HEK293T_simulation.metadata.tsv"
    )

    output_path = (
        f"{base}/"
        "RiboTxAlign/"
        "HEK293T_simulation/"
        "em/"
        "representative_tpm_comparison.tsv"
    )

    # --------------------------------------------------------
    # Read data
    # --------------------------------------------------------

    truth_df = read_simulation_truth(
        profile_path
    )

    abundance_df = read_simulation_abundance(
        abundance_path
    )

    ribotx_df = read_ribotxalign_representative(
        representative_tpm_path
    )

    cds_df = read_ribotxalign_cds_length(
        transcript_tpm_path=transcript_tpm_path,
        metadata_path=metadata_path,
    )

    # --------------------------------------------------------
    # Merge
    # --------------------------------------------------------

    df = truth_df.merge(
        abundance_df,
        on="transcript_id",
        how="left",
    )

    df = df.merge(
        ribotx_df[
            [
                "transcript_id",
                "ribotx_tpm",
                "ribotx_count",
                "source",
                "reference_id",
            ]
        ],
        on="transcript_id",
        how="left",
    )

    df = df.merge(
        cds_df,
        on="transcript_id",
        how="left",
    )

    # --------------------------------------------------------
    # Calculate truth TPM from footprint count and CDS length
    # --------------------------------------------------------

    df["truth_rpk_cds"] = np.nan

    valid_length = (
        df["cds_length"].notna()
        &
        (
            df["cds_length"] > 0
        )
    )

    df.loc[
        valid_length,
        "truth_rpk_cds",
    ] = (
        df.loc[
            valid_length,
            "truth_count",
        ]
        /
        df.loc[
            valid_length,
            "cds_length",
        ]
    )

    total_rpk = (
        df.loc[
            valid_length,
            "truth_rpk_cds",
        ]
        .sum()
    )

    df["truth_tpm_cds"] = np.nan

    df.loc[
        valid_length,
        "truth_tpm_cds",
    ] = (
        df.loc[
            valid_length,
            "truth_rpk_cds",
        ]
        /
        total_rpk
        *
        1_000_000
    )

    # --------------------------------------------------------
    # Matching summary
    # --------------------------------------------------------

    matched = (
        df["ribotx_tpm"]
        .notna()
    )

    print()
    print("=" * 70)
    print("Matching")
    print("=" * 70)

    print(
        "Truth transcripts:",
        len(df),
    )

    print(
        "Matched transcripts:",
        matched.sum(),
        "/",
        len(df),
        f"({matched.mean() * 100:.2f}%)",
    )

    # --------------------------------------------------------
    # Count correlation
    # --------------------------------------------------------

    report_correlation(
        df,
        truth_column="truth_count",
        estimate_column="ribotx_count",
        label="Truth count vs RiboTxAlign expected count",
    )

    # --------------------------------------------------------
    # TPM from simulator input
    # --------------------------------------------------------

    report_correlation(
        df,
        truth_column="truth_tpm_file",
        estimate_column="ribotx_tpm",
        label="Simulator input TPM vs RiboTxAlign TPM",
    )

    # --------------------------------------------------------
    # TPM calculated from simulated footprint counts
    # --------------------------------------------------------

    report_correlation(
        df,
        truth_column="truth_tpm_cds",
        estimate_column="ribotx_tpm",
        label="Footprint count/CDS truth TPM vs RiboTxAlign TPM",
    )

    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    output_columns = [
        "transcript_id",
        "truth_count",
        "ribotx_count",
        "cds_length",
        "truth_tpm_file",
        "truth_tpm_cds",
        "ribotx_tpm",
        "source",
        "reference_id",
    ]

    df[
        output_columns
    ].to_csv(
        output_path,
        sep="\t",
        index=False,
    )

    print()
    print(
        "Written:",
        output_path,
    )


compare()