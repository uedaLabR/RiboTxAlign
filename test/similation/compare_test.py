#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr


BASE = Path(
    "/mnt/share/project/riboML/"
    "validation_expression/similation/HEK293T"
)

TRUTH_PROFILE = BASE / "HEK293T.simulation.profile"
LONGREAD_TSV = BASE / "RiboTxAlign_LongRead.source_rescue_matching.tsv"

OUTPUT_DIR = BASE / "plots_longread_vs_truth"
OUTPUT_PNG = OUTPUT_DIR / "RiboTxAlign_LongRead_vs_Truth.png"
OUTPUT_PDF = OUTPUT_DIR / "RiboTxAlign_LongRead_vs_Truth.pdf"
OUTPUT_TSV = OUTPUT_DIR / "RiboTxAlign_LongRead_vs_Truth.tsv"


def norm_id(value):
    if value is None or pd.isna(value):
        return None

    value = str(value).strip()

    if not value:
        return None

    return value.split(".")[0]


def read_truth(path):
    records = []

    transcript_id = None
    profile_length = None

    with open(path, encoding="utf-8") as input_file:
        for line in input_file:
            line = line.strip()

            if line.startswith("tid:"):
                transcript_id = norm_id(
                    line.split(":", 1)[1]
                )

            elif line.startswith("plen:"):
                profile_length = int(
                    line.split(":", 1)[1]
                )

            elif line.startswith("rprofile:"):
                counts = [
                    int(x)
                    for x in line.split(":", 1)[1].split()
                ]

                records.append({
                    "transcript_id": transcript_id,
                    "truth_length": profile_length,
                    "truth_count": sum(counts),
                })

    df = (
        pd.DataFrame(records)
        .groupby(
            "transcript_id",
            as_index=False,
        )
        .agg(
            truth_length=("truth_length", "max"),
            truth_count=("truth_count", "sum"),
        )
    )

    rate = np.where(
        df["truth_length"] > 0,
        df["truth_count"] / df["truth_length"],
        0.0,
    )

    total = rate.sum()

    df["Truth_tpm"] = (
        rate / total * 1_000_000.0
        if total > 0
        else 0.0
    )

    return df


def read_longread(path):
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "transcript_id",
        "RiboTxAlign_LongRead_tpm",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            "Missing columns: "
            + ", ".join(sorted(missing))
        )

    df["transcript_id"] = (
        df["transcript_id"].map(norm_id)
    )

    df["RiboTxAlign_LongRead_tpm"] = (
        pd.to_numeric(
            df["RiboTxAlign_LongRead_tpm"],
            errors="coerce",
        )
        .fillna(0.0)
    )

    return df


def calc_corr(df):
    x = np.log10(
        df["Truth_tpm"].to_numpy(float) + 1.0
    )

    y = np.log10(
        df[
            "RiboTxAlign_LongRead_tpm"
        ].to_numpy(float) + 1.0
    )

    if len(df) < 2:
        return np.nan, np.nan

    return (
        pearsonr(x, y).statistic,
        spearmanr(x, y).statistic,
    )


def make_plot(df):
    plot_df = df[
        (df["Truth_tpm"] > 0)
        & (
            df["RiboTxAlign_LongRead_tpm"]
            > 0
        )
    ].copy()

    pearson, spearman = calc_corr(
        plot_df
    )

    x = np.log10(
        plot_df["Truth_tpm"].to_numpy(float)
        + 1.0
    )

    y = np.log10(
        plot_df[
            "RiboTxAlign_LongRead_tpm"
        ].to_numpy(float)
        + 1.0
    )

    fig, ax = plt.subplots(
        figsize=(8, 8)
    )

    ax.scatter(
        x,
        y,
        s=8,
        alpha=0.45,
    )

    max_value = max(
        x.max() if len(x) else 1.0,
        y.max() if len(y) else 1.0,
    )

    ax.plot(
        [0, max_value],
        [0, max_value],
        linestyle="--",
        linewidth=1.5,
    )

    ax.set_xlim(
        -0.05,
        max_value + 0.1,
    )
    ax.set_ylim(
        -0.05,
        max_value + 0.1,
    )

    ax.set_aspect(
        "equal",
        adjustable="box",
    )

    ax.set_xlabel(
        "log10(Simulated transcript TPM + 1)"
    )
    ax.set_ylabel(
        "log10(RiboTxAlign_LongRead transcript TPM + 1)"
    )

    ax.set_title(
        "RiboTxAlign LongRead vs Truth\n"
        "Truth > 0 and LongRead > 0"
    )

    ax.text(
        0.05,
        0.95,
        (
            f"n = {len(plot_df):,} transcripts\n"
            f"Pearson r = {pearson:.3f}\n"
            f"Spearman rho = {spearman:.3f}"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_PNG,
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        OUTPUT_PDF,
        bbox_inches="tight",
    )

    plt.close(fig)

    return plot_df, pearson, spearman


def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    truth = read_truth(
        TRUTH_PROFILE
    )
    longread = read_longread(
        LONGREAD_TSV
    )

    comparison = truth.merge(
        longread,
        on="transcript_id",
        how="left",
    )

    comparison[
        "RiboTxAlign_LongRead_tpm"
    ] = (
        comparison[
            "RiboTxAlign_LongRead_tpm"
        ]
        .fillna(0.0)
    )

    and_df = comparison[
        (comparison["Truth_tpm"] > 0)
        & (
            comparison[
                "RiboTxAlign_LongRead_tpm"
            ] > 0
        )
    ].copy()

    all_truth_df = comparison[
        comparison["Truth_tpm"] > 0
    ].copy()

    and_pearson, and_spearman = (
        calc_corr(and_df)
    )

    all_pearson, all_spearman = (
        calc_corr(all_truth_df)
    )

    make_plot(
        comparison
    )

    comparison[
        "log10_Truth_tpm_plus1"
    ] = np.log10(
        comparison["Truth_tpm"] + 1.0
    )

    comparison[
        "log10_RiboTxAlign_LongRead_tpm_plus1"
    ] = np.log10(
        comparison[
            "RiboTxAlign_LongRead_tpm"
        ] + 1.0
    )

    comparison.to_csv(
        OUTPUT_TSV,
        sep="\t",
        index=False,
    )

    print()
    print("=" * 72)
    print("RiboTxAlign LongRead vs Truth")
    print("=" * 72)

    print()
    print("AND: Truth > 0 and LongRead > 0")
    print(f"transcripts : {len(and_df):,}")
    print(f"Pearson r   : {and_pearson:.6f}")
    print(f"Spearman rho: {and_spearman:.6f}")

    print()
    print("ALL TRUTH-POSITIVE: LongRead zero included")
    print(f"transcripts : {len(all_truth_df):,}")
    print(f"Pearson r   : {all_pearson:.6f}")
    print(f"Spearman rho: {all_spearman:.6f}")

    print()
    print("Output:")
    print(OUTPUT_PNG)
    print(OUTPUT_PDF)
    print(OUTPUT_TSV)


if __name__ == "__main__":
    main()
