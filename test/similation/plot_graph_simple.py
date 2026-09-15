from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages


# ============================================================
# Configuration
# ============================================================

BASE = Path(
    "/mnt/share/project/riboML/validation_expression/similation"
)

SUMMARY_FILE = (
    BASE
    / "simulation_mapping_accuracy_summary.tsv"
)

RIBOTX_NH_FILE = (
    BASE
    / "ribotx_multihit_accuracy_by_nh.tsv"
)

STAR_MULTI_NH_FILE = (
    BASE
    / "star_multi_ribotx_accuracy_by_nh.tsv"
)

OUTPUT_DIR = (
    BASE
    / "plots"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

COMBINED_PDF = (
    OUTPUT_DIR
    / "simulation_mapping_all_figures.pdf"
)


DATASET_ORDER = [
    "HEK293T",
    "mouse_brain",
]

MAPPER_ORDER = [
    "STAR",
    "RiboTxAlign",
]


# ============================================================
# Load data
# ============================================================

def load_summary():

    if not SUMMARY_FILE.exists():
        raise FileNotFoundError(
            SUMMARY_FILE
        )

    return pd.read_csv(
        SUMMARY_FILE,
        sep="\t",
    )


def load_nh_table(
    path,
):

    if not path.exists():
        raise FileNotFoundError(
            path
        )

    return pd.read_csv(
        path,
        sep="\t",
    )


# ============================================================
# Utility
# ============================================================

def get_metric(
    df,
    dataset,
    mapper,
    metric,
    column="percent_of_evaluated_reads",
):

    selected = df[
        (
            df["dataset"]
            == dataset
        )
        & (
            df["mapper"]
            == mapper
        )
        & (
            df["metric"]
            == metric
        )
    ]

    if selected.empty:
        return 0.0

    value = selected.iloc[
        0
    ][
        column
    ]

    if pd.isna(
        value
    ):
        return 0.0

    return float(
        value
    )


def save_figure(
    fig,
    output_stem,
    pdf,
):

    output_png = (
        OUTPUT_DIR
        / f"{output_stem}.png"
    )

    output_pdf = (
        OUTPUT_DIR
        / f"{output_stem}.pdf"
    )

    fig.savefig(
        output_png,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        output_pdf,
        bbox_inches="tight",
    )

    if pdf is not None:

        pdf.savefig(
            fig,
            bbox_inches="tight",
        )

    plt.close(
        fig
    )

    print(
        f"Saved: {output_png}"
    )


def clean_axis(
    ax,
):

    ax.spines[
        "top"
    ].set_visible(
        False
    )

    ax.spines[
        "right"
    ].set_visible(
        False
    )


# ============================================================
# Graph 1
#
# Overall mapping result
#
# Correct:
#   unique hit
#   multi hit
#
# Error:
#   splice error
#   alignment error
#   unmapped
# ============================================================

def plot_overall_accuracy(
    summary_df,
    pdf=None,
):

    metrics = [
        (
            "correct_unique",
            "Correct\nunique",
        ),
        (
            "correct_multi",
            "Correct\nmulti",
        ),
        (
            "splice_error",
            "Splice\nerror",
        ),
        (
            "alignment_error",
            "Alignment\nerror",
        ),
        (
            "unmapped",
            "Unmapped",
        ),
    ]

    for dataset in DATASET_ORDER:

        fig, ax = plt.subplots(
            figsize=(
                8.5,
                4.8,
            )
        )

        x = np.arange(
            len(
                metrics
            )
        )

        width = 0.36

        for (
            mapper_index,
            mapper,
        ) in enumerate(
            MAPPER_ORDER
        ):

            values = []

            for (
                metric,
                _,
            ) in metrics:

                value = get_metric(
                    df=summary_df,
                    dataset=dataset,
                    mapper=mapper,
                    metric=metric,
                )

                values.append(
                    value
                )

            offset = (
                mapper_index
                - (
                    len(
                        MAPPER_ORDER
                    )
                    - 1
                )
                / 2
            ) * width

            ax.bar(
                x + offset,
                values,
                width=width,
                label=mapper,
            )

        ax.set_xticks(
            x
        )

        ax.set_xticklabels(
            [
                label
                for (
                    _,
                    label,
                )
                in metrics
            ]
        )

        ax.set_ylabel(
            "Reads (%)"
        )

        ax.set_title(
            dataset
        )

        ax.set_ylim(
            bottom=0
        )

        ax.legend(
            frameon=False
        )

        clean_axis(
            ax
        )

        plt.tight_layout()

        save_figure(
            fig=fig,
            output_stem=(
                f"{dataset}_"
                "overall_mapping_accuracy"
            ),
            pdf=pdf,
        )


# ============================================================
# Graph 2
#
# Correct mapping comparison
#
# STAR vs RiboTxAlign
#
# Unique and multi are shown separately.
# ============================================================

def plot_correct_mapping_comparison(
    summary_df,
    pdf=None,
):

    categories = []
    unique_values = []
    multi_values = []

    for dataset in DATASET_ORDER:
        for mapper in MAPPER_ORDER:

            categories.append(
                f"{dataset}\n{mapper}"
            )

            unique_values.append(
                get_metric(
                    df=summary_df,
                    dataset=dataset,
                    mapper=mapper,
                    metric="correct_unique",
                )
            )

            multi_values.append(
                get_metric(
                    df=summary_df,
                    dataset=dataset,
                    mapper=mapper,
                    metric="correct_multi",
                )
            )

    fig, ax = plt.subplots(
        figsize=(8.5, 4.8)
    )

    x = np.arange(
        len(categories)
    )

    width = 0.6

    # Correct unique
    ax.bar(
        x,
        unique_values,
        width=width,
        label="Unique hit",
    )

    # Correct multi stacked on unique
    ax.bar(
        x,
        multi_values,
        width=width,
        bottom=unique_values,
        label="Multi hit",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(categories)

    ax.set_ylabel(
        "Correctly mapped reads (%)"
    )

    ax.set_ylim(
        0,
        100,
    )

    ax.legend(
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0,
    )

    clean_axis(ax)

    plt.tight_layout(
        rect=[0, 0, 0.82, 1]
    )

    save_figure(
        fig=fig,
        output_stem="correct_mapping_comparison",
        pdf=pdf,
    )

# ============================================================
# Graph 3
#
# RiboTxAlign:
# correct assignment accuracy by NH
#
# Observed accuracy vs random expectation
# ============================================================

def plot_ribotx_nh_accuracy(
    nh_df,
    pdf=None,
):

    for dataset in DATASET_ORDER:

        df = nh_df[
            nh_df[
                "dataset"
            ]
            == dataset
        ].copy()

        if df.empty:
            continue

        df = df.sort_values(
            "nh"
        )

        fig, ax = plt.subplots(
            figsize=(
                6.5,
                4.8,
            )
        )

        ax.plot(
            df[
                "nh"
            ],
            df[
                "observed_accuracy_percent"
            ],
            marker="o",
            linewidth=1.5,
            label="RiboTxAlign",
        )

        ax.plot(
            df[
                "nh"
            ],
            df[
                "random_expectation_percent"
            ],
            marker="o",
            linewidth=1.5,
            linestyle="--",
            label="Random (1/N)",
        )

        ax.set_xlabel(
            "Number of candidate mappings (NH)"
        )

        ax.set_ylabel(
            "Correct assignment (%)"
        )

        ax.set_title(
            dataset
        )

        ax.set_ylim(
            0,
            100,
        )

        ax.set_xticks(
            df[
                "nh"
            ]
        )

        ax.legend(
            frameon=False
        )

        clean_axis(
            ax
        )

        plt.tight_layout()

        save_figure(
            fig=fig,
            output_stem=(
                f"{dataset}_"
                "ribotx_accuracy_by_NH"
            ),
            pdf=pdf,
        )


# ============================================================
# Graph 4
#
# RiboTxAlign:
# number of reads by NH
#
# This is the denominator used in Graph 3.
# ============================================================

def plot_ribotx_nh_counts(
    nh_df,
    pdf=None,
):

    for dataset in DATASET_ORDER:

        df = nh_df[
            nh_df[
                "dataset"
            ]
            == dataset
        ].copy()

        if df.empty:
            continue

        df = df.sort_values(
            "nh"
        )

        fig, ax = plt.subplots(
            figsize=(
                6.5,
                4.8,
            )
        )

        ax.bar(
            df[
                "nh"
            ],
            df[
                "total"
            ],
            width=0.7,
        )

        ax.set_xlabel(
            "Number of candidate mappings (NH)"
        )

        ax.set_ylabel(
            "Number of reads"
        )

        ax.set_title(
            dataset
        )

        ax.set_xticks(
            df[
                "nh"
            ]
        )

        ax.set_ylim(
            bottom=0
        )

        clean_axis(
            ax
        )

        plt.tight_layout()

        save_figure(
            fig=fig,
            output_stem=(
                f"{dataset}_"
                "ribotx_read_count_by_NH"
            ),
            pdf=pdf,
        )


# ============================================================
# Graph 5
#
# STAR multi-mapped reads:
#
# RiboTxAlign correct assignment accuracy
# by RiboTxAlign NH
# ============================================================

def plot_star_multi_nh_accuracy(
    nh_df,
    pdf=None,
):

    for dataset in DATASET_ORDER:

        df = nh_df[
            nh_df[
                "dataset"
            ]
            == dataset
        ].copy()

        if df.empty:
            continue

        df = df.sort_values(
            "ribotx_nh"
        )

        fig, ax = plt.subplots(
            figsize=(
                6.5,
                4.8,
            )
        )

        ax.plot(
            df[
                "ribotx_nh"
            ],
            df[
                "observed_accuracy_percent"
            ],
            marker="o",
            linewidth=1.5,
            label="RiboTxAlign",
        )

        ax.plot(
            df[
                "ribotx_nh"
            ],
            df[
                "random_expectation_percent"
            ],
            marker="o",
            linewidth=1.5,
            linestyle="--",
            label="Random (1/N)",
        )

        ax.set_xlabel(
            "Number of RiboTxAlign candidates (NH)"
        )

        ax.set_ylabel(
            "Correct assignment (%)"
        )

        ax.set_title(
            f"{dataset}: STAR multi-mapped reads"
        )

        ax.set_ylim(
            0,
            100,
        )

        ax.set_xticks(
            df[
                "ribotx_nh"
            ]
        )

        ax.legend(
            frameon=False
        )

        clean_axis(
            ax
        )

        plt.tight_layout()

        save_figure(
            fig=fig,
            output_stem=(
                f"{dataset}_"
                "STAR_multi_RiboTx_accuracy_by_NH"
            ),
            pdf=pdf,
        )


# ============================================================
# Graph 6
#
# STAR multi-mapped reads:
#
# Number of reads in each RiboTxAlign NH group.
# ============================================================

def plot_star_multi_nh_counts(
    nh_df,
    pdf=None,
):

    for dataset in DATASET_ORDER:

        df = nh_df[
            nh_df[
                "dataset"
            ]
            == dataset
        ].copy()

        if df.empty:
            continue

        df = df.sort_values(
            "ribotx_nh"
        )

        fig, ax = plt.subplots(
            figsize=(
                6.5,
                4.8,
            )
        )

        ax.bar(
            df[
                "ribotx_nh"
            ],
            df[
                "total"
            ],
            width=0.7,
        )

        ax.set_xlabel(
            "Number of RiboTxAlign candidates (NH)"
        )

        ax.set_ylabel(
            "Number of reads"
        )

        ax.set_title(
            f"{dataset}: STAR multi-mapped reads"
        )

        ax.set_xticks(
            df[
                "ribotx_nh"
            ]
        )

        ax.set_ylim(
            bottom=0
        )

        clean_axis(
            ax
        )

        plt.tight_layout()

        save_figure(
            fig=fig,
            output_stem=(
                f"{dataset}_"
                "STAR_multi_RiboTx_read_count_by_NH"
            ),
            pdf=pdf,
        )


# ============================================================
# Save combined plotting data
# ============================================================

def save_plot_tables(
    summary_df,
    ribotx_nh_df,
    star_multi_nh_df,
):

    summary_df.to_csv(
        (
            OUTPUT_DIR
            / "plot_mapping_summary.tsv"
        ),
        sep="\t",
        index=False,
    )

    ribotx_nh_df.to_csv(
        (
            OUTPUT_DIR
            / "plot_ribotx_accuracy_by_NH.tsv"
        ),
        sep="\t",
        index=False,
    )

    star_multi_nh_df.to_csv(
        (
            OUTPUT_DIR
            / "plot_STAR_multi_RiboTx_accuracy_by_NH.tsv"
        ),
        sep="\t",
        index=False,
    )


# ============================================================
# Main
# ============================================================

def main():

    # --------------------------------------------------------
    # Load
    # --------------------------------------------------------

    summary_df = (
        load_summary()
    )

    ribotx_nh_df = (
        load_nh_table(
            RIBOTX_NH_FILE
        )
    )

    star_multi_nh_df = (
        load_nh_table(
            STAR_MULTI_NH_FILE
        )
    )

    # --------------------------------------------------------
    # Check expected metrics
    # --------------------------------------------------------

    required_metrics = {
        "correct_unique",
        "correct_multi",
    }

    available_metrics = set(
        summary_df[
            "metric"
        ]
    )

    missing_metrics = (
        required_metrics
        - available_metrics
    )

    if missing_metrics:

        raise RuntimeError(
            "The summary file does not contain: "
            + ", ".join(
                sorted(
                    missing_metrics
                )
            )
            + "\n"
            + "Please rerun the simulation evaluation "
            + "after adding correct_unique / correct_multi."
        )

    # --------------------------------------------------------
    # Save tables
    # --------------------------------------------------------

    save_plot_tables(
        summary_df=summary_df,
        ribotx_nh_df=ribotx_nh_df,
        star_multi_nh_df=star_multi_nh_df,
    )

    # --------------------------------------------------------
    # Create all figures and combined PDF
    # --------------------------------------------------------

    with PdfPages(
        COMBINED_PDF
    ) as pdf:

        # Graph 1
        plot_overall_accuracy(
            summary_df=summary_df,
            pdf=pdf,
        )

        # Graph 2
        plot_correct_mapping_comparison(
            summary_df=summary_df,
            pdf=pdf,
        )

        # Graph 3
        plot_ribotx_nh_accuracy(
            nh_df=ribotx_nh_df,
            pdf=pdf,
        )

        # Graph 4
        plot_ribotx_nh_counts(
            nh_df=ribotx_nh_df,
            pdf=pdf,
        )

        # Graph 5
        plot_star_multi_nh_accuracy(
            nh_df=star_multi_nh_df,
            pdf=pdf,
        )

        # Graph 6
        plot_star_multi_nh_counts(
            nh_df=star_multi_nh_df,
            pdf=pdf,
        )

    print()
    print(
        "========================================"
    )

    print(
        "Finished"
    )

    print(
        "========================================"
    )

    print(
        f"Output directory: "
        f"{OUTPUT_DIR}"
    )

    print(
        f"Combined PDF: "
        f"{COMBINED_PDF}"
    )


if __name__ == "__main__":
    main()