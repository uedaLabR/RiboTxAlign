from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuration
# ============================================================

BASE = Path(
    "/mnt/share/project/riboML/validation"
)

COMPARISON_ROOT = (
    BASE
    / "bam_comparison_four_mappers"
)

OUTPUT_DIR = (
    COMPARISON_ROOT
    / "summary"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

EXCEL_PATH = (
    OUTPUT_DIR
    / "mapping_comparison_summary.xlsx"
)


# ============================================================
# Dataset configuration
# ============================================================

DATASETS = {
    "HEK293T": {
        "samples": [
            "SRR26942051",
            "SRR26942052",
            "SRR26942053",
            "SRR26942054",
        ],
    },
    "mouse_brain": {
        "samples": [
            "SRR2163083",
            "SRR2163093",
            "SRR2163103",
        ],
    },
}


# ============================================================
# Mapping methods
# ============================================================

MAPPERS = [
    "STAR_genome",
    "STAR_transcript",
    "RiboTxAlign_longread_stringtie",
    "RiboTxAlign_shortread_stringtie",
]


MAPPER_LABELS = {
    "STAR_genome": (
        "STAR\nGenome"
    ),

    "STAR_transcript": (
        "STAR\nTranscript"
    ),

    "RiboTxAlign_longread_stringtie": (
        "RiboTxAlign\nLong-read StringTie"
    ),

    "RiboTxAlign_shortread_stringtie": (
        "RiboTxAlign\nShort-read StringTie"
    ),
}


NH_MAX = 10


# ============================================================
# Utility
# ============================================================

def select_all_scope(
    df,
    input_path,
):
    """
    Use all mapper reads.

    common_to_all is intentionally not used.
    """

    if "scope" not in df.columns:
        return df.copy()

    scopes = (
        df["scope"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    target_scope = "all_mapper_reads"

    if target_scope not in scopes:
        raise RuntimeError(
            f"'{target_scope}' scope is not found in:\n"
            f"{input_path}\n"
            f"Available scopes: {scopes}"
        )

    return (
        df[
            df["scope"] == target_scope
        ]
        .copy()
    )

def safe_percent(
    numerator,
    denominator,
):

    if denominator <= 0:
        return np.nan

    return (
        numerator
        / denominator
        * 100.0
    )


# ============================================================
# Load raw mapping summaries
# ============================================================

def load_raw_mapping_summary():

    rows = []

    for (
        dataset,
        config,
    ) in DATASETS.items():

        dataset_dir = (
            COMPARISON_ROOT
            / dataset
        )

        for sample in config["samples"]:

            input_path = (
                dataset_dir
                / sample
                / "mapping_summary.tsv"
            )

            if not input_path.exists():

                print(
                    f"WARNING: not found: "
                    f"{input_path}"
                )

                continue

            df = pd.read_csv(
                input_path,
                sep="\t",
            )

            df = select_all_scope(
                df=df,
                input_path=input_path,
            )

            df["dataset"] = dataset
            df["sample"] = sample

            rows.append(
                df
            )

    if not rows:

        raise RuntimeError(
            "No mapping_summary.tsv files found."
        )

    return pd.concat(
        rows,
        ignore_index=True,
    )


# ============================================================
# Recalculate mapping statistics
# ============================================================

def recalculate_mapping_statistics(
    df,
):
    """
    For each sample:

        true total reads
            = maximum reported total_reads among all mappers

    Then:

        mapped_reads
            = unique_reads + multi_reads

        unmapped_reads
            = total_reads - mapped_reads

    All percentages use the same sample-level total.
    """

    required_columns = [
        "total_reads",
        "unique_reads",
        "multi_reads",
    ]

    for column in required_columns:

        if column not in df.columns:

            raise RuntimeError(
                f"Required column is missing: "
                f"{column}"
            )

    df = df.copy()

    # --------------------------------------------------------
    # Keep mapper-reported total for reference.
    # --------------------------------------------------------

    df[
        "reported_total_reads"
    ] = (
        df[
            "total_reads"
        ]
    )

    # --------------------------------------------------------
    # Determine true total number of reads for each sample.
    #
    # STAR BAM may contain mapped reads only.
    # Therefore use the maximum total observed among methods.
    # --------------------------------------------------------

    sample_totals = (
        df
        .groupby(
            [
                "dataset",
                "sample",
            ]
        )[
            "reported_total_reads"
        ]
        .max()
        .rename(
            "sample_total_reads"
        )
        .reset_index()
    )

    df = df.merge(
        sample_totals,
        on=[
            "dataset",
            "sample",
        ],
        how="left",
    )

    df[
        "total_reads"
    ] = (
        df[
            "sample_total_reads"
        ]
        .astype(int)
    )

    # --------------------------------------------------------
    # Standardized definitions
    # --------------------------------------------------------

    df[
        "mapped_reads"
    ] = (
        df[
            "unique_reads"
        ]
        + df[
            "multi_reads"
        ]
    )

    df[
        "unmapped_reads"
    ] = (
        df[
            "total_reads"
        ]
        - df[
            "mapped_reads"
        ]
    )

    # --------------------------------------------------------
    # Check for impossible values
    # --------------------------------------------------------

    negative_mask = (
        df[
            "unmapped_reads"
        ]
        < 0
    )

    if negative_mask.any():

        print()
        print(
            "WARNING:"
            " mapped reads exceed sample total."
        )

        print(
            df.loc[
                negative_mask,
                [
                    "dataset",
                    "sample",
                    "mapper",
                    "reported_total_reads",
                    "total_reads",
                    "unique_reads",
                    "multi_reads",
                    "mapped_reads",
                ],
            ].to_string(
                index=False
            )
        )

    # --------------------------------------------------------
    # Percent of total reads
    # --------------------------------------------------------

    df[
        "unique_percent_of_total"
    ] = (
        df[
            "unique_reads"
        ]
        / df[
            "total_reads"
        ]
        * 100.0
    )

    df[
        "multi_percent_of_total"
    ] = (
        df[
            "multi_reads"
        ]
        / df[
            "total_reads"
        ]
        * 100.0
    )

    df[
        "mapped_percent_of_total"
    ] = (
        df[
            "mapped_reads"
        ]
        / df[
            "total_reads"
        ]
        * 100.0
    )

    df[
        "unmapped_percent_of_total"
    ] = (
        df[
            "unmapped_reads"
        ]
        / df[
            "total_reads"
        ]
        * 100.0
    )

    # --------------------------------------------------------
    # Percent of mapped reads
    # --------------------------------------------------------

    df[
        "unique_percent_of_mapped"
    ] = np.where(
        df[
            "mapped_reads"
        ]
        > 0,

        df[
            "unique_reads"
        ]
        / df[
            "mapped_reads"
        ]
        * 100.0,

        np.nan,
    )

    df[
        "multi_percent_of_mapped"
    ] = np.where(
        df[
            "mapped_reads"
        ]
        > 0,

        df[
            "multi_reads"
        ]
        / df[
            "mapped_reads"
        ]
        * 100.0,

        np.nan,
    )

    # --------------------------------------------------------
    # Spliced reads
    # --------------------------------------------------------

    if "splice_reads" in df.columns:

        df[
            "splice_percent_of_total"
        ] = (
            df[
                "splice_reads"
            ]
            / df[
                "total_reads"
            ]
            * 100.0
        )

        df[
            "splice_percent_of_mapped"
        ] = np.where(
            df[
                "mapped_reads"
            ]
            > 0,

            df[
                "splice_reads"
            ]
            / df[
                "mapped_reads"
            ]
            * 100.0,

            np.nan,
        )

    return df


# ============================================================
# Load mapping summary
# ============================================================

def load_mapping_summary():

    raw_df = (
        load_raw_mapping_summary()
    )

    result = (
        recalculate_mapping_statistics(
            raw_df
        )
    )

    return result


# ============================================================
# Load NH distribution
# ============================================================

def load_nh_distribution(
    mapping_summary,
):

    rows = []

    for (
        dataset,
        config,
    ) in DATASETS.items():

        dataset_dir = (
            COMPARISON_ROOT
            / dataset
        )

        for sample in config["samples"]:

            nh_path = (
                dataset_dir
                / sample
                / "multi_reads_by_nh.tsv"
            )

            if not nh_path.exists():

                print(
                    f"WARNING: not found: "
                    f"{nh_path}"
                )

                continue

            multi_df = pd.read_csv(
                nh_path,
                sep="\t",
            )

            multi_df = select_all_scope(
                df=multi_df,
                input_path=nh_path,
            )

            for mapper in MAPPERS:

                mapper_summary = (
                    mapping_summary[
                        (
                            mapping_summary[
                                "dataset"
                            ]
                            == dataset
                        )
                        & (
                            mapping_summary[
                                "sample"
                            ]
                            == sample
                        )
                        & (
                            mapping_summary[
                                "mapper"
                            ]
                            == mapper
                        )
                    ]
                )

                if mapper_summary.empty:
                    continue

                total_reads = int(
                    mapper_summary.iloc[
                        0
                    ][
                        "total_reads"
                    ]
                )

                unique_reads = int(
                    mapper_summary.iloc[
                        0
                    ][
                        "unique_reads"
                    ]
                )

                # ------------------------------------------------
                # NH = 1
                # ------------------------------------------------

                rows.append(
                    {
                        "dataset": dataset,
                        "sample": sample,
                        "mapper": mapper,
                        "nh_group": 1,
                        "read_count": (
                            unique_reads
                        ),
                        "percent_of_total": (
                            safe_percent(
                                unique_reads,
                                total_reads,
                            )
                        ),
                    }
                )

                # ------------------------------------------------
                # NH >= 2
                # ------------------------------------------------

                mapper_multi = (
                    multi_df[
                        multi_df[
                            "mapper"
                        ]
                        == mapper
                    ]
                )

                for (
                    _,
                    row,
                ) in mapper_multi.iterrows():

                    nh = int(
                        row[
                            "nh"
                        ]
                    )

                    count = int(
                        row[
                            "read_count"
                        ]
                    )

                    if nh >= NH_MAX:

                        nh_group = (
                            NH_MAX
                        )

                    else:

                        nh_group = (
                            nh
                        )

                    rows.append(
                        {
                            "dataset": dataset,
                            "sample": sample,
                            "mapper": mapper,
                            "nh_group": (
                                nh_group
                            ),
                            "read_count": (
                                count
                            ),
                            "percent_of_total": (
                                safe_percent(
                                    count,
                                    total_reads,
                                )
                            ),
                        }
                    )

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "No NH distribution data found."
        )

    # --------------------------------------------------------
    # Combine NH >= NH_MAX
    # --------------------------------------------------------

    result = (
        result
        .groupby(
            [
                "dataset",
                "sample",
                "mapper",
                "nh_group",
            ],
            as_index=False,
        )
        .agg(
            read_count=(
                "read_count",
                "sum",
            ),

            percent_of_total=(
                "percent_of_total",
                "sum",
            ),
        )
    )

    return result


# ============================================================
# Graph 1
#
# Mapping composition
#
# Mean stacked bars:
#
#   Unique
#   Multi-hit
#   Unmapped
#
# Individual replicate points are plotted at:
#
#   1. Unique boundary
#   2. Mapped boundary (= Unique + Multi)
#
# Error bars represent mean +/- SD of those boundaries.
# ============================================================

def plot_mapping_composition(
    mapping_summary,
    dataset,
):

    df = (
        mapping_summary[
            mapping_summary[
                "dataset"
            ]
            == dataset
        ]
        .copy()
    )

    # --------------------------------------------------------
    # Calculate cumulative mapped boundary.
    # --------------------------------------------------------

    df[
        "mapped_boundary"
    ] = (
        df[
            "unique_percent_of_total"
        ]
        + df[
            "multi_percent_of_total"
        ]
    )

    # --------------------------------------------------------
    # Mean values
    # --------------------------------------------------------

    summary = (
        df
        .groupby(
            "mapper"
        )
        .agg(
            unique_mean=(
                "unique_percent_of_total",
                "mean",
            ),

            unique_sd=(
                "unique_percent_of_total",
                "std",
            ),

            multi_mean=(
                "multi_percent_of_total",
                "mean",
            ),

            unmapped_mean=(
                "unmapped_percent_of_total",
                "mean",
            ),

            mapped_mean=(
                "mapped_boundary",
                "mean",
            ),

            mapped_sd=(
                "mapped_boundary",
                "std",
            ),
        )
        .reindex(
            MAPPERS
        )
    )

    summary = (
        summary
        .fillna(
            0.0
        )
    )

    unique_mean = (
        summary[
            "unique_mean"
        ]
        .to_numpy()
    )

    multi_mean = (
        summary[
            "multi_mean"
        ]
        .to_numpy()
    )

    unmapped_mean = (
        summary[
            "unmapped_mean"
        ]
        .to_numpy()
    )

    unique_sd = (
        summary[
            "unique_sd"
        ]
        .to_numpy()
    )

    mapped_mean = (
        summary[
            "mapped_mean"
        ]
        .to_numpy()
    )

    mapped_sd = (
        summary[
            "mapped_sd"
        ]
        .to_numpy()
    )

    x = np.arange(
        len(
            MAPPERS
        )
    )

    fig, ax = plt.subplots(
        figsize=(
            8.5,
            5.5,
        )
    )

    # --------------------------------------------------------
    # Mean stacked bars
    # --------------------------------------------------------

    ax.bar(
        x,
        unique_mean,
        width=0.62,
        label="Unique",
    )

    ax.bar(
        x,
        multi_mean,
        width=0.62,
        bottom=unique_mean,
        label="Multi-hit",
    )

    ax.bar(
        x,
        unmapped_mean,
        width=0.62,
        bottom=(
            unique_mean
            + multi_mean
        ),
        label="Unmapped",
    )

    # --------------------------------------------------------
    # Error bars
    #
    # Boundary 1:
    #   Unique
    #
    # Boundary 2:
    #   Unique + Multi = mapped
    # --------------------------------------------------------

    ax.errorbar(
        x,
        unique_mean,
        yerr=unique_sd,
        fmt="none",
        capsize=4,
        linewidth=1.2,
        zorder=6,
    )

    ax.errorbar(
        x,
        mapped_mean,
        yerr=mapped_sd,
        fmt="none",
        capsize=4,
        linewidth=1.2,
        zorder=6,
    )

    # --------------------------------------------------------
    # Individual replicate points
    #
    # Points represent the actual segment boundaries.
    # --------------------------------------------------------

    samples = (
        DATASETS[
            dataset
        ][
            "samples"
        ]
    )

    offsets = np.linspace(
        -0.16,
        0.16,
        len(
            samples
        ),
    )

    for (
        mapper_index,
        mapper,
    ) in enumerate(
        MAPPERS
    ):

        mapper_df = (
            df[
                df[
                    "mapper"
                ]
                == mapper
            ]
        )

        for (
            sample_index,
            sample,
        ) in enumerate(
            samples
        ):

            sample_df = (
                mapper_df[
                    mapper_df[
                        "sample"
                    ]
                    == sample
                ]
            )

            if sample_df.empty:
                continue

            row = (
                sample_df.iloc[
                    0
                ]
            )

            unique_boundary = float(
                row[
                    "unique_percent_of_total"
                ]
            )

            mapped_boundary = float(
                row[
                    "mapped_boundary"
                ]
            )

            point_x = (
                mapper_index
                + offsets[
                    sample_index
                ]
            )

            # Unique boundary
            ax.scatter(
                point_x,
                unique_boundary,
                s=26,
                marker="o",
                facecolors="none",
                edgecolors="black",
                linewidths=1.0,
                zorder=7,
            )

            # Unique + Multi boundary
            ax.scatter(
                point_x,
                mapped_boundary,
                s=30,
                marker="^",
                facecolors="none",
                edgecolors="black",
                linewidths=1.0,
                zorder=7,
            )

    # --------------------------------------------------------
    # Axis
    # --------------------------------------------------------

    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        [
            MAPPER_LABELS[
                mapper
            ]
            for mapper
            in MAPPERS
        ],
        rotation=20,
        ha="right",
    )

    ax.set_ylabel(
        "Reads (% of total)"
    )

    ax.set_ylim(
        0,
        100,
    )

    ax.set_title(
        dataset
    )

    ax.legend(
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(
            1.02,
            1.0,
        ),
    )

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

    plt.tight_layout()

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    output_png = (
        OUTPUT_DIR
        / (
            f"{dataset}_"
            "mapping_composition.png"
        )
    )

    output_pdf = (
        OUTPUT_DIR
        / (
            f"{dataset}_"
            "mapping_composition.pdf"
        )
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

    plt.close(
        fig
    )

    print(
        f"Saved: "
        f"{output_png}"
    )


# ============================================================
# Graph 2
#
# NH distribution
# ============================================================

def plot_nh_distribution(
    nh_df,
    dataset,
):

    df = (
        nh_df[
            nh_df[
                "dataset"
            ]
            == dataset
        ]
        .copy()
    )

    fig, axes = plt.subplots(
        1,
        4,
        figsize=(
            17,
            4.4,
        ),
        sharey=True,
    )

    nh_values = list(
        range(
            1,
            NH_MAX + 1,
        )
    )

    x = np.arange(
        len(
            nh_values
        )
    )

    for (
        ax,
        mapper,
    ) in zip(
        axes,
        MAPPERS,
    ):

        mapper_df = (
            df[
                df[
                    "mapper"
                ]
                == mapper
            ]
        )

        summary = (
            mapper_df
            .groupby(
                "nh_group"
            )[
                "percent_of_total"
            ]
            .agg(
                mean="mean",
                sd="std",
            )
        )

        means = []
        sds = []

        for nh in nh_values:

            if nh in summary.index:

                mean_value = (
                    summary.loc[
                        nh,
                        "mean",
                    ]
                )

                sd_value = (
                    summary.loc[
                        nh,
                        "sd",
                    ]
                )

                if pd.isna(
                    sd_value
                ):
                    sd_value = 0.0

                means.append(
                    mean_value
                )

                sds.append(
                    sd_value
                )

            else:

                means.append(
                    0.0
                )

                sds.append(
                    0.0
                )

        ax.bar(
            x,
            means,
            yerr=sds,
            capsize=3,
            width=0.7,
        )

        # ----------------------------------------------------
        # Individual replicate points
        # ----------------------------------------------------

        samples = (
            DATASETS[
                dataset
            ][
                "samples"
            ]
        )

        offsets = np.linspace(
            -0.16,
            0.16,
            len(
                samples
            ),
        )

        for (
            nh_index,
            nh,
        ) in enumerate(
            nh_values
        ):

            for (
                sample_index,
                sample,
            ) in enumerate(
                samples
            ):

                values = (
                    mapper_df[
                        (
                            mapper_df[
                                "nh_group"
                            ]
                            == nh
                        )
                        & (
                            mapper_df[
                                "sample"
                            ]
                            == sample
                        )
                    ][
                        "percent_of_total"
                    ]
                )

                if values.empty:
                    continue

                value = float(
                    values.iloc[
                        0
                    ]
                )

                ax.scatter(
                    nh_index
                    + offsets[
                        sample_index
                    ],
                    value,
                    s=20,
                    facecolors="none",
                    edgecolors="black",
                    zorder=5,
                )

        labels = [
            str(
                nh
            )
            for nh in range(
                1,
                NH_MAX,
            )
        ]

        labels.append(
            f">={NH_MAX}"
        )

        ax.set_xticks(
            x
        )

        ax.set_xticklabels(
            labels
        )

        ax.set_xlabel(
            "Number of mapping loci (NH)"
        )

        ax.set_title(
            MAPPER_LABELS[
                mapper
            ]
        )

        ax.set_ylim(
            bottom=0
        )

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

    axes[
        0
    ].set_ylabel(
        "Reads (% of total)"
    )

    fig.suptitle(
        dataset,
        y=1.02,
    )

    plt.tight_layout()

    output_png = (
        OUTPUT_DIR
        / (
            f"{dataset}_"
            "NH_distribution.png"
        )
    )

    output_pdf = (
        OUTPUT_DIR
        / (
            f"{dataset}_"
            "NH_distribution.pdf"
        )
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

    plt.close(
        fig
    )

    print(
        f"Saved: "
        f"{output_png}"
    )


# ============================================================
# Mapping counts table
# ============================================================

def make_mapping_counts_table(
    mapping_summary,
):

    columns = [
        "dataset",
        "sample",
        "mapper",

        "reported_total_reads",
        "total_reads",

        "mapped_reads",
        "unmapped_reads",

        "unique_reads",
        "multi_reads",

        "splice_reads",
    ]

    columns = [
        column
        for column in columns
        if column
        in mapping_summary.columns
    ]

    return (
        mapping_summary[
            columns
        ]
        .copy()
    )


# ============================================================
# Mapping percentage table
# ============================================================

def make_mapping_percent_table(
    mapping_summary,
):

    columns = [
        "dataset",
        "sample",
        "mapper",

        "unique_percent_of_total",
        "multi_percent_of_total",
        "unmapped_percent_of_total",
        "mapped_percent_of_total",

        "unique_percent_of_mapped",
        "multi_percent_of_mapped",

        "splice_percent_of_total",
        "splice_percent_of_mapped",
    ]

    columns = [
        column
        for column in columns
        if column
        in mapping_summary.columns
    ]

    return (
        mapping_summary[
            columns
        ]
        .copy()
    )


# ============================================================
# NH counts table
# ============================================================

def make_nh_counts_table(
    nh_df,
):

    table = (
        nh_df
        .pivot_table(
            index=[
                "dataset",
                "sample",
                "mapper",
            ],

            columns=(
                "nh_group"
            ),

            values=(
                "read_count"
            ),

            aggfunc="sum",

            fill_value=0,
        )
        .reset_index()
    )

    rename_columns = {}

    for column in table.columns:

        if isinstance(
            column,
            (int, np.integer),
        ):

            if column == NH_MAX:

                rename_columns[
                    column
                ] = (
                    f"NH>={NH_MAX}"
                )

            else:

                rename_columns[
                    column
                ] = (
                    f"NH{column}"
                )

    return table.rename(
        columns=rename_columns
    )


# ============================================================
# NH percentage table
# ============================================================

def make_nh_percent_table(
    nh_df,
):

    table = (
        nh_df
        .pivot_table(
            index=[
                "dataset",
                "sample",
                "mapper",
            ],

            columns=(
                "nh_group"
            ),

            values=(
                "percent_of_total"
            ),

            aggfunc="sum",

            fill_value=0.0,
        )
        .reset_index()
    )

    rename_columns = {}

    for column in table.columns:

        if isinstance(
            column,
            (int, np.integer),
        ):

            if column == NH_MAX:

                rename_columns[
                    column
                ] = (
                    f"NH>={NH_MAX}"
                )

            else:

                rename_columns[
                    column
                ] = (
                    f"NH{column}"
                )

    return table.rename(
        columns=rename_columns
    )


# ============================================================
# Save Excel
# ============================================================

def save_excel(
    mapping_summary,
    nh_df,
):

    mapping_counts = (
        make_mapping_counts_table(
            mapping_summary
        )
    )

    mapping_percent = (
        make_mapping_percent_table(
            mapping_summary
        )
    )

    nh_counts = (
        make_nh_counts_table(
            nh_df
        )
    )

    nh_percent = (
        make_nh_percent_table(
            nh_df
        )
    )

    with pd.ExcelWriter(
        EXCEL_PATH,
        engine="openpyxl",
    ) as writer:

        mapping_counts.to_excel(
            writer,
            sheet_name="Mapping_counts",
            index=False,
        )

        mapping_percent.to_excel(
            writer,
            sheet_name="Mapping_percent",
            index=False,
        )

        nh_counts.to_excel(
            writer,
            sheet_name="NH_counts",
            index=False,
        )

        nh_percent.to_excel(
            writer,
            sheet_name="NH_percent",
            index=False,
        )

    print(
        f"Saved Excel: "
        f"{EXCEL_PATH}"
    )


# ============================================================
# Save combined TSV
# ============================================================

def save_tables(
    mapping_summary,
    nh_df,
):

    mapping_summary.to_csv(
        (
            OUTPUT_DIR
            / "mapping_summary_all_datasets.tsv"
        ),
        sep="\t",
        index=False,
    )

    nh_df.to_csv(
        (
            OUTPUT_DIR
            / "NH_distribution_all_datasets.tsv"
        ),
        sep="\t",
        index=False,
    )


# ============================================================
# Consistency check
# ============================================================

def print_consistency_check(
    mapping_summary,
):

    print()
    print(
        "========================================"
    )

    print(
        "Consistency check"
    )

    print(
        "========================================"
    )

    for (
        dataset,
        sample,
    ), sample_df in (
        mapping_summary
        .groupby(
            [
                "dataset",
                "sample",
            ]
        )
    ):

        total_reads = int(
            sample_df[
                "total_reads"
            ]
            .iloc[
                0
            ]
        )

        print()
        print(
            f"{dataset} / {sample}"
        )

        print(
            f"Total reads: "
            f"{total_reads:,}"
        )

        for _, row in sample_df.iterrows():

            unique = int(
                row[
                    "unique_reads"
                ]
            )

            multi = int(
                row[
                    "multi_reads"
                ]
            )

            unmapped = int(
                row[
                    "unmapped_reads"
                ]
            )

            total_check = (
                unique
                + multi
                + unmapped
            )

            if (
                total_check
                == total_reads
            ):

                status = "OK"

            else:

                status = "WARNING"

            print(
                f"  {row['mapper']}: "
                f"unique={unique:,}, "
                f"multi={multi:,}, "
                f"unmapped={unmapped:,}, "
                f"sum={total_check:,} "
                f"[{status}]"
            )


# ============================================================
# Main
# ============================================================

def main():

    # ========================================================
    # Load
    # ========================================================

    mapping_summary = (
        load_mapping_summary()
    )

    nh_df = (
        load_nh_distribution(
            mapping_summary
        )
    )

    # ========================================================
    # Print mapping summary
    # ========================================================

    print()
    print(
        "========================================"
    )

    print(
        "Mapping summary"
    )

    print(
        "========================================"
    )

    columns_to_show = [
        "dataset",
        "sample",
        "mapper",

        "reported_total_reads",
        "total_reads",

        "unique_reads",
        "multi_reads",
        "unmapped_reads",

        "unique_percent_of_total",
        "multi_percent_of_total",
        "unmapped_percent_of_total",
    ]

    columns_to_show = [
        column
        for column in columns_to_show
        if column
        in mapping_summary.columns
    ]

    print(
        mapping_summary[
            columns_to_show
        ].to_string(
            index=False
        )
    )

    # ========================================================
    # Check totals
    # ========================================================

    print_consistency_check(
        mapping_summary
    )

    # ========================================================
    # Plots
    # ========================================================

    for dataset in DATASETS:

        plot_mapping_composition(
            mapping_summary=(
                mapping_summary
            ),
            dataset=(
                dataset
            ),
        )

        plot_nh_distribution(
            nh_df=(
                nh_df
            ),
            dataset=(
                dataset
            ),
        )

    # ========================================================
    # Save tables
    # ========================================================

    save_tables(
        mapping_summary=(
            mapping_summary
        ),
        nh_df=(
            nh_df
        ),
    )

    # ========================================================
    # Save Excel
    # ========================================================

    save_excel(
        mapping_summary=(
            mapping_summary
        ),
        nh_df=(
            nh_df
        ),
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
        f"Excel: "
        f"{EXCEL_PATH}"
    )


if __name__ == "__main__":
    main()