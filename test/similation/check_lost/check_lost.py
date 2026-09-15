#!/usr/bin/env python3

import re
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# Configuration
# ============================================================

BASE = Path(
    "/mnt/share/project/riboML/"
    "validation_expression/similation/HEK293T"
)

METHOD = "RiboTxAlign_LongRead"

PAIRWISE_TSV = (
    BASE
    / "plots_transcript"
    / f"transcript_pairwise_values.{METHOD}.tsv"
)

GENCODE_GFF3 = Path(
    "/mnt/share/project/riboML/validation/reference/"
    "gencode.v44.annotation.gff3"
)

OUTPUT_TSV = (
    BASE
    / "plots_transcript"
    / f"{METHOD}.outliers_for_IGV.tsv"
)

TOP_N = 10
IGV_PADDING = 500


# ============================================================
# Utility
# ============================================================

def norm_id(value):
    """Remove common prefixes and version suffixes."""

    if value is None or pd.isna(value):
        return None

    value = str(value).strip()

    if not value:
        return None

    value = re.sub(
        r"^(transcript:|gene:|rna:)",
        "",
        value,
    )

    return value.split(".")[0]


def parse_gff3_attributes(attribute_text):
    """Parse GFF3 attributes into a dictionary."""

    attributes = {}

    for item in attribute_text.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts

        attributes[key.strip()] = value.strip().strip('"')

    return attributes


# ============================================================
# GENCODE positions
# ============================================================

def read_gencode_positions(gff3_path):
    """Read transcript and gene genomic positions from GENCODE GFF3."""

    transcript_records = []
    gene_records = []

    with open(gff3_path, encoding="utf-8") as input_file:
        for line in input_file:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            (
                chromosome,
                source,
                feature,
                start,
                end,
                score,
                strand,
                phase,
                attribute_text,
            ) = fields

            attributes = parse_gff3_attributes(attribute_text)

            if feature == "gene":
                gene_id = norm_id(
                    attributes.get("gene_id")
                    or attributes.get("ID")
                )

                gene_name = (
                    attributes.get("gene_name")
                    or attributes.get("gene")
                    or attributes.get("Name")
                    or ""
                )

                if gene_id is None:
                    continue

                gene_records.append({
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "gene_chromosome": chromosome,
                    "gene_start": int(start),
                    "gene_end": int(end),
                    "gene_strand": strand,
                })

            elif feature in {"transcript", "mRNA"}:
                transcript_id = norm_id(
                    attributes.get("transcript_id")
                    or attributes.get("ID")
                )

                gene_id = norm_id(
                    attributes.get("gene_id")
                    or attributes.get("Parent")
                )

                gene_name = (
                    attributes.get("gene_name")
                    or attributes.get("gene")
                    or ""
                )

                if transcript_id is None:
                    continue

                transcript_records.append({
                    "transcript_id": transcript_id,
                    "gene_id": gene_id,
                    "gene_name_transcript": gene_name,
                    "transcript_chromosome": chromosome,
                    "transcript_start": int(start),
                    "transcript_end": int(end),
                    "transcript_strand": strand,
                })

    transcript_df = (
        pd.DataFrame(transcript_records)
        .drop_duplicates(subset=["transcript_id"])
    )

    gene_df = (
        pd.DataFrame(gene_records)
        .drop_duplicates(subset=["gene_id"])
    )

    result = transcript_df.merge(
        gene_df,
        on="gene_id",
        how="left",
    )

    result["gene_name"] = result["gene_name"].fillna(
        result["gene_name_transcript"]
    )

    result["gene_chromosome"] = result["gene_chromosome"].fillna(
        result["transcript_chromosome"]
    )

    result["gene_start"] = result["gene_start"].fillna(
        result["transcript_start"]
    )

    result["gene_end"] = result["gene_end"].fillna(
        result["transcript_end"]
    )

    result["gene_strand"] = result["gene_strand"].fillna(
        result["transcript_strand"]
    )

    return result


# ============================================================
# Pairwise values and error calculation
# ============================================================

def read_pairwise_table(path, method):
    """Read transcript-level truth-versus-method table."""

    df = pd.read_csv(path, sep="\t")
    method_tpm_column = f"{method}_tpm"

    required = {
        "transcript_id",
        "Truth_tpm",
        method_tpm_column,
    }

    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            "Missing columns: "
            + ", ".join(sorted(missing))
        )

    df["transcript_id"] = df["transcript_id"].map(norm_id)

    return df


def calculate_outlier_scores(pairwise, method):
    """Calculate signed log10 error for transcripts positive in both sets."""

    method_tpm_column = f"{method}_tpm"

    df = pairwise[
        (pairwise["Truth_tpm"] > 0)
        & (pairwise[method_tpm_column] > 0)
    ].copy()

    df["log10_truth"] = np.log10(df["Truth_tpm"] + 1.0)
    df["log10_method"] = np.log10(df[method_tpm_column] + 1.0)

    df["log_error"] = (
        df["log10_method"]
        - df["log10_truth"]
    )

    df["abs_log_error"] = df["log_error"].abs()
    df["fold_ratio"] = df[method_tpm_column] / df["Truth_tpm"]

    return df


def add_genomic_positions(dataframe, annotation):
    """Attach transcript/gene genomic coordinates and IGV locus."""

    result = dataframe.merge(
        annotation,
        on="transcript_id",
        how="left",
    )

    result["IGV_gene_start"] = (
        pd.to_numeric(
            result["gene_start"],
            errors="coerce",
        )
        - IGV_PADDING
    ).clip(lower=1)

    result["IGV_gene_end"] = (
        pd.to_numeric(
            result["gene_end"],
            errors="coerce",
        )
        + IGV_PADDING
    )

    result["IGV_locus"] = ""

    valid = (
        result["gene_chromosome"].notna()
        & result["IGV_gene_start"].notna()
        & result["IGV_gene_end"].notna()
    )

    result.loc[valid, "IGV_locus"] = (
        result.loc[valid, "gene_chromosome"].astype(str)
        + ":"
        + result.loc[valid, "IGV_gene_start"].astype(int).astype(str)
        + "-"
        + result.loc[valid, "IGV_gene_end"].astype(int).astype(str)
    )

    return result


# ============================================================
# Outlier selection
# ============================================================

def select_overestimated(dataframe, top_n):
    return (
        dataframe[dataframe["log_error"] > 0]
        .sort_values("log_error", ascending=False)
        .head(top_n)
        .copy()
    )


def select_underestimated(dataframe, top_n):
    return (
        dataframe[dataframe["log_error"] < 0]
        .sort_values("log_error", ascending=True)
        .head(top_n)
        .copy()
    )


def select_largest_errors(dataframe, top_n):
    return (
        dataframe
        .sort_values("abs_log_error", ascending=False)
        .head(top_n)
        .copy()
    )


# ============================================================
# Printing
# ============================================================

def print_outlier_table(title, dataframe, method):
    """Print an IGV-friendly summary table."""

    method_tpm_column = f"{method}_tpm"
    merge_status_column = f"{method}_merge_status"

    print()
    print("=" * 120)
    print(title)
    print("=" * 120)

    if dataframe.empty:
        print("No transcripts found.")
        return

    columns = [
        "transcript_id",
        "gene_name",
        "gene_id",
        "Truth_tpm",
        method_tpm_column,
        "log_error",
        "fold_ratio",
        merge_status_column,
        "transcript_chromosome",
        "transcript_start",
        "transcript_end",
        "transcript_strand",
        "gene_chromosome",
        "gene_start",
        "gene_end",
        "gene_strand",
        "IGV_locus",
    ]

    columns = [
        column
        for column in columns
        if column in dataframe.columns
    ]

    display = dataframe[columns].copy()

    for column in [
        "Truth_tpm",
        method_tpm_column,
        "log_error",
        "fold_ratio",
    ]:
        if column in display.columns:
            display[column] = display[column].map(
                lambda value: f"{value:.4f}"
                if pd.notna(value)
                else ""
            )

    print(display.to_string(index=False))


def print_igv_shortlist(over, under, method):
    """Print five over- and five under-estimated loci for manual IGV review."""

    method_tpm_column = f"{method}_tpm"

    print()
    print("=" * 120)
    print("IGV SHORTLIST")
    print("=" * 120)

    selected = pd.concat(
        [
            over.head(5).assign(direction="OVER"),
            under.head(5).assign(direction="UNDER"),
        ],
        ignore_index=True,
    )

    for _, row in selected.iterrows():

        transcript_chr = str(
            row.get(
                "transcript_chromosome",
                "",
            )
        )

        transcript_start = row.get(
            "transcript_start",
            np.nan,
        )

        transcript_end = row.get(
            "transcript_end",
            np.nan,
        )

        gene_chr = str(
            row.get(
                "gene_chromosome",
                "",
            )
        )

        gene_start = row.get(
            "gene_start",
            np.nan,
        )

        gene_end = row.get(
            "gene_end",
            np.nan,
        )

        transcript_locus = (
            f"{transcript_chr}:"
            f"{int(transcript_start)}-"
            f"{int(transcript_end)}"
            if (
                transcript_chr
                and pd.notna(transcript_start)
                and pd.notna(transcript_end)
            )
            else ""
        )

        gene_locus = (
            f"{gene_chr}:"
            f"{int(gene_start)}-"
            f"{int(gene_end)}"
            if (
                gene_chr
                and pd.notna(gene_start)
                and pd.notna(gene_end)
            )
            else ""
        )

        print(
            f"{row['direction']:5s}  "
            f"{str(row.get('transcript_id', '')):18s}  "
            f"{str(row.get('gene_name', '')):16s}  "
            f"Truth={float(row['Truth_tpm']):10.3f}  "
            f"RiboTx={float(row[method_tpm_column]):10.3f}  "
            f"logErr={float(row['log_error']):+7.3f}"
        )

        print(
            f"       transcript locus : {transcript_locus}"
        )

        print(
            f"       gene locus       : {gene_locus}"
        )

        print(
            f"       IGV locus        : "
            f"{str(row.get('IGV_locus', ''))}"
        )

        print()


# ============================================================
# Main
# ============================================================

def main():

    print("Pairwise TSV :", PAIRWISE_TSV)
    print("GENCODE GFF3 :", GENCODE_GFF3)

    if not PAIRWISE_TSV.exists():
        raise FileNotFoundError(PAIRWISE_TSV)

    if not GENCODE_GFF3.exists():
        raise FileNotFoundError(GENCODE_GFF3)

    pairwise = read_pairwise_table(
        PAIRWISE_TSV,
        METHOD,
    )

    annotation = read_gencode_positions(
        GENCODE_GFF3
    )

    scored = calculate_outlier_scores(
        pairwise,
        METHOD,
    )

    scored = add_genomic_positions(
        scored,
        annotation,
    )

    method_tpm_column = f"{METHOD}_tpm"
    merge_status_column = f"{METHOD}_merge_status"

    over = select_overestimated(
        scored,
        TOP_N,
    )

    under = select_underestimated(
        scored,
        TOP_N,
    )

    print()
    print("Transcripts used:", len(scored))

    print_outlier_table(
        "TOP OVER-ESTIMATED TRANSCRIPTS: RiboTxAlign > Truth",
        over,
        METHOD,
    )

    print_outlier_table(
        "TOP UNDER-ESTIMATED TRANSCRIPTS: RiboTxAlign < Truth",
        under,
        METHOD,
    )

    print_igv_shortlist(
        over,
        under,
        METHOD,
    )

    # Never-merged outliers are useful for separating mapping/EM errors
    # from errors caused by transcript collapsing.
    if merge_status_column in scored.columns:
        never_merged = scored[
            scored[merge_status_column]
            .fillna("never_merged")
            .eq("never_merged")
        ].copy()

        never_merged = select_largest_errors(
            never_merged,
            TOP_N,
        )

        print_outlier_table(
            "LARGEST ERRORS AMONG NEVER-MERGED TRANSCRIPTS",
            never_merged,
            METHOD,
        )

    preferred_columns = [
        "transcript_id",
        "gene_name",
        "gene_id",
        "Truth_tpm",
        method_tpm_column,
        "log10_truth",
        "log10_method",
        "log_error",
        "abs_log_error",
        "fold_ratio",
        merge_status_column,
        f"{METHOD}_match_type",
        f"{METHOD}_reference_merged",
        f"{METHOD}_em_merged",
        "transcript_chromosome",
        "transcript_start",
        "transcript_end",
        "transcript_strand",
        "gene_chromosome",
        "gene_start",
        "gene_end",
        "gene_strand",
        "IGV_locus",
    ]

    output_columns = [
        column
        for column in preferred_columns
        if column in scored.columns
    ]

    remaining_columns = [
        column
        for column in scored.columns
        if column not in output_columns
    ]

    output = scored[
        output_columns + remaining_columns
    ].sort_values(
        "abs_log_error",
        ascending=False,
    )

    output.to_csv(
        OUTPUT_TSV,
        sep="\t",
        index=False,
    )

    print()
    print("=" * 120)
    print("OUTPUT")
    print("=" * 120)
    print(OUTPUT_TSV)

    print()
    print(
        "Load the RiboTxAlign genome BAM and BAI in IGV, "
        "then paste IGV_locus into the IGV search box."
    )


if __name__ == "__main__":
    main()