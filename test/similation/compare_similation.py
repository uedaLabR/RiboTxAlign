#!/usr/bin/env python3

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error


# ============================================================
# Configuration
# ============================================================

BASE = Path(
    "/mnt/share/project/riboML/"
    "validation_expression/similation/HEK293T"
)

ANNOTATION = Path(
    "/mnt/share/project/riboML/validation/reference/"
    "gencode.v44.annotation.gff3"
)

RIBOMAP_STATS = Path(
    "/mnt/share/project/riboML/"
    "validation_expression/riboMap/HEK293T/ribomap.stats"
)

DEEPSHAPE_LOG = Path(
    "/mnt/share/project/riboML/"
    "validation_expression/deepShape/HEK293T/"
    "DeepShapePrimeOutputs/runlog_199.txt"
)

# Each directory is searched recursively for exactly one
# *transcript_tpm_with_source.tsv file.
RIBOTX_METHOD_DIRS = {
    "RiboTxAlign_LongRead": BASE / "RiboTxAlign_LongRead",
    "RiboTxAlign_RNASeq": BASE / "RiboTxAlign_RNASeq",
}

# If the old RiboTxAlign result should also be included,
# uncomment the following line:
#
# RIBOTX_METHOD_DIRS["RiboTxAlign_Old"] = BASE / "RiboTxAlign"


# ============================================================
# Utility
# ============================================================

def norm_id(x):
    if x is None or pd.isna(x):
        return None

    x = str(x).strip()

    if not x:
        return None

    x = re.sub(
        r"^(transcript:|gene:|rna:)",
        "",
        x,
    )

    return x.split(".")[0]


def split_ids(x):
    if x is None or pd.isna(x):
        return []

    return [
        norm_id(v)
        for v in re.split(r"[,;]", str(x))
        if norm_id(v)
    ]


def calc_tpm(
    df,
    count_col,
    length_col,
    output_col,
):
    df = df.copy()

    rate = np.where(
        df[length_col].fillna(0) > 0,
        df[count_col] / df[length_col],
        0.0,
    )

    total = rate.sum()

    df[output_col] = (
        rate / total * 1_000_000
        if total > 0
        else 0.0
    )

    return df


def accuracy(
    df,
    column,
    method,
    subset,
):
    y = df["Truth_tpm"].to_numpy(float)
    p = df[column].to_numpy(float)

    if len(df) < 2:
        return {
            "subset": subset,
            "method": method,
            "transcripts": len(df),
            "pearson": np.nan,
            "spearman": np.nan,
            "mae": np.nan,
            "rmse": np.nan,
        }

    return {
        "subset": subset,
        "method": method,
        "transcripts": len(df),
        "pearson": pearsonr(y, p).statistic,
        "spearman": spearmanr(y, p).statistic,
        "mae": mean_absolute_error(y, p),
        "rmse": np.sqrt(mean_squared_error(y, p)),
    }


def find_ribotx_result(directory):
    """
    Find one RiboTxAlign transcript TPM output below a method directory.

    This accepts both naming styles, for example:

      em/transcript_tpm_with_source.tsv

    and

      em/HEK293T_simulation_LongRead.transcript_tpm_with_source.tsv
    """

    hits = sorted(
        directory.rglob("*transcript_tpm_with_source.tsv")
    )

    if len(hits) == 0:
        raise FileNotFoundError(
            f"No *transcript_tpm_with_source.tsv found under: "
            f"{directory}"
        )

    if len(hits) > 1:
        raise RuntimeError(
            "Multiple transcript TPM files found under "
            f"{directory}:\n"
            + "\n".join(str(x) for x in hits)
        )

    return hits[0]


def find_ribotx_metadata(directory):
    """
    Find the reference metadata table used before EM transcript merging.
    """

    # Diagnostic exports such as diagnosis_*/reference_metadata.tsv are
    # subsets, not the reference table used by the pipeline.
    candidates = sorted(directory.rglob("*.metadata.tsv"))
    reference_hits = [path for path in candidates if path.parent.name == "reference"]
    hits = reference_hits or candidates

    if len(hits) == 0:
        raise FileNotFoundError(
            f"No pipeline *.metadata.tsv found under: {directory}"
        )

    if len(hits) > 1:
        raise RuntimeError(
            "Multiple metadata files found under "
            f"{directory}:\n"
            + "\n".join(str(x) for x in hits)
        )

    return hits[0]


def find_ribotx_merge(directory):
    """
    Find the EM transcript merge table.
    """

    hits = sorted(
        directory.rglob("*transcript_merge.tsv")
    )

    if len(hits) == 0:
        raise FileNotFoundError(
            f"No *transcript_merge.tsv found under: {directory}"
        )

    if len(hits) > 1:
        raise RuntimeError(
            "Multiple transcript merge files found under "
            f"{directory}:\n"
            + "\n".join(str(x) for x in hits)
        )

    return hits[0]


# ============================================================
# Truth
# ============================================================

def read_truth(path):
    records = []

    tid = None
    plen = None

    with open(path) as f:
        for line in f:
            line = line.strip()

            if line.startswith("tid:"):
                tid = norm_id(
                    line.split(":", 1)[1]
                )

            elif line.startswith("plen:"):
                plen = int(
                    line.split(":", 1)[1]
                )

            elif line.startswith("rprofile:"):
                counts = [
                    int(x)
                    for x in (
                        line.split(":", 1)[1]
                        .split()
                    )
                ]

                records.append({
                    "transcript_id": tid,
                    "truth_length": plen,
                    "truth_count": sum(counts),
                })

    df = (
        pd.DataFrame(records)
        .groupby(
            "transcript_id",
            as_index=False,
        )
        .agg(
            truth_length=(
                "truth_length",
                "max",
            ),
            truth_count=(
                "truth_count",
                "sum",
            ),
        )
    )

    return calc_tpm(
        df,
        "truth_count",
        "truth_length",
        "Truth_tpm",
    )


def read_rna(path):
    df = pd.read_csv(
        path,
        sep="\t",
    )

    df = df.rename(
        columns={
            "Name": "transcript_id",
            "TPM": "RNA_input_tpm",
        }
    )

    df["transcript_id"] = (
        df["transcript_id"]
        .map(norm_id)
    )

    return (
        df.groupby(
            "transcript_id",
            as_index=False,
        )["RNA_input_tpm"]
        .sum()
    )


# ============================================================
# Annotation
# ============================================================

def parse_attributes(text):
    result = {}

    for item in text.split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" in item:
            key, value = item.split(
                "=",
                1,
            )
        else:
            parts = item.split(
                None,
                1,
            )

            if len(parts) != 2:
                continue

            key, value = parts

        result[key.strip()] = (
            value.strip()
            .strip('"')
        )

    return result


def read_transcript_gene(path):
    records = []

    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = (
                line.rstrip()
                .split("\t")
            )

            if len(fields) < 9:
                continue

            feature = fields[2]
            attr = parse_attributes(
                fields[8]
            )

            transcript_id = attr.get(
                "transcript_id"
            )

            gene_id = attr.get(
                "gene_id"
            )

            if (
                transcript_id is None
                and feature in {
                    "transcript",
                    "mRNA",
                    "lnc_RNA",
                    "ncRNA",
                }
            ):
                transcript_id = attr.get("ID")

            if (
                gene_id is None
                and feature in {
                    "transcript",
                    "mRNA",
                    "lnc_RNA",
                    "ncRNA",
                }
            ):
                gene_id = attr.get("Parent")

            transcript_id = norm_id(
                transcript_id
            )

            gene_id = norm_id(
                gene_id
            )

            if (
                transcript_id
                and gene_id
                and transcript_id.startswith("ENST")
                and gene_id.startswith("ENSG")
            ):
                records.append({
                    "transcript_id":
                        transcript_id,
                    "truth_gene_id":
                        gene_id,
                })

    df = (
        pd.DataFrame(records)
        .drop_duplicates(
            "transcript_id"
        )
    )

    print()
    print("=" * 72)
    print("ANNOTATION")
    print("=" * 72)
    print(
        "Transcript-gene mappings:",
        len(df),
    )
    print(
        "Genes:",
        df["truth_gene_id"].nunique(),
    )

    return df


# ============================================================
# RiboTxAlign
# ============================================================

def read_ribotx(path):
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    for col in [
        "tpm",
        "expected_count",
    ]:
        df[col] = (
            pd.to_numeric(
                df[col],
                errors="coerce",
            )
            .fillna(0.0)
        )

    for col in [
        "source_reference_id",
        "representative_transcript_id",
        "original_transcript_id",
        "reference_gene_id",
    ]:
        if col not in df.columns:
            df[col] = ""

        df[col + "_norm"] = (
            df[col]
            .map(norm_id)
        )

    # Preserve the source column: a GENCODE fallback is an explicit native
    # identity, not a StringTie source or an annotation-alias rescue.
    df["evaluation_transcript_id_norm"] = df.apply(
        evaluation_transcript_id, axis=1,
    )
    return df


def evaluation_transcript_id(row):
    source = norm_id(row.get("source_reference_id"))
    if source:
        return source
    # Both reference-construction paths carry an explicit GENCODE representative.
    # Resolve only that representative; do not expand the same-CDS member list
    # into multiple copies of the final abundance.
    if row.get("resolution_type") not in {
        "gencode_representative_uncovered",
        "gencode_representative_same_cds",
    }:
        return None
    original = norm_id(row.get("original_transcript_id"))
    representative = norm_id(row.get("representative_transcript_id"))
    if original and representative and original != representative:
        raise ValueError(
            f"Conflicting GENCODE identity for {row.get('reference_id')}: "
            f"{original} versus {representative}"
        )
    if original and original == representative and original.startswith("ENST"):
        return original
    return None


def read_ribotx_metadata(path):
    """
    Read the pre-EM reference metadata.

    This table is required because a source reference can disappear from
    the final transcript TPM table after EM transcript merging.
    """

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    for col in [
        "reference_id",
        "source_reference_id",
        "representative_transcript_id",
        "original_transcript_id",
        "gencode_transcript_ids",
        "source_transcript_ids",
        "overlap_group_id",
        "reference_gene_id",
        "resolution_type",
    ]:
        if col not in df.columns:
            df[col] = ""

    for col in [
        "source_reference_id",
        "representative_transcript_id",
        "original_transcript_id",
        "reference_gene_id",
    ]:
        df[col + "_norm"] = (
            df[col]
            .map(norm_id)
        )

    return df


def read_ribotx_merge(path):
    """
    Read source-reference -> target-reference merges performed before /
    during transcript-level EM.
    """

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "source_transcript_id",
        "target_transcript_id",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"{path} is missing columns: {sorted(missing)}"
        )

    return df


def build_merge_map(merge_df):
    """
    Return source_reference -> target_reference.
    """

    merge_map = {}
    for _, row in merge_df.iterrows():
        source = row["source_transcript_id"]
        target = row["target_transcript_id"]
        if not source or not target:
            continue
        if source in merge_map and merge_map[source] != target:
            raise ValueError(
                f"Conflicting merge targets for {source}: "
                f"{merge_map[source]} versus {target}. "
                "Regenerate the merge log from a single transcript EM run."
            )
        merge_map[source] = target
    # Validate even chains that are not reached by the truth-positive IDs.
    for source in merge_map:
        follow_merge_target(source, merge_map)
    return merge_map


def follow_merge_target(reference_id, merge_map):
    """
    Follow a merge chain until the terminal target reference.

    Example:
        A -> B -> C
    returns:
        C, [A, B, C]
    """

    chain = [
        reference_id
    ]

    seen = {
        reference_id
    }

    current = reference_id

    while current in merge_map:

        nxt = merge_map[
            current
        ]

        chain.append(
            nxt
        )

        if nxt in seen:
            raise ValueError("Cycle in transcript merge log: " + " -> ".join(chain))

        seen.add(
            nxt
        )

        current = nxt

    return (
        chain[-1],
        chain,
    )


def metadata_truth_reference_map(
    metadata,
    truth_ids,
):
    """
    Map truth ENST IDs to pre-EM RiboTxAlign reference IDs.

    The metadata can represent the same Gencode transcript in several
    columns, so all relevant fields are checked.
    """

    mapping = {}

    for _, row in metadata.iterrows():

        candidates = set()

        for col in [
            "source_reference_id_norm",
            "representative_transcript_id_norm",
            "original_transcript_id_norm",
        ]:
            tid = row.get(
                col
            )

            if tid:
                candidates.add(
                    tid
                )

        for col in [
            "gencode_transcript_ids",
            "source_transcript_ids",
        ]:
            candidates.update(
                split_ids(
                    row.get(
                        col,
                        "",
                    )
                )
            )

        for tid in candidates:

            if tid not in truth_ids:
                continue

            mapping.setdefault(
                tid,
                [],
            ).append({
                "reference_id":
                    row["reference_id"],

                "overlap_group_id":
                    row.get(
                        "overlap_group_id",
                        "",
                    ),

                "resolution_type":
                    row.get(
                        "resolution_type",
                        "",
                    ),
            })

    return mapping


def make_assignment(
    truth_id,
    rows,
    match_type,
    prefix,
):
    return {
        "transcript_id":
            truth_id,

        f"{prefix}_tpm":
            rows["tpm"].sum(),

        f"{prefix}_expected_count":
            rows["expected_count"].sum(),

        f"{prefix}_match_type":
            match_type,

        f"{prefix}_reference_ids":
            ";".join(
                rows["reference_id"]
            ),

        f"{prefix}_overlap_group_id":
            ";".join(
                sorted(
                    set(
                        rows["overlap_group_id"]
                    )
                )
            ),

        f"{prefix}_merge_chain":
            "",
    }


def match_ribotx(
    truth,
    ribotx,
    transcript_gene,
    metadata,
    merge_df,
    prefix,
):
    """
    Matching order:

    1. Exact source match in the final TPM table
    2. EM-merge target rescue
       - If a truth transcript's pre-EM reference was merged into another
         reference, use the FINAL TPM of the terminal merge target as-is.
       - No redistribution is performed.
    3. Unique overlap-group rescue
    4. Gene-level TPM redistribution
    5. Unresolved
    """

    truth_ids = set(
        truth["transcript_id"]
    )

    used_refs = set()
    assignments = []

    # --------------------------------------------------------
    # 1. Exact source
    # --------------------------------------------------------

    exact = ribotx[
        ribotx[
            "evaluation_transcript_id_norm"
        ].isin(
            truth_ids
        )
    ]

    for truth_id, rows in (
        exact.groupby(
            "evaluation_transcript_id_norm"
        )
    ):
        assignments.append(
            make_assignment(
                truth_id,
                rows,
                "exact_source",
                prefix,
            )
        )

        used_refs.update(
            rows["reference_id"]
        )

    assigned = {
        x["transcript_id"]
        for x in assignments
    }

    unmatched = (
        truth_ids
        -
        assigned
    )

    # --------------------------------------------------------
    # 2. EM-merge target rescue
    #
    # A source reference may disappear from the final TPM table because
    # transcript-level EM merged it into another reference.
    #
    # For evaluation, keep the official RiboTxAlign output unchanged:
    # simply assign the terminal target's FINAL TPM to the truth
    # transcript represented by the merged source.
    #
    # IMPORTANT:
    #   This is NOT redistribution. If multiple source references were
    #   merged into the same target, each source truth transcript sees
    #   the same terminal target TPM for evaluation.
    # --------------------------------------------------------

    merge_map = build_merge_map(
        merge_df
    )

    truth_to_metadata_refs = (
        metadata_truth_reference_map(
            metadata,
            truth_ids,
        )
    )

    final_by_reference = {
        ref:
            rows.copy()

        for ref, rows in (
            ribotx.groupby(
                "reference_id"
            )
        )
    }

    merge_rescue_records = []

    for truth_id in sorted(
        unmatched
    ):

        metadata_refs = (
            truth_to_metadata_refs.get(
                truth_id,
                [],
            )
        )

        rescued_targets = []
        chains = []

        for item in metadata_refs:

            source_ref = item[
                "reference_id"
            ]

            if source_ref not in merge_map:
                continue

            target_ref, chain = (
                follow_merge_target(
                    source_ref,
                    merge_map,
                )
            )

            chains.append(
                " -> ".join(
                    chain
                )
            )

            if target_ref in final_by_reference:
                rescued_targets.append(
                    target_ref
                )

        rescued_targets = list(
            dict.fromkeys(
                rescued_targets
            )
        )

        if not rescued_targets:
            continue

        rows = pd.concat(
            [
                final_by_reference[
                    ref
                ]
                for ref in rescued_targets
            ],
            ignore_index=True,
        )

        assignments.append({
            "transcript_id":
                truth_id,

            f"{prefix}_tpm":
                rows["tpm"].sum(),

            f"{prefix}_expected_count":
                rows[
                    "expected_count"
                ].sum(),

            f"{prefix}_match_type":
                "em_merge_target",

            f"{prefix}_reference_ids":
                ";".join(
                    rescued_targets
                ),

            f"{prefix}_overlap_group_id":
                ";".join(
                    sorted(
                        set(
                            rows[
                                "overlap_group_id"
                            ]
                        )
                    )
                ),

            f"{prefix}_merge_chain":
                " | ".join(
                    chains
                ),
        })

        merge_rescue_records.append({
            "transcript_id":
                truth_id,

            "target_reference_ids":
                ";".join(
                    rescued_targets
                ),

            "merge_chains":
                " | ".join(
                    chains
                ),
        })

        # Do not mark the target as unavailable for another merged source.
        # Several source transcripts may legitimately point to the same
        # terminal target in this evaluation.
        used_refs.update(
            rescued_targets
        )

    assigned = {
        x["transcript_id"]
        for x in assignments
    }

    unmatched = (
        truth_ids
        -
        assigned
    )

    # --------------------------------------------------------
    # 3. Overlap-group rescue
    # --------------------------------------------------------

    truth_to_groups = {}

    for _, row in ribotx.iterrows():

        candidate_ids = {
            row[
                "source_reference_id_norm"
            ],
            row[
                "representative_transcript_id_norm"
            ],
        }

        candidate_ids.update(
            split_ids(
                row.get(
                    "gencode_transcript_ids",
                    "",
                )
            )
        )

        for tid in candidate_ids:
            if tid in unmatched:
                truth_to_groups.setdefault(
                    tid,
                    set(),
                ).add(
                    row["overlap_group_id"]
                )

    group_to_truth = {}

    for tid, groups in (
        truth_to_groups.items()
    ):
        if len(groups) == 1:
            group_id = next(
                iter(groups)
            )

            group_to_truth.setdefault(
                group_id,
                [],
            ).append(
                tid
            )

    unused = ribotx[
        ~ribotx[
            "reference_id"
        ].isin(
            used_refs
        )
    ]

    for group_id, truth_candidates in (
        group_to_truth.items()
    ):
        rows = unused[
            unused[
                "overlap_group_id"
            ]
            ==
            group_id
        ]

        if (
            len(truth_candidates) == 1
            and len(rows) == 1
        ):
            tid = truth_candidates[0]

            assignments.append(
                make_assignment(
                    tid,
                    rows,
                    "unique_group_rescue",
                    prefix,
                )
            )

            used_refs.update(
                rows["reference_id"]
            )

    # --------------------------------------------------------
    # 3. Gene-level redistribution
    # --------------------------------------------------------

    assigned = {
        x["transcript_id"]
        for x in assignments
    }

    unmatched = (
        truth_ids
        -
        assigned
    )

    truth_gene = (
        transcript_gene[
            transcript_gene[
                "transcript_id"
            ].isin(
                unmatched
            )
        ]
    )

    unused_positive = ribotx[
        (
            ~ribotx[
                "reference_id"
            ].isin(
                used_refs
            )
        )
        &
        (
            ribotx["tpm"] > 0
        )
    ]

    gene_to_truth = {
        gene:
            group[
                "transcript_id"
            ].tolist()
        for gene, group in (
            truth_gene.groupby(
                "truth_gene_id"
            )
        )
    }

    for gene_id, truth_candidates in (
        gene_to_truth.items()
    ):
        rows = unused_positive[
            unused_positive[
                "reference_gene_id_norm"
            ]
            ==
            gene_id
        ]

        if len(rows) == 0:
            continue

        total_tpm = rows["tpm"].sum()
        total_count = (
            rows["expected_count"].sum()
        )

        n_truth = len(
            truth_candidates
        )

        refs = ";".join(
            rows["reference_id"]
        )

        groups = ";".join(
            sorted(
                set(
                    rows["overlap_group_id"]
                )
            )
        )

        for tid in truth_candidates:
            assignments.append({
                "transcript_id":
                    tid,

                f"{prefix}_tpm":
                    total_tpm / n_truth,

                f"{prefix}_expected_count":
                    total_count / n_truth,

                f"{prefix}_match_type":
                    "gene_average_rescue",

                f"{prefix}_reference_ids":
                    refs,

                f"{prefix}_overlap_group_id":
                    groups,

                f"{prefix}_merge_chain":
                    "",
            })

        used_refs.update(
            rows["reference_id"]
        )

    # --------------------------------------------------------
    # 4. Unresolved
    # --------------------------------------------------------

    assigned = {
        x["transcript_id"]
        for x in assignments
    }

    unresolved_ids = (
        truth_ids
        -
        assigned
    )

    unresolved = []

    for tid in sorted(
        unresolved_ids
    ):
        gene_id_series = (
            transcript_gene.loc[
                transcript_gene[
                    "transcript_id"
                ]
                ==
                tid,
                "truth_gene_id",
            ]
        )

        gene_id = (
            gene_id_series.iloc[0]
            if len(gene_id_series)
            else None
        )

        if gene_id is None:
            category = (
                "gene_id_not_found"
            )

        elif not (
            ribotx[
                "reference_gene_id_norm"
            ]
            ==
            gene_id
        ).any():
            category = (
                "gene_not_present"
            )

        else:
            category = (
                "same_gene_no_unused_reference"
            )

        unresolved.append({
            "transcript_id":
                tid,

            "truth_gene_id":
                gene_id,

            "category":
                category,
        })

        assignments.append({
            "transcript_id":
                tid,

            f"{prefix}_tpm":
                0.0,

            f"{prefix}_expected_count":
                0.0,

            f"{prefix}_match_type":
                "unresolved",

            f"{prefix}_reference_ids":
                "",

            f"{prefix}_overlap_group_id":
                "",

            f"{prefix}_merge_chain":
                "",
        })

    result = pd.DataFrame(
        assignments
    )

    unresolved_df = pd.DataFrame(
        unresolved
    )

    match_col = f"{prefix}_match_type"

    counts = (
        result[match_col]
        .value_counts()
    )

    matched = (
        len(truth_ids)
        -
        len(unresolved_ids)
    )

    print()
    print("=" * 72)
    print(f"{prefix} MATCHING")
    print("=" * 72)

    print(
        "Truth transcripts:",
        len(truth_ids),
    )

    for name in [
        "exact_source",
        "em_merge_target",
        "unique_group_rescue",
        "gene_average_rescue",
        "unresolved",
    ]:
        print(
            f"{name:24s}:",
            int(
                counts.get(
                    name,
                    0,
                )
            ),
        )

    print(
        "Total matched:",
        matched,
    )

    print(
        "Recovery:",
        f"{matched / len(truth_ids) * 100:.2f}%",
    )

    print(
        "RiboTxAlign references:",
        len(ribotx),
    )

    print(
        "Used references:",
        len(used_refs),
    )

    print(
        "Unused references:",
        len(ribotx) - len(used_refs),
    )

    if len(unresolved_df):
        print()
        print("Unresolved categories:")
        print(
            unresolved_df[
                "category"
            ]
            .value_counts()
            .to_string()
        )

    return (
        result,
        unresolved_df,
    )


# ============================================================
# RiboMap
# ============================================================

def read_cds(path):
    records = []

    with open(path) as f:
        for line in f:
            if (
                not line.strip()
                or line.startswith("#")
            ):
                continue

            fields = line.split()

            if len(fields) < 3:
                continue

            try:
                length = (
                    int(fields[2])
                    -
                    int(fields[1])
                )
            except ValueError:
                continue

            if length > 0:
                records.append({
                    "transcript_id":
                        norm_id(
                            fields[0]
                        ),

                    "cds_length":
                        length,
                })

    return (
        pd.DataFrame(records)
        .groupby(
            "transcript_id",
            as_index=False,
        )["cds_length"]
        .max()
    )


def read_ribomap(
    stats_path,
    cds,
):
    records = []
    tid = None

    with open(stats_path) as f:
        for line in f:
            line = line.strip()

            if line.startswith("tid:"):
                tid = norm_id(
                    line.split(
                        ":",
                        1,
                    )[1]
                )

            elif (
                line.startswith("rabd:")
                and tid
            ):
                records.append({
                    "transcript_id":
                        tid,

                    "count":
                        float(
                            line.split(
                                ":",
                                1,
                            )[1]
                        ),
                })

    df = (
        pd.DataFrame(records)
        .groupby(
            "transcript_id",
            as_index=False,
        )["count"]
        .sum()
        .merge(
            cds,
            on="transcript_id",
            how="left",
        )
    )

    print(
        "RiboMap transcripts without CDS length:",
        df["cds_length"]
        .isna()
        .sum(),
    )

    return calc_tpm(
        df,
        "count",
        "cds_length",
        "RiboMap_tpm",
    )


# ============================================================
# DeepShape
# ============================================================

def read_deepshape(path):
    records = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if (
                not line
                or line.startswith("Loop ")
            ):
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                continue

            try:
                tpm = float(fields[1])
            except ValueError:
                continue

            records.append({
                "transcript_id":
                    norm_id(
                        fields[0]
                    ),

                "DeepShape_tpm":
                    tpm,
            })

    return (
        pd.DataFrame(records)
        .groupby(
            "transcript_id",
            as_index=False,
        )["DeepShape_tpm"]
        .sum()
    )


# ============================================================
# Evaluation
# ============================================================

def evaluate(
    truth,
    rna,
    ribotx_results,
    ribomap,
    deepshape,
):
    df = truth.merge(
        rna,
        on="transcript_id",
        how="left",
    )

    for method_name, method_df in (
        ribotx_results.items()
    ):
        df = df.merge(
            method_df,
            on="transcript_id",
            how="left",
        )

    df = (
        df
        .merge(
            ribomap[
                [
                    "transcript_id",
                    "RiboMap_tpm",
                ]
            ],
            on="transcript_id",
            how="left",
        )
        .merge(
            deepshape,
            on="transcript_id",
            how="left",
        )
    )

    methods = {
        method_name:
            f"{method_name}_tpm"
        for method_name in ribotx_results
    }

    methods.update({
        "RiboMap":
            "RiboMap_tpm",

        "DeepShape":
            "DeepShape_tpm",
    })

    for col in methods.values():
        if col not in df.columns:
            df[col] = 0.0

        df[col] = (
            df[col]
            .fillna(0.0)
        )

    df = df[
        df["truth_count"] > 0
    ].copy()

    # --------------------------------------------------------
    # Recovery
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("RECOVERY")
    print("=" * 72)

    recovery_records = []

    for name, col in methods.items():
        n = int(
            (df[col] > 0).sum()
        )

        recovery_records.append({
            "method": name,
            "recovered_transcripts": n,
            "truth_transcripts": len(df),
            "recovery_percent":
                n / len(df) * 100,
        })

        print(
            f"{name:24s}: "
            f"{n} / {len(df)} "
            f"({n / len(df) * 100:.2f}%)"
        )

    # --------------------------------------------------------
    # Method-specific positive accuracy
    # --------------------------------------------------------

    accuracy_records = []

    print()
    print("=" * 72)
    print("METHOD-SPECIFIC ACCURACY")
    print("=" * 72)

    results = []

    for name, col in methods.items():
        subset = df[
            df[col] > 0
        ]

        result = accuracy(
            subset,
            col,
            name,
            f"{name}_positive",
        )

        results.append(result)
        accuracy_records.append(result)

    print(
        pd.DataFrame(results)
        .to_string(index=False)
    )

    # --------------------------------------------------------
    # Common-positive accuracy
    # --------------------------------------------------------

    positive_mask = pd.Series(
        True,
        index=df.index,
    )

    for col in methods.values():
        positive_mask &= (
            df[col] > 0
        )

    common = df[
        positive_mask
    ].copy()

    print()
    print("=" * 72)
    print("COMMON POSITIVE ACCURACY")
    print("=" * 72)

    print(
        "Common transcripts:",
        len(common),
    )

    results = []

    for name, col in methods.items():
        result = accuracy(
            common,
            col,
            name,
            "Common_positive",
        )

        results.append(result)
        accuracy_records.append(result)

    print(
        pd.DataFrame(results)
        .to_string(index=False)
    )

    # --------------------------------------------------------
    # End-to-end accuracy
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("END-TO-END ACCURACY")
    print("=" * 72)

    results = []

    for name, col in methods.items():
        result = accuracy(
            df,
            col,
            name,
            "All_truth_positive",
        )

        results.append(result)
        accuracy_records.append(result)

    print(
        pd.DataFrame(results)
        .to_string(index=False)
    )

    # --------------------------------------------------------
    # TPM mass on truth IDs
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("TPM MASS ON TRUTH IDS")
    print("=" * 72)

    print(
        "Truth                   :",
        df["Truth_tpm"].sum(),
    )

    tpm_mass_records = [{
        "method": "Truth",
        "tpm_mass_on_truth_ids":
            df["Truth_tpm"].sum(),
    }]

    for name, col in methods.items():
        value = df[col].sum()

        print(
            f"{name:24s}:",
            value,
        )

        tpm_mass_records.append({
            "method": name,
            "tpm_mass_on_truth_ids": value,
        })

    return (
        df,
        pd.DataFrame(
            accuracy_records
        ),
        pd.DataFrame(
            recovery_records
        ),
        pd.DataFrame(
            tpm_mass_records
        ),
    )


# ============================================================
# Gene-level evaluation
# ============================================================

def aggregate_truth_to_gene(
    truth,
    transcript_gene,
):
    """
    Aggregate transcript-level truth TPM to gene-level TPM.
    """

    df = (
        truth
        .merge(
            transcript_gene,
            on="transcript_id",
            how="left",
        )
    )

    missing = int(
        df["truth_gene_id"]
        .isna()
        .sum()
    )

    if missing:
        print(
            "Truth transcripts without gene ID:",
            missing,
        )

    gene = (
        df.dropna(
            subset=["truth_gene_id"]
        )
        .groupby(
            "truth_gene_id",
            as_index=False,
        )
        .agg(
            Truth_gene_tpm=(
                "Truth_tpm",
                "sum",
            ),
            Truth_gene_count=(
                "truth_count",
                "sum",
            ),
        )
        .rename(
            columns={
                "truth_gene_id":
                    "gene_id",
            }
        )
    )

    return gene


def make_transcript_to_gene_map(
    transcript_gene,
):
    return dict(
        zip(
            transcript_gene[
                "transcript_id"
            ],
            transcript_gene[
                "truth_gene_id"
            ],
        )
    )


def get_ribotx_gene_id(
    row,
    transcript_to_gene,
):
    """
    Determine a gene ID for one RiboTxAlign reference.

    Priority:
      1. reference_gene_id
      2. source_reference_id -> annotation
      3. representative_transcript_id -> annotation
      4. gencode_transcript_ids -> annotation
    """

    gene_id = row.get(
        "reference_gene_id_norm"
    )

    if (
        gene_id
        and isinstance(
            gene_id,
            str,
        )
        and gene_id.startswith("ENSG")
    ):
        return gene_id

    for col in [
        "source_reference_id_norm",
        "representative_transcript_id_norm",
    ]:
        tid = row.get(col)

        if tid in transcript_to_gene:
            return transcript_to_gene[tid]

    for tid in split_ids(
        row.get(
            "gencode_transcript_ids",
            "",
        )
    ):
        if tid in transcript_to_gene:
            return transcript_to_gene[tid]

    return None


def aggregate_ribotx_raw_to_gene(
    ribotx,
    transcript_gene,
    method_name,
):
    """
    Aggregate the RAW RiboTxAlign output to gene level.

    This is intentionally based on the method's own reported output,
    rather than only on truth-matched transcripts. Therefore genes
    reported outside the truth-positive set are retained and can be
    counted in the OR analysis.
    """

    df = ribotx.copy()

    transcript_to_gene = (
        make_transcript_to_gene_map(
            transcript_gene
        )
    )

    df["gene_id"] = df.apply(
        lambda row: get_ribotx_gene_id(
            row,
            transcript_to_gene,
        ),
        axis=1,
    )

    positive = df[
        df["tpm"] > 0
    ].copy()

    unmapped_positive_refs = int(
        positive["gene_id"]
        .isna()
        .sum()
    )

    gene = (
        positive
        .dropna(
            subset=["gene_id"]
        )
        .groupby(
            "gene_id",
            as_index=False,
        )["tpm"]
        .sum()
        .rename(
            columns={
                "tpm":
                    f"{method_name}_gene_tpm",
            }
        )
    )

    return (
        gene,
        unmapped_positive_refs,
    )


def aggregate_transcript_method_to_gene(
    method_df,
    transcript_gene,
    transcript_tpm_col,
    gene_tpm_col,
):
    """
    Aggregate a transcript-level method output to genes.
    Only method-positive transcripts are considered reported.
    """

    df = (
        method_df[
            [
                "transcript_id",
                transcript_tpm_col,
            ]
        ]
        .copy()
    )

    df[transcript_tpm_col] = (
        pd.to_numeric(
            df[transcript_tpm_col],
            errors="coerce",
        )
        .fillna(0.0)
    )

    df = df[
        df[transcript_tpm_col] > 0
    ]

    df = df.merge(
        transcript_gene,
        on="transcript_id",
        how="left",
    )

    unmapped_positive_transcripts = int(
        df["truth_gene_id"]
        .isna()
        .sum()
    )

    gene = (
        df.dropna(
            subset=["truth_gene_id"]
        )
        .groupby(
            "truth_gene_id",
            as_index=False,
        )[transcript_tpm_col]
        .sum()
        .rename(
            columns={
                "truth_gene_id":
                    "gene_id",
                transcript_tpm_col:
                    gene_tpm_col,
            }
        )
    )

    return (
        gene,
        unmapped_positive_transcripts,
    )


def correlation_record(
    df,
    method_col,
    method_name,
    set_name,
):
    """
    Calculate correlation and error metrics for a gene set.
    """

    n = len(df)

    if n < 2:
        return {
            "method":
                method_name,
            "set":
                set_name,
            "genes":
                n,
            "pearson":
                np.nan,
            "spearman":
                np.nan,
            "mae":
                np.nan,
            "rmse":
                np.nan,
        }

    y = df[
        "Truth_gene_tpm"
    ].to_numpy(float)

    p = df[
        method_col
    ].to_numpy(float)

    return {
        "method":
            method_name,
        "set":
            set_name,
        "genes":
            n,
        "pearson":
            pearsonr(
                y,
                p,
            ).statistic,
        "spearman":
            spearmanr(
                y,
                p,
            ).statistic,
        "mae":
            mean_absolute_error(
                y,
                p,
            ),
        "rmse":
            np.sqrt(
                mean_squared_error(
                    y,
                    p,
                )
            ),
    }


def evaluate_gene_level(
    truth,
    transcript_gene,
    ribotx_raw_results,
    ribomap,
    deepshape,
):
    """
    Gene-level comparison.

    For each method:
      reported genes:
          method gene TPM > 0

      truth-positive genes:
          truth gene TPM > 0

      AND:
          truth-positive AND method-positive

      OR:
          truth-positive OR method-positive
          Missing expression on either side is filled with 0.

    In addition, an all-method AND and all-method OR set are generated.
    """

    truth_gene = aggregate_truth_to_gene(
        truth,
        transcript_gene,
    )

    truth_gene = truth_gene[
        truth_gene[
            "Truth_gene_tpm"
        ] > 0
    ].copy()

    truth_positive_ids = set(
        truth_gene["gene_id"]
    )

    method_gene_tables = {}
    mapping_records = []

    # --------------------------------------------------------
    # RiboTxAlign methods
    # --------------------------------------------------------

    for method_name, raw in (
        ribotx_raw_results.items()
    ):
        (
            gene_df,
            unmapped,
        ) = aggregate_ribotx_raw_to_gene(
            raw,
            transcript_gene,
            method_name,
        )

        method_gene_tables[
            method_name
        ] = gene_df

        mapping_records.append({
            "method":
                method_name,
            "positive_records_without_gene_mapping":
                unmapped,
        })

    # --------------------------------------------------------
    # RiboMap
    # --------------------------------------------------------

    (
        ribomap_gene,
        ribomap_unmapped,
    ) = aggregate_transcript_method_to_gene(
        ribomap,
        transcript_gene,
        "RiboMap_tpm",
        "RiboMap_gene_tpm",
    )

    method_gene_tables[
        "RiboMap"
    ] = ribomap_gene

    mapping_records.append({
        "method":
            "RiboMap",
        "positive_records_without_gene_mapping":
            ribomap_unmapped,
    })

    # --------------------------------------------------------
    # DeepShape
    # --------------------------------------------------------

    (
        deepshape_gene,
        deepshape_unmapped,
    ) = aggregate_transcript_method_to_gene(
        deepshape,
        transcript_gene,
        "DeepShape_tpm",
        "DeepShape_gene_tpm",
    )

    method_gene_tables[
        "DeepShape"
    ] = deepshape_gene

    mapping_records.append({
        "method":
            "DeepShape",
        "positive_records_without_gene_mapping":
            deepshape_unmapped,
    })

    # --------------------------------------------------------
    # Counts and pairwise AND / OR correlations
    # --------------------------------------------------------

    count_records = []
    correlation_records = []
    pairwise_tables = {}

    print()
    print("=" * 88)
    print("GENE-LEVEL REPORTED / TRUTH OVERLAP")
    print("=" * 88)

    print(
        "Truth TPM > 0 genes:",
        len(truth_positive_ids),
    )

    for method_name, method_gene in (
        method_gene_tables.items()
    ):
        method_col = (
            f"{method_name}_gene_tpm"
        )

        method_positive_ids = set(
            method_gene[
                method_gene[
                    method_col
                ] > 0
            ]["gene_id"]
        )

        overlap_ids = (
            truth_positive_ids
            &
            method_positive_ids
        )

        union_ids = (
            truth_positive_ids
            |
            method_positive_ids
        )

        truth_only_ids = (
            truth_positive_ids
            -
            method_positive_ids
        )

        method_only_ids = (
            method_positive_ids
            -
            truth_positive_ids
        )

        count_records.append({
            "method":
                method_name,
            "truth_positive_genes":
                len(truth_positive_ids),
            "reported_genes":
                len(method_positive_ids),
            "overlap_genes":
                len(overlap_ids),
            "truth_only_genes":
                len(truth_only_ids),
            "method_only_genes":
                len(method_only_ids),
            "union_genes":
                len(union_ids),
            "overlap_vs_truth_percent":
                (
                    len(overlap_ids)
                    / len(truth_positive_ids)
                    * 100
                    if truth_positive_ids
                    else np.nan
                ),
            "overlap_vs_reported_percent":
                (
                    len(overlap_ids)
                    / len(method_positive_ids)
                    * 100
                    if method_positive_ids
                    else np.nan
                ),
        })

        print()
        print(method_name)
        print(
            f"  reported genes : "
            f"{len(method_positive_ids)}"
        )
        print(
            f"  overlap genes  : "
            f"{len(overlap_ids)} / "
            f"{len(truth_positive_ids)} truth-positive"
        )
        print(
            f"  method only    : "
            f"{len(method_only_ids)}"
        )
        print(
            f"  truth only     : "
            f"{len(truth_only_ids)}"
        )
        print(
            f"  OR / union     : "
            f"{len(union_ids)}"
        )

        # ----------------------------------------------------
        # AND = overlap only
        # ----------------------------------------------------

        and_df = (
            truth_gene[
                truth_gene[
                    "gene_id"
                ].isin(
                    overlap_ids
                )
            ]
            .merge(
                method_gene,
                on="gene_id",
                how="inner",
            )
        )

        and_record = correlation_record(
            and_df,
            method_col,
            method_name,
            "AND_overlap",
        )

        correlation_records.append(
            and_record
        )

        # ----------------------------------------------------
        # OR = union, missing side = 0
        # ----------------------------------------------------

        union_index = pd.DataFrame({
            "gene_id":
                sorted(
                    union_ids
                )
        })

        or_df = (
            union_index
            .merge(
                truth_gene[
                    [
                        "gene_id",
                        "Truth_gene_tpm",
                    ]
                ],
                on="gene_id",
                how="left",
            )
            .merge(
                method_gene,
                on="gene_id",
                how="left",
            )
        )

        or_df[
            "Truth_gene_tpm"
        ] = (
            or_df[
                "Truth_gene_tpm"
            ]
            .fillna(0.0)
        )

        or_df[
            method_col
        ] = (
            or_df[
                method_col
            ]
            .fillna(0.0)
        )

        or_record = correlation_record(
            or_df,
            method_col,
            method_name,
            "OR_union_missing_zero",
        )

        correlation_records.append(
            or_record
        )

        pairwise = (
            or_df
            .copy()
        )

        pairwise[
            "in_truth"
        ] = (
            pairwise[
                "gene_id"
            ].isin(
                truth_positive_ids
            )
        )

        pairwise[
            "in_method"
        ] = (
            pairwise[
                "gene_id"
            ].isin(
                method_positive_ids
            )
        )

        pairwise[
            "set_membership"
        ] = np.select(
            [
                (
                    pairwise["in_truth"]
                    &
                    pairwise["in_method"]
                ),
                (
                    pairwise["in_truth"]
                    &
                    ~pairwise["in_method"]
                ),
                (
                    ~pairwise["in_truth"]
                    &
                    pairwise["in_method"]
                ),
            ],
            [
                "overlap",
                "truth_only",
                "method_only",
            ],
            default="neither",
        )

        pairwise_tables[
            method_name
        ] = pairwise

        print(
            "  AND correlation: "
            f"Pearson={and_record['pearson']:.4f}, "
            f"Spearman={and_record['spearman']:.4f}"
        )

        print(
            "  OR correlation : "
            f"Pearson={or_record['pearson']:.4f}, "
            f"Spearman={or_record['spearman']:.4f}"
        )

    # --------------------------------------------------------
    # Global AND / OR across all methods + truth
    # --------------------------------------------------------

    method_positive_sets = {}

    for method_name, method_gene in (
        method_gene_tables.items()
    ):
        method_col = (
            f"{method_name}_gene_tpm"
        )

        method_positive_sets[
            method_name
        ] = set(
            method_gene[
                method_gene[
                    method_col
                ] > 0
            ]["gene_id"]
        )

    all_and_ids = set(
        truth_positive_ids
    )

    all_or_ids = set(
        truth_positive_ids
    )

    for ids in (
        method_positive_sets.values()
    ):
        all_and_ids &= ids
        all_or_ids |= ids

    print()
    print("=" * 88)
    print("ALL-METHOD GENE SETS")
    print("=" * 88)
    print(
        "ALL methods AND truth:",
        len(all_and_ids),
    )
    print(
        "ALL methods OR truth :",
        len(all_or_ids),
    )

    global_set_records = [{
        "set":
            "ALL_METHODS_AND_TRUTH",
        "genes":
            len(all_and_ids),
    }, {
        "set":
            "ALL_METHODS_OR_TRUTH",
        "genes":
            len(all_or_ids),
    }]

    global_correlation_records = []

    for method_name, method_gene in (
        method_gene_tables.items()
    ):
        method_col = (
            f"{method_name}_gene_tpm"
        )

        # All-method AND
        global_and = (
            pd.DataFrame({
                "gene_id":
                    sorted(
                        all_and_ids
                    )
            })
            .merge(
                truth_gene[
                    [
                        "gene_id",
                        "Truth_gene_tpm",
                    ]
                ],
                on="gene_id",
                how="left",
            )
            .merge(
                method_gene,
                on="gene_id",
                how="left",
            )
        )

        global_correlation_records.append(
            correlation_record(
                global_and,
                method_col,
                method_name,
                "ALL_METHODS_AND_TRUTH",
            )
        )

        # All-method OR, missing = zero
        global_or = (
            pd.DataFrame({
                "gene_id":
                    sorted(
                        all_or_ids
                    )
            })
            .merge(
                truth_gene[
                    [
                        "gene_id",
                        "Truth_gene_tpm",
                    ]
                ],
                on="gene_id",
                how="left",
            )
            .merge(
                method_gene,
                on="gene_id",
                how="left",
            )
        )

        global_or[
            "Truth_gene_tpm"
        ] = (
            global_or[
                "Truth_gene_tpm"
            ]
            .fillna(0.0)
        )

        global_or[
            method_col
        ] = (
            global_or[
                method_col
            ]
            .fillna(0.0)
        )

        global_correlation_records.append(
            correlation_record(
                global_or,
                method_col,
                method_name,
                "ALL_METHODS_OR_TRUTH_missing_zero",
            )
        )

    return {
        "gene_truth":
            truth_gene,

        "gene_method_counts":
            pd.DataFrame(
                count_records
            ),

        "gene_pairwise_correlations":
            pd.DataFrame(
                correlation_records
            ),

        "gene_global_set_counts":
            pd.DataFrame(
                global_set_records
            ),

        "gene_global_correlations":
            pd.DataFrame(
                global_correlation_records
            ),

        "gene_mapping_diagnostics":
            pd.DataFrame(
                mapping_records
            ),

        "pairwise_tables":
            pairwise_tables,
    }


# ============================================================
# Main
# ============================================================

def main():

    truth = read_truth(
        BASE / "HEK293T.simulation.profile"
    )

    rna = read_rna(
        BASE / "HEK293T.simulation_abundance.sf"
    )

    transcript_gene = (
        read_transcript_gene(
            ANNOTATION
        )
    )

    # --------------------------------------------------------
    # RiboTxAlign methods
    # --------------------------------------------------------

    ribotx_results = {}
    ribotx_raw_results = {}
    unresolved_results = {}

    print()
    print("=" * 72)
    print("RiboTxAlign INPUT FILES")
    print("=" * 72)

    for method_name, method_dir in (
        RIBOTX_METHOD_DIRS.items()
    ):
        result_path = find_ribotx_result(
            method_dir
        )

        print(
            f"{method_name:24s}: "
            f"{result_path}"
        )

        metadata_path = (
            find_ribotx_metadata(
                method_dir
            )
        )

        merge_path = (
            find_ribotx_merge(
                method_dir
            )
        )

        print(
            f"{'metadata':24s}: "
            f"{metadata_path}"
        )

        print(
            f"{'merge':24s}: "
            f"{merge_path}"
        )

        raw = read_ribotx(
            result_path
        )

        metadata = (
            read_ribotx_metadata(
                metadata_path
            )
        )

        merge_df = (
            read_ribotx_merge(
                merge_path
            )
        )

        ribotx_raw_results[
            method_name
        ] = raw

        matched, unresolved = (
            match_ribotx(
                truth,
                raw,
                transcript_gene,
                metadata,
                merge_df,
                method_name,
            )
        )

        ribotx_results[
            method_name
        ] = matched

        unresolved_results[
            method_name
        ] = unresolved

    # --------------------------------------------------------
    # RiboMap / DeepShape
    # --------------------------------------------------------

    ribomap = read_ribomap(
        RIBOMAP_STATS,
        read_cds(
            BASE
            / "HEK293T.simulation_cds.txt"
        ),
    )

    deepshape = read_deepshape(
        DEEPSHAPE_LOG
    )

    # --------------------------------------------------------
    # Gene-level evaluation
    # --------------------------------------------------------

    gene_results = evaluate_gene_level(
        truth,
        transcript_gene,
        ribotx_raw_results,
        ribomap,
        deepshape,
    )

    # --------------------------------------------------------
    # Transcript-level evaluation
    # --------------------------------------------------------

    (
        comparison,
        accuracy_df,
        recovery_df,
        tpm_mass_df,
    ) = evaluate(
        truth,
        rna,
        ribotx_results,
        ribomap,
        deepshape,
    )

    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    outputs = {
        "gene_truth_tpm.LongRead_RNASeq.tsv":
            gene_results[
                "gene_truth"
            ],

        "gene_report_overlap_counts.LongRead_RNASeq.tsv":
            gene_results[
                "gene_method_counts"
            ],

        "gene_pairwise_correlations.LongRead_RNASeq.tsv":
            gene_results[
                "gene_pairwise_correlations"
            ],

        "gene_global_set_counts.LongRead_RNASeq.tsv":
            gene_results[
                "gene_global_set_counts"
            ],

        "gene_global_correlations.LongRead_RNASeq.tsv":
            gene_results[
                "gene_global_correlations"
            ],

        "gene_mapping_diagnostics.LongRead_RNASeq.tsv":
            gene_results[
                "gene_mapping_diagnostics"
            ],

        "transcript_tpm_comparison.LongRead_RNASeq.tsv":
            comparison,

        "transcript_tpm_accuracy.LongRead_RNASeq.tsv":
            accuracy_df,

        "transcript_recovery.LongRead_RNASeq.tsv":
            recovery_df,

        "transcript_tpm_mass.LongRead_RNASeq.tsv":
            tpm_mass_df,
    }

    for method_name, matched in (
        ribotx_results.items()
    ):
        outputs[
            f"{method_name}.source_rescue_matching.tsv"
        ] = matched

    for method_name, unresolved in (
        unresolved_results.items()
    ):
        unresolved = (
            unresolved
            .merge(
                truth[
                    [
                        "transcript_id",
                        "truth_count",
                        "Truth_tpm",
                    ]
                ],
                on="transcript_id",
                how="left",
            )
            .merge(
                rna,
                on="transcript_id",
                how="left",
            )
        )

        outputs[
            f"{method_name}.unresolved_diagnosis.tsv"
        ] = unresolved

    for method_name, pairwise_df in (
        gene_results[
            "pairwise_tables"
        ].items()
    ):
        outputs[
            f"gene_pairwise_values.{method_name}.tsv"
        ] = pairwise_df

    print()
    print("=" * 72)
    print("OUTPUT FILES")
    print("=" * 72)

    for filename, df in outputs.items():
        path = BASE / filename

        df.to_csv(
            path,
            sep="\t",
            index=False,
        )

        print(path)


if __name__ == "__main__":
    main()
