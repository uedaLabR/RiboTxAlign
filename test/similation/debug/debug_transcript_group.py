"""Read-only, one-group diagnostic. Run with the pipeline's Linux Python.

No truth information is used to change mapping or EM. CDS overlap below is
alignment-span overlap, not P-site evidence of translation.
"""
from __future__ import annotations

import csv
import itertools
import json
from collections import defaultdict
from pathlib import Path


# Edit these paths for the Linux run, then execute this file without arguments.
SAMPLE_DIR = Path(
    "/mnt/share/project/riboML/validation_expression/similation/HEK293T/"
    "RiboTxAlign_LongRead/HEK293T_simulation_LongRead"
)
SAMPLE_NAME = "HEK293T_simulation_LongRead"
GROUP_ID = "chr16_plus_G00536"  # PPP4C
METADATA_TSV = SAMPLE_DIR / "reference" / f"{SAMPLE_NAME}.metadata.tsv"
READS_INFO_PARQUET = SAMPLE_DIR / "em" / f"{SAMPLE_NAME}.reads_info.parquet"
TRANSCRIPT_TPM_TSV = SAMPLE_DIR / "em" / f"{SAMPLE_NAME}.transcript_tpm.tsv"
SORTED_MAPPING_DIR = SAMPLE_DIR / "mapping" / "sorted"  # None if deleted
OUTPUT_DIR = SAMPLE_DIR / "diagnosis_ppp4c_01"  # Use a new directory each run
SORTED_MAPPING_DIR = None

def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def unique_file(root, pattern, explicit=None):
    if explicit:
        path = Path(explicit)
        if not path.is_file():
            raise FileNotFoundError(path)
        return path
    found = sorted(root.rglob(pattern))
    if len(found) != 1:
        raise ValueError(f"Expected one {pattern} under {root}, found {found}. "
                         "Set its explicit path in the file's configuration section.")
    return found[0]


def group_rows(paths, group):
    import pyarrow.dataset as ds
    for path in paths:
        dataset = ds.dataset(str(path), format="parquet")
        if "group_id" not in dataset.schema.names:
            raise ValueError(f"Missing group_id: {path}")
        scanner = dataset.scanner(filter=ds.field("group_id") == group,
                                  batch_size=65536)
        for batch in scanner.to_batches():
            yield from batch.to_pylist()


def span_class(start, end, metadata):
    """Coordinates are zero-based, half-open on the mapped reference."""
    try:
        left = int(metadata["cds_reference_start"])
        right = int(metadata["cds_reference_end"])
    except (KeyError, ValueError, TypeError):
        return "cds_unknown"
    if right <= left or start is None or end is None or end <= start:
        return "cds_unknown"
    if end <= left or start >= right:
        return "outside_cds"
    if start >= left and end <= right:
        return "inside_cds"
    return "crosses_cds_boundary"


def summarize(rows, metadata, stage):
    candidates = defaultdict(set)
    support = defaultdict(set)
    classes = defaultdict(lambda: defaultdict(set))
    positions = defaultdict(list)
    weights = defaultdict(float)
    for row in rows:
        tid, rid = row["transcript_id"], row["read_id"]
        candidates[rid].add(tid)
        support[tid].add(rid)
        alignment_tid = row.get("alignment_transcript_id") or tid
        start = row.get("reference_start", row.get("start"))
        end = row.get("reference_end", row.get("end"))
        category = span_class(start, end, metadata.get(alignment_tid, {}))
        classes[tid][category].add(rid)
        # Do not pool coordinates from references absorbed during merging.
        if alignment_tid == tid and start is not None and end is not None:
            positions[tid].append((start, end))
        if stage == "post_em":
            weights[tid] += float(row["prior_probability"]) * float(row["transcript_probability"])
    summaries = []
    for tid in sorted(support):
        spans = positions[tid]
        summaries.append(dict(
            transcript_id=tid, candidate_reads=len(support[tid]),
            unique_within_group_reads=sum(len(candidates[r]) == 1 for r in support[tid]),
            mapped_start_same_reference=min((s for s, e in spans), default=""),
            mapped_end_same_reference=max((e for s, e in spans), default=""),
            inside_cds_reads=len(classes[tid]["inside_cds"]),
            outside_cds_reads=len(classes[tid]["outside_cds"]),
            crosses_cds_boundary_reads=len(classes[tid]["crosses_cds_boundary"]),
            cds_unknown_reads=len(classes[tid]["cds_unknown"]),
            posterior_count=weights[tid] if stage == "post_em" else ""))
    pairs = []
    for a, b in itertools.combinations(sorted(support), 2):
        pairs.append(dict(transcript_a=a, transcript_b=b,
                          shared_reads=len(support[a] & support[b]),
                          a_only_reads=len(support[a] - support[b]),
                          b_only_reads=len(support[b] - support[a]),
                          identical_read_sets=support[a] == support[b]))
    return summaries, pairs


def main():
    root = SAMPLE_DIR
    if not root.is_dir():
        raise FileNotFoundError(f"Missing sample directory: {root}")
    paths = {
        "metadata": unique_file(root, "*metadata.tsv", METADATA_TSV),
        "reads_info": unique_file(root, "*.reads_info.parquet", READS_INFO_PARQUET),
        "transcript_tpm": unique_file(root, "*.transcript_tpm.tsv", TRANSCRIPT_TPM_TSV),
    }
    metadata = {r["reference_id"]: r for r in read_tsv(paths["metadata"])}
    post = list(group_rows([paths["reads_info"]], GROUP_ID))
    if not post:
        raise ValueError(f"No post-EM records for {GROUP_ID}")
    stages = {"post_em": post}
    warnings = ["Use all inputs from the same run; timestamps cannot prove provenance.",
                "CDS classification uses alignment spans, not CIGAR-resolved bases or P-sites.",
                "Unique means unique within this group, not globally unique.",
                "Read-set identity after EM cannot establish pre-EM merge eligibility."]
    if SORTED_MAPPING_DIR:
        raw_paths = sorted(SORTED_MAPPING_DIR.glob("*.parquet"))
        if not raw_paths:
            raise ValueError(f"No parquet files: {SORTED_MAPPING_DIR}")
    raw_paths = (
        sorted(SORTED_MAPPING_DIR.glob("*.parquet"))
        if SORTED_MAPPING_DIR else []
    )
    if raw_paths:
        pre = list(group_rows(raw_paths, GROUP_ID))
        if not pre:
            raise ValueError(f"No pre-EM records for {GROUP_ID}")
        stages["pre_em"] = pre
        paths.update({f"mapping_{i}": p for i, p in enumerate(raw_paths)})
    else:
        warnings.append("Pre-EM mapping not supplied: original unique-read support is unknown.")
        message = (
            f"Pre-EM mapping unavailable ({SORTED_MAPPING_DIR}): "
            "continuing with post-EM data only. Original unique-read support "
            "is unknown; pre_em outputs will not be generated."
        )
        warnings.append(message)
        print(f"NOTE: {message}")
    selected = {r["transcript_id"] for rows in stages.values() for r in rows}
    selected.update(r.get("alignment_transcript_id") for r in post)
    native = [r for r in read_tsv(paths["transcript_tpm"]) if r["transcript_id"] in selected]
    meta = [metadata[t] for t in sorted(selected - {None}) if t in metadata]
    if any(span_class(0, 1, r) == "cds_unknown" for r in meta):
        warnings.append("Some reference CDS coordinates unavailable; inspect metadata fields.")
    # Fail rather than silently mix an earlier diagnostic with a new run.
    OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    for stage, rows in stages.items():
        summaries, pairs = summarize(rows, metadata, stage)
        write_tsv(OUTPUT_DIR / f"{stage}_summary.tsv", summaries, list(summaries[0]))
        write_tsv(OUTPUT_DIR / f"{stage}_pairs.tsv", pairs,
                  ["transcript_a", "transcript_b", "shared_reads", "a_only_reads",
                   "b_only_reads", "identical_read_sets"])
        write_tsv(OUTPUT_DIR / f"{stage}_reads.tsv", rows, list(rows[0]))
    write_tsv(OUTPUT_DIR / "reference_metadata.tsv", meta,
              list(meta[0]) if meta else ["reference_id"])
    write_tsv(OUTPUT_DIR / "native_tpm.tsv", native,
              list(native[0]) if native else ["transcript_id"])
    manifest = dict(group_id=GROUP_ID, warnings=warnings, inputs={
        k: dict(path=str(p.resolve()), bytes=p.stat().st_size, mtime=p.stat().st_mtime)
        for k, p in paths.items()})
    (OUTPUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Output: {OUTPUT_DIR}")
    for warning in warnings:
        print(f"NOTE: {warning}")


if __name__ == "__main__":
    main()