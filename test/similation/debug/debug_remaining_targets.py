"""Inspect two residual cases using same-run outputs. No CLI arguments."""
import csv
import json
import re
from datetime import datetime
from pathlib import Path
from collections import Counter, defaultdict


SAMPLE_DIR = Path(
    "/mnt/share/project/riboML/validation_expression/similation/HEK293T/"
    "RiboTxAlign_LongRead/HEK293T_simulation_LongRead"
)
SAMPLE_NAME = "HEK293T_simulation_LongRead"
TARGETS = ("ENST00000647841", "ENST00000643777")
# Optional: annotation actually used to build the simulation, if available.
TRUTH_GTF = Path(
    "/mnt/share/project/riboML/validation_expression/similation/HEK293T/"
    "HEK293T.simulation_transcripts.gtf"
)
SORTED_MAPPING_DIR = SAMPLE_DIR / "mapping" / "sorted"
SIMULATION_INFO = TRUTH_GTF.parent / "HEK293T.simulation_transcript_info.tsv"
OUTPUT_DIR = SAMPLE_DIR / ("diagnosis_remaining_" + datetime.now().strftime("%Y%m%d_%H%M%S"))


def contains_id(row, target):
    return any(re.search(r"(?<![A-Za-z0-9])" + re.escape(target) + r"(?:\.\d+)?(?![A-Za-z0-9])",
                         str(value)) for value in row.values())


def group_key(value):
    match = re.match(r"(.+_(?:plus|minus))_(?:T\d+_)?(G\d+)(?:_|$)", value or "")
    return f"{match[1]}_{match[2]}" if match else value


def write_table(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def attributes(text):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', text))


def read_truth_loci(path, output, aliases=None):
    aliases = aliases or {target: {target} for target in TARGETS}
    loci = defaultdict(list)
    with path.open(encoding="utf-8") as handle, output.open("w", encoding="utf-8") as out:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            attrs = attributes(cols[8])
            tid = re.sub(r"^(ENST\d+)\.\d+$", r"\1", attrs.get("transcript_id", ""))
            matches = [target for target in TARGETS
                       if tid in aliases[target] or contains_id(attrs, target)]
            if matches:
                out.write(line)
                for target in matches:
                    loci[target].append(dict(chromosome=cols[0], start=int(cols[3])-1,
                                             end=int(cols[4]), strand=cols[6], feature=cols[2],
                                             gene_id=attrs.get("gene_id", "")))
    return loci


def overlaps_truth(row, loci):
    return any(row["chromosome"] == locus["chromosome"]
               and row["strand"] == locus["strand"]
               and int(row["genome_start"]) < locus["end"]
               and int(row["genome_end"]) > locus["start"] for locus in loci)


def inspect_group_competitors(group_ids, tables, manifest):
    paths = sorted(SORTED_MAPPING_DIR.glob("*.parquet")) if SORTED_MAPPING_DIR else []
    if not paths:
        manifest["competitor_scan"] = {"status": "mapping/sorted unavailable; competitors not reconstructed"}
        return
    import pyarrow.dataset as ds
    hits = defaultdict(dict)
    for path in paths:
        dataset = ds.dataset(str(path), format="parquet")
        scanner = dataset.scanner(columns=["group_id", "read_id", "other_group_ids"],
                                  filter=ds.field("group_id").isin(sorted(group_ids)), batch_size=65536)
        for batch in scanner.to_batches():
            for row in batch.to_pylist():
                candidates = {row["group_id"], *(row["other_group_ids"] or [])}
                hits[row["group_id"]].setdefault(row["read_id"], set()).update(candidates)
    stats = {r["group_id"]: r for r in tables["group_tpm"]}
    metadata_by_group = defaultdict(list)
    for row in tables["metadata"]:
        metadata_by_group[group_key(row["overlap_group_id"])].append(row)
    for group in sorted(group_ids):
        shared = Counter()
        mass = Counter()
        patterns = Counter()
        for candidates in hits[group].values():
            patterns[tuple(sorted(candidates))] += 1
            missing = candidates - stats.keys()
            if missing:
                raise ValueError(f"Candidate groups missing from group TPM: {missing}")
            total = sum(float(stats[g]["final_tpm"]) for g in candidates)
            for candidate in candidates:
                shared[candidate] += 1
                mass[candidate] += (float(stats[candidate]["final_tpm"]) / total
                                    if total > 0 else 1.0 / len(candidates))
        rows = [dict(candidate_group=g, shared_target_group_reads=shared[g],
                     posterior_on_target_group_reads=mass[g],
                     whole_group_expected_count=stats[g]["final_expected_count"],
                     whole_group_unique_reads=stats[g]["unique_read_count"],
                     final_tpm=stats[g]["final_tpm"],
                     gene_names=";".join(sorted({r.get("gene_name", "") for r in metadata_by_group[g]})))
                for g in sorted(shared, key=lambda g: -mass[g])]
        write_table(OUTPUT_DIR / f"{group}.competitors.tsv", rows,
                    ["candidate_group", "shared_target_group_reads", "posterior_on_target_group_reads",
                     "whole_group_expected_count", "whole_group_unique_reads", "final_tpm", "gene_names"])
        write_table(OUTPUT_DIR / f"{group}.candidate_patterns.tsv",
                    [dict(candidate_groups=";".join(k), read_count=v) for k, v in patterns.most_common()],
                    ["candidate_groups", "read_count"])
        selected = [r for g in shared for r in metadata_by_group[g]]
        write_table(OUTPUT_DIR / f"{group}.competitor_metadata.tsv", selected, list(tables["metadata"][0]))
    manifest["competitor_scan"] = {
        "status": "complete", "target_groups": sorted(group_ids),
        "read_counts": {g: len(hits[g]) for g in group_ids},
        "note": "Posterior recalculated from final group TPM on reads hitting target group; not truth-origin tracking.",
        "inputs": [{"path": str(p), "bytes": p.stat().st_size, "mtime": p.stat().st_mtime} for p in paths],
    }


def main():
    inputs = {
        "metadata": SAMPLE_DIR / "reference" / f"{SAMPLE_NAME}.metadata.tsv",
        "native_tpm": SAMPLE_DIR / "em" / f"{SAMPLE_NAME}.transcript_tpm.tsv",
        "source_tpm": SAMPLE_DIR / "em" / f"{SAMPLE_NAME}.transcript_tpm_with_source.tsv",
        "group_tpm": SAMPLE_DIR / "em" / f"{SAMPLE_NAME}.group_tpm.tsv",
        "merge": SAMPLE_DIR / "em" / f"{SAMPLE_NAME}.transcript_merge.tsv",
        "cds_regions": SAMPLE_DIR / "em" / f"{SAMPLE_NAME}.transcript_tpm.cds_regions.tsv",
        "cds_filter": SAMPLE_DIR / "em" / f"{SAMPLE_NAME}.transcript_tpm.cds_filter.tsv",
    }
    tables, fields, manifest = {}, {}, {}
    for name, path in inputs.items():
        if not path.exists():
            if name in {"metadata", "native_tpm", "source_tpm", "group_tpm"}:
                raise FileNotFoundError(path)
            tables[name], fields[name] = [], ["status"]
            manifest[name] = {"path": str(path), "status": "missing"}
            continue
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields[name], tables[name] = reader.fieldnames, list(reader)
        manifest[name] = {"path": str(path), "bytes": path.stat().st_size,
                          "mtime": path.stat().st_mtime}
    OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    aliases = {target: {target} for target in TARGETS}
    if SIMULATION_INFO.exists():
        with SIMULATION_INFO.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            info_fields, info_rows = reader.fieldnames, list(reader)
        selected_info = []
        for row in info_rows:
            matches = [target for target in TARGETS if contains_id(row, target)]
            if matches:
                selected_info.append(row)
                for target in matches:
                    for key, value in row.items():
                        if "transcript_id" in key and value:
                            # STRG.3744.1 is an ID, not an Ensembl version.
                            for token in re.split(r"[,;|\s]+", value):
                                aliases[target].add(re.sub(r"^(ENST\d+)\.\d+$", r"\1", token))
        write_table(OUTPUT_DIR / "simulation_target_info.tsv", selected_info, info_fields)
        manifest["simulation_info"] = {"path": str(SIMULATION_INFO),
                                       "target_aliases": {k: sorted(v) for k, v in aliases.items()}}
    else:
        manifest["simulation_info"] = {"path": str(SIMULATION_INFO), "status": "missing"}
    loci = read_truth_loci(TRUTH_GTF, OUTPUT_DIR / "truth_targets.gtf", aliases) if TRUTH_GTF else {}
    if TRUTH_GTF:
        manifest["truth_gtf"] = {"path": str(TRUTH_GTF), "bytes": TRUTH_GTF.stat().st_size,
                                 "mtime": TRUTH_GTF.stat().st_mtime,
                                 "matching_lines": {t: len(loci.get(t, [])) for t in TARGETS}}
    all_groups = set()
    selected_reference_ids = set()
    summary = []
    for target in TARGETS:
        direct = [r for r in tables["metadata"] if contains_id(r, target)]
        source = [r for r in tables["source_tpm"] if contains_id(r, target)]
        overlapping = [r for r in tables["metadata"] if overlaps_truth(r, loci.get(target, []))]
        groups = {group_key(r.get("overlap_group_id", "")) for r in direct + source + overlapping}
        groups.discard("")
        if target == "ENST00000647841":
            groups.add("chr15_minus_G00926")  # Previously observed group; verify against current metadata.
        all_groups.update(groups)
        selected_reference_ids.update(r["reference_id"] for r in tables["metadata"]
                                      if group_key(r["overlap_group_id"]) in groups)
        note = "reference ID found" if direct else "ID absent from metadata; structural absence NOT established"
        summary.append(dict(target=target, metadata_rows=len(direct), final_source_rows=len(source),
                            truth_gtf_rows=len(loci.get(target, [])),
                            locus_overlap_rows=len(overlapping), groups=";".join(sorted(groups)), status=note))
        write_table(OUTPUT_DIR / f"{target}.locus_overlap_metadata.tsv", overlapping, fields["metadata"])
        write_table(OUTPUT_DIR / f"{target}.direct_metadata.tsv", direct, fields["metadata"])
        for name, rows in tables.items():
            selected = [r for r in rows if contains_id(r, target) or any(
                group_key(r.get(key, "")) in groups
                for key in ("group_id", "overlap_group_id", "reference_id", "transcript_id"))]
            write_table(OUTPUT_DIR / f"{target}.{name}.tsv", selected, fields[name])
    reference_gtf = SAMPLE_DIR / "reference" / f"{SAMPLE_NAME}.gtf"
    if reference_gtf.exists():
        with reference_gtf.open(encoding="utf-8") as handle, (OUTPUT_DIR / "reference_candidates.gtf").open("w", encoding="utf-8") as out:
            for line in handle:
                cols = line.rstrip("\n").split("\t")
                if len(cols) == 9 and attributes(cols[8]).get("transcript_id") in selected_reference_ids:
                    out.write(line)
        manifest["reference_gtf"] = {"path": str(reference_gtf), "mtime": reference_gtf.stat().st_mtime}
    else:
        manifest["reference_gtf"] = {"path": str(reference_gtf), "status": "missing"}
    inspect_group_competitors(all_groups, tables, manifest)
    write_table(OUTPUT_DIR / "summary.tsv", summary, list(summary[0]))
    (OUTPUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Output: {OUTPUT_DIR}")
    for row in summary:
        print(row)
    print("Share this output directory. A missing ID alone does not prove a missing reference sequence.")


if __name__ == "__main__":
    main()
