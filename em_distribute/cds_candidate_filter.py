"""Resolve CDS-like regions before transcript classification and merging.

All intervals are zero-based and half-open. Group fallback projects the
group's genomic CDS start/end envelope onto the target's own exons. No P-site
assumption is made. Read fallback is used only for EM length, never as CDS
evidence. Non-CDS hits are dropped only when the read has a CDS hit.
"""
from collections import defaultdict


def valid_interval(start, end):
    return start is not None and end is not None and 0 <= start < end


def project_genomic_span(definition, start, end):
    intervals = []
    offset = 0
    reverse = definition["strand"] == "-"
    for left, right in sorted(set(definition["exons"]), reverse=reverse):
        a, b = max(left, start), min(right, end)
        if a < b:
            intervals.append((offset + right - b, offset + right - a)
                             if reverse else (offset + a - left, offset + b - left))
        offset += right - left
    return intervals


def own_region(definition):
    start, end = definition.get("cds_reference_start"), definition.get("cds_reference_end")
    if valid_interval(start, end):
        if end > definition["transcript_length"]:
            raise ValueError(f"CDS exceeds transcript length: {definition['reference_id']}")
        return [(start, end)]
    start, end = definition.get("cds_start"), definition.get("cds_end")
    if valid_interval(start, end):
        return project_genomic_span(definition, start, end)
    return []


def genomic_region(definition):
    """Recover genomic CDS bounds also for metadata with only spliced bounds."""
    region = own_region(definition)
    offset = 0
    genomic = []
    reverse = definition["strand"] == "-"
    for left, right in sorted(set(definition["exons"]), reverse=reverse):
        for start, end in region:
            a, b = max(start, offset), min(end, offset + right - left)
            if a < b:
                genomic.append((right - (b - offset), right - (a - offset))
                               if reverse else (left + a - offset, left + b - offset))
        offset += right - left
    return genomic


def build_group_cds_spans(definitions):
    spans = defaultdict(list)
    for definition in definitions.values():
        key = (definition["group_id"], definition["chromosome"], definition["strand"])
        spans[key].extend(genomic_region(definition))
    return {key: (min(a for a, b in values), max(b for a, b in values))
            for key, values in spans.items() if values}


def resolve_cds_region(definition, group_spans):
    """Prefer a transcript's own CDS; otherwise project its group's bounds."""
    region = own_region(definition)
    if region:
        return region, "own_cds"
    key = (definition["group_id"], definition["chromosome"], definition["strand"])
    if key in group_spans:
        return project_genomic_span(definition, *group_spans[key]), "group_cds"
    return [], "unresolved"


def filter_cds_candidates(records, definitions, group_spans):
    by_transcript = defaultdict(list)
    for record in records:
        by_transcript[record["transcript_id"]].append(record)
    diagnostics = []
    cds_regions = {}
    for tid, hits in by_transcript.items():
        definition = definitions[tid]
        region, source = resolve_cds_region(definition, group_spans)
        cds_regions[tid] = region
        if not region:
            spans = [(max(0, int(r["reference_start"])),
                      min(definition["transcript_length"], int(r["reference_end"]))) for r in hits]
            spans = [(a, b) for a, b in spans if a < b]
            region = [(min(a for a, b in spans), max(b for a, b in spans))] if spans else []
            source = "read_span" if region else "unresolved"
        definition["em_region_length"] = sum(b - a for a, b in region)
        diagnostics.append(dict(transcript_id=tid, region_source=source,
                                region_intervals=";".join(f"{a}-{b}" for a, b in region),
                                em_region_length=definition["em_region_length"]))
    cds_hits = [r for r in records if any(
        int(r["reference_start"]) < b and int(r["reference_end"]) > a
        for a, b in cds_regions[r["transcript_id"]])]
    cds_reads = {r["read_id"] for r in cds_hits}
    cds_hit_ids = {id(r) for r in cds_hits}
    retained = [r for r in records if r["read_id"] not in cds_reads or id(r) in cds_hit_ids]
    retained_reads = {r["read_id"] for r in retained}
    lost = {r["read_id"] for r in records} - retained_reads
    return retained, diagnostics, lost
