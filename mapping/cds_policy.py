"""Annotate and prefer CDS-overlapping best-score mapping candidates."""
from pathlib import Path
import warnings

from em_distribute.cds_candidate_filter import build_group_cds_spans, resolve_cds_region
from em_distribute.reference_utils import load_reference_definitions


def load_mapping_cds_regions(reference_path, metadata_tsv=None, reference_gtf=None):
    """Load explicit annotation paths, or discover the reference's sidecars."""
    prefix = Path(reference_path).with_suffix("")
    explicit = metadata_tsv is not None or reference_gtf is not None
    metadata = Path(metadata_tsv) if metadata_tsv is not None else Path(f"{prefix}.metadata.tsv")
    gtf = Path(reference_gtf) if reference_gtf is not None else Path(f"{prefix}.gtf")
    if not metadata.is_file() or not gtf.is_file():
        if explicit:
            raise FileNotFoundError(f"CDS annotation requires both {metadata} and {gtf}")
        warnings.warn(f"CDS preference disabled: missing {metadata} or {gtf}", stacklevel=2)
        return {}
    definitions = load_reference_definitions(metadata, gtf)
    group_spans = build_group_cds_spans(definitions)
    return {tid: resolve_cds_region(definition, group_spans)
            for tid, definition in definitions.items()}


def prefer_cds_hits(alignments, cds_regions):
    """Keep all candidates if none overlaps CDS; otherwise keep only CDS hits."""
    for alignment in alignments:
        intervals, source = cds_regions.get(alignment.reference_name, ([], "unresolved"))
        alignment.cds_source = source
        alignment.is_cds = any(alignment.reference_start < end and alignment.reference_end > start
                               for start, end in intervals)
    cds_hits = [alignment for alignment in alignments if alignment.is_cds]
    return cds_hits or alignments
