"""Prefer primary chromosomes over equal-best ALT hits as a reporting policy."""
import re
from functools import lru_cache


def chromosome(reference_name):
    return re.split(r"_(?:plus|minus)_", reference_name, maxsplit=1)[0]


@lru_cache(maxsize=65_536)
def _reference_flags(reference_name):
    """Cache only immutable reference properties, never alignment scores."""
    chrom = chromosome(reference_name)
    return bool(re.fullmatch(r"chr(?:[1-9][0-9]*|X|Y|M|MT)", chrom)), chrom.endswith("_alt")


def prefer_primary_alt_ties(alignments, enabled=True, limit_reached=False):
    """Return retained alignments and removed ALT hit count.

    The mapper normally supplies best-score ties only. Check scores explicitly
    so a superior ALT is never removed. A truncated candidate set is left alone.
    Non-ALT scaffolds remain candidates; this is not an ALT/primary locus liftover.
    """
    if not enabled or limit_reached or not alignments:
        return alignments, 0
    best_score = max(a.score for a in alignments)
    has_primary = any(
        a.score == best_score
        and _reference_flags(a.reference_name)[0]
        for a in alignments
    )
    if not has_primary:
        return alignments, 0
    retained = [a for a in alignments
                if not (a.score == best_score and _reference_flags(a.reference_name)[1])]
    return retained, len(alignments) - len(retained)
