from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_mapping_stats(path: str | Path) -> dict[str, Any]:
    """Load mapping statistics from JSON."""

    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(
            f"Mapping statistics file was not found: {source}"
        )

    with source.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError(
            f"Mapping statistics must be a JSON object: {source}"
        )

    return data


def write_mapping_stats(
    path: str | Path,
    stats: dict[str, Any],
) -> Path:
    """Write mapping statistics as formatted JSON."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)

    temporary = destination.with_suffix(
        destination.suffix + ".tmp"
    )

    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(
            stats,
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")

    temporary.replace(destination)
    return destination


def percentage(
    numerator: int,
    denominator: int,
) -> float:
    """Return a percentage, using zero when the denominator is zero."""

    if denominator <= 0:
        return 0.0

    return 100.0 * numerator / denominator


def validate_assignment_bam_counts(stats, assigned_reads, main_reads, rescue_reads, final_reads):
    """Validate BAM against post-EM assignments, preserving mapping-stage totals."""
    mapped = int(stats["transcript_mapped"])
    rescued = int(stats["genome_rescued_unique"])
    if not 0 <= assigned_reads <= mapped:
        raise ValueError(f"Assignment count outside mapping total: assignments={assigned_reads}, mapping={mapped}")
    if main_reads != assigned_reads:
        raise ValueError(f"Main BAM count does not match final assignments: BAM={main_reads}, assignments={assigned_reads}")
    if rescue_reads != rescued:
        raise ValueError(f"Rescue BAM count does not match rescue statistics: BAM={rescue_reads}, statistics={rescued}")
    if int(stats["final_mapped"]) != mapped + rescued:
        raise ValueError("Final mapping statistics are inconsistent with transcript mapping and rescue counts")
    if final_reads != assigned_reads + rescued:
        raise ValueError(f"Final BAM count does not match assignments plus rescue: BAM={final_reads}, expected={assigned_reads + rescued}")
    return {
        "transcript_assigned_reads": assigned_reads,
        "transcript_mapped_without_final_assignment": mapped - assigned_reads,
        "expected_final_bam_reads": assigned_reads + rescued,
    }