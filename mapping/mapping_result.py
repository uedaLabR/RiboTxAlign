from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class RawAlignment:
    """One filtered alignment against the transcript reference."""

    read_id: str
    reference_name: str
    reference_start: int
    reference_end: int
    query_start: int
    query_end: int
    strand: int
    matching_bases: int
    alignment_length: int
    edit_distance: int
    cigartuples: list[tuple[int, int]]
    is_cds: bool = False
    cds_source: str = "unresolved"

    @property
    def score(self) -> int:
        """Return the internal score used only to select best ties."""

        return self.matching_bases - self.edit_distance

    @property
    def identity(self) -> float:
        if self.alignment_length <= 0:
            return 0.0
        return self.matching_bases / self.alignment_length

    def query_coverage(self, read_length: int) -> float:
        if read_length <= 0:
            return 0.0
        return (self.query_end - self.query_start) / read_length
