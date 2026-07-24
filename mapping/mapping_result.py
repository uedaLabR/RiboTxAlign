from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class RawAlignment:
    read_id: str
    source: str
    reference_name: str
    reference_start: int
    reference_end: int
    query_start: int
    query_end: int
    strand: int
    mapq: int
    matching_bases: int
    alignment_length: int
    edit_distance: int
    is_primary: bool
    cigartuples: list[tuple[int, int]]

    @property
    def score(self) -> int:
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

    def is_full_length(self, read_length: int) -> bool:
        return self.query_start == 0 and self.query_end == read_length


@dataclass
class ProjectedAlignment:
    read_id: str
    chrom: str
    reference_start: int
    cigartuples: list[tuple[int, int]]
    is_reverse: bool
    mapq: int
    edit_distance: int
    score: int
    source_reference: str
    source_reference_start: int
    source_reference_end: int
    gene_id: str = ''
    record_type: str = ''
    support_sources: set[str] = field(default_factory=set)
    assignment_weight: float = 1.0

    @property
    def genomic_blocks(self) -> tuple[tuple[int, int], ...]:
        """Return normalized genomic blocks covered by this alignment.

        Match, mismatch, deletion and padding operations remain in the same
        genomic block. An intron/skipped-reference operation (N) starts a new
        block. Query-only clipping and insertion operations do not consume the
        reference. This representation makes alignments from different
        synthetic transcript references comparable after genomic projection.
        """

        reference_position = self.reference_start
        block_start = reference_position
        blocks: list[tuple[int, int]] = []

        for operation, length in self.cigartuples:
            if length <= 0:
                continue

            if operation == 3:  # N: skipped reference / intron
                if reference_position > block_start:
                    blocks.append((block_start, reference_position))
                reference_position += length
                block_start = reference_position
                continue

            if operation in {0, 2, 7, 8}:  # M, D, =, X
                reference_position += length

        if reference_position > block_start:
            blocks.append((block_start, reference_position))

        return tuple(blocks)

    @property
    def genomic_placement_key(self) -> tuple:
        """Return a key for one genomic start position.

        Different transcript references and CIGAR representations are
        collapsed when they project to the same chromosome, genomic start
        position and strand.
        """

        return (
            self.chrom,
            int(self.reference_start),
            bool(self.is_reverse),
        )

    @property
    def deduplication_key(self) -> tuple:
        # Backward-compatible alias used by the mapping modules.
        return self.genomic_placement_key
