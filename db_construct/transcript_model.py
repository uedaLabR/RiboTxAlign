from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ExonModel:
    """One exon interval in one-based inclusive genomic coordinates."""

    start: int
    end: int
    exon_number: int | None = None


@dataclass
class TranscriptModel:
    """Transcript annotation loaded from GENCODE or StringTie2."""

    transcript_id: str
    gene_id: str
    chromosome: str
    strand: str
    exons: list[ExonModel]
    source: str
    gene_name: str | None = None
    cds_intervals: list[tuple[int, int]] = field(default_factory=list)

    @property
    def genomic_start(self) -> int:
        return min(exon.start for exon in self.exons)

    @property
    def genomic_end(self) -> int:
        return max(exon.end for exon in self.exons)


@dataclass
class OverlapGroup:
    """Connected group of transcripts sharing exonic regions."""

    overlap_group_id: str
    transcripts: list[TranscriptModel]


@dataclass
class ReferenceTranscript:
    """Resolved transcript model written to the final reference."""

    reference_id: str
    overlap_group_id: str
    gene_id: str
    chromosome: str
    strand: str
    exons: list[ExonModel]
    representative_transcript_id: str
    source_transcript_ids: list[str]
    gene_name: str | None = None
    cds_intervals: list[tuple[int, int]] = field(default_factory=list)
    gencode_transcript_ids: list[str] = field(default_factory=list)
    stringtie_transcript_ids: list[str] = field(default_factory=list)
    resolution_type: str = "not_implemented"

    @property
    def transcript_length(self) -> int:
        return sum(exon.end - exon.start + 1 for exon in self.exons)
