from __future__ import annotations

from dataclasses import dataclass, field


Interval = tuple[int, int]
TranscriptPath = tuple[int, ...]
Junction = tuple[int, int]


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str | None
    chrom: str
    strand: str
    source: str
    start: int = 0
    end: int = 0
    tpm: float = 0.0
    reference_id: str | None = None
    reference_gene_id: str | None = None
    gene_name: str | None = None
    exons: list[Interval] = field(default_factory=list)
    cds_regions: list[Interval] = field(default_factory=list)

    @property
    def exonic_length(self) -> int:
        return sum(end - start for start, end in self.exons)

    @property
    def cds_length(self) -> int:
        return sum(end - start for start, end in self.cds_regions)

    def ordered_exons(self) -> list[Interval]:
        exons = sorted(self.exons)
        if self.strand == "-":
            exons.reverse()
        return exons


@dataclass
class ExonBlock:
    block_id: int
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class GeneGroup:
    gene_id: str
    chrom: str
    strand: str
    gene_name: str | None = None
    transcripts: dict[str, Transcript] = field(default_factory=dict)
    exon_blocks: list[ExonBlock] = field(default_factory=list)
    transcript_paths: dict[str, TranscriptPath] = field(default_factory=dict)
    transcript_exon_blocks: dict[str, tuple[tuple[int, ...], ...]] = field(default_factory=dict)
    transcript_junctions: dict[str, tuple[Junction, ...]] = field(default_factory=dict)
    path_tpm: dict[TranscriptPath, float] = field(default_factory=dict)
    gene_tpm: float = 0.0
    representative_path: TranscriptPath | None = None
    representative_transcript_id: str | None = None
    representative_path_tpm: float = 0.0
    alternative_block_ids: set[int] = field(default_factory=set)
    alternative_junctions: set[Junction] = field(default_factory=set)

    def add_transcript(self, transcript: Transcript) -> None:
        self.transcripts[transcript.transcript_id] = transcript
        if self.gene_name is None and transcript.gene_name:
            self.gene_name = transcript.gene_name

    @property
    def start(self) -> int:
        starts = [t.start for t in self.transcripts.values() if t.end > t.start]
        return min(starts) if starts else 0

    @property
    def end(self) -> int:
        ends = [t.end for t in self.transcripts.values() if t.end > t.start]
        return max(ends) if ends else 0

    def get_block(self, block_id: int) -> ExonBlock:
        for block in self.exon_blocks:
            if block.block_id == block_id:
                return block
        raise KeyError(f"Block {block_id} was not found in gene {self.gene_id}")


@dataclass
class ReferenceMapRecord:
    reference_id: str
    gene_id: str
    record_type: str
    segment_index: int
    reference_start: int
    reference_end: int
    chrom: str
    genome_start: int
    genome_end: int
    strand: str
    transcript_id: str | None = None
    block_id: int | None = None
    junction_left_block: int | None = None
    junction_right_block: int | None = None
