from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Transcript:
    """Transcript annotation used during reference construction."""

    transcript_id: str
    gene_id: str | None
    chrom: str
    strand: str
    source: str

    start: int | None = None
    end: int | None = None
    gene_name: str | None = None

    exons: list[tuple[int, int]] = field(
        default_factory=list
    )
    cds_regions: list[tuple[int, int]] = field(
        default_factory=list
    )

    reference_id: str | None = None
    reference_gene_id: str | None = None
    tpm: float = 0.0

    overlap_group_id: str | None = None
    representative_transcript_id: str | None = None
    source_transcript_ids: list[str] = field(
        default_factory=list
    )
    gencode_transcript_ids: list[str] = field(
        default_factory=list
    )
    stringtie_transcript_ids: list[str] = field(
        default_factory=list
    )
    resolution_type: str = "unresolved"
    # Annotation aliases with identical genomic CDS, not observed source members.
    cds_equivalent_transcript_ids: list[str] = field(default_factory=list)

    @property
    def transcript_length(self) -> int:
        """Return the total spliced exon length."""

        return sum(
            end - start
            for start, end in self.exons
        )
