from __future__ import annotations

from db_construct.data_class import Transcript


def sort_transcripts(
    transcripts: dict[str, Transcript],
) -> list[Transcript]:
    """Sort transcripts by chromosome, strand, and genomic position."""

    return sorted(
        transcripts.values(),
        key=lambda transcript: (
            transcript.chrom,
            transcript.strand,
            transcript.start
            if transcript.start is not None
            else -1,
            transcript.end
            if transcript.end is not None
            else -1,
            transcript.transcript_id,
        ),
    )
