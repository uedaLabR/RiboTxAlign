from __future__ import annotations

from pathlib import Path

import pysam

from db_construct.data_class import Transcript
from db_construct.decoy_fasta import (
    write_decoy_fasta_records,
)
from db_construct.fasta_utils import (
    write_fasta_record,
)
from db_construct.gtf_writer import (
    get_output_transcript_id,
    normalize_exons,
)


_COMPLEMENT_TABLE = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def reverse_complement(
    sequence: str,
) -> str:
    """Return the reverse-complemented nucleotide sequence."""

    return sequence.translate(
        _COMPLEMENT_TABLE
    )[::-1]


def get_exons_in_transcript_order(
    transcript: Transcript,
) -> list[tuple[int, int]]:
    """Return valid exons in biological 5'-to-3' order."""

    exons = normalize_exons(
        transcript
    )

    if transcript.strand == "-":
        exons.reverse()

    return exons


def fetch_transcript_sequence(
    transcript: Transcript,
    genome: pysam.FastaFile,
) -> str:
    """Fetch and splice one transcript sequence from the genome FASTA."""

    sequence_parts: list[str] = []

    for exon_start, exon_end in get_exons_in_transcript_order(
        transcript
    ):
        sequence = genome.fetch(
            transcript.chrom,
            exon_start,
            exon_end,
        ).upper()

        if transcript.strand == "-":
            sequence = reverse_complement(
                sequence
            )

        sequence_parts.append(
            sequence
        )

    return "".join(
        sequence_parts
    )


def write_reference_fasta(
    reference_transcripts: list[Transcript],
    genome_fasta: str | Path,
    output_path: str | Path,
    decoy_directory: str | Path | None = None,
) -> None:
    """Write spliced transcript sequences followed by decoy sequences."""

    genome_fasta = Path(
        genome_fasta
    )
    output_path = Path(
        output_path
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not genome_fasta.exists():
        raise FileNotFoundError(
            f"Genome FASTA was not found: {genome_fasta}"
        )

    written_transcript_count = 0
    written_decoy_count = 0
    skipped_without_exons = 0
    transcript_bases = 0
    decoy_bases = 0
    seen_reference_ids: set[str] = set()

    sorted_transcripts = sorted(
        reference_transcripts,
        key=lambda transcript: (
            transcript.chrom,
            transcript.strand,
            transcript.start
            if transcript.start is not None
            else -1,
            transcript.end
            if transcript.end is not None
            else -1,
            get_output_transcript_id(
                transcript
            ),
        ),
    )

    fasta_index = Path(
        f"{genome_fasta}.fai"
    )

    if not fasta_index.exists():
        print(
            "Create genome FASTA index:",
            fasta_index,
        )
        pysam.faidx(
            str(genome_fasta)
        )

    with pysam.FastaFile(
        str(genome_fasta)
    ) as genome, output_path.open(
        "w",
        encoding="utf-8",
    ) as output_file:

        available_chromosomes = set(
            genome.references
        )

        for transcript in sorted_transcripts:
            reference_id = (
                get_output_transcript_id(
                    transcript
                )
            )

            if reference_id in seen_reference_ids:
                raise ValueError(
                    "Duplicate reference ID in FASTA output: "
                    f"{reference_id}"
                )

            seen_reference_ids.add(
                reference_id
            )

            exons = get_exons_in_transcript_order(
                transcript
            )

            if not exons:
                skipped_without_exons += 1
                continue

            if transcript.chrom not in available_chromosomes:
                raise KeyError(
                    "Chromosome is absent from genome FASTA: "
                    f"{transcript.chrom}"
                )

            sequence = fetch_transcript_sequence(
                transcript=transcript,
                genome=genome,
            )

            expected_length = sum(
                exon_end - exon_start
                for exon_start, exon_end in exons
            )

            if len(sequence) != expected_length:
                raise RuntimeError(
                    "Transcript sequence length mismatch: "
                    f"{reference_id}; expected={expected_length}; "
                    f"observed={len(sequence)}"
                )

            write_fasta_record(
                output_file=output_file,
                reference_id=reference_id,
                sequence=sequence,
            )

            written_transcript_count += 1
            transcript_bases += len(sequence)

        (
            written_decoy_count,
            decoy_bases,
        ) = write_decoy_fasta_records(
            output_file=output_file,
            decoy_directory=decoy_directory,
            seen_reference_ids=seen_reference_ids,
        )

    total_record_count = (
        written_transcript_count
        + written_decoy_count
    )
    total_bases = (
        transcript_bases
        + decoy_bases
    )

    print(
        "Reference FASTA written:",
        output_path,
    )
    print(
        "FASTA transcripts:",
        written_transcript_count,
        "transcript_bases:",
        transcript_bases,
        "decoy_records:",
        written_decoy_count,
        "decoy_bases:",
        decoy_bases,
        "total_records:",
        total_record_count,
        "total_bases:",
        total_bases,
        "skipped_without_exons:",
        skipped_without_exons,
        "bytes:",
        output_path.stat().st_size,
    )

    if written_transcript_count == 0:
        raise RuntimeError(
            "The output FASTA contains no transcript sequences."
        )
