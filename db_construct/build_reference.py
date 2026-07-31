from __future__ import annotations

from pathlib import Path

from db_construct.coordinate_writer import (
    write_coordinate_table,
)
from db_construct.fasta_writer import (
    write_reference_fasta,
)
from db_construct.gtf_reader import (
    load_gencode_gff3,
    load_stringtie_transcripts,
)
from db_construct.gtf_writer import (
    write_reference_gtf,
)
from db_construct.metadata_writer import (
    write_reference_metadata,
)
from db_construct.overlap_group import (
    build_reference_transcripts,
)
from db_construct.transcript_sort import (
    sort_transcripts,
)


def build_transcript_reference(
    gencode_gff3: str,
    stringtie_gtf: str,
    genome_fasta: str,
    output_prefix: str,
):
    """Build a transcript-only reference for RiboTxAlign."""

    output_prefix_path = Path(
        output_prefix
    )
    output_prefix_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("1. Read GENCODE GFF3.")

    gencode_transcript_dict = (
        load_gencode_gff3(
            gencode_gff3
        )
    )

    print(
        "GENCODE transcripts:",
        len(gencode_transcript_dict),
    )

    print("2. Read StringTie2 GTF.")

    stringtie_transcript_dict = (
        load_stringtie_transcripts(
            stringtie_gtf
        )
    )

    print(
        "StringTie2 transcripts:",
        len(stringtie_transcript_dict),
    )

    print("3. Sort transcript annotations.")

    gencode_transcripts = sort_transcripts(
        gencode_transcript_dict
    )

    stringtie_transcripts = sort_transcripts(
        stringtie_transcript_dict
    )

    print(
        "4. Build resolved reference transcripts."
    )

    reference_transcripts = (
        build_reference_transcripts(
            gencode_transcripts=(
                gencode_transcripts
            ),
            stringtie_transcripts=(
                stringtie_transcripts
            ),
        )
    )

    output_gtf = Path(
        f"{output_prefix}.gtf"
    )
    output_coordinate_tsv = Path(
        f"{output_prefix}.coordinates.tsv"
    )
    output_fasta = Path(
        f"{output_prefix}.fa"
    )
    output_metadata_tsv = Path(
        f"{output_prefix}.metadata.tsv"
    )

    print("5. Write reference GTF.")

    write_reference_gtf(
        reference_transcripts=(
            reference_transcripts
        ),
        output_path=output_gtf,
    )

    print("6. Write coordinate table.")

    write_coordinate_table(
        reference_transcripts=(
            reference_transcripts
        ),
        output_path=output_coordinate_tsv,
    )

    print("7. Write reference FASTA.")

    write_reference_fasta(
        reference_transcripts=(
            reference_transcripts
        ),
        genome_fasta=genome_fasta,
        output_path=output_fasta,
    )

    print("8. Write metadata table.")

    write_reference_metadata(
        reference_transcripts=(
            reference_transcripts
        ),
        output_path=output_metadata_tsv,
    )

    print(
        "Reference construction flow completed."
    )

    return reference_transcripts
