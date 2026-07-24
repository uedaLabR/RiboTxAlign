from __future__ import annotations

from pathlib import Path

from db_construct.create_fasta import (
    write_coordinate_table,
    write_reference_fasta,
    write_transcript_path_table,
)
from db_construct.gtf_reader import load_gencode_gff3, load_stringtie_transcripts
from db_construct.merge_transcript import merge_transcripts_by_gene
from db_construct.resolve_each_gene import resolve_all_genes


def build_common_gene_references(
    gencode_gff3: str,
    stringtie_gtf: str | None,
    genome_fasta: str,
    output_prefix: str,
    flank_length: int = 40,
):
    output_prefix_path = Path(output_prefix)
    output_prefix_path.parent.mkdir(parents=True, exist_ok=True)

    gencode_transcripts = load_gencode_gff3(gencode_gff3)
    stringtie_transcripts = (
        load_stringtie_transcripts(stringtie_gtf)
        if stringtie_gtf
        else {}
    )

    gene_groups, novel_transcripts = merge_transcripts_by_gene(
        gencode_transcripts=gencode_transcripts,
        stringtie_transcripts=stringtie_transcripts,
    )

    resolve_all_genes(gene_groups)

    output_fasta = f"{output_prefix}.fa"
    output_coordinate_tsv = f"{output_prefix}.coordinates.tsv"
    output_transcript_path_tsv = f"{output_prefix}.transcript_paths.tsv"

    coordinate_records = write_reference_fasta(
        gene_groups=gene_groups,
        genome_fasta=genome_fasta,
        output_fasta=output_fasta,
        flank_length=flank_length,
    )
    write_coordinate_table(coordinate_records, output_coordinate_tsv)
    write_transcript_path_table(gene_groups, output_transcript_path_tsv)

    return gene_groups, novel_transcripts
