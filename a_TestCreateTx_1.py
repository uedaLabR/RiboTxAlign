from db_construct.build_reference import build_transcript_reference


GENCODE_GFF3 = (
    "/mnt/ssdnas07/pipeline/rna_modscore_v2/source/"
    "gencode.v44.annotation.gff3"
)

STRINGTIE_GTF = (
    "/mnt/share/project/riboML/stringtie2/"
    "HEK293T.transcripts.gtf"
)

GENOME_FASTA = (
    "/mnt/ssdnas07/pipeline/rna_modscore_v2/source/"
    "genome/hg38.fa"
)

OUTPUT_PREFIX = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign"
)


reference_transcripts = build_transcript_reference(
    gencode_gff3=GENCODE_GFF3,
    stringtie_gtf=STRINGTIE_GTF,
    genome_fasta=GENOME_FASTA,
    output_prefix=OUTPUT_PREFIX,
)

print(
    "Reference transcripts:",
    len(reference_transcripts),
)
