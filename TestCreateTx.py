from db_construct.build_reference import build_common_gene_references
from db_construct.calculate_mappability import calculate_reference_mappability


GENCODE_GFF3 = (
    "/mnt/ssdnas07/pipeline/rna_modscore_v2/source/"
    "gencode.v44.annotation.gff3"
)

# Optional. Set to None to build the reference from GENCODE only.
STRINGTIE_GTF = None

# Example:
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
    "HEK293T_RiboTxAlgin"
)


gene_groups, novel_transcripts = build_common_gene_references(
    gencode_gff3=GENCODE_GFF3,
    stringtie_gtf=STRINGTIE_GTF,
    genome_fasta=GENOME_FASTA,
    output_prefix=OUTPUT_PREFIX,
    flank_length=40,
)

print("Gene groups:", len(gene_groups))
print(
    "StringTie transcripts not linked to GENCODE:",
    len(novel_transcripts),
)


# Calculate exact-match mappability after the reference FASTA is created.
MAPPABILITY_READ_LENGTH = 30
MAPPABILITY_STEP = 1
MAPPABILITY_MAXIMUM_HITS = 100
MAPPABILITY_PROGRESS_GENE_INTERVAL = 100

summary_df, _ = calculate_reference_mappability(
    reference_fasta=f"{OUTPUT_PREFIX}.fa",
    coordinate_tsv=f"{OUTPUT_PREFIX}.coordinates.tsv",
    output_prefix=OUTPUT_PREFIX,
    read_length=MAPPABILITY_READ_LENGTH,
    step=MAPPABILITY_STEP,
    maximum_hits=MAPPABILITY_MAXIMUM_HITS,
    progress_gene_interval=MAPPABILITY_PROGRESS_GENE_INTERVAL,
    write_position_table=False,
)

print("Mappability references:", len(summary_df))
