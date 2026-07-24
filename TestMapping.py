from mapping.map_read_mt import run_dual_mapping_to_bam_mt


RIBO_FASTQ = (
    "/mnt/share/project/riboML/riboseq_raw/trimmed_26_34/"
    "SRR26942051.26_34.fastq.gz"
)
GENOME_FASTA = (
    "/mnt/ssdnas07/pipeline/rna_modscore_v2/source/"
    "genome/hg38.fa"
)
TRANSCRIPT_FASTA = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlgin.fa"
)
COORDINATE_TSV = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlgin.coordinates.tsv"
)

UNIQUE_BAM = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlgin.unique.bam"
)
MULTIPLEMAP_BAM = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlgin.multiplemap.bam"
)
MAPPING_PARQUET = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlgin.multiplemap.parquet"
)


if __name__ == "__main__":
    run_dual_mapping_to_bam_mt(
        fastq_path=RIBO_FASTQ,
        genome_fasta=GENOME_FASTA,
        transcript_fasta=TRANSCRIPT_FASTA,
        coordinate_tsv=COORDINATE_TSV,
        unique_bam=UNIQUE_BAM,
        multiplemap_bam=MULTIPLEMAP_BAM,
        mapping_parquet=MAPPING_PARQUET,
        maximum_reads=None,
        transcript_candidate_count=8,
        kmer_size=11,
        window_size=5,
        minimum_minimizers=2,
        minimum_chain_score=15,
        minimum_dp_score=15,
        minimum_query_coverage=0.85,
        minimum_identity=0.90,
        worker_count=6,
        reads_per_worker=50000,
        sort_threads=8,
        progress_interval=100000,
        parquet_batch_size=100000,
        log_file=(
            "/mnt/share/project/riboML/sample_ref/"
            "HEK293T_RiboTxAlgin.mapping_mt.log"
        ),
    )
