from mapping.map_read_mt import run_transcript_mapping_to_parquet_mt
# For single-process debugging, use:
# from mapping.map_reads import run_transcript_mapping_to_parquet


RIBO_FASTQ = (
    "/mnt/share/project/riboML/riboseq_raw/trimmed_26_34/deduplicated/SRR26942051.26_34.deduplicated.fastq.gz"
)

TRANSCRIPT_FASTA = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.fa"
)
MAPPING_PARTS_DIRECTORY = (
    "/mnt/share/project/riboML/sample_ref/"
    "mapping_parts"
)
UNMAPPED_FASTQ = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.unmapped.fastq.gz"
)

STATS_JSON = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.mapping_stats.json"
)

LOG_FILE = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.mapping_mt.log"
)

if __name__ == "__main__":
    run_transcript_mapping_to_parquet_mt(
        fastq_path=RIBO_FASTQ,
        transcript_fasta=TRANSCRIPT_FASTA,
        mapping_parts_directory=(
            MAPPING_PARTS_DIRECTORY
        ),
        maximum_reads=None,
        maximum_candidates=16,
        kmer_size=11,
        window_size=5,
        minimum_minimizers=2,
        minimum_chain_score=15,
        minimum_dp_score=15,
        minimum_query_coverage=0.85,
        minimum_identity=0.90,
        worker_count=6,
        reads_per_worker=500000,
        progress_interval=100000,
        parquet_row_group_size=100000,
        log_file=LOG_FILE,
        unmapped_fastq_path=UNMAPPED_FASTQ,
        stats_json_path=STATS_JSON,
    )

