from em_distribute.run_redistribute import run_redistribution


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
MAPPABILITY_PARQUET = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlgin.mappability.parquet"
)
OUTPUT_PREFIX = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlgin.redistributed"
)


if __name__ == "__main__":
    result = run_redistribution(
        coordinate_tsv=COORDINATE_TSV,
        unique_bam=UNIQUE_BAM,
        multiplemap_bam=MULTIPLEMAP_BAM,
        mapping_parquet=MAPPING_PARQUET,
        mappability_parquet=MAPPABILITY_PARQUET,
        output_prefix=OUTPUT_PREFIX,
        pseudocount=0.01,
        maximum_iterations=100,
        convergence_threshold=1e-7,
        random_seed=13,
        sort_threads=8,
        bam_progress_interval=1_000_000,
        parquet_batch_size=100_000,
    )

    print("Final TPM:", result.final_tpm_parquet)
    print("Assignments:", result.assignments_parquet)
    print("Final BAM:", result.final_bam)
