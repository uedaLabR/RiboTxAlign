from em_distribute.transcript_em import (
    em_transcript,
)


GROUP_EM_TSV = (
    "/mnt/share/project/riboML/validation_expression/similation/HEK293T/RiboTxAlign/HEK293T_simulation/em/"
    "HEK293T_simulation.group_tpm.tsv"
)

SORTED_MAPPING_DIRECTORY = (
    "/mnt/share/project/riboML/validation_expression/similation/HEK293T/RiboTxAlign/HEK293T_simulation/mapping/sorted"
)


if __name__ == "__main__":


    em_transcript(
        group_em_tsv=GROUP_EM_TSV,
        sorted_mapping_dir=SORTED_MAPPING_DIRECTORY,
        reads_info_parquet="/mnt/share/project/riboML/validation_expression/similation/HEK293T/RiboTxAlign/HEK293T_simulation/em/transcript_reads_info.parquet",
        transcript_tpm_tsv="/mnt/share/project/riboML/validation_expression/similation/HEK293T/RiboTxAlign/HEK293T_simulation/em/transcript_tpm.tsv",
        transcript_merge_log_tsv="/mnt/share/project/riboML/validation_expression/similation/HEK293T/RiboTxAlign/HEK293T_simulation/em/transcript_merge_log_tsv",
        reference_metadata_tsv="/mnt/share/project/riboML/validation_expression/similation/HEK293T/RiboTxAlign/HEK293T_simulation/reference/HEK293T_simulation.metadata.tsv",
        reference_gtf = "/mnt/share/project/riboML/validation_expression/similation/HEK293T/RiboTxAlign/HEK293T_simulation/reference/HEK293T_simulation.gtf"
    )