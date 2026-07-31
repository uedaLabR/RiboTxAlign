from em_distribute.transcript_em import (
    em_transcript,
)


GROUP_EM_TSV = (
    "/mnt/share/project/riboML/"
    "sample_ref/group_em_tpm.tsv"
)

SORTED_MAPPING_DIRECTORY = (
    "/mnt/share/project/riboML/"
    "sample_ref/sorted_mapping"
)


if __name__ == "__main__":


    em_transcript(
        group_em_tsv=GROUP_EM_TSV,
        sorted_mapping_dir=SORTED_MAPPING_DIRECTORY,
        reads_info_parquet="/mnt/share/project/riboML/sample_ref/transcript_reads_info.parquet",
        transcript_tpm_tsv="/mnt/share/project/riboML/sample_ref/transcript_tpm.tsv",
    )