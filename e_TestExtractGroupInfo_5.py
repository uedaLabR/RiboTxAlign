from em_distribute.group_cds_definition import (
    create_group_cds_definitions,
)
from em_distribute.extract_group_em_input import (
    extract_em_input,
)

def main() -> None:

    """Create distributed EM input Parquet files."""
    print("start reading def")
    group_definitions = (
        create_group_cds_definitions(
            gencode_gtf_path=(
                "/mnt/ssdnas07/pipeline/"
                "rna_modscore_v2/source/"
                "gencode.v44.annotation.gff3"
            ),
            reference_metadata_path=(
                "/mnt/share/project/riboML/"
                "sample_ref/"
                "HEK293T_RiboTxAlign.metadata.tsv"
            ),
        )
    )
    print("finish reading def")

    extract_em_input(
        sorted_mapping_directory=(
            "/mnt/share/project/riboML/"
            "sample_ref/sorted_mapping"
        ),
        output_directory=(
            "/mnt/share/project/riboML/"
            "sample_ref/em_input"
        ),
        group_definitions=group_definitions,
        batch_size=100000,
        write_batch_size=100000,
        worker_count=4,
        overwrite=True,
    )

if __name__ == "__main__":
    main()