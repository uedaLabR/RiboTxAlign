from parquet_work.sort_mapping_by_group import (
    sort_mapping_by_group,
)


MAPPING_PARTS_DIRECTORY = (
    "/mnt/share/project/riboML/sample_ref/"
    "mapping_parts"
)

SORTED_MAPPING_DIRECTORY = (
    "/mnt/share/project/riboML/sample_ref/"
    "sorted_mapping"
)


if __name__ == "__main__":

    from parquet_work.sort_mapping_by_group_mt import (
        sort_mapping_by_group,
    )

    output_files = sort_mapping_by_group(
        mapping_parts_directory=MAPPING_PARTS_DIRECTORY,
        output_directory=SORTED_MAPPING_DIRECTORY,
        read_batch_size=100000,
        row_group_size=100000,
        worker_count=4,
        overwrite=True,
        delete_intermediate_parquets=False,
    )

    # output_files = sort_mapping_by_group(
    #     mapping_parts_directory=(
    #         MAPPING_PARTS_DIRECTORY
    #     ),
    #     output_directory=(
    #         SORTED_MAPPING_DIRECTORY
    #     ),
    #     row_group_size=100000,
    #     overwrite=True,
    #     delete_intermediate_parquets=False,
    # )
    #
    # print(
    #     "Created sorted Parquet files:",
    #     len(
    #         output_files
    #     ),
    # )
    #
    # for output_file in output_files:
    #     print(
    #         output_file
    #     )