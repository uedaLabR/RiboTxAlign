from em_distribute.group_effective_length import (
    build_group_effective_length_parquet,
)


GROUP_UNIQUE_PARQUET = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.group_unique.parquet"
)
GROUP_MULTI_PARQUET = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.group_multi.parquet"
)

COORDINATE_TSV = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.coordinates.tsv"
)

OUTPUT_PARQUET = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.group_length.parquet"
)

build_group_effective_length_parquet(
    group_unique_parquet=GROUP_UNIQUE_PARQUET,
    group_multi_parquet=GROUP_MULTI_PARQUET,
    coordinate_tsv=COORDINATE_TSV,
    output_parquet=OUTPUT_PARQUET,
    batch_size=100000,
    strict_coordinates=True,
)