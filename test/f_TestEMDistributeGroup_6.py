from em_distribute.group_em import (
    run_group_em,
)

result = run_group_em(
    em_input_directory=(
        "/mnt/share/project/riboML/"
        "sample_ref/em_input"
    ),
    mapping_log_path=(
        "/mnt/share/project/riboML/"
        "sample_ref/HEK293T_RiboTxAlign.mapping_mt.log"
    ),
    output_tsv_path=(
        "/mnt/share/project/riboML/"
        "sample_ref/group_em_tpm.tsv"
    ),
    output_log_path=(
        "/mnt/share/project/riboML/"
        "sample_ref/group_em.log"
    ),
    trial_count=300,
    convergence_threshold=60.0,
    read_length=32,
    batch_size=100000,
)