"""Run mapping-free UMI deduplication from PyCharm."""

from pathlib import Path

from de_duplication.umi_deduplication import (
    process_buckets,
    split_fastq_into_buckets,
    write_stats,
)


RIBO_FASTQ = (
    "/mnt/share/project/riboML/validation/HEK293T/riboseq/processed/clean/"
    "SRR26942051.trimmed.clean.fastq.gz"
)



OUTPUT_DIR = Path(
    "/mnt/share/project/riboML/riboseq_raw/trimmed_26_34/deduplicated"
)

BUCKET_DIR = OUTPUT_DIR / "SRR26942051_umi_buckets"

OUTPUT_FASTQ = (
    OUTPUT_DIR
    / "SRR26942051.26_34.deduplicated.fastq.gz"
)

STATS_PATH = (
    OUTPUT_DIR
    / "SRR26942051.26_34.deduplicated.stats.tsv"
)


def run_test_deduplication() -> None:
    """Create UMI buckets and collapse duplicate reads sequentially."""

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    bucket_stats = split_fastq_into_buckets(
        input_fastq=RIBO_FASTQ,
        output_dir=BUCKET_DIR,
    )

    collapse_stats = process_buckets(
        bucket_dir=BUCKET_DIR,
        output_fastq_gz=OUTPUT_FASTQ,
    )

    write_stats(
        path=STATS_PATH,
        bucket_stats=bucket_stats,
        collapse_stats=collapse_stats,
    )

    print(
        f"Input reads: "
        f"{bucket_stats.total_reads:,}"
    )
    print(
        f"Written bucket records: "
        f"{bucket_stats.written_records:,}"
    )
    print(
        f"Exact sequence groups: "
        f"{collapse_stats.exact_sequence_groups:,}"
    )
    print(
        f"Collapsed groups: "
        f"{collapse_stats.collapsed_sequence_groups:,}"
    )
    print(
        f"Output reads: "
        f"{collapse_stats.output_reads:,}"
    )
    print(
        f"Output FASTQ: {OUTPUT_FASTQ}"
    )
    print(
        f"Statistics: {STATS_PATH}"
    )


run_test_deduplication()

