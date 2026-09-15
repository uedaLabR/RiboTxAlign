from __future__ import annotations

from pathlib import Path

import click

from RiboTxAlign_pipeline import (
    build_paths,
    cleanup_deduplication_intermediates,
    ensure_parent_directories,
    load_yaml,
    required,
    resolve_cli_path,
    resolve_config_path,
    run_pipeline,
    stage,
)


CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
    "show_default": True,
}


def _sample_option(function):
    """Add the shared sample-name option."""

    return click.option(
        "--sample-name",
        required=True,
        type=str,
        help="Sample name used consistently in output paths and files.",
    )(function)


def _output_option(function):
    """Add the shared output-directory option."""

    return click.option(
        "--output-dir",
        required=True,
        type=click.Path(file_okay=False, path_type=Path),
        help="Root output directory.",
    )(function)


def _fastq_option(function):
    """Add the shared input-FASTQ option."""

    return click.option(
        "--fastq",
        required=True,
        type=click.Path(
            exists=True,
            dir_okay=False,
            path_type=Path,
        ),
        help="Input single-end FASTQ or FASTQ.GZ file.",
    )(function)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version="0.5.0", prog_name="RiboTxAlign")
def cli() -> None:
    """Long-read-guided mapping and EM assignment of Ribo-seq reads."""


@cli.command("create-reference")
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        path_type=Path,
    ),
    help="Pipeline YAML file containing reference settings.",
)
@_sample_option
@_output_option
@click.option(
    "--stringtie-gtf",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        path_type=Path,
    ),
    help="Sample-specific StringTie2 transcript GTF.",
)
def create_reference_command(
    config_path: Path,
    sample_name: str,
    output_dir: Path,
    stringtie_gtf: Path,
) -> None:
    """Run step a: create the transcript reference."""

    from db_construct.build_reference import build_transcript_reference

    config = load_yaml(config_path)
    paths = build_paths(
        sample_name=sample_name,
        output_dir=output_dir,
    )
    ensure_parent_directories(paths)

    gencode_gff3 = resolve_config_path(
        config,
        required(config, "reference", "gencode_gff3"),
    )
    genome_fasta = resolve_config_path(
        config,
        required(config, "reference", "genome_fasta"),
    )
    stringtie_gtf = resolve_cli_path(stringtie_gtf)

    stage("a", "Creating the transcript reference")
    reference_transcripts = build_transcript_reference(
        gencode_gff3=str(gencode_gff3),
        stringtie_gtf=str(stringtie_gtf),
        genome_fasta=str(genome_fasta),
        output_prefix=str(paths.reference_prefix),
    )

    click.echo(f"Reference transcripts: {len(reference_transcripts):,}")
    click.echo(f"Reference prefix: {paths.reference_prefix}")




@cli.command("deduplicate")
@_fastq_option
@click.option(
    "--output-fastq",
    required=True,
    type=click.Path(
        dir_okay=False,
        path_type=Path,
    ),
    help="Output deduplicated FASTQ.GZ file.",
)
def deduplicate_command(
    fastq: Path,
    output_fastq: Path,
) -> None:
    """Perform mapping-free UMI deduplication."""

    import shutil

    from de_duplication.umi_deduplication import (
        process_buckets,
        split_fastq_into_buckets,
        write_stats,
    )

    fastq = fastq.resolve()
    output_fastq = output_fastq.resolve()

    output_fastq.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    bucket_dir = Path(
        f"{output_fastq}.umi_buckets"
    )

    stats_path = Path(
        f"{output_fastq}.stats.tsv"
    )

    stage(
        "b",
        "Deduplicating the input FASTQ",
    )

    bucket_stats = split_fastq_into_buckets(
        input_fastq=fastq,
        output_dir=bucket_dir,
    )

    collapse_stats = process_buckets(
        bucket_dir=bucket_dir,
        output_fastq_gz=output_fastq,
    )

    write_stats(
        path=stats_path,
        bucket_stats=bucket_stats,
        collapse_stats=collapse_stats,
    )

    shutil.rmtree(
        bucket_dir,
        ignore_errors=True,
    )

    click.echo(
        f"Input reads: {bucket_stats.total_reads:,}"
    )
    click.echo(
        f"Output reads: {collapse_stats.output_reads:,}"
    )
    click.echo(
        f"Deduplicated FASTQ: {output_fastq}"
    )
    click.echo(
        f"Statistics: {stats_path}"
    )

@cli.command("run")
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        path_type=Path,
    ),
    help="Pipeline YAML file containing reference and analysis settings.",
)
@_sample_option
@_output_option
@_fastq_option
def run_command(
    config_path: Path,
    sample_name: str,
    output_dir: Path,
    fastq: Path,
) -> None:
    """Run stages c through l with parallel rescue and depth branches."""

    run_pipeline(
        config_path=config_path,
        sample_name=sample_name,
        output_dir=output_dir,
        fastq=fastq,
    )


if __name__ == "__main__":
    cli()