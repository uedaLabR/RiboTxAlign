from __future__ import annotations

from pipeline_profiling import measure, timed, profile_run, profile_child

import json
import logging
import multiprocessing as multiprocessing_module
import shutil
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import click
import yaml


@dataclass(frozen=True, slots=True)
class PipelinePaths:
    """Canonical paths shared by all RiboTxAlign commands."""

    sample_name: str
    work_dir: Path
    reference_prefix: Path
    input_fastq: Path
    deduplicated_fastq: Path
    dedup_bucket_dir: Path
    dedup_stats: Path
    mapping_parts_dir: Path
    mapping_log: Path
    mapping_stats: Path
    unmapped_fastq: Path
    sorted_mapping_dir: Path
    em_input_dir: Path
    group_tpm: Path
    group_em_log: Path
    transcript_em_log: Path
    reads_info: Path
    transcript_tpm: Path
    transcript_merge_log: Path
    depth_prefix: Path
    bam_dir: Path
    final_assignment: Path
    final_assignment_stats: Path
    assignment_parts_dir: Path
    main_unsorted_bam: Path
    rescue_unsorted_bam: Path
    rescue_stats: Path
    merged_unsorted_bam: Path
    final_bam: Path
    child_log_dir: Path

    @property
    def reference_fasta(self) -> Path:
        return Path(f"{self.reference_prefix}.fa")

    @property
    def reference_gtf(self) -> Path:
        return Path(f"{self.reference_prefix}.gtf")

    @property
    def coordinate_tsv(self) -> Path:
        return Path(f"{self.reference_prefix}.coordinates.tsv")

    @property
    def metadata_tsv(self) -> Path:
        return Path(f"{self.reference_prefix}.metadata.tsv")


@dataclass(frozen=True, slots=True)
class PipelineArguments:
    """Serializable arguments used to rebuild one pipeline context."""

    config_path: Path
    sample_name: str
    output_dir: Path
    fastq: Path


@dataclass(frozen=True, slots=True)
class PipelineContext:
    """Resolved configuration and file paths for one pipeline run."""

    arguments: PipelineArguments
    config: dict[str, Any]
    paths: PipelinePaths
    fastq_path: Path
    genome_fasta: Path
    genome_fai: Path
    gencode_gff3: Path


@dataclass(slots=True)
class ChildProcessHandle:
    """One child process and its persisted traceback file."""

    stage_name: str
    process: multiprocessing_module.Process
    error_path: Path


def load_yaml(config_path: str | Path) -> dict[str, Any]:
    """Load and validate one YAML configuration file."""

    path = Path(config_path).expanduser().resolve()
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}

    if not isinstance(config, dict):
        raise click.ClickException("The YAML root must be a mapping.")

    config["_config_dir"] = str(path.parent)
    return config


def section(config: dict[str, Any], name: str) -> dict[str, Any]:
    """Return one YAML section as a mapping."""

    value = config.get(name, {})
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise click.ClickException(f"YAML section '{name}' must be a mapping.")
    return value


def required(
    config: dict[str, Any],
    section_name: str,
    key: str,
) -> Any:
    """Return one required YAML value."""

    value = section(config, section_name).get(key)
    if value is None or value == "":
        raise click.ClickException(
            f"Missing required YAML value: {section_name}.{key}"
        )
    return value


def parameter(
    config: dict[str, Any],
    section_name: str,
    key: str,
    default: Any,
) -> Any:
    """Return one optional YAML value."""

    return section(config, section_name).get(key, default)


def resolve_config_path(
    config: dict[str, Any],
    value: str | Path,
) -> Path:
    """Resolve one path relative to the YAML directory."""

    path = Path(value).expanduser()
    if not path.is_absolute():
        path = Path(config["_config_dir"]) / path
    return path.resolve()


def resolve_cli_path(value: str | Path) -> Path:
    """Resolve one command-line path."""

    return Path(value).expanduser().resolve()


def build_paths(
    sample_name: str,
    output_dir: str | Path,
    input_fastq: str | Path | None = None,
) -> PipelinePaths:
    """Build all standard paths from command-line inputs."""

    sample_name = sample_name.strip()
    if not sample_name:
        raise click.ClickException("--sample-name must not be empty.")

    output_dir = resolve_cli_path(output_dir)
    work_dir = output_dir / sample_name
    reference_prefix = work_dir / "reference" / sample_name
    resolved_fastq = (
        resolve_cli_path(input_fastq)
        if input_fastq is not None
        else Path()
    )

    reads_dir = work_dir / "reads"
    mapping_dir = work_dir / "mapping"
    em_dir = work_dir / "em"
    depth_dir = work_dir / "depth"
    bam_dir = work_dir / "bam"
    logs_dir = work_dir / "logs"

    return PipelinePaths(
        sample_name=sample_name,
        work_dir=work_dir,
        reference_prefix=reference_prefix,
        input_fastq=resolved_fastq,
        deduplicated_fastq=reads_dir / f"{sample_name}.deduplicated.fastq.gz",
        dedup_bucket_dir=reads_dir / "umi_buckets",
        dedup_stats=reads_dir / f"{sample_name}.deduplication.stats.tsv",
        mapping_parts_dir=mapping_dir / "parts",
        mapping_log=mapping_dir / f"{sample_name}.mapping.log",
        mapping_stats=mapping_dir / f"{sample_name}.mapping_stats.json",
        unmapped_fastq=mapping_dir / f"{sample_name}.unmapped.fastq.gz",
        sorted_mapping_dir=mapping_dir / "sorted",
        em_input_dir=em_dir / "group_input",
        group_tpm=em_dir / f"{sample_name}.group_tpm.tsv",
        group_em_log=em_dir / f"{sample_name}.group_em.log",
        transcript_em_log=em_dir / f"{sample_name}.transcript_em.log",
        reads_info=em_dir / f"{sample_name}.reads_info.parquet",
        transcript_tpm=em_dir / f"{sample_name}.transcript_tpm.tsv",
        transcript_merge_log=em_dir / f"{sample_name}.transcript_merge.tsv",
        depth_prefix=depth_dir / sample_name,
        bam_dir=bam_dir,
        final_assignment=bam_dir / f"{sample_name}.final_assignments.parquet",
        final_assignment_stats=(
            bam_dir / f"{sample_name}.final_assignments.stats.tsv"
        ),
        assignment_parts_dir=bam_dir / "assignment_parts",
        main_unsorted_bam=(
            bam_dir / f"{sample_name}.transcript_supported.unsorted.bam"
        ),
        rescue_unsorted_bam=(
            bam_dir / f"{sample_name}.genome_rescued_unique.unsorted.bam"
        ),
        rescue_stats=mapping_dir / f"{sample_name}.genome_rescue_stats.json",
        merged_unsorted_bam=bam_dir / f"{sample_name}.merged.unsorted.bam",
        final_bam=bam_dir / f"{sample_name}.genome.sorted.bam",
        child_log_dir=logs_dir / "child_processes",
    )


def ensure_parent_directories(paths: PipelinePaths) -> None:
    """Create all standard output directories."""

    for path in (
        paths.reference_prefix.parent,
        paths.deduplicated_fastq.parent,
        paths.mapping_parts_dir,
        paths.sorted_mapping_dir,
        paths.em_input_dir,
        paths.group_tpm.parent,
        paths.depth_prefix.parent,
        paths.bam_dir,
        paths.assignment_parts_dir,
        paths.child_log_dir,
    ):
        path.mkdir(parents=True, exist_ok=True)


def stage(number: str, message: str) -> None:
    """Print one concise pipeline progress line."""

    click.echo(f"[{number}] {message}", err=True)


def create_logger(name: str) -> logging.Logger:
    """Create one shared console logger."""

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s")
        )
        logger.addHandler(handler)

    return logger


def remove_path(path: Path) -> None:
    """Remove one file or directory when it exists."""

    if not path.exists():
        return

    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def cleanup_deduplication_intermediates(
    config: dict[str, Any],
    paths: PipelinePaths,
) -> None:
    """Remove temporary UMI bucket files."""

    if not bool(
        parameter(config, "pipeline", "delete_intermediate_files", True)
    ):
        return

    remove_path(paths.dedup_bucket_dir)


def cleanup_run_intermediates(context: PipelineContext) -> None:
    """Remove intermediates after all output branches finish."""

    if not bool(
        parameter(
            context.config,
            "pipeline",
            "delete_intermediate_files",
            True,
        )
    ):
        return

    for path in (
        context.paths.mapping_parts_dir,
        context.paths.sorted_mapping_dir,
        context.paths.em_input_dir,
        context.paths.assignment_parts_dir,
        context.paths.final_assignment,
        context.paths.main_unsorted_bam,
        context.paths.rescue_unsorted_bam,
        context.paths.merged_unsorted_bam,
    ):
        remove_path(path)


def prepare_context(arguments: PipelineArguments) -> PipelineContext:
    """Resolve all configuration and input paths for one run."""

    config = load_yaml(arguments.config_path)
    paths = build_paths(
        sample_name=arguments.sample_name,
        output_dir=arguments.output_dir,
        input_fastq=arguments.fastq,
    )
    ensure_parent_directories(paths)

    genome_fasta = resolve_config_path(
        config,
        required(config, "reference", "genome_fasta"),
    )
    genome_fai_value = section(config, "reference").get("genome_fai")
    genome_fai = (
        resolve_config_path(config, genome_fai_value)
        if genome_fai_value
        else Path(f"{genome_fasta}.fai")
    )
    gencode_gff3 = resolve_config_path(
        config,
        required(config, "reference", "gencode_gff3"),
    )

    return PipelineContext(
        arguments=arguments,
        config=config,
        paths=paths,
        fastq_path=arguments.fastq,
        genome_fasta=genome_fasta,
        genome_fai=genome_fai,
        gencode_gff3=gencode_gff3,
    )


@timed("c")
def run_stage_c(context: PipelineContext) -> None:
    """Map Ribo-seq reads to the transcript reference."""

    from mapping.map_read_mt import run_transcript_mapping_to_parquet_mt

    config = context.config
    paths = context.paths

    stage("c", "Mapping reads to the transcript reference")
    run_transcript_mapping_to_parquet_mt(
        fastq_path=str(context.fastq_path),
        transcript_fasta=str(paths.reference_fasta),
        reference_metadata_tsv=str(paths.metadata_tsv),
        reference_gtf=str(paths.reference_gtf),
        mapping_parts_directory=str(paths.mapping_parts_dir),
        overwrite=bool(parameter(config, "mapping", "overwrite", True)),
        maximum_reads=parameter(config, "mapping", "maximum_reads", None),
        maximum_candidates=int(
            parameter(config, "mapping", "maximum_candidates", 16)
        ),
        kmer_size=int(parameter(config, "mapping", "kmer_size", 11)),
        window_size=int(parameter(config, "mapping", "window_size", 5)),
        minimum_minimizers=int(
            parameter(config, "mapping", "minimum_minimizers", 2)
        ),
        minimum_chain_score=int(
            parameter(config, "mapping", "minimum_chain_score", 15)
        ),
        minimum_dp_score=int(
            parameter(config, "mapping", "minimum_dp_score", 15)
        ),
        minimum_query_coverage=float(
            parameter(config, "mapping", "minimum_query_coverage", 0.85)
        ),
        minimum_identity=float(
            parameter(config, "mapping", "minimum_identity", 0.90)
        ),
        worker_count=int(
            parameter(config, "mapping", "worker_count", 6)
        ),
        reads_per_worker=int(
            parameter(config, "mapping", "reads_per_worker", 500000)
        ),
        progress_interval=int(
            parameter(config, "mapping", "progress_interval", 100000)
        ),
        parquet_row_group_size=int(
            parameter(
                config,
                "mapping",
                "parquet_row_group_size",
                100000,
            )
        ),
        exclude_decoy_reads=bool(
            parameter(config, "mapping", "exclude_decoy_reads", True)
        ),
        prefer_primary_on_alt_tie=bool(
            parameter(config, "mapping", "prefer_primary_on_alt_tie", False)
        ),
        decoy_reference_patterns=str(
            parameter(
                config,
                "mapping",
                "decoy_reference_patterns",
                "chrUn_*,*_random",
            )
        ),
        log_file=str(paths.mapping_log),
        unmapped_fastq_path=str(paths.unmapped_fastq),
        stats_json_path=str(paths.mapping_stats),
        multiprocessing_start_method=str(
            parameter(config, "mapping", "start_method", "spawn")
        ),
    )


@timed("d")
def run_stage_d(context: PipelineContext) -> None:
    """Sort mapping Parquet files by genome-equivalent group."""

    from parquet_work.sort_mapping_by_group_mt import sort_mapping_by_group

    config = context.config
    paths = context.paths

    stage("d", "Sorting mapping Parquet files by group")
    sort_mapping_by_group(
        mapping_parts_directory=str(paths.mapping_parts_dir),
        output_directory=str(paths.sorted_mapping_dir),
        read_batch_size=int(
            parameter(config, "sorting", "read_batch_size", 100000)
        ),
        row_group_size=int(
            parameter(config, "sorting", "row_group_size", 100000)
        ),
        worker_count=int(
            parameter(config, "sorting", "worker_count", 4)
        ),
        overwrite=bool(
            parameter(config, "sorting", "overwrite", True)
        ),
        delete_intermediate_parquets=bool(
            parameter(
                config,
                "sorting",
                "delete_intermediate_parquets",
                False,
            )
        ),
        multiprocessing_start_method=str(
            parameter(config, "sorting", "start_method", "spawn")
        ),
    )


@timed("e")
def run_stage_e(context: PipelineContext) -> None:
    """Prepare group-level EM input."""

    from em_distribute.extract_group_em_input import extract_em_input
    from em_distribute.group_cds_definition import create_group_cds_definitions

    config = context.config
    paths = context.paths

    stage("e", "Preparing group-level EM input")
    group_definitions = create_group_cds_definitions(
        gencode_gtf_path=context.gencode_gff3,
        reference_metadata_path=paths.metadata_tsv,
    )
    extract_em_input(
        sorted_mapping_directory=str(paths.sorted_mapping_dir),
        output_directory=str(paths.em_input_dir),
        group_definitions=group_definitions,
        batch_size=int(
            parameter(config, "group_input", "batch_size", 100000)
        ),
        write_batch_size=int(
            parameter(config, "group_input", "write_batch_size", 100000)
        ),
        worker_count=int(
            parameter(config, "group_input", "worker_count", 4)
        ),
        overwrite=bool(
            parameter(config, "group_input", "overwrite", True)
        ),
        minimum_mapped_reads=int(
            parameter(config, "group_input", "minimum_mapped_reads", 10)
        ),
        minimum_mapped_span=int(
            parameter(config, "group_input", "minimum_mapped_span", 100)
        ),
        minimum_transcript_fraction=float(
            parameter(
                config,
                "group_input",
                "minimum_transcript_fraction",
                1.0 / 3.0,
            )
        ),
    )


@timed("f")
def run_stage_f(context: PipelineContext) -> None:
    """Estimate group-level TPM by EM."""

    from em_distribute.group_em import run_group_em

    config = context.config
    paths = context.paths

    stage("f", "Estimating group-level TPM by EM")
    run_group_em(
        em_input_directory=paths.em_input_dir,
        mapping_log_path=paths.mapping_log,
        output_tsv_path=paths.group_tpm,
        output_log_path=paths.group_em_log,
        trial_count=int(
            parameter(config, "group_em", "trial_count", 300)
        ),
        convergence_threshold=float(
            parameter(
                config,
                "group_em",
                "convergence_threshold",
                60.0,
            )
        ),
        read_length=int(
            parameter(config, "group_em", "read_length", 32)
        ),
        batch_size=int(
            parameter(config, "group_em", "batch_size", 100000)
        ),
        engine=str(parameter(config, "group_em", "engine", "numpy")),
        keep_iteration_history=bool(
            parameter(config, "group_em", "keep_iteration_history", False)
        ),
    )


@timed("g")
def run_stage_g(context: PipelineContext) -> None:
    """Estimate transcript probabilities within groups."""

    from em_distribute.transcript_em import em_transcript

    paths = context.paths

    stage("g", "Estimating transcript probabilities within groups")
    em_transcript(
        group_em_tsv=paths.group_tpm,
        sorted_mapping_dir=paths.sorted_mapping_dir,
        reads_info_parquet=paths.reads_info,
        transcript_tpm_tsv=paths.transcript_tpm,
        transcript_merge_log_tsv=paths.transcript_merge_log,
        reference_metadata_tsv=paths.metadata_tsv,
        reference_gtf=paths.reference_gtf,
        output_log_path=paths.transcript_em_log,
    )


@timed("h")
def run_stage_h(context: PipelineContext) -> None:
    """Create transcript and genome depth tracks."""

    from depth_output.transcript_depth import create_transcript_depth_tracks

    config = context.config
    paths = context.paths

    stage("h", "Creating transcript and genome depth tracks")
    create_transcript_depth_tracks(
        reads_info_parquet=paths.reads_info,
        reference_metadata_tsv=paths.metadata_tsv,
        reference_gtf=paths.reference_gtf,
        transcript_fasta=paths.reference_fasta,
        coordinate_table=paths.coordinate_tsv,
        genome_fai=context.genome_fai,
        output_prefix=paths.depth_prefix,
        profile_minimum=int(
            parameter(config, "depth", "profile_minimum", -40)
        ),
        profile_maximum=int(
            parameter(config, "depth", "profile_maximum", 20)
        ),
        minimum_offset=int(
            parameter(config, "depth", "minimum_offset", 8)
        ),
        maximum_offset=int(
            parameter(config, "depth", "maximum_offset", 20)
        ),
        fallback_offset=int(
            parameter(config, "depth", "fallback_offset", 12)
        ),
        minimum_profile_weight=float(
            parameter(config, "depth", "minimum_profile_weight", 10.0)
        ),
        metagene_minimum=int(
            parameter(config, "depth", "metagene_minimum", -50)
        ),
        metagene_maximum=int(
            parameter(config, "depth", "metagene_maximum", 100)
        ),
        batch_size=int(
            parameter(config, "depth", "batch_size", 100000)
        ),
        create_bigwig=bool(
            parameter(config, "depth", "create_bigwig", True)
        ),
        create_genome_bigwig=bool(
            parameter(config, "depth", "create_genome_bigwig", True)
        ),
    )


@timed("i")
def run_stage_i(context: PipelineContext) -> None:
    """Assign one genomic alignment to each transcript-supported read."""

    from bam_output.final_alignment_assignment import (
        create_final_genome_assignments,
    )

    config = context.config
    paths = context.paths

    stage("i", "Assigning one genomic alignment to each mapped read")
    create_final_genome_assignments(
        reads_info_parquet=paths.reads_info,
        coordinate_table=paths.coordinate_tsv,
        output_parquet=paths.final_assignment,
        output_stats_tsv=paths.final_assignment_stats,
        batch_size=int(
            parameter(config, "assignment", "batch_size", 250000)
        ),
        local_bin_size=int(
            parameter(config, "assignment", "local_bin_size", 30)
        ),
        logger=create_logger("RiboTxAlignAssignment"),
    )


@timed("j")
def run_stage_j(context: PipelineContext) -> None:
    """Write the unsorted transcript-supported genome BAM."""

    from bam_output.bam_writer import write_bam_from_fastq_order_partitions
    from bam_output.fastq_order_partition import (
        partition_assignments_by_fastq_order,
    )

    config = context.config
    paths = context.paths

    stage("j", "Writing the unsorted transcript-supported genome BAM")
    partition_count = int(
        parameter(config, "bam", "partition_count", 10)
    )
    partition_prefix = f"{paths.sample_name}.final_assignments"

    with measure("j.partition_assignments"):
        partition_assignments_by_fastq_order(
            assignment_parquet=paths.final_assignment,
            fastq_path=context.fastq_path,
            output_dir=paths.assignment_parts_dir,
            output_prefix=partition_prefix,
            partition_count=partition_count,
            batch_size=int(
                parameter(config, "bam", "batch_size", 250000)
            ),
            compression=str(
                parameter(config, "bam", "compression", "zstd")
            ),
            logger=create_logger("RiboTxAlignBam"),
        )

    with measure("j.write_bam"):
        write_bam_from_fastq_order_partitions(
            fastq_path=context.fastq_path,
            partition_dir=paths.assignment_parts_dir,
            partition_prefix=partition_prefix,
            partition_count=partition_count,
            genome_fai=context.genome_fai,
            output_bam=paths.main_unsorted_bam,
            sort_output=False,
            index_output=False,
            logger=create_logger("RiboTxAlignBam"),
        )


def _write_json_atomic(path: Path, data: dict[str, Any]) -> None:
    """Write one JSON file atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = Path(f"{path}.tmp")
    temporary_path.write_text(
        json.dumps(data, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary_path.replace(path)


@timed("k")
def run_stage_k(context: PipelineContext) -> None:
    """Rescue unique genome alignments and persist rescue counters."""

    from mapping.genome_rescue import rescue_unique_genome_alignments

    config = context.config
    paths = context.paths

    stage("k", "Rescuing unique genome alignments from unmapped reads")
    temporary_bam = Path(f"{paths.rescue_unsorted_bam}.tmp")
    temporary_bam.unlink(missing_ok=True)
    paths.rescue_stats.unlink(missing_ok=True)

    rescue_counters = rescue_unique_genome_alignments(
        unmapped_fastq=paths.unmapped_fastq,
        genome_fasta=context.genome_fasta,
        genome_fai=context.genome_fai,
        output_bam=temporary_bam,
        remaining_unmapped_fastq=None,
        contiguous_preset=str(
            parameter(config, "rescue", "contiguous_preset", "sr")
        ),
        splice_preset=str(
            parameter(config, "rescue", "splice_preset", "splice")
        ),
        kmer_size=int(
            parameter(config, "rescue", "kmer_size", 11)
        ),
        window_size=int(
            parameter(config, "rescue", "window_size", 5)
        ),
        minimum_minimizers=int(
            parameter(config, "rescue", "minimum_minimizers", 2)
        ),
        minimum_chain_score=int(
            parameter(config, "rescue", "minimum_chain_score", 15)
        ),
        minimum_dp_score=int(
            parameter(config, "rescue", "minimum_dp_score", 15)
        ),
        minimum_query_coverage=float(
            parameter(config, "rescue", "minimum_query_coverage", 0.85)
        ),
        minimum_identity=float(
            parameter(config, "rescue", "minimum_identity", 0.90)
        ),
        enable_splice_rescue=bool(
            parameter(config, "rescue", "enable_splice_rescue", True)
        ),
        minimum_anchor_length=int(
            parameter(config, "rescue", "minimum_anchor_length", 11)
        ),
        minimum_intron_length=int(
            parameter(config, "rescue", "minimum_intron_length", 20)
        ),
        maximum_intron_length=int(
            parameter(config, "rescue", "maximum_intron_length", 2000)
        ),
        canonical_motifs=tuple(
            parameter(
                config,
                "rescue",
                "canonical_motifs",
                ["GT-AG", "GC-AG", "AT-AC"],
            )
        ),
        progress_interval=int(
            parameter(config, "rescue", "progress_interval", 100000)
        ),
        logger=create_logger("RiboTxAlignGenomeRescue"),
    )

    temporary_bam.replace(paths.rescue_unsorted_bam)
    _write_json_atomic(
        paths.rescue_stats,
        {
            "processed_reads": rescue_counters.processed_reads,
            "contiguous_rescued_unique_reads": (
                rescue_counters.contiguous_rescued_unique_reads
            ),
            "splice_rescued_unique_reads": (
                rescue_counters.splice_rescued_unique_reads
            ),
            "rescued_unique_reads": rescue_counters.rescued_unique_reads,
            "remaining_unmapped_reads": (
                rescue_counters.remaining_unmapped_reads
            ),
        },
    )


def load_rescue_counters(context: PipelineContext) -> dict[str, int]:
    """Load and validate persisted rescue counters."""

    if not context.paths.rescue_stats.exists():
        raise click.ClickException(
            f"Genome rescue statistics were not found: "
            f"{context.paths.rescue_stats}"
        )

    data = json.loads(
        context.paths.rescue_stats.read_text(encoding="utf-8")
    )
    required_keys = {
        "processed_reads",
        "contiguous_rescued_unique_reads",
        "splice_rescued_unique_reads",
        "rescued_unique_reads",
        "remaining_unmapped_reads",
    }
    missing = required_keys - set(data)
    if missing:
        raise click.ClickException(
            "Genome rescue statistics are incomplete: "
            + ", ".join(sorted(missing))
        )

    return {key: int(data[key]) for key in required_keys}


def update_rescue_mapping_statistics(
    context: PipelineContext,
    rescue_counters: dict[str, int],
) -> dict[str, Any]:
    """Add genome-rescue and final mapping values to the mapping JSON."""

    from mapping.mapping_stats import (
        load_mapping_stats,
        percentage,
        write_mapping_stats,
    )

    stats = load_mapping_stats(context.paths.mapping_stats)
    transcript_unmapped = int(stats["transcript_unmapped"])
    processed_reads = rescue_counters["processed_reads"]

    if processed_reads != transcript_unmapped:
        raise click.ClickException(
            "Genome rescue FASTQ count does not match transcript-unmapped "
            f"count: FASTQ={processed_reads}, "
            f"statistics={transcript_unmapped}"
        )

    input_reads = int(stats["input_reads"])
    transcript_mapped = int(stats["transcript_mapped"])
    group_unique = int(stats["group_unique"])
    group_multi = int(stats["group_multi"])
    decoy_filtered = int(stats["decoy_filtered"])
    analyzable_reads = input_reads - decoy_filtered

    contiguous_rescued = (
        rescue_counters["contiguous_rescued_unique_reads"]
    )
    splice_rescued = rescue_counters["splice_rescued_unique_reads"]
    rescued_unique = rescue_counters["rescued_unique_reads"]
    final_unmapped = rescue_counters["remaining_unmapped_reads"]
    final_mapped = transcript_mapped + rescued_unique
    final_unique = group_unique + rescued_unique
    final_multi = group_multi

    stats.update({
        "genome_rescue_input": processed_reads,
        "genome_contiguous_rescued_unique": contiguous_rescued,
        "genome_splice_rescued_unique": splice_rescued,
        "genome_rescued_unique": rescued_unique,
        "genome_rescue_rate_percent": percentage(
            rescued_unique,
            processed_reads,
        ),
        "final_mapped": final_mapped,
        "final_unique": final_unique,
        "final_multi": final_multi,
        "final_unmapped": final_unmapped,
        "final_mapped_rate_percent": percentage(
            final_mapped,
            input_reads,
        ),
        "final_unique_rate_percent": percentage(
            final_unique,
            input_reads,
        ),
        "final_multi_rate_percent": percentage(
            final_multi,
            input_reads,
        ),
        "final_unmapped_rate_percent": percentage(
            final_unmapped,
            input_reads,
        ),
        "final_mapped_non_decoy_rate_percent": percentage(
            final_mapped,
            analyzable_reads,
        ),
    })
    write_mapping_stats(context.paths.mapping_stats, stats)
    return stats

def _count_primary_bam_records(
    path: Path,
    pysam_module,
) -> int:
    """Count primary mapped BAM records."""

    count = 0

    with pysam_module.AlignmentFile(
        path,
        "rb",
    ) as bam_handle:

        for record in bam_handle.fetch(
            until_eof=True
        ):
            if (
                record.is_secondary
                or record.is_supplementary
                or record.is_unmapped
            ):
                continue

            count += 1

    return count


@timed("l")
def run_stage_l(
    context: PipelineContext,
    stats: dict[str, Any],
) -> None:
    """Merge, sort, index, and validate the final genome BAM."""

    from mapping.mapping_stats import write_mapping_stats

    try:
        import pysam
    except ImportError as error:
        raise click.ClickException(
            "Final BAM processing requires pysam."
        ) from error

    config = context.config
    paths = context.paths

    stage("l", "Merging, sorting, indexing, and validating BAM")

    with measure("l.count.main_unsorted_bam"):
        main_bam_reads = _count_primary_bam_records(
            paths.main_unsorted_bam,
            pysam,
        )
    with measure("l.count.rescue_unsorted_bam"):
        rescue_bam_reads = _count_primary_bam_records(
            paths.rescue_unsorted_bam,
            pysam,
        )

    paths.merged_unsorted_bam.unlink(missing_ok=True)
    with measure("l.merge"):
        pysam.merge(
            "-f",
            str(paths.merged_unsorted_bam),
            str(paths.main_unsorted_bam),
            str(paths.rescue_unsorted_bam),
        )

    with measure("l.count.merged_unsorted_bam"):
        merged_bam_reads = _count_primary_bam_records(
            paths.merged_unsorted_bam,
            pysam,
        )
    expected_final_reads = main_bam_reads + rescue_bam_reads

    if merged_bam_reads != expected_final_reads:
        raise click.ClickException(
            "Merged BAM count is inconsistent: "
            f"main={main_bam_reads}, rescue={rescue_bam_reads}, "
            f"merged={merged_bam_reads}"
        )

    sort_threads = max(
        int(parameter(config, "bam", "sort_threads", 8)),
        1,
    )
    sort_memory = str(
        parameter(config, "bam", "sort_memory", "2G")
    )

    with measure("l.sort"):
        pysam.sort(
            "-@",
            str(sort_threads),
            "-m",
            sort_memory,
            "-o",
            str(paths.final_bam),
            str(paths.merged_unsorted_bam),
        )
    with measure("l.index"):
        pysam.index(
            "-@",
            str(sort_threads),
            str(paths.final_bam),
        )

    with measure("l.count.final_bam"):
        final_bam_reads = _count_primary_bam_records(
            paths.final_bam,
            pysam,
        )
    import pyarrow.parquet as pq
    from mapping.mapping_stats import validate_assignment_bam_counts

    assigned_reads = pq.ParquetFile(paths.final_assignment).metadata.num_rows
    try:
        with measure("l.validate_counts"):
            assignment_stats = validate_assignment_bam_counts(
                stats, assigned_reads, main_bam_reads, rescue_bam_reads, final_bam_reads
            )
    except ValueError as error:
        raise click.ClickException(str(error)) from error
    stats.update(assignment_stats)

    stats.update({
        "main_bam_primary_reads": main_bam_reads,
        "rescue_bam_primary_reads": rescue_bam_reads,
        "final_bam_primary_reads": final_bam_reads,
        "bam_merge_validated": True,
        "final_bam": str(paths.final_bam),
        "final_bam_index": str(paths.final_bam) + ".bai",
    })
    write_mapping_stats(paths.mapping_stats, stats)
    paths.merged_unsorted_bam.unlink(missing_ok=True)


def _write_child_error(
    error_path: Path,
    stage_name: str,
    error: BaseException,
) -> None:
    """Persist one child-process traceback."""

    error_path.parent.mkdir(parents=True, exist_ok=True)
    error_path.write_text(
        (
            f"Stage: {stage_name}\n"
            f"Error: {type(error).__name__}: {error}\n\n"
            f"{traceback.format_exc()}"
        ),
        encoding="utf-8",
    )


@profile_child
def _child_stage_entry(
    stage_name: str,
    arguments: PipelineArguments,
    error_path: Path,
) -> None:
    """Rebuild context and run one named child stage."""

    try:
        error_path.unlink(missing_ok=True)
        context = prepare_context(arguments)

        if stage_name == "h":
            run_stage_h(context)
        elif stage_name == "k":
            run_stage_k(context)
        else:
            raise ValueError(f"Unsupported child stage: {stage_name}")

    except BaseException as error:
        _write_child_error(
            error_path=error_path,
            stage_name=stage_name,
            error=error,
        )
        raise


def _start_child(
    multiprocessing_context,
    stage_name: str,
    arguments: PipelineArguments,
    error_path: Path,
    logger: logging.Logger,
) -> ChildProcessHandle:
    """Start one child stage."""

    error_path.unlink(missing_ok=True)
    process = multiprocessing_context.Process(
        target=_child_stage_entry,
        args=(stage_name, arguments, error_path),
        name=f"RiboTxAlign-stage-{stage_name}",
    )
    process.start()

    logger.info(
        "Started child stage %s: pid=%s",
        stage_name,
        process.pid,
    )
    return ChildProcessHandle(
        stage_name=stage_name,
        process=process,
        error_path=error_path,
    )


def _child_failure_message(handle: ChildProcessHandle) -> str:
    """Return the persisted child traceback when available."""

    if handle.error_path.exists():
        return handle.error_path.read_text(
            encoding="utf-8",
            errors="replace",
        )

    return (
        f"Stage {handle.stage_name} failed with exit code "
        f"{handle.process.exitcode}, but no traceback was produced."
    )


def _raise_if_child_failed(handle: ChildProcessHandle) -> None:
    """Raise when one child has already failed."""

    exit_code = handle.process.exitcode
    if exit_code is None or exit_code == 0:
        return

    raise click.ClickException(_child_failure_message(handle))


def _wait_for_child(
    handle: ChildProcessHandle,
    logger: logging.Logger,
) -> None:
    """Wait for one child and validate completion."""

    if handle.process.is_alive():
        logger.info(
            "Waiting for child stage %s: pid=%s",
            handle.stage_name,
            handle.process.pid,
        )

    handle.process.join()
    _raise_if_child_failed(handle)
    logger.info(
        "Child stage %s completed successfully.",
        handle.stage_name,
    )


def _terminate_child(
    handle: ChildProcessHandle | None,
    logger: logging.Logger,
    timeout_seconds: float = 30.0,
) -> None:
    """Terminate one child after another branch fails."""

    if handle is None or not handle.process.is_alive():
        return

    logger.error(
        "Terminating child stage %s: pid=%s",
        handle.stage_name,
        handle.process.pid,
    )
    handle.process.terminate()
    handle.process.join(timeout=timeout_seconds)

    if handle.process.is_alive():
        logger.error(
            "Killing child stage %s: pid=%s",
            handle.stage_name,
            handle.process.pid,
        )
        handle.process.kill()
        handle.process.join()


def report_pipeline_outputs(context: PipelineContext) -> None:
    """Print final output paths."""

    click.echo("RiboTxAlign workflow completed.")
    click.echo(f"Final BAM: {context.paths.final_bam}")
    click.echo(f"BAM index: {context.paths.final_bam}.bai")
    click.echo(
        f"Transcript-unmapped FASTQ: {context.paths.unmapped_fastq}"
    )
    click.echo(f"Mapping statistics: {context.paths.mapping_stats}")


@profile_run
def run_pipeline(
    config_path: Path,
    sample_name: str,
    output_dir: Path,
    fastq: Path,
) -> None:
    """Run stages c-l with independent rescue and depth branches."""

    logger = create_logger("RiboTxAlignPipeline")
    arguments = PipelineArguments(
        config_path=resolve_cli_path(config_path),
        sample_name=sample_name,
        output_dir=resolve_cli_path(output_dir),
        fastq=resolve_cli_path(fastq),
    )
    context = prepare_context(arguments)

    start_method = str(
        parameter(
            context.config,
            "pipeline",
            "start_method",
            "spawn",
        )
    )
    multiprocessing_context = multiprocessing_module.get_context(
        start_method
    )

    rescue_handle: ChildProcessHandle | None = None
    depth_handle: ChildProcessHandle | None = None

    try:
        run_stage_c(context)

        rescue_handle = _start_child(
            multiprocessing_context=multiprocessing_context,
            stage_name="k",
            arguments=arguments,
            error_path=(
                context.paths.child_log_dir
                / f"{sample_name}.stage_k.error.txt"
            ),
            logger=logger,
        )

        run_stage_d(context)
        _raise_if_child_failed(rescue_handle)

        run_stage_e(context)
        _raise_if_child_failed(rescue_handle)

        run_stage_f(context)
        _raise_if_child_failed(rescue_handle)

        run_stage_g(context)
        _raise_if_child_failed(rescue_handle)

        depth_handle = _start_child(
            multiprocessing_context=multiprocessing_context,
            stage_name="h",
            arguments=arguments,
            error_path=(
                context.paths.child_log_dir
                / f"{sample_name}.stage_h.error.txt"
            ),
            logger=logger,
        )

        run_stage_i(context)
        _raise_if_child_failed(rescue_handle)
        _raise_if_child_failed(depth_handle)

        run_stage_j(context)
        _raise_if_child_failed(rescue_handle)
        _raise_if_child_failed(depth_handle)

        with measure("wait.k"):
            _wait_for_child(rescue_handle, logger)
        rescue_counters = load_rescue_counters(context)
        stats = update_rescue_mapping_statistics(
            context,
            rescue_counters,
        )

        run_stage_l(context, stats)

        with measure("wait.h"):
            _wait_for_child(depth_handle, logger)

        stage("cleanup", "Removing intermediate files")
        with measure("cleanup"):
            cleanup_run_intermediates(context)
        report_pipeline_outputs(context)

    except BaseException:
        _terminate_child(depth_handle, logger)
        _terminate_child(rescue_handle, logger)
        raise
