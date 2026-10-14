from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass, field
from pathlib import Path

import pyarrow.parquet as pq

from em_distribute.group_em_numpy import NumpyMultiHits


@dataclass(slots=True)
class MultiHitRecord:
    """One multi-hit read associated with one candidate group."""

    read_id: str
    other_group_ids: tuple[str, ...]


@dataclass(slots=True)
class GroupEmInfo:
    """All EM information retained for one mapping group."""

    group_id: str
    cds_length: int
    unique_read_count: int

    multi_hit_records: list[MultiHitRecord] = field(
        default_factory=list
    )

    tpm_by_iteration: list[float] = field(
        default_factory=list
    )

    expected_count_by_iteration: list[float] = field(
        default_factory=list
    )

    # The NumPy backend stores candidate rows compactly, outside this object.
    multi_hit_record_count: int | None = None

    def effective_length(
        self,
        read_length: int,
    ) -> int:
        """Return an effective CDS length for fixed-length reads."""

        return max(
            self.cds_length
            - read_length
            + 1,
            1,
        )


@dataclass(slots=True)
class EmRunResult:
    """Final result returned by the EM workflow."""

    groups: dict[str, GroupEmInfo]
    final_tpm_by_group: dict[str, float]
    mapped_read_count: int
    completed_iterations: int
    converged: bool
    final_delta: float


def create_logger(
    log_path: str | Path,
) -> logging.Logger:
    """Create a file and console logger."""

    path = Path(
        log_path
    ).resolve()

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    logger = logging.getLogger(
        f"RiboTxAlignEM:{path}"
    )
    logger.setLevel(
        logging.INFO
    )
    logger.propagate = False
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s\t%(levelname)s\t%(message)s"
    )

    file_handler = logging.FileHandler(
        path,
        mode="w",
        encoding="utf-8",
    )
    file_handler.setFormatter(
        formatter
    )

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(
        formatter
    )

    logger.addHandler(
        file_handler
    )
    logger.addHandler(
        stream_handler
    )

    return logger


def read_mapped_read_count(
    mapping_log_path: str | Path,
) -> int:
    """
    Read the final mapped-read count from a mapping log.

    The parser first looks for a final line such as:

        Mapped: 123456

    If unavailable, it uses the last progress value such as:

        Mapped=123456
    """

    path = Path(
        mapping_log_path
    ).resolve()

    if not path.is_file():
        raise FileNotFoundError(
            "Mapping log was not found: "
            f"{path}"
        )

    mapped_values: list[int] = []

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            final_match = re.search(
                r"\bMapped:\s*([0-9]+)\b",
                line,
            )

            if final_match is not None:
                mapped_values.append(
                    int(
                        final_match.group(
                            1
                        )
                    )
                )
                continue

            progress_match = re.search(
                r"\bMapped=([0-9]+)\b",
                line,
            )

            if progress_match is not None:
                mapped_values.append(
                    int(
                        progress_match.group(
                            1
                        )
                    )
                )

    if not mapped_values:
        raise ValueError(
            "No mapped-read count was found in "
            f"{path}"
        )

    return mapped_values[-1]


def find_em_input_files(
    em_input_directory: str | Path,
) -> tuple[
    list[Path],
    list[Path],
]:
    """Find group-info and multi-read Parquet shards."""

    input_directory = Path(
        em_input_directory
    ).resolve()

    if not input_directory.is_dir():
        raise NotADirectoryError(
            "EM-input directory was not found: "
            f"{input_directory}"
        )

    group_info_paths = sorted(
        input_directory.glob(
            "*.group_info.parquet"
        )
    )

    multi_read_paths = sorted(
        input_directory.glob(
            "*.multi_reads.parquet"
        )
    )

    if not group_info_paths:
        raise FileNotFoundError(
            "No group-info Parquet files were found in "
            f"{input_directory}"
        )

    if not multi_read_paths:
        raise FileNotFoundError(
            "No multi-read Parquet files were found in "
            f"{input_directory}"
        )

    return (
        group_info_paths,
        multi_read_paths,
    )


def load_group_information(
    group_info_paths: list[Path],
    batch_size: int,
) -> dict[str, GroupEmInfo]:
    """Load group ID, CDS length, and unique-read count."""

    groups: dict[
        str,
        GroupEmInfo,
    ] = {}

    columns = [
        "group_id",
        "cds_length",
        "group_unique_read_count",
    ]

    for parquet_path in group_info_paths:
        parquet_file = pq.ParquetFile(
            parquet_path
        )

        for batch in parquet_file.iter_batches(
            batch_size=batch_size,
            columns=columns,
        ):
            data = batch.to_pydict()

            for (
                group_id,
                cds_length,
                unique_read_count,
            ) in zip(
                data["group_id"],
                data["cds_length"],
                data[
                    "group_unique_read_count"
                ],
            ):
                if group_id in groups:
                    raise ValueError(
                        "Duplicate group_id was found: "
                        f"{group_id}"
                    )

                groups[group_id] = GroupEmInfo(
                    group_id=group_id,
                    cds_length=int(
                        cds_length
                    ),
                    unique_read_count=int(
                        unique_read_count
                    ),
                )

    return groups


def load_multi_hit_records(
    groups: dict[str, GroupEmInfo],
    multi_read_paths: list[Path],
    batch_size: int,
) -> int:
    """
    Attach multi-hit records to their candidate groups.

    One read is retained once per candidate group. Multiple transcript
    hits within the same group have already been collapsed upstream.
    """

    columns = [
        "group_id",
        "read_id",
        "other_group_ids",
    ]

    record_count = 0

    for parquet_path in multi_read_paths:
        parquet_file = pq.ParquetFile(
            parquet_path
        )

        for batch in parquet_file.iter_batches(
            batch_size=batch_size,
            columns=columns,
        ):
            data = batch.to_pydict()

            for (
                group_id,
                read_id,
                other_group_ids,
            ) in zip(
                data["group_id"],
                data["read_id"],
                data["other_group_ids"],
            ):
                group = groups.get(
                    group_id
                )

                if group is None:
                    raise KeyError(
                        "A multi-hit record refers to an "
                        "unknown group_id: "
                        f"{group_id}"
                    )

                group.multi_hit_records.append(
                    MultiHitRecord(
                        read_id=read_id,
                        other_group_ids=tuple(
                            other_group_ids
                            or ()
                        ),
                    )
                )

                record_count += 1

    return record_count


def calculate_tpm(
    count_by_group: dict[str, float],
    groups: dict[str, GroupEmInfo],
    read_length: int,
) -> dict[str, float]:
    """Calculate group TPM values from expected read counts."""

    rate_by_group: dict[
        str,
        float,
    ] = {}

    for group_id, group in groups.items():
        count = count_by_group.get(
            group_id,
            0.0,
        )

        rate_by_group[group_id] = (
            count
            / group.effective_length(
                read_length
            )
        )

    total_rate = sum(
        rate_by_group.values()
    )

    if total_rate <= 0.0:
        return {
            group_id: 0.0
            for group_id in groups
        }

    scale = (
        1_000_000.0
        / total_rate
    )

    return {
        group_id: rate * scale
        for group_id, rate
        in rate_by_group.items()
    }


def initialize_iteration_zero(
    groups: dict[str, GroupEmInfo],
    read_length: int,
) -> dict[str, float]:
    """
    Calculate iteration-0 TPM from unique reads only.

    The returned dictionary is the current TPM lookup table used by the
    first EM iteration.
    """

    unique_count_by_group = {
        group_id: float(
            group.unique_read_count
        )
        for group_id, group
        in groups.items()
    }

    tpm_by_group = calculate_tpm(
        count_by_group=(
            unique_count_by_group
        ),
        groups=groups,
        read_length=read_length,
    )

    for group_id, group in groups.items():
        group.expected_count_by_iteration.append(
            unique_count_by_group[
                group_id
            ]
        )

        group.tpm_by_iteration.append(
            tpm_by_group[
                group_id
            ]
        )

    return tpm_by_group


def calculate_multi_hit_contribution(
    group_id: str,
    record: MultiHitRecord,
    current_tpm_by_group: dict[str, float],
) -> float:
    """
    Calculate the posterior contribution of one read to one group.

    The same read has one record for every candidate group. The
    contributions across all candidate groups sum to one.
    """

    candidate_group_ids = (
        group_id,
        *record.other_group_ids,
    )

    candidate_tpm_values = [
        current_tpm_by_group.get(
            candidate_group_id,
            0.0,
        )
        for candidate_group_id
        in candidate_group_ids
    ]

    denominator = sum(
        candidate_tpm_values
    )

    if denominator > 0.0:
        return (
            current_tpm_by_group.get(
                group_id,
                0.0,
            )
            / denominator
        )

    candidate_count = len(
        candidate_group_ids
    )

    if candidate_count == 0:
        return 0.0

    return (
        1.0
        / candidate_count
    )


def run_em_iteration(
    groups: dict[str, GroupEmInfo],
    current_tpm_by_group: dict[str, float],
    read_length: int,
) -> tuple[
    dict[str, float],
    dict[str, float],
    float,
]:
    """Run one EM iteration and return new counts, TPM, and delta."""

    expected_count_by_group: dict[
        str,
        float,
    ] = {}

    for group_id, group in groups.items():
        expected_count = float(
            group.unique_read_count
        )

        # Avoid cache overhead when the first records have no repeated pattern.
        # This only selects a computation path; every read is still included.
        sample = group.multi_hit_records[:64]
        if len({record.other_group_ids for record in sample}) == len(sample):
            for record in group.multi_hit_records:
                expected_count += calculate_multi_hit_contribution(
                    group_id=group_id,
                    record=record,
                    current_tpm_by_group=current_tpm_by_group,
                )
        else:
            # Candidate order is part of the key: preserve floating-point sums
            # within each posterior AND across reads. TPM changes at every
            # M-step, so cache only within this group/iteration.
            contributions = {}
            for record in group.multi_hit_records:
                candidates = record.other_group_ids
                contribution = contributions.get(candidates)
                if contribution is None:
                    contribution = calculate_multi_hit_contribution(
                        group_id=group_id,
                        record=record,
                        current_tpm_by_group=current_tpm_by_group,
                    )
                    # Bound memory for unusually diverse candidate sets.
                    if len(contributions) >= 4096:
                        contributions.clear()
                    contributions[candidates] = contribution
                expected_count += contribution

        expected_count_by_group[
            group_id
        ] = expected_count

    new_tpm_by_group = calculate_tpm(
        count_by_group=(
            expected_count_by_group
        ),
        groups=groups,
        read_length=read_length,
    )

    delta = max(
        (
            abs(
                new_tpm_by_group[
                    group_id
                ]
                - current_tpm_by_group.get(
                    group_id,
                    0.0,
                )
            )
            for group_id in groups
        ),
        default=0.0,
    )

    return (
        expected_count_by_group,
        new_tpm_by_group,
        delta,
    )


def append_iteration_results(
    groups: dict[str, GroupEmInfo],
    expected_count_by_group: dict[
        str,
        float,
    ],
    tpm_by_group: dict[str, float],
    keep_iteration_history: bool = True,
) -> None:
    """Store all iterations, or only iteration zero and the latest result."""

    for group_id, group in groups.items():
        if not keep_iteration_history:
            del group.expected_count_by_iteration[1:]
            del group.tpm_by_iteration[1:]
        group.expected_count_by_iteration.append(
            expected_count_by_group[
                group_id
            ]
        )

        group.tpm_by_iteration.append(
            tpm_by_group[
                group_id
            ]
        )


def _load_em_data(
    em_input_directory: str | Path,
    mapping_log_path: str | Path,
    read_length: int,
    batch_size: int,
    logger: logging.Logger,
    engine: str,
) -> tuple[
    dict[str, GroupEmInfo],
    dict[str, float],
    int,
    int,
    NumpyMultiHits | None,
]:
    """
    Load EM input and create the iteration-0 TPM dictionary.

    This function performs all initialization required before the EM loop.
    """

    mapped_read_count = (
        read_mapped_read_count(
            mapping_log_path
        )
    )

    (
        group_info_paths,
        multi_read_paths,
    ) = find_em_input_files(
        em_input_directory
    )

    groups = load_group_information(
        group_info_paths=group_info_paths,
        batch_size=batch_size,
    )

    numpy_data = None
    if engine == "numpy":
        numpy_data = NumpyMultiHits.from_parquet(groups, multi_read_paths, batch_size)
        multi_record_count = numpy_data.record_count
        logger.info("Compact multi-hit arrays: %d bytes", numpy_data.nbytes)
    else:
        multi_record_count = load_multi_hit_records(
            groups=groups,
            multi_read_paths=(
                multi_read_paths
            ),
            batch_size=batch_size,
        )

    initial_tpm_by_group = (
        initialize_iteration_zero(
            groups=groups,
            read_length=read_length,
        )
    )

    logger.info(
        "Mapped read count: %d",
        mapped_read_count,
    )
    logger.info(
        "Loaded groups: %d",
        len(groups),
    )
    logger.info(
        "Loaded multi-hit records: %d",
        multi_record_count,
    )
    logger.info(
        "Initial TPM sum: %.6f",
        sum(
            initial_tpm_by_group.values()
        ),
    )

    return (
        groups,
        initial_tpm_by_group,
        mapped_read_count,
        multi_record_count,
        numpy_data,
    )


def load_data(em_input_directory, mapping_log_path, read_length, batch_size, logger):
    """Legacy loader retaining MultiHitRecord objects for Python API callers."""
    return _load_em_data(em_input_directory, mapping_log_path, read_length,
                         batch_size, logger, engine="python")[:4]


def write_group_tpm_table(
    groups: dict[str, GroupEmInfo],
    output_path: str | Path,
) -> None:
    """Write initial and final TPM values for each group."""

    path = Path(
        output_path
    ).resolve()

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        header = [
            "group_id",
            "cds_length",
            "unique_read_count",
            "multi_hit_record_count",
            "final_expected_count",
            "initial_tpm",
            "final_tpm",
        ]

        handle.write(
            "\t".join(
                header
            )
            + "\n"
        )

        for group_id in sorted(
            groups
        ):
            group = groups[
                group_id
            ]

            initial_tpm = (
                group.tpm_by_iteration[0]
                if group.tpm_by_iteration
                else 0.0
            )

            final_tpm = (
                group.tpm_by_iteration[-1]
                if group.tpm_by_iteration
                else 0.0
            )

            final_expected_count = (
                group.expected_count_by_iteration[-1]
                if group.expected_count_by_iteration
                else 0.0
            )

            row = [
                group.group_id,
                str(
                    group.cds_length
                ),
                str(
                    group.unique_read_count
                ),
                str(
                    group.multi_hit_record_count
                    if group.multi_hit_record_count is not None
                    else len(group.multi_hit_records)
                ),
                f"{final_expected_count:.12g}",
                f"{initial_tpm:.12g}",
                f"{final_tpm:.12g}",
            ]

            handle.write(
                "\t".join(
                    row
                )
                + "\n"
            )

def run_group_em(
    em_input_directory: str | Path,
    mapping_log_path: str | Path,
    output_tsv_path: str | Path,
    output_log_path: str | Path,
    trial_count: int = 300,
    convergence_threshold: float = 10,
    read_length: int = 30,
    batch_size: int = 100000,
    engine: str = "numpy",
    keep_iteration_history: bool = False,
) -> EmRunResult:
    """
    Run the complete group-level EM workflow.

    The default NumPy backend keeps compact candidates and per-group record
    counts, not MultiHitRecord objects. Use engine="python" for those objects
    or to reproduce legacy arithmetic. By default histories contain iteration
    zero and the final iteration; keep_iteration_history=True retains all.

    Flow:

        load_data()
        for iteration in 1..trial_count:
            run_em_iteration()
            store results
            break when delta reaches the threshold
        write group TPM table
        return final result
    """

    if engine not in {"numpy", "python"}:
        raise ValueError("engine must be 'numpy' or 'python'.")
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1.")
    if trial_count < 1:
        raise ValueError(
            "trial_count must be at least 1."
        )

    if convergence_threshold < 0.0:
        raise ValueError(
            "convergence_threshold must be non-negative."
        )

    if read_length < 1:
        raise ValueError(
            "read_length must be at least 1."
        )

    logger = create_logger(
        output_log_path
    )
    workflow_started_at = time.perf_counter()
    logger.info("EM engine: %s; keep iteration history: %s", engine, keep_iteration_history)

    logger.info(
        "Group EM started"
    )
    logger.info(
        "EM-input directory: %s",
        Path(
            em_input_directory
        ).resolve(),
    )
    logger.info(
        "Mapping log: %s",
        Path(
            mapping_log_path
        ).resolve(),
    )
    logger.info(
        "Trial count: %d",
        trial_count,
    )
    logger.info(
        "Convergence threshold: %.12g",
        convergence_threshold,
    )
    logger.info(
        "Fixed read length: %d",
        read_length,
    )

    load_started_at = time.perf_counter()
    (
        groups,
        current_tpm_by_group,
        mapped_read_count,
        multi_record_count,
        numpy_data,
    ) = _load_em_data(
        em_input_directory=(
            em_input_directory
        ),
        mapping_log_path=(
            mapping_log_path
        ),
        read_length=read_length,
        batch_size=batch_size,
        logger=logger,
        engine=engine,
    )
    logger.info(
        "Group EM input loading elapsed: %.3f s",
        time.perf_counter() - load_started_at,
    )

    completed_iterations = 0
    converged = False
    final_delta = float(
        "inf"
    )

    for iteration in range(
        1,
        trial_count + 1,
    ):
        iteration_started_at = time.perf_counter()
        (
            expected_count_by_group,
            new_tpm_by_group,
            delta,
        ) = (numpy_data.iteration if numpy_data is not None else run_em_iteration)(
            groups=groups,
            current_tpm_by_group=(
                current_tpm_by_group
            ),
            read_length=read_length,
        )

        append_iteration_results(
            groups=groups,
            expected_count_by_group=(
                expected_count_by_group
            ),
            tpm_by_group=(
                new_tpm_by_group
            ),
            keep_iteration_history=keep_iteration_history,
        )

        current_tpm_by_group = (
            new_tpm_by_group
        )

        completed_iterations = (
            iteration
        )
        final_delta = delta

        logger.info(
            "Iteration=%d DeltaTPMMax=%.12g TPMsum=%.6f Elapsed=%.3fs",
            iteration,
            delta,
            sum(
                current_tpm_by_group.values()
            ),
            time.perf_counter() - iteration_started_at,
        )

        if delta <= convergence_threshold:
            converged = True

            logger.info(
                "Converged at iteration %d",
                iteration,
            )
            break

    output_started_at = time.perf_counter()
    write_group_tpm_table(
        groups=groups,
        output_path=output_tsv_path,
    )
    logger.info(
        "Group EM output writing elapsed: %.3f s",
        time.perf_counter() - output_started_at,
    )

    logger.info(
        "Group EM finished"
    )
    logger.info(
        "Completed iterations: %d",
        completed_iterations,
    )
    logger.info(
        "Converged: %s",
        converged,
    )
    logger.info(
        "Final delta: %.12g",
        final_delta,
    )
    logger.info(
        "Mapped read count: %d",
        mapped_read_count,
    )
    logger.info(
        "Multi-hit record count: %d",
        multi_record_count,
    )
    logger.info(
        "Output table: %s",
        Path(
            output_tsv_path
        ).resolve(),
    )
    logger.info(
        "Group EM total elapsed: %.3f s",
        time.perf_counter() - workflow_started_at,
    )

    return EmRunResult(
        groups=groups,
        final_tpm_by_group=(
            current_tpm_by_group
        ),
        mapped_read_count=(
            mapped_read_count
        ),
        completed_iterations=(
            completed_iterations
        ),
        converged=converged,
        final_delta=final_delta,
    )
