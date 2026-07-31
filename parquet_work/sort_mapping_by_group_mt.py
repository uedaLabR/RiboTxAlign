from __future__ import annotations

import argparse
import re
import time
from collections import defaultdict
from concurrent.futures import (
    ProcessPoolExecutor,
    as_completed,
)
from dataclasses import dataclass
from multiprocessing import (
    get_context,
)
from pathlib import Path


_PART_FILE_PATTERN = re.compile(
    r"^(?P<partition_name>.+)"
    r"\.part(?P<part_index>\d+)\.parquet$"
)


def _print_status(
    message: str,
) -> None:
    """Print one timestamped status message immediately."""

    timestamp = time.strftime(
        "%Y-%m-%d %H:%M:%S"
    )

    print(
        f"{timestamp} {message}",
        flush=True,
    )


@dataclass(slots=True)
class PartitionTask:
    """One chromosome/partition sorting task."""

    partition_name: str
    part_paths: list[str]
    output_path: str
    read_batch_size: int
    row_group_size: int
    overwrite: bool


@dataclass(slots=True)
class PartitionResult:
    """Summary returned from one chromosome worker."""

    partition_name: str
    output_path: str
    total_rows: int
    group_count: int
    elapsed_seconds: float


def find_mapping_parts(
    mapping_parts_directory: str,
) -> dict[
    str,
    list[Path],
]:
    """
    Find mapping part files and group them by partition name.

    This initial file-grouping stage runs in the parent process.

    Examples:

        chr1.part000000.parquet
        chr1.part000001.parquet
        chr2.part000000.parquet
    """

    input_directory = Path(
        mapping_parts_directory
    ).resolve()

    if not input_directory.is_dir():
        raise NotADirectoryError(
            "Mapping-parts directory was not found: "
            f"{input_directory}"
        )

    indexed_parts: dict[
        str,
        list[
            tuple[
                int,
                Path,
            ]
        ],
    ] = defaultdict(
        list
    )

    _print_status(
        "Scanning mapping-parts directory: "
        f"{input_directory}"
    )

    for part_path in input_directory.glob(
        "*.part*.parquet"
    ):
        match = _PART_FILE_PATTERN.fullmatch(
            part_path.name
        )

        if match is None:
            _print_status(
                "Ignored unrecognized file: "
                f"{part_path.name}"
            )
            continue

        partition_name = match.group(
            "partition_name"
        )

        part_index = int(
            match.group(
                "part_index"
            )
        )

        indexed_parts[
            partition_name
        ].append(
            (
                part_index,
                part_path.resolve(),
            )
        )

    if not indexed_parts:
        raise FileNotFoundError(
            "No mapping part Parquet files were found in "
            f"{input_directory}"
        )

    parts_by_partition: dict[
        str,
        list[Path],
    ] = {}

    for partition_name, indexed_paths in (
        indexed_parts.items()
    ):
        indexed_paths.sort(
            key=lambda item: (
                item[0],
                item[1].name,
            )
        )

        parts_by_partition[
            partition_name
        ] = [
            path
            for _, path
            in indexed_paths
        ]

        _print_status(
            f"Detected {partition_name}: "
            f"{len(indexed_paths)} part files"
        )

    return parts_by_partition


def _collect_group_fragments(
    part_paths: list[Path],
    partition_name: str,
    read_batch_size: int,
) -> tuple[
    dict[
        str,
        list,
    ],
    object,
    int,
]:
    """
    Read one chromosome and store Arrow fragments by group ID.

    The dictionary structure is:

        group_id -> list[pyarrow.Table]

    Records are not converted to Python dictionaries.
    """

    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as error:
        raise ImportError(
            "Sorting requires pyarrow. "
            "Install it with: pip install pyarrow"
        ) from error

    group_fragments: dict[
        str,
        list,
    ] = defaultdict(
        list
    )

    schema = None
    total_rows = 0
    partition_start = time.time()

    _print_status(
        f"[{partition_name}] Start reading "
        f"{len(part_paths)} part files"
    )

    for part_number, part_path in enumerate(
        part_paths,
        start=1,
    ):
        part_start = time.time()

        _print_status(
            f"[{partition_name}] Reading part "
            f"{part_number}/{len(part_paths)}: "
            f"{part_path.name}"
        )

        parquet_file = pq.ParquetFile(
            part_path
        )

        if schema is None:
            schema = parquet_file.schema_arrow

        elif parquet_file.schema_arrow != schema:
            raise ValueError(
                "Parquet schemas differ within partition "
                f"{partition_name}: {part_path}"
            )

        part_rows = 0

        for batch_number, batch in enumerate(
            parquet_file.iter_batches(
                batch_size=read_batch_size,
            ),
            start=1,
        ):
            table = pa.Table.from_batches(
                [
                    batch
                ]
            )

            group_ids = table[
                "group_id"
            ].to_pylist()

            row_indices_by_group: dict[
                str,
                list[int],
            ] = defaultdict(
                list
            )

            for row_index, group_id in enumerate(
                group_ids
            ):
                row_indices_by_group[
                    group_id
                ].append(
                    row_index
                )

            for group_id, row_indices in (
                row_indices_by_group.items()
            ):
                fragment = table.take(
                    pa.array(
                        row_indices,
                        type=pa.int64(),
                    )
                )

                group_fragments[
                    group_id
                ].append(
                    fragment
                )

            batch_rows = table.num_rows
            part_rows += batch_rows
            total_rows += batch_rows

            _print_status(
                f"[{partition_name}] "
                f"Part {part_number}/{len(part_paths)}, "
                f"batch {batch_number}: "
                f"rows={batch_rows:,}, "
                f"groups_in_batch={len(row_indices_by_group):,}, "
                f"groups_held={len(group_fragments):,}, "
                f"total_rows={total_rows:,}"
            )

        _print_status(
            f"[{partition_name}] Finished part "
            f"{part_number}/{len(part_paths)}: "
            f"rows={part_rows:,}, "
            f"elapsed={time.time() - part_start:.1f}s"
        )

    if schema is None:
        raise ValueError(
            "No readable Parquet schema found for "
            f"{partition_name}"
        )

    _print_status(
        f"[{partition_name}] Finished reading: "
        f"rows={total_rows:,}, "
        f"group_ids={len(group_fragments):,}, "
        f"elapsed={time.time() - partition_start:.1f}s"
    )

    return (
        dict(
            group_fragments
        ),
        schema,
        total_rows,
    )


def _write_partition_by_group(
    group_fragments: dict[
        str,
        list,
    ],
    schema,
    output_path: Path,
    partition_name: str,
    row_group_size: int,
) -> None:
    """
    Write one chromosome with group IDs in sorted order.

    Records inside each group are not sorted again. Fragments are written
    in their original part/batch order.

    Entire groups are never concatenated.
    """

    try:
        import pyarrow.parquet as pq
    except ImportError as error:
        raise ImportError(
            "Sorting requires pyarrow. "
            "Install it with: pip install pyarrow"
        ) from error

    temporary_output = output_path.with_name(
        f".{output_path.name}.tmp"
    )

    if temporary_output.exists():
        temporary_output.unlink()

    group_ids = sorted(
        group_fragments
    )

    total_group_count = len(
        group_ids
    )

    writer = None
    written_rows = 0
    write_start = time.time()

    _print_status(
        f"[{partition_name}] Opening output: "
        f"{temporary_output}"
    )

    try:
        writer = pq.ParquetWriter(
            temporary_output,
            schema=schema,
            compression="zstd",
        )

        for processed_group_count, group_id in enumerate(
            group_ids,
            start=1,
        ):
            fragments = group_fragments[
                group_id
            ]

            group_row_count = sum(
                fragment.num_rows
                for fragment in fragments
            )

            _print_status(
                f"[{partition_name}] Processing group_id "
                f"{processed_group_count:,}/"
                f"{total_group_count:,}: "
                f"{group_id}, "
                f"fragments={len(fragments):,}, "
                f"records={group_row_count:,}"
            )

            group_start = time.time()

            for fragment_number, fragment in enumerate(
                fragments,
                start=1,
            ):
                writer.write_table(
                    fragment,
                    row_group_size=(
                        row_group_size
                    ),
                )

                written_rows += (
                    fragment.num_rows
                )

                if (
                    len(fragments) > 20
                    and (
                        fragment_number == 1
                        or fragment_number % 20 == 0
                        or fragment_number
                        == len(fragments)
                    )
                ):
                    _print_status(
                        f"[{partition_name}] "
                        f"  group_id={group_id}: "
                        f"fragment "
                        f"{fragment_number:,}/"
                        f"{len(fragments):,}, "
                        f"written_rows={written_rows:,}"
                    )

            _print_status(
                f"[{partition_name}] Finished group_id "
                f"{processed_group_count:,}/"
                f"{total_group_count:,}: "
                f"{group_id}, "
                f"records={group_row_count:,}, "
                f"total_written={written_rows:,}, "
                f"elapsed={time.time() - group_start:.1f}s"
            )

    except Exception:
        _print_status(
            f"[{partition_name}] Error while writing output"
        )

        if writer is not None:
            writer.close()

        temporary_output.unlink(
            missing_ok=True
        )

        raise

    else:
        _print_status(
            f"[{partition_name}] Closing output writer"
        )

        writer.close()

        temporary_output.replace(
            output_path
        )

    _print_status(
        f"[{partition_name}] Finished output: "
        f"{output_path}, "
        f"group_ids={total_group_count:,}, "
        f"rows={written_rows:,}, "
        f"elapsed={time.time() - write_start:.1f}s"
    )


def _process_partition(
    task: PartitionTask,
) -> PartitionResult:
    """
    Collect and write one chromosome in one worker process.
    """

    start_time = time.time()

    output_path = Path(
        task.output_path
    )

    if output_path.exists():
        if task.overwrite:
            output_path.unlink()
        else:
            raise FileExistsError(
                "Sorted output already exists: "
                f"{output_path}"
            )

    part_paths = [
        Path(
            path
        )
        for path in task.part_paths
    ]

    _print_status(
        f"[{task.partition_name}] Worker started"
    )

    group_fragments, schema, total_rows = (
        _collect_group_fragments(
            part_paths=part_paths,
            partition_name=(
                task.partition_name
            ),
            read_batch_size=(
                task.read_batch_size
            ),
        )
    )

    group_count = len(
        group_fragments
    )

    _write_partition_by_group(
        group_fragments=group_fragments,
        schema=schema,
        output_path=output_path,
        partition_name=(
            task.partition_name
        ),
        row_group_size=(
            task.row_group_size
        ),
    )

    del group_fragments

    elapsed = (
        time.time()
        - start_time
    )

    _print_status(
        f"[{task.partition_name}] Worker finished: "
        f"rows={total_rows:,}, "
        f"group_ids={group_count:,}, "
        f"elapsed={elapsed:.1f}s"
    )

    return PartitionResult(
        partition_name=(
            task.partition_name
        ),
        output_path=str(
            output_path
        ),
        total_rows=total_rows,
        group_count=group_count,
        elapsed_seconds=elapsed,
    )


def _delete_intermediate_parquets(
    part_paths: list[Path],
    input_directory: Path,
) -> None:
    """Delete part Parquets after every output succeeds."""

    for file_number, part_path in enumerate(
        part_paths,
        start=1,
    ):
        part_path.unlink()

        if (
            file_number == 1
            or file_number % 100 == 0
            or file_number == len(part_paths)
        ):
            _print_status(
                "Deleted intermediate Parquet "
                f"{file_number:,}/{len(part_paths):,}"
            )

    try:
        input_directory.rmdir()

    except OSError:
        # Keep the directory when logs or unrelated files remain.
        pass


def sort_mapping_by_group(
    mapping_parts_directory: str,
    output_directory: str,
    read_batch_size: int = 100000,
    row_group_size: int = 100000,
    worker_count: int = 4,
    overwrite: bool = False,
    delete_intermediate_parquets: bool = False,
    multiprocessing_start_method: str = "spawn",
) -> list[Path]:
    """
    Group input part files in the parent, then process chromosomes in
    parallel.

    Each worker handles one chromosome:

        read parts
        -> hold fragments in group_id dictionary
        -> write group IDs in sorted order

    SQL and DuckDB are not used.
    """

    if read_batch_size < 1:
        raise ValueError(
            "read_batch_size must be at least 1."
        )

    if row_group_size < 1:
        raise ValueError(
            "row_group_size must be at least 1."
        )

    if worker_count < 1:
        raise ValueError(
            "worker_count must be at least 1."
        )

    input_directory = Path(
        mapping_parts_directory
    ).resolve()

    output_path = Path(
        output_directory
    ).resolve()

    parts_by_partition = find_mapping_parts(
        str(
            input_directory
        )
    )

    output_path.mkdir(
        parents=True,
        exist_ok=True,
    )

    all_part_paths = [
        part_path
        for part_paths
        in parts_by_partition.values()
        for part_path
        in part_paths
    ]

    tasks = [
        PartitionTask(
            partition_name=partition_name,
            part_paths=[
                str(
                    part_path
                )
                for part_path
                in parts_by_partition[
                    partition_name
                ]
            ],
            output_path=str(
                output_path
                / f"{partition_name}.parquet"
            ),
            read_batch_size=read_batch_size,
            row_group_size=row_group_size,
            overwrite=overwrite,
        )
        for partition_name
        in sorted(
            parts_by_partition
        )
    ]

    actual_worker_count = min(
        worker_count,
        len(tasks),
    )

    _print_status(
        f"Started multiprocess group ordering: "
        f"partitions={len(tasks):,}, "
        f"part_files={len(all_part_paths):,}, "
        f"workers={actual_worker_count:,}"
    )

    context = get_context(
        multiprocessing_start_method
    )

    results: list[
        PartitionResult
    ] = []

    with ProcessPoolExecutor(
        max_workers=actual_worker_count,
        mp_context=context,
    ) as executor:
        future_to_partition = {
            executor.submit(
                _process_partition,
                task,
            ): task.partition_name
            for task in tasks
        }

        for completed_count, future in enumerate(
            as_completed(
                future_to_partition
            ),
            start=1,
        ):
            partition_name = (
                future_to_partition[
                    future
                ]
            )

            try:
                result = future.result()

            except Exception:
                _print_status(
                    f"[{partition_name}] Worker failed"
                )

                for pending_future in (
                    future_to_partition
                ):
                    pending_future.cancel()

                raise

            results.append(
                result
            )

            _print_status(
                f"Completed partition "
                f"{completed_count:,}/{len(tasks):,}: "
                f"{result.partition_name}, "
                f"rows={result.total_rows:,}, "
                f"group_ids={result.group_count:,}, "
                f"elapsed={result.elapsed_seconds:.1f}s"
            )

    results.sort(
        key=lambda result: (
            result.partition_name
        )
    )

    output_files = [
        Path(
            result.output_path
        )
        for result in results
    ]

    if delete_intermediate_parquets:
        _print_status(
            "Deleting intermediate mapping Parquets"
        )

        _delete_intermediate_parquets(
            part_paths=all_part_paths,
            input_directory=(
                input_directory
            ),
        )

    _print_status(
        f"Completed multiprocess group ordering: "
        f"outputs={len(output_files):,}"
    )

    return output_files