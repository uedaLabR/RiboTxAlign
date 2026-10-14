from __future__ import annotations

import multiprocessing as mp
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from em_distribute.group_cds_definition import (
    GroupCdsDefinitions,
)


_WORKER_GROUP_DEFINITIONS: (
    GroupCdsDefinitions
    | None
) = None


MULTI_SCHEMA = pa.schema(
    [
        ("group_id", pa.string()),
        ("read_id", pa.string()),
        ("group_index", pa.int32()),
        (
            "other_group_ids",
            pa.list_(pa.string()),
        ),
    ]
)


GROUP_INFO_SCHEMA = pa.schema(
    [
        ("group_id", pa.string()),
        ("cds_length", pa.int64()),
        (
            "group_unique_read_count",
            pa.int64(),
        ),
    ]
)


def initialize_worker(
    group_definitions: GroupCdsDefinitions,
) -> None:
    """Store group definitions in the current process."""

    global _WORKER_GROUP_DEFINITIONS

    _WORKER_GROUP_DEFINITIONS = (
        group_definitions
    )


def read_groups(
    parquet_path: str | Path,
    batch_size: int = 100000,
):
    """
    Yield one complete group from a group-sorted Parquet file.

    The final group in an Arrow batch is retained and continued in the
    following batch.
    """

    parquet_file = pq.ParquetFile(
        parquet_path
    )

    columns = [
        "read_id",
        "transcript_id",
        "group_id",
        "group_index",
        "group_count",
        "other_group_ids",
        "reference_start",
        "reference_end",
    ]

    current_group_id: str | None = None
    current_rows: list[dict] = []

    for batch in parquet_file.iter_batches(
        batch_size=batch_size,
        columns=columns,
    ):
        for row in batch.to_pylist():
            group_id = row[
                "group_id"
            ]

            if current_group_id is None:
                current_group_id = group_id

            elif group_id != current_group_id:
                yield (
                    current_group_id,
                    current_rows,
                )

                current_group_id = group_id
                current_rows = []

            current_rows.append(row)

    if current_group_id is not None:
        yield (
            current_group_id,
            current_rows,
        )


def union_interval_length(
    intervals: set[
        tuple[
            str,
            str,
            int,
            int,
        ]
    ],
) -> int:
    """
    Calculate non-redundant inclusive CDS length.

    Intervals on different chromosomes or strands are merged separately.
    """

    intervals_by_location: dict[
        tuple[
            str,
            str,
        ],
        list[
            tuple[
                int,
                int,
            ]
        ],
    ] = {}

    for (
        chromosome,
        strand,
        start,
        end,
    ) in intervals:
        intervals_by_location.setdefault(
            (
                chromosome,
                strand,
            ),
            [],
        ).append(
            (
                start,
                end,
            )
        )

    total_length = 0

    for location_intervals in (
        intervals_by_location.values()
    ):
        location_intervals.sort(
            key=lambda interval: (
                interval[0],
                interval[1],
            )
        )

        current_start, current_end = (
            location_intervals[0]
        )

        for start, end in (
            location_intervals[1:]
        ):
            if start <= current_end + 1:
                current_end = max(
                    current_end,
                    end,
                )
                continue

            total_length += (
                current_end
                - current_start
                + 1
            )

            current_start = start
            current_end = end

        total_length += (
            current_end
            - current_start
            + 1
        )

    return total_length



def determine_group_length(
    rows: list[dict],
    reference_definitions: dict,
    hit_transcript_ids: set[str],
    minimum_mapped_reads: int,
    minimum_mapped_span: int,
    minimum_transcript_fraction: float,
) -> int:
    """
    Determine the group length from CDS annotation or mapped-read span.

    Priority:

        1. Non-redundant annotated CDS length.
        2. Mapped span when annotation is absent and support is sufficient.
        3. Maximum transcript length among hit references.
    """

    cds_intervals: set[
        tuple[str, str, int, int]
    ] = set()

    maximum_transcript_length = 0

    for transcript_id in hit_transcript_ids:
        definition = reference_definitions.get(
            transcript_id
        )

        if definition is None:
            continue

        cds_intervals.update(
            definition.cds_intervals
        )

        maximum_transcript_length = max(
            maximum_transcript_length,
            int(
                definition.transcript_length
            ),
        )

    annotated_cds_length = (
        union_interval_length(
            cds_intervals
        )
    )

    if annotated_cds_length > 0:
        return annotated_cds_length

    distinct_read_count = len(
        {
            row["read_id"]
            for row in rows
        }
    )

    if rows:
        mapped_start = min(
            int(
                row["reference_start"]
            )
            for row in rows
        )

        mapped_end = max(
            int(
                row["reference_end"]
            )
            for row in rows
        )

        mapped_span = max(
            mapped_end - mapped_start,
            1,
        )
    else:
        mapped_span = 0

    comparison_length = max(
        maximum_transcript_length,
        1,
    )

    use_mapped_span = (
        distinct_read_count
        > minimum_mapped_reads
        and mapped_span
        >= minimum_mapped_span
        and mapped_span
        >= comparison_length
        * minimum_transcript_fraction
    )

    if use_mapped_span:
        return mapped_span

    return comparison_length

def extract_group(
    group_id: str,
    rows: list[dict],
    minimum_mapped_reads: int,
    minimum_mapped_span: int,
    minimum_transcript_fraction: float,
) -> tuple[
    dict,
    list[dict],
]:
    """
    Extract final group information and multi-read records.

    Only mapping references actually hit in this group contribute CDS
    intervals.
    """

    if _WORKER_GROUP_DEFINITIONS is None:
        raise RuntimeError(
            "Worker group definitions were not initialized."
        )

    unique_read_ids: set[str] = set()
    hit_transcript_ids: set[str] = set()

    multi_records_by_candidate: dict[
        tuple[
            str,
            int,
        ],
        dict,
    ] = {}

    for row in rows:
        read_id = row["read_id"]
        transcript_id = row["transcript_id"]
        group_index = int(
            row["group_index"]
        )
        group_count = int(
            row["group_count"]
        )

        hit_transcript_ids.add(
            transcript_id
        )

        if group_count == 1:
            unique_read_ids.add(
                read_id
            )
            continue

        candidate_key = (
            read_id,
            group_index,
        )

        if candidate_key in (
            multi_records_by_candidate
        ):
            continue

        other_group_ids = row.get(
            "other_group_ids"
        )

        if other_group_ids is None:
            other_group_ids = []

        multi_records_by_candidate[
            candidate_key
        ] = {
            "group_id": group_id,
            "read_id": read_id,
            "group_index": group_index,
            "other_group_ids": list(
                other_group_ids
            ),
        }

    reference_definitions = (
        _WORKER_GROUP_DEFINITIONS.get(
            group_id,
            {},
        )
    )

    group_length = determine_group_length(
        rows=rows,
        reference_definitions=(
            reference_definitions
        ),
        hit_transcript_ids=(
            hit_transcript_ids
        ),
        minimum_mapped_reads=(
            minimum_mapped_reads
        ),
        minimum_mapped_span=(
            minimum_mapped_span
        ),
        minimum_transcript_fraction=(
            minimum_transcript_fraction
        ),
    )

    group_record = {
        "group_id": group_id,
        "cds_length": group_length,
        "group_unique_read_count": len(
            unique_read_ids
        ),
    }

    return (
        group_record,
        list(
            multi_records_by_candidate.values()
        ),
    )


def write_buffer(
    writer: pq.ParquetWriter,
    rows: list[dict],
    schema: pa.Schema,
) -> None:
    """Write and clear one row buffer."""

    if not rows:
        return

    table = pa.Table.from_pylist(
        rows,
        schema=schema,
    )

    writer.write_table(table)
    rows.clear()


def process_file_worker(
    parquet_path: str,
    output_directory: str,
    batch_size: int,
    write_batch_size: int,
    overwrite: bool,
    minimum_mapped_reads: int,
    minimum_mapped_span: int,
    minimum_transcript_fraction: float,
) -> tuple[
    str,
    str,
]:
    """
    Process one sorted mapping Parquet.

    Final output shards:

        <stem>.multi_reads.parquet
        <stem>.group_info.parquet
    """

    input_path = Path(parquet_path)
    output_path = Path(output_directory)

    multi_output_path = (
        output_path
        / f"{input_path.stem}.multi_reads.parquet"
    )

    group_output_path = (
        output_path
        / f"{input_path.stem}.group_info.parquet"
    )

    for file_path in (
        multi_output_path,
        group_output_path,
    ):
        if file_path.exists():
            if overwrite:
                file_path.unlink()
            else:
                raise FileExistsError(
                    f"Output already exists: {file_path}"
                )

    multi_writer = pq.ParquetWriter(
        multi_output_path,
        schema=MULTI_SCHEMA,
        compression="zstd",
    )

    group_writer = pq.ParquetWriter(
        group_output_path,
        schema=GROUP_INFO_SCHEMA,
        compression="zstd",
    )

    multi_buffer: list[dict] = []
    group_buffer: list[dict] = []

    try:
        for group_id, rows in read_groups(
            parquet_path=input_path,
            batch_size=batch_size,
        ):
            (
                group_record,
                multi_records,
            ) = extract_group(
                group_id=group_id,
                rows=rows,
                minimum_mapped_reads=(
                    minimum_mapped_reads
                ),
                minimum_mapped_span=(
                    minimum_mapped_span
                ),
                minimum_transcript_fraction=(
                    minimum_transcript_fraction
                ),
            )

            group_buffer.append(
                group_record
            )

            multi_buffer.extend(
                multi_records
            )

            if (
                len(group_buffer)
                >= write_batch_size
            ):
                write_buffer(
                    writer=group_writer,
                    rows=group_buffer,
                    schema=GROUP_INFO_SCHEMA,
                )

            if (
                len(multi_buffer)
                >= write_batch_size
            ):
                write_buffer(
                    writer=multi_writer,
                    rows=multi_buffer,
                    schema=MULTI_SCHEMA,
                )

        write_buffer(
            writer=group_writer,
            rows=group_buffer,
            schema=GROUP_INFO_SCHEMA,
        )

        write_buffer(
            writer=multi_writer,
            rows=multi_buffer,
            schema=MULTI_SCHEMA,
        )

    except Exception:
        group_writer.close()
        multi_writer.close()

        group_output_path.unlink(
            missing_ok=True
        )
        multi_output_path.unlink(
            missing_ok=True
        )

        raise

    else:
        group_writer.close()
        multi_writer.close()

    return (
        str(multi_output_path),
        str(group_output_path),
    )


def _run_one_file_process(
    parquet_path: str,
    output_directory: str,
    batch_size: int,
    write_batch_size: int,
    overwrite: bool,
    minimum_mapped_reads: int,
    minimum_mapped_span: int,
    minimum_transcript_fraction: float,
    status_directory: str,
) -> None:
    """
    Process one input file and terminate without Python cleanup.

    A status file is written before os._exit() so the parent can
    distinguish successful completion from failure.
    """

    import os
    import traceback

    input_path = Path(
        parquet_path
    )

    status_path = (
        Path(status_directory)
        / f"{input_path.stem}.status"
    )

    try:
        process_file_worker(
            parquet_path=parquet_path,
            output_directory=output_directory,
            batch_size=batch_size,
            write_batch_size=write_batch_size,
            overwrite=overwrite,
            minimum_mapped_reads=(
                minimum_mapped_reads
            ),
            minimum_mapped_span=(
                minimum_mapped_span
            ),
            minimum_transcript_fraction=(
                minimum_transcript_fraction
            ),
        )

        status_path.write_text(
            "OK\n",
            encoding="utf-8",
        )

        # Parquet writers are already closed. Skip Python/PyArrow
        # interpreter cleanup, which can hang in worker processes.
        os._exit(0)

    except BaseException:
        status_path.write_text(
            traceback.format_exc(),
            encoding="utf-8",
        )

        os._exit(1)


def extract_em_input(
    sorted_mapping_directory: str,
    output_directory: str,
    group_definitions: GroupCdsDefinitions,
    batch_size: int = 100000,
    write_batch_size: int = 100000,
    worker_count: int = 4,
    overwrite: bool = False,
    minimum_mapped_reads: int = 10,
    minimum_mapped_span: int = 100,
    minimum_transcript_fraction: float = (
        1.0 / 3.0
    ),
) -> None:
    """
    Create final EM-input shards with bounded multiprocessing.

    Each child processes exactly one input Parquet and calls os._exit()
    immediately after closing its output writers. No Pool or executor is
    used, and no worker result is returned.
    """

    import multiprocessing as mp
    import shutil
    import tempfile
    import time

    sorted_directory = Path(
        sorted_mapping_directory
    ).resolve()

    output_path = Path(
        output_directory
    ).resolve()

    if output_path == sorted_directory:
        raise ValueError("Sorted-mapping and EM-input directories must be different.")

    if not sorted_directory.is_dir():
        raise NotADirectoryError(
            "Sorted mapping directory was not found: "
            f"{sorted_directory}"
        )

    parquet_paths = sorted(
        sorted_directory.glob(
            "*.parquet"
        )
    )

    if not parquet_paths:
        raise FileNotFoundError(
            "No sorted mapping Parquet files were found in "
            f"{sorted_directory}"
        )

    if worker_count < 1:
        raise ValueError(
            "worker_count must be at least 1."
        )

    output_path.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Remove all shards owned by this stage, including chromosomes absent
    # from the new sorted inputs. Per-worker replacement alone leaves these
    # old shards visible to group EM on a rerun.
    previous_outputs = sorted([
        *output_path.glob("*.group_info.parquet"),
        *output_path.glob("*.multi_reads.parquet"),
    ])
    if previous_outputs and not overwrite:
        raise FileExistsError(
            f"EM-input directory already contains {len(previous_outputs)} shard(s), "
            f"including: {previous_outputs[0]}. Use a fresh output directory or set "
            "group_input.overwrite: true (overwrite=True in the Python API) for a full restart."
        )
    if previous_outputs:
        print(f"Removing {len(previous_outputs):,} previous EM-input shards for a fresh run", flush=True)
        for previous_output in previous_outputs:
            previous_output.unlink()

    # Make the large annotation object available before fork.
    initialize_worker(
        group_definitions
    )

    context = mp.get_context(
        "fork"
    )

    status_directory = Path(
        tempfile.mkdtemp(
            prefix="em_input_status_",
            dir=output_path,
        )
    )

    pending_paths = list(
        parquet_paths
    )

    running: dict[
        mp.Process,
        Path,
    ] = {}

    completed_count = 0

    try:
        while pending_paths or running:
            while (
                pending_paths
                and len(running) < worker_count
            ):
                parquet_path = pending_paths.pop(
                    0
                )

                process = context.Process(
                    target=_run_one_file_process,
                    args=(
                        str(parquet_path),
                        str(output_path),
                        batch_size,
                        write_batch_size,
                        overwrite,
                        minimum_mapped_reads,
                        minimum_mapped_span,
                        minimum_transcript_fraction,
                        str(status_directory),
                    ),
                )

                process.start()

                running[
                    process
                ] = parquet_path

                print(
                    "Started:",
                    parquet_path.name,
                    f"pid={process.pid}",
                    flush=True,
                )

            finished_processes = [
                process
                for process
                in running
                if not process.is_alive()
            ]

            if not finished_processes:
                time.sleep(
                    0.2
                )
                continue

            for process in finished_processes:
                parquet_path = running.pop(
                    process
                )

                process.join(
                    timeout=1
                )

                status_path = (
                    status_directory
                    / f"{parquet_path.stem}.status"
                )

                status_text = (
                    status_path.read_text(
                        encoding="utf-8",
                    )
                    if status_path.exists()
                    else ""
                )

                if (
                    process.exitcode != 0
                    or not status_text.startswith(
                        "OK"
                    )
                ):
                    raise RuntimeError(
                        "Worker failed for "
                        f"{parquet_path}\n"
                        f"exitcode={process.exitcode}\n"
                        f"{status_text}"
                    )

                completed_count += 1

                print(
                    "Finished:",
                    parquet_path.name,
                    f"({completed_count}/{len(parquet_paths)})",
                    flush=True,
                )

        print(
            "All files finished.",
            flush=True,
        )

    except Exception:
        for process in running:
            if process.is_alive():
                process.terminate()

        for process in running:
            process.join(
                timeout=5
            )

        raise

    finally:
        shutil.rmtree(
            status_directory,
            ignore_errors=True,
        )
