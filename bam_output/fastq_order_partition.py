from __future__ import annotations

import gzip
import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


@dataclass(frozen=True, slots=True)
class PartitionInfo:
    """Describe one contiguous FASTQ-order assignment partition."""

    partition_id: int
    first_fastq_index: int
    last_fastq_index: int
    fastq_read_count: int
    assignment_count: int
    parquet_path: Path


def open_text_auto(path: str | Path) -> TextIO:
    """Open plain-text or gzip-compressed input in text mode."""

    input_path = Path(path)

    if input_path.suffix.lower() == ".gz":
        return gzip.open(input_path, "rt", encoding="utf-8", newline="")

    return input_path.open("r", encoding="utf-8", newline="")


def normalize_fastq_read_id(header: str) -> str:
    """Return the first whitespace-delimited FASTQ header token."""

    value = header.rstrip("\r\n")

    if not value.startswith("@"):
        raise ValueError(f"Invalid FASTQ header: {value[:100]}")

    read_id = value[1:].split(maxsplit=1)[0]

    if not read_id:
        raise ValueError("FASTQ header contains an empty read ID.")

    return read_id


def iter_fastq_read_ids(fastq_path: str | Path) -> Iterator[str]:
    """Yield FASTQ read IDs in file order."""

    with open_text_auto(fastq_path) as handle:
        record_index = 0

        while True:
            header = handle.readline()

            if not header:
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise ValueError(
                    "Truncated FASTQ record at zero-based record index "
                    f"{record_index}."
                )

            if not plus.startswith("+"):
                raise ValueError(
                    "Invalid FASTQ plus line at zero-based record index "
                    f"{record_index}."
                )

            yield normalize_fastq_read_id(header)
            record_index += 1


def build_fastq_read_order(
    fastq_path: str | Path,
    logger: logging.Logger | None = None,
) -> dict[str, int]:
    """Read FASTQ once and return read ID to zero-based FASTQ index."""

    read_order: dict[str, int] = {}

    for fastq_index, read_id in enumerate(iter_fastq_read_ids(fastq_path)):
        if read_id in read_order:
            raise ValueError(f"Duplicate FASTQ read ID: {read_id}")

        read_order[read_id] = fastq_index

        if logger is not None and (fastq_index + 1) % 1_000_000 == 0:
            logger.info("FASTQ read IDs loaded: %d", fastq_index + 1)

    if not read_order:
        raise RuntimeError("FASTQ contains no reads.")

    if logger is not None:
        logger.info("FASTQ read IDs loaded in total: %d", len(read_order))

    return read_order


def _get_partition_id(
    fastq_index: int,
    reads_per_partition: int,
    partition_count: int,
) -> int:
    """Return a contiguous partition ID from one FASTQ index."""

    return min(
        fastq_index // reads_per_partition,
        partition_count - 1,
    )


def partition_assignments_by_fastq_order(
    assignment_parquet: str | Path,
    fastq_path: str | Path,
    output_dir: str | Path,
    output_prefix: str = "final_assignments",
    partition_count: int = 10,
    batch_size: int = 250_000,
    compression: str = "zstd",
    logger: logging.Logger | None = None,
) -> list[PartitionInfo]:
    """Sort assignment rows by FASTQ order and write contiguous Parquet parts.

    The input assignment table must contain exactly one selected genomic
    alignment row per read and a ``read_id`` column. All existing columns are
    preserved. ``fastq_index`` and ``partition_id`` are appended.
    """

    try:
        import pyarrow as pa
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
    except ImportError as error:
        raise RuntimeError(
            "FASTQ-order Parquet partitioning requires pyarrow."
        ) from error

    if partition_count <= 0:
        raise ValueError("partition_count must be positive.")

    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")

    assignment_path = Path(assignment_parquet)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)

    read_order = build_fastq_read_order(
        fastq_path=fastq_path,
        logger=logger,
    )
    total_fastq_reads = len(read_order)
    reads_per_partition = math.ceil(total_fastq_reads / partition_count)

    parquet_file = pq.ParquetFile(assignment_path)

    if "read_id" not in parquet_file.schema_arrow.names:
        raise ValueError("Assignment Parquet does not contain read_id.")

    temp_paths = [
        destination / f".{output_prefix}.part{part_id:02d}.unsorted.parquet"
        for part_id in range(partition_count)
    ]
    final_paths = [
        destination / f"{output_prefix}.part{part_id:02d}.parquet"
        for part_id in range(partition_count)
    ]

    writers: list[object | None] = [None] * partition_count
    assignment_counts = [0] * partition_count
    missing_read_ids: list[str] = []
    output_schema = None

    try:
        for batch_number, batch in enumerate(
            parquet_file.iter_batches(batch_size=batch_size),
            start=1,
        ):
            table = pa.Table.from_batches([batch])
            read_ids = table.column("read_id").to_pylist()

            fastq_indices: list[int] = []
            partition_ids: list[int] = []
            keep_indices: list[int] = []

            for row_index, read_id in enumerate(read_ids):
                fastq_index = read_order.get(read_id)

                if fastq_index is None:
                    if len(missing_read_ids) < 20:
                        missing_read_ids.append(str(read_id))
                    continue

                partition_id = _get_partition_id(
                    fastq_index=fastq_index,
                    reads_per_partition=reads_per_partition,
                    partition_count=partition_count,
                )
                keep_indices.append(row_index)
                fastq_indices.append(fastq_index)
                partition_ids.append(partition_id)

            if not keep_indices:
                continue

            filtered = table.take(pa.array(keep_indices, type=pa.int64()))
            filtered = filtered.append_column(
                "fastq_index",
                pa.array(fastq_indices, type=pa.int64()),
            )
            filtered = filtered.append_column(
                "partition_id",
                pa.array(partition_ids, type=pa.int16()),
            )

            if output_schema is None:
                output_schema = filtered.schema

            for partition_id in sorted(set(partition_ids)):
                mask = pc.equal(
                    filtered.column("partition_id"),
                    pa.scalar(partition_id, type=pa.int16()),
                )
                partition_table = filtered.filter(mask)

                if partition_table.num_rows == 0:
                    continue

                if writers[partition_id] is None:
                    writers[partition_id] = pq.ParquetWriter(
                        temp_paths[partition_id],
                        output_schema,
                        compression=compression,
                    )

                writers[partition_id].write_table(partition_table)
                assignment_counts[partition_id] += partition_table.num_rows

            if logger is not None:
                logger.info(
                    "Assignment batches partitioned: %d; rows processed: %d",
                    batch_number,
                    batch_number * batch_size,
                )
    finally:
        for writer in writers:
            if writer is not None:
                writer.close()

    if missing_read_ids:
        raise RuntimeError(
            "Assignment read IDs were not found in FASTQ. Examples: "
            + ", ".join(missing_read_ids)
        )

    partition_info: list[PartitionInfo] = []

    for partition_id in range(partition_count):
        first_fastq_index = partition_id * reads_per_partition
        last_fastq_index = min(
            (partition_id + 1) * reads_per_partition,
            total_fastq_reads,
        )
        fastq_read_count = max(last_fastq_index - first_fastq_index, 0)

        temp_path = temp_paths[partition_id]
        final_path = final_paths[partition_id]

        if temp_path.exists():
            table = pq.read_table(temp_path)
            sort_indices = pc.sort_indices(
                table,
                sort_keys=[("fastq_index", "ascending")],
            )
            sorted_table = table.take(sort_indices)
            pq.write_table(
                sorted_table,
                final_path,
                compression=compression,
            )
            temp_path.unlink()
        else:
            if output_schema is None:
                raise RuntimeError("No assignment rows were written.")

            empty_arrays = [
                pa.array([], type=field.type)
                for field in output_schema
            ]
            empty_table = pa.Table.from_arrays(
                empty_arrays,
                schema=output_schema,
            )
            pq.write_table(
                empty_table,
                final_path,
                compression=compression,
            )

        partition_info.append(
            PartitionInfo(
                partition_id=partition_id,
                first_fastq_index=first_fastq_index,
                last_fastq_index=last_fastq_index,
                fastq_read_count=fastq_read_count,
                assignment_count=assignment_counts[partition_id],
                parquet_path=final_path,
            )
        )

        if logger is not None:
            logger.info(
                "Partition %02d: FASTQ [%d, %d), assignments=%d, path=%s",
                partition_id,
                first_fastq_index,
                last_fastq_index,
                assignment_counts[partition_id],
                final_path,
            )

    manifest_path = destination / f"{output_prefix}.partitions.tsv"

    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(
            "partition_id\tfirst_fastq_index\tlast_fastq_index\t"
            "fastq_read_count\tassignment_count\tparquet_path\n"
        )

        for info in partition_info:
            handle.write(
                f"{info.partition_id}\t"
                f"{info.first_fastq_index}\t"
                f"{info.last_fastq_index}\t"
                f"{info.fastq_read_count}\t"
                f"{info.assignment_count}\t"
                f"{info.parquet_path}\n"
            )

    return partition_info