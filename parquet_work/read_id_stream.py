"""Disk-backed grouping of candidates without retaining a sample's read IDs.

Hash partitions keep each read together. Bounded sorted runs and a bounded
merge fan-in also handle uneven partitions. Rows of a read retain input order,
including across input files, so posterior sums and representative hits agree
with the in-memory implementation.
"""
from contextlib import ExitStack, closing, contextmanager
import heapq
from pathlib import Path
import tempfile
import zlib

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from parquet_work.buffered_writer import BufferedParquetWriter


INPUT_ORDER = "_assignment_input_order"


def _partition_batches(batches, directory, partition_count, logger):
    paths = {}
    ordinal = 0
    next_progress = 1_000_000
    with ExitStack() as stack:
        writers = {}
        for batch in batches:
            table = pa.Table.from_batches([batch])
            size = table.num_rows
            if not size:
                continue
            table = table.append_column(INPUT_ORDER, pa.array(
                np.arange(ordinal, ordinal + size, dtype=np.int64)))
            ordinal += size
            buckets = np.fromiter(
                (zlib.crc32(str(r).encode("utf-8")) % partition_count
                 for r in table["read_id"].to_pylist()), dtype=np.int32, count=size)
            order = np.argsort(buckets, kind="stable")
            ordered_buckets = buckets[order]
            boundaries = np.flatnonzero(ordered_buckets[1:] != ordered_buckets[:-1]) + 1
            for indices in np.split(order, boundaries):
                bucket = int(buckets[indices[0]])
                part = table.take(pa.array(indices))
                if bucket not in writers:
                    paths[bucket] = directory / f"bucket-{bucket:03d}.parquet"
                    writers[bucket] = stack.enter_context(
                        pq.ParquetWriter(paths[bucket], table.schema, compression="zstd"))
                writers[bucket].write_table(part)
            if ordinal >= next_progress:
                logger.info("Read grouping spill progress: rows=%d", ordinal)
                next_progress = (ordinal // 1_000_000 + 1) * 1_000_000
    logger.info("Read grouping spill completed: rows=%d partitions=%d", ordinal, len(paths))
    return [paths[bucket] for bucket in sorted(paths)]


def _rows(path, batch_size):
    with pq.ParquetFile(path) as parquet:
        for batch in parquet.iter_batches(batch_size=batch_size):
            yield from batch.to_pylist()


def _key(row):
    return str(row["read_id"]), row[INPUT_ORDER]


def _merge_rows(paths, batch_size):
    with ExitStack() as stack:
        streams = [stack.enter_context(closing(_rows(path, batch_size))) for path in paths]
        yield from heapq.merge(*streams, key=_key)


def _sorted_rows(path, directory, run_rows, merge_fan_in):
    # At most run_rows Arrow rows are sorted at once, including skewed buckets.
    runs = []
    with pq.ParquetFile(path) as parquet:
        schema = parquet.schema_arrow
        if parquet.metadata.num_rows <= run_rows:
            # Most small partitions fit within the existing row budget.
            # Sort in Arrow and stream once, avoiding a write/read round trip
            # and a Python heap merge. Ties still follow original input order.
            table = parquet.read()
            indices = pc.sort_indices(table, sort_keys=[("read_id", "ascending"),
                                                        (INPUT_ORDER, "ascending")])
            ordered = table.take(indices)
            del table, indices
            for batch in ordered.to_batches(max_chunksize=min(2048, run_rows)):
                yield from batch.to_pylist()
            return
        for index, batch in enumerate(parquet.iter_batches(batch_size=run_rows)):
            table = pa.Table.from_batches([batch])
            indices = pc.sort_indices(table, sort_keys=[("read_id", "ascending"),
                                                        (INPUT_ORDER, "ascending")])
            run = directory / f"{path.stem}.run-{index:06d}.parquet"
            pq.write_table(table.take(indices), run, compression="zstd")
            runs.append(run)
        # Do not keep the final sort buffers during the merge.
        if runs:
            del batch, table, indices
    path.unlink()
    merge_batch = max(1, min(2048, run_rows // merge_fan_in))
    level = 0
    while len(runs) > merge_fan_in:
        next_runs = []
        for start in range(0, len(runs), merge_fan_in):
            group = runs[start:start + merge_fan_in]
            output = directory / f"{path.stem}.merge-{level}-{start}.parquet"
            with BufferedParquetWriter(output, schema, batch_size=run_rows) as writer:
                with closing(_merge_rows(group, merge_batch)) as rows:
                    for row in rows:
                        writer.append(row)
            for run in group:
                run.unlink()
            next_runs.append(output)
        runs = next_runs
        level += 1
    with closing(_merge_rows(runs, merge_batch)) as rows:
        yield from rows
    for run in runs:
        run.unlink()


@contextmanager
def read_id_sorted_rows(batches, temporary_parent, batch_size, logger,
                        partition_count=64, merge_fan_in=32):
    """Yield a bounded iterator with each read contiguous, then remove spills.

Different partitions need not be globally sorted by ID. All occurrences of
one ID enter the same partition, even across input files and record batches.
"""
    if batch_size < 1 or partition_count < 1 or merge_fan_in < 2:
        raise ValueError("Invalid read grouping buffer/partition/merge limits")
    with tempfile.TemporaryDirectory(prefix=".assignment-read-groups-",
                                     dir=temporary_parent) as temporary:
        directory = Path(temporary)
        paths = _partition_batches(batches, directory, partition_count, logger)

        def rows():
            for index, path in enumerate(paths):
                with closing(_sorted_rows(path, directory, batch_size, merge_fan_in)) as stream:
                    yield from stream
                logger.info("Read grouping partition completed: %d/%d", index + 1, len(paths))

        with closing(rows()) as stream:
            yield stream
