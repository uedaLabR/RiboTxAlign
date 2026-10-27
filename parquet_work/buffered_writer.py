"""Bounded, atomic Parquet output for append-only read records."""

from pathlib import Path
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq


class BufferedParquetWriter:
    """Keep at most one batch of Python rows; publish on successful exit.

    Row order and the supplied schema are preserved. The temporary file is
    beside the destination so a failed stage leaves any previous output intact.
    """

    def __init__(self, output_path, schema, batch_size=100_000, defaults=None):
        if batch_size < 1:
            raise ValueError("batch_size must be at least 1")
        self.output_path = Path(output_path)
        self.schema = schema
        self.batch_size = batch_size
        self.defaults = defaults or {}
        self.row_count = 0
        self._rows = []
        self._writer = None
        self._temporary_path = None

    @property
    def buffered_row_count(self):
        return len(self._rows)

    def __enter__(self):
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            dir=self.output_path.parent,
            prefix=f".{self.output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            self._temporary_path = Path(handle.name)
        try:
            self._writer = pq.ParquetWriter(
                self._temporary_path, self.schema, compression="zstd"
            )
        except BaseException:
            self._temporary_path.unlink(missing_ok=True)
            raise
        return self

    def append(self, row):
        if self._writer is None:
            raise RuntimeError("Use BufferedParquetWriter as a context manager")
        normalized = dict(row)
        for name, value in self.defaults.items():
            normalized.setdefault(name, value)
        self._rows.append(normalized)
        self.row_count += 1
        if len(self._rows) >= self.batch_size:
            self.flush()

    def flush(self):
        if self._rows:
            table = pa.Table.from_pylist(self._rows, schema=self.schema)
            self._writer.write_table(table)
            self._rows.clear()

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            try:
                if exc_type is None:
                    self.flush()
            finally:
                if self._writer is not None:
                    self._writer.close()
                    self._writer = None
            if exc_type is None:
                self._temporary_path.replace(self.output_path)
        finally:
            self._rows.clear()
            self._temporary_path.unlink(missing_ok=True)
        return False
