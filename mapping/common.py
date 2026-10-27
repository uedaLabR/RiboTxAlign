from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import mappy as mp


@dataclass(slots=True)
class ReadRecord:
    read_id: str
    sequence: str
    quality: str | None = None


@dataclass(slots=True)
class MappingCounters:
    total_reads: int = 0
    mapped_reads: int = 0
    group_unique_reads: int = 0
    group_multi_reads: int = 0
    unmapped_reads: int = 0
    best_hit_overflow_reads: int = 0
    primary_preferred_alt_tie_reads: int = 0
    removed_alt_tie_hits: int = 0

    def add(self, other: "MappingCounters") -> None:
        self.total_reads += other.total_reads
        self.mapped_reads += other.mapped_reads
        self.group_unique_reads += other.group_unique_reads
        self.group_multi_reads += other.group_multi_reads
        self.unmapped_reads += other.unmapped_reads
        self.best_hit_overflow_reads += other.best_hit_overflow_reads
        self.primary_preferred_alt_tie_reads += other.primary_preferred_alt_tie_reads
        self.removed_alt_tie_hits += other.removed_alt_tie_hits


def read_fastq_batches(
    fastq_path: str,
    batch_size: int,
    maximum_reads: int | None,
) -> Iterator[list[ReadRecord]]:
    """Read FASTQ records in bounded batches."""

    batch: list[ReadRecord] = []
    total = 0

    for read_id, sequence, quality in mp.fastx_read(fastq_path):
        batch.append(
            ReadRecord(
                read_id=read_id,
                sequence=sequence,
                quality=quality,
            )
        )
        total += 1

        if len(batch) >= batch_size:
            yield batch
            batch = []

        if maximum_reads is not None and total >= maximum_reads:
            break

    if batch:
        yield batch


def create_logger(name: str, log_file: str) -> logging.Logger:
    """Create a file and console logger."""

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()

    Path(log_file).parent.mkdir(parents=True, exist_ok=True)

    formatter = logging.Formatter(
        "%(asctime)s\t%(levelname)s\t%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger
