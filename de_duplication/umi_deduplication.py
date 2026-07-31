"""Mapping-free UMI deduplication for short Ribo-seq reads.

FASTQ parsing is performed with pysam.FastxFile.
UMIs are detected automatically from the trailing A/C/G/T sequence in each
FASTQ read ID. Base qualities are stored as raw Phred integer bytes in
intermediate files and converted back to FASTQ ASCII only when representative
reads are written.
"""

from __future__ import annotations

import gzip
import re
import shutil
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterator, TextIO

import pysam


DNA_TO_BITS = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
}

BITS_TO_DNA = "ACGT"

BLOCK_SIZE = 8
NUM_BLOCKS = 5
MAX_SEQUENCE_LENGTH = BLOCK_SIZE * NUM_BLOCKS

BUCKET_PREFIX_LENGTH = 5
FLUSH_INTERVAL = 5_000_000

UMI_LENGTH_FILE = "umi_length.txt"
UMI_SUFFIX_PATTERN = re.compile(r"([ACGT]+)$")

# Binary record layout:
# read_id_length : uint16
# umi_code       : uint64
# sequence_length: uint8
# block_0..4     : five uint16 values
# quality_sum    : uint32
# quality_length : uint8
#
# Variable data following the header:
# read_id bytes
# raw Phred quality bytes
RECORD_HEADER = struct.Struct("<HQB5HIB")


@dataclass
class SequenceRecord:
    """Aggregated information for one exact sequence within one UMI."""

    sequence_length: int
    blocks: tuple[int, int, int, int, int]
    count: int
    max_quality_sum: int
    quality_scores: bytes
    read_id: str


@dataclass
class ProcessingStats:
    """Processing counters."""

    total_reads: int = 0
    written_records: int = 0
    skipped_short: int = 0
    skipped_long: int = 0
    skipped_non_acgt: int = 0
    exact_sequence_groups: int = 0
    collapsed_sequence_groups: int = 0
    output_reads: int = 0


def count_full_block_mismatches(
    block_a: int,
    block_b: int,
) -> int:
    """Count mismatches between two 8-nt two-bit encoded blocks."""

    xor_value = block_a ^ block_b

    return (
        (
            xor_value
            | (xor_value >> 1)
        )
        & 0x5555
    ).bit_count()


def count_sequence_mismatches(
    record_a: SequenceRecord,
    record_b: SequenceRecord,
    maximum_mismatches: int = 3,
) -> int:
    """Count sequence mismatches with early termination."""

    if record_a.sequence_length != record_b.sequence_length:
        return maximum_mismatches + 1

    mismatch_count = 0
    full_block_count = (
        record_a.sequence_length
        // BLOCK_SIZE
    )
    partial_block_length = (
        record_a.sequence_length
        % BLOCK_SIZE
    )

    for block_index in range(full_block_count):
        block_a = record_a.blocks[block_index]
        block_b = record_b.blocks[block_index]

        if block_a == block_b:
            continue

        mismatch_count += count_full_block_mismatches(
            block_a,
            block_b,
        )

        if mismatch_count > maximum_mismatches:
            return mismatch_count

    if partial_block_length > 0:
        block_a = record_a.blocks[full_block_count]
        block_b = record_b.blocks[full_block_count]

        if block_a != block_b:
            xor_value = block_a ^ block_b

            mismatch_bits = (
                xor_value
                | (xor_value >> 1)
            ) & 0x5555

            valid_mask = (
                1 << (
                    partial_block_length * 2
                )
            ) - 1

            mismatch_count += (
                mismatch_bits
                & valid_mask
            ).bit_count()

    return mismatch_count


def read_fastq(
    path: str | Path,
) -> Iterator[tuple[str, str, bytes]]:
    """Yield read ID, sequence, and raw integer Phred scores."""

    with pysam.FastxFile(str(path)) as handle:
        for record in handle:
            if record.quality is None:
                raise ValueError(
                    f"Quality is missing for read: {record.name}"
                )

            yield (
                record.name,
                record.sequence.upper(),
                bytes(record.get_quality_array(offset=33)),
            )


def extract_umi_from_read_id(
    read_id: str,
) -> tuple[str, str]:
    """Extract a trailing A/C/G/T UMI from a FASTQ read ID."""

    match = UMI_SUFFIX_PATTERN.search(read_id)

    if match is None:
        raise ValueError(
            f"No trailing UMI was found in read ID: {read_id}"
        )

    umi = match.group(1)

    output_read_id = read_id[:match.start()].rstrip("_:-")

    return output_read_id or read_id, umi


def encode_dna(sequence: str) -> int:
    """Encode an A/C/G/T sequence as a two-bit integer."""

    code = 0

    for base in sequence:
        code = (code << 2) | DNA_TO_BITS[base]

    return code


def decode_dna(code: int, length: int) -> str:
    """Decode a two-bit integer into an A/C/G/T sequence."""

    bases = ["A"] * length

    for index in range(length - 1, -1, -1):
        bases[index] = BITS_TO_DNA[code & 0b11]
        code >>= 2

    return "".join(bases)


def encode_sequence_blocks(
    sequence: str,
) -> tuple[int, int, int, int, int]:
    """Encode up to 40 nt as five 8-mer integer blocks."""

    if len(sequence) > MAX_SEQUENCE_LENGTH:
        raise ValueError(
            f"Sequence is longer than {MAX_SEQUENCE_LENGTH} nt: "
            f"{len(sequence)}"
        )

    blocks = [0] * NUM_BLOCKS

    for block_index in range(NUM_BLOCKS):
        start = block_index * BLOCK_SIZE
        block_sequence = sequence[start : start + BLOCK_SIZE]

        if block_sequence:
            blocks[block_index] = encode_dna(block_sequence)

    return (
        blocks[0],
        blocks[1],
        blocks[2],
        blocks[3],
        blocks[4],
    )


def decode_sequence_blocks(
    blocks: tuple[int, int, int, int, int],
    sequence_length: int,
) -> str:
    """Restore a sequence from five encoded blocks."""

    remaining = sequence_length
    sequence_parts: list[str] = []

    for block in blocks:
        block_length = min(BLOCK_SIZE, remaining)

        if block_length <= 0:
            break

        sequence_parts.append(
            decode_dna(
                block,
                block_length,
            )
        )
        remaining -= block_length

    return "".join(sequence_parts)


def quality_scores_to_fastq(
    quality_scores: bytes,
) -> str:
    """Convert raw integer Phred scores to a Phred+33 FASTQ string."""

    return bytes(
        score + 33
        for score in quality_scores
    ).decode("ascii")


def pack_record(
    read_id: str,
    umi_code: int,
    sequence_length: int,
    blocks: tuple[int, int, int, int, int],
    quality_sum: int,
    quality_scores: bytes,
) -> bytes:
    """Pack one intermediate record into binary format."""

    read_id_bytes = read_id.encode("utf-8")
    quality_length = len(quality_scores)

    if len(read_id_bytes) > 65535:
        raise ValueError(
            f"Read ID is too long: {read_id}"
        )

    if quality_length > 255:
        raise ValueError(
            f"Quality array is too long: {quality_length}"
        )

    header = RECORD_HEADER.pack(
        len(read_id_bytes),
        umi_code,
        sequence_length,
        *blocks,
        quality_sum,
        quality_length,
    )

    return (
        header
        + read_id_bytes
        + quality_scores
    )


def read_exact(
    handle: BinaryIO,
    size: int,
) -> bytes:
    """Read exactly the requested number of bytes."""

    data = handle.read(size)

    if len(data) != size:
        raise ValueError(
            "Truncated intermediate bucket file"
        )

    return data


def read_bucket_records(
    path: str | Path,
) -> Iterator[
    tuple[
        str,
        int,
        int,
        tuple[int, int, int, int, int],
        int,
        bytes,
    ]
]:
    """Read records from one binary bucket file."""

    with open(path, "rb") as handle:
        while True:
            header = handle.read(RECORD_HEADER.size)

            if not header:
                break

            if len(header) != RECORD_HEADER.size:
                raise ValueError(
                    f"Truncated record header in {path}"
                )

            values = RECORD_HEADER.unpack(header)

            read_id_length = values[0]
            umi_code = values[1]
            sequence_length = values[2]
            blocks = (
                values[3],
                values[4],
                values[5],
                values[6],
                values[7],
            )
            quality_sum = values[8]
            quality_length = values[9]

            read_id = read_exact(
                handle,
                read_id_length,
            ).decode("utf-8")

            quality_scores = read_exact(
                handle,
                quality_length,
            )

            yield (
                read_id,
                umi_code,
                sequence_length,
                blocks,
                quality_sum,
                quality_scores,
            )


def flush_all_buffers(
    buffers: dict[str, bytearray],
    output_dir: Path,
) -> None:
    """Append all non-empty bucket buffers to disk."""

    for filename, buffer in buffers.items():
        if not buffer:
            continue

        with open(
            output_dir / filename,
            "ab",
        ) as output_handle:
            output_handle.write(buffer)

        buffer.clear()


def split_fastq_into_buckets(
    input_fastq: str | Path,
    output_dir: str | Path,
) -> ProcessingStats:
    """Create UMI bucket files and return processing statistics."""

    stats = ProcessingStats()

    output_dir = Path(output_dir)
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for old_path in output_dir.glob("*.data"):
        old_path.unlink()

    umi_length_path = output_dir / UMI_LENGTH_FILE

    if umi_length_path.exists():
        umi_length_path.unlink()

    buffers: dict[str, bytearray] = {}
    next_flush = FLUSH_INTERVAL
    detected_umi_length: int | None = None

    print(
        "Start reading FASTQ",
        file=sys.stderr,
        flush=True,
    )

    for (
        read_id,
        sequence,
        quality_scores,
    ) in read_fastq(input_fastq):
        stats.total_reads += 1

        try:
            output_read_id, umi = extract_umi_from_read_id(
                read_id
            )
        except ValueError:
            stats.skipped_short += 1
            continue

        if detected_umi_length is None:
            detected_umi_length = len(umi)

            if detected_umi_length < BUCKET_PREFIX_LENGTH:
                raise ValueError(
                    "Detected UMI is shorter than the "
                    f"{BUCKET_PREFIX_LENGTH}-nt bucket prefix: {umi}"
                )

            print(
                f"Detected UMI length: "
                f"{detected_umi_length} nt",
                file=sys.stderr,
                flush=True,
            )

        elif len(umi) != detected_umi_length:
            raise ValueError(
                f"Inconsistent UMI length in {read_id}: "
                f"expected {detected_umi_length}, "
                f"found {len(umi)}"
            )

        sequence_length = len(sequence)

        if sequence_length > MAX_SEQUENCE_LENGTH:
            stats.skipped_long += 1
            continue

        try:
            umi_code = encode_dna(umi)
            blocks = encode_sequence_blocks(sequence)
        except KeyError:
            stats.skipped_non_acgt += 1
            continue

        filename = (
            f"{umi[:BUCKET_PREFIX_LENGTH]}.data"
        )

        packed_record = pack_record(
            read_id=output_read_id,
            umi_code=umi_code,
            sequence_length=sequence_length,
            blocks=blocks,
            quality_sum=sum(quality_scores),
            quality_scores=quality_scores,
        )

        buffer = buffers.get(filename)

        if buffer is None:
            buffer = bytearray()
            buffers[filename] = buffer

        buffer.extend(packed_record)
        stats.written_records += 1

        if stats.total_reads % 10000 == 0:
            print(stats.total_reads," reads processed")
        if stats.total_reads == next_flush:
            flush_all_buffers(
                buffers,
                output_dir,
            )

            print(
                f"[bucket] processed={stats.total_reads:,} "
                f"written={stats.written_records:,}",
                file=sys.stderr,
                flush=True,
            )

            next_flush += FLUSH_INTERVAL

    flush_all_buffers(
        buffers,
        output_dir,
    )

    if detected_umi_length is None:
        raise ValueError(
            "No valid UMI was detected"
        )

    umi_length_path.write_text(
        f"{detected_umi_length}\n",
        encoding="utf-8",
    )

    print(
        f"[bucket] finished "
        f"processed={stats.total_reads:,} "
        f"written={stats.written_records:,}",
        file=sys.stderr,
        flush=True,
    )

    return stats


def aggregate_bucket(
    path: str | Path,
    stats: ProcessingStats,
) -> dict[
    int,
    dict[
        tuple[int, int, int, int, int, int],
        SequenceRecord,
    ],
]:
    """Aggregate exact sequences within each UMI."""

    umi_dictionary: dict[
        int,
        dict[
            tuple[int, int, int, int, int, int],
            SequenceRecord,
        ],
    ] = {}

    for (
        read_id,
        umi_code,
        sequence_length,
        blocks,
        quality_sum,
        quality_scores,
    ) in read_bucket_records(path):
        sequence_key = (
            sequence_length,
            *blocks,
        )

        sequence_dictionary = umi_dictionary.get(
            umi_code
        )

        if sequence_dictionary is None:
            sequence_dictionary = {}
            umi_dictionary[umi_code] = sequence_dictionary

        existing = sequence_dictionary.get(
            sequence_key
        )

        if existing is None:
            sequence_dictionary[sequence_key] = SequenceRecord(
                sequence_length=sequence_length,
                blocks=blocks,
                count=1,
                max_quality_sum=quality_sum,
                quality_scores=quality_scores,
                read_id=read_id,
            )
            stats.exact_sequence_groups += 1
            continue

        existing.count += 1

        if quality_sum > existing.max_quality_sum:
            existing.max_quality_sum = quality_sum
            existing.quality_scores = quality_scores
            existing.read_id = read_id

    return umi_dictionary


def get_two_block_keys(
    record: SequenceRecord,
) -> tuple[
    tuple[int, int, int, int],
    tuple[int, int, int, int],
    tuple[int, int, int, int],
    tuple[int, int, int, int],
    tuple[int, int, int, int],
    tuple[int, int, int, int],
]:
    """Return all two-block keys from the first four blocks."""

    block_0, block_1, block_2, block_3 = record.blocks[:4]

    return (
        (0, 1, block_0, block_1),
        (0, 2, block_0, block_2),
        (0, 3, block_0, block_3),
        (1, 2, block_1, block_2),
        (1, 3, block_1, block_3),
        (2, 3, block_2, block_3),
    )


def collapse_one_umi(
    sequence_records: list[SequenceRecord],
) -> list[SequenceRecord]:
    """Collapse sequences using two-block keys and at most three mismatches."""

    sorted_records = sorted(
        sequence_records,
        key=lambda record: (
            -record.count,
            -record.max_quality_sum,
        ),
    )

    deduplicated: list[SequenceRecord] = []

    key_to_indices: dict[
        tuple[int, int, int, int],
        list[int],
    ] = {}

    for candidate in sorted_records:
        candidate_indices: set[int] = set()

        for key in get_two_block_keys(candidate):
            candidate_indices.update(
                key_to_indices.get(
                    key,
                    [],
                )
            )

        collapsed = False

        for representative_index in sorted(
            candidate_indices
        ):
            representative = deduplicated[
                representative_index
            ]

            mismatch_count = count_sequence_mismatches(
                representative,
                candidate,
                maximum_mismatches=3,
            )

            if mismatch_count <= 3:
                representative.count += candidate.count
                collapsed = True
                break

        if collapsed:
            continue

        representative_index = len(
            deduplicated
        )
        deduplicated.append(candidate)

        for key in get_two_block_keys(candidate):
            key_to_indices.setdefault(
                key,
                [],
            ).append(
                representative_index
            )

    return deduplicated


def write_fastq_record(
    handle: TextIO,
    record: SequenceRecord,
    umi_code: int,
    umi_length: int,
) -> None:
    """Write one representative sequence in FASTQ format."""

    sequence = decode_sequence_blocks(
        record.blocks,
        record.sequence_length,
    )

    umi = decode_dna(
        umi_code,
        umi_length,
    )

    if (
        len(record.quality_scores)
        != len(sequence)
    ):
        raise ValueError(
            "Quality length differs from sequence length "
            f"for {record.read_id}"
        )

    output_id = (
        f"{record.read_id}"
        f"|umi={umi}"
        f"|family_size={record.count}"
    )

    handle.write(
        f"@{output_id}\n"
        f"{sequence}\n"
        "+\n"
        f"{quality_scores_to_fastq(record.quality_scores)}\n"
    )


def process_buckets(
    bucket_dir: str | Path,
    output_fastq_gz: str | Path,
) -> ProcessingStats:
    """Collapse bucket records and return processing statistics."""

    stats = ProcessingStats()

    bucket_dir = Path(bucket_dir)

    umi_length_path = (
        bucket_dir / UMI_LENGTH_FILE
    )

    if not umi_length_path.exists():
        raise FileNotFoundError(
            f"UMI length file was not found: "
            f"{umi_length_path}"
        )

    umi_length = int(
        umi_length_path.read_text(
            encoding="utf-8"
        ).strip()
    )

    bucket_paths = sorted(
        bucket_dir.glob("*.data")
    )

    if not bucket_paths:
        raise FileNotFoundError(
            f"No bucket files were found in {bucket_dir}"
        )

    with gzip.open(
        output_fastq_gz,
        "wt",
        encoding="utf-8",
        compresslevel=6,
    ) as output_handle:
        for (
            bucket_index,
            bucket_path,
        ) in enumerate(
            bucket_paths,
            start=1,
        ):
            umi_dictionary = aggregate_bucket(
                bucket_path,
                stats,
            )

            for (
                umi_code,
                sequence_dictionary,
            ) in umi_dictionary.items():
                retained = collapse_one_umi(
                    list(
                        sequence_dictionary.values()
                    )
                )

                stats.collapsed_sequence_groups += (
                    len(sequence_dictionary)
                    - len(retained)
                )

                for record in retained:
                    write_fastq_record(
                        output_handle,
                        record,
                        umi_code,
                        umi_length,
                    )

                    stats.output_reads += 1

            print(
                f"[collapse] "
                f"bucket={bucket_index}/{len(bucket_paths)} "
                f"name={bucket_path.name} "
                f"output_reads={stats.output_reads:,}",
                file=sys.stderr,
                flush=True,
            )

    shutil.rmtree(
        bucket_dir
    )
    print(
        f"[cleanup] removed bucket directory: "
        f"{bucket_dir}",
        file=sys.stderr,
        flush=True,
    )
    return stats


def write_stats(
    path: str | Path,
    bucket_stats: ProcessingStats,
    collapse_stats: ProcessingStats,
) -> None:
    """Write bucket and collapse statistics as TSV."""

    values = {
        "total_reads": bucket_stats.total_reads,
        "written_records": bucket_stats.written_records,
        "skipped_short": bucket_stats.skipped_short,
        "skipped_long": bucket_stats.skipped_long,
        "skipped_non_acgt": bucket_stats.skipped_non_acgt,
        "exact_sequence_groups": collapse_stats.exact_sequence_groups,
        "collapsed_sequence_groups": (
            collapse_stats.collapsed_sequence_groups
        ),
        "output_reads": collapse_stats.output_reads,
    }

    with open(
        path,
        "wt",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "metric\tvalue\n"
        )

        for key, value in values.items():
            handle.write(
                f"{key}\t{value}\n"
            )