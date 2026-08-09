from __future__ import annotations

import gzip
import logging
import math
import re
from pathlib import Path
from typing import Iterator, TextIO


_CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
_CIGAR_OPERATION_CODES = {
    "M": 0,
    "I": 1,
    "D": 2,
    "N": 3,
    "S": 4,
    "H": 5,
    "P": 6,
    "=": 7,
    "X": 8,
}


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

    return value[1:].split(maxsplit=1)[0]


def iter_fastq_records(
    fastq_path: str | Path,
) -> Iterator[tuple[int, str, str, str]]:
    """Yield FASTQ index, read ID, sequence, and quality."""

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

            yield (
                record_index,
                normalize_fastq_read_id(header),
                sequence.rstrip("\r\n"),
                quality.rstrip("\r\n"),
            )
            record_index += 1


def parse_cigar(cigar: str) -> list[tuple[int, int]]:
    """Parse one compact CIGAR string into pysam CIGAR tuples."""

    if not cigar or cigar == "*":
        raise ValueError("Mapped assignment requires a non-empty CIGAR.")

    matches = _CIGAR_PATTERN.findall(cigar)

    if not matches or "".join(f"{length}{op}" for length, op in matches) != cigar:
        raise ValueError(f"Invalid CIGAR string: {cigar}")

    return [
        (_CIGAR_OPERATION_CODES[operation], int(length))
        for length, operation in matches
    ]


def cigar_query_length(
    cigar_tuples: list[tuple[int, int]],
) -> int:
    """Return query length consumed by pysam CIGAR tuples."""

    query_consuming_codes = {0, 1, 4, 7, 8}

    return sum(
        length
        for operation, length in cigar_tuples
        if operation in query_consuming_codes
    )


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a nucleotide sequence."""

    translation = str.maketrans(
        "ACGTRYMKBDHVacgtrymkbdhvNn",
        "TGCAYRKMVHDBtgcayrkmvhdbNn",
    )
    return sequence.translate(translation)[::-1]


def read_genome_fai(genome_fai: str | Path) -> list[tuple[str, int]]:
    """Read BAM reference names and lengths from a FASTA index."""

    references: list[tuple[str, int]] = []

    with Path(genome_fai).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 2:
                raise ValueError(f"Invalid FAI line {line_number}.")

            references.append((fields[0], int(fields[1])))

    if not references:
        raise RuntimeError("Genome FAI contains no references.")

    return references


def _mapq_from_posterior(value: float | None) -> int:
    """Convert assignment posterior to a capped informative MAPQ."""

    if value is None or not math.isfinite(value):
        return 0

    posterior = min(max(float(value), 0.0), 1.0)

    if posterior >= 1.0:
        return 60

    if posterior <= 0.0:
        return 0

    return min(60, max(0, round(-10.0 * math.log10(1.0 - posterior))))


def _optional_value(row: dict, *names: str):
    """Return the first available non-null row value."""

    for name in names:
        value = row.get(name)

        if value is not None:
            return value

    return None


def _load_partition_rows(
    parquet_path: str | Path,
) -> list[dict]:
    """Load one FASTQ-order assignment partition."""

    try:
        import pyarrow.parquet as pq
    except ImportError as error:
        raise RuntimeError("BAM writing requires pyarrow.") from error

    table = pq.read_table(parquet_path)
    required = {
        "read_id",
        "fastq_index",
        "chrom",
        "genome_start",
        "strand",
        "genome_cigar",
    }
    missing = sorted(required - set(table.column_names))

    if missing:
        raise ValueError(
            "Assignment partition is missing columns: " + ", ".join(missing)
        )

    rows = table.to_pylist()
    rows.sort(key=lambda row: int(row["fastq_index"]))
    return rows


def _make_aligned_segment(
    pysam_module,
    row: dict,
    sequence: str,
    quality: str,
    reference_id_by_name: dict[str, int],
):
    """Construct one pysam AlignedSegment from one selected assignment."""

    read_id = str(row["read_id"])
    chrom = str(row["chrom"])
    strand_value = row["strand"]
    is_reverse = (
        strand_value == "-"
        or strand_value == -1
        or strand_value == "-1"
    )

    if chrom not in reference_id_by_name:
        raise ValueError(f"Chromosome is absent from genome FAI: {chrom}")

    record = pysam_module.AlignedSegment()
    record.query_name = read_id
    record.flag = 16 if is_reverse else 0
    record.reference_id = reference_id_by_name[chrom]
    record.reference_start = int(row["genome_start"])

    cigar_tuples = parse_cigar(str(row["genome_cigar"]))
    consumed_query_length = cigar_query_length(cigar_tuples)

    if consumed_query_length != len(sequence):
        raise ValueError(
            "Genome CIGAR and FASTQ sequence lengths differ before BAM "
            f"writing: read_id={read_id}, "
            f"sequence_length={len(sequence)}, "
            f"cigar={row['genome_cigar']}, "
            f"cigar_query_length={consumed_query_length}"
        )

    record.cigartuples = cigar_tuples

    posterior_value = _optional_value(
        row,
        "posterior_probability",
        "assignment_probability",
        "final_probability",
    )
    posterior = float(posterior_value) if posterior_value is not None else None
    record.mapping_quality = _mapq_from_posterior(posterior)

    if is_reverse:
        record.query_sequence = reverse_complement(sequence)
        record.query_qualities = pysam_module.qualitystring_to_array(quality[::-1])
    else:
        record.query_sequence = sequence
        record.query_qualities = pysam_module.qualitystring_to_array(quality)

    selected_transcript = _optional_value(
        row,
        "selected_transcript_id",
        "transcript_id",
    )
    group_id = _optional_value(row, "group_id")
    assignment_type = _optional_value(row, "assignment_type")
    candidate_count = _optional_value(
        row,
        "original_candidate_count",
        "candidate_count",
        "genome_candidate_count",
    )

    if selected_transcript is not None:
        record.set_tag("XT", str(selected_transcript), value_type="Z")

    if group_id is not None:
        record.set_tag("XG", str(group_id), value_type="Z")

    if assignment_type is not None:
        record.set_tag("XA", str(assignment_type), value_type="Z")

    if posterior is not None:
        record.set_tag("XP", posterior, value_type="f")

    if candidate_count is not None:
        record.set_tag("NH", int(candidate_count), value_type="i")

    return record


def write_bam_from_fastq_order_partitions(
    fastq_path: str | Path,
    partition_dir: str | Path,
    partition_prefix: str,
    partition_count: int,
    genome_fai: str | Path,
    output_bam: str | Path,
    sort_threads: int = 4,
    sort_memory: str = "1G",
    sort_output: bool = True,
    index_output: bool = True,
    logger: logging.Logger | None = None,
) -> Path:
    """Stream FASTQ once and write an optional sorted/indexed BAM."""

    try:
        import pysam
    except ImportError as error:
        raise RuntimeError(
            "BAM output requires pysam. Install it with: python -m pip install pysam"
        ) from error

    if partition_count <= 0:
        raise ValueError("partition_count must be positive.")

    destination = Path(output_bam)
    destination.parent.mkdir(parents=True, exist_ok=True)
    unsorted_bam = (
        destination.with_suffix(".unsorted.bam")
        if sort_output
        else destination
    )

    references = read_genome_fai(genome_fai)
    reference_id_by_name = {
        name: reference_id
        for reference_id, (name, _) in enumerate(references)
    }
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [
            {"SN": name, "LN": length}
            for name, length in references
        ],
        "PG": [
            {
                "ID": "RiboTxAlign",
                "PN": "RiboTxAlign",
                "VN": "0.1.0",
                "CL": "FASTQ-order partition BAM writer",
            }
        ],
    }

    current_partition_id = -1
    current_rows: list[dict] = []
    current_row_index = 0
    written_records = 0
    missing_assignments = 0

    def load_partition(partition_id: int) -> None:
        nonlocal current_partition_id, current_rows, current_row_index

        parquet_path = (
            Path(partition_dir)
            / f"{partition_prefix}.part{partition_id:02d}.parquet"
        )

        if not parquet_path.exists():
            raise FileNotFoundError(parquet_path)

        current_partition_id = partition_id
        current_rows = _load_partition_rows(parquet_path)
        current_row_index = 0

        if logger is not None:
            logger.info(
                "Loaded assignment partition %02d: %d rows",
                partition_id,
                len(current_rows),
            )

    with pysam.AlignmentFile(unsorted_bam, "wb", header=header) as bam_handle:
        for fastq_index, read_id, sequence, quality in iter_fastq_records(
            fastq_path
        ):
            while (
                current_partition_id + 1 < partition_count
                and (
                    current_partition_id < 0
                    or not current_rows
                    or fastq_index
                    > int(current_rows[-1]["fastq_index"])
                )
            ):
                load_partition(current_partition_id + 1)

            while (
                current_row_index < len(current_rows)
                and int(current_rows[current_row_index]["fastq_index"])
                < fastq_index
            ):
                missing_assignments += 1
                current_row_index += 1

            if current_row_index >= len(current_rows):
                continue

            row = current_rows[current_row_index]

            if int(row["fastq_index"]) != fastq_index:
                continue

            if str(row["read_id"]) != read_id:
                raise RuntimeError(
                    "FASTQ-order mismatch at index "
                    f"{fastq_index}: FASTQ={read_id}, "
                    f"assignment={row['read_id']}"
                )

            record = _make_aligned_segment(
                pysam_module=pysam,
                row=row,
                sequence=sequence,
                quality=quality,
                reference_id_by_name=reference_id_by_name,
            )
            bam_handle.write(record)
            written_records += 1
            current_row_index += 1

            if logger is not None and written_records % 1_000_000 == 0:
                logger.info("Unsorted BAM records written: %d", written_records)

    if missing_assignments:
        raise RuntimeError(
            "Assignment rows were skipped before matching FASTQ records: "
            f"{missing_assignments}"
        )

    if written_records == 0:
        raise RuntimeError("No BAM records were written.")

    if not sort_output:
        if index_output:
            raise ValueError(
                "index_output=True requires sort_output=True."
            )

        if logger is not None:
            logger.info("Unsorted BAM written: %s", destination)
            logger.info("BAM records written: %d", written_records)

        return destination

    if logger is not None:
        logger.info("Sorting BAM by genome coordinate")

    pysam.sort(
        "-@",
        str(max(sort_threads, 1)),
        "-m",
        sort_memory,
        "-o",
        str(destination),
        str(unsorted_bam),
    )
    unsorted_bam.unlink(missing_ok=True)

    if index_output:
        if logger is not None:
            logger.info("Indexing sorted BAM")
        pysam.index("-@", str(max(sort_threads, 1)), str(destination))

    if logger is not None:
        logger.info("Sorted BAM written: %s", destination)
        if index_output:
            logger.info("BAM index written: %s.bai", destination)
        logger.info("BAM records written: %d", written_records)

    return destination
