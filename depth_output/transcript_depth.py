from __future__ import annotations

import csv
import logging
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pyarrow.parquet as pq

from db_construct.gtf_reader import parse_gtf_attributes
from depth_output.transcript_genome_coordinate import (
    TranscriptGenomeCoordinateConverter,
    read_genome_sizes,
)


@dataclass(slots=True)
class TranscriptDefinition:
    """Transcript definition used for depth-track output."""

    transcript_id: str
    track_name: str
    transcript_length: int
    cds_start: int | None
    cds_end: int | None
    gencode_transcript_ids: str
    stringtie_transcript_ids: str


@dataclass(slots=True)
class PsiteOffset:
    """Estimated P-site offset for one read length."""

    read_length: int
    offset: int
    peak_weight: float
    total_profile_weight: float
    used_fallback: bool


def create_logger(log_path: str | Path) -> logging.Logger:
    """Create a file and console logger."""

    path = Path(log_path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(f"RiboTxAlignDepth:{path}")
    logger.setLevel(logging.INFO)
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
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


def read_fasta_lengths(
    fasta_path: str | Path,
) -> dict[str, int]:
    """Read sequence lengths without loading complete sequences."""

    lengths: dict[str, int] = {}
    current_id: str | None = None
    current_length = 0

    with Path(fasta_path).open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_length

                current_id = line[1:].strip().split()[0]
                current_length = 0
            else:
                current_length += len(line.strip())

    if current_id is not None:
        lengths[current_id] = current_length

    return lengths


def load_gtf_features(
    reference_gtf: str | Path,
) -> tuple[
    dict[str, str],
    dict[str, list[tuple[int, int]]],
    dict[str, list[tuple[int, int]]],
]:
    """Load transcript strand, exon intervals, and CDS intervals."""

    strand_by_transcript: dict[str, str] = {}
    exons_by_transcript: dict[
        str,
        list[tuple[int, int]],
    ] = defaultdict(list)
    cds_by_transcript: dict[
        str,
        list[tuple[int, int]],
    ] = defaultdict(list)

    with Path(reference_gtf).open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            feature = fields[2]
            if feature not in {"transcript", "exon", "CDS"}:
                continue

            attributes = parse_gtf_attributes(fields[8])
            transcript_id = attributes.get("transcript_id")
            if not transcript_id:
                continue

            strand_by_transcript[transcript_id] = fields[6]

            if feature == "transcript":
                continue

            # Convert 1-based closed GTF coordinates to 0-based half-open.
            start = int(fields[3]) - 1
            end = int(fields[4])

            if feature == "exon":
                exons_by_transcript[transcript_id].append((start, end))
            else:
                cds_by_transcript[transcript_id].append((start, end))

    return (
        strand_by_transcript,
        dict(exons_by_transcript),
        dict(cds_by_transcript),
    )


def genomic_interval_to_transcript_interval(
    exons: list[tuple[int, int]],
    cds_intervals: list[tuple[int, int]],
    strand: str,
) -> tuple[int, int] | None:
    """Project genomic CDS intervals onto 5-prime-to-3-prime transcript coordinates."""

    if not exons or not cds_intervals:
        return None

    ordered_exons = sorted(
        exons,
        reverse=(strand == "-"),
    )

    transcript_cds_intervals: list[tuple[int, int]] = []
    transcript_offset = 0

    for exon_start, exon_end in ordered_exons:
        exon_length = exon_end - exon_start

        for cds_start, cds_end in cds_intervals:
            overlap_start = max(exon_start, cds_start)
            overlap_end = min(exon_end, cds_end)

            if overlap_start >= overlap_end:
                continue

            if strand == "+":
                relative_start = overlap_start - exon_start
                relative_end = overlap_end - exon_start
            else:
                relative_start = exon_end - overlap_end
                relative_end = exon_end - overlap_start

            transcript_cds_intervals.append(
                (
                    transcript_offset + relative_start,
                    transcript_offset + relative_end,
                )
            )

        transcript_offset += exon_length

    if not transcript_cds_intervals:
        return None

    return (
        min(start for start, _ in transcript_cds_intervals),
        max(end for _, end in transcript_cds_intervals),
    )



def _sanitize_track_component(value: str) -> str:
    """Return a BigWig-safe, whitespace-free name component."""

    return "_".join(value.strip().split())


def build_track_name(
    reference_id: str,
    metadata_row: dict[str, str],
) -> str:
    """Build a unique output name containing source transcript IDs."""

    parts = [_sanitize_track_component(reference_id)]

    gencode_ids = metadata_row.get("gencode_transcript_ids", "").strip()
    stringtie_ids = metadata_row.get("stringtie_transcript_ids", "").strip()

    if gencode_ids:
        parts.append(
            "GENCODE=" + _sanitize_track_component(gencode_ids)
        )

    if stringtie_ids:
        parts.append(
            "StringTie=" + _sanitize_track_component(stringtie_ids)
        )

    return "|".join(parts)

def _parse_optional_int(value: str | None) -> int | None:
    """Parse an optional integer field."""

    if value is None:
        return None

    value = value.strip()
    if not value:
        return None

    return int(float(value))


def load_transcript_definitions(
    reference_metadata_tsv: str | Path,
    reference_gtf: str | Path,
    transcript_fasta: str | Path,
    logger: logging.Logger | None = None,
) -> dict[str, TranscriptDefinition]:
    """Load transcript lengths and transcript-coordinate CDS intervals.

    Explicit ``cds_reference_start`` and ``cds_reference_end`` columns in
    reference metadata are preferred. GTF projection is used only as a
    fallback for transcripts without valid metadata CDS coordinates.
    """

    fasta_lengths = read_fasta_lengths(transcript_fasta)

    (
        strand_by_transcript,
        exons_by_transcript,
        cds_by_transcript,
    ) = load_gtf_features(reference_gtf)

    metadata_rows: dict[str, dict[str, str]] = {}
    metadata_fieldnames: set[str] = set()

    with Path(reference_metadata_tsv).open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        metadata_fieldnames = set(reader.fieldnames or ())

        if "reference_id" not in metadata_fieldnames:
            raise ValueError(
                "Reference metadata does not contain a reference_id column."
            )

        for row in reader:
            reference_id = row.get("reference_id", "").strip()
            if reference_id:
                metadata_rows[reference_id] = row

    has_explicit_cds_columns = {
        "cds_reference_start",
        "cds_reference_end",
    }.issubset(metadata_fieldnames)

    if logger is not None and not has_explicit_cds_columns:
        logger.warning(
            "Reference metadata does not contain both "
            "cds_reference_start and cds_reference_end; "
            "GTF projection will be used as fallback."
        )

    definitions: dict[str, TranscriptDefinition] = {}
    matched_metadata_count = 0
    metadata_cds_count = 0
    gtf_cds_count = 0
    invalid_cds_count = 0

    for transcript_id, transcript_length in fasta_lengths.items():
        row = metadata_rows.get(transcript_id, {})

        if row:
            matched_metadata_count += 1

        cds_start = None
        cds_end = None

        if has_explicit_cds_columns and row:
            candidate_start = _parse_optional_int(
                row.get("cds_reference_start")
            )
            candidate_end = _parse_optional_int(
                row.get("cds_reference_end")
            )

            if candidate_start is not None and candidate_end is not None:
                if (
                    0 <= candidate_start
                    < candidate_end
                    <= transcript_length
                ):
                    cds_start = candidate_start
                    cds_end = candidate_end
                    metadata_cds_count += 1
                else:
                    invalid_cds_count += 1

        # Retain compatibility with older metadata column names.
        if cds_start is None or cds_end is None:
            for start_column, end_column in (
                ("transcript_cds_start", "transcript_cds_end"),
            ):
                if start_column in row and end_column in row:
                    candidate_start = _parse_optional_int(
                        row.get(start_column)
                    )
                    candidate_end = _parse_optional_int(
                        row.get(end_column)
                    )

                    if (
                        candidate_start is not None
                        and candidate_end is not None
                        and 0 <= candidate_start
                        < candidate_end
                        <= transcript_length
                    ):
                        cds_start = candidate_start
                        cds_end = candidate_end
                        metadata_cds_count += 1
                        break

        if cds_start is None or cds_end is None:
            projected = genomic_interval_to_transcript_interval(
                exons=exons_by_transcript.get(transcript_id, []),
                cds_intervals=cds_by_transcript.get(transcript_id, []),
                strand=strand_by_transcript.get(transcript_id, "+"),
            )

            if projected is not None:
                projected_start, projected_end = projected

                if (
                    0 <= projected_start
                    < projected_end
                    <= transcript_length
                ):
                    cds_start = projected_start
                    cds_end = projected_end
                    gtf_cds_count += 1
                else:
                    invalid_cds_count += 1

        definitions[transcript_id] = TranscriptDefinition(
            transcript_id=transcript_id,
            track_name=build_track_name(transcript_id, row),
            transcript_length=transcript_length,
            cds_start=cds_start,
            cds_end=cds_end,
            gencode_transcript_ids=row.get(
                "gencode_transcript_ids",
                "",
            ),
            stringtie_transcript_ids=row.get(
                "stringtie_transcript_ids",
                "",
            ),
        )

    if logger is not None:
        logger.info(
            "Metadata rows matched to FASTA: %d",
            matched_metadata_count,
        )
        logger.info(
            "CDS loaded directly from metadata: %d",
            metadata_cds_count,
        )
        logger.info(
            "CDS recovered from GTF fallback: %d",
            gtf_cds_count,
        )
        logger.info(
            "Invalid CDS coordinates: %d",
            invalid_cds_count,
        )

    return definitions


def iter_read_assignments(
    reads_info_parquet: str | Path,
    batch_size: int = 100000,
):
    """Yield weighted transcript-alignment records."""

    parquet_file = pq.ParquetFile(reads_info_parquet)
    available_columns = set(parquet_file.schema.names)

    required_columns = {
        "transcript_id",
        "start",
        "end",
        "prior_probability",
        "transcript_probability",
    }

    missing_columns = required_columns - available_columns
    if missing_columns:
        raise ValueError(
            "Missing reads_info columns: "
            + ",".join(sorted(missing_columns))
        )

    columns = [
        "transcript_id",
        "start",
        "end",
        "prior_probability",
        "transcript_probability",
    ]

    has_read_length = "read_length" in available_columns
    if has_read_length:
        columns.append("read_length")

    for batch in parquet_file.iter_batches(
        batch_size=batch_size,
        columns=columns,
    ):
        data = batch.to_pydict()
        row_count = batch.num_rows

        for index in range(row_count):
            start = int(data["start"][index])
            end = int(data["end"][index])

            read_length = (
                int(data["read_length"][index])
                if has_read_length
                else end - start
            )

            weight = (
                float(data["prior_probability"][index])
                * float(data["transcript_probability"][index])
            )

            yield (
                data["transcript_id"][index],
                start,
                end,
                read_length,
                weight,
            )


def build_cds_start_profiles(
    reads_info_parquet: str | Path,
    transcript_definitions: dict[str, TranscriptDefinition],
    profile_minimum: int = -40,
    profile_maximum: int = 20,
    batch_size: int = 100000,
) -> dict[int, dict[int, float]]:
    """Build weighted 5-prime-end profiles around CDS starts."""

    profiles: dict[int, dict[int, float]] = defaultdict(
        lambda: defaultdict(float)
    )

    for (
        transcript_id,
        start,
        _end,
        read_length,
        weight,
    ) in iter_read_assignments(
        reads_info_parquet,
        batch_size=batch_size,
    ):
        definition = transcript_definitions.get(transcript_id)
        if definition is None or definition.cds_start is None:
            continue

        relative_start = start - definition.cds_start

        if profile_minimum <= relative_start <= profile_maximum:
            profiles[read_length][relative_start] += weight

    return {
        read_length: dict(profile)
        for read_length, profile in profiles.items()
    }


def estimate_psite_offsets(
    profiles: dict[int, dict[int, float]],
    minimum_offset: int = 8,
    maximum_offset: int = 20,
    fallback_offset: int = 12,
    minimum_profile_weight: float = 10.0,
) -> dict[int, PsiteOffset]:
    """Estimate one P-site offset from the CDS-start profile of each read length."""

    offsets: dict[int, PsiteOffset] = {}

    for read_length in sorted(profiles):
        profile = profiles[read_length]
        total_weight = sum(profile.values())

        candidate_offsets = range(
            minimum_offset,
            maximum_offset + 1,
        )

        best_offset = max(
            candidate_offsets,
            key=lambda offset: profile.get(-offset, 0.0),
        )
        peak_weight = profile.get(-best_offset, 0.0)

        used_fallback = (
            total_weight < minimum_profile_weight
            or peak_weight <= 0.0
        )

        if used_fallback:
            best_offset = fallback_offset
            peak_weight = profile.get(-best_offset, 0.0)

        offsets[read_length] = PsiteOffset(
            read_length=read_length,
            offset=best_offset,
            peak_weight=peak_weight,
            total_profile_weight=total_weight,
            used_fallback=used_fallback,
        )

    return offsets


def write_cds_start_profile(
    profiles: dict[int, dict[int, float]],
    output_path: str | Path,
) -> None:
    """Write read-length-specific CDS-start profiles."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["read_length", "relative_start", "weighted_count"])

        for read_length in sorted(profiles):
            for relative_start in sorted(profiles[read_length]):
                writer.writerow(
                    [
                        read_length,
                        relative_start,
                        f"{profiles[read_length][relative_start]:.12g}",
                    ]
                )


def write_psite_offsets(
    offsets: dict[int, PsiteOffset],
    output_path: str | Path,
) -> None:
    """Write estimated P-site offsets."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "read_length",
                "psite_offset",
                "peak_weight",
                "total_profile_weight",
                "used_fallback",
            ]
        )

        for read_length in sorted(offsets):
            offset = offsets[read_length]
            writer.writerow(
                [
                    offset.read_length,
                    offset.offset,
                    f"{offset.peak_weight:.12g}",
                    f"{offset.total_profile_weight:.12g}",
                    int(offset.used_fallback),
                ]
            )



def write_cds_start_stats(
    profiles: dict[int, dict[int, float]],
    offsets: dict[int, PsiteOffset],
    output_path: str | Path,
) -> None:
    """Write summary statistics for each read-length CDS-start profile."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "read_length",
                "total_weight",
                "nonzero_position_count",
                "profile_peak_position",
                "profile_peak_weight",
                "weighted_mean_position",
                "weighted_sd_position",
                "psite_offset",
                "offset_peak_weight",
                "used_fallback",
            ]
        )

        for read_length in sorted(profiles):
            profile = profiles[read_length]
            positive_items = [
                (position, weight)
                for position, weight in profile.items()
                if weight > 0.0
            ]
            total_weight = sum(weight for _, weight in positive_items)

            if positive_items and total_weight > 0.0:
                peak_position, peak_weight = max(
                    positive_items,
                    key=lambda item: item[1],
                )
                weighted_mean = sum(
                    position * weight
                    for position, weight in positive_items
                ) / total_weight
                weighted_variance = sum(
                    weight * ((position - weighted_mean) ** 2)
                    for position, weight in positive_items
                ) / total_weight
                weighted_sd = weighted_variance ** 0.5
            else:
                peak_position = ""
                peak_weight = 0.0
                weighted_mean = 0.0
                weighted_sd = 0.0

            offset = offsets.get(read_length)

            writer.writerow(
                [
                    read_length,
                    f"{total_weight:.12g}",
                    len(positive_items),
                    peak_position,
                    f"{peak_weight:.12g}",
                    f"{weighted_mean:.12g}",
                    f"{weighted_sd:.12g}",
                    offset.offset if offset is not None else "",
                    (
                        f"{offset.peak_weight:.12g}"
                        if offset is not None
                        else ""
                    ),
                    (
                        int(offset.used_fallback)
                        if offset is not None
                        else ""
                    ),
                ]
            )


def plot_cds_start_profiles(
    profiles: dict[int, dict[int, float]],
    output_path: str | Path,
) -> None:
    """Plot weighted CDS-start profiles by read length."""

    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError(
            "CDS-start plotting requires matplotlib."
        ) from error

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(10, 6))

    for read_length in sorted(profiles):
        profile = profiles[read_length]
        positions = sorted(profile)
        values = [profile[position] for position in positions]

        axis.plot(
            positions,
            values,
            linewidth=1.0,
            label=str(read_length),
        )

    axis.axvline(0, linewidth=0.8)
    axis.set_xlabel("Read 5-prime start relative to CDS start")
    axis.set_ylabel("Weighted read count")
    axis.set_title("CDS-start profiles by read length")
    axis.legend(
        title="Read length",
        fontsize="small",
        ncol=2,
    )
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)


def plot_psite_offsets(
    offsets: dict[int, PsiteOffset],
    output_path: str | Path,
) -> None:
    """Plot inferred P-site offsets by read length."""

    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError(
            "P-site offset plotting requires matplotlib."
        ) from error

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    read_lengths = sorted(offsets)
    offset_values = [
        offsets[read_length].offset
        for read_length in read_lengths
    ]

    figure, axis = plt.subplots(figsize=(8, 5))
    axis.bar(
        [str(read_length) for read_length in read_lengths],
        offset_values,
    )
    axis.set_xlabel("Read length")
    axis.set_ylabel("P-site offset from 5-prime end")
    axis.set_title("Estimated P-site offsets")
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)

def build_depth_events(
    reads_info_parquet: str | Path,
    transcript_definitions: dict[str, TranscriptDefinition],
    offsets: dict[int, PsiteOffset],
    fallback_offset: int = 12,
    batch_size: int = 100000,
) -> tuple[
    dict[str, dict[int, float]],
    dict[str, dict[int, float]],
]:
    """Build sparse difference arrays for coverage and sparse P-site depths."""

    coverage_events: dict[
        str,
        dict[int, float],
    ] = defaultdict(lambda: defaultdict(float))

    psite_depth: dict[
        str,
        dict[int, float],
    ] = defaultdict(lambda: defaultdict(float))

    for (
        transcript_id,
        start,
        end,
        read_length,
        weight,
    ) in iter_read_assignments(
        reads_info_parquet,
        batch_size=batch_size,
    ):
        definition = transcript_definitions.get(transcript_id)
        if definition is None:
            continue

        clipped_start = max(start, 0)
        clipped_end = min(end, definition.transcript_length)

        if clipped_start >= clipped_end or weight <= 0.0:
            continue

        coverage_events[transcript_id][clipped_start] += weight
        coverage_events[transcript_id][clipped_end] -= weight

        psite_offset = offsets.get(read_length)
        offset_value = (
            psite_offset.offset
            if psite_offset is not None
            else fallback_offset
        )

        psite_position = start + offset_value

        if 0 <= psite_position < definition.transcript_length:
            psite_depth[transcript_id][psite_position] += weight

    return (
        {key: dict(value) for key, value in coverage_events.items()},
        {key: dict(value) for key, value in psite_depth.items()},
    )



def build_psite_metagene(
    psite_depth: dict[str, dict[int, float]],
    transcript_definitions: dict[str, TranscriptDefinition],
    anchor: str,
    minimum_position: int = -50,
    maximum_position: int = 100,
) -> dict[int, dict[str, float | int]]:
    """Build raw and per-transcript-normalized P-site metagene profiles."""

    if anchor not in {"start", "stop"}:
        raise ValueError("anchor must be 'start' or 'stop'.")

    raw_sum: dict[int, float] = defaultdict(float)
    normalized_sum: dict[int, float] = defaultdict(float)
    contributing_transcripts: dict[int, int] = defaultdict(int)
    normalized_transcript_count = 0

    for transcript_id, positions in psite_depth.items():
        definition = transcript_definitions.get(transcript_id)
        if definition is None:
            continue

        anchor_position = (
            definition.cds_start
            if anchor == "start"
            else definition.cds_end
        )
        if anchor_position is None:
            continue

        # cds_end is stored as a 0-based half-open boundary. The stop-codon
        # anchor is the first base after the CDS, which keeps relative
        # positions consistent with the stored annotation boundary.
        transcript_profile: dict[int, float] = {}

        for relative_position in range(
            minimum_position,
            maximum_position + 1,
        ):
            absolute_position = anchor_position + relative_position
            value = positions.get(absolute_position, 0.0)
            transcript_profile[relative_position] = value
            raw_sum[relative_position] += value

            if 0 <= absolute_position < definition.transcript_length:
                contributing_transcripts[relative_position] += 1

        window_total = sum(transcript_profile.values())
        if window_total <= 0.0:
            continue

        normalized_transcript_count += 1
        for relative_position, value in transcript_profile.items():
            normalized_sum[relative_position] += value / window_total

    result: dict[int, dict[str, float | int]] = {}

    for relative_position in range(
        minimum_position,
        maximum_position + 1,
    ):
        normalized_mean = (
            normalized_sum[relative_position]
            / normalized_transcript_count
            if normalized_transcript_count > 0
            else 0.0
        )

        result[relative_position] = {
            "raw_weighted_sum": raw_sum[relative_position],
            "normalized_mean": normalized_mean,
            "available_transcript_count": contributing_transcripts[
                relative_position
            ],
            "normalized_transcript_count": normalized_transcript_count,
        }

    return result


def write_psite_metagene(
    profile: dict[int, dict[str, float | int]],
    output_path: str | Path,
) -> None:
    """Write raw and normalized P-site metagene values."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "relative_position",
                "raw_weighted_sum",
                "normalized_mean",
                "available_transcript_count",
                "normalized_transcript_count",
            ]
        )

        for relative_position in sorted(profile):
            values = profile[relative_position]
            writer.writerow(
                [
                    relative_position,
                    f"{float(values['raw_weighted_sum']):.12g}",
                    f"{float(values['normalized_mean']):.12g}",
                    int(values["available_transcript_count"]),
                    int(values["normalized_transcript_count"]),
                ]
            )


def plot_psite_metagene(
    profile: dict[int, dict[str, float | int]],
    output_path: str | Path,
    title: str,
) -> None:
    """Plot normalized metagene with raw weighted signal for reference."""

    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError(
            "Metagene plotting requires matplotlib."
        ) from error

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    positions = sorted(profile)
    normalized_values = [
        float(profile[position]["normalized_mean"])
        for position in positions
    ]
    raw_values = [
        float(profile[position]["raw_weighted_sum"])
        for position in positions
    ]

    figure, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    axes[0].plot(positions, normalized_values, linewidth=1.2)
    axes[0].axvline(
        0,
        color="black",
        linestyle="--",
        linewidth=0.8,
    )
    axes[0].set_ylabel("Normalized mean")
    axes[0].set_title(title)

    axes[1].plot(positions, raw_values, linewidth=1.0)
    axes[1].axvline(
        0,
        color="black",
        linestyle="--",
        linewidth=0.8,
    )
    axes[1].set_xlabel("Position relative to CDS boundary")
    axes[1].set_ylabel("Raw weighted sum")

    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)

def write_depth_bigwig(
    coverage_events: dict[str, dict[int, float]],
    psite_depth: dict[str, dict[int, float]],
    transcript_order: list[str],
    transcript_definitions: dict[str, TranscriptDefinition],
    coverage_output_path: str | Path,
    psite_output_path: str | Path,
) -> None:
    """Write transcript coverage and P-site depth directly to BigWig."""

    try:
        import pyBigWig  # type: ignore
    except ImportError as error:
        raise RuntimeError(
            "Direct BigWig output requires pyBigWig. "
            "Install it with: python -m pip install pyBigWig"
        ) from error

    coverage_output = Path(coverage_output_path)
    psite_output = Path(psite_output_path)
    coverage_output.parent.mkdir(parents=True, exist_ok=True)
    psite_output.parent.mkdir(parents=True, exist_ok=True)

    header = [
        (
            transcript_definitions[transcript_id].track_name,
            transcript_definitions[transcript_id].transcript_length,
        )
        for transcript_id in transcript_order
    ]

    coverage_bigwig = pyBigWig.open(str(coverage_output), "w")
    psite_bigwig = pyBigWig.open(str(psite_output), "w")

    try:
        coverage_bigwig.addHeader(header)
        psite_bigwig.addHeader(header)

        for transcript_id in transcript_order:
            definition = transcript_definitions[transcript_id]
            track_name = definition.track_name

            events = coverage_events.get(transcript_id, {})
            if events:
                starts: list[int] = []
                ends: list[int] = []
                values: list[float] = []
                current_depth = 0.0
                previous_position: int | None = None

                for position in sorted(events):
                    if (
                        previous_position is not None
                        and position > previous_position
                        and current_depth != 0.0
                    ):
                        starts.append(previous_position)
                        ends.append(position)
                        values.append(float(current_depth))

                    current_depth += events[position]
                    previous_position = position

                if starts:
                    coverage_bigwig.addEntries(
                        [track_name] * len(starts),
                        starts,
                        ends=ends,
                        values=values,
                    )

            positions = psite_depth.get(transcript_id, {})
            if positions:
                sorted_positions = [
                    position
                    for position in sorted(positions)
                    if positions[position] != 0.0
                ]

                if sorted_positions:
                    psite_bigwig.addEntries(
                        [track_name] * len(sorted_positions),
                        sorted_positions,
                        ends=[
                            position + 1
                            for position in sorted_positions
                        ],
                        values=[
                            float(positions[position])
                            for position in sorted_positions
                        ],
                    )

    finally:
        coverage_bigwig.close()
        psite_bigwig.close()




def iter_coverage_intervals(
    events: dict[int, float],
):
    """Yield non-zero half-open intervals from one sparse difference map."""

    current_depth = 0.0
    previous_position: int | None = None

    for position in sorted(events):
        if (
            previous_position is not None
            and position > previous_position
            and current_depth != 0.0
        ):
            yield (
                previous_position,
                position,
                float(current_depth),
            )

        current_depth += events[position]
        previous_position = position


def project_depth_to_genome(
    coverage_events: dict[str, dict[int, float]],
    psite_depth: dict[str, dict[int, float]],
    converter: TranscriptGenomeCoordinateConverter,
) -> tuple[
    dict[str, dict[str, dict[int, float]]],
    dict[str, dict[str, dict[int, float]]],
]:
    """Project transcript coverage and P-site depth to genomic coordinates."""

    genome_coverage_events: dict[
        str,
        dict[str, dict[int, float]],
    ] = {
        "+": defaultdict(lambda: defaultdict(float)),
        "-": defaultdict(lambda: defaultdict(float)),
    }

    genome_psite_depth: dict[
        str,
        dict[str, dict[int, float]],
    ] = {
        "+": defaultdict(lambda: defaultdict(float)),
        "-": defaultdict(lambda: defaultdict(float)),
    }

    for transcript_id, events in coverage_events.items():
        for reference_start, reference_end, value in iter_coverage_intervals(
            events
        ):
            projected_intervals = converter.project_interval(
                reference_id=transcript_id,
                reference_start=reference_start,
                reference_end=reference_end,
            )

            for chrom, genome_start, genome_end, strand in projected_intervals:
                genome_coverage_events[strand][chrom][genome_start] += value
                genome_coverage_events[strand][chrom][genome_end] -= value

    for transcript_id, positions in psite_depth.items():
        for reference_position, value in positions.items():
            if value == 0.0:
                continue

            projected = converter.project_position(
                reference_id=transcript_id,
                reference_position=reference_position,
            )

            if projected is None:
                continue

            chrom, genome_position, strand = projected
            genome_psite_depth[strand][chrom][genome_position] += value

    return (
        {
            strand: {
                chrom: dict(events)
                for chrom, events in chrom_map.items()
            }
            for strand, chrom_map in genome_coverage_events.items()
        },
        {
            strand: {
                chrom: dict(positions)
                for chrom, positions in chrom_map.items()
            }
            for strand, chrom_map in genome_psite_depth.items()
        },
    )


def write_genome_bigwig(
    genome_coverage_events: dict[
        str,
        dict[str, dict[int, float]],
    ],
    genome_psite_depth: dict[
        str,
        dict[str, dict[int, float]],
    ],
    genome_sizes: list[tuple[str, int]],
    coverage_output_paths: dict[str, str | Path],
    psite_output_paths: dict[str, str | Path],
) -> None:
    """Write strand-specific genomic coverage and P-site BigWig files."""

    try:
        import pyBigWig  # type: ignore
    except ImportError as error:
        raise RuntimeError(
            "Genome BigWig output requires pyBigWig. "
            "Install it with: python -m pip install pyBigWig"
        ) from error

    genome_size_map = dict(genome_sizes)

    for strand in ("+", "-"):
        coverage_path = Path(coverage_output_paths[strand])
        psite_path = Path(psite_output_paths[strand])
        coverage_path.parent.mkdir(parents=True, exist_ok=True)
        psite_path.parent.mkdir(parents=True, exist_ok=True)

        coverage_bigwig = pyBigWig.open(str(coverage_path), "w")
        psite_bigwig = pyBigWig.open(str(psite_path), "w")

        try:
            coverage_bigwig.addHeader(genome_sizes)
            psite_bigwig.addHeader(genome_sizes)

            for chrom, chrom_length in genome_sizes:
                events = genome_coverage_events[strand].get(chrom, {})

                starts: list[int] = []
                ends: list[int] = []
                values: list[float] = []

                for start, end, value in iter_coverage_intervals(events):
                    clipped_start = max(start, 0)
                    clipped_end = min(end, chrom_length)

                    if clipped_start >= clipped_end or value == 0.0:
                        continue

                    starts.append(clipped_start)
                    ends.append(clipped_end)
                    values.append(float(value))

                if starts:
                    coverage_bigwig.addEntries(
                        [chrom] * len(starts),
                        starts,
                        ends=ends,
                        values=values,
                    )

                positions = genome_psite_depth[strand].get(chrom, {})
                sorted_positions = [
                    position
                    for position in sorted(positions)
                    if (
                        0 <= position < chrom_length
                        and positions[position] != 0.0
                    )
                ]

                if sorted_positions:
                    psite_bigwig.addEntries(
                        [chrom] * len(sorted_positions),
                        sorted_positions,
                        ends=[
                            position + 1
                            for position in sorted_positions
                        ],
                        values=[
                            float(positions[position])
                            for position in sorted_positions
                        ],
                    )

        finally:
            coverage_bigwig.close()
            psite_bigwig.close()



def write_chrom_sizes(
    transcript_definitions: dict[str, TranscriptDefinition],
    transcript_order: list[str],
    output_path: str | Path,
) -> None:
    """Write only transcripts having non-zero weighted coverage."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as handle:
        for transcript_id in transcript_order:
            definition = transcript_definitions[transcript_id]
            handle.write(
                f"{definition.track_name}\t"
                f"{definition.transcript_length}\n"
            )


def write_transcript_name_map(
    transcript_definitions: dict[str, TranscriptDefinition],
    transcript_order: list[str],
    output_path: str | Path,
) -> None:
    """Write transcript IDs and appended source metadata."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(
            [
                "reference_id",
                "track_name",
                "transcript_length",
                "gencode_transcript_ids",
                "stringtie_transcript_ids",
            ]
        )

        for transcript_id in transcript_order:
            definition = transcript_definitions[transcript_id]
            writer.writerow(
                [
                    definition.transcript_id,
                    definition.track_name,
                    definition.transcript_length,
                    definition.gencode_transcript_ids,
                    definition.stringtie_transcript_ids,
                ]
            )


def _format_gtf_attributes(
    attributes: dict[str, str],
) -> str:
    """Format GTF attributes deterministically."""

    preferred_order = [
        "gene_id",
        "transcript_id",
        "gene_name",
        "gencode_transcript_ids",
        "stringtie_transcript_ids",
        "source_transcript_ids",
    ]

    ordered_keys = [
        key
        for key in preferred_order
        if key in attributes and attributes[key] != ""
    ]
    ordered_keys.extend(
        sorted(
            key
            for key, value in attributes.items()
            if key not in preferred_order and value != ""
        )
    )

    return " ".join(
        f'{key} "{str(attributes[key]).replace(chr(34), chr(92) + chr(34))}";'
        for key in ordered_keys
    )


def load_reference_gtf_attributes(
    reference_gtf: str | Path,
) -> dict[str, dict[str, str]]:
    """Load transcript attributes from the reference-definition GTF."""

    attributes_by_transcript: dict[str, dict[str, str]] = {}

    with Path(reference_gtf).open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            attributes = parse_gtf_attributes(fields[8])
            transcript_id = attributes.get("transcript_id")
            if not transcript_id:
                continue

            current = attributes_by_transcript.setdefault(
                transcript_id,
                {},
            )

            for key, value in attributes.items():
                if value and key not in current:
                    current[key] = value

    return attributes_by_transcript


def write_transcript_annotation_gtf(
    transcript_definitions: dict[str, TranscriptDefinition],
    transcript_order: list[str],
    reference_gtf: str | Path,
    output_path: str | Path,
) -> None:
    """
    Write transcript-coordinate GTF for transcripts with weighted depth.

    The sequence name remains the unchanged reference transcript ID so the
    annotation can be loaded together with transcript-coordinate BigWig files.
    Original GTF attributes are retained when available. Source transcript IDs
    from metadata are appended as attributes.
    """

    source_attributes = load_reference_gtf_attributes(reference_gtf)

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as handle:
        for transcript_id in transcript_order:
            definition = transcript_definitions[transcript_id]
            attributes = dict(
                source_attributes.get(transcript_id, {})
            )

            attributes["transcript_id"] = transcript_id
            attributes.setdefault("gene_id", transcript_id)

            if definition.gencode_transcript_ids:
                attributes["gencode_transcript_ids"] = (
                    definition.gencode_transcript_ids
                )

            if definition.stringtie_transcript_ids:
                attributes["stringtie_transcript_ids"] = (
                    definition.stringtie_transcript_ids
                )

            attribute_text = _format_gtf_attributes(attributes)
            transcript_end = definition.transcript_length

            for feature, start, end, frame in (
                ("transcript", 1, transcript_end, "."),
                ("exon", 1, transcript_end, "."),
            ):
                handle.write(
                    "\t".join(
                        [
                            definition.track_name,
                            "RiboTxAlign",
                            feature,
                            str(start),
                            str(end),
                            ".",
                            "+",
                            frame,
                            attribute_text,
                        ]
                    )
                    + "\n"
                )

            if (
                definition.cds_start is not None
                and definition.cds_end is not None
                and definition.cds_start < definition.cds_end
            ):
                handle.write(
                    "\t".join(
                        [
                            definition.track_name,
                            "RiboTxAlign",
                            "CDS",
                            str(definition.cds_start + 1),
                            str(definition.cds_end),
                            ".",
                            "+",
                            "0",
                            attribute_text,
                        ]
                    )
                    + "\n"
                )


def create_transcript_depth_tracks(
    reads_info_parquet: str | Path,
    reference_metadata_tsv: str | Path,
    reference_gtf: str | Path,
    transcript_fasta: str | Path,
    output_prefix: str | Path,
    coordinate_table: str | Path | None = None,
    genome_fai: str | Path | None = None,
    profile_minimum: int = -40,
    profile_maximum: int = 20,
    minimum_offset: int = 8,
    maximum_offset: int = 20,
    fallback_offset: int = 12,
    minimum_profile_weight: float = 10.0,
    metagene_minimum: int = -50,
    metagene_maximum: int = 100,
    batch_size: int = 100000,
    create_bigwig: bool = True,
    create_genome_bigwig: bool = True,
) -> dict[str, Path]:
    """Run CDS-start profiling, P-site inference, and depth-track generation."""

    prefix = Path(output_prefix).resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)

    paths = {
        "profile": Path(f"{prefix}.cds_start_profile.tsv"),
        "offsets": Path(f"{prefix}.psite_offsets.tsv"),
        "profile_stats": Path(f"{prefix}.cds_start_stats.tsv"),
        "profile_plot": Path(f"{prefix}.cds_start_profile.png"),
        "offset_plot": Path(f"{prefix}.psite_offsets.png"),
        "start_metagene": Path(
            f"{prefix}.cds_start_psite_metagene.tsv"
        ),
        "start_metagene_plot": Path(
            f"{prefix}.cds_start_psite_metagene.png"
        ),
        "stop_metagene": Path(
            f"{prefix}.cds_stop_psite_metagene.tsv"
        ),
        "stop_metagene_plot": Path(
            f"{prefix}.cds_stop_psite_metagene.png"
        ),
        "coverage_bigwig": Path(f"{prefix}.transcript_depth.bw"),
        "psite_bigwig": Path(f"{prefix}.psite_depth.bw"),
        "genome_coverage_plus_bigwig": Path(
            f"{prefix}.genome_transcript_depth.plus.bw"
        ),
        "genome_coverage_minus_bigwig": Path(
            f"{prefix}.genome_transcript_depth.minus.bw"
        ),
        "genome_psite_plus_bigwig": Path(
            f"{prefix}.genome_psite_depth.plus.bw"
        ),
        "genome_psite_minus_bigwig": Path(
            f"{prefix}.genome_psite_depth.minus.bw"
        ),
        "chrom_sizes": Path(f"{prefix}.transcript_sizes.tsv"),
        "name_map": Path(f"{prefix}.transcript_name_map.tsv"),
        "annotation_gtf": Path(f"{prefix}.depth_transcripts.gtf"),
        "log": Path(f"{prefix}.depth.log"),
    }

    logger = create_logger(paths["log"])
    logger.info("Transcript depth workflow started")

    definitions = load_transcript_definitions(
        reference_metadata_tsv=reference_metadata_tsv,
        reference_gtf=reference_gtf,
        transcript_fasta=transcript_fasta,
        logger=logger,
    )

    cds_count = sum(
        definition.cds_start is not None
        for definition in definitions.values()
    )

    logger.info("Loaded transcripts: %d", len(definitions))
    logger.info("Transcripts with CDS: %d", cds_count)

    if cds_count == 0:
        raise RuntimeError(
            "No transcript-coordinate CDS intervals were loaded. "
            "Confirm that reference_id values match the FASTA headers and "
            "that cds_reference_start/cds_reference_end are populated."
        )

    profiles = build_cds_start_profiles(
        reads_info_parquet=reads_info_parquet,
        transcript_definitions=definitions,
        profile_minimum=profile_minimum,
        profile_maximum=profile_maximum,
        batch_size=batch_size,
    )

    offsets = estimate_psite_offsets(
        profiles=profiles,
        minimum_offset=minimum_offset,
        maximum_offset=maximum_offset,
        fallback_offset=fallback_offset,
        minimum_profile_weight=minimum_profile_weight,
    )

    write_cds_start_profile(profiles, paths["profile"])
    write_psite_offsets(offsets, paths["offsets"])
    write_cds_start_stats(
        profiles,
        offsets,
        paths["profile_stats"],
    )
    plot_cds_start_profiles(
        profiles,
        paths["profile_plot"],
    )
    plot_psite_offsets(
        offsets,
        paths["offset_plot"],
    )

    for read_length in sorted(offsets):
        result = offsets[read_length]
        logger.info(
            "ReadLength=%d PsiteOffset=%d PeakWeight=%.6f "
            "ProfileWeight=%.6f Fallback=%s",
            read_length,
            result.offset,
            result.peak_weight,
            result.total_profile_weight,
            result.used_fallback,
        )

    coverage_events, psite_depth = build_depth_events(
        reads_info_parquet=reads_info_parquet,
        transcript_definitions=definitions,
        offsets=offsets,
        fallback_offset=fallback_offset,
        batch_size=batch_size,
    )

    start_metagene = build_psite_metagene(
        psite_depth=psite_depth,
        transcript_definitions=definitions,
        anchor="start",
        minimum_position=metagene_minimum,
        maximum_position=metagene_maximum,
    )
    stop_metagene = build_psite_metagene(
        psite_depth=psite_depth,
        transcript_definitions=definitions,
        anchor="stop",
        minimum_position=metagene_minimum,
        maximum_position=metagene_maximum,
    )

    write_psite_metagene(
        start_metagene,
        paths["start_metagene"],
    )
    plot_psite_metagene(
        start_metagene,
        paths["start_metagene_plot"],
        "CDS-start P-site metagene",
    )
    write_psite_metagene(
        stop_metagene,
        paths["stop_metagene"],
    )
    plot_psite_metagene(
        stop_metagene,
        paths["stop_metagene_plot"],
        "CDS-stop P-site metagene",
    )

    transcript_order = [
        transcript_id
        for transcript_id in definitions
        if coverage_events.get(transcript_id)
    ]

    logger.info(
        "Transcripts with weighted depth: %d",
        len(transcript_order),
    )
    logger.info(
        "Transcripts removed for zero depth: %d",
        len(definitions) - len(transcript_order),
    )

    if create_bigwig:
        write_depth_bigwig(
            coverage_events=coverage_events,
            psite_depth=psite_depth,
            transcript_order=transcript_order,
            transcript_definitions=definitions,
            coverage_output_path=paths["coverage_bigwig"],
            psite_output_path=paths["psite_bigwig"],
        )

    if create_genome_bigwig:
        if coordinate_table is None:
            raise ValueError(
                "coordinate_table is required when "
                "create_genome_bigwig=True."
            )

        if genome_fai is None:
            raise ValueError(
                "genome_fai is required when create_genome_bigwig=True."
            )

        logger.info("Loading transcript-to-genome coordinate table")

        converter = TranscriptGenomeCoordinateConverter.from_coordinate_table(
            coordinate_table
        )
        genome_sizes = read_genome_sizes(genome_fai)

        logger.info("Projecting transcript depth to genome coordinates")

        (
            genome_coverage_events,
            genome_psite_depth,
        ) = project_depth_to_genome(
            coverage_events=coverage_events,
            psite_depth=psite_depth,
            converter=converter,
        )

        write_genome_bigwig(
            genome_coverage_events=genome_coverage_events,
            genome_psite_depth=genome_psite_depth,
            genome_sizes=genome_sizes,
            coverage_output_paths={
                "+": paths["genome_coverage_plus_bigwig"],
                "-": paths["genome_coverage_minus_bigwig"],
            },
            psite_output_paths={
                "+": paths["genome_psite_plus_bigwig"],
                "-": paths["genome_psite_minus_bigwig"],
            },
        )

    write_chrom_sizes(
        definitions,
        transcript_order,
        paths["chrom_sizes"],
    )
    write_transcript_name_map(
        definitions,
        transcript_order,
        paths["name_map"],
    )
    write_transcript_annotation_gtf(
        transcript_definitions=definitions,
        transcript_order=transcript_order,
        reference_gtf=reference_gtf,
        output_path=paths["annotation_gtf"],
    )

    logger.info("Transcript depth workflow finished")
    logger.info("CDS-start profile: %s", paths["profile"])
    logger.info("P-site offsets: %s", paths["offsets"])
    logger.info(
        "CDS-start normalized metagene: %s",
        paths["start_metagene"],
    )
    logger.info(
        "CDS-stop normalized metagene: %s",
        paths["stop_metagene"],
    )
    logger.info("Transcript name map: %s", paths["name_map"])
    logger.info(
        "Transcript annotation GTF: %s",
        paths["annotation_gtf"],
    )

    if create_bigwig:
        logger.info("Transcript depth BigWig: %s", paths["coverage_bigwig"])
        logger.info("P-site depth BigWig: %s", paths["psite_bigwig"])

    if create_genome_bigwig:
        logger.info(
            "Genome transcript depth plus BigWig: %s",
            paths["genome_coverage_plus_bigwig"],
        )
        logger.info(
            "Genome transcript depth minus BigWig: %s",
            paths["genome_coverage_minus_bigwig"],
        )
        logger.info(
            "Genome P-site depth plus BigWig: %s",
            paths["genome_psite_plus_bigwig"],
        )
        logger.info(
            "Genome P-site depth minus BigWig: %s",
            paths["genome_psite_minus_bigwig"],
        )

    return paths