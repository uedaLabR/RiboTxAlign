from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class CoordinateSegment:
    """One transcript-reference interval mapped to one genomic exon interval."""

    reference_id: str
    reference_start: int
    reference_end: int
    chrom: str
    genome_start: int
    genome_end: int
    strand: str

    @property
    def length(self) -> int:
        """Return the segment length."""

        return self.reference_end - self.reference_start


class TranscriptGenomeCoordinateConverter:
    """Convert transcript-reference coordinates to genomic coordinates."""

    def __init__(
        self,
        segments_by_reference: dict[str, tuple[CoordinateSegment, ...]],
    ) -> None:
        self._segments_by_reference = segments_by_reference

    @classmethod
    def from_coordinate_table(
        cls,
        coordinate_table: str | Path,
    ) -> "TranscriptGenomeCoordinateConverter":
        """Load a reference coordinate table written by coordinate_writer.py."""

        required_columns = {
            "reference_id",
            "reference_start",
            "reference_end",
            "chrom",
            "genome_start",
            "genome_end",
            "strand",
        }

        segments: dict[str, list[CoordinateSegment]] = defaultdict(list)

        with Path(coordinate_table).open(
            "r",
            encoding="utf-8",
            newline="",
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = set(reader.fieldnames or ())

            missing_columns = sorted(required_columns - fieldnames)
            if missing_columns:
                raise ValueError(
                    "Coordinate table is missing columns: "
                    + ", ".join(missing_columns)
                )

            for row in reader:
                if row.get("record_type", "TRANSCRIPT") != "TRANSCRIPT":
                    continue

                reference_id = row["reference_id"].strip()
                if not reference_id:
                    continue

                segment = CoordinateSegment(
                    reference_id=reference_id,
                    reference_start=int(row["reference_start"]),
                    reference_end=int(row["reference_end"]),
                    chrom=row["chrom"].strip(),
                    genome_start=int(row["genome_start"]),
                    genome_end=int(row["genome_end"]),
                    strand=row["strand"].strip(),
                )

                cls._validate_segment(segment)
                segments[reference_id].append(segment)

        finalized: dict[str, tuple[CoordinateSegment, ...]] = {}

        for reference_id, reference_segments in segments.items():
            ordered = sorted(
                reference_segments,
                key=lambda segment: (
                    segment.reference_start,
                    segment.reference_end,
                ),
            )
            cls._validate_reference_segments(reference_id, ordered)
            finalized[reference_id] = tuple(ordered)

        if not finalized:
            raise RuntimeError(
                "No transcript coordinate segments were loaded."
            )

        return cls(finalized)

    @staticmethod
    def _validate_segment(segment: CoordinateSegment) -> None:
        """Validate one coordinate segment."""

        if segment.strand not in {"+", "-"}:
            raise ValueError(
                f"Invalid strand for {segment.reference_id}: "
                f"{segment.strand}"
            )

        if segment.reference_end <= segment.reference_start:
            raise ValueError(
                f"Invalid reference interval for {segment.reference_id}: "
                f"{segment.reference_start}-{segment.reference_end}"
            )

        if segment.genome_end <= segment.genome_start:
            raise ValueError(
                f"Invalid genome interval for {segment.reference_id}: "
                f"{segment.genome_start}-{segment.genome_end}"
            )

        if (
            segment.reference_end - segment.reference_start
            != segment.genome_end - segment.genome_start
        ):
            raise ValueError(
                f"Length mismatch for {segment.reference_id}: "
                f"reference={segment.reference_start}-"
                f"{segment.reference_end}, "
                f"genome={segment.genome_start}-{segment.genome_end}"
            )

    @staticmethod
    def _validate_reference_segments(
        reference_id: str,
        segments: list[CoordinateSegment],
    ) -> None:
        """Validate segment continuity and common chromosome/strand."""

        chroms = {segment.chrom for segment in segments}
        strands = {segment.strand for segment in segments}

        if len(chroms) != 1 or len(strands) != 1:
            raise ValueError(
                f"Inconsistent chromosome or strand for {reference_id}."
            )

        previous_end = 0

        for index, segment in enumerate(segments):
            if index == 0 and segment.reference_start != 0:
                raise ValueError(
                    f"First reference segment does not start at zero for "
                    f"{reference_id}."
                )

            if segment.reference_start != previous_end:
                raise ValueError(
                    f"Non-contiguous reference segments for {reference_id}: "
                    f"expected {previous_end}, "
                    f"observed {segment.reference_start}."
                )

            previous_end = segment.reference_end

    def get_segments(
        self,
        reference_id: str,
    ) -> tuple[CoordinateSegment, ...]:
        """Return coordinate segments for one reference."""

        return self._segments_by_reference.get(reference_id, ())

    def project_interval(
        self,
        reference_id: str,
        reference_start: int,
        reference_end: int,
    ) -> list[tuple[str, int, int, str]]:
        """Project one half-open transcript interval onto genomic exons."""

        if reference_end <= reference_start:
            return []

        projected: list[tuple[str, int, int, str]] = []

        for segment in self.get_segments(reference_id):
            overlap_start = max(reference_start, segment.reference_start)
            overlap_end = min(reference_end, segment.reference_end)

            if overlap_start >= overlap_end:
                continue

            relative_start = overlap_start - segment.reference_start
            relative_end = overlap_end - segment.reference_start

            if segment.strand == "+":
                genome_start = segment.genome_start + relative_start
                genome_end = segment.genome_start + relative_end
            else:
                genome_start = segment.genome_end - relative_end
                genome_end = segment.genome_end - relative_start

            projected.append(
                (
                    segment.chrom,
                    genome_start,
                    genome_end,
                    segment.strand,
                )
            )

        return projected

    def project_position(
        self,
        reference_id: str,
        reference_position: int,
    ) -> tuple[str, int, str] | None:
        """Project one transcript position to one genomic base position."""

        for segment in self.get_segments(reference_id):
            if not (
                segment.reference_start
                <= reference_position
                < segment.reference_end
            ):
                continue

            relative_position = (
                reference_position - segment.reference_start
            )

            if segment.strand == "+":
                genome_position = (
                    segment.genome_start + relative_position
                )
            else:
                genome_position = (
                    segment.genome_end - 1 - relative_position
                )

            return (
                segment.chrom,
                genome_position,
                segment.strand,
            )

        return None


def read_genome_sizes(
    genome_fai: str | Path,
) -> list[tuple[str, int]]:
    """Read chromosome names and lengths from a FASTA index."""

    sizes: list[tuple[str, int]] = []

    with Path(genome_fai).open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 2:
                raise ValueError(
                    f"Invalid FASTA index line {line_number}: {line.rstrip()}"
                )

            chrom = fields[0]
            length = int(fields[1])

            if length <= 0:
                raise ValueError(
                    f"Invalid chromosome length for {chrom}: {length}"
                )

            sizes.append((chrom, length))

    if not sizes:
        raise RuntimeError("Genome FASTA index contains no sequences.")

    return sizes