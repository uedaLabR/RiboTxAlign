from __future__ import annotations

import csv
from dataclasses import dataclass

import pysam

from mapping.cigar_utils import (
    REFERENCE_CONSUMING_OPERATIONS,
    add_soft_clips,
    merge_adjacent_cigar_operations,
)
from mapping.mapping_result import (
    ProjectedAlignment,
    RawAlignment,
)


@dataclass(frozen=True)
class ReferenceSegment:
    """
    One section of a synthetic reference sequence.

    All coordinates are 0-based and half-open.
    """

    reference_id: str
    gene_id: str
    record_type: str
    segment_index: int

    reference_start: int
    reference_end: int

    chrom: str
    genome_start: int
    genome_end: int
    strand: str

    transcript_id: str = ""
    block_id: str = ""


@dataclass
class ProjectedChunk:
    """
    One part of an alignment mapped to one genomic segment.
    """

    segment: ReferenceSegment
    synthetic_start: int
    synthetic_end: int
    cigartuples: list[tuple[int, int]]

    @property
    def genome_start(self) -> int:
        if self.segment.strand == "+":
            return (
                self.segment.genome_start
                + (
                    self.synthetic_start
                    - self.segment.reference_start
                )
            )

        return (
            self.segment.genome_end
            - (
                self.synthetic_end
                - self.segment.reference_start
            )
        )

    @property
    def genome_end(self) -> int:
        if self.segment.strand == "+":
            return (
                self.segment.genome_start
                + (
                    self.synthetic_end
                    - self.segment.reference_start
                )
            )

        return (
            self.segment.genome_end
            - (
                self.synthetic_start
                - self.segment.reference_start
            )
        )


class CoordinateConverter:
    """
    Convert synthetic transcript-reference alignments
    to genomic coordinates.
    """

    def __init__(
        self,
        coordinate_tsv: str,
    ):
        self.segments_by_reference = (
            self._load_coordinate_table(
                coordinate_tsv
            )
        )

    @staticmethod
    def _load_coordinate_table(
        coordinate_tsv: str,
    ) -> dict[str, list[ReferenceSegment]]:
        segments_by_reference: dict[
            str,
            list[ReferenceSegment],
        ] = {}

        with open(coordinate_tsv) as input_file:
            reader = csv.DictReader(
                input_file,
                delimiter="\t",
            )

            required_columns = {
                "reference_id",
                "gene_id",
                "record_type",
                "segment_index",
                "reference_start",
                "reference_end",
                "chrom",
                "genome_start",
                "genome_end",
                "strand",
            }

            fieldnames = set(
                reader.fieldnames or []
            )

            missing_columns = (
                required_columns
                - fieldnames
            )

            if missing_columns:
                raise ValueError(
                    "Missing coordinate-table columns: "
                    + ", ".join(
                        sorted(
                            missing_columns
                        )
                    )
                )

            for row in reader:
                segment = ReferenceSegment(
                    reference_id=row[
                        "reference_id"
                    ],
                    gene_id=row["gene_id"],
                    record_type=row[
                        "record_type"
                    ],
                    segment_index=int(
                        row["segment_index"]
                    ),
                    reference_start=int(
                        row["reference_start"]
                    ),
                    reference_end=int(
                        row["reference_end"]
                    ),
                    chrom=row["chrom"],
                    genome_start=int(
                        row["genome_start"]
                    ),
                    genome_end=int(
                        row["genome_end"]
                    ),
                    strand=row["strand"],
                    transcript_id=row.get(
                        "transcript_id",
                        "",
                    ),
                    block_id=row.get(
                        "block_id",
                        "",
                    ),
                )

                segments_by_reference.setdefault(
                    segment.reference_id,
                    [],
                ).append(
                    segment
                )

        for reference_id, segments in (
            segments_by_reference.items()
        ):
            segments.sort(
                key=lambda item: (
                    item.reference_start,
                    item.reference_end,
                )
            )

            CoordinateConverter._validate_reference(
                reference_id,
                segments,
            )

        return segments_by_reference

    @staticmethod
    def _validate_reference(
        reference_id: str,
        segments: list[ReferenceSegment],
    ) -> None:
        if not segments:
            return

        expected_chrom = segments[0].chrom
        expected_strand = segments[0].strand
        expected_gene_id = segments[0].gene_id

        previous_end = None

        for segment in segments:
            if segment.chrom != expected_chrom:
                raise ValueError(
                    "Synthetic reference maps to "
                    "multiple chromosomes: "
                    f"{reference_id}"
                )

            if segment.strand != expected_strand:
                raise ValueError(
                    "Synthetic reference has mixed "
                    f"strands: {reference_id}"
                )

            if segment.gene_id != expected_gene_id:
                raise ValueError(
                    "Synthetic reference maps to "
                    f"multiple genes: {reference_id}"
                )

            if (
                previous_end is not None
                and segment.reference_start
                < previous_end
            ):
                raise ValueError(
                    "Synthetic coordinate segments "
                    f"overlap: {reference_id}"
                )

            previous_end = (
                segment.reference_end
            )

    def has_reference(
        self,
        reference_id: str,
    ) -> bool:
        return (
            reference_id
            in self.segments_by_reference
        )

    def convert_position(
        self,
        reference_id: str,
        reference_position: int,
    ) -> tuple[str, int, str]:
        """
        Convert one synthetic reference position.
        """

        segment = self._find_segment(
            reference_id=reference_id,
            reference_position=reference_position,
        )

        offset = (
            reference_position
            - segment.reference_start
        )

        if segment.strand == "+":
            genome_position = (
                segment.genome_start
                + offset
            )
        else:
            genome_position = (
                segment.genome_end
                - 1
                - offset
            )

        return (
            segment.chrom,
            genome_position,
            segment.strand,
        )

    def project_alignment(
        self,
        alignment: RawAlignment,
        read_length: int,
    ) -> ProjectedAlignment:
        """
        Project one synthetic-reference alignment
        to genomic coordinates.

        Segment boundaries are converted to genomic
        intron skips using CIGAR N operations.
        """

        segments = self.segments_by_reference.get(
            alignment.reference_name
        )

        if not segments:
            raise KeyError(
                "Reference absent from coordinate table: "
                f"{alignment.reference_name}"
            )

        chunks = self._split_alignment(
            alignment=alignment,
            segments=segments,
        )

        if not chunks:
            raise ValueError(
                "No projected chunks were produced"
            )

        synthetic_strand = (
            alignment.strand
        )
        gene_strand = (
            chunks[0].segment.strand
        )

        genomic_strand = (
            synthetic_strand
            * (
                1
                if gene_strand == "+"
                else -1
            )
        )

        is_reverse = (
            genomic_strand == -1
        )

        if gene_strand == "+":
            ordered_chunks = chunks
        else:
            ordered_chunks = list(
                reversed(chunks)
            )

            for chunk in ordered_chunks:
                chunk.cigartuples = list(
                    reversed(
                        chunk.cigartuples
                    )
                )

        genomic_cigar: list[
            tuple[int, int]
        ] = []

        for index, chunk in enumerate(
            ordered_chunks
        ):
            genomic_cigar.extend(
                chunk.cigartuples
            )

            if index + 1 >= len(
                ordered_chunks
            ):
                continue

            next_chunk = (
                ordered_chunks[
                    index + 1
                ]
            )

            intron_length = (
                next_chunk.genome_start
                - chunk.genome_end
            )

            if intron_length < 0:
                raise ValueError(
                    "Projected chunks overlap in "
                    "genomic coordinates"
                )

            if intron_length > 0:
                genomic_cigar.append(
                    (
                        pysam.CREF_SKIP,
                        intron_length,
                    )
                )

        genomic_cigar = (
            add_soft_clips(
                cigartuples=genomic_cigar,
                query_start=alignment.query_start,
                query_end=alignment.query_end,
                read_length=read_length,
                is_reverse=is_reverse,
            )
        )

        first_segment = segments[0]

        return ProjectedAlignment(
            read_id=alignment.read_id,
            chrom=first_segment.chrom,
            reference_start=min(
                chunk.genome_start
                for chunk in ordered_chunks
            ),
            cigartuples=genomic_cigar,
            is_reverse=is_reverse,
            mapq=alignment.mapq,
            edit_distance=alignment.edit_distance,
            score=alignment.score,
            source_reference=alignment.reference_name,
            source_reference_start=alignment.reference_start,
            source_reference_end=alignment.reference_end,
            gene_id=first_segment.gene_id,
            record_type=first_segment.record_type,
            support_sources={
                "transcript"
            },
        )

    def _split_alignment(
        self,
        alignment: RawAlignment,
        segments: list[ReferenceSegment],
    ) -> list[ProjectedChunk]:
        """
        Split a synthetic-reference CIGAR at
        coordinate-table segment boundaries.
        """

        chunks: list[ProjectedChunk] = []

        synthetic_position = (
            alignment.reference_start
        )

        current_chunk: (
            ProjectedChunk | None
        ) = None

        for operation, length in (
            alignment.cigartuples
        ):
            if (
                operation
                not in REFERENCE_CONSUMING_OPERATIONS
            ):
                if current_chunk is None:
                    segment = self._find_segment(
                        reference_id=(
                            alignment.reference_name
                        ),
                        reference_position=(
                            synthetic_position
                        ),
                    )

                    current_chunk = (
                        ProjectedChunk(
                            segment=segment,
                            synthetic_start=(
                                synthetic_position
                            ),
                            synthetic_end=(
                                synthetic_position
                            ),
                            cigartuples=[],
                        )
                    )

                    chunks.append(
                        current_chunk
                    )

                current_chunk.cigartuples.append(
                    (
                        operation,
                        length,
                    )
                )

                continue

            remaining_length = length

            while remaining_length > 0:
                segment = self._find_segment(
                    reference_id=(
                        alignment.reference_name
                    ),
                    reference_position=(
                        synthetic_position
                    ),
                )

                segment_available = (
                    segment.reference_end
                    - synthetic_position
                )

                take_length = min(
                    remaining_length,
                    segment_available,
                )

                if (
                    current_chunk is None
                    or current_chunk.segment
                    != segment
                ):
                    current_chunk = (
                        ProjectedChunk(
                            segment=segment,
                            synthetic_start=(
                                synthetic_position
                            ),
                            synthetic_end=(
                                synthetic_position
                            ),
                            cigartuples=[],
                        )
                    )

                    chunks.append(
                        current_chunk
                    )

                current_chunk.cigartuples.append(
                    (
                        operation,
                        take_length,
                    )
                )

                synthetic_position += (
                    take_length
                )

                current_chunk.synthetic_end = (
                    synthetic_position
                )

                remaining_length -= (
                    take_length
                )

        if (
            synthetic_position
            != alignment.reference_end
        ):
            raise ValueError(
                "CIGAR/reference length mismatch: "
                f"{synthetic_position} != "
                f"{alignment.reference_end}"
            )

        for chunk in chunks:
            chunk.cigartuples = (
                merge_adjacent_cigar_operations(
                    chunk.cigartuples
                )
            )

        return chunks

    def _find_segment(
        self,
        reference_id: str,
        reference_position: int,
    ) -> ReferenceSegment:
        segments = self.segments_by_reference.get(
            reference_id
        )

        if not segments:
            raise KeyError(
                "Reference absent from coordinate table: "
                f"{reference_id}"
            )

        for segment in segments:
            if (
                segment.reference_start
                <= reference_position
                < segment.reference_end
            ):
                return segment

        raise ValueError(
            "Position is not covered by the "
            "coordinate table: "
            f"{reference_id}:"
            f"{reference_position}"
        )
