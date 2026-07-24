from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass(frozen=True)
class SelectedAlignment:
    """One selected genomic placement for a multimap read."""

    read_id: str
    chrom: str
    reference_start: int
    is_reverse: bool
    cigar: str
    genomic_hit_count: int
    genomic_hit_index: int
    candidate_posterior: float
    group_posterior: float
    gene_group_id: str
    gene_id: str
    record_type: str
    source_reference: str
    support_sources: str
    edit_distance: int
    score: int
    sequence_group_size: int
    allocation_error: float
    assignment_method: str


@dataclass(frozen=True)
class BamRedistributionResult:
    """Summary of rebuilding and merging redistributed BAM records."""

    selected_alignment_count: int
    written_mapped_count: int
    written_unmapped_count: int
    missing_multiple_bam_read_count: int
    redistributed_multi_bam: str
    final_bam: str


_CIGAR_CODE_BY_CHARACTER = {
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


def _require_pyarrow_parquet():
    try:
        import pyarrow.parquet as pq
    except ImportError as error:
        raise ImportError(
            "BAM redistribution requires pyarrow. Install it with: "
            "pip install pyarrow"
        ) from error
    return pq


def _parse_cigar(cigar: str) -> list[tuple[int, int]]:
    """Convert a SAM CIGAR string into pysam CIGAR tuples."""

    if not cigar:
        raise ValueError("Selected alignment has an empty CIGAR")

    tuples: list[tuple[int, int]] = []
    number_text = ""

    for character in cigar:
        if character.isdigit():
            number_text += character
            continue

        if not number_text:
            raise ValueError(f"Invalid CIGAR: {cigar}")
        if character not in _CIGAR_CODE_BY_CHARACTER:
            raise ValueError(f"Unsupported CIGAR operation: {character}")

        tuples.append(
            (
                _CIGAR_CODE_BY_CHARACTER[character],
                int(number_text),
            )
        )
        number_text = ""

    if number_text:
        raise ValueError(f"Invalid CIGAR: {cigar}")

    return tuples


def _reverse_complement(sequence: str) -> str:
    translation = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan",
    )
    return sequence.translate(translation)[::-1]


def load_selected_alignments(
    assignments_parquet: str,
) -> dict[str, SelectedAlignment]:
    """Load the single selected candidate row for every assigned read."""

    assignment_path = Path(assignments_parquet)
    if not assignment_path.exists():
        raise FileNotFoundError(assignments_parquet)

    pq = _require_pyarrow_parquet()
    table = pq.read_table(
        str(assignment_path),
        filters=[("selected", "=", True)],
        columns=[
            "read_id",
            "chrom",
            "reference_start",
            "is_reverse",
            "cigar",
            "genomic_hit_count",
            "genomic_hit_index",
            "candidate_posterior",
            "group_posterior",
            "gene_group_id",
            "gene_id",
            "record_type",
            "source_reference",
            "support_sources",
            "edit_distance",
            "score",
            "sequence_group_size",
            "allocation_error",
            "assignment_method",
        ],
    )

    selected: dict[str, SelectedAlignment] = {}

    for row in table.to_pylist():
        alignment = SelectedAlignment(
            read_id=str(row["read_id"]),
            chrom=str(row["chrom"]),
            reference_start=int(row["reference_start"]),
            is_reverse=bool(row["is_reverse"]),
            cigar=str(row["cigar"]),
            genomic_hit_count=int(row["genomic_hit_count"]),
            genomic_hit_index=int(row["genomic_hit_index"]),
            candidate_posterior=float(row["candidate_posterior"]),
            group_posterior=float(row["group_posterior"]),
            gene_group_id=str(row["gene_group_id"]),
            gene_id=str(row.get("gene_id") or ""),
            record_type=str(row.get("record_type") or ""),
            source_reference=str(row.get("source_reference") or ""),
            support_sources=str(row.get("support_sources") or ""),
            edit_distance=int(row.get("edit_distance") or 0),
            score=int(row.get("score") or 0),
            sequence_group_size=int(row.get("sequence_group_size") or 0),
            allocation_error=float(row.get("allocation_error") or 0.0),
            assignment_method=str(row.get("assignment_method") or ""),
        )

        if alignment.read_id in selected:
            raise ValueError(
                "Multiple selected candidates for read: "
                f"{alignment.read_id}"
            )
        selected[alignment.read_id] = alignment

    return selected


def _create_mapped_record(
    source_record: pysam.AlignedSegment,
    selected: SelectedAlignment,
    output_header: pysam.AlignmentHeader,
) -> pysam.AlignedSegment:
    """Create one mapped BAM record from an unmapped multimap record."""

    record = pysam.AlignedSegment(output_header)
    record.query_name = source_record.query_name

    sequence = source_record.query_sequence or ""
    qualities = (
        list(source_record.query_qualities)
        if source_record.query_qualities is not None
        else None
    )

    if selected.is_reverse:
        record.query_sequence = _reverse_complement(sequence)
        if qualities is not None:
            record.query_qualities = qualities[::-1]
        record.flag = 16
    else:
        record.query_sequence = sequence
        if qualities is not None:
            record.query_qualities = qualities
        record.flag = 0

    reference_id = output_header.get_tid(selected.chrom)
    if reference_id < 0:
        raise ValueError(
            f"Chromosome absent from BAM header: {selected.chrom}"
        )

    record.reference_id = reference_id
    record.reference_start = selected.reference_start
    record.mapping_quality = 0
    record.cigartuples = _parse_cigar(selected.cigar)

    record.set_tag("NM", selected.edit_distance)
    record.set_tag("AS", selected.score)
    record.set_tag("NH", selected.genomic_hit_count)
    record.set_tag("HI", selected.genomic_hit_index)
    record.set_tag("ZW", selected.candidate_posterior)
    record.set_tag("ZP", selected.group_posterior)
    record.set_tag("ZG", selected.gene_group_id)
    record.set_tag(
        "ZA",
        selected.assignment_method or "SEQUENCE_GROUP_EM",
    )
    record.set_tag("ZQ", selected.sequence_group_size)
    record.set_tag("ZX", selected.allocation_error)

    if selected.gene_id:
        record.set_tag("ZE", selected.gene_id)
    if selected.record_type:
        record.set_tag("ZT", selected.record_type)
    if selected.source_reference:
        record.set_tag("ZR", selected.source_reference)
    if selected.support_sources:
        record.set_tag("ZS", selected.support_sources)

    return record


def rebuild_multiplemap_bam(
    unique_bam: str,
    multiplemap_bam: str,
    assignments_parquet: str,
    output_bam: str,
    sort_threads: int = 8,
    logger: logging.Logger | None = None,
) -> tuple[int, int, int, int]:
    """Add selected genomic alignments to the stored multimap BAM reads."""

    selected_alignments = load_selected_alignments(assignments_parquet)
    output_path = Path(output_bam)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_bam = Path(f"{output_bam}.unsorted.bam")

    written_mapped_count = 0
    written_unmapped_count = 0
    seen_read_ids: set[str] = set()

    with pysam.AlignmentFile(unique_bam, "rb") as unique_input:
        output_header = unique_input.header

        with (
            pysam.AlignmentFile(multiplemap_bam, "rb") as multi_input,
            pysam.AlignmentFile(
                str(temporary_bam),
                "wb",
                header=output_header,
            ) as output,
        ):
            for source_record in multi_input.fetch(until_eof=True):
                read_id = source_record.query_name
                selected = selected_alignments.get(read_id)

                if selected is None:
                    unmapped_record = pysam.AlignedSegment.fromstring(
                        source_record.to_string(),
                        output_header,
                    )
                    output.write(unmapped_record)
                    written_unmapped_count += 1
                    continue

                output.write(
                    _create_mapped_record(
                        source_record=source_record,
                        selected=selected,
                        output_header=output_header,
                    )
                )
                written_mapped_count += 1
                seen_read_ids.add(read_id)

    missing_multiple_bam_read_count = len(
        set(selected_alignments) - seen_read_ids
    )

    pysam.sort(
        "-@",
        str(sort_threads),
        "-o",
        str(output_path),
        str(temporary_bam),
    )
    pysam.index(
        "-@",
        str(sort_threads),
        str(output_path),
    )
    temporary_bam.unlink(missing_ok=True)

    if logger is not None:
        logger.info(
            "Redistributed multi BAM: mapped=%d unmapped=%d "
            "selected_missing_from_bam=%d output=%s",
            written_mapped_count,
            written_unmapped_count,
            missing_multiple_bam_read_count,
            output_path,
        )

    return (
        len(selected_alignments),
        written_mapped_count,
        written_unmapped_count,
        missing_multiple_bam_read_count,
    )


def merge_unique_and_redistributed_bams(
    unique_bam: str,
    redistributed_multi_bam: str,
    output_bam: str,
    sort_threads: int = 8,
) -> None:
    """Merge unique and redistributed records, then sort and index."""

    output_path = Path(output_bam)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_bam = Path(f"{output_bam}.unsorted.bam")

    with pysam.AlignmentFile(unique_bam, "rb") as unique_input:
        header = unique_input.header

        with pysam.AlignmentFile(
            str(temporary_bam),
            "wb",
            header=header,
        ) as output:
            for record in unique_input.fetch(until_eof=True):
                output.write(record)

            with pysam.AlignmentFile(
                redistributed_multi_bam,
                "rb",
            ) as multi_input:
                for record in multi_input.fetch(until_eof=True):
                    output.write(record)

    pysam.sort(
        "-@",
        str(sort_threads),
        "-o",
        str(output_path),
        str(temporary_bam),
    )
    pysam.index(
        "-@",
        str(sort_threads),
        str(output_path),
    )
    temporary_bam.unlink(missing_ok=True)


def redistribute_and_merge_bams(
    unique_bam: str,
    multiplemap_bam: str,
    assignments_parquet: str,
    redistributed_multi_bam: str,
    final_bam: str,
    sort_threads: int = 8,
    logger: logging.Logger | None = None,
) -> BamRedistributionResult:
    """Rebuild multimap alignments and merge them with the unique BAM."""

    (
        selected_alignment_count,
        written_mapped_count,
        written_unmapped_count,
        missing_multiple_bam_read_count,
    ) = rebuild_multiplemap_bam(
        unique_bam=unique_bam,
        multiplemap_bam=multiplemap_bam,
        assignments_parquet=assignments_parquet,
        output_bam=redistributed_multi_bam,
        sort_threads=sort_threads,
        logger=logger,
    )

    merge_unique_and_redistributed_bams(
        unique_bam=unique_bam,
        redistributed_multi_bam=redistributed_multi_bam,
        output_bam=final_bam,
        sort_threads=sort_threads,
    )

    if logger is not None:
        logger.info("Final merged BAM: %s", final_bam)

    return BamRedistributionResult(
        selected_alignment_count=selected_alignment_count,
        written_mapped_count=written_mapped_count,
        written_unmapped_count=written_unmapped_count,
        missing_multiple_bam_read_count=missing_multiple_bam_read_count,
        redistributed_multi_bam=redistributed_multi_bam,
        final_bam=final_bam,
    )
