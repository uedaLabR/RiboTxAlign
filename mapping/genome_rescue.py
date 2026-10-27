from __future__ import annotations

import gzip
import logging
import time
from dataclasses import dataclass
from pathlib import Path

import mappy as mp

from bam_output.bam_writer import (
    normalize_fastq_read_id,
    open_text_auto,
    read_genome_fai,
    reverse_complement,
)


_MATCH_OPERATIONS = {0, 7, 8}
_INTRON_OPERATION = 3


@dataclass(slots=True)
class RescueCounters:
    """Counters produced by contiguous and splice genome rescue."""

    processed_reads: int = 0
    contiguous_rescued_unique_reads: int = 0
    splice_rescued_unique_reads: int = 0

    @property
    def rescued_unique_reads(self) -> int:
        """Return the total number of uniquely rescued reads."""

        return (
            self.contiguous_rescued_unique_reads
            + self.splice_rescued_unique_reads
        )

    @property
    def remaining_unmapped_reads(self) -> int:
        """Return reads that were not uniquely rescued."""

        return self.processed_reads - self.rescued_unique_reads


def _log_rescue_progress(
    logger: logging.Logger,
    counters: RescueCounters,
    started_at: float,
) -> None:
    """Log rescue counts, elapsed time, and throughput."""

    elapsed_seconds = max(time.monotonic() - started_at, 1e-9)
    reads_per_second = counters.processed_reads / elapsed_seconds
    logger.info(
        (
            "Genome rescue progress: processed=%d contiguous=%d "
            "splice=%d unresolved=%d elapsed=%.1fs rate=%.1f reads/s"
        ),
        counters.processed_reads,
        counters.contiguous_rescued_unique_reads,
        counters.splice_rescued_unique_reads,
        counters.remaining_unmapped_reads,
        elapsed_seconds,
        reads_per_second,
    )


def _iter_fastq_records_with_header(fastq_path: str | Path):
    """Yield the original FASTQ header together with parsed fields."""

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
                header.rstrip("\r\n"),
                normalize_fastq_read_id(header),
                sequence.rstrip("\r\n"),
                quality.rstrip("\r\n"),
            )
            record_index += 1


def _query_coverage(hit, read_length: int) -> float:
    """Return the fraction of the query covered by one alignment."""

    if read_length <= 0:
        return 0.0
    return (hit.q_en - hit.q_st) / read_length


def _identity(hit) -> float:
    """Return alignment identity reported by mappy."""

    if hit.blen <= 0:
        return 0.0
    return hit.mlen / hit.blen


def _score(hit) -> int:
    """Return the score used to identify equal best genome hits."""

    return int(hit.mlen) - int(hit.NM)


def _has_intron(hit) -> bool:
    """Return True when one alignment CIGAR contains an intron."""

    return any(operation == _INTRON_OPERATION for _, operation in hit.cigar)


def _splice_geometry(
    hit,
    minimum_anchor_length: int,
    minimum_intron_length: int,
    maximum_intron_length: int,
) -> tuple[int, int, int] | None:
    """Validate one-junction geometry and return anchor/intron lengths."""

    intron_indices = [
        index
        for index, (_, operation) in enumerate(hit.cigar)
        if operation == _INTRON_OPERATION
    ]
    if len(intron_indices) != 1:
        return None

    intron_index = intron_indices[0]
    if intron_index == 0 or intron_index == len(hit.cigar) - 1:
        return None

    left_length, left_operation = hit.cigar[intron_index - 1]
    intron_length, _ = hit.cigar[intron_index]
    right_length, right_operation = hit.cigar[intron_index + 1]

    if left_operation not in _MATCH_OPERATIONS:
        return None
    if right_operation not in _MATCH_OPERATIONS:
        return None
    if left_length < minimum_anchor_length:
        return None
    if right_length < minimum_anchor_length:
        return None
    if not minimum_intron_length <= intron_length <= maximum_intron_length:
        return None

    return int(left_length), int(right_length), int(intron_length)


def _intron_reference_interval(hit) -> tuple[int, int] | None:
    """Return the zero-based reference interval occupied by one intron."""

    reference_position = int(hit.r_st)
    intron_interval = None

    for length, operation in hit.cigar:
        if operation == _INTRON_OPERATION:
            if intron_interval is not None:
                return None
            intron_interval = (
                reference_position,
                reference_position + int(length),
            )
            reference_position += int(length)
        elif operation in {0, 2, 7, 8}:
            reference_position += int(length)

    return intron_interval


def _is_canonical_junction(
    aligner: mp.Aligner,
    hit,
    allowed_motifs: set[tuple[str, str]],
) -> bool:
    """Return True when the genomic splice motif is canonical."""

    interval = _intron_reference_interval(hit)
    if interval is None:
        return False

    intron_start, intron_end = interval
    if intron_end - intron_start < 4:
        return False

    donor = aligner.seq(hit.ctg, intron_start, intron_start + 2)
    acceptor = aligner.seq(hit.ctg, intron_end - 2, intron_end)
    if donor is None or acceptor is None:
        return False

    donor = donor.upper()
    acceptor = acceptor.upper()

    if hit.strand == -1:
        donor, acceptor = (
            reverse_complement(acceptor),
            reverse_complement(donor),
        )

    return (donor, acceptor) in allowed_motifs


def _best_unique_hit(
    aligner: mp.Aligner,
    thread_buffer: mp.ThreadBuffer,
    sequence: str,
    minimum_query_coverage: float,
    minimum_identity: float,
    require_splice: bool,
    minimum_anchor_length: int,
    minimum_intron_length: int,
    maximum_intron_length: int,
    allowed_motifs: set[tuple[str, str]],
):
    """Return one unique best accepted hit or None when tied/unsupported."""

    filtered_hits = []

    for hit in aligner.map(sequence, buf=thread_buffer, cs=False, MD=False):
        if _query_coverage(hit, len(sequence)) < minimum_query_coverage:
            continue
        if _identity(hit) < minimum_identity:
            continue

        if require_splice:
            if not _has_intron(hit):
                continue
            if _splice_geometry(
                hit,
                minimum_anchor_length=minimum_anchor_length,
                minimum_intron_length=minimum_intron_length,
                maximum_intron_length=maximum_intron_length,
            ) is None:
                continue
            if not _is_canonical_junction(
                aligner=aligner,
                hit=hit,
                allowed_motifs=allowed_motifs,
            ):
                continue
        elif _has_intron(hit):
            continue

        filtered_hits.append(hit)

    if not filtered_hits:
        return None

    best_score = max(_score(hit) for hit in filtered_hits)
    best_hits = [hit for hit in filtered_hits if _score(hit) == best_score]
    if len(best_hits) != 1:
        return None

    return best_hits[0]


def _create_bam_record(
    pysam_module,
    read_id: str,
    sequence: str,
    quality: str,
    hit,
    reference_id_by_name: dict[str, int],
    assignment_tag: str,
):
    """Create one BAM record from a mappy alignment."""

    if hit.ctg not in reference_id_by_name:
        raise ValueError(f"Genome hit is absent from FAI: {hit.ctg}")

    read_length = len(sequence)

    if not 0 <= hit.q_st <= hit.q_en <= read_length:
        raise ValueError(
            f"Invalid query interval for {read_id}: "
            f"q_st={hit.q_st}, q_en={hit.q_en}, "
            f"read_length={read_length}"
        )

    aligned_cigar = [
        (operation, int(length))
        for length, operation in hit.cigar
    ]

    # mappy reports only the aligned query interval.
    # Add terminal soft clips so that the BAM CIGAR consumes
    # the complete stored query sequence.
    if hit.strand == -1:
        left_soft_clip = read_length - int(hit.q_en)
        right_soft_clip = int(hit.q_st)
    else:
        left_soft_clip = int(hit.q_st)
        right_soft_clip = read_length - int(hit.q_en)

    cigar_tuples = []

    if left_soft_clip > 0:
        cigar_tuples.append((4, left_soft_clip))

    cigar_tuples.extend(aligned_cigar)

    if right_soft_clip > 0:
        cigar_tuples.append((4, right_soft_clip))

    record = pysam_module.AlignedSegment()
    record.query_name = read_id
    record.flag = 16 if hit.strand == -1 else 0
    record.reference_id = reference_id_by_name[hit.ctg]
    record.reference_start = int(hit.r_st)
    record.mapping_quality = 60
    record.cigartuples = cigar_tuples

    if hit.strand == -1:
        record.query_sequence = reverse_complement(sequence)
        record.query_qualities = pysam_module.qualitystring_to_array(
            quality[::-1]
        )
    else:
        record.query_sequence = sequence
        record.query_qualities = pysam_module.qualitystring_to_array(quality)

    record.set_tag("XA", assignment_tag, value_type="Z")
    record.set_tag("NH", 1, value_type="i")
    record.set_tag("NM", int(hit.NM), value_type="i")
    record.set_tag("AS", _score(hit), value_type="i")

    return record


def _create_unmapped_bam_record(
    pysam_module,
    read_id: str,
    sequence: str,
    quality: str,
):
    """Create one unmapped BAM record."""

    record = pysam_module.AlignedSegment()
    record.query_name = read_id
    record.flag = 4
    record.reference_id = -1
    record.reference_start = -1
    record.mapping_quality = 0
    record.cigartuples = None
    record.query_sequence = sequence
    record.query_qualities = pysam_module.qualitystring_to_array(quality)

    record.set_tag("XA", "genome_unmapped", value_type="Z")
    record.set_tag("NH", 0, value_type="i")

    return record


def rescue_unique_genome_alignments(
    unmapped_fastq: str | Path,
    genome_fasta: str | Path,
    genome_fai: str | Path,
    output_bam: str | Path,
    remaining_unmapped_fastq: str | Path | None = None,
    contiguous_preset: str = "sr",
    splice_preset: str = "splice",
    kmer_size: int = 11,
    window_size: int = 5,
    minimum_minimizers: int = 2,
    minimum_chain_score: int = 15,
    minimum_dp_score: int = 15,
    minimum_query_coverage: float = 0.85,
    minimum_identity: float = 0.90,
    enable_splice_rescue: bool = True,
    minimum_anchor_length: int = 11,
    minimum_intron_length: int = 20,
    maximum_intron_length: int = 2000,
    canonical_motifs: tuple[str, ...] = ("GT-AG", "GC-AG", "AT-AC"),
    progress_interval: int = 100_000,
    logger: logging.Logger | None = None,
) -> RescueCounters:
    """Rescue unique contiguous hits, then conservative splice hits."""

    try:
        import pysam
    except ImportError as error:
        raise RuntimeError("Genome rescue requires pysam.") from error

    if logger is not None:
        logger.info(
            "Loading contiguous genome aligner: fasta=%s preset=%s k=%d w=%d",
            genome_fasta,
            contiguous_preset,
            kmer_size,
            window_size,
        )

    contiguous_aligner = mp.Aligner(
        str(genome_fasta),
        preset=contiguous_preset,
        k=kmer_size,
        w=window_size,
        min_cnt=minimum_minimizers,
        min_chain_score=minimum_chain_score,
        min_dp_score=minimum_dp_score,
        best_n=2,
        n_threads=1,
    )
    if not contiguous_aligner:
        raise RuntimeError(f"Failed to load genome reference: {genome_fasta}")

    if logger is not None:
        logger.info("Contiguous genome aligner loaded.")

    splice_aligner = None
    if enable_splice_rescue:
        if logger is not None:
            logger.info(
                (
                    "Loading splice genome aligner: fasta=%s preset=%s "
                    "anchor>=%d intron=%d-%d"
                ),
                genome_fasta,
                splice_preset,
                minimum_anchor_length,
                minimum_intron_length,
                maximum_intron_length,
            )

        splice_aligner = mp.Aligner(
            str(genome_fasta),
            preset=splice_preset,
            k=kmer_size,
            w=window_size,
            min_cnt=minimum_minimizers,
            min_chain_score=minimum_chain_score,
            min_dp_score=minimum_dp_score,
            best_n=2,
            n_threads=1,
        )
        if not splice_aligner:
            raise RuntimeError(
                f"Failed to load splice genome reference: {genome_fasta}"
            )

        if logger is not None:
            logger.info("Splice genome aligner loaded.")
    elif logger is not None:
        logger.info("Splice genome rescue is disabled.")

    contiguous_buffer = mp.ThreadBuffer()
    splice_buffer = mp.ThreadBuffer() if splice_aligner is not None else None
    references = read_genome_fai(genome_fai)
    reference_id_by_name = {
        name: reference_id
        for reference_id, (name, _) in enumerate(references)
    }
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": name, "LN": length} for name, length in references],
        "PG": [{
            "ID": "RiboTxAlignGenomeRescue",
            "PN": "RiboTxAlign",
            "VN": "0.5.0",
            "CL": "mappy contiguous and conservative splice rescue",
        }],
    }

    allowed_motifs = set()
    for motif in canonical_motifs:
        fields = motif.upper().split("-")
        if len(fields) != 2 or len(fields[0]) != 2 or len(fields[1]) != 2:
            raise ValueError(f"Invalid canonical splice motif: {motif}")
        allowed_motifs.add((fields[0], fields[1]))

    destination = Path(output_bam)
    destination.parent.mkdir(parents=True, exist_ok=True)
    remaining_path = (
        Path(remaining_unmapped_fastq)
        if remaining_unmapped_fastq is not None
        else None
    )
    if remaining_path is not None:
        remaining_path.parent.mkdir(parents=True, exist_ok=True)

    counters = RescueCounters()
    started_at = time.monotonic()

    if logger is not None:
        logger.info("Genome rescue started.")
        logger.info("Input unmapped FASTQ: %s", unmapped_fastq)
        logger.info("Output rescue BAM: %s", destination)
        logger.info("Progress interval: %d reads", progress_interval)

    unresolved_handle = (
        gzip.open(remaining_path, "wt", encoding="utf-8")
        if remaining_path is not None
        else None
    )

    try:
        with pysam.AlignmentFile(destination, "wb", header=header) as bam_handle:
            for header_line, read_id, sequence, quality in (
                _iter_fastq_records_with_header(unmapped_fastq)
            ):
                counters.processed_reads += 1

                hit = _best_unique_hit(
                    aligner=contiguous_aligner,
                    thread_buffer=contiguous_buffer,
                    sequence=sequence,
                    minimum_query_coverage=minimum_query_coverage,
                    minimum_identity=minimum_identity,
                    require_splice=False,
                    minimum_anchor_length=minimum_anchor_length,
                    minimum_intron_length=minimum_intron_length,
                    maximum_intron_length=maximum_intron_length,
                    allowed_motifs=allowed_motifs,
                )
                assignment_tag = "genome_rescued_unique"

                if hit is not None:
                    counters.contiguous_rescued_unique_reads += 1
                elif splice_aligner is not None and splice_buffer is not None:
                    hit = _best_unique_hit(
                        aligner=splice_aligner,
                        thread_buffer=splice_buffer,
                        sequence=sequence,
                        minimum_query_coverage=minimum_query_coverage,
                        minimum_identity=minimum_identity,
                        require_splice=True,
                        minimum_anchor_length=minimum_anchor_length,
                        minimum_intron_length=minimum_intron_length,
                        maximum_intron_length=maximum_intron_length,
                        allowed_motifs=allowed_motifs,
                    )
                    assignment_tag = "genome_splice_rescued_unique"
                    if hit is not None:
                        counters.splice_rescued_unique_reads += 1

                if hit is None:
                    bam_handle.write(
                        _create_unmapped_bam_record(
                            pysam_module=pysam,
                            read_id=read_id,
                            sequence=sequence,
                            quality=quality,
                        )
                    )

                    if unresolved_handle is not None:
                        unresolved_handle.write(header_line.rstrip("\n") + "\n")
                        unresolved_handle.write(sequence + "\n+\n")
                        unresolved_handle.write(quality + "\n")
                else:
                    bam_handle.write(_create_bam_record(
                        pysam_module=pysam,
                        read_id=read_id,
                        sequence=sequence,
                        quality=quality,
                        hit=hit,
                        reference_id_by_name=reference_id_by_name,
                        assignment_tag=assignment_tag,
                    ))

                if (
                    logger is not None
                    and progress_interval > 0
                    and counters.processed_reads % progress_interval == 0
                ):
                    _log_rescue_progress(
                        logger=logger,
                        counters=counters,
                        started_at=started_at,
                    )
    finally:
        if unresolved_handle is not None:
            unresolved_handle.close()

    if logger is not None:
        _log_rescue_progress(
            logger=logger,
            counters=counters,
            started_at=started_at,
        )
        logger.info("Genome rescue completed.")
        logger.info("Genome rescue FASTQ: %s", unmapped_fastq)
        logger.info("Genome rescue BAM: %s", destination)
        logger.info("Genome rescue processed reads: %d", counters.processed_reads)
        logger.info(
            "Contiguous genome-rescued unique reads: %d",
            counters.contiguous_rescued_unique_reads,
        )
        logger.info(
            "Splice genome-rescued unique reads: %d",
            counters.splice_rescued_unique_reads,
        )
        logger.info(
            "Genome rescue remaining unmapped reads: %d",
            counters.remaining_unmapped_reads,
        )
        if remaining_path is not None:
            logger.info("Final unmapped FASTQ: %s", remaining_path)

    return counters