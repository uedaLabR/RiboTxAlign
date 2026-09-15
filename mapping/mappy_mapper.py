from __future__ import annotations

import mappy as mp

from mapping.mapping_result import RawAlignment
from mapping.cds_policy import load_mapping_cds_regions, prefer_cds_hits


class MappyMapper:
    """Map single-end Ribo-seq reads against a transcript reference."""

    def __init__(
        self,
        reference_path: str,
        preset: str = "sr",
        kmer_size: int = 11,
        window_size: int = 5,
        minimum_minimizers: int = 2,
        minimum_chain_score: int = 15,
        minimum_dp_score: int = 15,
        maximum_candidates: int = 16,
        index_threads: int = 1,
        minimum_query_coverage: float = 0.85,
        minimum_identity: float = 0.90,
        reference_metadata_tsv: str | None = None,
        reference_gtf: str | None = None,
    ) -> None:
        if maximum_candidates < 1:
            raise ValueError("maximum_candidates must be at least 1.")

        self.maximum_candidates = maximum_candidates
        self.minimum_query_coverage = minimum_query_coverage
        self.minimum_identity = minimum_identity
        self.cds_regions = load_mapping_cds_regions(
            reference_path, reference_metadata_tsv, reference_gtf,
        )

        # Request one extra candidate so a best-score overflow can be
        # detected when more than maximum_candidates equal best hits exist.
        self.aligner = mp.Aligner(
            reference_path,
            preset=preset,
            k=kmer_size,
            w=window_size,
            min_cnt=minimum_minimizers,
            min_chain_score=minimum_chain_score,
            min_dp_score=minimum_dp_score,
            best_n=maximum_candidates + 1,
            n_threads=index_threads,
        )

        if not self.aligner:
            raise RuntimeError(
                "Failed to load mappy reference: "
                f"{reference_path}"
            )

        self.thread_buffer = mp.ThreadBuffer()

    def map_read(
        self,
        read_id: str,
        sequence: str,
    ) -> tuple[list[RawAlignment], bool]:
        """Return best-score ties and whether they exceeded the limit.

        Query coverage and identity filters are applied first. Among the
        surviving alignments, only alignments with the maximum internal score
        are retained. Lower-scoring alignments are discarded. Within these
        ties, CDS hits take priority before the candidate limit is applied.
        """

        read_length = len(sequence)
        filtered: list[RawAlignment] = []

        for hit in self.aligner.map(
            sequence,
            buf=self.thread_buffer,
            cs=False,
            MD=False,
        ):

            # Transcript references are stored in the 5' to 3' direction.
            # Retain only forward-strand alignments.
            if hit.strand != 1:
                continue

            alignment = RawAlignment(
                read_id=read_id,
                reference_name=hit.ctg,
                reference_start=hit.r_st,
                reference_end=hit.r_en,
                query_start=hit.q_st,
                query_end=hit.q_en,
                strand=hit.strand,
                matching_bases=hit.mlen,
                alignment_length=hit.blen,
                edit_distance=hit.NM,
                cigartuples=[
                    (operation, length)
                    for length, operation in hit.cigar
                ],
            )

            if alignment.query_coverage(read_length) < self.minimum_query_coverage:
                continue

            if alignment.identity < self.minimum_identity:
                continue

            filtered.append(alignment)

        if not filtered:
            return [], False

        best_score = max(alignment.score for alignment in filtered)
        best_hits = [
            alignment
            for alignment in filtered
            if alignment.score == best_score
        ]

        best_hits = prefer_cds_hits(best_hits, self.cds_regions)
        best_hits.sort(
            key=lambda alignment: (
                alignment.reference_name,
                alignment.reference_start,
                alignment.reference_end,
                alignment.strand,
            )
        )

        overflow = len(best_hits) > self.maximum_candidates
        return best_hits[: self.maximum_candidates], overflow
