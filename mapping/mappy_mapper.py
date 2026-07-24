from __future__ import annotations

import mappy as mp

from mapping.mapping_result import RawAlignment


class MappyMapper:
    """
    Map single-end Ribo-seq reads with mappy.
    """

    def __init__(
        self,
        reference_path: str,
        source: str,
        preset: str = "sr",
        kmer_size: int = 9,
        window_size: int = 5,
        minimum_minimizers: int = 2,
        minimum_chain_score: int = 15,
        minimum_dp_score: int = 15,
        maximum_hits: int = 5,
        index_threads: int = 4,
        minimum_query_coverage: float = 0.85,
        minimum_identity: float = 0.90,
    ):
        self.source = source
        self.minimum_query_coverage = (
            minimum_query_coverage
        )
        self.minimum_identity = (
            minimum_identity
        )

        self.aligner = mp.Aligner(
            reference_path,
            preset=preset,
            k=kmer_size,
            w=window_size,
            min_cnt=minimum_minimizers,
            min_chain_score=minimum_chain_score,
            min_dp_score=minimum_dp_score,
            best_n=maximum_hits,
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
    ) -> list[RawAlignment]:
        """
        Map one read and return filtered candidates.
        """

        read_length = len(sequence)
        alignments: list[RawAlignment] = []

        for hit in self.aligner.map(
            sequence,
            buf=self.thread_buffer,
            cs=False,
            MD=False,
        ):
            alignment = RawAlignment(
                read_id=read_id,
                source=self.source,
                reference_name=hit.ctg,
                reference_start=hit.r_st,
                reference_end=hit.r_en,
                query_start=hit.q_st,
                query_end=hit.q_en,
                strand=hit.strand,
                mapq=hit.mapq,
                matching_bases=hit.mlen,
                alignment_length=hit.blen,
                edit_distance=hit.NM,
                is_primary=hit.is_primary,

                # mappy:
                #     (length, operation)
                #
                # pysam:
                #     (operation, length)
                cigartuples=[
                    (
                        operation,
                        length,
                    )
                    for length, operation
                    in hit.cigar
                ],
            )

            if (
                alignment.query_coverage(
                    read_length
                )
                < self.minimum_query_coverage
            ):
                continue

            if (
                alignment.identity
                < self.minimum_identity
            ):
                continue

            alignments.append(
                alignment
            )

        alignments.sort(
            key=lambda item: (
                item.score,
                item.mapq,
                item.matching_bases,
                item.is_primary,
            ),
            reverse=True,
        )

        return alignments
