from __future__ import annotations

from pathlib import Path

import mappy as mp
import pysam

from mapping.mapping_result import (
    ProjectedAlignment,
)


class GenomeBamWriter:
    """
    Write projected genomic alignments to BAM.

    The output BAM is coordinate sorted and indexed.
    """

    def __init__(
        self,
        genome_fasta: str,
        output_bam: str,
        sort_threads: int = 4,
    ):
        self.genome_fasta = genome_fasta
        self.output_bam = Path(
            output_bam
        )
        self.sort_threads = (
            sort_threads
        )

        self.output_bam.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        fasta_index = Path(
            f"{genome_fasta}.fai"
        )

        if not fasta_index.exists():
            pysam.faidx(
                genome_fasta
            )

        with pysam.FastaFile(
            genome_fasta
        ) as genome:
            header_dict = {
                "HD": {
                    "VN": "1.6",
                    "SO": "unsorted",
                },
                "SQ": [
                    {
                        "SN": reference_name,
                        "LN": reference_length,
                    }
                    for (
                        reference_name,
                        reference_length,
                    ) in zip(
                        genome.references,
                        genome.lengths,
                    )
                ],
                "PG": [
                    {
                        "ID": "RiboTxAlign",
                        "PN": "RiboTxAlign",
                        "VN": "0.1.0",
                    }
                ],
            }

        self.header = (
            pysam.AlignmentHeader.from_dict(
                header_dict
            )
        )

        self.temporary_bam = Path(
            f"{self.output_bam}.unsorted.bam"
        )

        self.bam = pysam.AlignmentFile(
            self.temporary_bam,
            "wb",
            header=self.header,
        )

    def write_read(
        self,
        read_id: str,
        sequence: str,
        quality: str | None,
        alignments: list[
            ProjectedAlignment
        ],
        score_tolerance: int = 2,
    ) -> None:
        """
        Write one read.

        Best alignment:
            primary

        Remaining alignments:
            secondary
        """

        if not alignments:
            self._write_unmapped(
                read_id=read_id,
                sequence=sequence,
                quality=quality,
            )

            return

        alignments = sorted(
            alignments,
            key=lambda item: (
                item.score,
                item.mapq,
            ),
            reverse=True,
        )

        best_score = (
            alignments[0].score
        )

        near_best_count = sum(
            1
            for alignment
            in alignments
            if (
                best_score
                - alignment.score
                <= score_tolerance
            )
        )

        hit_count = len(
            alignments
        )

        for hit_index, alignment in enumerate(
            alignments,
            start=1,
        ):
            is_secondary = (
                hit_index > 1
            )

            if (
                near_best_count > 1
                or is_secondary
            ):
                mapq = 0
            else:
                mapq = (
                    alignment.mapq
                )

            record = self._create_record(
                read_id=read_id,
                sequence=sequence,
                quality=quality,
                alignment=alignment,
                mapq=mapq,
                is_secondary=is_secondary,
                hit_count=hit_count,
                hit_index=hit_index,
            )

            self.bam.write(
                record
            )

    def _create_record(
        self,
        read_id: str,
        sequence: str,
        quality: str | None,
        alignment: ProjectedAlignment,
        mapq: int,
        is_secondary: bool,
        hit_count: int,
        hit_index: int,
    ) -> pysam.AlignedSegment:
        record = pysam.AlignedSegment(
            self.header
        )

        record.query_name = (
            read_id
        )

        if alignment.is_reverse:
            stored_sequence = (
                mp.revcomp(
                    sequence
                )
            )

            stored_quality = (
                quality[::-1]
                if quality
                else None
            )
        else:
            stored_sequence = (
                sequence
            )

            stored_quality = (
                quality
            )

        record.query_sequence = (
            stored_sequence
        )

        if stored_quality:
            record.query_qualities = (
                pysam.qualitystring_to_array(
                    stored_quality
                )
            )

        flag = 0

        if alignment.is_reverse:
            flag |= 16

        if is_secondary:
            flag |= 256

        record.flag = flag

        reference_id = (
            self.header.get_tid(
                alignment.chrom
            )
        )

        if reference_id < 0:
            raise ValueError(
                "Chromosome absent from BAM header: "
                f"{alignment.chrom}"
            )

        record.reference_id = (
            reference_id
        )

        record.reference_start = (
            alignment.reference_start
        )

        record.mapping_quality = max(
            0,
            min(
                60,
                int(mapq),
            ),
        )

        record.cigartuples = (
            alignment.cigartuples
        )

        record.set_tag(
            "NM",
            int(
                alignment.edit_distance
            ),
        )

        record.set_tag(
            "AS",
            int(
                alignment.score
            ),
        )

        record.set_tag(
            "NH",
            int(
                hit_count
            ),
        )

        record.set_tag(
            "HI",
            int(
                hit_index
            ),
        )

        record.set_tag(
            "ZS",
            ",".join(
                sorted(
                    alignment.support_sources
                )
            ),
        )

        record.set_tag(
            "ZR",
            alignment.source_reference,
        )

        record.set_tag(
            "ZB",
            int(
                alignment.source_reference_start
            ),
        )

        record.set_tag(
            "ZE",
            int(
                alignment.source_reference_end
            ),
        )

        if alignment.gene_id:
            record.set_tag(
                "ZG",
                alignment.gene_id,
            )

        if alignment.record_type:
            record.set_tag(
                "ZT",
                alignment.record_type,
            )

        return record

    def _write_unmapped(
        self,
        read_id: str,
        sequence: str,
        quality: str | None,
    ) -> None:
        record = pysam.AlignedSegment(
            self.header
        )

        record.query_name = (
            read_id
        )

        record.query_sequence = (
            sequence
        )

        if quality:
            record.query_qualities = (
                pysam.qualitystring_to_array(
                    quality
                )
            )

        record.flag = 4
        record.reference_id = -1
        record.reference_start = -1
        record.mapping_quality = 0

        self.bam.write(
            record
        )

    def close(self) -> None:
        if not self.bam.closed:
            self.bam.close()

        pysam.sort(
            "-@",
            str(
                self.sort_threads
            ),
            "-o",
            str(
                self.output_bam
            ),
            str(
                self.temporary_bam
            ),
        )

        pysam.index(
            "-@",
            str(
                self.sort_threads
            ),
            str(
                self.output_bam
            ),
        )

        self.temporary_bam.unlink(
            missing_ok=True
        )

    def __enter__(
        self,
    ) -> "GenomeBamWriter":
        return self

    def __exit__(
        self,
        exc_type,
        exc_value,
        traceback,
    ) -> None:
        if not self.bam.closed:
            self.bam.close()

        if exc_type is not None:
            return

        pysam.sort(
            "-@",
            str(
                self.sort_threads
            ),
            "-o",
            str(
                self.output_bam
            ),
            str(
                self.temporary_bam
            ),
        )

        pysam.index(
            "-@",
            str(
                self.sort_threads
            ),
            str(
                self.output_bam
            ),
        )

        self.temporary_bam.unlink(
            missing_ok=True
        )
