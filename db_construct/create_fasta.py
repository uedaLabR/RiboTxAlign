from __future__ import annotations

import csv
from pathlib import Path

import pysam

from db_construct.data_class import GeneGroup, Interval, ReferenceMapRecord


_COMPLEMENT_TABLE = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT_TABLE)[::-1]


def fetch_genomic_sequence(
    genome: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> str:
    sequence = genome.fetch(chrom, start, end).upper()
    if strand == "-":
        sequence = reverse_complement(sequence)
    return sequence


def _write_fasta_record(
    output_file,
    reference_id: str,
    sequence: str,
    line_width: int = 80,
) -> None:
    output_file.write(f">{reference_id}\n")
    for index in range(0, len(sequence), line_width):
        output_file.write(sequence[index:index + line_width] + "\n")


def create_representative_sequence(
    gene_group: GeneGroup,
    genome: pysam.FastaFile,
) -> tuple[str, str, list[ReferenceMapRecord]] | None:
    transcript_id = gene_group.representative_transcript_id
    if transcript_id is None:
        return None

    transcript = gene_group.transcripts[transcript_id]
    ordered_exons = transcript.ordered_exons()
    if not ordered_exons:
        return None

    reference_id = f"REP|{gene_group.gene_id}|{transcript_id}"
    sequence_parts: list[str] = []
    records: list[ReferenceMapRecord] = []
    reference_position = 0

    for segment_index, (start, end) in enumerate(ordered_exons, start=1):
        sequence = fetch_genomic_sequence(
            genome=genome,
            chrom=gene_group.chrom,
            start=start,
            end=end,
            strand=gene_group.strand,
        )
        segment_end = reference_position + len(sequence)
        records.append(
            ReferenceMapRecord(
                reference_id=reference_id,
                gene_id=gene_group.gene_id,
                record_type="REP",
                segment_index=segment_index,
                reference_start=reference_position,
                reference_end=segment_end,
                chrom=gene_group.chrom,
                genome_start=start,
                genome_end=end,
                strand=gene_group.strand,
                transcript_id=transcript_id,
            )
        )
        sequence_parts.append(sequence)
        reference_position = segment_end

    return reference_id, "".join(sequence_parts), records


def create_alternative_block_sequences(
    gene_group: GeneGroup,
    genome: pysam.FastaFile,
) -> list[tuple[str, str, list[ReferenceMapRecord]]]:
    entries: list[tuple[str, str, list[ReferenceMapRecord]]] = []

    for block_id in sorted(gene_group.alternative_block_ids):
        block = gene_group.get_block(block_id)
        sequence = fetch_genomic_sequence(
            genome=genome,
            chrom=gene_group.chrom,
            start=block.start,
            end=block.end,
            strand=gene_group.strand,
        )
        reference_id = f"ALT_BLOCK|{gene_group.gene_id}|B{block_id}"
        record = ReferenceMapRecord(
            reference_id=reference_id,
            gene_id=gene_group.gene_id,
            record_type="ALT_BLOCK",
            segment_index=1,
            reference_start=0,
            reference_end=len(sequence),
            chrom=gene_group.chrom,
            genome_start=block.start,
            genome_end=block.end,
            strand=gene_group.strand,
            block_id=block_id,
        )
        entries.append((reference_id, sequence, [record]))

    return entries


def _find_supporting_exon_pair(
    gene_group: GeneGroup,
    junction: tuple[int, int],
) -> tuple[str, Interval, Interval] | None:
    for transcript_id, junctions in gene_group.transcript_junctions.items():
        ordered_exons = gene_group.transcripts[transcript_id].ordered_exons()
        for junction_index, current_junction in enumerate(junctions):
            if current_junction != junction:
                continue
            if junction_index + 1 >= len(ordered_exons):
                continue
            return (
                transcript_id,
                ordered_exons[junction_index],
                ordered_exons[junction_index + 1],
            )
    return None


def _upstream_suffix_interval(
    exon: Interval,
    strand: str,
    length: int,
) -> Interval:
    start, end = exon
    length = min(length, end - start)
    if strand == "+":
        return end - length, end
    return start, start + length


def _downstream_prefix_interval(
    exon: Interval,
    strand: str,
    length: int,
) -> Interval:
    start, end = exon
    length = min(length, end - start)
    if strand == "+":
        return start, start + length
    return end - length, end


def create_alternative_junction_sequences(
    gene_group: GeneGroup,
    genome: pysam.FastaFile,
    flank_length: int = 40,
) -> list[tuple[str, str, list[ReferenceMapRecord]]]:
    entries: list[tuple[str, str, list[ReferenceMapRecord]]] = []

    for left_block, right_block in sorted(gene_group.alternative_junctions):
        support = _find_supporting_exon_pair(
            gene_group=gene_group,
            junction=(left_block, right_block),
        )
        if support is None:
            continue

        transcript_id, upstream_exon, downstream_exon = support
        upstream_interval = _upstream_suffix_interval(
            exon=upstream_exon,
            strand=gene_group.strand,
            length=flank_length,
        )
        downstream_interval = _downstream_prefix_interval(
            exon=downstream_exon,
            strand=gene_group.strand,
            length=flank_length,
        )

        upstream_sequence = fetch_genomic_sequence(
            genome=genome,
            chrom=gene_group.chrom,
            start=upstream_interval[0],
            end=upstream_interval[1],
            strand=gene_group.strand,
        )
        downstream_sequence = fetch_genomic_sequence(
            genome=genome,
            chrom=gene_group.chrom,
            start=downstream_interval[0],
            end=downstream_interval[1],
            strand=gene_group.strand,
        )

        junction_position = len(upstream_sequence)
        sequence = upstream_sequence + downstream_sequence
        reference_id = (
            f"ALT_JUNC|{gene_group.gene_id}|"
            f"B{left_block}-B{right_block}|JPOS={junction_position}"
        )

        records = [
            ReferenceMapRecord(
                reference_id=reference_id,
                gene_id=gene_group.gene_id,
                record_type="ALT_JUNC",
                segment_index=1,
                reference_start=0,
                reference_end=len(upstream_sequence),
                chrom=gene_group.chrom,
                genome_start=upstream_interval[0],
                genome_end=upstream_interval[1],
                strand=gene_group.strand,
                transcript_id=transcript_id,
                junction_left_block=left_block,
                junction_right_block=right_block,
            ),
            ReferenceMapRecord(
                reference_id=reference_id,
                gene_id=gene_group.gene_id,
                record_type="ALT_JUNC",
                segment_index=2,
                reference_start=len(upstream_sequence),
                reference_end=len(sequence),
                chrom=gene_group.chrom,
                genome_start=downstream_interval[0],
                genome_end=downstream_interval[1],
                strand=gene_group.strand,
                transcript_id=transcript_id,
                junction_left_block=left_block,
                junction_right_block=right_block,
            ),
        ]
        entries.append((reference_id, sequence, records))

    return entries


def write_reference_fasta(
    gene_groups: dict[str, GeneGroup],
    genome_fasta: str,
    output_fasta: str,
    flank_length: int = 40,
) -> list[ReferenceMapRecord]:
    output_path = Path(output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    coordinate_records: list[ReferenceMapRecord] = []

    with pysam.FastaFile(genome_fasta) as genome:
        with open(output_path, "w", encoding="utf-8") as output_file:
            ordered_genes = sorted(
                gene_groups.values(),
                key=lambda gene: (gene.chrom, gene.start, gene.end, gene.gene_id),
            )

            for gene_group in ordered_genes:
                representative = create_representative_sequence(
                    gene_group=gene_group,
                    genome=genome,
                )
                if representative is not None:
                    reference_id, sequence, records = representative
                    _write_fasta_record(output_file, reference_id, sequence)
                    coordinate_records.extend(records)

                for reference_id, sequence, records in create_alternative_block_sequences(
                    gene_group=gene_group,
                    genome=genome,
                ):
                    _write_fasta_record(output_file, reference_id, sequence)
                    coordinate_records.extend(records)

                for reference_id, sequence, records in create_alternative_junction_sequences(
                    gene_group=gene_group,
                    genome=genome,
                    flank_length=flank_length,
                ):
                    _write_fasta_record(output_file, reference_id, sequence)
                    coordinate_records.extend(records)

    return coordinate_records


def write_coordinate_table(
    records: list[ReferenceMapRecord],
    output_tsv: str,
) -> None:
    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
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
        "transcript_id",
        "block_id",
        "junction_left_block",
        "junction_right_block",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    column: getattr(record, column) if getattr(record, column) is not None else ""
                    for column in columns
                }
            )


def write_gene_tpm_table(
    gene_groups: dict[str, GeneGroup],
    output_tsv: str,
) -> None:
    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "gene_id",
        "gene_name",
        "chrom",
        "strand",
        "gene_tpm",
        "representative_transcript_id",
        "representative_path_tpm",
        "representative_path",
        "transcript_count",
        "alternative_block_count",
        "alternative_junction_count",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for gene_group in sorted(
            gene_groups.values(),
            key=lambda gene: (gene.chrom, gene.start, gene.gene_id),
        ):
            representative_path = ",".join(
                str(block_id) for block_id in (gene_group.representative_path or ())
            )
            writer.writerow(
                {
                    "gene_id": gene_group.gene_id,
                    "gene_name": gene_group.gene_name or "",
                    "chrom": gene_group.chrom,
                    "strand": gene_group.strand,
                    "gene_tpm": gene_group.gene_tpm,
                    "representative_transcript_id": gene_group.representative_transcript_id or "",
                    "representative_path_tpm": gene_group.representative_path_tpm,
                    "representative_path": representative_path,
                    "transcript_count": len(gene_group.transcripts),
                    "alternative_block_count": len(gene_group.alternative_block_ids),
                    "alternative_junction_count": len(gene_group.alternative_junctions),
                }
            )


def write_transcript_path_table(
    gene_groups: dict[str, GeneGroup],
    output_tsv: str,
) -> None:
    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "gene_id",
        "transcript_id",
        "source",
        "reference_id",
        "tpm",
        "path",
        "junctions",
        "is_representative",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for gene_group in sorted(
            gene_groups.values(),
            key=lambda gene: (gene.chrom, gene.start, gene.gene_id),
        ):
            for transcript_id, transcript in gene_group.transcripts.items():
                path = gene_group.transcript_paths.get(transcript_id, ())
                junctions = gene_group.transcript_junctions.get(transcript_id, ())
                writer.writerow(
                    {
                        "gene_id": gene_group.gene_id,
                        "transcript_id": transcript_id,
                        "source": transcript.source,
                        "reference_id": transcript.reference_id or "",
                        "tpm": transcript.tpm,
                        "path": ",".join(str(block_id) for block_id in path),
                        "junctions": ",".join(f"{left}-{right}" for left, right in junctions),
                        "is_representative": int(
                            transcript_id == gene_group.representative_transcript_id
                        ),
                    }
                )
