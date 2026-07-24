from __future__ import annotations

import mappy as mp

from mapping.bam_writer import GenomeBamWriter
from mapping.cigar_utils import add_soft_clips
from mapping.coordinate_converter import CoordinateConverter
from mapping.mappy_mapper import MappyMapper
from mapping.mapping_result import ProjectedAlignment, RawAlignment
from mapping.tpm_assigner import assign_weights_by_gene_tpm, load_gene_tpm
import logging
import time
from pathlib import Path

def project_genome_alignment(alignment: RawAlignment, read_length: int) -> ProjectedAlignment:
    is_reverse = alignment.strand == -1
    return ProjectedAlignment(
        read_id=alignment.read_id,
        chrom=alignment.reference_name,
        reference_start=alignment.reference_start,
        cigartuples=add_soft_clips(
            cigartuples=alignment.cigartuples,
            query_start=alignment.query_start,
            query_end=alignment.query_end,
            read_length=read_length,
            is_reverse=is_reverse,
        ),
        is_reverse=is_reverse,
        mapq=alignment.mapq,
        edit_distance=alignment.edit_distance,
        score=alignment.score,
        source_reference=alignment.reference_name,
        source_reference_start=alignment.reference_start,
        source_reference_end=alignment.reference_end,
        support_sources={'genome'},
        assignment_weight=1.0,
    )


def classify_genome_hits(
    hits: list[RawAlignment],
    read_length: int,
) -> tuple[str, RawAlignment | None]:
    if not hits:
        return 'unmapped', None

    best_hit = hits[0]

    if not best_hit.is_full_length(read_length):
        return 'soft_clipped', None

    if len(hits) == 1:
        return 'unique', best_hit

    if best_hit.score == hits[1].score:
        return 'multi_equal_score', None

    return 'unique', best_hit


def select_equal_score_transcript_hits(
    hits: list[RawAlignment],
    read_length: int,
    maximum_candidates: int | None = None,
) -> list[RawAlignment]:
    full_length_hits = [hit for hit in hits if hit.is_full_length(read_length)]
    if not full_length_hits:
        return []

    best_score = full_length_hits[0].score
    equal_score_hits = [hit for hit in full_length_hits if hit.score == best_score]
    if maximum_candidates is None:
        return equal_score_hits

    return equal_score_hits[:maximum_candidates]


def deduplicate_projected_alignments(
    alignments: list[ProjectedAlignment],
) -> list[ProjectedAlignment]:
    deduplicated: dict[tuple, ProjectedAlignment] = {}
    for alignment in alignments:
        key = alignment.deduplication_key
        if key not in deduplicated:
            deduplicated[key] = alignment
        else:
            deduplicated[key].support_sources.update(alignment.support_sources)
    return list(deduplicated.values())

def run_dual_mapping_to_bam(
        fastq_path: str,
        genome_fasta: str,
        transcript_fasta: str,
        coordinate_tsv: str,
        gene_tpm_tsv: str,
        output_bam: str,
        maximum_reads: int | None = None,
        transcript_candidate_count: int = 8,
        kmer_size: int = 9,
        window_size: int = 5,
        minimum_minimizers: int = 1,
        minimum_chain_score: int = 15,
        minimum_dp_score: int = 15,
        minimum_query_coverage: float = 0.85,
        minimum_identity: float = 0.90,
        tpm_pseudocount: float = 0.01,
        index_threads: int = 4,
        sort_threads: int = 4,
        progress_interval: int = 100000,
        log_file: str | None = None,
) -> None:

    genome_mapper = MappyMapper(
        reference_path=genome_fasta,
        source='genome',
        preset='sr',
        kmer_size=kmer_size,
        window_size=window_size,
        minimum_minimizers=minimum_minimizers,
        minimum_chain_score=minimum_chain_score,
        minimum_dp_score=minimum_dp_score,
        maximum_hits=2,
        index_threads=index_threads,
        minimum_query_coverage=minimum_query_coverage,
        minimum_identity=minimum_identity,
    )

    transcript_mapper = MappyMapper(
        reference_path=transcript_fasta,
        source='transcript',
        preset='sr',
        kmer_size=kmer_size,
        window_size=window_size,
        minimum_minimizers=minimum_minimizers,
        minimum_chain_score=minimum_chain_score,
        minimum_dp_score=minimum_dp_score,
        maximum_hits=max(32, transcript_candidate_count * 16),
        index_threads=index_threads,
        minimum_query_coverage=minimum_query_coverage,
        minimum_identity=minimum_identity,
    )

    converter = CoordinateConverter(coordinate_tsv)
    gene_tpm = load_gene_tpm(gene_tpm_tsv)

    read_count = 0
    genome_unique_count = 0
    transcript_single_count = 0
    transcript_multi_count = 0
    transcript_multihit_2_count = 0
    transcript_multihit_3_count = 0
    transcript_multihit_4_count = 0
    transcript_multihit_5_count = 0
    transcript_multihit_6_count = 0
    transcript_multihit_7_count = 0
    transcript_multihit_8_count = 0
    transcript_multihit_over8_count = 0
    unmapped_count = 0
    conversion_errors = 0
    genome_multi_transcript_unmapped_count = 0
    genome_softclip_transcript_unmapped_count = 0
    genome_unmapped_transcript_unmapped_count = 0
    genome_no_hit_count = 0
    genome_single_hit_count = 0
    genome_equal_score_count = 0
    genome_unequal_score_count = 0
    genome_best_not_full_length_count = 0

    with GenomeBamWriter(
        genome_fasta=genome_fasta,
        output_bam=output_bam,
        sort_threads=sort_threads,
    ) as bam_writer:
        for read_id, sequence, quality in mp.fastx_read(fastq_path):
            read_length = len(sequence)
            genome_hits = genome_mapper.map_read(read_id=read_id, sequence=sequence)
            genome_status, genome_best_hit = classify_genome_hits(
                hits=genome_hits,
                read_length=read_length,
            )
            # print(read_id,"genome_hits",len(genome_hits))
            if not genome_hits:
                genome_no_hit_count += 1

            else:
                best_hit = genome_hits[0]

                if not best_hit.is_full_length(
                        read_length
                ):
                    genome_best_not_full_length_count += 1

                elif len(genome_hits) == 1:
                    genome_single_hit_count += 1

                elif (
                        genome_hits[0].score
                        == genome_hits[1].score
                ):
                    genome_equal_score_count += 1

                else:
                    genome_unequal_score_count += 1

            if genome_status == 'unique' and genome_best_hit is not None:
                genome_unique_count += 1
                output_alignments = [
                    project_genome_alignment(
                        alignment=genome_best_hit,
                        read_length=read_length,
                    )
                ]
                unresolved_multi = False

            else:
                transcript_hits = transcript_mapper.map_read(
                    read_id=read_id,
                    sequence=sequence,
                )

                transcript_candidates = select_equal_score_transcript_hits(
                    hits=transcript_hits,
                    read_length=read_length,
                    maximum_candidates=None,
                )
                # print(read_id, "transcript_candidates", len(transcript_candidates))

                output_alignments: list[ProjectedAlignment] = []
                for transcript_hit in transcript_candidates:
                    try:
                        output_alignments.append(
                            converter.project_alignment(
                                alignment=transcript_hit,
                                read_length=read_length,
                            )
                        )
                    except (KeyError, ValueError) as error:
                        conversion_errors += 1
                        if conversion_errors <= 10:
                            print('Coordinate conversion skipped:', read_id, transcript_hit.reference_name, error)

                deduplicated_alignments = deduplicate_projected_alignments(
                    output_alignments
                )
                candidate_count = len(deduplicated_alignments)

                if candidate_count == 0:
                    unmapped_count += 1
                    unresolved_multi = False
                    output_alignments = []

                    if genome_status == "multi_equal_score":
                        genome_multi_transcript_unmapped_count += 1

                    elif genome_status == "soft_clipped":
                        genome_softclip_transcript_unmapped_count += 1

                    elif genome_status == "unmapped":
                        genome_unmapped_transcript_unmapped_count += 1

                elif candidate_count == 1:
                    transcript_single_count += 1
                    output_alignments = deduplicated_alignments
                    output_alignments[0].assignment_weight = 1.0
                    unresolved_multi = False

                elif candidate_count > transcript_candidate_count:
                    transcript_multi_count += 1
                    transcript_multihit_over8_count += 1
                    output_alignments = []
                    unresolved_multi = True

                else:
                    transcript_multi_count += 1
                    output_alignments = deduplicated_alignments

                    if candidate_count == 2:
                        transcript_multihit_2_count += 1
                    elif candidate_count == 3:
                        transcript_multihit_3_count += 1
                    elif candidate_count == 4:
                        transcript_multihit_4_count += 1
                    elif candidate_count == 5:
                        transcript_multihit_5_count += 1
                    elif candidate_count == 6:
                        transcript_multihit_6_count += 1
                    elif candidate_count == 7:
                        transcript_multihit_7_count += 1
                    elif candidate_count == 8:
                        transcript_multihit_8_count += 1

                    assign_weights_by_gene_tpm(
                        alignments=output_alignments,
                        gene_tpm=gene_tpm,
                        pseudocount=tpm_pseudocount,
                    )
                    unresolved_multi = True

            bam_writer.write_read(
                read_id=read_id,
                sequence=sequence,
                quality=quality,
                alignments=output_alignments
            )

            if read_count < 20:
                print(
                    read_id,
                    'length',
                    read_length,
                    'genome_status',
                    genome_status,
                    'genome_hits',
                    len(genome_hits),
                    'output_hits',
                    len(output_alignments),
                )
                for alignment in output_alignments:
                    print(
                        '  ',
                        alignment.source_reference,
                        alignment.gene_id,
                        'score',
                        alignment.score,
                        'weight',
                        round(alignment.assignment_weight, 6),
                    )

            read_count += 1
            if maximum_reads is not None and read_count >= maximum_reads:
                break
            if read_count % 100000 == 0:
                print('Processed reads:', read_count)

    print('Reads:', read_count)
    print('Genome unique:', genome_unique_count)
    print('Transcript single candidate:', transcript_single_count)
    print('Transcript TPM-weighted multi:', transcript_multi_count)
    print('Transcript multi-hit 2:', transcript_multihit_2_count)
    print('Transcript multi-hit 3:', transcript_multihit_3_count)
    print('Transcript multi-hit 4:', transcript_multihit_4_count)
    print('Transcript multi-hit 5:', transcript_multihit_5_count)
    print('Transcript multi-hit 6:', transcript_multihit_6_count)
    print('Transcript multi-hit 7:', transcript_multihit_7_count)
    print('Transcript multi-hit 8:', transcript_multihit_8_count)
    print('Transcript multi-hit over 8:', transcript_multihit_over8_count)
    print('Unmapped:', unmapped_count)
    print('Coordinate conversion errors:', conversion_errors)
    print('BAM:', output_bam)
    print('BAI:', f'{output_bam}.bai')
    print(
        "Genome multi, transcript unmapped:",
        genome_multi_transcript_unmapped_count,
    )

    print(
        "Genome soft-clipped, transcript unmapped:",
        genome_softclip_transcript_unmapped_count,
    )

    print(
        "Genome unmapped, transcript unmapped:",
        genome_unmapped_transcript_unmapped_count,
    )

    print(
        "Genome no hit:",
        genome_no_hit_count,
    )

    print(
        "Genome one hit:",
        genome_single_hit_count,
    )

    print(
        "Genome equal-score top hits:",
        genome_equal_score_count,
    )

    print(
        "Genome unequal-score top hits:",
        genome_unequal_score_count,
    )

    print(
        "Genome best not full length:",
        genome_best_not_full_length_count,
    )