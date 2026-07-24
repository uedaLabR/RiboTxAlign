from __future__ import annotations

from collections import defaultdict

from db_construct.data_class import ExonBlock, GeneGroup, Junction, Transcript, TranscriptPath


def build_exon_blocks(gene_group: GeneGroup) -> None:
    events: dict[int, int] = defaultdict(int)
    for transcript in gene_group.transcripts.values():
        for start, end in transcript.exons:
            if end <= start:
                continue
            events[start] += 1
            events[end] -= 1

    if not events:
        gene_group.exon_blocks = []
        return

    positions = sorted(events)
    active = 0
    previous_position = positions[0]
    genomic_blocks: list[tuple[int, int]] = []

    for position in positions:
        if active > 0 and previous_position < position:
            genomic_blocks.append((previous_position, position))
        active += events[position]
        previous_position = position

    if gene_group.strand == "-":
        genomic_blocks.reverse()

    gene_group.exon_blocks = [
        ExonBlock(block_id=index, start=start, end=end)
        for index, (start, end) in enumerate(genomic_blocks, start=1)
    ]


def _blocks_for_exon(
    gene_group: GeneGroup,
    exon_start: int,
    exon_end: int,
) -> tuple[int, ...]:
    blocks = [
        block
        for block in gene_group.exon_blocks
        if block.start >= exon_start and block.end <= exon_end
    ]
    blocks.sort(
        key=lambda block: block.start,
        reverse=(gene_group.strand == "-"),
    )
    return tuple(block.block_id for block in blocks)


def build_transcript_paths(gene_group: GeneGroup) -> None:
    gene_group.transcript_paths.clear()
    gene_group.transcript_exon_blocks.clear()
    gene_group.transcript_junctions.clear()

    for transcript_id, transcript in gene_group.transcripts.items():
        exon_block_groups: list[tuple[int, ...]] = []

        for exon_start, exon_end in transcript.ordered_exons():
            block_ids = _blocks_for_exon(
                gene_group=gene_group,
                exon_start=exon_start,
                exon_end=exon_end,
            )
            if block_ids:
                exon_block_groups.append(block_ids)

        path = tuple(
            block_id
            for exon_blocks in exon_block_groups
            for block_id in exon_blocks
        )

        junctions: list[Junction] = []
        for first_exon, second_exon in zip(exon_block_groups, exon_block_groups[1:]):
            junctions.append((first_exon[-1], second_exon[0]))

        gene_group.transcript_paths[transcript_id] = path
        gene_group.transcript_exon_blocks[transcript_id] = tuple(exon_block_groups)
        gene_group.transcript_junctions[transcript_id] = tuple(junctions)


def calculate_gene_tpm(gene_group: GeneGroup) -> None:
    gene_group.gene_tpm = sum(
        transcript.tpm
        for transcript in gene_group.transcripts.values()
        if transcript.source == "stringtie"
    )


def calculate_path_tpm(gene_group: GeneGroup) -> None:
    path_tpm: dict[TranscriptPath, float] = defaultdict(float)
    for transcript_id, path in gene_group.transcript_paths.items():
        if path:
            path_tpm[path] += gene_group.transcripts[transcript_id].tpm
    gene_group.path_tpm = dict(path_tpm)


def _path_exonic_length(gene_group: GeneGroup, path: TranscriptPath) -> int:
    return sum(gene_group.get_block(block_id).length for block_id in path)


def _path_has_gencode_transcript(gene_group: GeneGroup, path: TranscriptPath) -> bool:
    return any(
        gene_group.transcripts[transcript_id].source == "gencode"
        for transcript_id, transcript_path in gene_group.transcript_paths.items()
        if transcript_path == path
    )


def _path_max_cds_length(gene_group: GeneGroup, path: TranscriptPath) -> int:
    lengths = [
        gene_group.transcripts[transcript_id].cds_length
        for transcript_id, transcript_path in gene_group.transcript_paths.items()
        if transcript_path == path
    ]
    return max(lengths, default=0)


def select_representative_path(gene_group: GeneGroup) -> None:
    available_paths = [path for path in gene_group.path_tpm if path]
    if not available_paths:
        gene_group.representative_path = None
        gene_group.representative_transcript_id = None
        gene_group.representative_path_tpm = 0.0
        return

    representative_path = max(
        available_paths,
        key=lambda path: (
            gene_group.path_tpm.get(path, 0.0),
            _path_has_gencode_transcript(gene_group, path),
            _path_max_cds_length(gene_group, path),
            _path_exonic_length(gene_group, path),
        ),
    )

    gene_group.representative_path = representative_path
    gene_group.representative_path_tpm = gene_group.path_tpm.get(representative_path, 0.0)
    candidate_ids = [
        transcript_id
        for transcript_id, path in gene_group.transcript_paths.items()
        if path == representative_path
    ]

    if not candidate_ids:
        gene_group.representative_transcript_id = None
        return

    def transcript_score(transcript: Transcript) -> tuple[int, int, float, int, int]:
        return (
            int(transcript.source == "gencode" and transcript.cds_length > 0),
            int(transcript.source == "gencode"),
            transcript.tpm,
            transcript.cds_length,
            transcript.exonic_length,
        )

    gene_group.representative_transcript_id = max(
        candidate_ids,
        key=lambda transcript_id: transcript_score(gene_group.transcripts[transcript_id]),
    )


def find_alternative_blocks(gene_group: GeneGroup) -> None:
    representative_blocks = set(gene_group.representative_path or ())
    all_blocks = {
        block_id
        for path in gene_group.transcript_paths.values()
        for block_id in path
    }
    gene_group.alternative_block_ids = all_blocks - representative_blocks


def find_alternative_junctions(gene_group: GeneGroup) -> None:
    all_junctions = {
        junction
        for junctions in gene_group.transcript_junctions.values()
        for junction in junctions
    }
    representative_junctions: set[Junction] = set()
    representative_transcript_id = gene_group.representative_transcript_id
    if representative_transcript_id is not None:
        representative_junctions = set(
            gene_group.transcript_junctions.get(representative_transcript_id, ())
        )
    gene_group.alternative_junctions = all_junctions - representative_junctions


def resolve_gene(gene_group: GeneGroup) -> None:
    calculate_gene_tpm(gene_group)
    build_exon_blocks(gene_group)
    build_transcript_paths(gene_group)
    calculate_path_tpm(gene_group)
    select_representative_path(gene_group)
    find_alternative_blocks(gene_group)
    find_alternative_junctions(gene_group)


def resolve_all_genes(gene_groups: dict[str, GeneGroup]) -> None:
    for gene_group in gene_groups.values():
        resolve_gene(gene_group)
