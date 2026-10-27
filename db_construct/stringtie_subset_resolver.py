from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass

from db_construct.data_class import Transcript


@dataclass(frozen=True)
class GeneCdsModel:
    """Outermost GENCODE CDS window and exonic support for one gene."""

    gene_key: str
    gene_id: str
    chrom: str
    strand: str
    cds_start: int
    cds_end: int
    exons: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class TranscriptCollapseAnnotation:
    """Information used only to decide whether StringTie models may collapse."""

    gene_key: str | None
    cds_signature: tuple[tuple[int, int], ...]
    cds_resolved: bool
    candidate_gene_count: int


class _UnionFind:
    """Union-find that prevents collapse across conflicting genes/CDS models."""

    def __init__(
        self,
        transcripts: list[Transcript],
        annotations: dict[str, TranscriptCollapseAnnotation],
    ) -> None:
        self.parent = list(range(len(transcripts)))
        self.rank = [0] * len(transcripts)
        self.gene_keys: list[set[str]] = []
        self.cds_keys: list[
            set[tuple[str, tuple[tuple[int, int], ...]]]
        ] = []

        for transcript in transcripts:
            annotation = annotations[transcript.transcript_id]

            gene_keys: set[str] = set()
            if annotation.gene_key is not None:
                gene_keys.add(annotation.gene_key)

            cds_keys: set[
                tuple[str, tuple[tuple[int, int], ...]]
            ] = set()
            if (
                annotation.cds_resolved
                and annotation.gene_key is not None
                and annotation.cds_signature
            ):
                cds_keys.add(
                    (
                        annotation.gene_key,
                        annotation.cds_signature,
                    )
                )

            self.gene_keys.append(gene_keys)
            self.cds_keys.append(cds_keys)

    def find(self, index: int) -> int:
        if self.parent[index] != index:
            self.parent[index] = self.find(self.parent[index])
        return self.parent[index]

    def union(self, first: int, second: int) -> bool:
        first_root = self.find(first)
        second_root = self.find(second)

        if first_root == second_root:
            return True

        combined_gene_keys = (
            self.gene_keys[first_root]
            | self.gene_keys[second_root]
        )
        if len(combined_gene_keys) > 1:
            return False

        combined_cds_keys = (
            self.cds_keys[first_root]
            | self.cds_keys[second_root]
        )
        if len(combined_cds_keys) > 1:
            return False

        if self.rank[first_root] < self.rank[second_root]:
            first_root, second_root = second_root, first_root

        self.parent[second_root] = first_root
        self.gene_keys[first_root] = combined_gene_keys
        self.cds_keys[first_root] = combined_cds_keys

        if self.rank[first_root] == self.rank[second_root]:
            self.rank[first_root] += 1

        return True


def normalize_gene_id(value: str | None) -> str | None:
    """Normalize an Ensembl-like gene ID while retaining non-Ensembl IDs."""

    if value is None:
        return None

    text = str(value).strip()
    if not text:
        return None

    if text.startswith("gene:"):
        text = text[len("gene:"):]

    return text.split(".", 1)[0]


def get_tpm(transcript: Transcript) -> float:
    """Return transcript TPM as a finite numeric value when possible."""

    value = getattr(transcript, "tpm", 0.0)

    if value is None:
        return 0.0

    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def get_transcript_length(transcript: Transcript) -> int:
    """Return total spliced exon length."""

    return sum(
        end - start
        for start, end in getattr(
            transcript,
            "exons",
            [],
        )
        if end > start
    )


def get_exons_in_genomic_order(
    transcript: Transcript,
) -> tuple[tuple[int, int], ...]:
    """Return unique valid exons in genomic coordinate order."""

    return tuple(
        sorted(
            {
                (int(start), int(end))
                for start, end in getattr(
                    transcript,
                    "exons",
                    [],
                )
                if int(end) > int(start)
            }
        )
    )


def get_exons_in_transcript_order(
    transcript: Transcript,
) -> tuple[tuple[int, int], ...]:
    """Return exons in biological 5'-to-3' order."""

    exons = list(
        get_exons_in_genomic_order(
            transcript
        )
    )

    if transcript.strand == "-":
        exons.reverse()

    return tuple(exons)


def intervals_overlap(
    first: tuple[int, int],
    second: tuple[int, int],
) -> bool:
    """Return True when two half-open intervals overlap."""

    return (
        first[0] < second[1]
        and second[0] < first[1]
    )


def exon_chains_overlap(
    first_exons: tuple[tuple[int, int], ...],
    second_exons: tuple[tuple[int, int], ...],
) -> bool:
    """Return True when two exon chains share at least one genomic base."""

    first_index = 0
    second_index = 0

    while (
        first_index < len(first_exons)
        and second_index < len(second_exons)
    ):
        first = first_exons[first_index]
        second = second_exons[second_index]

        if intervals_overlap(first, second):
            return True

        if first[1] <= second[0]:
            first_index += 1
        else:
            second_index += 1

    return False


def build_gencode_gene_cds_models(
    gencode_transcripts: list[Transcript],
) -> dict[str, GeneCdsModel]:
    """
    Build one conservative gene-level CDS window from GENCODE.

    The requested rule is used here: CDS start/end are the outermost
    coordinates observed across all coding transcripts of the gene.
    """

    grouped: dict[str, list[Transcript]] = defaultdict(list)

    for transcript in gencode_transcripts:
        gene_key = normalize_gene_id(
            getattr(transcript, "gene_id", None)
        )
        if gene_key is None:
            continue

        if not getattr(transcript, "cds_regions", None):
            continue

        grouped[gene_key].append(transcript)

    models: dict[str, GeneCdsModel] = {}

    for gene_key, transcripts in grouped.items():
        chromosomes = {
            transcript.chrom
            for transcript in transcripts
        }
        strands = {
            transcript.strand
            for transcript in transcripts
        }

        if len(chromosomes) != 1 or len(strands) != 1:
            continue

        cds_regions = [
            (int(start), int(end))
            for transcript in transcripts
            for start, end in getattr(
                transcript,
                "cds_regions",
                [],
            )
            if int(end) > int(start)
        ]
        if not cds_regions:
            continue

        gene_exons = tuple(
            sorted(
                {
                    (int(start), int(end))
                    for transcript in transcripts
                    for start, end in getattr(
                        transcript,
                        "exons",
                        [],
                    )
                    if int(end) > int(start)
                }
            )
        )

        representative_gene_id = next(
            (
                str(transcript.gene_id)
                for transcript in transcripts
                if getattr(transcript, "gene_id", None)
            ),
            gene_key,
        )

        models[gene_key] = GeneCdsModel(
            gene_key=gene_key,
            gene_id=representative_gene_id,
            chrom=next(iter(chromosomes)),
            strand=next(iter(strands)),
            cds_start=min(start for start, _ in cds_regions),
            cds_end=max(end for _, end in cds_regions),
            exons=gene_exons,
        )

    return models


def clip_exons_to_cds_window(
    transcript: Transcript,
    cds_start: int,
    cds_end: int,
) -> tuple[tuple[int, int], ...]:
    """
    Return the transcript exon structure inside an outermost CDS window.

    Differences outside this window are treated as UTR-only differences.
    """

    clipped: list[tuple[int, int]] = []

    for exon_start, exon_end in get_exons_in_genomic_order(
        transcript
    ):
        start = max(exon_start, cds_start)
        end = min(exon_end, cds_end)

        if start < end:
            clipped.append((start, end))

    return tuple(clipped)


def get_gene_exons_inside_cds_window(
    model: GeneCdsModel,
) -> tuple[tuple[int, int], ...]:
    """Return GENCODE gene exonic support only inside its outermost CDS window.

    The outermost CDS window remains the gene-level definition used by this
    resolver. UTR-only exons outside that window must not contribute to gene
    assignment or transcript uniqueness for reference collapsing.
    """

    clipped: list[
        tuple[int, int]
    ] = []

    for exon_start, exon_end in model.exons:
        start = max(
            int(exon_start),
            int(model.cds_start),
        )
        end = min(
            int(exon_end),
            int(model.cds_end),
        )

        if start < end:
            clipped.append(
                (
                    start,
                    end,
                )
            )

    return tuple(
        sorted(
            set(
                clipped
            )
        )
    )


def infer_gencode_gene_for_stringtie(
    transcript: Transcript,
    gene_models: dict[str, GeneCdsModel],
) -> tuple[str | None, int]:
    """
    Infer one GENCODE gene for a StringTie transcript.

    A valid StringTie ref_gene_id has priority.

    Otherwise, candidate genes are evaluated only with exon sequence inside
    each gene's outermost CDS window. UTR-only exons outside that window are
    deliberately ignored. This prevents a transcript-specific UTR exon from
    making a coding-gene assignment ambiguous and subsequently blocking the
    same-CDS reference collapse.
    """

    direct_gene_key = normalize_gene_id(
        getattr(
            transcript,
            "reference_gene_id",
            None,
        )
    )

    if direct_gene_key is not None:
        model = gene_models.get(
            direct_gene_key
        )
        if (
            model is not None
            and model.chrom
            == transcript.chrom
            and model.strand
            == transcript.strand
        ):
            return direct_gene_key, 1

    candidate_gene_keys: list[
        str
    ] = []

    for gene_key, model in (
        gene_models.items()
    ):
        if (
            model.chrom
            != transcript.chrom
            or model.strand
            != transcript.strand
        ):
            continue

        transcript_cds_exons = (
            clip_exons_to_cds_window(
                transcript=transcript,
                cds_start=model.cds_start,
                cds_end=model.cds_end,
            )
        )

        if not transcript_cds_exons:
            continue

        gene_cds_exons = (
            get_gene_exons_inside_cds_window(
                model
            )
        )

        if not gene_cds_exons:
            continue

        if exon_chains_overlap(
            transcript_cds_exons,
            gene_cds_exons,
        ):
            candidate_gene_keys.append(
                gene_key
            )

    candidate_gene_keys = sorted(
        set(
            candidate_gene_keys
        )
    )

    if len(
        candidate_gene_keys
    ) == 1:
        return (
            candidate_gene_keys[0],
            1,
        )

    return (
        None,
        len(
            candidate_gene_keys
        ),
    )


def annotate_stringtie_transcripts(
    stringtie_transcripts: list[Transcript],
    gencode_transcripts: list[Transcript],
) -> tuple[
    list[Transcript],
    dict[str, TranscriptCollapseAnnotation],
]:
    """Attach conservative GENCODE-derived CDS windows to StringTie copies.

    Both gene inference and CDS-equivalence are evaluated within the gene-level
    outermost CDS window. Exon differences outside that window are UTR-only
    for the purpose of reference collapsing.
    """

    gene_models = build_gencode_gene_cds_models(
        gencode_transcripts
    )

    annotated_transcripts: list[Transcript] = []
    annotations: dict[
        str,
        TranscriptCollapseAnnotation,
    ] = {}

    for original in stringtie_transcripts:
        transcript = deepcopy(original)

        gene_key, candidate_gene_count = (
            infer_gencode_gene_for_stringtie(
                transcript=transcript,
                gene_models=gene_models,
            )
        )

        cds_signature: tuple[
            tuple[int, int], ...
        ] = ()
        cds_resolved = False

        if gene_key is not None:
            model = gene_models[gene_key]
            cds_signature = clip_exons_to_cds_window(
                transcript=transcript,
                cds_start=model.cds_start,
                cds_end=model.cds_end,
            )

            if cds_signature:
                cds_resolved = True
                transcript.cds_regions = list(
                    cds_signature
                )

                if not getattr(
                    transcript,
                    "reference_gene_id",
                    None,
                ):
                    transcript.reference_gene_id = (
                        model.gene_id
                    )

        annotations[transcript.transcript_id] = (
            TranscriptCollapseAnnotation(
                gene_key=gene_key,
                cds_signature=cds_signature,
                cds_resolved=cds_resolved,
                candidate_gene_count=(
                    candidate_gene_count
                ),
            )
        )
        annotated_transcripts.append(transcript)

    return annotated_transcripts, annotations


def exact_exon_structure_matches(
    first: Transcript,
    second: Transcript,
) -> bool:
    """Return True when the complete exon structures are identical."""

    return (
        first.chrom == second.chrom
        and first.strand == second.strand
        and get_exons_in_genomic_order(first)
        == get_exons_in_genomic_order(second)
    )


def terminal_exons_only_differ(
    first: Transcript,
    second: Transcript,
) -> bool:
    """
    Return True when only first/last exon outer boundaries may differ.

    Internal exons must be identical. The splice-connected boundary of the
    biological first and last exons must remain identical. This deliberately
    does not collapse transcripts that gain or lose an entire terminal exon.
    """

    if first.chrom != second.chrom:
        return False
    if first.strand != second.strand:
        return False

    first_exons = get_exons_in_transcript_order(
        first
    )
    second_exons = get_exons_in_transcript_order(
        second
    )

    if len(first_exons) != len(second_exons):
        return False

    # Without CDS information, single-exon models are kept separate unless
    # they were already caught by exact_exon_structure_matches().
    if len(first_exons) < 2:
        return False

    if first_exons[1:-1] != second_exons[1:-1]:
        return False

    first_first_exon = first_exons[0]
    second_first_exon = second_exons[0]
    first_last_exon = first_exons[-1]
    second_last_exon = second_exons[-1]

    if first.strand == "+":
        first_exon_splice_boundary_matches = (
            first_first_exon[1]
            == second_first_exon[1]
        )
        last_exon_splice_boundary_matches = (
            first_last_exon[0]
            == second_last_exon[0]
        )
    else:
        first_exon_splice_boundary_matches = (
            first_first_exon[0]
            == second_first_exon[0]
        )
        last_exon_splice_boundary_matches = (
            first_last_exon[1]
            == second_last_exon[1]
        )

    return (
        first_exon_splice_boundary_matches
        and last_exon_splice_boundary_matches
    )


def known_genes_are_compatible(
    first_annotation: TranscriptCollapseAnnotation,
    second_annotation: TranscriptCollapseAnnotation,
) -> bool:
    """Reject a pair only when both have different confidently inferred genes."""

    if (
        first_annotation.gene_key is not None
        and second_annotation.gene_key is not None
        and first_annotation.gene_key
        != second_annotation.gene_key
    ):
        return False

    return True


def classify_collapse_pair(
    first: Transcript,
    second: Transcript,
    first_annotation: TranscriptCollapseAnnotation,
    second_annotation: TranscriptCollapseAnnotation,
) -> str | None:
    """Return the allowed collapse reason, or None when models must stay separate."""

    if first.chrom != second.chrom:
        return None
    if first.strand != second.strand:
        return None

    if not known_genes_are_compatible(
        first_annotation,
        second_annotation,
    ):
        return None

    if exact_exon_structure_matches(
        first,
        second,
    ):
        return "stringtie_exact_exon"

    if (
        first_annotation.cds_resolved
        and second_annotation.cds_resolved
    ):
        if (
            first_annotation.gene_key
            == second_annotation.gene_key
            and first_annotation.cds_signature
            == second_annotation.cds_signature
        ):
            return "stringtie_same_inferred_cds"

        # CDS could be assigned to both transcripts, but differs. Do not use
        # the terminal-exon fallback in this case.
        return None

    # CDS is not uniquely/fully inferable for at least one transcript.
    # Collapse only when the difference is restricted to biological terminal
    # exon outer boundaries.
    if terminal_exons_only_differ(
        first,
        second,
    ):
        return "stringtie_terminal_utr_like"

    return None


def build_locus_groups(
    transcripts: list[Transcript],
) -> list[list[Transcript]]:
    """Build broad StringTie locus groups from overlapping genomic spans."""

    sorted_transcripts = sorted(
        transcripts,
        key=lambda transcript: (
            transcript.chrom,
            transcript.strand,
            transcript.start
            if transcript.start is not None
            else -1,
            transcript.end
            if transcript.end is not None
            else -1,
            transcript.transcript_id,
        ),
    )

    locus_groups: list[list[Transcript]] = []
    current_group: list[Transcript] = []
    current_chrom: str | None = None
    current_strand: str | None = None
    current_end: int | None = None

    for transcript in sorted_transcripts:
        if (
            transcript.start is None
            or transcript.end is None
        ):
            continue

        same_reference = (
            transcript.chrom == current_chrom
            and transcript.strand == current_strand
        )
        overlaps_current_span = (
            current_end is not None
            and transcript.start < current_end
        )

        if (
            current_group
            and same_reference
            and overlaps_current_span
        ):
            current_group.append(transcript)
            current_end = max(
                current_end,
                transcript.end,
            )
            continue

        if current_group:
            locus_groups.append(current_group)

        current_group = [transcript]
        current_chrom = transcript.chrom
        current_strand = transcript.strand
        current_end = transcript.end

    if current_group:
        locus_groups.append(current_group)

    return locus_groups


def choose_representative(
    members: list[Transcript],
) -> Transcript:
    """
    Keep the longest structure so Ribo-seq reads are not lost.

    TPM is only the second criterion, so a short major isoform cannot replace
    a longer compatible structure. The merged abundance is calculated
    separately as the sum of member TPMs.
    """

    return max(
        members,
        key=lambda transcript: (
            get_transcript_length(transcript),
            get_tpm(transcript),
            transcript.transcript_id,
        ),
    )


def flatten_source_ids(
    members: list[Transcript],
) -> list[str]:
    """Return all original StringTie transcript IDs represented by a group."""

    source_ids: set[str] = set()

    for transcript in members:
        existing = getattr(
            transcript,
            "source_transcript_ids",
            None,
        )

        if existing:
            source_ids.update(
                str(value)
                for value in existing
                if value is not None
            )
        else:
            source_ids.add(
                str(transcript.transcript_id)
            )

    return sorted(source_ids)


def determine_group_resolution_type(
    members: list[Transcript],
    annotations: dict[str, TranscriptCollapseAnnotation],
) -> str:
    """Assign one readable reason to a completed collapse group."""

    exon_structures = {
        get_exons_in_genomic_order(transcript)
        for transcript in members
    }
    if len(exon_structures) == 1:
        return "stringtie_exact_exon"

    resolved_cds_keys = {
        (
            annotation.gene_key,
            annotation.cds_signature,
        )
        for transcript in members
        for annotation in [
            annotations[transcript.transcript_id]
        ]
        if annotation.cds_resolved
        and annotation.gene_key is not None
        and annotation.cds_signature
    }

    if (
        len(resolved_cds_keys) == 1
        and all(
            annotations[
                transcript.transcript_id
            ].cds_resolved
            for transcript in members
        )
    ):
        return "stringtie_same_inferred_cds"

    return "stringtie_terminal_utr_like"


def merge_source_information(
    members: list[Transcript],
    annotations: dict[str, TranscriptCollapseAnnotation],
) -> Transcript:
    """Keep the longest structure and attach summed abundance/merge metadata."""

    representative = choose_representative(
        members
    )
    merged = deepcopy(representative)

    source_ids = flatten_source_ids(
        members
    )

    merged.source_transcript_ids = source_ids
    merged.stringtie_transcript_ids = source_ids
    merged.collapsed_transcript_count = len(
        source_ids
    )
    merged.representative_transcript_id = (
        representative.transcript_id
    )
    merged.resolution_type = (
        determine_group_resolution_type(
            members=members,
            annotations=annotations,
        )
    )

    # Structural representative and abundance prior are intentionally
    # separated: keep the longest sequence, but sum expression across all
    # compatible StringTie models.
    merged.tpm = sum(
        get_tpm(transcript)
        for transcript in members
    )

    return merged


def collapse_one_locus(
    transcripts: list[Transcript],
    annotations: dict[str, TranscriptCollapseAnnotation],
) -> list[Transcript]:
    """Collapse only exact/CDS-equivalent/terminal-UTR-like StringTie models."""

    if len(transcripts) <= 1:
        return [
            deepcopy(transcript)
            for transcript in transcripts
        ]

    union_find = _UnionFind(
        transcripts=transcripts,
        annotations=annotations,
    )

    pair_candidates: dict[
        str,
        list[tuple[int, int]],
    ] = {
        "stringtie_exact_exon": [],
        "stringtie_same_inferred_cds": [],
        "stringtie_terminal_utr_like": [],
    }

    for first_index in range(len(transcripts)):
        for second_index in range(
            first_index + 1,
            len(transcripts),
        ):
            first = transcripts[first_index]
            second = transcripts[second_index]

            reason = classify_collapse_pair(
                first=first,
                second=second,
                first_annotation=(
                    annotations[
                        first.transcript_id
                    ]
                ),
                second_annotation=(
                    annotations[
                        second.transcript_id
                    ]
                ),
            )

            if reason is not None:
                pair_candidates[reason].append(
                    (
                        first_index,
                        second_index,
                    )
                )

    # Strongest evidence first. The union-find metadata prevents a terminal
    # fallback from bridging two incompatible known CDS structures.
    for reason in (
        "stringtie_exact_exon",
        "stringtie_same_inferred_cds",
        "stringtie_terminal_utr_like",
    ):
        for first_index, second_index in (
            pair_candidates[reason]
        ):
            union_find.union(
                first_index,
                second_index,
            )

    members_by_root: dict[
        int,
        list[Transcript],
    ] = defaultdict(list)

    for index, transcript in enumerate(
        transcripts
    ):
        members_by_root[
            union_find.find(index)
        ].append(transcript)

    collapsed: list[Transcript] = []

    for members in members_by_root.values():
        if len(members) == 1:
            collapsed.append(
                deepcopy(members[0])
            )
        else:
            collapsed.append(
                merge_source_information(
                    members=members,
                    annotations=annotations,
                )
            )

    return collapsed


def collapse_stringtie_subsets(
    stringtie_transcripts: list[Transcript],
    gencode_transcripts: list[Transcript] | None = None,
) -> list[Transcript]:
    """
    Conservatively collapse redundant StringTie transcript models.

    Collapse rules, in order:
      1. Exact exon structure.
      2. Same inferred CDS structure using the outermost GENCODE CDS window
         of one uniquely assigned gene; UTR-only differences are ignored.
      3. If CDS cannot be uniquely inferred for at least one model, allow
         only first/last-exon outer-boundary differences with identical
         internal exon structure.

    For every collapsed group, the longest transcript structure is retained
    and member TPMs are summed. This preserves the broadest mapping target
    without letting a low-TPM long isoform incorrectly lower the abundance
    prior of a major compatible isoform.
    """

    if gencode_transcripts is None:
        gencode_transcripts = []

    (
        annotated_transcripts,
        annotations,
    ) = annotate_stringtie_transcripts(
        stringtie_transcripts=(
            stringtie_transcripts
        ),
        gencode_transcripts=(
            gencode_transcripts
        ),
    )

    locus_groups = build_locus_groups(
        annotated_transcripts
    )

    collapsed_transcripts: list[Transcript] = []

    for locus_group in locus_groups:
        collapsed_transcripts.extend(
            collapse_one_locus(
                transcripts=locus_group,
                annotations=annotations,
            )
        )

    collapsed_transcripts.sort(
        key=lambda transcript: (
            transcript.chrom,
            transcript.strand,
            transcript.start
            if transcript.start is not None
            else -1,
            transcript.end
            if transcript.end is not None
            else -1,
            transcript.transcript_id,
        )
    )

    resolution_counts: dict[str, int] = (
        defaultdict(int)
    )
    for transcript in collapsed_transcripts:
        if (
            getattr(
                transcript,
                "collapsed_transcript_count",
                1,
            )
            > 1
        ):
            resolution_counts[
                str(
                    getattr(
                        transcript,
                        "resolution_type",
                        "collapsed",
                    )
                )
            ] += 1

    resolved_cds_count = sum(
        annotation.cds_resolved
        for annotation in annotations.values()
    )
    ambiguous_gene_count = sum(
        annotation.candidate_gene_count > 1
        for annotation in annotations.values()
    )

    print(
        "StringTie conservative collapse:",
        "before=",
        len(stringtie_transcripts),
        "after=",
        len(collapsed_transcripts),
        "removed=",
        len(stringtie_transcripts)
        - len(collapsed_transcripts),
        "loci=",
        len(locus_groups),
    )
    print(
        "StringTie GENCODE CDS inference:",
        "resolved=",
        resolved_cds_count,
        "ambiguous_gene=",
        ambiguous_gene_count,
    )

    if resolution_counts:
        print(
            "StringTie collapse groups by rule:"
        )
        for reason, count in sorted(
            resolution_counts.items()
        ):
            print(
                "  ",
                reason,
                count,
            )

    return collapsed_transcripts
