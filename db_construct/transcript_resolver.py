from __future__ import annotations

from copy import deepcopy

from db_construct.data_class import Transcript

UTR_PADDING = 0  # nt added outside the GENCODE UTR boundary (no extension)
MIN_EXON_LENGTH = 10  # exons at or below this length (nt) are "micro-exons"
MAX_GENCODE_EXONS_SPANNED = 8  # a StringTie exon overlapping >= this many
# GENCODE exons signals a large structural mismatch
LOW_TPM_FRACTION = 0.05  # drop transcripts at or below this fraction of the
# max TPM within their gene group


def is_stringtie_transcript(transcript: Transcript) -> bool:
    """Return True when the transcript was derived from StringTie."""
    source = str(getattr(transcript, "source", "")).lower()
    return "stringtie" in source or transcript.transcript_id.startswith(
        "MSTRG."
    )


def is_gencode_transcript(transcript: Transcript) -> bool:
    """Return True when the transcript is a GENCODE annotation record."""
    source = str(getattr(transcript, "source", "")).lower()
    return "gencode" in source or transcript.transcript_id.startswith("ENST")


def _has_no_micro_exons(transcript: Transcript) -> bool:
    """False when the transcript contains any exon <= MIN_EXON_LENGTH nt.

    Transcripts containing even one such micro-exon are discarded
    outright, as they are likely assembly artifacts.
    """
    return all(
        (end - start) > MIN_EXON_LENGTH for start, end in transcript.exons
    )


def _junctions_five_to_three(
    transcript: Transcript,
) -> tuple[tuple[int, int], ...]:
    """Return internal splice junctions ordered 5' -> 3' along the transcript.

    ``_internal_junctions`` returns junctions in ascending genomic order;
    on the - strand that is 3' -> 5', so it is reversed here to get
    transcription-direction order.
    """
    junctions = _internal_junctions(transcript)
    if transcript.strand == "-":
        return tuple(reversed(junctions))
    return junctions


def _is_three_prime_anchored_subset(
    candidate: Transcript,
    reference: Transcript,
) -> bool:
    """True when ``candidate`` is a 3'-anchored internal subset of ``reference``.

    "3'-anchored subset" means: candidate's internal splice junctions,
    read 5' -> 3', exactly match a run of reference's junctions that
    reaches all the way to reference's 3' end. This lets candidate be
    truncated at the 5' end (the unreliable end of a StringTie assembly)
    while still requiring its splicing to agree with reference from
    wherever it starts through to the shared 3' end.

    Single-exon candidates have no junctions to compare, so they are
    matched purely on genomic containment plus reaching reference's 3'
    end coordinate.
    """
    if candidate.chrom != reference.chrom or candidate.strand != reference.strand:
        return False
    if (
        candidate.start is None
        or candidate.end is None
        or reference.start is None
        or reference.end is None
    ):
        return False
    if not (reference.start <= candidate.start and candidate.end <= reference.end):
        return False

    candidate_junctions = _junctions_five_to_three(candidate)
    reference_junctions = _junctions_five_to_three(reference)

    three_prime_coord = (
        reference.end if reference.strand != "-" else reference.start
    )
    candidate_three_prime_coord = (
        candidate.end if candidate.strand != "-" else candidate.start
    )

    if not candidate_junctions:
        return candidate_three_prime_coord == three_prime_coord

    if len(candidate_junctions) > len(reference_junctions):
        return False

    return candidate_junctions == reference_junctions[-len(candidate_junctions):]


def _find_gencode_by_structure(
    transcript: Transcript,
    gencode_candidates: list[Transcript],
) -> Transcript | None:
    """Find a GENCODE transcript that 3'-structurally contains ``transcript``.

    Used as a fallback when no ID-based match (``reference_id`` /
    ``gencode_transcript_ids``) is available. Among GENCODE transcripts in
    the same overlap group, picks the one for which ``transcript`` is a
    3'-anchored internal subset (see ``_is_three_prime_anchored_subset``),
    preferring the match with the most shared internal junctions (i.e.
    the tightest match) and breaking ties by transcript_id.
    """
    best: Transcript | None = None
    best_junction_count = -1

    for gencode_transcript in gencode_candidates:
        if not _is_three_prime_anchored_subset(transcript, gencode_transcript):
            continue

        junction_count = len(_internal_junctions(gencode_transcript))
        if junction_count > best_junction_count or (
            junction_count == best_junction_count
            and (
                best is None
                or gencode_transcript.transcript_id < best.transcript_id
            )
        ):
            best = gencode_transcript
            best_junction_count = junction_count

    return best


def _find_matched_gencode(
    transcript: Transcript,
    gencode_by_id: dict[str, Transcript],
) -> tuple[Transcript | None, bool]:
    """Find the GENCODE transcript ``transcript`` was matched against.

    Looks first at ``reference_id``, then falls back to the first entry
    of ``gencode_transcript_ids``. If neither yields a match, falls back
    to structural matching: a GENCODE transcript in the same overlap
    group whose 3' internal splice structure contains ``transcript`` as a
    3'-anchored subset (see ``_find_gencode_by_structure``).

    Returns ``(matched_transcript, is_structural_match)``. A structural
    match only guarantees agreement from the matched region through to
    the shared 3' end; it says nothing about UTR length or 5' extent, so
    callers should treat it as a full replacement by the GENCODE
    definition rather than a partial boundary patch.
    """
    candidate_ids: list[str] = []
    if transcript.reference_id:
        candidate_ids.append(transcript.reference_id)
    candidate_ids.extend(transcript.gencode_transcript_ids)

    for candidate_id in candidate_ids:
        matched = gencode_by_id.get(candidate_id)
        if matched is not None:
            return matched, False

    structural_match = _find_gencode_by_structure(
        transcript, list(gencode_by_id.values())
    )
    return structural_match, structural_match is not None


def _stringtie_reaches_five_prime(
    stringtie_transcript: Transcript,
    gencode_transcript: Transcript,
) -> bool:
    """True when the StringTie model extends to the GENCODE 5' end.

    On the + strand the 5' end is the lower coordinate (transcript
    start); on the - strand it is the higher coordinate (transcript
    end).
    """
    if (
        stringtie_transcript.start is None
        or stringtie_transcript.end is None
        or gencode_transcript.start is None
        or gencode_transcript.end is None
    ):
        return False

    if stringtie_transcript.strand == "+":
        return stringtie_transcript.start <= gencode_transcript.start
    return stringtie_transcript.end >= gencode_transcript.end


def _apply_gencode_boundaries(
    stringtie_transcript: Transcript,
    gencode_transcript: Transcript,
    force_full_replacement: bool = False,
) -> Transcript:
    """Return a copy with GENCODE-derived CDS/UTR boundaries applied.

    - CDS regions come verbatim from the GENCODE transcript.
    - The outer UTR edge comes from the GENCODE transcript, padded by
      ``UTR_PADDING`` nt on the outside.
    - For an ID-based match, the whole exon structure is only replaced by
      the GENCODE definition when the StringTie assembly does not reach
      the GENCODE 5' end; otherwise the StringTie exon structure is kept
      and only the outer CDS/UTR boundaries are patched in.
    - For a structural (3'-anchored) match, ``force_full_replacement``
      should be True: such a match only guarantees agreement from the
      matched junctions through to the shared 3' end and says nothing
      about UTR length or 5' extent, so the whole exon structure is
      always replaced by the GENCODE definition, ignoring any UTR-side
      difference in the StringTie assembly.
    """
    resolved = deepcopy(stringtie_transcript)

    padded_start = gencode_transcript.start - UTR_PADDING
    padded_end = gencode_transcript.end + UTR_PADDING

    reaches_five_prime = not force_full_replacement and (
        _stringtie_reaches_five_prime(stringtie_transcript, gencode_transcript)
    )

    exons = (
        sorted(stringtie_transcript.exons)
        if reaches_five_prime
        else sorted(gencode_transcript.exons)
    )

    if exons:
        exons[0] = (padded_start, exons[0][1])
        exons[-1] = (exons[-1][0], padded_end)

    resolved.exons = exons
    resolved.cds_regions = deepcopy(gencode_transcript.cds_regions)
    resolved.start = padded_start
    resolved.end = padded_end
    resolved.reference_id = gencode_transcript.transcript_id
    resolved.reference_gene_id = gencode_transcript.gene_id

    return resolved


def _max_gencode_exons_overlapped_by_one_exon(
    stringtie_transcript: Transcript,
    gencode_transcript: Transcript,
) -> int:
    """Return the largest number of GENCODE exons any single StringTie
    exon overlaps.

    A single StringTie exon spanning across several real GENCODE exon
    boundaries indicates the assembly likely merged separate exons
    together (e.g. from unspliced intronic read-through), rather than a
    genuine structural difference.
    """
    gencode_exons = sorted(gencode_transcript.exons)
    max_overlap = 0
    for exon_start, exon_end in stringtie_transcript.exons:
        overlap_count = sum(
            1
            for gencode_start, gencode_end in gencode_exons
            if gencode_start < exon_end and gencode_end > exon_start
        )
        max_overlap = max(max_overlap, overlap_count)
    return max_overlap


def _spans_too_many_gencode_exons(
    stringtie_transcript: Transcript,
    gencode_transcript: Transcript,
) -> bool:
    """True when a single StringTie exon overlaps too many GENCODE exons."""
    return (
        _max_gencode_exons_overlapped_by_one_exon(
            stringtie_transcript, gencode_transcript
        )
        >= MAX_GENCODE_EXONS_SPANNED
    )


def _cds_key(transcript: Transcript) -> tuple[tuple[int, int], ...]:
    """Key used to detect transcripts that share the same CDS content."""
    return tuple(sorted(transcript.cds_regions))


def _span(transcript: Transcript) -> int:
    if transcript.start is None or transcript.end is None:
        return 0
    return transcript.end - transcript.start


def _collapse_cds_duplicates(
    candidates: list[Transcript],
    prefer_gencode: bool = False,
) -> list[Transcript]:
    """Collapse transcripts that share the same CDS into one representative.

    Keep the largest genomic span by default. In the final pass, prefer
    an actual GENCODE model (MANE/canonical/length ordering) when present.
    Preserve member identities and abundance across successive collapses.
    """
    # The final pass also sees added GENCODE representatives. Chromosome and
    # strand are part of identity even when this helper is called directly.
    groups: dict[tuple, list[Transcript]] = {}

    for transcript in candidates:
        key = (transcript.chrom, transcript.strand, _cds_key(transcript))
        groups.setdefault(key, []).append(transcript)

    collapsed: list[Transcript] = []

    for key, group in groups.items():
        if not key[2]:
            # No CDS information to compare on; keep each independently.
            collapsed.extend(group)
            continue

        group.sort(key=lambda t: (-_span(t), t.transcript_id))
        representative = group[0]
        if prefer_gencode:
            representative = _select_representative_gencode(
                [t for t in group if is_gencode_transcript(t)]
            ) or representative

        # A GENCODE model stores its gene in gene_id, whereas the output
        # correspondence table uses reference_gene_id. Preserve that identity
        # when it replaces a StringTie representative.
        representative = deepcopy(representative)
        if is_gencode_transcript(representative):
            representative.reference_gene_id = representative.gene_id
        elif not representative.reference_gene_id:
            genes = {t.reference_gene_id for t in group if t.reference_gene_id}
            if len(genes) == 1:
                representative.reference_gene_id = next(iter(genes))
        representative.cds_equivalent_transcript_ids = sorted({
            tid for t in group for tid in t.cds_equivalent_transcript_ids
        })

        if len(group) > 1:
            representative = deepcopy(representative)
            source_ids, gencode_ids, stringtie_ids = set(), set(), set()
            for member in group:
                source_ids.update(member.source_transcript_ids)
                source_ids.add(member.transcript_id)
                gencode_ids.update(member.gencode_transcript_ids)
                stringtie_ids.update(member.stringtie_transcript_ids)
                if member.reference_id:
                    source_ids.add(member.reference_id)
                    gencode_ids.add(member.reference_id)
                if is_gencode_transcript(member):
                    gencode_ids.add(member.transcript_id)
                elif is_stringtie_transcript(member):
                    stringtie_ids.add(member.transcript_id)
            source_ids.update(gencode_ids)
            source_ids.update(stringtie_ids)
            representative.source_transcript_ids = sorted(source_ids)
            representative.gencode_transcript_ids = sorted(gencode_ids)
            representative.stringtie_transcript_ids = sorted(stringtie_ids)
            representative.tpm = sum(t.tpm for t in group)
            if prefer_gencode and is_gencode_transcript(representative):
                representative.representative_transcript_id = representative.transcript_id
                representative.resolution_type = "gencode_representative_same_cds"

        collapsed.append(representative)

    return collapsed


def _internal_junctions(transcript: Transcript) -> tuple[tuple[int, int], ...]:
    """Return the internal splice junctions (exon boundary pairs).

    A single-exon transcript has no internal junctions and returns an
    empty tuple.
    """
    exons = sorted(transcript.exons)
    return tuple(
        (exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)
    )


def _is_exon_subset(candidate: Transcript, other: Transcript) -> bool:
    """True when ``candidate``'s exon structure is contained in ``other``'s.

    "Contained" means: genomically spanned by ``other``, and, if
    ``candidate`` has internal splice junctions, that its junction chain
    appears as a contiguous run within ``other``'s junction chain. This
    covers exact subsets as well as transcripts that start or end
    partway through ``other`` but otherwise splice the same way.
    """
    if candidate.chrom != other.chrom or candidate.strand != other.strand:
        return False
    if (
        candidate.start is None
        or candidate.end is None
        or other.start is None
        or other.end is None
    ):
        return False
    if not (other.start <= candidate.start and candidate.end <= other.end):
        return False

    candidate_junctions = _internal_junctions(candidate)
    if not candidate_junctions:
        return True

    other_junctions = _internal_junctions(other)
    window = len(candidate_junctions)
    for i in range(len(other_junctions) - window + 1):
        if other_junctions[i : i + window] == candidate_junctions:
            return True
    return False


def _dominates(a: Transcript, b: Transcript) -> bool:
    """True when ``a`` should absorb ``b`` (``b`` is a subset of ``a``).

    Ties (identical exon structure/span) are broken by transcript_id so
    that two mutually-subset transcripts don't both get dropped.
    """
    if not _is_exon_subset(b, a):
        return False
    if _span(a) != _span(b):
        return True
    return a.transcript_id < b.transcript_id


def _collapse_exon_subset_duplicates(
    transcripts: list[Transcript],
) -> list[Transcript]:
    """Collapse GENCODE-unmatched transcripts by exon-structure containment.

    Any transcript whose exon structure is a subset of, or splice-
    compatible but truncated within, another transcript in the list is
    dropped; only the maximal transcript is kept, with the dropped
    transcript_ids recorded on ``source_transcript_ids``.
    """
    dominators: dict[str, Transcript] = {}
    for candidate in transcripts:
        best_dominator = None
        for other in transcripts:
            if other is candidate:
                continue
            if _dominates(other, candidate):
                if best_dominator is None or _span(other) > _span(
                    best_dominator
                ):
                    best_dominator = other
        if best_dominator is not None:
            dominators[candidate.transcript_id] = best_dominator

    collapsed: list[Transcript] = []
    for transcript in transcripts:
        if transcript.transcript_id in dominators:
            continue  # absorbed into a larger transcript

        absorbed_ids = [
            t.transcript_id
            for t in transcripts
            if dominators.get(t.transcript_id) is transcript
        ]
        if absorbed_ids:
            transcript = deepcopy(transcript)
            transcript.source_transcript_ids = (
                transcript.source_transcript_ids + absorbed_ids
            )
        collapsed.append(transcript)

    return collapsed


def _filter_low_tpm_within_gene(
    transcripts: list[Transcript],
) -> list[Transcript]:
    """Drop transcripts at or below LOW_TPM_FRACTION of their gene's max TPM.

    Transcripts are grouped by gene (``gene_id``, falling back to
    ``reference_gene_id`` for GENCODE-boundary-resolved transcripts).
    Within each group, any transcript whose TPM is <= LOW_TPM_FRACTION of
    the group's maximum TPM is dropped as a minor/noise isoform.
    Transcripts with no resolvable gene key are kept as-is.

    A transcript with ``tpm == 0.0`` is treated as having unmeasured TPM
    (the dataclass default) rather than genuinely zero expression, and is
    kept unconditionally, since a real transcript essentially never has
    exactly zero TPM. Such transcripts are also excluded from the
    per-gene max-TPM calculation so a single unmeasured entry can't
    suppress the max used to filter its siblings.
    """
    max_tpm_by_gene: dict[str, float] = {}
    for transcript in transcripts:
        if transcript.tpm == 0.0:
            continue  # treat as unmeasured, not a real zero
        gene_key = _gencode_gene_key(transcript)
        if gene_key is None:
            continue
        max_tpm_by_gene[gene_key] = max(
            max_tpm_by_gene.get(gene_key, 0.0), transcript.tpm
        )

    kept: list[Transcript] = []
    for transcript in transcripts:
        if transcript.tpm == 0.0:
            kept.append(transcript)  # unmeasured TPM: keep unconditionally
            continue

        gene_key = _gencode_gene_key(transcript)
        if gene_key is None:
            kept.append(transcript)
            continue

        max_tpm = max_tpm_by_gene[gene_key]
        if max_tpm <= 0:
            kept.append(transcript)
            continue

        if transcript.tpm <= LOW_TPM_FRACTION * max_tpm:
            continue

        kept.append(transcript)

    return kept


def _gencode_gene_key(transcript: Transcript) -> str | None:
    """Return the gene identity key used to relate a transcript to GENCODE.

    GENCODE transcripts are keyed by their own ``gene_id`` (e.g.
    ``ENSG...``). StringTie transcripts are keyed by ``reference_gene_id``
    (the GENCODE gene they were matched against via the GTF's
    ``ref_gene_id`` attribute), since their own ``gene_id`` is normally
    an internal StringTie grouping ID (e.g. ``MSTRG.123``) that has no
    relation to any GENCODE gene ID. Falls back to ``gene_id`` when no
    ``reference_gene_id`` is set (e.g. a StringTie transcript at a novel
    locus with no GENCODE overlap at all).
    """
    if is_gencode_transcript(transcript):
        return transcript.gene_id
    return transcript.reference_gene_id or transcript.gene_id


def _select_representative_gencode(
    candidates: list[Transcript],
) -> Transcript | None:
    """Pick the single most representative GENCODE transcript for a gene.

    Preference order: MANE_Select tag, then Ensembl_canonical tag, then
    longest CDS, then longest genomic span, then lowest transcript_id as
    a final tiebreak. (If ``tags`` is absent/empty on every candidate ?
    e.g. the GFF3 loader hasn't been updated to populate it yet ?
    selection reduces to longest CDS, then span, then id.)
    """
    if not candidates:
        return None

    def sort_key(transcript: Transcript):
        tags = set(getattr(transcript, "tags", None) or [])
        cds_length = sum(end - start for start, end in transcript.cds_regions)
        return (
            "MANE_Select" not in tags,
            "Ensembl_canonical" not in tags,
            -cds_length,
            -_span(transcript),
            transcript.transcript_id,
        )

    return sorted(candidates, key=sort_key)[0]


def _transcripts_overlap(
    first: Transcript,
    second: Transcript,
) -> bool:
    """Return True when two transcripts overlap on the same strand."""

    if first.chrom != second.chrom or first.strand != second.strand:
        return False

    if (
        first.start is None
        or first.end is None
        or second.start is None
        or second.end is None
    ):
        return False

    return (
        first.start < second.end
        and second.start < first.end
    )


def _add_gencode_representative_when_uncovered(
    resolved_transcripts: list[Transcript],
    overlap_group: list[Transcript],
    overlap_group_id: str,
    added_gencode_ids: set[str],
) -> list[Transcript]:
    """
    Add the representative GENCODE transcript when it is not represented.

    The decision is made per representative transcript ID, not merely per
    gene. Therefore, when StringTie supports another isoform of the same
    gene but does not support the selected representative transcript, the
    representative GENCODE transcript is still added.

    A representative is added only when at least one StringTie transcript
    overlaps that representative on the same chromosome and strand.
    """

    gencode_by_gene: dict[str, list[Transcript]] = {}

    stringtie_transcripts = [
        transcript
        for transcript in overlap_group
        if is_stringtie_transcript(transcript)
    ]

    for transcript in overlap_group:
        if not is_gencode_transcript(transcript):
            continue

        gene_key = _gencode_gene_key(transcript)

        if gene_key is None:
            continue

        gencode_by_gene.setdefault(
            gene_key,
            [],
        ).append(
            transcript
        )

    represented_gencode_ids = {
        transcript.reference_id
        for transcript in resolved_transcripts
        if transcript.reference_id
    }

    represented_gencode_ids.update(
        transcript.transcript_id
        for transcript in resolved_transcripts
        if is_gencode_transcript(transcript)
    )

    additions: list[Transcript] = []

    for _, candidates in gencode_by_gene.items():
        representative = _select_representative_gencode(
            candidates
        )

        if representative is None:
            continue

        representative_id = representative.transcript_id

        if representative_id in represented_gencode_ids:
            continue

        if representative_id in added_gencode_ids:
            continue

        has_overlapping_stringtie = any(
            _transcripts_overlap(
                stringtie_transcript,
                representative,
            )
            for stringtie_transcript in stringtie_transcripts
        )

        if not has_overlapping_stringtie:
            continue

        added = deepcopy(representative)
        added.reference_gene_id = representative.gene_id
        added.overlap_group_id = overlap_group_id
        added.resolution_type = "gencode_representative_uncovered"
        added.source_transcript_ids = [representative_id]

        additions.append(added)
        added_gencode_ids.add(representative_id)

    return resolved_transcripts + additions


def resolve_one_overlap_group(
    overlap_group: list[Transcript],
    overlap_group_id: str,
    added_gencode_ids: set[str] | None = None,
) -> list[Transcript]:
    """
    Return the StringTie-supported transcripts kept for one overlap group.

    Resolution rules:
    1. GENCODE transcripts in the group are used only as an annotation
       reference and are never returned directly.
    2. StringTie transcripts containing any micro-exon (length <=
       MIN_EXON_LENGTH nt) are discarded outright, as if they never
       existed.
    3. Each StringTie transcript is matched to a GENCODE transcript by
       ``reference_id``/``gencode_transcript_ids`` if available, otherwise
       by 3'-anchored intron-chain matching (its splice junctions, read
       5' -> 3', form a subset of the GENCODE transcript's junctions that
       reaches the GENCODE 3' end). CDS regions and the outer UTR
       boundary (padded by ``UTR_PADDING`` nt) always come from the
       GENCODE transcript. For an ID-based match, the StringTie exon
       structure is kept as-is unless it fails to reach the GENCODE 5'
       end, in which case it is replaced wholesale by the GENCODE exon
       structure. For a structural match, the exon structure is always
       replaced wholesale by the GENCODE structure, since such a match
       only guarantees 3'-side agreement.
    4. Among GENCODE-matched transcripts, those sharing the same resolved
       CDS are treated as duplicates; only the one with the largest
       genomic span is kept.
    5. Among transcripts with no GENCODE match, exon structure subsets
       (or splice-compatible but truncated transcripts) are collapsed
       into the largest enclosing transcript.
    6. GENCODE-matched StringTie transcripts where a single exon overlaps
       MAX_GENCODE_EXONS_SPANNED or more real GENCODE exons are dropped
       (likely unspliced/merged-exon assembly artifacts).
    7. Within each gene group, transcripts at or below LOW_TPM_FRACTION of
       the group's maximum TPM are dropped as minor isoforms.
    8. Add the gene's selected GENCODE representative when its ID is not
       represented and an overlapping StringTie transcript exists.
       ``added_gencode_ids`` is shared across calls (by the caller,
       across all overlap groups) so the same GENCODE transcript is
       never added more than once even if a gene's loci are split across
       multiple overlap groups.
    9. Collapse identical CDS again after the additions, preferring the
       GENCODE representative and preserving all absorbed identities.
    """
    if added_gencode_ids is None:
        added_gencode_ids = set()

    gencode_by_id = {
        transcript.transcript_id: transcript
        for transcript in overlap_group
        if is_gencode_transcript(transcript)
    }

    matched_candidates: list[Transcript] = []
    unmatched_candidates: list[Transcript] = []

    for transcript in overlap_group:
        if not is_stringtie_transcript(transcript):
            continue
        if not _has_no_micro_exons(transcript):
            continue

        matched_gencode, _is_structural_match = _find_matched_gencode(
            transcript, gencode_by_id
        )

        if matched_gencode is not None and _spans_too_many_gencode_exons(
            transcript, matched_gencode
        ):
            continue  # one StringTie exon merges too many real GENCODE exons

        resolved = deepcopy(transcript)
        resolved.overlap_group_id = overlap_group_id
        resolved.resolution_type = "stringtie_supported"
        resolved.source_transcript_ids = sorted(set(
            transcript.source_transcript_ids + [transcript.transcript_id]
        ))

        if matched_gencode is not None:
            resolved = _apply_gencode_boundaries(
                resolved,
                matched_gencode,
                force_full_replacement=_is_structural_match,
            )
            resolved.overlap_group_id = overlap_group_id
            resolved.resolution_type = "stringtie_supported"
            resolved.source_transcript_ids = sorted(set(
                transcript.source_transcript_ids + [transcript.transcript_id]
            ))
            matched_candidates.append(resolved)
        else:
            unmatched_candidates.append(resolved)

    resolved_transcripts = _collapse_cds_duplicates(
        matched_candidates
    ) + _collapse_exon_subset_duplicates(unmatched_candidates)

    resolved_transcripts = _filter_low_tpm_within_gene(resolved_transcripts)

    resolved_transcripts = _add_gencode_representative_when_uncovered(
        resolved_transcripts, overlap_group, overlap_group_id, added_gencode_ids
    )

    # Fallback is ID-based and can add a representative whose CDS already
    # exists under a different ENST ID. Collapse again before assigning final
    # reference IDs, retaining the GENCODE model and all absorbed identities.
    resolved_transcripts = _collapse_cds_duplicates(
        resolved_transcripts, prefer_gencode=True,
    )

    # Preserve annotation equivalence separately from the IDs actually
    # observed/absorbed. UTR-only alternatives can be absent from the surviving
    # StringTie candidates yet name the same translated unit in simulation.
    by_cds: dict[tuple, list[Transcript]] = {}
    for annotation in gencode_by_id.values():
        if annotation.cds_regions:
            key = (annotation.chrom, annotation.strand, _cds_key(annotation))
            by_cds.setdefault(key, []).append(annotation)
    for transcript in resolved_transcripts:
        key = (transcript.chrom, transcript.strand, _cds_key(transcript))
        aliases = by_cds.get(key, []) if transcript.cds_regions else []
        if transcript.reference_gene_id:
            aliases = [a for a in aliases if a.gene_id == transcript.reference_gene_id]
        genes = {a.gene_id for a in aliases if a.gene_id}
        if not transcript.reference_gene_id and len(genes) == 1:
            transcript.reference_gene_id = next(iter(genes))
        # Do not infer an alias across ambiguous gene annotations.
        if len(genes) == 1:
            transcript.cds_equivalent_transcript_ids = sorted({
                a.transcript_id for a in aliases
            })

    return sorted(
        resolved_transcripts,
        key=lambda transcript: (
            transcript.start if transcript.start is not None else -1,
            transcript.end if transcript.end is not None else -1,
            transcript.transcript_id,
        ),
    )
