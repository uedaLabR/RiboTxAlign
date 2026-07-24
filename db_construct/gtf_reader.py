from __future__ import annotations

from urllib.parse import unquote

from db_construct.data_class import Transcript


def parse_gtf_attributes(attribute_text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in attribute_text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(None, 1)
        if len(parts) != 2:
            continue
        key, value = parts
        attributes[key] = value.strip().strip('"')
    return attributes


def parse_gff3_attributes(attribute_text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in attribute_text.strip().split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[key] = unquote(value)
    return attributes


def remove_prefix(value: str | None, prefix: str) -> str | None:
    if value is None:
        return None
    if value.startswith(prefix):
        return value[len(prefix):]
    return value


def safe_float(value: str | None, default: float = 0.0) -> float:
    if value is None:
        return default
    try:
        return float(value)
    except ValueError:
        return default


def _get_or_create_transcript(
    transcripts: dict[str, Transcript],
    transcript_id: str,
    gene_id: str | None,
    chrom: str,
    strand: str,
    source: str,
) -> Transcript:
    if transcript_id not in transcripts:
        transcripts[transcript_id] = Transcript(
            transcript_id=transcript_id,
            gene_id=gene_id,
            chrom=chrom,
            strand=strand,
            source=source,
        )
    return transcripts[transcript_id]


def _finalize_transcripts(transcripts: dict[str, Transcript]) -> None:
    for transcript in transcripts.values():
        transcript.exons = sorted(set(transcript.exons))
        transcript.cds_regions = sorted(set(transcript.cds_regions))
        if transcript.exons:
            transcript.start = min(start for start, _ in transcript.exons)
            transcript.end = max(end for _, end in transcript.exons)


def load_gencode_gff3(gff3_path: str) -> dict[str, Transcript]:
    transcripts: dict[str, Transcript] = {}
    with open(gff3_path, encoding="utf-8") as input_file:
        for line in input_file:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                continue
            chrom = columns[0]
            feature = columns[2]
            start = int(columns[3]) - 1
            end = int(columns[4])
            strand = columns[6]
            attributes = parse_gff3_attributes(columns[8])

            if feature in {"transcript", "mRNA"}:
                transcript_id = remove_prefix(attributes.get("ID"), "transcript:")
                gene_id = remove_prefix(attributes.get("Parent"), "gene:")
                if transcript_id is None:
                    continue
                transcript = _get_or_create_transcript(
                    transcripts, transcript_id, gene_id, chrom, strand, "gencode"
                )
                transcript.gene_id = gene_id
                transcript.chrom = chrom
                transcript.strand = strand
                transcript.start = start
                transcript.end = end
                transcript.gene_name = attributes.get("gene_name")

            elif feature in {"exon", "CDS"}:
                parent_text = attributes.get("Parent")
                if parent_text is None:
                    continue
                for parent_id in parent_text.split(","):
                    transcript_id = remove_prefix(parent_id, "transcript:")
                    if transcript_id is None:
                        continue
                    transcript = _get_or_create_transcript(
                        transcripts,
                        transcript_id,
                        remove_prefix(attributes.get("gene_id"), "gene:"),
                        chrom,
                        strand,
                        "gencode",
                    )
                    if feature == "exon":
                        transcript.exons.append((start, end))
                    else:
                        transcript.cds_regions.append((start, end))
                    if transcript.gene_name is None:
                        transcript.gene_name = attributes.get("gene_name")

    _finalize_transcripts(transcripts)
    return transcripts


def load_stringtie_transcripts(gtf_path: str) -> dict[str, Transcript]:
    transcripts: dict[str, Transcript] = {}
    with open(gtf_path, encoding="utf-8") as input_file:
        for line in input_file:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                continue
            chrom = columns[0]
            feature = columns[2]
            start = int(columns[3]) - 1
            end = int(columns[4])
            strand = columns[6]
            attributes = parse_gtf_attributes(columns[8])
            transcript_id = attributes.get("transcript_id")
            if transcript_id is None:
                continue

            transcript = _get_or_create_transcript(
                transcripts,
                transcript_id,
                attributes.get("gene_id"),
                chrom,
                strand,
                "stringtie",
            )
            transcript.gene_id = attributes.get("gene_id") or transcript.gene_id
            transcript.chrom = chrom
            transcript.strand = strand
            transcript.reference_id = attributes.get("reference_id") or transcript.reference_id
            transcript.reference_gene_id = attributes.get("ref_gene_id") or transcript.reference_gene_id
            transcript.gene_name = (
                attributes.get("ref_gene_name")
                or attributes.get("gene_name")
                or transcript.gene_name
            )

            if feature == "transcript":
                transcript.start = start
                transcript.end = end
                transcript.tpm = safe_float(attributes.get("TPM"))
            elif feature == "exon":
                transcript.exons.append((start, end))

    _finalize_transcripts(transcripts)
    return transcripts
