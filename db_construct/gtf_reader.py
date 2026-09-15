from __future__ import annotations

from urllib.parse import unquote

from db_construct.data_class import Transcript


def parse_gtf_attributes(
    attribute_text: str,
) -> dict[str, str]:
    """Parse standard GTF attributes."""

    attributes: dict[str, str] = {}

    for item in attribute_text.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        parts = item.split(None, 1)

        if len(parts) != 2:
            continue

        key, value = parts

        attributes[key] = (
            value.strip().strip('"')
        )

    return attributes


def parse_gff3_attributes(
    attribute_text: str,
) -> dict[str, str]:
    """Parse standard GFF3 key=value attributes."""

    attributes: dict[str, str] = {}

    for item in attribute_text.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" not in item:
            continue

        key, value = item.split("=", 1)

        attributes[key] = unquote(value)

    return attributes


def detect_attribute_format(
    attribute_text: str,
) -> str:
    """
    Detect whether one annotation attribute field is GTF or GFF3 style.

    GTF:
        transcript_id "ENST...";
        gene_id "ENSG...";

    GFF3:
        ID=ENST...;
        Parent=ENST...;
        geneID=ENSG...;
    """

    for item in attribute_text.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" in item:
            return "gff3"

        return "gtf"

    return "gtf"


def remove_prefix(
    value: str | None,
    prefix: str,
) -> str | None:
    """Remove a prefix when present."""

    if value is None:
        return None

    if value.startswith(prefix):
        return value[len(prefix):]

    return value


def safe_float(
    value: str | None,
    default: float = 0.0,
) -> float:
    """Convert text to float while tolerating missing values."""

    if value is None:
        return default

    try:
        return float(value)

    except ValueError:
        return default


def parse_tags(
    tag_text: str | None,
) -> list[str]:
    """
    Parse a GFF3 'tag' attribute into individual tag names.

    Example:
        tag=basic,Ensembl_canonical,MANE_Select
    """

    if not tag_text:
        return []

    return [
        tag
        for tag in tag_text.split(",")
        if tag
    ]


def _get_or_create_transcript(
    transcripts: dict[str, Transcript],
    transcript_id: str,
    gene_id: str | None,
    chrom: str,
    strand: str,
    source: str,
) -> Transcript:
    """Get an existing transcript or create a new one."""

    if transcript_id not in transcripts:
        transcripts[transcript_id] = Transcript(
            transcript_id=transcript_id,
            gene_id=gene_id,
            chrom=chrom,
            strand=strand,
            source=source,
        )

    return transcripts[transcript_id]


def _finalize_transcripts(
    transcripts: dict[str, Transcript],
) -> None:
    """Sort coordinates and define transcript genomic boundaries."""

    for transcript in transcripts.values():

        transcript.exons = sorted(
            set(transcript.exons)
        )

        transcript.cds_regions = sorted(
            set(transcript.cds_regions)
        )

        if transcript.exons:
            transcript.start = min(
                start
                for start, _ in transcript.exons
            )

            transcript.end = max(
                end
                for _, end in transcript.exons
            )


def load_gencode_gff3(
    gff3_path: str,
) -> dict[str, Transcript]:
    """Load GENCODE transcripts from GFF3."""

    transcripts: dict[
        str,
        Transcript,
    ] = {}

    with open(
        gff3_path,
        encoding="utf-8",
    ) as input_file:

        for line in input_file:

            if line.startswith("#"):
                continue

            columns = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            if len(columns) != 9:
                continue

            chrom = columns[0]
            feature = columns[2]

            start = int(
                columns[3]
            ) - 1

            end = int(
                columns[4]
            )

            strand = columns[6]

            attributes = parse_gff3_attributes(
                columns[8]
            )

            if feature in {
                "transcript",
                "mRNA",
            }:

                transcript_id = remove_prefix(
                    attributes.get("ID"),
                    "transcript:",
                )

                gene_id = remove_prefix(
                    attributes.get("Parent"),
                    "gene:",
                )

                if transcript_id is None:
                    continue

                transcript = (
                    _get_or_create_transcript(
                        transcripts,
                        transcript_id,
                        gene_id,
                        chrom,
                        strand,
                        "gencode",
                    )
                )

                transcript.gene_id = gene_id
                transcript.chrom = chrom
                transcript.strand = strand
                transcript.start = start
                transcript.end = end

                transcript.gene_name = (
                    attributes.get(
                        "gene_name"
                    )
                )

                transcript.tags = parse_tags(
                    attributes.get(
                        "tag"
                    )
                )

            elif feature in {
                "exon",
                "CDS",
            }:

                parent_text = attributes.get(
                    "Parent"
                )

                if parent_text is None:
                    continue

                for parent_id in parent_text.split(
                    ","
                ):

                    transcript_id = remove_prefix(
                        parent_id,
                        "transcript:",
                    )

                    if transcript_id is None:
                        continue

                    transcript = (
                        _get_or_create_transcript(
                            transcripts,
                            transcript_id,
                            remove_prefix(
                                attributes.get(
                                    "gene_id"
                                ),
                                "gene:",
                            ),
                            chrom,
                            strand,
                            "gencode",
                        )
                    )

                    if feature == "exon":
                        transcript.exons.append(
                            (
                                start,
                                end,
                            )
                        )

                    else:
                        transcript.cds_regions.append(
                            (
                                start,
                                end,
                            )
                        )

                    if transcript.gene_name is None:
                        transcript.gene_name = (
                            attributes.get(
                                "gene_name"
                            )
                        )

                    if not transcript.tags:
                        parsed_tags = parse_tags(
                            attributes.get(
                                "tag"
                            )
                        )

                        if parsed_tags:
                            transcript.tags = (
                                parsed_tags
                            )

    _finalize_transcripts(
        transcripts
    )

    return transcripts


def _get_stringtie_transcript_id(
    feature: str,
    attributes: dict[str, str],
    attribute_format: str,
) -> str | None:
    """
    Obtain transcript ID from either StringTie GTF or GFF3-style output.
    """

    if attribute_format == "gtf":
        return attributes.get(
            "transcript_id"
        )

    if feature in {
        "transcript",
        "mRNA",
    }:
        return remove_prefix(
            attributes.get("ID"),
            "transcript:",
        )

    parent_text = attributes.get(
        "Parent"
    )

    if not parent_text:
        return None

    # A StringTie exon should normally have one transcript parent.
    parent_id = parent_text.split(
        ",",
        1,
    )[0]

    return remove_prefix(
        parent_id,
        "transcript:",
    )


def _get_stringtie_gene_id(
    attributes: dict[str, str],
) -> str | None:
    """
    Obtain gene ID from either conventional GTF or GFF3-style
    StringTie output.
    """

    return (
        attributes.get("gene_id")
        or attributes.get("geneID")
        or attributes.get("ref_gene_id")
    )


def _get_stringtie_expression_value(
    attributes: dict[str, str],
) -> float:
    """
    Return the expression value used internally as transcript.tpm.

    TPM is preferred when available. If TPM is unavailable, coverage is
    used as a fallback. RiboTxAlign currently only requires this value
    for expression filtering, so coverage can be used for StringTie
    outputs that do not report TPM.
    """

    tpm = attributes.get(
        "TPM"
    )

    if tpm is not None:
        return safe_float(
            tpm
        )

    coverage = (
        attributes.get("coverage")
        or attributes.get("cov")
    )

    return safe_float(
        coverage
    )


def load_stringtie_transcripts(
    gtf_path: str,
) -> dict[str, Transcript]:
    """
    Load transcripts from StringTie output.

    Both standard StringTie GTF attributes and GFF3-style attributes
    are supported.

    Examples
    --------
    Standard GTF:
        transcript_id "MSTRG.1.1";
        gene_id "MSTRG.1";
        TPM "5.23";

    GFF3-style output:
        ID=ENST00000623083.4;
        geneID=ENSG00000279457.4;
        gene_name=WASH9P;
        coverage=191.34

    Transcripts with an unknown strand are excluded.
    Only '+' and '-' strands are accepted.

    The Transcript.source field always records "stringtie" because
    these records originate from the StringTie output file, regardless
    of whether the second GTF/GFF column says StringTie, ENSEMBL,
    HAVANA, or another annotation source.

    If TPM is unavailable, coverage is stored in transcript.tpm as a
    fallback because downstream RiboTxAlign processing currently uses
    this value only for expression filtering.
    """

    transcripts: dict[
        str,
        Transcript,
    ] = {}

    total_transcript_count = 0
    excluded_transcript_count = 0

    gtf_attribute_count = 0
    gff3_attribute_count = 0

    expression_from_tpm_count = 0
    expression_from_coverage_count = 0
    expression_missing_count = 0

    with open(
        gtf_path,
        encoding="utf-8",
    ) as input_file:

        for line in input_file:

            if line.startswith("#"):
                continue

            columns = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            if len(columns) != 9:
                continue

            chrom = columns[0]
            feature = columns[2]

            start = int(
                columns[3]
            ) - 1

            end = int(
                columns[4]
            )

            strand = columns[6]

            attribute_text = columns[8]

            attribute_format = (
                detect_attribute_format(
                    attribute_text
                )
            )

            if attribute_format == "gff3":
                attributes = (
                    parse_gff3_attributes(
                        attribute_text
                    )
                )

                gff3_attribute_count += 1

            else:
                attributes = (
                    parse_gtf_attributes(
                        attribute_text
                    )
                )

                gtf_attribute_count += 1

            transcript_id = (
                _get_stringtie_transcript_id(
                    feature=feature,
                    attributes=attributes,
                    attribute_format=(
                        attribute_format
                    ),
                )
            )

            if transcript_id is None:
                continue

            if feature in {
                "transcript",
                "mRNA",
            }:
                total_transcript_count += 1

                if strand not in {
                    "+",
                    "-",
                }:
                    excluded_transcript_count += 1
                    continue

            # Exclude exon/CDS records with an unknown strand as well.
            if strand not in {
                "+",
                "-",
            }:
                continue

            gene_id = (
                _get_stringtie_gene_id(
                    attributes
                )
            )

            transcript = (
                _get_or_create_transcript(
                    transcripts,
                    transcript_id,
                    gene_id,
                    chrom,
                    strand,
                    "stringtie",
                )
            )

            # Always identify records loaded from this file as
            # StringTie-derived. The original annotation source
            # (ENSEMBL/HAVANA/etc.) is not used as Transcript.source.
            transcript.source = (
                "stringtie"
            )

            transcript.gene_id = (
                gene_id
                or transcript.gene_id
            )

            transcript.chrom = chrom
            transcript.strand = strand

            transcript.reference_id = (
                attributes.get(
                    "reference_id"
                )
                or attributes.get(
                    "reference_transcript_id"
                )
                or transcript.reference_id
            )

            transcript.reference_gene_id = (
                attributes.get(
                    "ref_gene_id"
                )
                or attributes.get(
                    "geneID"
                )
                or transcript.reference_gene_id
            )

            transcript.gene_name = (
                attributes.get(
                    "ref_gene_name"
                )
                or attributes.get(
                    "gene_name"
                )
                or transcript.gene_name
            )

            if feature in {
                "transcript",
                "mRNA",
            }:

                transcript.start = start
                transcript.end = end

                if attributes.get(
                    "TPM"
                ) is not None:

                    transcript.tpm = (
                        _get_stringtie_expression_value(
                            attributes
                        )
                    )

                    expression_from_tpm_count += 1

                elif (
                    attributes.get(
                        "coverage"
                    ) is not None
                    or attributes.get(
                        "cov"
                    ) is not None
                ):

                    transcript.tpm = (
                        _get_stringtie_expression_value(
                            attributes
                        )
                    )

                    expression_from_coverage_count += 1

                else:
                    transcript.tpm = 0.0
                    expression_missing_count += 1

            elif feature == "exon":

                transcript.exons.append(
                    (
                        start,
                        end,
                    )
                )

            elif feature == "CDS":

                transcript.cds_regions.append(
                    (
                        start,
                        end,
                    )
                )

    _finalize_transcripts(
        transcripts
    )

    excluded_percentage = (
        100.0
        * excluded_transcript_count
        / total_transcript_count
        if total_transcript_count > 0
        else 0.0
    )

    print(
        "StringTie transcripts loaded:",
        len(transcripts),
    )

    print(
        "StringTie transcript records:",
        total_transcript_count,
    )

    print(
        "Unknown-strand transcripts excluded:",
        excluded_transcript_count,
        f"({excluded_percentage:.4f}%)",
    )

    print(
        "StringTie attribute formats:",
        "GTF=",
        gtf_attribute_count,
        "GFF3=",
        gff3_attribute_count,
    )

    print(
        "StringTie expression values:",
        "TPM=",
        expression_from_tpm_count,
        "coverage=",
        expression_from_coverage_count,
        "missing=",
        expression_missing_count,
    )

    return transcripts