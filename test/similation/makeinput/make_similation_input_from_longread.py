#!/usr/bin/env python3

from collections import defaultdict
from pathlib import Path

import pysam


# ============================================================
# Configuration
# ============================================================

BASE = Path(
    "/mnt/share/project/riboML/validation"
)

BASE2 = Path(
    "/mnt/share/project/riboML/validation_expression"
)

OUTPUT_BASE = (
    BASE2
    / "similation"
)

MIN_RELATIVE_TPM = 0.03


primary_chroms = {
    f"chr{i}" for i in range(1, 23)
} | {
    "chrX",
    "chrY",
    "chrM",
}

DATASETS = {
    "HEK293T": {
        "stringtie_gtf": (
            BASE
            / "HEK293T"
            / "longread"
            / "mapping"
            / "HEK293T.longread.stringtie.gtf"
        ),
        "gencode_gff3": (
            BASE
            / "reference"
            / "gencode.v44.annotation.gff3"
        ),
        "genome_fasta": (
            BASE
            / "reference"
            / "hg38.fa"
        ),
    },

    "mouse_brain": {
        "stringtie_gtf": (
            BASE
            / "mouse_brain"
            / "longread"
            / "mapping"
            / "mouse_brain.longread.stringtie.gtf"
        ),
        "gencode_gff3": (
            BASE
            / "reference"
            / "gencode.vM25.annotation.gff3"
        ),
        "genome_fasta": (
            BASE
            / "reference"
            / "mm10.fa"
        ),
    },
}


# ============================================================
# Utility
# ============================================================

def parse_gtf_attributes(text):
    """Parse GTF attributes."""

    attrs = {}

    for item in text.strip().strip(";").split(";"):

        item = item.strip()

        if not item:
            continue

        parts = item.split(
            None,
            1,
        )

        if len(parts) != 2:
            continue

        key = parts[0]
        value = parts[1].strip().strip('"')

        attrs[key] = value

    return attrs


def parse_gff3_attributes(text):
    """Parse GFF3 attributes."""

    attrs = {}

    for item in text.strip().split(";"):

        if "=" not in item:
            continue

        key, value = item.split(
            "=",
            1,
        )

        attrs[key] = value

    return attrs


def strip_version(transcript_id):
    """Remove a GENCODE transcript version suffix."""

    return transcript_id.split(".")[0]


def reverse_complement(seq):
    """Return reverse-complement sequence."""

    complement = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan",
    )

    return seq.translate(
        complement
    )[::-1]


# ============================================================
# Read StringTie GTF
# ============================================================

def load_stringtie_gtf(gtf_path):
    """Read StringTie transcript structures, TPM, and reference IDs."""

    transcripts = {}

    with open(gtf_path) as input_file:

        for line in input_file:

            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")

            if len(cols) != 9:
                continue

            chrom = cols[0]
            feature = cols[2]
            start = int(cols[3]) - 1
            end = int(cols[4])
            strand = cols[6]

            attrs = parse_gtf_attributes(
                cols[8]
            )

            transcript_id = attrs.get(
                "transcript_id"
            )

            if not transcript_id:
                continue

            gene_id = (
                attrs.get("ref_gene_id")
                or attrs.get("gene_id")
                or ""
            )

            gene_name = (
                attrs.get("ref_gene_name")
                or attrs.get("gene_name")
                or ""
            )

            if transcript_id not in transcripts:

                try:
                    tpm = float(
                        attrs.get(
                            "TPM",
                            "0",
                        )
                    )
                except ValueError:
                    tpm = 0.0

                reference_id = (
                    attrs.get("reference_id")
                    or attrs.get("ref_id")
                    or attrs.get("reference_transcript_id")
                    or ""
                )

                if (
                    not reference_id
                    and (
                        transcript_id.startswith("ENST")
                        or transcript_id.startswith("ENSMUST")
                    )
                ):
                    reference_id = transcript_id

                transcripts[
                    transcript_id
                ] = {
                    "transcript_id": transcript_id,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "chrom": chrom,
                    "strand": strand,
                    "tpm": tpm,
                    "reference_id": reference_id,
                    "exons": [],
                }

            if feature == "transcript":

                if "TPM" in attrs:

                    try:
                        transcripts[
                            transcript_id
                        ]["tpm"] = float(
                            attrs["TPM"]
                        )
                    except ValueError:
                        pass

                reference_id = (
                    attrs.get("reference_id")
                    or attrs.get("ref_id")
                    or attrs.get("reference_transcript_id")
                )

                if reference_id:
                    transcripts[
                        transcript_id
                    ]["reference_id"] = reference_id

                reference_gene_id = attrs.get(
                    "ref_gene_id"
                )

                if reference_gene_id:
                    transcripts[
                        transcript_id
                    ]["gene_id"] = reference_gene_id

                reference_gene_name = attrs.get(
                    "ref_gene_name"
                )

                if reference_gene_name:
                    transcripts[
                        transcript_id
                    ]["gene_name"] = reference_gene_name

            elif feature == "exon":

                transcripts[
                    transcript_id
                ]["exons"].append(
                    (
                        start,
                        end,
                    )
                )

    for tx in transcripts.values():

        tx["exons"].sort(
            key=lambda x: x[0],
            reverse=(
                tx["strand"] == "-"
            ),
        )

    return transcripts


# ============================================================
# Filter by relative TPM
# ============================================================

def filter_by_relative_tpm(
    transcripts,
    relative_cutoff,
):
    """
    Keep transcripts with TPM >= relative_cutoff
    times the maximum TPM of the same gene.
    """

    gene_to_transcripts = defaultdict(
        list
    )

    for tx in transcripts.values():

        gene_id = tx["gene_id"]

        if not gene_id:
            gene_id = (
                "__NO_GENE__:"
                + tx["transcript_id"]
            )

        gene_to_transcripts[
            gene_id
        ].append(tx)

    selected = {}

    for tx_list in gene_to_transcripts.values():

        max_tpm = max(
            tx["tpm"]
            for tx in tx_list
        )

        if max_tpm <= 0.0:
            continue

        threshold = (
            max_tpm
            * relative_cutoff
        )

        for tx in tx_list:

            if tx["tpm"] >= threshold:

                selected[
                    tx["transcript_id"]
                ] = tx

    return selected


# ============================================================
# Read GENCODE protein-coding transcript/CDS information
# ============================================================

def load_gencode_cds(
    gff3_path,
):
    """Load GENCODE transcript metadata and CDS intervals."""

    transcript_info = {}
    protein_coding = set()

    cds_intervals = defaultdict(
        list
    )

    with open(gff3_path) as input_file:

        for line in input_file:

            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")

            if len(cols) != 9:
                continue

            chrom = cols[0]
            feature = cols[2]
            start = int(cols[3]) - 1
            end = int(cols[4])
            strand = cols[6]

            attrs = parse_gff3_attributes(
                cols[8]
            )

            if feature in {
                "transcript",
                "mRNA",
            }:

                transcript_id = (
                    attrs.get("transcript_id")
                    or attrs.get("ID", "")
                )

                transcript_id = transcript_id.replace(
                    "transcript:",
                    "",
                )

                if not transcript_id:
                    continue

                tid_base = strip_version(
                    transcript_id
                )

                gene_id = attrs.get(
                    "gene_id",
                    "",
                )

                gene_name = attrs.get(
                    "gene_name",
                    "",
                )

                transcript_type = (
                    attrs.get("transcript_type")
                    or attrs.get("transcript_biotype")
                    or ""
                )

                transcript_info[
                    tid_base
                ] = {
                    "transcript_id": transcript_id,
                    "chrom": chrom,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "transcript_type": transcript_type,
                }

                if transcript_type == "protein_coding":

                    protein_coding.add(
                        tid_base
                    )

            elif feature == "CDS":

                parent = attrs.get(
                    "Parent",
                    "",
                )

                if not parent:
                    continue

                for parent_id in parent.split(","):

                    parent_id = parent_id.replace(
                        "transcript:",
                        "",
                    )

                    tid_base = strip_version(
                        parent_id
                    )

                    cds_intervals[
                        tid_base
                    ].append(
                        (
                            start,
                            end,
                        )
                    )

    for transcript_id in cds_intervals:

        cds_intervals[
            transcript_id
        ] = sorted(
            set(
                cds_intervals[
                    transcript_id
                ]
            )
        )

    return (
        transcript_info,
        protein_coding,
        cds_intervals,
    )


# ============================================================
# Build transcript sequence
# ============================================================

def build_transcript_sequence(
    genome,
    tx,
):
    """Construct transcript sequence from genomic exons."""

    sequence_parts = []

    for start, end in tx["exons"]:

        seq = genome.fetch(
            tx["chrom"],
            start,
            end,
        )

        if tx["strand"] == "-":

            seq = reverse_complement(
                seq
            )

        sequence_parts.append(
            seq
        )

    return "".join(
        sequence_parts
    ).upper()


# ============================================================
# Convert genomic CDS to transcript coordinates
# ============================================================

def genomic_cds_to_transcript(
    tx,
    cds_intervals,
):
    """
    Convert GENCODE genomic CDS intervals to StringTie transcript coordinates.

    Returns:
        cds_start, cds_stop

    Coordinates are 0-based transcript coordinates.
    cds_stop is exclusive.

    All annotated CDS bases must be covered by the StringTie transcript,
    and the projected CDS must be contiguous in transcript coordinates.
    """

    exon_positions = []

    transcript_pos = 0

    for exon_start, exon_end in tx["exons"]:

        exon_length = (
            exon_end
            - exon_start
        )

        exon_positions.append(
            {
                "genome_start": exon_start,
                "genome_end": exon_end,
                "tx_start": transcript_pos,
                "tx_end": (
                    transcript_pos
                    + exon_length
                ),
            }
        )

        transcript_pos += exon_length

    projected_intervals = []
    covered_cds_bases = 0

    expected_cds_bases = sum(
        cds_end - cds_start
        for cds_start, cds_end
        in cds_intervals
    )

    for cds_start, cds_end in cds_intervals:

        for exon in exon_positions:

            overlap_start = max(
                cds_start,
                exon["genome_start"],
            )

            overlap_end = min(
                cds_end,
                exon["genome_end"],
            )

            if overlap_start >= overlap_end:
                continue

            covered_cds_bases += (
                overlap_end
                - overlap_start
            )

            if tx["strand"] == "+":

                local_start = (
                    exon["tx_start"]
                    + (
                        overlap_start
                        - exon["genome_start"]
                    )
                )

                local_end = (
                    exon["tx_start"]
                    + (
                        overlap_end
                        - exon["genome_start"]
                    )
                )

            else:

                local_start = (
                    exon["tx_start"]
                    + (
                        exon["genome_end"]
                        - overlap_end
                    )
                )

                local_end = (
                    exon["tx_start"]
                    + (
                        exon["genome_end"]
                        - overlap_start
                    )
                )

            projected_intervals.append(
                (
                    local_start,
                    local_end,
                )
            )

    if (
        not projected_intervals
        or covered_cds_bases
        != expected_cds_bases
    ):
        return None

    projected_intervals.sort()

    for previous, current in zip(
        projected_intervals,
        projected_intervals[1:],
    ):

        if previous[1] != current[0]:
            return None

    cds_start = projected_intervals[0][0]
    cds_stop = projected_intervals[-1][1]

    if cds_stop <= cds_start:
        return None

    return (
        cds_start,
        cds_stop,
    )


# ============================================================
# CDS validation for footprint_generator
# ============================================================

def prepare_cds_for_simulator(
    transcript_sequence,
    cds_start,
    cds_stop,
):
    """
    Validate CDS for RiboMap footprint_generator.

    CDS must:
      - be within transcript
      - have length divisible by 3
      - start with ATG
      - end with TAA/TAG/TGA
      - contain no ambiguous bases
      - contain no internal stop codon

    cds_stop remains exclusive and INCLUDES the terminal stop codon.
    """

    if cds_start < 0:
        return None

    if cds_stop > len(transcript_sequence):
        return None

    if cds_stop <= cds_start:
        return None

    cds_sequence = transcript_sequence[
        cds_start:cds_stop
    ]

    if len(cds_sequence) < 9:
        return None

    if len(cds_sequence) % 3 != 0:
        return None

    if any(
        base not in "ACGT"
        for base in cds_sequence
    ):
        return None

    # RiboMap expects a complete CDS
    if cds_sequence[:3] != "ATG":
        return None

    if cds_sequence[-3:] not in {
        "TAA",
        "TAG",
        "TGA",
    }:
        return None

    # Reject internal stops
    for i in range(
        3,
        len(cds_sequence) - 3,
        3,
    ):
        codon = cds_sequence[i:i + 3]

        if codon in {
            "TAA",
            "TAG",
            "TGA",
        }:
            return None

    return (
        cds_start,
        cds_stop,
        cds_sequence,
    )

# ============================================================
# Build candidates
# ============================================================

def build_simulation_candidates(
    transcripts,
    gencode_info,
    protein_coding,
    gencode_cds,
    genome,
):
    """
    Build simulation candidates after the relative TPM filter.

    Biological selection:
      1. GENCODE reference transcript exists.
      2. GENCODE transcript is protein coding.
      3. GENCODE CDS annotation exists.

    Projection and CDS checks below are technical validity checks.
    """

    candidates = []

    counters = {
        "no_reference": 0,
        "reference_not_in_gencode": 0,
        "not_protein_coding": 0,
        "no_cds": 0,
        "cds_projection_failed": 0,
        "invalid_cds": 0,
        'non_primary_contig':0
    }

    for stringtie_id, tx in sorted(
        transcripts.items()
    ):

        if tx["chrom"] not in primary_chroms:
            counters["non_primary_contig"] += 1
            continue

        reference_id = tx[
            "reference_id"
        ]

        if not reference_id:

            counters[
                "no_reference"
            ] += 1

            continue

        reference_base = strip_version(
            reference_id
        )

        if reference_base not in gencode_info:

            counters[
                "reference_not_in_gencode"
            ] += 1

            continue

        if reference_base not in protein_coding:

            counters[
                "not_protein_coding"
            ] += 1

            continue

        if reference_base not in gencode_cds:

            counters[
                "no_cds"
            ] += 1

            continue

        transcript_sequence = (
            build_transcript_sequence(
                genome,
                tx,
            )
        )

        cds_result = (
            genomic_cds_to_transcript(
                tx,
                gencode_cds[
                    reference_base
                ],
            )
        )

        if cds_result is None:

            counters[
                "cds_projection_failed"
            ] += 1

            continue

        (
            cds_start,
            cds_stop,
        ) = cds_result

        prepared_cds = (
            prepare_cds_for_simulator(
                transcript_sequence,
                cds_start,
                cds_stop,
            )
        )

        if prepared_cds is None:

            counters[
                "invalid_cds"
            ] += 1

            continue

        (
            cds_start,
            cds_stop,
            cds_sequence,
        ) = prepared_cds

        gencode_tx = gencode_info[
            reference_base
        ]

        gene_id = (
            gencode_tx.get("gene_id")
            or tx.get("gene_id")
            or ""
        )

        gene_name = (
            gencode_tx.get("gene_name")
            or tx.get("gene_name")
            or ""
        )

        candidates.append(
            {
                "stringtie_transcript_id": stringtie_id,
                "gencode_transcript_id": reference_base,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "chrom": tx["chrom"],
                "strand": tx["strand"],
                "tpm": float(
                    tx["tpm"]
                ),
                "transcript_sequence": (
                    transcript_sequence
                ),
                "transcript_length": len(
                    transcript_sequence
                ),
                "cds_start": cds_start,
                "cds_stop": cds_stop,
                "cds_sequence": cds_sequence,
                "cds_length": len(
                    cds_sequence
                ),
            }
        )

    return (
        candidates,
        counters,
    )


# ============================================================
# Collapse identical CDS
# ============================================================

def collapse_identical_cds(
    candidates,
):
    """
    Collapse identical spliced CDS sequences within the same GENCODE gene.

    The representative is chosen by:
      1. highest StringTie TPM,
      2. longest transcript,
      3. GENCODE transcript ID,
      4. StringTie transcript ID.

    TPM is summed across all members.
    """

    groups = defaultdict(
        list
    )

    for candidate in candidates:

        gene_key = (
            candidate["gene_id"]
            or (
                candidate["chrom"]
                + ":"
                + candidate["strand"]
            )
        )

        key = (
            gene_key,
            candidate["cds_sequence"],
        )

        groups[
            key
        ].append(
            candidate
        )

    collapsed = []

    for members in groups.values():

        members_sorted = sorted(
            members,
            key=lambda x: (
                -x["tpm"],
                -x["transcript_length"],
                x["gencode_transcript_id"],
                x["stringtie_transcript_id"],
            ),
        )

        representative = dict(
            members_sorted[0]
        )

        representative[
            "class_tpm"
        ] = sum(
            member["tpm"]
            for member in members
        )

        representative[
            "class_size"
        ] = len(
            members
        )

        representative[
            "member_stringtie_transcript_ids"
        ] = sorted(
            {
                member[
                    "stringtie_transcript_id"
                ]
                for member in members
            }
        )

        representative[
            "member_gencode_transcript_ids"
        ] = sorted(
            {
                member[
                    "gencode_transcript_id"
                ]
                for member in members
            }
        )

        collapsed.append(
            representative
        )

    collapsed.sort(
        key=lambda x: (
            x["gene_id"],
            x["gencode_transcript_id"],
            x["stringtie_transcript_id"],
        )
    )

    return collapsed


# ============================================================
# Write filtered GTF
# ============================================================

def write_filtered_gtf(
    input_gtf,
    output_gtf,
    selected_ids,
):
    """Write representative StringTie transcripts to GTF."""

    with open(input_gtf) as src, open(
        output_gtf,
        "w",
    ) as dst:

        for line in src:

            if line.startswith("#"):

                dst.write(line)
                continue

            cols = line.rstrip().split("\t")

            if len(cols) != 9:
                continue

            attrs = parse_gtf_attributes(
                cols[8]
            )

            transcript_id = attrs.get(
                "transcript_id"
            )

            if transcript_id in selected_ids:

                dst.write(line)


# ============================================================
# Process one dataset
# ============================================================

def process_dataset(
    dataset_name,
    config,
):
    """Create simulation input files."""

    print()
    print(
        "============================================================"
    )
    print(
        f"Dataset: {dataset_name}"
    )
    print(
        "============================================================"
    )

    output_dir = (
        OUTPUT_BASE
        / dataset_name
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    transcripts = load_stringtie_gtf(
        config["stringtie_gtf"]
    )

    print(
        f"StringTie transcripts: "
        f"{len(transcripts):,}"
    )

    transcripts = (
        filter_by_relative_tpm(
            transcripts,
            MIN_RELATIVE_TPM,
        )
    )

    print(
        f"After >=3% relative TPM filtering: "
        f"{len(transcripts):,}"
    )

    (
        gencode_info,
        protein_coding,
        gencode_cds,
    ) = load_gencode_cds(
        config["gencode_gff3"]
    )

    genome = pysam.FastaFile(
        str(
            config["genome_fasta"]
        )
    )

    (
        candidates,
        counters,
    ) = build_simulation_candidates(
        transcripts=transcripts,
        gencode_info=gencode_info,
        protein_coding=protein_coding,
        gencode_cds=gencode_cds,
        genome=genome,
    )

    print(
        f"Eligible candidates before CDS deduplication: "
        f"{len(candidates):,}"
    )

    collapsed = collapse_identical_cds(
        candidates
    )

    duplicate_count = (
        len(candidates)
        - len(collapsed)
    )

    duplicate_classes = sum(
        item["class_size"] > 1
        for item in collapsed
    )

    print(
        f"After identical-CDS deduplication: "
        f"{len(collapsed):,}"
    )

    print(
        f"Removed duplicate CDS transcripts: "
        f"{duplicate_count:,}"
    )

    print(
        f"CDS classes with >1 member: "
        f"{duplicate_classes:,}"
    )

    fasta_path = (
        output_dir
        / (
            f"{dataset_name}."
            "simulation_transcripts.fa"
        )
    )

    cds_path = (
        output_dir
        / (
            f"{dataset_name}."
            "simulation_cds.txt"
        )
    )

    abundance_path = (
        output_dir
        / (
            f"{dataset_name}."
            "simulation_abundance.sf"
        )
    )

    info_path = (
        output_dir
        / (
            f"{dataset_name}."
            "simulation_transcript_info.tsv"
        )
    )

    filtered_gtf_path = (
        output_dir
        / (
            f"{dataset_name}."
            "simulation_transcripts.gtf"
        )
    )

    selected_ids = set()

    with (
        open(
            fasta_path,
            "w",
        ) as fasta_out,
        open(
            cds_path,
            "w",
        ) as cds_out,
        open(
            abundance_path,
            "w",
        ) as abundance_out,
        open(
            info_path,
            "w",
        ) as info_out,
    ):

        abundance_out.write(
            "Name\tLength\tEffectiveLength\tTPM\n"
        )

        info_out.write(
            "\t".join(
                [
                    "simulation_transcript_id",
                    "representative_stringtie_transcript_id",
                    "representative_gencode_transcript_id",
                    "member_stringtie_transcript_ids",
                    "member_gencode_transcript_ids",
                    "class_size",
                    "gene_id",
                    "gene_name",
                    "chrom",
                    "strand",
                    "tpm_sum",
                    "representative_tpm",
                    "transcript_length",
                    "cds_start",
                    "cds_stop",
                    "cds_length",
                ]
            )
            + "\n"
        )

        simulation_id_seen = set()

        for candidate in collapsed:

            simulation_id = candidate[
                "gencode_transcript_id"
            ]

            if simulation_id in simulation_id_seen:

                raise RuntimeError(
                    "Duplicate simulation transcript ID after CDS "
                    f"deduplication: {simulation_id}"
                )

            transcript_sequence = candidate[
                "transcript_sequence"
            ]

            transcript_length = candidate[
                "transcript_length"
            ]

            cds_start = candidate[
                "cds_start"
            ]

            cds_stop = candidate[
                "cds_stop"
            ]

            class_tpm = candidate[
                "class_tpm"
            ]

            fasta_out.write(
                f">{simulation_id}\n"
            )

            for i in range(
                0,
                len(transcript_sequence),
                80,
            ):

                fasta_out.write(
                    transcript_sequence[
                        i:i + 80
                    ]
                    + "\n"
                )

            cds_out.write(
                f"{simulation_id}\t"
                f"{cds_start}\t"
                f"{cds_stop}\n"
            )

            effective_length = max(
                transcript_length - 29,
                1,
            )

            abundance_out.write(
                f"{simulation_id}\t"
                f"{transcript_length}\t"
                f"{effective_length}\t"
                f"{class_tpm:.8f}\n"
            )

            info_out.write(
                "\t".join(
                    [
                        simulation_id,
                        candidate[
                            "stringtie_transcript_id"
                        ],
                        candidate[
                            "gencode_transcript_id"
                        ],
                        ";".join(
                            candidate[
                                "member_stringtie_transcript_ids"
                            ]
                        ),
                        ";".join(
                            candidate[
                                "member_gencode_transcript_ids"
                            ]
                        ),
                        str(
                            candidate[
                                "class_size"
                            ]
                        ),
                        candidate[
                            "gene_id"
                        ],
                        candidate[
                            "gene_name"
                        ],
                        candidate[
                            "chrom"
                        ],
                        candidate[
                            "strand"
                        ],
                        f"{class_tpm:.8f}",
                        f"{candidate['tpm']:.8f}",
                        str(
                            transcript_length
                        ),
                        str(
                            cds_start
                        ),
                        str(
                            cds_stop
                        ),
                        str(
                            candidate[
                                "cds_length"
                            ]
                        ),
                    ]
                )
                + "\n"
            )

            simulation_id_seen.add(
                simulation_id
            )

            selected_ids.add(
                candidate[
                    "stringtie_transcript_id"
                ]
            )

    genome.close()

    write_filtered_gtf(
        config["stringtie_gtf"],
        filtered_gtf_path,
        selected_ids,
    )

    print()
    print(
        f"Simulation transcripts: "
        f"{len(collapsed):,}"
    )

    print(
        f"No GENCODE reference ID: "
        f"{counters['no_reference']:,}"
    )

    print(
        f"Reference ID absent from GENCODE: "
        f"{counters['reference_not_in_gencode']:,}"
    )

    print(
        f"Not protein coding: "
        f"{counters['not_protein_coding']:,}"
    )

    print(
        f"No CDS annotation: "
        f"{counters['no_cds']:,}"
    )

    print(
        f"CDS projection failed: "
        f"{counters['cds_projection_failed']:,}"
    )

    print(
        f"Simulator-incompatible CDS: "
        f"{counters['invalid_cds']:,}"
    )

    print(
        f"Identical-CDS transcripts removed: "
        f"{duplicate_count:,}"
    )

    print()
    print("Output:")

    print(
        f"  FASTA     : {fasta_path}"
    )

    print(
        f"  CDS       : {cds_path}"
    )

    print(
        f"  Abundance : {abundance_path}"
    )

    print(
        f"  GTF       : {filtered_gtf_path}"
    )

    print(
        f"  Info      : {info_path}"
    )


# ============================================================
# Main
# ============================================================

def main():

    OUTPUT_BASE.mkdir(
        parents=True,
        exist_ok=True,
    )

    for dataset_name, config in DATASETS.items():

        process_dataset(
            dataset_name,
            config,
        )


if __name__ == "__main__":
    main()
