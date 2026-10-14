#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import itertools
import re
import shutil
import subprocess
import tempfile

from collections import Counter
from pathlib import Path

import pysam


# ============================================================
# Configuration
# ============================================================

BASE = Path(
    "/mnt/share/project/riboML/validation_expression/similation"
)


DATASETS = [
    {
        "name": "mouse_brain",

        "fastq": (
            BASE
            / "mouse_brain"
            / "mouse_brain.simulation.fixed.fastq"
        ),

        "transcript_info": (
            BASE
            / "mouse_brain"
            / "mouse_brain.simulation_transcript_info.tsv"
        ),

        "stringtie_gtf": (
            BASE
            / "mouse_brain"
            / "mouse_brain.simulation_transcripts.gtf"
        ),

        "mappers": [
            {
                "name": "STAR",

                "bam": (
                    BASE
                    / "mouse_brain"
                    / "STAR"
                    / "Aligned.sortedByCoord.out.bam"
                ),
            },

            {
                "name": "RiboTxAlign_LongRead",

                "bam": (
                    BASE
                    / "mouse_brain"
                    / "RiboTxAlign_LongRead"
                    / "mouse_brain_simulation_LongRead"
                    / "bam"
                    / "mouse_brain_simulation_LongRead.genome.sorted.bam"
                ),
            },
            {
                "name": "RiboTxAlign_RNASeq",
                "bam": (
                    BASE / "mouse_brain" / "RiboTxAlign_RNASeq"
                    / "mouse_brain_simulation_RNASeq" / "bam"
                    / "mouse_brain_simulation_RNASeq.genome.sorted.bam"
                ),
            },
        ],
    },

    {
        "name": "HEK293T",

        "fastq": (
            BASE
            / "HEK293T"
            / "HEK293T.simulation.fastq"
        ),

        "transcript_info": (
            BASE
            / "HEK293T"
            / "HEK293T.simulation_transcript_info.tsv"
        ),

        "stringtie_gtf": (
            BASE
            / "HEK293T"
            / "HEK293T.simulation_transcripts.gtf"
        ),

        "mappers": [
            {
                "name": "STAR",

                "bam": (
                    BASE
                    / "HEK293T"
                    / "STAR"
                    / "Aligned.sortedByCoord.out.bam"
                ),
            },

            {
                "name": "RiboTxAlign_LongRead",

                "bam": (
                    BASE
                    / "HEK293T"
                    / "RiboTxAlign_LongRead"
                    / "HEK293T_simulation_LongRead"
                    / "bam"
                    / "HEK293T_simulation_LongRead.genome.sorted.bam"
                ),
            },
            {
                "name": "RiboTxAlign_RNASeq",
                "bam": (
                    BASE / "HEK293T" / "RiboTxAlign_RNASeq"
                    / "HEK293T_simulation_RNASeq" / "bam"
                    / "HEK293T_simulation_RNASeq.genome.sorted.bam"
                ),
            },
        ],
    },
]


OUTPUT = (
    BASE
    / "simulation_mapping_accuracy_summary.tsv"
)


DEBUG_MAX = 100


def is_ribotx_mapper(name: str) -> bool:
    """Recognize legacy and separately named StringTie reference runs."""
    return name == "RiboTxAlign" or name.startswith("RiboTxAlign_")


HEADER_PATTERN = re.compile(
    r"^@(?P<read_id>\S+)\s+"
    r"(?P<simulation_transcript_id>\S+)\s+"
    r"\S+\s+"
    r"codon:\s*-?\d+\s+"
    r"base:\s*(?P<start>\d+)-(?P<end>\d+)"
)


# ============================================================
# GTF
# ============================================================

def parse_gtf_attributes(
    text: str,
) -> dict[str, str]:

    attrs = {}

    for item in (
        text
        .strip()
        .strip(";")
        .split(";")
    ):

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

        value = (
            parts[1]
            .strip()
            .strip('"')
        )

        attrs[
            key
        ] = value

    return attrs


def load_stringtie_transcripts(
    gtf_path: Path,
) -> dict[str, dict]:

    """
    Load exon structures from the simulation StringTie GTF.

    Exons are stored in transcript 5' -> 3' order.

    GTF coordinates are converted to 0-based half-open.
    """

    transcripts: dict[
        str,
        dict,
    ] = {}

    with gtf_path.open(
        "r",
        encoding="utf-8",
    ) as handle:

        for line in handle:

            if line.startswith(
                "#"
            ):
                continue

            cols = (
                line
                .rstrip()
                .split("\t")
            )

            if len(cols) != 9:
                continue

            if cols[2] != "exon":
                continue

            chrom = cols[0]

            start = (
                int(
                    cols[3]
                )
                - 1
            )

            end = int(
                cols[4]
            )

            strand = cols[6]

            attrs = (
                parse_gtf_attributes(
                    cols[8]
                )
            )

            transcript_id = (
                attrs.get(
                    "transcript_id"
                )
            )

            if not transcript_id:
                continue

            transcript = (
                transcripts.setdefault(
                    transcript_id,
                    {
                        "chrom": chrom,
                        "strand": strand,
                        "exons": [],
                    },
                )
            )

            transcript[
                "exons"
            ].append(
                (
                    start,
                    end,
                )
            )

    for transcript in (
        transcripts.values()
    ):

        transcript[
            "exons"
        ].sort(
            key=lambda exon: exon[0],

            reverse=(
                transcript[
                    "strand"
                ]
                == "-"
            ),
        )

    return transcripts


# ============================================================
# Transcript ID mapping
# ============================================================

def load_transcript_id_map(
    transcript_info_path: Path,
) -> dict[str, str]:

    """
    FASTQ header uses simulation_transcript_id.

    Simulation GTF uses StringTie IDs. CDS-collapsed simulation inputs use
    representative_stringtie_transcript_id; older inputs use stringtie_transcript_id.
    Member IDs are not interchangeable with the representative's coordinates.
    """

    result = {}

    with Path(transcript_info_path).open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        fields = reader.fieldnames or []
        target_columns = [name for name in (
            "representative_stringtie_transcript_id", "stringtie_transcript_id"
        ) if name in fields]
        if "simulation_transcript_id" not in fields or not target_columns:
            raise ValueError(
                f"Invalid simulation transcript info: {transcript_info_path}. "
                "Required: simulation_transcript_id and representative_stringtie_transcript_id "
                f"(or legacy stringtie_transcript_id). Found columns: {fields}"
            )
        for row in reader:
            simulation_id = (row.get("simulation_transcript_id") or "").strip()
            target_id = (row.get(target_columns[0]) or "").strip()
            if simulation_id in ("", ".") or target_id in ("", "."):
                raise ValueError(
                    f"Missing transcript ID in {transcript_info_path}, line {reader.line_num}"
                )
            if len(target_columns) == 2:
                legacy_id = (row.get("stringtie_transcript_id") or "").strip()
                if legacy_id not in ("", ".", target_id):
                    raise ValueError(
                        f"Conflicting representative/legacy StringTie IDs for {simulation_id} "
                        f"in {transcript_info_path}, line {reader.line_num}"
                    )
            if simulation_id in result and result[simulation_id] != target_id:
                raise ValueError(
                    f"Conflicting StringTie IDs for {simulation_id} in {transcript_info_path}, "
                    f"line {reader.line_num}"
                )
            result[simulation_id] = target_id

    return result


# ============================================================
# Transcript -> genome projection
# ============================================================

def project_transcript_interval(
    transcript: dict,
    tx_start: int,
    tx_end: int,
):

    """
    Convert transcript-coordinate interval to genomic blocks.

    Input:
        tx_start, tx_end
        0-based half-open transcript coordinates.

    Output:
        {
            chrom,
            strand,
            blocks
        }
    """

    if tx_end <= tx_start:
        return None

    chrom = transcript[
        "chrom"
    ]

    strand = transcript[
        "strand"
    ]

    exons = transcript[
        "exons"
    ]

    transcript_position = 0

    genome_blocks = []

    for (
        exon_start,
        exon_end,
    ) in exons:

        exon_length = (
            exon_end
            - exon_start
        )

        exon_tx_start = (
            transcript_position
        )

        exon_tx_end = (
            transcript_position
            + exon_length
        )

        overlap_tx_start = max(
            tx_start,
            exon_tx_start,
        )

        overlap_tx_end = min(
            tx_end,
            exon_tx_end,
        )

        if (
            overlap_tx_start
            < overlap_tx_end
        ):

            local_start = (
                overlap_tx_start
                - exon_tx_start
            )

            local_end = (
                overlap_tx_end
                - exon_tx_start
            )

            if strand == "+":

                genome_start = (
                    exon_start
                    + local_start
                )

                genome_end = (
                    exon_start
                    + local_end
                )

            else:

                genome_start = (
                    exon_end
                    - local_end
                )

                genome_end = (
                    exon_end
                    - local_start
                )

            genome_blocks.append(
                (
                    genome_start,
                    genome_end,
                )
            )

        transcript_position = (
            exon_tx_end
        )

        if (
            transcript_position
            >= tx_end
        ):
            break

    covered_length = sum(
        end - start
        for (
            start,
            end,
        ) in genome_blocks
    )

    expected_length = (
        tx_end
        - tx_start
    )

    if (
        covered_length
        != expected_length
    ):
        return None

    genome_blocks.sort(
        key=lambda block: block[0]
    )

    return {
        "chrom": chrom,
        "strand": strand,
        "blocks": tuple(
            genome_blocks
        ),
    }


# ============================================================
# FASTQ truth
# ============================================================

def open_fastq(
    path: Path,
):

    if path.suffix == ".gz":

        return gzip.open(
            path,
            "rt",
        )

    return path.open(
        "r",
        encoding="utf-8",
    )


def write_truth_tsv(
    fastq_path: Path,
    transcript_info_path: Path,
    output_path: Path,
) -> Counter:

    """
    Write compact simulation truth table.
    """

    transcript_map = (
        load_transcript_id_map(
            transcript_info_path
        )
    )

    counts = Counter()

    with (
        open_fastq(
            fastq_path
        ) as fastq,

        output_path.open(
            "w",
            encoding="utf-8",
        ) as output,
    ):

        while True:

            header = fastq.readline()

            if not header:
                break

            sequence = (
                fastq.readline()
            )

            plus = (
                fastq.readline()
            )

            quality = (
                fastq.readline()
            )

            if not quality:

                counts[
                    "truncated_fastq"
                ] += 1

                break

            counts[
                "fastq_reads"
            ] += 1

            match = (
                HEADER_PATTERN.match(
                    header.rstrip()
                )
            )

            if match is None:

                counts[
                    "header_error"
                ] += 1

                continue

            simulation_id = (
                match.group(
                    "simulation_transcript_id"
                )
            )

            stringtie_id = (
                transcript_map.get(
                    simulation_id
                )
            )

            if (
                stringtie_id is None
            ):

                counts[
                    "transcript_id_conversion_error"
                ] += 1

                continue

            read_id = (
                match.group(
                    "read_id"
                )
            )

            start = int(
                match.group(
                    "start"
                )
            )

            end = int(
                match.group(
                    "end"
                )
            )

            if end <= start:

                counts[
                    "invalid_truth_interval"
                ] += 1

                continue

            output.write(
                f"{read_id}\t"
                f"{stringtie_id}\t"
                f"{start}\t"
                f"{end}\n"
            )

            counts[
                "truth_reads"
            ] += 1

    return counts


def sort_truth_tsv(
    input_path: Path,
    output_path: Path,
) -> None:

    with output_path.open(
        "wb"
    ) as output:

        subprocess.run(
            [
                "sort",
                "-t",
                "\t",
                "-k1,1n",
                str(
                    input_path
                ),
            ],
            stdout=output,
            check=True,
        )


def iter_truth(
    path: Path,
):

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:

        for line in handle:

            (
                read_id,
                transcript_id,
                start,
                end,
            ) = (
                line
                .rstrip()
                .split("\t")
            )

            yield (
                read_id,
                transcript_id,
                int(
                    start
                ),
                int(
                    end
                ),
            )


# ============================================================
# BAM utilities
# ============================================================

def name_sort_bam(
    input_bam: Path,
    output_bam: Path,
) -> None:

    pysam.sort(
        "-n",
        "-@",
        "1",
        "-o",
        str(
            output_bam
        ),
        str(
            input_bam
        ),
    )


def iter_bam_groups(
    bam_path: Path,
):

    with pysam.AlignmentFile(
        str(
            bam_path
        ),
        "rb",
    ) as bam:

        for (
            read_id,
            records,
        ) in itertools.groupby(
            bam.fetch(
                until_eof=True
            ),

            key=lambda record: (
                record.query_name
            ),
        ):

            yield (
                read_id,
                list(
                    records
                ),
            )


def primary_alignment(
    records,
):

    for record in records:

        if (
            not record.is_unmapped
            and not record.is_secondary
            and not record.is_supplementary
        ):

            return record

    return None


def observed_blocks(
    record,
):

    return tuple(
        (
            int(
                start
            ),
            int(
                end
            ),
        )
        for (
            start,
            end,
        ) in record.get_blocks()
    )


def get_ribotx_assignment_type(
    record,
) -> str:

    """
    Return RiboTxAlign assignment type stored in XA.
    """

    if record.has_tag(
        "XA"
    ):

        return str(
            record.get_tag(
                "XA"
            )
        )

    return "unknown"


def get_original_hit_count(
    record,
) -> int:

    """
    Return original candidate count stored in NH.
    """

    if record.has_tag(
        "NH"
    ):

        return int(
            record.get_tag(
                "NH"
            )
        )

    return 1


def load_star_multi_read_ids(
    bam_path: Path,
) -> set[str]:

    """
    STAR multi-mapped reads are defined as NH > 1.
    """

    multi_ids = set()

    with pysam.AlignmentFile(
        str(
            bam_path
        ),
        "rb",
    ) as bam:

        for record in bam.fetch(
            until_eof=True
        ):

            if record.is_unmapped:
                continue

            if record.is_secondary:
                continue

            if record.is_supplementary:
                continue

            nh = (
                int(
                    record.get_tag(
                        "NH"
                    )
                )
                if record.has_tag(
                    "NH"
                )
                else 1
            )

            if nh > 1:

                multi_ids.add(
                    record.query_name
                )

    return multi_ids


# ============================================================
# Alignment classification
# ============================================================

def is_correct(
    truth,
    record,
) -> bool:

    if record.is_unmapped:
        return False

    strand = (
        "-"
        if record.is_reverse
        else "+"
    )

    return (
        record.reference_name
        == truth[
            "chrom"
        ]
        and strand
        == truth[
            "strand"
        ]
        and observed_blocks(
            record
        )
        == truth[
            "blocks"
        ]
    )


def classify(
    truth,
    record,
) -> str:

    if is_correct(
        truth,
        record,
    ):

        return "correct"

    strand = (
        "-"
        if record.is_reverse
        else "+"
    )

    blocks = (
        observed_blocks(
            record
        )
    )

    if (
        record.reference_name
        == truth[
            "chrom"
        ]
        and strand
        == truth[
            "strand"
        ]
        and blocks
    ):

        truth_start = min(
            start
            for (
                start,
                _,
            ) in truth[
                "blocks"
            ]
        )

        truth_end = max(
            end
            for (
                _,
                end,
            ) in truth[
                "blocks"
            ]
        )

        observed_start = min(
            start
            for (
                start,
                _,
            ) in blocks
        )

        observed_end = max(
            end
            for (
                _,
                end,
            ) in blocks
        )

        overlap = max(
            0,

            min(
                truth_end,
                observed_end,
            )
            - max(
                truth_start,
                observed_start,
            ),
        )

        if (
            overlap > 0
            and (
                len(
                    truth[
                        "blocks"
                    ]
                )
                > 1
                or len(
                    blocks
                )
                > 1
            )
        ):

            return "splice_error"

    return "alignment_error"


# ============================================================
# Evaluation
# ============================================================

def evaluate_bam(
    truth_path: Path,
    bam_path: Path,
    transcripts: dict,
    star_multi_ids: set[str] | None = None,
) -> Counter:

    counts = Counter()

    debug_printed = 0

    temp_bam = Path(
        tempfile.mktemp(
            suffix=(
                ".name_sorted.bam"
            )
        )
    )

    try:

        print(
            "  Name-sort BAM..."
        )

        name_sort_bam(
            bam_path,
            temp_bam,
        )

        bam_iterator = (
            iter_bam_groups(
                temp_bam
            )
        )

        try:

            (
                bam_read_id,
                records,
            ) = next(
                bam_iterator
            )

        except StopIteration:

            bam_read_id = None
            records = []

        print(
            "  Compare..."
        )

        for (
            read_id,
            transcript_id,
            tx_start,
            tx_end,
        ) in iter_truth(
            truth_path
        ):

            counts[
                "evaluated_reads"
            ] += 1

            # ------------------------------------------------
            # Advance BAM until BAM ID >= truth ID
            # ------------------------------------------------

            while (
                bam_read_id is not None
                and int(
                    bam_read_id
                )
                < int(
                    read_id
                )
            ):

                try:

                    (
                        bam_read_id,
                        records,
                    ) = next(
                        bam_iterator
                    )

                except StopIteration:

                    bam_read_id = None
                    records = []

                    break

            transcript = (
                transcripts.get(
                    transcript_id
                )
            )

            if transcript is None:

                counts[
                    "truth_transcript_missing"
                ] += 1

                continue

            truth = (
                project_transcript_interval(
                    transcript,
                    tx_start,
                    tx_end,
                )
            )

            if truth is None:

                counts[
                    "truth_projection_error"
                ] += 1

                if (
                    debug_printed
                    < DEBUG_MAX
                ):

                    print()
                    print(
                        "DEBUG TRUTH PROJECTION ERROR"
                    )
                    print(
                        "read_id      :",
                        read_id,
                    )
                    print(
                        "transcript_id:",
                        transcript_id,
                    )
                    print(
                        "tx_start-end :",
                        tx_start,
                        tx_end,
                    )
                    print(
                        "strand       :",
                        transcript[
                            "strand"
                        ],
                    )
                    print(
                        "exons        :",
                        transcript[
                            "exons"
                        ],
                    )

                    debug_printed += 1

                continue

            if (
                len(
                    truth[
                        "blocks"
                    ]
                )
                > 1
            ):

                counts[
                    "truth_spliced"
                ] += 1

            # ------------------------------------------------
            # No BAM record
            # ------------------------------------------------

            if (
                bam_read_id is None
                or bam_read_id
                != read_id
            ):

                counts[
                    "unmapped"
                ] += 1

                if (
                    star_multi_ids
                    is not None
                    and read_id
                    in star_multi_ids
                ):

                    counts[
                        "star_multi_total"
                    ] += 1

                    counts[
                        "star_multi_unmapped"
                    ] += 1

                continue

            primary = (
                primary_alignment(
                    records
                )
            )

            # ------------------------------------------------
            # No primary mapped alignment
            # ------------------------------------------------

            if primary is None:

                counts[
                    "unmapped"
                ] += 1

                if (
                    star_multi_ids
                    is not None
                    and read_id
                    in star_multi_ids
                ):

                    counts[
                        "star_multi_total"
                    ] += 1

                    counts[
                        "star_multi_unmapped"
                    ] += 1

            else:

                status = (
                    classify(
                        truth,
                        primary,
                    )
                )

                counts[
                    status
                ] += 1

                assignment_type = (
                    get_ribotx_assignment_type(
                        primary
                    )
                )

                original_nh = (
                    get_original_hit_count(
                        primary
                    )
                )

                # ============================================
                # NEW:
                # Correct mapping split by unique / multi
                # ============================================

                if status == "correct":

                    if original_nh == 1:

                        counts[
                            "correct_unique"
                        ] += 1

                    else:

                        counts[
                            "correct_multi"
                        ] += 1

                # ============================================
                # Optional diagnostic:
                # classify all mapped reads by NH
                # ============================================

                if original_nh == 1:

                    counts[
                        "mapped_unique"
                    ] += 1

                else:

                    counts[
                        "mapped_multi"
                    ] += 1

                # ============================================
                # RiboTxAlign multimap assignment
                # ============================================

                if (
                    assignment_type
                    == "multimap_assigned"
                ):

                    counts[
                        "multi_assigned_total"
                    ] += 1

                    counts[
                        f"multi_assigned_{status}"
                    ] += 1

                    counts[
                        f"multi_assigned_nh_{original_nh}"
                    ] += 1

                    counts[
                        f"multi_assigned_nh_{original_nh}_total"
                    ] += 1

                    if status == "correct":

                        counts[
                            f"multi_assigned_nh_{original_nh}_correct"
                        ] += 1

                # ============================================
                # STAR multi -> RiboTxAlign outcome
                # ============================================

                if (
                    star_multi_ids
                    is not None
                    and read_id
                    in star_multi_ids
                ):

                    counts[
                        "star_multi_total"
                    ] += 1

                    counts[
                        f"star_multi_{status}"
                    ] += 1

                    if status == "correct":

                        if (
                            assignment_type
                            == "multimap_assigned"
                        ):

                            counts[
                                "star_multi_correct_via_multimap_assigned"
                            ] += 1

                        else:

                            counts[
                                "star_multi_correct_as_non_multimap_assigned"
                            ] += 1

                    if (
                        assignment_type
                        == "multimap_assigned"
                    ):

                        counts[
                            f"star_multi_ribotx_nh_{original_nh}_total"
                        ] += 1

                        if status == "correct":

                            counts[
                                f"star_multi_ribotx_nh_{original_nh}_correct"
                            ] += 1

                # ============================================
                # Debug incorrect alignments
                # ============================================

                if (
                    status != "correct"
                    and debug_printed
                    < DEBUG_MAX
                ):

                    print()
                    print(
                        "========================================"
                    )
                    print(
                        "DEBUG",
                        status.upper(),
                    )
                    print(
                        "========================================"
                    )

                    print(
                        "read_id      :",
                        read_id,
                    )

                    print(
                        "transcript_id:",
                        transcript_id,
                    )

                    print(
                        "tx_start-end :",
                        tx_start,
                        tx_end,
                    )

                    print(
                        "truth chrom  :",
                        truth[
                            "chrom"
                        ],
                    )

                    print(
                        "truth strand :",
                        truth[
                            "strand"
                        ],
                    )

                    print(
                        "truth blocks :",
                        truth[
                            "blocks"
                        ],
                    )

                    print(
                        "BAM chrom    :",
                        primary.reference_name,
                    )

                    print(
                        "BAM strand   :",
                        (
                            "-"
                            if primary.is_reverse
                            else "+"
                        ),
                    )

                    print(
                        "BAM cigar    :",
                        primary.cigarstring,
                    )

                    print(
                        "BAM blocks   :",
                        observed_blocks(
                            primary
                        ),
                    )

                    print(
                        "BAM MAPQ     :",
                        primary.mapping_quality,
                    )

                    print(
                        "BAM XA       :",
                        assignment_type,
                    )

                    print(
                        "BAM NH       :",
                        original_nh,
                    )

                    debug_printed += 1

            # ------------------------------------------------
            # Check whether any candidate alignment is correct
            # ------------------------------------------------

            if any(
                is_correct(
                    truth,
                    record,
                )
                for record
                in records
                if not record.is_unmapped
            ):

                counts[
                    "any_candidate_correct"
                ] += 1

            # ------------------------------------------------
            # Advance BAM group
            # ------------------------------------------------

            try:

                (
                    bam_read_id,
                    records,
                ) = next(
                    bam_iterator
                )

            except StopIteration:

                bam_read_id = None
                records = []

    finally:

        if temp_bam.exists():

            temp_bam.unlink()

    counts[
        "incorrect"
    ] = (
        counts[
            "splice_error"
        ]
        + counts[
            "alignment_error"
        ]
    )

    # --------------------------------------------------------
    # Consistency check
    # --------------------------------------------------------

    correct_split_total = (
        counts[
            "correct_unique"
        ]
        + counts[
            "correct_multi"
        ]
    )

    if (
        correct_split_total
        != counts[
            "correct"
        ]
    ):

        print()
        print(
            "WARNING: correct count mismatch"
        )

        print(
            "correct:",
            counts[
                "correct"
            ],
        )

        print(
            "correct_unique + correct_multi:",
            correct_split_total,
        )

    return counts


# ============================================================
# Summary output
# ============================================================

def write_summary(
    results,
):

    metrics = [
        "fastq_reads",
        "truth_reads",
        "evaluated_reads",

        # Overall correct mapping
        "correct",

        # Correct mapping split
        "correct_unique",
        "correct_multi",

        # All mapped reads split
        "mapped_unique",
        "mapped_multi",

        # Errors
        "incorrect",
        "splice_error",
        "alignment_error",
        "unmapped",

        # Candidate diagnostic
        "any_candidate_correct",

        # RiboTxAlign multi-assignment
        "multi_assigned_total",
        "multi_assigned_correct",
        "multi_assigned_splice_error",
        "multi_assigned_alignment_error",

        # STAR multi subset
        "star_multi_total",
        "star_multi_correct",
        "star_multi_splice_error",
        "star_multi_alignment_error",
        "star_multi_unmapped",

        "star_multi_correct_via_multimap_assigned",
        "star_multi_correct_as_non_multimap_assigned",

        # Truth diagnostics
        "truth_spliced",
        "truth_transcript_missing",
        "truth_projection_error",
        "header_error",
        "transcript_id_conversion_error",
        "invalid_truth_interval",
        "truncated_fastq",
    ]

    with OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "dataset",
                "mapper",
                "metric",
                "count",
                "percent_of_evaluated_reads",
                "percent_of_multi_assigned",
                "percent_of_star_multi",
            ]
        )

        for (
            dataset_name,
            mapper_name,
            counts,
        ) in results:

            evaluated_reads = (
                counts[
                    "evaluated_reads"
                ]
            )

            multi_assigned_total = (
                counts[
                    "multi_assigned_total"
                ]
            )

            star_multi_total = (
                counts[
                    "star_multi_total"
                ]
            )

            for metric in metrics:

                count = counts[
                    metric
                ]

                percent_evaluated = (
                    100.0
                    * count
                    / evaluated_reads
                    if evaluated_reads
                    else 0.0
                )

                if metric.startswith(
                    "multi_assigned_"
                ):

                    percent_multi = (
                        100.0
                        * count
                        / multi_assigned_total
                        if multi_assigned_total
                        else 0.0
                    )

                else:

                    percent_multi = ""

                if metric.startswith(
                    "star_multi_"
                ):

                    percent_star_multi = (
                        100.0
                        * count
                        / star_multi_total
                        if star_multi_total
                        else 0.0
                    )

                else:

                    percent_star_multi = ""

                writer.writerow(
                    [
                        dataset_name,
                        mapper_name,
                        metric,
                        count,

                        f"{percent_evaluated:.4f}",

                        (
                            f"{percent_multi:.4f}"
                            if (
                                percent_multi
                                != ""
                            )
                            else ""
                        ),

                        (
                            f"{percent_star_multi:.4f}"
                            if (
                                percent_star_multi
                                != ""
                            )
                            else ""
                        ),
                    ]
                )


# ============================================================
# NH accuracy output
# ============================================================

def write_nh_accuracy_tables(
    results,
):

    """
    Write:
        observed accuracy
        random expectation = 1/N
        read count

    by NH.
    """

    overall_output = (
        BASE
        / "ribotx_multihit_accuracy_by_nh.tsv"
    )

    star_multi_output = (
        BASE
        / "star_multi_ribotx_accuracy_by_nh.tsv"
    )

    def write_table(
        output_path,
        prefix,
        nh_index,
        nh_column_name,
    ):

        with output_path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:

            writer = csv.writer(
                handle,
                delimiter="\t",
            )

            writer.writerow(
                [
                    "dataset",
                    nh_column_name,
                    "mapper",
                    "total",
                    "correct",
                    "observed_accuracy_percent",
                    "random_expectation_percent",
                    "fold_over_random",
                    "absolute_gain_percent_point",
                ]
            )

            for (
                dataset_name,
                mapper_name,
                counts,
            ) in results:

                if not is_ribotx_mapper(mapper_name):
                    continue

                nh_values = sorted(
                    {
                        int(
                            key.split(
                                "_"
                            )[
                                nh_index
                            ]
                        )
                        for key in counts
                        if (
                            key.startswith(
                                prefix
                            )
                            and key.endswith(
                                "_total"
                            )
                        )
                    }
                )

                for nh in nh_values:

                    total_key = (
                        f"{prefix}"
                        f"{nh}_total"
                    )

                    correct_key = (
                        f"{prefix}"
                        f"{nh}_correct"
                    )

                    total = (
                        counts[
                            total_key
                        ]
                    )

                    correct = (
                        counts[
                            correct_key
                        ]
                    )

                    observed = (
                        100.0
                        * correct
                        / total
                        if total
                        else 0.0
                    )

                    random_expectation = (
                        100.0
                        / nh
                        if nh > 0
                        else 0.0
                    )

                    fold_over_random = (
                        observed
                        / random_expectation
                        if random_expectation
                        else 0.0
                    )

                    gain = (
                        observed
                        - random_expectation
                    )

                    writer.writerow(
                        [
                            dataset_name,
                            nh,
                            mapper_name,
                            total,
                            correct,
                            f"{observed:.4f}",
                            f"{random_expectation:.4f}",
                            f"{fold_over_random:.4f}",
                            f"{gain:.4f}",
                        ]
                    )

    # --------------------------------------------------------
    # All RiboTxAlign multimap-assigned reads
    # --------------------------------------------------------

    write_table(
        overall_output,
        "multi_assigned_nh_",
        3,
        "nh",
    )

    # --------------------------------------------------------
    # STAR multi-mapped reads only
    # --------------------------------------------------------

    write_table(
        star_multi_output,
        "star_multi_ribotx_nh_",
        4,
        "ribotx_nh",
    )

    print()
    print(
        f"NH accuracy table: "
        f"{overall_output}"
    )

    print(
        f"STAR-multi NH table: "
        f"{star_multi_output}"
    )


# ============================================================
# Console summary
# ============================================================

def print_mapper_summary(
    mapper_name,
    counts,
):

    evaluated = (
        counts[
            "evaluated_reads"
        ]
    )

    correct = (
        counts[
            "correct"
        ]
    )

    correct_unique = (
        counts[
            "correct_unique"
        ]
    )

    correct_multi = (
        counts[
            "correct_multi"
        ]
    )

    correct_rate = (
        100.0
        * correct
        / evaluated
        if evaluated
        else 0.0
    )

    unique_rate = (
        100.0
        * correct_unique
        / evaluated
        if evaluated
        else 0.0
    )

    multi_rate = (
        100.0
        * correct_multi
        / evaluated
        if evaluated
        else 0.0
    )

    print()
    print(
        f"  Mapper          : "
        f"{mapper_name}"
    )

    print(
        f"  Evaluated       : "
        f"{evaluated:,}"
    )

    print(
        f"  Correct         : "
        f"{correct:,} "
        f"({correct_rate:.2f}%)"
    )

    print(
        f"    Unique hit    : "
        f"{correct_unique:,} "
        f"({unique_rate:.2f}%)"
    )

    print(
        f"    Multi hit     : "
        f"{correct_multi:,} "
        f"({multi_rate:.2f}%)"
    )

    print(
        f"  Splice error    : "
        f"{counts['splice_error']:,}"
    )

    print(
        f"  Alignment error : "
        f"{counts['alignment_error']:,}"
    )

    print(
        f"  Unmapped        : "
        f"{counts['unmapped']:,}"
    )


# ============================================================
# Main
# ============================================================

def main():

    results = []

    for dataset in DATASETS:

        print()
        print(
            "============================================================"
        )

        print(
            dataset[
                "name"
            ]
        )

        print(
            "============================================================"
        )

        # ----------------------------------------------------
        # Check files
        # ----------------------------------------------------

        for key in (
            "fastq",
            "transcript_info",
            "stringtie_gtf",
        ):

            path = Path(
                dataset[
                    key
                ]
            )

            if not path.exists():

                raise FileNotFoundError(
                    f"{key}: "
                    f"{path}"
                )

        # ----------------------------------------------------
        # Load transcript structures
        # ----------------------------------------------------

        transcripts = (
            load_stringtie_transcripts(
                dataset[
                    "stringtie_gtf"
                ]
            )
        )

        print(
            "StringTie transcripts:",
            len(
                transcripts
            ),
        )

        temp_dir = Path(
            tempfile.mkdtemp(
                prefix=(
                    "simulation_truth_"
                    f"{dataset['name']}_"
                )
            )
        )

        truth_unsorted = (
            temp_dir
            / "truth.unsorted.tsv"
        )

        truth_sorted = (
            temp_dir
            / "truth.sorted.tsv"
        )

        try:

            # ------------------------------------------------
            # Prepare truth
            # ------------------------------------------------

            prep_counts = (
                write_truth_tsv(
                    fastq_path=(
                        dataset[
                            "fastq"
                        ]
                    ),
                    transcript_info_path=(
                        dataset[
                            "transcript_info"
                        ]
                    ),
                    output_path=(
                        truth_unsorted
                    ),
                )
            )

            sort_truth_tsv(
                truth_unsorted,
                truth_sorted,
            )

            # ------------------------------------------------
            # STAR multi read IDs
            # ------------------------------------------------

            star_bam_path = Path(
                next(
                    mapper[
                        "bam"
                    ]
                    for mapper
                    in dataset[
                        "mappers"
                    ]
                    if (
                        mapper[
                            "name"
                        ]
                        == "STAR"
                    )
                )
            )

            print()
            print(
                "Load STAR multi read IDs..."
            )

            star_multi_ids = (
                load_star_multi_read_ids(
                    star_bam_path
                )
            )

            print(
                "STAR multi reads:",
                f"{len(star_multi_ids):,}",
            )

            # ------------------------------------------------
            # Evaluate each mapper
            # ------------------------------------------------

            for mapper in dataset[
                "mappers"
            ]:

                bam_path = Path(
                    mapper[
                        "bam"
                    ]
                )

                if not bam_path.exists():

                    raise FileNotFoundError(
                        f"bam: "
                        f"{bam_path}"
                    )

                print()
                print(
                    "------------------------------------------------------------"
                )

                print(
                    mapper[
                        "name"
                    ]
                )

                print(
                    "------------------------------------------------------------"
                )

                counts = (
                    evaluate_bam(
                        truth_path=(
                            truth_sorted
                        ),
                        bam_path=(
                            bam_path
                        ),
                        transcripts=(
                            transcripts
                        ),
                        star_multi_ids=(
                            star_multi_ids
                            if is_ribotx_mapper(mapper["name"])
                            else None
                        ),
                    )
                )

                counts.update(
                    prep_counts
                )

                results.append(
                    (
                        dataset[
                            "name"
                        ],
                        mapper[
                            "name"
                        ],
                        counts,
                    )
                )

                print_mapper_summary(
                    mapper_name=(
                        mapper[
                            "name"
                        ]
                    ),
                    counts=(
                        counts
                    ),
                )

                # --------------------------------------------
                # RiboTxAlign-specific summary
                # --------------------------------------------

                if is_ribotx_mapper(mapper["name"]):

                    multi_total = (
                        counts[
                            "multi_assigned_total"
                        ]
                    )

                    multi_correct = (
                        counts[
                            "multi_assigned_correct"
                        ]
                    )

                    multi_assignment_rate = (
                        100.0
                        * multi_correct
                        / multi_total
                        if multi_total
                        else 0.0
                    )

                    print()
                    print(
                        "  Multi assigned  :",
                        f"{multi_total:,}",
                    )

                    print(
                        "  Multi -> correct:",
                        f"{multi_correct:,}",
                    )

                    print(
                        "  Multi -> splice error:",
                        f"{counts['multi_assigned_splice_error']:,}",
                    )

                    print(
                        "  Multi -> alignment error:",
                        f"{counts['multi_assigned_alignment_error']:,}",
                    )

                    print(
                        "  Multi assignment accuracy:",
                        f"{multi_assignment_rate:.2f}%",
                    )

                    star_multi_total = (
                        counts[
                            "star_multi_total"
                        ]
                    )

                    star_multi_correct = (
                        counts[
                            "star_multi_correct"
                        ]
                    )

                    star_multi_rate = (
                        100.0
                        * star_multi_correct
                        / star_multi_total
                        if star_multi_total
                        else 0.0
                    )

                    print()
                    print(
                        "  STAR multi total:",
                        f"{star_multi_total:,}",
                    )

                    print(
                        "  STAR multi -> RiboTx correct:",
                        f"{star_multi_correct:,}",
                    )

                    print(
                        "    via multimap_assigned:",
                        f"{counts['star_multi_correct_via_multimap_assigned']:,}",
                    )

                    print(
                        "    as non-multimap_assigned:",
                        f"{counts['star_multi_correct_as_non_multimap_assigned']:,}",
                    )

                    print(
                        "  STAR multi -> splice error:",
                        f"{counts['star_multi_splice_error']:,}",
                    )

                    print(
                        "  STAR multi -> alignment error:",
                        f"{counts['star_multi_alignment_error']:,}",
                    )

                    print(
                        "  STAR multi -> unmapped:",
                        f"{counts['star_multi_unmapped']:,}",
                    )

                    print(
                        "  STAR multi rescue accuracy:",
                        f"{star_multi_rate:.2f}%",
                    )

        finally:

            shutil.rmtree(
                temp_dir,
                ignore_errors=True,
            )

    # ========================================================
    # Write outputs
    # ========================================================

    write_summary(
        results
    )

    write_nh_accuracy_tables(
        results
    )

    print()
    print(
        "============================================================"
    )

    print(
        "Finished"
    )

    print(
        "============================================================"
    )

    print(
        f"Summary: "
        f"{OUTPUT}"
    )


if __name__ == "__main__":
    main()
