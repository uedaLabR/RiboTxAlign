from collections import defaultdict
from itertools import combinations
from pathlib import Path
import re

import pandas as pd
import pysam


# ============================================================
# Configuration
# ============================================================

BASE = Path(
    "/mnt/share/project/riboML/validation"
)

OUTPUT_DIR = (
    BASE
    / "bam_comparison_four_mappers"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Dataset configuration
#
# Four mapping methods:
#
# 1. STAR_genome
# 2. STAR_transcript
# 3. RiboTxAlign_longread_stringtie
# 4. RiboTxAlign_shortread_stringtie
# ============================================================

DATASETS = {
    # "HEK293T": {
    #     "samples": [
    #         "SRR26942051",
    #         "SRR26942052",
    #         "SRR26942053",
    #         "SRR26942054",
    #     ],
    #
    #     "star_genome_dir": (
    #         BASE
    #         / "star_validation"
    #     ),
    #
    #     "star_transcript_dir": (
    #         BASE
    #         / "star_transcript_validation_HEK293T"
    #     ),
    #
    #     # ----------------------------------------------------
    #     # CHANGE THESE TWO PATHS IF NECESSARY
    #     # ----------------------------------------------------
    #
    #     # RiboTxAlign built from long-read StringTie transcripts
    #     "ribotx_longread_stringtie_dir": (
    #         BASE
    #         / "ribotx_validation"
    #     ),
    #
    #     # RiboTxAlign built from short-read RNA-seq StringTie
    #     # transcripts
    #     "ribotx_shortread_stringtie_dir": (
    #         BASE
    #         / "ribotx_validation_rnaseq"
    #     ),
    # },

    "mouse_brain": {
        "samples": [
            "SRR2163083",
            "SRR2163093",
            "SRR2163103",
        ],

        "star_genome_dir": (
            BASE
            / "star_validation"
        ),

        "star_transcript_dir": (
            BASE
            / "star_transcript_validation_mouse_brain"
        ),

        # ----------------------------------------------------
        # CHANGE THESE TWO PATHS IF NECESSARY
        # ----------------------------------------------------

        "ribotx_longread_stringtie_dir": (
            BASE
            / "ribotx_validation"
        ),

        "ribotx_shortread_stringtie_dir": (
            BASE
            / "ribotx_validation_rnaseq"
        ),
    },
}


# ============================================================
# General utilities
# ============================================================

def safe_percent(
    numerator: int,
    denominator: int,
) -> float:

    if denominator == 0:
        return 0.0

    return (
        numerator
        / denominator
        * 100.0
    )


def has_splice_cigar(
    cigar,
) -> bool:

    if cigar is None:
        return False

    return "N" in cigar


# ============================================================
# Read ID normalization
# ============================================================

def normalize_read_id(
    read_id: str,
) -> str:
    """
    Normalize read IDs across mapping methods.

    Examples
    --------

    SRR26942051.34957075_CGACTTAGC

    SRR26942051.34957075|umi=CCGTCGCCC|family_size=1

    SRR26942051.34957075

    become

    SRR26942051.34957075
    """

    read_id = (
        read_id
        .split()[0]
    )

    match = re.search(
        r"(SRR\d+\.\d+)",
        read_id,
    )

    if match is not None:
        return match.group(1)

    return read_id


# ============================================================
# Mapping categories
# ============================================================

def make_nh_category(
    status: str,
    nh: int,
) -> str:

    if (
        status == "unmapped"
        or nh == 0
    ):
        return "unmapped"

    if status == "unique":
        return f"unique_NH{nh}"

    return f"multi_NH{nh}"


def mapping_category_sort_key(
    category: str,
):

    if category == "unmapped":
        return 0, 0

    if "_NH" not in category:
        return 99, 999999

    status, nh_text = (
        category.rsplit(
            "_NH",
            1,
        )
    )

    status_order = {
        "unique": 1,
        "multi": 2,
    }

    try:
        nh = int(
            nh_text
        )

    except ValueError:
        nh = 999999

    return (
        status_order.get(
            status,
            99,
        ),
        nh,
    )


def sort_mapping_transition(
    table: pd.DataFrame,
) -> pd.DataFrame:

    row_order = sorted(
        table.index,
        key=mapping_category_sort_key,
    )

    column_order = sorted(
        table.columns,
        key=mapping_category_sort_key,
    )

    return table.reindex(
        index=row_order,
        columns=column_order,
    )


# ============================================================
# BAM discovery
# ============================================================

def remove_bam_suffix(
    bam_path: Path,
) -> str:
    """
    Obtain sample name from a BAM filename.

    The longest suffixes must be tested first.
    """

    name = bam_path.name

    suffixes = [
        ".Aligned.sortedByCoord.out.dedup.bam",
        ".Aligned.sortedByCoord.out.deduplicated.bam",
        ".deduplicated.bam",
        ".dedup.bam",
        ".Aligned.sortedByCoord.out.bam",
        ".genome.sorted.bam",
        ".sorted.bam",
        ".bam",
    ]

    for suffix in suffixes:

        if name.endswith(
            suffix
        ):
            return name[
                :-len(suffix)
            ]

    return bam_path.stem


def choose_single_bam(
    candidates,
    description: str,
) -> Path:

    candidates = sorted(
        set(candidates)
    )

    if not candidates:

        raise FileNotFoundError(
            f"{description}: BAM not found"
        )

    if len(candidates) == 1:
        return candidates[0]

    # --------------------------------------------------------
    # Prefer UMI-deduplicated BAM
    # --------------------------------------------------------

    dedup_candidates = [
        path
        for path in candidates
        if (
            path.name.endswith(
                ".dedup.bam"
            )
            or path.name.endswith(
                ".deduplicated.bam"
            )
        )
    ]

    if len(dedup_candidates) == 1:
        return dedup_candidates[0]

    if len(dedup_candidates) > 1:

        text = "\n".join(
            f"    {path}"
            for path in dedup_candidates
        )

        raise RuntimeError(
            f"{description}: "
            f"multiple deduplicated BAM files found:\n"
            f"{text}"
        )

    # --------------------------------------------------------
    # Prefer STAR coordinate-sorted BAM
    # --------------------------------------------------------

    coordinate_sorted = [
        path
        for path in candidates
        if (
            "Aligned.sortedByCoord.out.bam"
            in path.name
        )
    ]

    if len(coordinate_sorted) == 1:
        return coordinate_sorted[0]

    # --------------------------------------------------------
    # Prefer RiboTxAlign genome.sorted.bam
    # --------------------------------------------------------

    genome_sorted = [
        path
        for path in candidates
        if path.name.endswith(
            ".genome.sorted.bam"
        )
    ]

    if len(genome_sorted) == 1:
        return genome_sorted[0]

    text = "\n".join(
        f"    {path}"
        for path in candidates
    )

    raise RuntimeError(
        f"{description}: "
        f"multiple BAM files found:\n"
        f"{text}"
    )


def find_star_genome_bam(
    directory: Path,
    srr_id: str,
) -> Path:
    """
    STAR genome BAMs are stored directly under star_validation.
    """

    candidates = [
        path
        for path in directory.glob(
            "*.bam"
        )
        if (
            path.is_file()
            and srr_id in path.name
        )
    ]

    return choose_single_bam(
        candidates,
        f"STAR genome {srr_id}",
    )


def find_star_transcript_bam(
    directory: Path,
    srr_id: str,
) -> Path:
    """
    STAR transcript BAMs are inside sample directories.

    The UMI-deduplicated BAM is explicitly preferred:

        SRRxxxx.Aligned.sortedByCoord.out.dedup.bam

    Example:

        star_transcript_validation_HEK293T/
            SRR26942051/
                SRR26942051.Aligned.sortedByCoord.out.bam
                SRR26942051.Aligned.sortedByCoord.out.dedup.bam
    """

    sample_dir = (
        directory
        / srr_id
    )

    # --------------------------------------------------------
    # First choice: exact deduplicated filename
    # --------------------------------------------------------

    dedup_bam = (
        sample_dir
        / (
            f"{srr_id}."
            "Aligned.sortedByCoord.out.dedup.bam"
        )
    )

    if dedup_bam.is_file():
        return dedup_bam

    # --------------------------------------------------------
    # Second choice: any dedup BAM in sample directory
    # --------------------------------------------------------

    dedup_candidates = [
        path
        for path in sample_dir.glob(
            "*.dedup.bam"
        )
        if path.is_file()
    ]

    if len(dedup_candidates) == 1:
        return dedup_candidates[0]

    if len(dedup_candidates) > 1:
        return choose_single_bam(
            dedup_candidates,
            f"STAR transcript dedup {srr_id}",
        )

    # --------------------------------------------------------
    # Fallback: original STAR transcript BAM
    # --------------------------------------------------------

    original_bam = (
        sample_dir
        / (
            f"{srr_id}."
            "Aligned.sortedByCoord.out.bam"
        )
    )

    if original_bam.is_file():

        print(
            "WARNING: deduplicated STAR transcript BAM "
            f"not found for {srr_id}."
        )

        print(
            f"Using original BAM instead: {original_bam}"
        )

        return original_bam

    # --------------------------------------------------------
    # Last fallback: search recursively
    # --------------------------------------------------------

    candidates = []

    if sample_dir.exists():

        candidates.extend(
            sample_dir.rglob(
                "*.bam"
            )
        )

    if not candidates:

        for path in directory.iterdir():

            if not path.is_dir():
                continue

            if srr_id not in path.name:
                continue

            candidates.extend(
                path.rglob(
                    "*.bam"
                )
            )

    return choose_single_bam(
        candidates,
        f"STAR transcript {srr_id}",
    )


def find_ribotx_bam(
    directory: Path,
    sample_name: str,
    srr_id: str,
    description: str,
) -> Path:
    """
    Locate RiboTxAlign genome BAM.

    Expected structure:

        directory/
            sample_name/
                bam/
                    sample_name.genome.sorted.bam
    """

    expected = (
        directory
        / sample_name
        / "bam"
        / (
            f"{sample_name}."
            "genome.sorted.bam"
        )
    )

    if expected.is_file():
        return expected

    # Some directories may use SRR ID directly.
    expected_srr = (
        directory
        / srr_id
        / "bam"
        / (
            f"{srr_id}."
            "genome.sorted.bam"
        )
    )

    if expected_srr.is_file():
        return expected_srr

    candidates = [
        path
        for path in directory.rglob(
            "*.genome.sorted.bam"
        )
        if (
            srr_id
            in str(path)
        )
    ]

    return choose_single_bam(
        candidates,
        f"{description} {srr_id}",
    )


def find_sample_bams(
    dataset_name: str,
    config: dict,
):

    results = []

    for srr_id in config[
        "samples"
    ]:

        print()
        print(
            "----------------------------------------"
        )
        print(
            f"Searching BAMs: "
            f"{dataset_name} / {srr_id}"
        )
        print(
            "----------------------------------------"
        )

        # ----------------------------------------------------
        # STAR genome
        # ----------------------------------------------------

        star_genome_bam = (
            find_star_genome_bam(
                config[
                    "star_genome_dir"
                ],
                srr_id,
            )
        )

        sample_name = (
            remove_bam_suffix(
                star_genome_bam
            )
        )

        # ----------------------------------------------------
        # STAR transcript
        # ----------------------------------------------------

        star_transcript_bam = (
            find_star_transcript_bam(
                config[
                    "star_transcript_dir"
                ],
                srr_id,
            )
        )

        # ----------------------------------------------------
        # RiboTxAlign long-read StringTie
        # ----------------------------------------------------

        ribotx_longread_bam = (
            find_ribotx_bam(
                directory=(
                    config[
                        "ribotx_longread_stringtie_dir"
                    ]
                ),
                sample_name=(
                    sample_name
                ),
                srr_id=(
                    srr_id
                ),
                description=(
                    "RiboTxAlign long-read StringTie"
                ),
            )
        )

        # ----------------------------------------------------
        # RiboTxAlign short-read StringTie
        # ----------------------------------------------------

        ribotx_shortread_bam = (
            find_ribotx_bam(
                directory=(
                    config[
                        "ribotx_shortread_stringtie_dir"
                    ]
                ),
                sample_name=(
                    sample_name
                ),
                srr_id=(
                    srr_id
                ),
                description=(
                    "RiboTxAlign short-read StringTie"
                ),
            )
        )

        print(
            f"Sample name                  : "
            f"{sample_name}"
        )

        print(
            f"STAR genome                  : "
            f"{star_genome_bam}"
        )

        print(
            f"STAR transcript              : "
            f"{star_transcript_bam}"
        )

        print(
            f"RiboTx long-read StringTie   : "
            f"{ribotx_longread_bam}"
        )

        print(
            f"RiboTx short-read StringTie  : "
            f"{ribotx_shortread_bam}"
        )

        results.append(
            {
                "dataset_name": (
                    dataset_name
                ),
                "sample_name": (
                    sample_name
                ),
                "srr_id": (
                    srr_id
                ),
                "star_genome_bam": (
                    star_genome_bam
                ),
                "star_transcript_bam": (
                    star_transcript_bam
                ),
                "ribotx_longread_bam": (
                    ribotx_longread_bam
                ),
                "ribotx_shortread_bam": (
                    ribotx_shortread_bam
                ),
            }
        )

    return results


# ============================================================
# STAR BAM loader
# ============================================================

def load_star_reads(
    bam_path: Path,
) -> dict:

    records_by_read = defaultdict(
        list
    )

    with pysam.AlignmentFile(
        bam_path,
        "rb",
    ) as bam:

        for record in bam.fetch(
            until_eof=True,
        ):

            if record.is_supplementary:
                continue

            read_id = (
                normalize_read_id(
                    record.query_name
                )
            )

            records_by_read[
                read_id
            ].append(
                record
            )

    reads = {}

    for (
        read_id,
        records,
    ) in records_by_read.items():

        mapped_records = [
            record
            for record in records
            if not record.is_unmapped
        ]

        # ----------------------------------------------------
        # Unmapped
        # ----------------------------------------------------

        if not mapped_records:

            sequence = None

            if records:
                sequence = (
                    records[
                        0
                    ].query_sequence
                )

            reads[
                read_id
            ] = {
                "status": (
                    "unmapped"
                ),
                "nh": 0,
                "sequence": (
                    sequence
                ),
                "alignments": [],
                "spliced": False,
            }

            continue

        # ----------------------------------------------------
        # Primary alignment
        # ----------------------------------------------------

        primary_records = [
            record
            for record in mapped_records
            if not record.is_secondary
        ]

        if primary_records:

            primary = (
                primary_records[
                    0
                ]
            )

        else:

            primary = (
                mapped_records[
                    0
                ]
            )

        # ----------------------------------------------------
        # NH
        # ----------------------------------------------------

        if primary.has_tag(
            "NH"
        ):

            nh = int(
                primary.get_tag(
                    "NH"
                )
            )

        else:

            nh = len(
                mapped_records
            )

        if nh == 1:

            status = "unique"

        else:

            status = "multi"

        # ----------------------------------------------------
        # Alignment details
        # ----------------------------------------------------

        alignments = []

        for record in mapped_records:

            if (
                record.reference_name
                is None
            ):
                continue

            alignments.append(
                {
                    "reference": (
                        record.reference_name
                    ),
                    "start": (
                        record.reference_start
                    ),
                    "end": (
                        record.reference_end
                    ),
                    "strand": (
                        "-"
                        if record.is_reverse
                        else "+"
                    ),
                    "mapq": (
                        record.mapping_quality
                    ),
                    "cigar": (
                        record.cigarstring
                    ),
                    "secondary": (
                        record.is_secondary
                    ),
                }
            )

        reads[
            read_id
        ] = {
            "status": (
                status
            ),
            "nh": (
                nh
            ),
            "sequence": (
                primary.query_sequence
            ),
            "alignments": (
                alignments
            ),
            "spliced": (
                any(
                    has_splice_cigar(
                        record.cigarstring
                    )
                    for record
                    in mapped_records
                )
            ),
        }

    return reads


# ============================================================
# RiboTxAlign BAM loader
# ============================================================

def load_ribotx_reads(
    bam_path: Path,
) -> dict:

    reads = {}

    with pysam.AlignmentFile(
        bam_path,
        "rb",
    ) as bam:

        for record in bam.fetch(
            until_eof=True,
        ):

            if record.is_secondary:
                continue

            if record.is_supplementary:
                continue

            read_id = (
                normalize_read_id(
                    record.query_name
                )
            )

            # ------------------------------------------------
            # Unmapped
            # ------------------------------------------------

            if record.is_unmapped:

                reads[
                    read_id
                ] = {
                    "status": (
                        "unmapped"
                    ),
                    "nh": 0,
                    "assignment_type": (
                        "genome_unmapped"
                    ),
                    "sequence": (
                        record.query_sequence
                    ),
                    "spliced": False,
                }

                continue

            # ------------------------------------------------
            # NH
            # ------------------------------------------------

            if record.has_tag(
                "NH"
            ):

                nh = int(
                    record.get_tag(
                        "NH"
                    )
                )

            else:

                nh = 1

            # ------------------------------------------------
            # XA assignment type
            # ------------------------------------------------

            if record.has_tag(
                "XA"
            ):

                assignment_type = str(
                    record.get_tag(
                        "XA"
                    )
                )

            else:

                assignment_type = (
                    "unknown"
                )

            # ------------------------------------------------
            # Mapping status
            # ------------------------------------------------

            if (
                assignment_type
                == "multimap_assigned"
            ):

                status = (
                    "multi_assigned"
                )

            else:

                status = (
                    "unique"
                )

            reads[
                read_id
            ] = {
                "status": (
                    status
                ),
                "nh": (
                    nh
                ),
                "assignment_type": (
                    assignment_type
                ),
                "sequence": (
                    record.query_sequence
                ),
                "spliced": (
                    has_splice_cigar(
                        record.cigarstring
                    )
                ),
            }

    return reads


# ============================================================
# Debug read IDs
# ============================================================

def print_read_id_examples(
    title: str,
    reads: dict,
    n: int = 5,
):

    print()
    print(
        title
    )

    for read_id in list(
        reads
    )[:n]:

        print(
            f"    {read_id}"
        )


# ============================================================
# Mapping summary
# ============================================================

def summarize_mapper_reads(
    mapper_name: str,
    reads: dict,
    read_ids,
    scope: str,
    splice_valid: bool = True,
):

    selected = [
        reads[
            read_id
        ]
        for read_id in read_ids
        if read_id in reads
    ]

    total_reads = len(
        selected
    )

    unmapped_reads = sum(
        info[
            "status"
        ]
        == "unmapped"
        for info in selected
    )

    unique_reads = sum(
        info[
            "status"
        ]
        == "unique"
        for info in selected
    )

    multi_reads = sum(
        info[
            "status"
        ]
        in {
            "multi",
            "multi_assigned",
        }
        for info in selected
    )

    mapped_reads = (
        unique_reads
        + multi_reads
    )

    if splice_valid:

        splice_reads = sum(
            (
                info[
                    "status"
                ]
                != "unmapped"
            )
            and info.get(
                "spliced",
                False,
            )
            for info in selected
        )

        splice_percent_of_mapped = (
            safe_percent(
                splice_reads,
                mapped_reads,
            )
        )

        splice_percent_of_total = (
            safe_percent(
                splice_reads,
                total_reads,
            )
        )

    else:

        splice_reads = None
        splice_percent_of_mapped = None
        splice_percent_of_total = None

    return {
        "scope": scope,
        "mapper": mapper_name,

        "total_reads": (
            total_reads
        ),

        "mapped_reads": (
            mapped_reads
        ),

        "unmapped_reads": (
            unmapped_reads
        ),

        "unique_reads": (
            unique_reads
        ),

        "multi_reads": (
            multi_reads
        ),

        "mapped_percent_of_total": (
            safe_percent(
                mapped_reads,
                total_reads,
            )
        ),

        "unmapped_percent_of_total": (
            safe_percent(
                unmapped_reads,
                total_reads,
            )
        ),

        "unique_percent_of_total": (
            safe_percent(
                unique_reads,
                total_reads,
            )
        ),

        "multi_percent_of_total": (
            safe_percent(
                multi_reads,
                total_reads,
            )
        ),

        "unique_percent_of_mapped": (
            safe_percent(
                unique_reads,
                mapped_reads,
            )
        ),

        "multi_percent_of_mapped": (
            safe_percent(
                multi_reads,
                mapped_reads,
            )
        ),

        "splice_reads": (
            splice_reads
        ),

        "splice_percent_of_mapped": (
            splice_percent_of_mapped
        ),

        "splice_percent_of_total": (
            splice_percent_of_total
        ),
    }


# ============================================================
# Multi-hit summary
# ============================================================

def summarize_multi_by_nh(
    mapper_name: str,
    reads: dict,
    read_ids,
    scope: str,
):

    counts = defaultdict(
        int
    )

    for read_id in read_ids:

        if read_id not in reads:
            continue

        info = (
            reads[
                read_id
            ]
        )

        if (
            info[
                "status"
            ]
            not in {
                "multi",
                "multi_assigned",
            }
        ):
            continue

        counts[
            int(
                info[
                    "nh"
                ]
            )
        ] += 1

    multi_total = sum(
        counts.values()
    )

    total_reads = len(
        read_ids
    )

    rows = []

    for nh in sorted(
        counts
    ):

        count = (
            counts[
                nh
            ]
        )

        rows.append(
            {
                "scope": (
                    scope
                ),
                "mapper": (
                    mapper_name
                ),
                "nh": (
                    nh
                ),
                "read_count": (
                    count
                ),
                "percent_of_multi": (
                    safe_percent(
                        count,
                        multi_total,
                    )
                ),
                "percent_of_total": (
                    safe_percent(
                        count,
                        total_reads,
                    )
                ),
            }
        )

    return rows


# ============================================================
# RiboTxAlign internal rescue
# ============================================================

def summarize_ribotx_rescue(
    mapper_name: str,
    ribotx_reads: dict,
    read_ids,
):

    selected = [
        ribotx_reads[
            read_id
        ]
        for read_id in read_ids
        if read_id in ribotx_reads
    ]

    contiguous_rescued = sum(
        info[
            "assignment_type"
        ]
        == "genome_rescued_unique"
        for info in selected
    )

    splice_rescued = sum(
        info[
            "assignment_type"
        ]
        == "genome_splice_rescued_unique"
        for info in selected
    )

    final_unmapped = sum(
        info[
            "status"
        ]
        == "unmapped"
        for info in selected
    )

    rescued_total = (
        contiguous_rescued
        + splice_rescued
    )

    rescue_input_reads = (
        rescued_total
        + final_unmapped
    )

    total_reads = len(
        selected
    )

    return {
        "mapper": (
            mapper_name
        ),

        "total_ribotx_reads": (
            total_reads
        ),

        "rescue_input_reads": (
            rescue_input_reads
        ),

        "contiguous_rescued_unique_reads": (
            contiguous_rescued
        ),

        "splice_rescued_unique_reads": (
            splice_rescued
        ),

        "rescued_unique_reads": (
            rescued_total
        ),

        "remaining_unmapped_reads": (
            final_unmapped
        ),

        "rescue_percent_of_rescue_input": (
            safe_percent(
                rescued_total,
                rescue_input_reads,
            )
        ),

        "contiguous_rescue_percent_of_rescue_input": (
            safe_percent(
                contiguous_rescued,
                rescue_input_reads,
            )
        ),

        "splice_rescue_percent_of_rescue_input": (
            safe_percent(
                splice_rescued,
                rescue_input_reads,
            )
        ),

        "rescue_percent_of_all_ribotx_reads": (
            safe_percent(
                rescued_total,
                total_reads,
            )
        ),
    }


# ============================================================
# Pairwise read dataframe
# ============================================================

def build_pairwise_dataframe(
    left_reads: dict,
    right_reads: dict,
    left_name: str,
    right_name: str,
):

    common_ids = (
        set(
            left_reads
        )
        & set(
            right_reads
        )
    )

    rows = []

    for read_id in sorted(
        common_ids
    ):

        left = (
            left_reads[
                read_id
            ]
        )

        right = (
            right_reads[
                read_id
            ]
        )

        rows.append(
            {
                "read_id": (
                    read_id
                ),

                f"{left_name}_status": (
                    left[
                        "status"
                    ]
                ),

                f"{left_name}_nh": (
                    left[
                        "nh"
                    ]
                ),

                f"{left_name}_nh_category": (
                    make_nh_category(
                        left[
                            "status"
                        ],
                        left[
                            "nh"
                        ],
                    )
                ),

                f"{left_name}_spliced": (
                    left.get(
                        "spliced",
                        False,
                    )
                ),

                f"{right_name}_status": (
                    right[
                        "status"
                    ]
                ),

                f"{right_name}_nh": (
                    right[
                        "nh"
                    ]
                ),

                f"{right_name}_nh_category": (
                    make_nh_category(
                        right[
                            "status"
                        ],
                        right[
                            "nh"
                        ],
                    )
                ),

                f"{right_name}_spliced": (
                    right.get(
                        "spliced",
                        False,
                    )
                ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# Transition matrices
# ============================================================

def save_transition_tables(
    df: pd.DataFrame,
    left_name: str,
    right_name: str,
    output_prefix: Path,
):

    if df.empty:

        print(
            f"WARNING: empty comparison: "
            f"{left_name} -> {right_name}"
        )

        return

    left_status = (
        f"{left_name}_status"
    )

    right_status = (
        f"{right_name}_status"
    )

    left_category = (
        f"{left_name}_nh_category"
    )

    right_category = (
        f"{right_name}_nh_category"
    )

    left_nh = (
        f"{left_name}_nh"
    )

    right_nh = (
        f"{right_name}_nh"
    )

    # --------------------------------------------------------
    # Status transition
    # --------------------------------------------------------

    transition = pd.crosstab(
        df[
            left_status
        ],
        df[
            right_status
        ],
    )

    transition.to_csv(
        Path(
            f"{output_prefix}."
            "transition.tsv"
        ),
        sep="\t",
    )

    transition_percent = (
        transition.div(
            transition.sum(
                axis=1
            ),
            axis=0,
        )
        * 100.0
    )

    transition_percent.to_csv(
        Path(
            f"{output_prefix}."
            "transition_percent.tsv"
        ),
        sep="\t",
        float_format="%.4f",
    )

    # --------------------------------------------------------
    # Detailed NH category transition
    # --------------------------------------------------------

    detailed = pd.crosstab(
        df[
            left_category
        ],
        df[
            right_category
        ],
    )

    detailed = (
        sort_mapping_transition(
            detailed
        )
    )

    detailed.to_csv(
        Path(
            f"{output_prefix}."
            "detailed_transition.tsv"
        ),
        sep="\t",
    )

    detailed_percent = (
        detailed.div(
            detailed.sum(
                axis=1
            ),
            axis=0,
        )
        * 100.0
    )

    detailed_percent.to_csv(
        Path(
            f"{output_prefix}."
            "detailed_transition_percent.tsv"
        ),
        sep="\t",
        float_format="%.4f",
    )

    # --------------------------------------------------------
    # Numeric NH transition
    # --------------------------------------------------------

    nh_transition = pd.crosstab(
        df[
            left_nh
        ],
        df[
            right_nh
        ],
    )

    nh_transition = (
        nh_transition
        .sort_index(
            axis=0
        )
        .sort_index(
            axis=1
        )
    )

    nh_transition.to_csv(
        Path(
            f"{output_prefix}."
            "nh_transition.tsv"
        ),
        sep="\t",
    )

    nh_transition_percent = (
        nh_transition.div(
            nh_transition.sum(
                axis=1
            ),
            axis=0,
        )
        * 100.0
    )

    nh_transition_percent.to_csv(
        Path(
            f"{output_prefix}."
            "nh_transition_percent.tsv"
        ),
        sep="\t",
        float_format="%.4f",
    )

    print()
    print(
        "========================================"
    )
    print(
        f"{left_name} -> "
        f"{right_name}"
    )
    print(
        "========================================"
    )

    print(
        transition
    )

    print()

    print(
        "Percent"
    )

    print(
        transition_percent.round(
            2
        )
    )


# ============================================================
# STAR -> RiboTx assignment
# ============================================================

def save_star_to_ribotx_assignment(
    star_reads: dict,
    ribotx_reads: dict,
    star_name: str,
    ribotx_name: str,
    output_prefix: Path,
):

    common_ids = (
        set(
            star_reads
        )
        & set(
            ribotx_reads
        )
    )

    rows = []

    for read_id in common_ids:

        star_info = (
            star_reads[
                read_id
            ]
        )

        ribotx_info = (
            ribotx_reads[
                read_id
            ]
        )

        rows.append(
            {
                "star_mapper": (
                    star_name
                ),

                "ribotx_mapper": (
                    ribotx_name
                ),

                "star_category": (
                    make_nh_category(
                        star_info[
                            "status"
                        ],
                        star_info[
                            "nh"
                        ],
                    )
                ),

                "ribotx_assignment_type": (
                    ribotx_info[
                        "assignment_type"
                    ]
                ),
            }
        )

    df = pd.DataFrame(
        rows
    )

    if df.empty:
        return

    transition = pd.crosstab(
        df[
            "star_category"
        ],
        df[
            "ribotx_assignment_type"
        ],
    )

    transition = (
        transition.reindex(
            index=sorted(
                transition.index,
                key=(
                    mapping_category_sort_key
                ),
            ),
            columns=sorted(
                transition.columns
            ),
        )
    )

    transition.to_csv(
        Path(
            f"{output_prefix}."
            "assignment.tsv"
        ),
        sep="\t",
    )

    transition_percent = (
        transition.div(
            transition.sum(
                axis=1
            ),
            axis=0,
        )
        * 100.0
    )

    transition_percent.to_csv(
        Path(
            f"{output_prefix}."
            "assignment_percent.tsv"
        ),
        sep="\t",
        float_format="%.4f",
    )


# ============================================================
# STAR multimapper -> RiboTx
# ============================================================

def save_star_multi_to_ribotx(
    star_reads: dict,
    ribotx_reads: dict,
    star_name: str,
    ribotx_name: str,
    output_path: Path,
):

    common_ids = (
        set(
            star_reads
        )
        & set(
            ribotx_reads
        )
    )

    rows = []

    for read_id in common_ids:

        star_info = (
            star_reads[
                read_id
            ]
        )

        if (
            star_info[
                "status"
            ]
            != "multi"
        ):
            continue

        ribotx_info = (
            ribotx_reads[
                read_id
            ]
        )

        rows.append(
            {
                "star_nh": (
                    star_info[
                        "nh"
                    ]
                ),

                "ribotx_status": (
                    ribotx_info[
                        "status"
                    ]
                ),

                "ribotx_nh": (
                    ribotx_info[
                        "nh"
                    ]
                ),

                "ribotx_assignment_type": (
                    ribotx_info[
                        "assignment_type"
                    ]
                ),
            }
        )

    df = pd.DataFrame(
        rows
    )

    if df.empty:

        df.to_csv(
            output_path,
            sep="\t",
            index=False,
        )

        return

    total = len(
        df
    )

    summary = (
        df.groupby(
            [
                "star_nh",
                "ribotx_status",
                "ribotx_nh",
                "ribotx_assignment_type",
            ],
            dropna=False,
        )
        .size()
        .rename(
            "read_count"
        )
        .reset_index()
    )

    summary[
        "percent_of_star_multi"
    ] = (
        summary[
            "read_count"
        ].map(
            lambda value: (
                safe_percent(
                    value,
                    total,
                )
            )
        )
    )

    summary = (
        summary.sort_values(
            [
                "star_nh",
                "ribotx_nh",
                "ribotx_status",
            ]
        )
        .reset_index(
            drop=True
        )
    )

    summary.insert(
        0,
        "ribotx_mapper",
        ribotx_name,
    )

    summary.insert(
        0,
        "star_mapper",
        star_name,
    )

    summary.to_csv(
        output_path,
        sep="\t",
        index=False,
        float_format="%.4f",
    )


# ============================================================
# Four-mapper read table
# ============================================================

def build_four_mapper_dataframe(
    mapper_dict: dict,
):

    all_ids = set()

    for reads in mapper_dict.values():
        all_ids |= set(
            reads
        )

    rows = []

    for read_id in sorted(
        all_ids
    ):

        row = {
            "read_id": (
                read_id
            ),
        }

        for (
            mapper_name,
            reads,
        ) in mapper_dict.items():

            if read_id not in reads:

                row[
                    f"{mapper_name}_present"
                ] = False

                row[
                    f"{mapper_name}_status"
                ] = "absent"

                row[
                    f"{mapper_name}_nh"
                ] = None

                row[
                    f"{mapper_name}_nh_category"
                ] = "absent"

                row[
                    f"{mapper_name}_spliced"
                ] = None

                if mapper_name.startswith(
                    "ribotx"
                ):

                    row[
                        f"{mapper_name}_assignment_type"
                    ] = "absent"

                continue

            info = (
                reads[
                    read_id
                ]
            )

            row[
                f"{mapper_name}_present"
            ] = True

            row[
                f"{mapper_name}_status"
            ] = (
                info[
                    "status"
                ]
            )

            row[
                f"{mapper_name}_nh"
            ] = (
                info[
                    "nh"
                ]
            )

            row[
                f"{mapper_name}_nh_category"
            ] = (
                make_nh_category(
                    info[
                        "status"
                    ],
                    info[
                        "nh"
                    ],
                )
            )

            row[
                f"{mapper_name}_spliced"
            ] = (
                info.get(
                    "spliced",
                    False,
                )
            )

            if mapper_name.startswith(
                "ribotx"
            ):

                row[
                    f"{mapper_name}_assignment_type"
                ] = (
                    info.get(
                        "assignment_type",
                        "unknown",
                    )
                )

        rows.append(
            row
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# RiboTx long-read vs short-read StringTie
# ============================================================

def save_ribotx_stringtie_comparison(
    longread_reads: dict,
    shortread_reads: dict,
    output_dir: Path,
):

    common_ids = (
        set(
            longread_reads
        )
        & set(
            shortread_reads
        )
    )

    rows = []

    for read_id in common_ids:

        long_info = (
            longread_reads[
                read_id
            ]
        )

        short_info = (
            shortread_reads[
                read_id
            ]
        )

        rows.append(
            {
                "read_id": (
                    read_id
                ),

                "longread_status": (
                    long_info[
                        "status"
                    ]
                ),

                "longread_nh": (
                    long_info[
                        "nh"
                    ]
                ),

                "longread_category": (
                    make_nh_category(
                        long_info[
                            "status"
                        ],
                        long_info[
                            "nh"
                        ],
                    )
                ),

                "longread_assignment_type": (
                    long_info[
                        "assignment_type"
                    ]
                ),

                "shortread_status": (
                    short_info[
                        "status"
                    ]
                ),

                "shortread_nh": (
                    short_info[
                        "nh"
                    ]
                ),

                "shortread_category": (
                    make_nh_category(
                        short_info[
                            "status"
                        ],
                        short_info[
                            "nh"
                        ],
                    )
                ),

                "shortread_assignment_type": (
                    short_info[
                        "assignment_type"
                    ]
                ),
            }
        )

    df = pd.DataFrame(
        rows
    )

    df.to_csv(
        (
            output_dir
            / "ribotx_longread_vs_shortread_reads.tsv.gz"
        ),
        sep="\t",
        index=False,
        compression="gzip",
    )

    if df.empty:
        return

    # --------------------------------------------------------
    # Assignment-type transition
    # --------------------------------------------------------

    assignment_transition = pd.crosstab(
        df[
            "longread_assignment_type"
        ],
        df[
            "shortread_assignment_type"
        ],
    )

    assignment_transition.to_csv(
        (
            output_dir
            / "ribotx_longread_to_shortread_assignment.tsv"
        ),
        sep="\t",
    )

    assignment_percent = (
        assignment_transition.div(
            assignment_transition.sum(
                axis=1
            ),
            axis=0,
        )
        * 100.0
    )

    assignment_percent.to_csv(
        (
            output_dir
            / "ribotx_longread_to_shortread_assignment_percent.tsv"
        ),
        sep="\t",
        float_format="%.4f",
    )

    # --------------------------------------------------------
    # Particularly useful summary:
    #
    # long-read multi -> short-read unique
    # short-read multi -> long-read unique
    # --------------------------------------------------------

    long_multi_short_unique = sum(
        (
            row[
                "longread_status"
            ]
            in {
                "multi",
                "multi_assigned",
            }
        )
        and (
            row[
                "shortread_status"
            ]
            == "unique"
        )
        for row in rows
    )

    short_multi_long_unique = sum(
        (
            row[
                "shortread_status"
            ]
            in {
                "multi",
                "multi_assigned",
            }
        )
        and (
            row[
                "longread_status"
            ]
            == "unique"
        )
        for row in rows
    )

    summary = pd.DataFrame(
        [
            {
                "comparison": (
                    "longread_multi_to_shortread_unique"
                ),
                "read_count": (
                    long_multi_short_unique
                ),
                "percent_of_common_reads": (
                    safe_percent(
                        long_multi_short_unique,
                        len(common_ids),
                    )
                ),
            },
            {
                "comparison": (
                    "shortread_multi_to_longread_unique"
                ),
                "read_count": (
                    short_multi_long_unique
                ),
                "percent_of_common_reads": (
                    safe_percent(
                        short_multi_long_unique,
                        len(common_ids),
                    )
                ),
            },
        ]
    )

    summary.to_csv(
        (
            output_dir
            / "ribotx_stringtie_unique_gain_summary.tsv"
        ),
        sep="\t",
        index=False,
        float_format="%.4f",
    )


# ============================================================
# Compare one sample
# ============================================================

def compare_sample(
    sample_info: dict,
):

    dataset_name = (
        sample_info[
            "dataset_name"
        ]
    )

    sample_name = (
        sample_info[
            "sample_name"
        ]
    )

    srr_id = (
        sample_info[
            "srr_id"
        ]
    )

    star_genome_bam = (
        sample_info[
            "star_genome_bam"
        ]
    )

    star_transcript_bam = (
        sample_info[
            "star_transcript_bam"
        ]
    )

    ribotx_longread_bam = (
        sample_info[
            "ribotx_longread_bam"
        ]
    )

    ribotx_shortread_bam = (
        sample_info[
            "ribotx_shortread_bam"
        ]
    )

    sample_output_dir = (
        OUTPUT_DIR
        / dataset_name
        / srr_id
    )

    sample_output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print()
    print(
        "============================================================"
    )

    print(
        f"Dataset : "
        f"{dataset_name}"
    )

    print(
        f"Sample  : "
        f"{sample_name}"
    )

    print(
        f"SRR     : "
        f"{srr_id}"
    )

    print(
        "============================================================"
    )

    # ========================================================
    # Load four BAMs
    # ========================================================

    print()
    print(
        "Loading STAR genome..."
    )

    star_genome_reads = (
        load_star_reads(
            star_genome_bam
        )
    )

    print(
        f"Reads: "
        f"{len(star_genome_reads):,}"
    )

    print()
    print(
        "Loading STAR transcript..."
    )

    print(
        f"BAM: {star_transcript_bam}"
    )

    star_transcript_reads = (
        load_star_reads(
            star_transcript_bam
        )
    )

    print(
        f"Reads: "
        f"{len(star_transcript_reads):,}"
    )

    print()
    print(
        "Loading RiboTxAlign long-read StringTie..."
    )

    ribotx_longread_reads = (
        load_ribotx_reads(
            ribotx_longread_bam
        )
    )

    print(
        f"Reads: "
        f"{len(ribotx_longread_reads):,}"
    )

    print()
    print(
        "Loading RiboTxAlign short-read StringTie..."
    )

    ribotx_shortread_reads = (
        load_ribotx_reads(
            ribotx_shortread_bam
        )
    )

    print(
        f"Reads: "
        f"{len(ribotx_shortread_reads):,}"
    )

    # ========================================================
    # Unified mapper dictionary
    # ========================================================

    mapper_dict = {
        "star_genome": (
            star_genome_reads
        ),
        "star_transcript": (
            star_transcript_reads
        ),
        "ribotx_longread_stringtie": (
            ribotx_longread_reads
        ),
        "ribotx_shortread_stringtie": (
            ribotx_shortread_reads
        ),
    }

    # ========================================================
    # Read ID examples
    # ========================================================

    print()
    print(
        "========================================"
    )
    print(
        "Normalized read ID examples"
    )
    print(
        "========================================"
    )

    for (
        mapper_name,
        reads,
    ) in mapper_dict.items():

        print_read_id_examples(
            f"{mapper_name}:",
            reads,
        )

    # ========================================================
    # Read ID overlap
    # ========================================================

    mapper_ids = {
        mapper_name: set(
            reads
        )
        for (
            mapper_name,
            reads,
        ) in mapper_dict.items()
    }

    common_all = set.intersection(
        *mapper_ids.values()
    )

    union_all = set.union(
        *mapper_ids.values()
    )

    print()
    print(
        "========================================"
    )
    print(
        "Read ID overlap"
    )
    print(
        "========================================"
    )

    for (
        mapper_name,
        ids,
    ) in mapper_ids.items():

        print(
            f"{mapper_name:32s}: "
            f"{len(ids):,}"
        )

    print()

    print(
        f"{'Union of all':32s}: "
        f"{len(union_all):,}"
    )

    print(
        f"{'Common to all':32s}: "
        f"{len(common_all):,}"
    )

    # Pairwise overlap
    overlap_rows = []

    for (
        left_name,
        right_name,
    ) in combinations(
        mapper_dict.keys(),
        2,
    ):

        common = (
            mapper_ids[
                left_name
            ]
            & mapper_ids[
                right_name
            ]
        )

        overlap_rows.append(
            {
                "left_mapper": (
                    left_name
                ),
                "right_mapper": (
                    right_name
                ),
                "left_reads": (
                    len(
                        mapper_ids[
                            left_name
                        ]
                    )
                ),
                "right_reads": (
                    len(
                        mapper_ids[
                            right_name
                        ]
                    )
                ),
                "common_reads": (
                    len(
                        common
                    )
                ),
            }
        )

    pd.DataFrame(
        overlap_rows
    ).to_csv(
        (
            sample_output_dir
            / "read_id_overlap.tsv"
        ),
        sep="\t",
        index=False,
    )

    # ========================================================
    # Mapping summary
    # ========================================================

    summary_rows = []

    mapper_configs = [
        (
            "STAR_genome",
            star_genome_reads,
            True,
        ),
        (
            "STAR_transcript",
            star_transcript_reads,
            False,
        ),
        (
            "RiboTxAlign_longread_stringtie",
            ribotx_longread_reads,
            True,
        ),
        (
            "RiboTxAlign_shortread_stringtie",
            ribotx_shortread_reads,
            True,
        ),
    ]

    for (
        mapper_name,
        reads,
        splice_valid,
    ) in mapper_configs:

        # ----------------------------------------------------
        # All reads present in that mapper
        # ----------------------------------------------------

        summary_rows.append(
            summarize_mapper_reads(
                mapper_name=(
                    mapper_name
                ),
                reads=(
                    reads
                ),
                read_ids=(
                    set(
                        reads
                    )
                ),
                scope=(
                    "all_mapper_reads"
                ),
                splice_valid=(
                    splice_valid
                ),
            )
        )

        # ----------------------------------------------------
        # Reads common to all four methods
        # ----------------------------------------------------

        summary_rows.append(
            summarize_mapper_reads(
                mapper_name=(
                    mapper_name
                ),
                reads=(
                    reads
                ),
                read_ids=(
                    common_all
                ),
                scope=(
                    "common_to_all"
                ),
                splice_valid=(
                    splice_valid
                ),
            )
        )

    mapping_summary = pd.DataFrame(
        summary_rows
    )

    mapping_summary.to_csv(
        (
            sample_output_dir
            / "mapping_summary.tsv"
        ),
        sep="\t",
        index=False,
        float_format="%.4f",
    )

    print()
    print(
        "========================================"
    )
    print(
        "Mapping summary"
    )
    print(
        "========================================"
    )

    print(
        mapping_summary.to_string(
            index=False
        )
    )

    # ========================================================
    # NH distribution
    # ========================================================

    multi_rows = []

    for (
        mapper_name,
        reads,
        _,
    ) in mapper_configs:

        # All mapper reads
        multi_rows.extend(
            summarize_multi_by_nh(
                mapper_name=(
                    mapper_name
                ),
                reads=(
                    reads
                ),
                read_ids=(
                    set(
                        reads
                    )
                ),
                scope=(
                    "all_mapper_reads"
                ),
            )
        )

        # Common-to-all reads
        multi_rows.extend(
            summarize_multi_by_nh(
                mapper_name=(
                    mapper_name
                ),
                reads=(
                    reads
                ),
                read_ids=(
                    common_all
                ),
                scope=(
                    "common_to_all"
                ),
            )
        )

    multi_summary = pd.DataFrame(
        multi_rows
    )

    if not multi_summary.empty:

        multi_summary = (
            multi_summary.sort_values(
                [
                    "scope",
                    "mapper",
                    "nh",
                ]
            )
            .reset_index(
                drop=True
            )
        )

    multi_summary.to_csv(
        (
            sample_output_dir
            / "multi_reads_by_nh.tsv"
        ),
        sep="\t",
        index=False,
        float_format="%.4f",
    )

    # ========================================================
    # RiboTxAlign internal rescue
    # ========================================================

    rescue_rows = [
        summarize_ribotx_rescue(
            mapper_name=(
                "RiboTxAlign_longread_stringtie"
            ),
            ribotx_reads=(
                ribotx_longread_reads
            ),
            read_ids=(
                set(
                    ribotx_longread_reads
                )
            ),
        ),

        summarize_ribotx_rescue(
            mapper_name=(
                "RiboTxAlign_shortread_stringtie"
            ),
            ribotx_reads=(
                ribotx_shortread_reads
            ),
            read_ids=(
                set(
                    ribotx_shortread_reads
                )
            ),
        ),
    ]

    rescue_summary = pd.DataFrame(
        rescue_rows
    )

    rescue_summary.to_csv(
        (
            sample_output_dir
            / "ribotx_rescue_summary.tsv"
        ),
        sep="\t",
        index=False,
        float_format="%.4f",
    )

    print()
    print(
        "========================================"
    )
    print(
        "RiboTxAlign internal rescue"
    )
    print(
        "========================================"
    )

    print(
        rescue_summary.to_string(
            index=False
        )
    )

    # ========================================================
    # Pairwise comparisons
    #
    # 4 methods -> 6 comparisons
    # ========================================================

    for (
        left_name,
        right_name,
    ) in combinations(
        mapper_dict.keys(),
        2,
    ):

        left_reads = (
            mapper_dict[
                left_name
            ]
        )

        right_reads = (
            mapper_dict[
                right_name
            ]
        )

        pair_df = (
            build_pairwise_dataframe(
                left_reads=(
                    left_reads
                ),
                right_reads=(
                    right_reads
                ),
                left_name=(
                    left_name
                ),
                right_name=(
                    right_name
                ),
            )
        )

        prefix = (
            sample_output_dir
            / (
                f"{left_name}"
                f"_to_"
                f"{right_name}"
            )
        )

        pair_df.to_csv(
            Path(
                f"{prefix}."
                "paired_reads.tsv.gz"
            ),
            sep="\t",
            index=False,
            compression="gzip",
        )

        save_transition_tables(
            df=(
                pair_df
            ),
            left_name=(
                left_name
            ),
            right_name=(
                right_name
            ),
            output_prefix=(
                prefix
            ),
        )

    # ========================================================
    # STAR -> RiboTx assignment tables
    #
    # 2 STAR mappings x 2 RiboTxAlign mappings
    # ========================================================

    star_mappers = {
        "star_genome": (
            star_genome_reads
        ),
        "star_transcript": (
            star_transcript_reads
        ),
    }

    ribotx_mappers = {
        "ribotx_longread_stringtie": (
            ribotx_longread_reads
        ),
        "ribotx_shortread_stringtie": (
            ribotx_shortread_reads
        ),
    }

    for (
        star_name,
        star_reads,
    ) in star_mappers.items():

        for (
            ribotx_name,
            ribotx_reads,
        ) in ribotx_mappers.items():

            prefix = (
                sample_output_dir
                / (
                    f"{star_name}"
                    f"_to_"
                    f"{ribotx_name}"
                )
            )

            save_star_to_ribotx_assignment(
                star_reads=(
                    star_reads
                ),
                ribotx_reads=(
                    ribotx_reads
                ),
                star_name=(
                    star_name
                ),
                ribotx_name=(
                    ribotx_name
                ),
                output_prefix=(
                    prefix
                ),
            )

            save_star_multi_to_ribotx(
                star_reads=(
                    star_reads
                ),
                ribotx_reads=(
                    ribotx_reads
                ),
                star_name=(
                    star_name
                ),
                ribotx_name=(
                    ribotx_name
                ),
                output_path=(
                    sample_output_dir
                    / (
                        f"{star_name}"
                        f"_multi_to_"
                        f"{ribotx_name}.tsv"
                    )
                ),
            )

    # ========================================================
    # Direct comparison:
    #
    # RiboTxAlign long-read StringTie
    # vs
    # RiboTxAlign short-read StringTie
    # ========================================================

    save_ribotx_stringtie_comparison(
        longread_reads=(
            ribotx_longread_reads
        ),
        shortread_reads=(
            ribotx_shortread_reads
        ),
        output_dir=(
            sample_output_dir
        ),
    )

    # ========================================================
    # Four-mapper read-level table
    # ========================================================

    # four_mapper_df = (
    #     build_four_mapper_dataframe(
    #         mapper_dict
    #     )
    # )
    #
    # four_mapper_df.to_csv(
    #     (
    #         sample_output_dir
    #         / "four_mapper_reads.tsv.gz"
    #     ),
    #     sep="\t",
    #     index=False,
    #     compression="gzip",
    # )

    print()
    print(
        "========================================"
    )
    print(
        "Finished"
    )
    print(
        "========================================"
    )

    print(
        f"{dataset_name} / "
        f"{srr_id}"
    )

    print(
        f"Output: "
        f"{sample_output_dir}"
    )


# ============================================================
# Main
# ============================================================

def main():

    all_samples = []

    # ========================================================
    # Discover BAM files
    # ========================================================

    for (
        dataset_name,
        config,
    ) in DATASETS.items():

        print()
        print(
            "############################################################"
        )
        print(
            f"Finding BAMs: "
            f"{dataset_name}"
        )
        print(
            "############################################################"
        )

        for srr_id in config[
            "samples"
        ]:

            single_config = {
                **config,
                "samples": [
                    srr_id
                ],
            }

            try:

                samples = (
                    find_sample_bams(
                        dataset_name=(
                            dataset_name
                        ),
                        config=(
                            single_config
                        ),
                    )
                )

                all_samples.extend(
                    samples
                )

            except Exception as exc:

                print()
                print(
                    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                )

                print(
                    f"BAM discovery failed: "
                    f"{dataset_name} / {srr_id}"
                )

                print(
                    exc
                )

                print(
                    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                )

    # ========================================================
    # Resolved BAM files
    # ========================================================

    print()
    print(
        "############################################################"
    )
    print(
        "Resolved samples"
    )
    print(
        "############################################################"
    )

    for sample_info in all_samples:

        print()

        print(
            f"{sample_info['dataset_name']} "
            f"{sample_info['srr_id']}"
        )

        print(
            f"  STAR genome                  : "
            f"{sample_info['star_genome_bam']}"
        )

        print(
            f"  STAR transcript              : "
            f"{sample_info['star_transcript_bam']}"
        )

        print(
            f"  RiboTx long-read StringTie   : "
            f"{sample_info['ribotx_longread_bam']}"
        )

        print(
            f"  RiboTx short-read StringTie  : "
            f"{sample_info['ribotx_shortread_bam']}"
        )

    print()

    print(
        f"Resolved samples: "
        f"{len(all_samples)}"
    )

    # ========================================================
    # Run comparisons
    # ========================================================

    for sample_info in all_samples:

        try:

            compare_sample(
                sample_info
            )

        except Exception as exc:

            print()
            print(
                "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            )

            print(
                f"Comparison failed: "
                f"{sample_info['dataset_name']} / "
                f"{sample_info['srr_id']}"
            )

            print(
                exc
            )

            print(
                "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            )

    print()
    print(
        "============================================================"
    )
    print(
        "All comparisons finished"
    )
    print(
        "============================================================"
    )

    print(
        f"Output: "
        f"{OUTPUT_DIR}"
    )


if __name__ == "__main__":
    main()