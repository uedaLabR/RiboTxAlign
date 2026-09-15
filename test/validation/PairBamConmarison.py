from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam


# ============================================================
# Directories
# ============================================================

BASE = Path(
    "/mnt/share/project/riboML/validation"
)

STAR_DIR = (
    BASE
    / "star_validation"
)

RIBOTX_DIR = (
    BASE
    / "ribotx_validation"
)

OUTPUT_DIR = (
    BASE
    / "bam_comparison"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Parameters
# ============================================================

N_EXAMPLES = 20


def has_splice_cigar(
    cigar,
) -> bool:
    """Return True when a CIGAR contains an N (splice) operation."""

    if cigar is None:
        return False

    return "N" in cigar


def make_nh_category(
    status: str,
    nh: int,
) -> str:
    """Create a detailed mapping category including hit multiplicity."""

    if status == "unmapped" or nh == 0:
        return "unmapped"

    if status == "unique":
        return f"unique_NH{nh}"

    return f"multi_NH{nh}"


def mapping_category_sort_key(
    category: str,
) -> tuple[int, int]:
    """Sort mapping categories by status and numeric NH."""

    if category == "unmapped":
        return 0, 0

    status, nh_text = category.rsplit(
        "_NH",
        1,
    )
    status_order = {
        "unique": 1,
        "multi": 2,
    }

    return (
        status_order.get(status, 99),
        int(nh_text),
    )


def sort_mapping_transition(
    table: pd.DataFrame,
) -> pd.DataFrame:
    """Sort both axes of a mapping-category transition matrix."""

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
# Read ID normalization
# ============================================================

def normalize_star_read_id(
    read_id: str,
) -> str:
    """
    Normalize STAR read ID.

    Example:
        SRR26942051.34957075_CGACTTAGC
        ->
        SRR26942051.34957075
    """

    return read_id.rsplit(
        "_",
        1,
    )[0]


def normalize_ribotx_read_id(
    read_id: str,
) -> str:
    """
    Normalize RiboTxAlign read ID.

    Example:
        SRR26942051.29119401|umi=CCGTCGCCC|family_size=1
        ->
        SRR26942051.29119401
    """

    return read_id.split(
        "|",
        1,
    )[0]


def write_star_unique_ribotx_unmapped_fastq(
    sample_name: str,
    star_bam: Path,
    star_reads: dict,
    ribotx_reads: dict,
):
    """
    Write FASTQ reads satisfying:

        1. Read ID exists in both STAR and RiboTxAlign BAMs
        2. STAR mapping is unique
        3. RiboTxAlign is unmapped

    STAR mapping information is added to the FASTQ header.
    """

    # --------------------------------------------------------
    # Select target read IDs
    # --------------------------------------------------------

    common_ids = (
        set(star_reads)
        & set(ribotx_reads)
    )

    target_ids = {
        read_id
        for read_id in common_ids
        if (
            star_reads[read_id]["status"]
            == "unique"
        )
        and (
            ribotx_reads[read_id]["status"]
            == "unmapped"
        )
    }

    print()
    print(
        "========================================"
    )
    print(
        "STAR unique -> RiboTxAlign unmapped"
    )
    print(
        "========================================"
    )

    print(
        f"Reads: {len(target_ids):,}"
    )

    # --------------------------------------------------------
    # Output FASTQ
    # --------------------------------------------------------

    output_fastq = (
        OUTPUT_DIR
        / (
            f"{sample_name}."
            "STAR_unique_RiboTx_unmapped.fastq"
        )
    )

    written_ids = set()

    with pysam.AlignmentFile(
        star_bam,
        "rb",
    ) as bam, open(
        output_fastq,
        "w",
    ) as out:

        for record in bam.fetch(
            until_eof=True,
        ):

            if record.is_unmapped:
                continue

            if record.is_secondary:
                continue

            if record.is_supplementary:
                continue

            read_id = normalize_star_read_id(
                record.query_name
            )

            if read_id not in target_ids:
                continue

            # Prevent duplicate output
            if read_id in written_ids:
                continue

            # ------------------------------------------------
            # Sequence
            # ------------------------------------------------

            sequence = (
                record.query_sequence
            )

            if sequence is None:
                continue

            # ------------------------------------------------
            # Quality
            # ------------------------------------------------

            if (
                record.query_qualities
                is not None
            ):

                quality = (
                    pysam.array_to_qualitystring(
                        record.query_qualities
                    )
                )

            else:

                # Fallback if BAM contains no base qualities
                quality = (
                    "I"
                    * len(sequence)
                )

            # ------------------------------------------------
            # STAR mapping information
            # ------------------------------------------------

            chrom = (
                record.reference_name
            )

            start = (
                record.reference_start
                + 1
            )

            end = (
                record.reference_end
            )

            strand = (
                "-"
                if record.is_reverse
                else "+"
            )

            cigar = (
                record.cigarstring
            )

            mapq = (
                record.mapping_quality
            )

            # ------------------------------------------------
            # FASTQ header
            # ------------------------------------------------

            header = (
                f"@{read_id}"
                f"|STAR={chrom}:{start}-{end}:{strand}"
                f"|CIGAR={cigar}"
                f"|MAPQ={mapq}"
            )

            # ------------------------------------------------
            # Write FASTQ
            # ------------------------------------------------

            out.write(
                header
                + "\n"
            )

            out.write(
                sequence
                + "\n"
            )

            out.write(
                "+\n"
            )

            out.write(
                quality
                + "\n"
            )

            written_ids.add(
                read_id
            )

    print(
        f"FASTQ written: "
        f"{len(written_ids):,}"
    )

    print(
        f"Output: "
        f"{output_fastq}"
    )

    return output_fastq

# ============================================================
# STAR
# ============================================================

def load_star_reads(
    bam_path: Path,
) -> dict:
    """
    Load STAR BAM.

    Returns one dictionary entry per normalized read ID.

    Stored information:
        status
        nh
        sequence
        alignments
    """

    records_by_read = defaultdict(list)

    with pysam.AlignmentFile(
        bam_path,
        "rb",
    ) as bam:

        for record in bam.fetch(
            until_eof=True,
        ):

            if record.is_supplementary:
                continue

            read_id = normalize_star_read_id(
                record.query_name
            )

            records_by_read[
                read_id
            ].append(record)

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
                    records[0].query_sequence
                )

            reads[read_id] = {
                "status": "unmapped",
                "nh": 0,
                "sequence": sequence,
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
            primary_record = (
                primary_records[0]
            )
        else:
            primary_record = (
                mapped_records[0]
            )

        # ----------------------------------------------------
        # NH
        # ----------------------------------------------------

        if primary_record.has_tag(
            "NH"
        ):

            nh = int(
                primary_record.get_tag(
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
        # Alignment information
        # ----------------------------------------------------

        alignment_list = []

        for record in mapped_records:

            if record.reference_name is None:
                continue

            alignment_list.append(
                {
                    "chrom": (
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

        reads[read_id] = {
            "status": status,
            "nh": nh,
            "sequence": (
                primary_record.query_sequence
            ),
            "alignments": alignment_list,
            "spliced": has_splice_cigar(
                primary_record.cigarstring
            ),
        }

    return reads


# ============================================================
# RiboTxAlign
# ============================================================

def load_ribotx_reads(
    bam_path: Path,
) -> dict:
    """
    Load RiboTxAlign genome BAM.
    """

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
                normalize_ribotx_read_id(
                    record.query_name
                )
            )

            # ------------------------------------------------
            # Unmapped
            # ------------------------------------------------

            if record.is_unmapped:

                reads[read_id] = {
                    "status": "unmapped",
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
            # Assignment type
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

            if (
                assignment_type
                == "multimap_assigned"
            ):

                status = (
                    "multi_assigned"
                )

            else:

                status = "unique"

            reads[read_id] = {
                "status": status,
                "nh": nh,
                "assignment_type": (
                    assignment_type
                ),
                "sequence": (
                    record.query_sequence
                ),
                "spliced": has_splice_cigar(
                    record.cigarstring
                ),
            }

    return reads


# ============================================================
# Find BAM pairs
# ============================================================

def find_bam_pairs():

    pairs = []

    # --------------------------------------------------------
    # HEK293T
    # --------------------------------------------------------

    for star_bam in sorted(
        STAR_DIR.glob(
            "HEK293T_*.deduplicated.bam"
        )
    ):

        sample_name = (
            star_bam.name
            .removesuffix(
                ".deduplicated.bam"
            )
        )

        ribotx_bam = (
            RIBOTX_DIR
            / sample_name
            / "bam"
            / (
                f"{sample_name}."
                "genome.sorted.bam"
            )
        )

        if not ribotx_bam.exists():

            print(
                f"WARNING: "
                f"RiboTxAlign BAM not found: "
                f"{ribotx_bam}"
            )

            continue

        pairs.append(
            (
                sample_name,
                star_bam,
                ribotx_bam,
            )
        )

    # --------------------------------------------------------
    # Hela / mouse
    # --------------------------------------------------------

    for dataset in [
        "Hela",
        "mouse_brain",
        "mouse_liver",
    ]:

        pattern = (
            f"{dataset}_*."
            "Aligned.sortedByCoord.out.bam"
        )

        for star_bam in sorted(
            STAR_DIR.glob(
                pattern
            )
        ):

            sample_name = (
                star_bam.name
                .removesuffix(
                    ".Aligned.sortedByCoord.out.bam"
                )
            )

            ribotx_bam = (
                RIBOTX_DIR
                / sample_name
                / "bam"
                / (
                    f"{sample_name}."
                    "genome.sorted.bam"
                )
            )

            if not ribotx_bam.exists():

                print(
                    f"WARNING: "
                    f"RiboTxAlign BAM not found: "
                    f"{ribotx_bam}"
                )

                continue

            pairs.append(
                (
                    sample_name,
                    star_bam,
                    ribotx_bam,
                )
            )

    return pairs


# ============================================================
# Convert STAR read information to rows
# ============================================================

def make_star_detail_rows(
    read_ids,
    star_reads,
    category,
):
    """
    Convert selected STAR reads into a flat DataFrame.
    """

    rows = []

    for read_id in read_ids:

        info = star_reads[
            read_id
        ]

        alignments = (
            info["alignments"]
        )

        # If no alignment exists, still make one row.
        if not alignments:

            rows.append(
                {
                    "category": category,
                    "read_id": read_id,
                    "star_status": (
                        info["status"]
                    ),
                    "star_nh": (
                        info["nh"]
                    ),
                    "sequence": (
                        info["sequence"]
                    ),
                    "read_length": (
                        len(info["sequence"])
                        if info["sequence"]
                        else 0
                    ),
                    "chrom": None,
                    "start": None,
                    "end": None,
                    "strand": None,
                    "mapq": None,
                    "cigar": None,
                    "secondary": None,
                }
            )

            continue

        for alignment in alignments:

            rows.append(
                {
                    "category": category,
                    "read_id": read_id,
                    "star_status": (
                        info["status"]
                    ),
                    "star_nh": (
                        info["nh"]
                    ),
                    "sequence": (
                        info["sequence"]
                    ),
                    "read_length": (
                        len(info["sequence"])
                        if info["sequence"]
                        else 0
                    ),
                    "chrom": (
                        alignment["chrom"]
                    ),
                    "start": (
                        alignment["start"]
                    ),
                    "end": (
                        alignment["end"]
                    ),
                    "strand": (
                        alignment["strand"]
                    ),
                    "mapq": (
                        alignment["mapq"]
                    ),
                    "cigar": (
                        alignment["cigar"]
                    ),
                    "secondary": (
                        alignment[
                            "secondary"
                        ]
                    ),
                }
            )

    return rows


# ============================================================
# Print example reads
# ============================================================

def print_examples(
    title,
    read_ids,
    star_reads,
    n_examples=N_EXAMPLES,
):

    print()
    print(
        "============================================================"
    )
    print(
        title
    )
    print(
        "============================================================"
    )

    print(
        f"Total reads: "
        f"{len(read_ids):,}"
    )

    selected_ids = list(
        sorted(read_ids)
    )[
        :n_examples
    ]

    print(
        f"Showing first "
        f"{len(selected_ids)} reads"
    )

    for read_id in selected_ids:

        info = star_reads[
            read_id
        ]

        print()
        print(
            "------------------------------------------------------------"
        )

        print(
            f"Read ID : {read_id}"
        )

        print(
            f"Status  : {info['status']}"
        )

        print(
            f"STAR NH : {info['nh']}"
        )

        print(
            f"Length  : "
            f"{len(info['sequence']) if info['sequence'] else 0}"
        )

        print(
            f"Sequence: {info['sequence']}"
        )

        for i, aln in enumerate(
            info["alignments"],
            start=1,
        ):

            print(
                f"  Alignment {i}: "
                f"{aln['chrom']}:"
                f"{aln['start']}-"
                f"{aln['end']} "
                f"{aln['strand']} "
                f"CIGAR={aln['cigar']} "
                f"MAPQ={aln['mapq']} "
                f"secondary="
                f"{aln['secondary']}"
            )


# ============================================================
# Summary helpers
# ============================================================

def safe_percent(
    numerator: int,
    denominator: int,
) -> float:
    """Return a percentage while avoiding division by zero."""

    if denominator == 0:
        return 0.0

    return numerator / denominator * 100.0


def summarize_mapper_reads(
    mapper_name: str,
    reads: dict,
    read_ids,
    scope: str,
) -> dict:
    """Summarize unmapped, unique, multi, and splice reads."""

    selected = [
        reads[read_id]
        for read_id in read_ids
    ]

    total_reads = len(selected)
    unmapped_reads = sum(
        info["status"] == "unmapped"
        for info in selected
    )
    unique_reads = sum(
        info["status"] == "unique"
        for info in selected
    )
    multi_reads = sum(
        info["status"] in {
            "multi",
            "multi_assigned",
        }
        for info in selected
    )
    mapped_reads = unique_reads + multi_reads
    splice_reads = sum(
        info["status"] != "unmapped"
        and info.get("spliced", False)
        for info in selected
    )

    return {
        "scope": scope,
        "mapper": mapper_name,
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "unmapped_reads": unmapped_reads,
        "unique_reads": unique_reads,
        "multi_reads": multi_reads,
        "splice_reads": splice_reads,
        "mapped_percent_of_total": safe_percent(
            mapped_reads,
            total_reads,
        ),
        "unmapped_percent_of_total": safe_percent(
            unmapped_reads,
            total_reads,
        ),
        "unique_percent_of_total": safe_percent(
            unique_reads,
            total_reads,
        ),
        "multi_percent_of_total": safe_percent(
            multi_reads,
            total_reads,
        ),
        "unique_percent_of_mapped": safe_percent(
            unique_reads,
            mapped_reads,
        ),
        "multi_percent_of_mapped": safe_percent(
            multi_reads,
            mapped_reads,
        ),
        "splice_percent_of_mapped": safe_percent(
            splice_reads,
            mapped_reads,
        ),
        "splice_percent_of_total": safe_percent(
            splice_reads,
            total_reads,
        ),
    }


def summarize_multi_by_nh(
    mapper_name: str,
    reads: dict,
    read_ids,
    scope: str,
) -> list[dict]:
    """Count multimapping reads separately for every NH value."""

    counts = defaultdict(int)

    for read_id in read_ids:
        info = reads[read_id]
        if info["status"] not in {
            "multi",
            "multi_assigned",
        }:
            continue

        counts[int(info["nh"])] += 1

    multi_total = sum(counts.values())
    total_reads = len(read_ids)

    return [
        {
            "scope": scope,
            "mapper": mapper_name,
            "nh": nh,
            "read_count": read_count,
            "percent_of_multi": safe_percent(
                read_count,
                multi_total,
            ),
            "percent_of_total": safe_percent(
                read_count,
                total_reads,
            ),
        }
        for nh, read_count in sorted(counts.items())
    ]


def summarize_ribotx_rescue(
    ribotx_reads: dict,
    read_ids,
) -> dict:
    """Summarize RiboTxAlign's internal unmapped-read rescue stage."""

    selected = [
        ribotx_reads[read_id]
        for read_id in read_ids
    ]

    contiguous_rescued = sum(
        info["assignment_type"]
        == "genome_rescued_unique"
        for info in selected
    )
    splice_rescued = sum(
        info["assignment_type"]
        == "genome_splice_rescued_unique"
        for info in selected
    )
    final_unmapped = sum(
        info["status"] == "unmapped"
        for info in selected
    )

    rescued_total = (
        contiguous_rescued
        + splice_rescued
    )

    # Every final unmapped read and every genome-rescued read entered
    # RiboTxAlign's internal rescue stage.
    rescue_input_reads = (
        rescued_total
        + final_unmapped
    )
    total_reads = len(selected)

    return {
        "total_ribotx_reads": total_reads,
        "rescue_input_reads": rescue_input_reads,
        "contiguous_rescued_unique_reads": contiguous_rescued,
        "splice_rescued_unique_reads": splice_rescued,
        "rescued_unique_reads": rescued_total,
        "remaining_unmapped_reads": final_unmapped,
        "rescue_percent_of_rescue_input": safe_percent(
            rescued_total,
            rescue_input_reads,
        ),
        "contiguous_rescue_percent_of_rescue_input": safe_percent(
            contiguous_rescued,
            rescue_input_reads,
        ),
        "splice_rescue_percent_of_rescue_input": safe_percent(
            splice_rescued,
            rescue_input_reads,
        ),
        "rescue_percent_of_all_ribotx_reads": safe_percent(
            rescued_total,
            total_reads,
        ),
    }


# ============================================================
# Compare one sample
# ============================================================

def compare_sample(
    sample_name: str,
    star_bam: Path,
    ribotx_bam: Path,
):

    print()
    print(
        "============================================================"
    )
    print(
        f"Sample: {sample_name}"
    )
    print(
        "============================================================"
    )

    print(
        f"STAR   : {star_bam}"
    )

    print(
        f"RiboTx : {ribotx_bam}"
    )

    # --------------------------------------------------------
    # Load
    # --------------------------------------------------------

    print()
    print(
        "Loading STAR..."
    )

    star_reads = load_star_reads(
        star_bam
    )

    print(
        f"STAR reads: "
        f"{len(star_reads):,}"
    )

    print(
        "Loading RiboTxAlign..."
    )

    ribotx_reads = load_ribotx_reads(
        ribotx_bam
    )

    write_star_unique_ribotx_unmapped_fastq(
        sample_name=sample_name,
        star_bam=star_bam,
        star_reads=star_reads,
        ribotx_reads=ribotx_reads,
    )

    print(
        f"RiboTxAlign reads: "
        f"{len(ribotx_reads):,}"
    )

    # --------------------------------------------------------
    # Read IDs
    # --------------------------------------------------------

    star_ids = set(
        star_reads
    )

    ribotx_ids = set(
        ribotx_reads
    )

    common_ids = (
        star_ids
        & ribotx_ids
    )

    star_only_ids = (
        star_ids
        - ribotx_ids
    )

    ribotx_only_ids = (
        ribotx_ids
        - star_ids
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

    print(
        f"Common reads : "
        f"{len(common_ids):,}"
    )

    print(
        f"STAR only    : "
        f"{len(star_only_ids):,}"
    )

    print(
        f"RiboTx only  : "
        f"{len(ribotx_only_ids):,}"
    )

    # ========================================================
    # Common-read comparison
    # ========================================================

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
                "read_id": read_id,
                "star_status": (
                    star_info["status"]
                ),
                "star_nh": (
                    star_info["nh"]
                ),
                "star_nh_category": (
                    make_nh_category(
                        star_info["status"],
                        star_info["nh"],
                    )
                ),
                "star_spliced": (
                    star_info["spliced"]
                ),
                "ribotx_status": (
                    ribotx_info["status"]
                ),
                "ribotx_nh": (
                    ribotx_info["nh"]
                ),
                "ribotx_nh_category": (
                    make_nh_category(
                        ribotx_info["status"],
                        ribotx_info["nh"],
                    )
                ),
                "ribotx_spliced": (
                    ribotx_info["spliced"]
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

    paired_output = (
        OUTPUT_DIR
        / (
            f"{sample_name}."
            "paired_reads.tsv.gz"
        )
    )

    df.to_csv(
        paired_output,
        sep="\t",
        index=False,
        compression="gzip",
    )

    # ========================================================
    # Mapping, multi-hit, splice, and rescue summaries
    # ========================================================

    mapping_summary_rows = []

    for scope, selected_star_ids, selected_ribotx_ids in [
        (
            "all_mapper_reads",
            star_ids,
            ribotx_ids,
        ),
        (
            "common_reads",
            common_ids,
            common_ids,
        ),
    ]:
        mapping_summary_rows.append(
            summarize_mapper_reads(
                mapper_name="STAR",
                reads=star_reads,
                read_ids=selected_star_ids,
                scope=scope,
            )
        )
        mapping_summary_rows.append(
            summarize_mapper_reads(
                mapper_name="RiboTxAlign",
                reads=ribotx_reads,
                read_ids=selected_ribotx_ids,
                scope=scope,
            )
        )

    mapping_summary = pd.DataFrame(
        mapping_summary_rows
    )
    mapping_summary_output = (
        OUTPUT_DIR
        / f"{sample_name}.mapping_summary.tsv"
    )
    mapping_summary.to_csv(
        mapping_summary_output,
        sep="\t",
        index=False,
        float_format="%.4f",
    )

    multi_by_nh_rows = []

    for scope, selected_star_ids, selected_ribotx_ids in [
        (
            "all_mapper_reads",
            star_ids,
            ribotx_ids,
        ),
        (
            "common_reads",
            common_ids,
            common_ids,
        ),
    ]:
        multi_by_nh_rows.extend(
            summarize_multi_by_nh(
                mapper_name="STAR",
                reads=star_reads,
                read_ids=selected_star_ids,
                scope=scope,
            )
        )
        multi_by_nh_rows.extend(
            summarize_multi_by_nh(
                mapper_name="RiboTxAlign",
                reads=ribotx_reads,
                read_ids=selected_ribotx_ids,
                scope=scope,
            )
        )

    multi_by_nh = pd.DataFrame(
        multi_by_nh_rows
    )
    multi_by_nh_output = (
        OUTPUT_DIR
        / f"{sample_name}.multi_reads_by_nh.tsv"
    )
    multi_by_nh.to_csv(
        multi_by_nh_output,
        sep="\t",
        index=False,
        float_format="%.4f",
    )

    rescue_summary = pd.DataFrame(
        [
            summarize_ribotx_rescue(
                ribotx_reads=ribotx_reads,
                read_ids=ribotx_ids,
            )
        ]
    )
    rescue_summary_output = (
        OUTPUT_DIR
        / f"{sample_name}.ribotx_rescue_summary.tsv"
    )
    rescue_summary.to_csv(
        rescue_summary_output,
        sep="\t",
        index=False,
        float_format="%.4f",
    )

    # ========================================================
    # STAR mapped -> RiboTx unmapped
    # ========================================================

    star_mapped_ribotx_unmapped = set(
        df.loc[
            (
                df["star_status"]
                != "unmapped"
            )
            &
            (
                df["ribotx_status"]
                == "unmapped"
            ),
            "read_id",
        ]
    )

    # ========================================================
    # STAR mapped -> read absent from RiboTx BAM
    # ========================================================

    star_mapped_ribotx_absent = {
        read_id
        for read_id in star_only_ids
        if (
            star_reads[
                read_id
            ]["status"]
            != "unmapped"
        )
    }

    # ========================================================
    # Print examples
    # ========================================================

    print_examples(
        title=(
            "STAR mapped -> "
            "RiboTxAlign unmapped"
        ),
        read_ids=(
            star_mapped_ribotx_unmapped
        ),
        star_reads=star_reads,
    )

    print_examples(
        title=(
            "STAR mapped -> "
            "read absent from RiboTxAlign BAM"
        ),
        read_ids=(
            star_mapped_ribotx_absent
        ),
        star_reads=star_reads,
    )

    # ========================================================
    # Save detailed reads
    # ========================================================

    detail_rows = []

    detail_rows.extend(
        make_star_detail_rows(
            read_ids=(
                star_mapped_ribotx_unmapped
            ),
            star_reads=star_reads,
            category=(
                "ribotx_unmapped"
            ),
        )
    )

    detail_rows.extend(
        make_star_detail_rows(
            read_ids=(
                star_mapped_ribotx_absent
            ),
            star_reads=star_reads,
            category=(
                "ribotx_absent"
            ),
        )
    )

    detail_df = pd.DataFrame(
        detail_rows
    )

    detail_output = (
        OUTPUT_DIR
        / (
            f"{sample_name}."
            "star_mapped_ribotx_failed.tsv.gz"
        )
    )

    detail_df.to_csv(
        detail_output,
        sep="\t",
        index=False,
        compression="gzip",
    )

    # ========================================================
    # Transition table
    # ========================================================

    transition = pd.crosstab(
        df["star_status"],
        df["ribotx_status"],
    )
    transition = transition.reindex(
        index=[
            status
            for status in [
                "unmapped",
                "unique",
                "multi",
            ]
            if status in transition.index
        ],
        columns=[
            status
            for status in [
                "unmapped",
                "unique",
                "multi_assigned",
            ]
            if status in transition.columns
        ],
    )

    transition_output = (
        OUTPUT_DIR
        / (
            f"{sample_name}."
            "transition.tsv"
        )
    )

    transition.to_csv(
        transition_output,
        sep="\t",
    )

    transition_percent = (
        transition.div(
            transition.sum(axis=1),
            axis=0,
        )
        * 100.0
    )
    transition_percent_output = (
        OUTPUT_DIR
        / f"{sample_name}.transition_percent.tsv"
    )
    transition_percent.to_csv(
        transition_percent_output,
        sep="\t",
        float_format="%.4f",
    )

    # Detailed semantic status + NH transition.
    detailed_transition = pd.crosstab(
        df["star_nh_category"],
        df["ribotx_nh_category"],
    )
    detailed_transition = sort_mapping_transition(
        detailed_transition
    )
    detailed_transition_output = (
        OUTPUT_DIR
        / f"{sample_name}.detailed_transition.tsv"
    )
    detailed_transition.to_csv(
        detailed_transition_output,
        sep="\t",
    )

    detailed_transition_percent = (
        detailed_transition.div(
            detailed_transition.sum(axis=1),
            axis=0,
        )
        * 100.0
    )
    detailed_transition_percent_output = (
        OUTPUT_DIR
        / f"{sample_name}.detailed_transition_percent.tsv"
    )
    detailed_transition_percent.to_csv(
        detailed_transition_percent_output,
        sep="\t",
        float_format="%.4f",
    )

    # Numeric NH-to-NH matrix retained for exact multiplicity analysis.
    nh_transition = pd.crosstab(
        df["star_nh"],
        df["ribotx_nh"],
    )
    nh_transition = (
        nh_transition
        .sort_index(axis=0)
        .sort_index(axis=1)
    )
    nh_transition_output = (
        OUTPUT_DIR
        / f"{sample_name}.star_nh_to_ribotx_nh.tsv"
    )
    nh_transition.to_csv(
        nh_transition_output,
        sep="\t",
    )

    nh_transition_percent = (
        nh_transition.div(
            nh_transition.sum(axis=1),
            axis=0,
        )
        * 100.0
    )
    nh_transition_percent_output = (
        OUTPUT_DIR
        / f"{sample_name}.star_nh_to_ribotx_nh_percent.tsv"
    )
    nh_transition_percent.to_csv(
        nh_transition_percent_output,
        sep="\t",
        float_format="%.4f",
    )

    # RiboTxAlign assignment types preserve internal rescue provenance.
    assignment_transition = pd.crosstab(
        df["star_nh_category"],
        df["ribotx_assignment_type"],
    )
    assignment_transition = assignment_transition.reindex(
        index=sorted(
            assignment_transition.index,
            key=mapping_category_sort_key,
        ),
        columns=sorted(
            assignment_transition.columns,
        ),
    )
    assignment_transition_output = (
        OUTPUT_DIR
        / f"{sample_name}.star_to_ribotx_assignment.tsv"
    )
    assignment_transition.to_csv(
        assignment_transition_output,
        sep="\t",
    )

    assignment_transition_percent = (
        assignment_transition.div(
            assignment_transition.sum(axis=1),
            axis=0,
        )
        * 100.0
    )
    assignment_transition_percent_output = (
        OUTPUT_DIR
        / f"{sample_name}.star_to_ribotx_assignment_percent.tsv"
    )
    assignment_transition_percent.to_csv(
        assignment_transition_percent_output,
        sep="\t",
        float_format="%.4f",
    )

    # Among STAR multimappers, show the RiboTxAlign NH distribution.
    star_multi = df.loc[
        df["star_status"] == "multi"
    ].copy()

    star_multi_to_ribotx_nh = (
        star_multi.groupby(
            [
                "ribotx_status",
                "ribotx_nh",
            ],
            dropna=False,
        )
        .size()
        .rename("read_count")
        .reset_index()
    )
    star_multi_to_ribotx_nh[
        "percent_of_star_multi"
    ] = star_multi_to_ribotx_nh[
        "read_count"
    ].map(
        lambda value: safe_percent(
            value,
            len(star_multi),
        )
    )
    star_multi_to_ribotx_nh = (
        star_multi_to_ribotx_nh.sort_values(
            [
                "ribotx_nh",
                "ribotx_status",
            ],
        )
        .reset_index(drop=True)
    )
    star_multi_to_ribotx_nh_output = (
        OUTPUT_DIR
        / f"{sample_name}.star_multi_to_ribotx_nh.tsv"
    )
    star_multi_to_ribotx_nh.to_csv(
        star_multi_to_ribotx_nh_output,
        sep="\t",
        index=False,
        float_format="%.4f",
    )

    print()
    print(
        "========================================"
    )
    print(
        "STAR -> RiboTxAlign"
    )
    print(
        "========================================"
    )

    print(
        transition
    )

    print()
    print(
        "STAR -> RiboTxAlign (%)"
    )
    print(
        transition_percent.round(2)
    )

    print()
    print(
        "Detailed STAR NH -> RiboTxAlign NH"
    )
    print(
        detailed_transition
    )

    print()
    print(
        "Mapping and splice summary"
    )
    print(
        mapping_summary.loc[
            mapping_summary["scope"]
            == "common_reads"
        ].to_string(index=False)
    )

    print()
    print(
        "RiboTxAlign internal rescue summary"
    )
    print(
        rescue_summary.to_string(index=False)
    )

    # ========================================================
    # Summary
    # ========================================================

    print()
    print(
        "========================================"
    )
    print(
        "Failure summary"
    )
    print(
        "========================================"
    )

    print(
        "STAR mapped -> "
        "RiboTx unmapped : "
        f"{len(star_mapped_ribotx_unmapped):,}"
    )

    print(
        "STAR mapped -> "
        "RiboTx absent   : "
        f"{len(star_mapped_ribotx_absent):,}"
    )

    print()
    print(
        f"Detailed output:"
    )

    print(
        detail_output
    )

    print(
        f"Mapping summary: {mapping_summary_output}"
    )
    print(
        f"Multi by NH: {multi_by_nh_output}"
    )
    print(
        f"Rescue summary: {rescue_summary_output}"
    )


# ============================================================
# Main
# ============================================================

bam_pairs = (
    find_bam_pairs()
)

print()
print(
    "============================================================"
)
print(
    f"Matched samples: "
    f"{len(bam_pairs)}"
)
print(
    "============================================================"
)

for (
    sample_name,
    star_bam,
    ribotx_bam,
) in bam_pairs:

    print(
        sample_name
    )

    print(
        f"  STAR   : "
        f"{star_bam.name}"
    )

    print(
        f"  RiboTx : "
        f"{ribotx_bam}"
    )


for (
    sample_name,
    star_bam,
    ribotx_bam,
) in bam_pairs:

    compare_sample(
        sample_name=sample_name,
        star_bam=star_bam,
        ribotx_bam=ribotx_bam,
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