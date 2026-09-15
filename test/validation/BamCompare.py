#!/usr/bin/env python3

from collections import defaultdict
from pathlib import Path
import gzip

import pandas as pd
import pysam


SAMPLES = [
    "SRR26942051",
    "SRR26942052",
    "SRR26942053",
    "SRR26942054",
]

FASTQ_DIR = Path("/mnt/share/project/riboML/validation/HEK293T/riboseq/fastq")
STAR_GENOME_DIR = Path("/mnt/share/project/riboML/validation/star_validation")
STAR_TRANSCRIPT_DIR = Path("/mnt/share/project/riboML/validation/star_transcript_validation")
RIBOTX_DIR = Path("/mnt/share/project/riboML/validation/ribotx_validation_rnaseq")
OUTPUT_DIR = Path("/mnt/share/project/riboML/validation/three_mapper_comparison")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def find_one_bam(root: Path, sample: str, preferred_names: list[str]) -> Path:
    for name in preferred_names:
        path = root / name
        if path.exists():
            return path

    candidates = sorted(
        p for p in root.rglob(f"*{sample}*.bam")
        if not p.name.endswith(".bai")
    )

    if len(candidates) == 1:
        return candidates[0]

    if len(candidates) == 0:
        raise FileNotFoundError(
            f"No BAM found for {sample} under {root}"
        )

    raise RuntimeError(
        f"Multiple BAMs found for {sample} under {root}:\n"
        + "\n".join(str(p) for p in candidates)
    )


def get_star_genome_bam(sample: str) -> Path:
    return find_one_bam(
        STAR_GENOME_DIR,
        sample,
        [
            f"HEK293T_{sample}.deduplicated.bam",
            f"{sample}.deduplicated.bam",
            f"HEK293T_{sample}.Aligned.sortedByCoord.out.bam",
            f"{sample}.Aligned.sortedByCoord.out.bam",
        ],
    )


def get_star_transcript_bam(sample: str) -> Path:
    return find_one_bam(
        STAR_TRANSCRIPT_DIR,
        sample,
        [
            f"{sample}/{sample}.Aligned.out.bam",
            f"{sample}/{sample}.Aligned.sortedByCoord.out.bam",
            f"{sample}.Aligned.out.bam",
            f"{sample}.Aligned.sortedByCoord.out.bam",
        ],
    )


def get_ribotx_bam(sample: str) -> Path:
    return find_one_bam(
        RIBOTX_DIR,
        sample,
        [
            f"HEK293T_{sample}/bam/HEK293T_{sample}.genome.sorted.bam",
            f"{sample}/bam/{sample}.genome.sorted.bam",
        ],
    )


def count_fastq_reads(fastq_path: Path) -> int:
    opener = gzip.open if fastq_path.suffix == ".gz" else open
    line_count = 0

    with opener(fastq_path, "rt") as handle:
        for _ in handle:
            line_count += 1

    if line_count % 4 != 0:
        raise ValueError(
            f"FASTQ line count is not divisible by 4: {fastq_path}"
        )

    return line_count // 4


def summarize_star_bam(
    bam_path: Path,
    total_reads: int,
    detect_splice: bool,
) -> tuple[dict, list[dict]]:

    mapped_reads = 0
    unique_reads = 0
    multi_reads = 0
    splice_reads = 0
    multi_by_nh = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for record in bam.fetch(until_eof=True):

            if (
                record.is_unmapped
                or record.is_secondary
                or record.is_supplementary
            ):
                continue

            mapped_reads += 1

            nh = int(record.get_tag("NH")) if record.has_tag("NH") else 1

            if nh == 1:
                unique_reads += 1
            else:
                multi_reads += 1
                multi_by_nh[nh] += 1

            if (
                detect_splice
                and record.cigarstring
                and "N" in record.cigarstring
            ):
                splice_reads += 1

    unmapped_reads = max(total_reads - mapped_reads, 0)

    summary = {
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "unique_reads": unique_reads,
        "multi_reads": multi_reads,
        "unmapped_reads": unmapped_reads,
        "unique_percent_total": unique_reads / total_reads * 100 if total_reads else 0.0,
        "multi_percent_total": multi_reads / total_reads * 100 if total_reads else 0.0,
        "unmapped_percent_total": unmapped_reads / total_reads * 100 if total_reads else 0.0,
        "unique_percent_mapped": unique_reads / mapped_reads * 100 if mapped_reads else 0.0,
        "multi_percent_mapped": multi_reads / mapped_reads * 100 if mapped_reads else 0.0,
    }

    if detect_splice:
        summary["splice_reads"] = splice_reads
        summary["splice_percent_mapped"] = splice_reads / mapped_reads * 100 if mapped_reads else 0.0
        summary["splice_percent_total"] = splice_reads / total_reads * 100 if total_reads else 0.0
    else:
        summary["splice_reads"] = "N.D."
        summary["splice_percent_mapped"] = "N.D."
        summary["splice_percent_total"] = "N.D."

    nh_rows = [
        {
            "nh": nh,
            "read_count": count,
            "percent_of_multi": count / multi_reads * 100 if multi_reads else 0.0,
            "percent_of_total": count / total_reads * 100 if total_reads else 0.0,
        }
        for nh, count in sorted(multi_by_nh.items())
    ]

    return summary, nh_rows


def summarize_ribotx_bam(
    bam_path: Path,
    total_reads: int,
) -> tuple[dict, list[dict]]:

    mapped_reads = 0
    unique_reads = 0
    multi_reads = 0
    splice_reads = 0
    multi_by_nh = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for record in bam.fetch(until_eof=True):

            if (
                record.is_unmapped
                or record.is_secondary
                or record.is_supplementary
            ):
                continue

            mapped_reads += 1

            nh = int(record.get_tag("NH")) if record.has_tag("NH") else 1

            assignment_type = (
                str(record.get_tag("XA"))
                if record.has_tag("XA")
                else "unknown"
            )

            if assignment_type == "multimap_assigned":
                multi_reads += 1
                multi_by_nh[nh] += 1
            else:
                unique_reads += 1

            if record.cigarstring and "N" in record.cigarstring:
                splice_reads += 1

    unmapped_reads = max(total_reads - mapped_reads, 0)

    summary = {
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "unique_reads": unique_reads,
        "multi_reads": multi_reads,
        "unmapped_reads": unmapped_reads,
        "unique_percent_total": unique_reads / total_reads * 100 if total_reads else 0.0,
        "multi_percent_total": multi_reads / total_reads * 100 if total_reads else 0.0,
        "unmapped_percent_total": unmapped_reads / total_reads * 100 if total_reads else 0.0,
        "unique_percent_mapped": unique_reads / mapped_reads * 100 if mapped_reads else 0.0,
        "multi_percent_mapped": multi_reads / mapped_reads * 100 if mapped_reads else 0.0,
        "splice_reads": splice_reads,
        "splice_percent_mapped": splice_reads / mapped_reads * 100 if mapped_reads else 0.0,
        "splice_percent_total": splice_reads / total_reads * 100 if total_reads else 0.0,
    }

    nh_rows = [
        {
            "nh": nh,
            "read_count": count,
            "percent_of_multi": count / multi_reads * 100 if multi_reads else 0.0,
            "percent_of_total": count / total_reads * 100 if total_reads else 0.0,
        }
        for nh, count in sorted(multi_by_nh.items())
    ]

    return summary, nh_rows


summary_rows = []
multi_nh_rows = []

for sample in SAMPLES:

    print()
    print("=" * 60)
    print(f"Sample: {sample}")
    print("=" * 60)

    fastq = FASTQ_DIR / f"{sample}.fastq.gz"
    star_genome_bam = get_star_genome_bam(sample)
    star_transcript_bam = get_star_transcript_bam(sample)
    ribotx_bam = get_ribotx_bam(sample)

    print(f"FASTQ           : {fastq}")
    print(f"STAR genome     : {star_genome_bam}")
    print(f"STAR transcript : {star_transcript_bam}")
    print(f"RiboTxAlign     : {ribotx_bam}")

    total_reads = count_fastq_reads(fastq)
    print(f"Total reads: {total_reads:,}")

    for mapper, func, bam_path, detect_splice in [
        ("STAR-genome", summarize_star_bam, star_genome_bam, True),
        ("STAR-transcript", summarize_star_bam, star_transcript_bam, False),
    ]:
        summary, nh_rows = func(
            bam_path=bam_path,
            total_reads=total_reads,
            detect_splice=detect_splice,
        )

        summary_rows.append(
            {
                "sample": sample,
                "mapper": mapper,
                **summary,
            }
        )

        for row in nh_rows:
            multi_nh_rows.append(
                {
                    "sample": sample,
                    "mapper": mapper,
                    **row,
                }
            )

    summary, nh_rows = summarize_ribotx_bam(
        bam_path=ribotx_bam,
        total_reads=total_reads,
    )

    summary_rows.append(
        {
            "sample": sample,
            "mapper": "RiboTxAlign",
            **summary,
        }
    )

    for row in nh_rows:
        multi_nh_rows.append(
            {
                "sample": sample,
                "mapper": "RiboTxAlign",
                **row,
            }
        )


summary_df = pd.DataFrame(summary_rows)
multi_nh_df = pd.DataFrame(multi_nh_rows)

summary_output = OUTPUT_DIR / "HEK293T.mapping_summary.tsv"
multi_nh_output = OUTPUT_DIR / "HEK293T.multi_by_nh.tsv"

summary_df.to_csv(
    summary_output,
    sep="\t",
    index=False,
    float_format="%.4f",
)

multi_nh_df.to_csv(
    multi_nh_output,
    sep="\t",
    index=False,
    float_format="%.4f",
)

print()
print("=" * 60)
print("Final summary")
print("=" * 60)

print(
    summary_df[
        [
            "sample",
            "mapper",
            "total_reads",
            "mapped_reads",
            "unique_reads",
            "multi_reads",
            "unmapped_reads",
            "splice_reads",
        ]
    ].to_string(index=False)
)

print()
print(f"Summary : {summary_output}")
print(f"Multi NH: {multi_nh_output}")