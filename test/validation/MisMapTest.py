from pathlib import Path

import mappy as mp


# ============================================================
# Input
# ============================================================

FASTQ = Path(
    "/mnt/share/project/riboML/validation/bam_comparison/"
    "HEK293T_SRR26942051.clean.STAR_unique_RiboTx_unmapped.fastq"
)

REFERENCE = Path(
    "/mnt/share/project/riboML/validation/reference/hg38.fa"
)

import time
# ============================================================
# Mappy parameters
# ============================================================
start_time = time.perf_counter()
aligner = mp.Aligner(
    str(REFERENCE),
    preset="sr",
    k=13,
    w=3,
    min_cnt=1,
    min_chain_score=10,
    min_dp_score=10,
    best_n=16,
)

if not aligner:
    raise RuntimeError(
        f"Failed to load reference: {REFERENCE}"
    )


total_reads = 0
mapped_reads = 0
unmapped_reads = 0


for name, seq, qual in mp.fastx_read(
    str(FASTQ)
):

    total_reads += 1

    hits = list(
        aligner.map(
            seq,
            MD=True,
        )
    )

    if not hits:

        unmapped_reads += 1
        continue

    mapped_reads += 1

elapsed = time.perf_counter() - start_time

print()
print(
    "========================================"
)
print(
    "Mappy mapping summary"
)
print(
    "========================================"
)

print(
    f"Total reads   : {total_reads:,}"
)

print(
    f"Mapped reads  : {mapped_reads:,}"
)

print(
    f"Unmapped reads: {unmapped_reads:,}"
)

if total_reads > 0:

    print(
        f"Mapped rate   : "
        f"{mapped_reads / total_reads * 100:.2f}%"
    )

    print(
        f"Unmapped rate : "
        f"{unmapped_reads / total_reads * 100:.2f}%"
    )

    print(
        f"Speed          : "
        f"{total_reads / elapsed:,.0f} reads/sec"
    )