from pathlib import Path

from Bio import SeqIO


input_fasta = Path(
    "/mnt/share/project/riboML/decoy_fasta/hg38-tRNAs.fa"
)

output_fasta = Path(
    "/mnt/share/project/riboML/decoy_fasta/hg38-tRNAs.nr.fa"
)



records = []

for record in SeqIO.parse(
    input_fasta,
    "fasta",
):
    sequence = str(
        record.seq
    ).upper()

    if not sequence:
        continue

    record.seq = record.seq.upper()
    record.description = ""

    records.append(
        (
            sequence,
            record,
        )
    )


# Process longer sequences first.
records.sort(
    key=lambda item: len(item[0]),
    reverse=True,
)


kept_records = []
kept_sequences: list[str] = []

exact_duplicate_count = 0
subset_count = 0


for sequence, record in records:

    if sequence in kept_sequences:
        exact_duplicate_count += 1
        continue

    # Remove the sequence when it is fully contained
    # within an already retained longer sequence.
    if any(
        sequence in kept_sequence
        for kept_sequence in kept_sequences
    ):
        subset_count += 1
        continue

    kept_sequences.append(
        sequence
    )

    kept_records.append(
        record
    )


SeqIO.write(
    kept_records,
    output_fasta,
    "fasta",
)


print(
    "Input records:",
    len(records),
)

print(
    "Retained records:",
    len(kept_records),
)

print(
    "Removed exact duplicates:",
    exact_duplicate_count,
)

print(
    "Removed subset sequences:",
    subset_count,
)

print(
    "Output:",
    output_fasta,
)

