from pathlib import Path

from Bio import SeqIO


input_path = Path(
    "/mnt/share/project/riboML/decoy_fasta/"
    "human_nuclear_rRNA.fa"
)

output_path = Path(
    "/mnt/share/project/riboML/decoy_fasta/"
    "human_nuclear_rRNA_named.fa"
)

accession_to_name = {
    "NR_003286.4": "18S_rRNA",
    "NR_003287.4": "28S_rRNA",
    "NR_003285.3": "5_8S_rRNA",
    "NR_023363.1": "5S_rRNA",
}

records = []

for record in SeqIO.parse(
    input_path,
    "fasta",
):
    accession = record.id

    if accession not in accession_to_name:
        raise ValueError(
            f"Unexpected accession: {accession}"
        )

    reference_id = accession_to_name[
        accession
    ]

    record.id = reference_id
    record.name = reference_id
    record.description = ""

    records.append(
        record
    )

SeqIO.write(
    records,
    output_path,
    "fasta",
)