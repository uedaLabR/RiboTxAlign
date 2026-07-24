from __future__ import annotations

import csv
from collections import Counter

from mapping.mapping_result import ProjectedAlignment


def load_gene_tpm(gene_tpm_tsv: str) -> dict[str, float]:
    gene_tpm: dict[str, float] = {}
    with open(gene_tpm_tsv) as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        required = {'gene_id', 'gene_tpm'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError('Missing gene TPM columns: ' + ', '.join(sorted(missing)))
        for row in reader:
            try:
                value = float(row['gene_tpm'])
            except (TypeError, ValueError):
                value = 0.0
            gene_tpm[row['gene_id']] = max(0.0, value)
    return gene_tpm


def assign_weights_by_gene_tpm(
    alignments: list[ProjectedAlignment],
    gene_tpm: dict[str, float],
    pseudocount: float = 0.01,
) -> None:
    if not alignments:
        return

    counts = Counter(alignment.gene_id for alignment in alignments)
    values = {gene_id: gene_tpm.get(gene_id, 0.0) + pseudocount for gene_id in counts}
    total = sum(values.values())

    if total <= 0:
        equal_weight = 1.0 / len(alignments)
        for alignment in alignments:
            alignment.assignment_weight = equal_weight
        return

    for alignment in alignments:
        gene_weight = values[alignment.gene_id] / total
        alignment.assignment_weight = gene_weight / counts[alignment.gene_id]
