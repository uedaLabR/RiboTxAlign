"""Count source FASTQ records directly. Edit paths below; no CLI.

Does not infer mapping input from BAM counts or label absent records unmapped.
Requires pysam. Missing FASTQs are reported, never replaced by zero.
"""
from pathlib import Path
import csv
import logging
import time
import pysam

BASE = Path('/mnt/share/project/riboML/validation')
OUTPUT_DIR = BASE / 'fastq_read_counts'
SAMPLES = {
    'HEK293T': ['SRR26942051', 'SRR26942052', 'SRR26942053', 'SRR26942054'],
    'mouse_brain': ['SRR2163083', 'SRR2163093', 'SRR2163103'],
}
# These are explicit files to count, not proof of which FASTQ a historical run used.
# Add/correct actual mapper inputs from the corresponding run command if different.
FASTQS = []
for dataset, samples in SAMPLES.items():
    for sample in samples:
        FASTQS.append((dataset, sample, 'clean_before_dedup',
            BASE / dataset / 'riboseq/processed/clean' / f'{sample}.clean.fastq.gz'))
        if dataset == 'HEK293T':
            for workflow, folder in [('RiboTxAlign_longread', 'ribotx_validation'),
                                     ('RiboTxAlign_shortread', 'ribotx_validation_rnaseq')]:
                name = f'{dataset}_{sample}.clean'
                FASTQS.append((dataset, sample, workflow + '_deduplicated',
                    BASE / folder / name / 'reads' / f'{name}.deduplicated.fastq.gz'))


def count_fastq(path):
    """Stream parsed FASTQ records, checking sequence/quality presence and length."""
    total = 0
    with pysam.FastxFile(str(path)) as records:
        for record in records:
            if record.quality is None or len(record.sequence) != len(record.quality):
                raise ValueError(f'Invalid FASTQ sequence/quality: {record.name}')
            total += 1
            if total % 5_000_000 == 0:
                logging.info('%s: %s records', path.name, f'{total:,}')
    return total


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
    rows = []
    failed = False
    # Reuse only identical paths within this execution; no stale on-disk cache.
    counts = {}
    for dataset, sample, stage, path in FASTQS:
        row = dict(dataset=dataset, sample=sample, stage=stage, fastq=str(path),
                   fastq_records='', status='ok', message='')
        try:
            started = time.monotonic()
            if path not in counts:
                counts[path] = count_fastq(path)
            row['fastq_records'] = counts[path]
            logging.info('%s / %s / %s: %s reads (%.1fs)', dataset, sample,
                         stage, counts[path], time.monotonic() - started)
        except Exception as error:
            failed = True
            row.update(status='error', message=str(error))
            logging.error('%s: %s', path, error)
        rows.append(row)
    output = OUTPUT_DIR / 'fastq_read_counts.tsv'
    temporary = output.with_suffix('.tsv.tmp')
    with temporary.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)
    print(f'Output: {output}')
    if failed:
        raise SystemExit('Some FASTQs could not be counted. Check error rows and edit FASTQS paths.')


if __name__ == '__main__':
    main()
