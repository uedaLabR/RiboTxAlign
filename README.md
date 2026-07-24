# RiboTxAlign

RiboTxAlign is a Python command-line pipeline for building a sample-specific transcript reference, mapping short Ribo-seq reads, and redistributing multi-mapping reads by expectation-maximization (EM).

The workflow has three main stages:

1. Build a merged transcript reference from GENCODE and, optionally, StringTie transcripts.
2. Map Ribo-seq reads to the genome and transcript reference, retaining transcript-derived genomic candidates with two to six hits.
3. Estimate transcript-group abundance and redistribute multi-mapping reads into a final coordinate-sorted BAM file.

## Requirements

- Linux
- Python 3.10 or later
- A genome FASTA
- A GENCODE GFF3 annotation
- Optional StringTie GTF output
- Single-end Ribo-seq FASTQ or FASTQ.GZ

Python dependencies are listed in `requirements.txt` and `pyproject.toml`:

- Click
- mappy
- NumPy
- pandas
- PyArrow
- pysam

## Installation

### Install from a cloned Git repository

```bash
git clone https://github.com/uedaLabR/RiboTxAlign.git
cd RiboTxAlign

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install .
```

The `pip install .` command installs RiboTxAlign and all Python dependencies declared in `pyproject.toml`.

Confirm the installation:

```bash
ribotxalign --help
ribotxalign --version
```

### Editable installation for development

Use an editable installation when modifying the source code:

```bash
git clone https://github.com/uedaLabR/RiboTxAlign.git
cd RiboTxAlign

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -e .
```

Changes to the Python source files will then be reflected without reinstalling the package.

### Install dependencies only

To run the source file directly without installing the command:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The CLI can then be started with:

```bash
python RiboTxAlign_Cmd.py --help
```

## Command overview

```text
ribotxalign create-reference   Build a transcript reference and mappability table.
ribotxalign map                Map reads and create unique/multi-mapping outputs.
ribotxalign redistribute       Run TPM estimation and EM redistribution.
ribotxalign run                Run mapping and redistribution consecutively.
```

Detailed help is available for every command:

```bash
ribotxalign create-reference --help
ribotxalign map --help
ribotxalign redistribute --help
ribotxalign run --help
```

## 1. Build the transcript reference

The `create-reference` command corresponds to the previous `TestCreateTx.py` workflow. It performs both reference construction and exact-match mappability calculation.

### GENCODE plus StringTie

```bash
ribotxalign create-reference \
  --gencode-gff3 /path/to/gencode.annotation.gff3 \
  --stringtie-gtf /path/to/sample.transcripts.gtf \
  --genome-fasta /path/to/hg38.fa \
  --output-prefix /path/to/output/HEK293T_RiboTxAlign \
  --flank-length 40 \
  --mappability-read-length 30 \
  --mappability-step 1
```

### GENCODE only

Omit `--stringtie-gtf`:

```bash
ribotxalign create-reference \
  --gencode-gff3 /path/to/gencode.annotation.gff3 \
  --genome-fasta /path/to/hg38.fa \
  --output-prefix /path/to/output/GENCODE_RiboTxAlign
```

### Reference outputs

For an output prefix named `HEK293T_RiboTxAlign`, the command creates:

```text
HEK293T_RiboTxAlign.fa
HEK293T_RiboTxAlign.coordinates.tsv
HEK293T_RiboTxAlign.transcript_paths.tsv
HEK293T_RiboTxAlign.mappability.parquet
```

The coordinate TSV is required for mapping and redistribution. The mappability Parquet file is required for initial TPM estimation and EM redistribution.

## 2. Map Ribo-seq reads

The `map` command corresponds to `TestMapping.py`.

```bash
ribotxalign map \
  --fastq /path/to/reads.fastq.gz \
  --genome-fasta /path/to/hg38.fa \
  --transcript-fasta /path/to/HEK293T_RiboTxAlign.fa \
  --coordinate-tsv /path/to/HEK293T_RiboTxAlign.coordinates.tsv \
  --output-prefix /path/to/output/SRR26942051 \
  --transcript-candidate-count 8 \
  --kmer-size 11 \
  --window-size 5 \
  --minimum-minimizers 2 \
  --minimum-chain-score 15 \
  --minimum-dp-score 15 \
  --minimum-query-coverage 0.85 \
  --minimum-identity 0.90 \
  --worker-count 6 \
  --reads-per-worker 50000 \
  --sort-threads 8 \
  --progress-interval 100000 \
  --parquet-batch-size 100000
```

### Mapping outputs

For an output prefix named `SRR26942051`, the command creates:

```text
SRR26942051.unique.bam
SRR26942051.multiplemap.bam
SRR26942051.multiplemap.parquet
SRR26942051.mapping_mt.log
```

The mapping log reports the cumulative number of reads in each class:

```text
GenomeUnique
TranscriptSingle
TranscriptMulti
MultiHit2
MultiHit3
MultiHit4
MultiHit5
MultiHit6
MultiHit7
MultiHit8
MultiHitOver8
Unmapped
```

Reads with two to eight projected genomic candidates are retained for EM redistribution. Reads with more than eight candidates are recorded as `MultiHitOver8`, preserved as unresolved reads in the multi-mapping BAM, and excluded from the EM candidate Parquet file.

When upgrading from a six-hit reference, rerun `create-reference` so that the mappability Parquet includes `hit7_mappable_length` and `hit8_mappable_length`.

### Important mapping options

`--maximum-reads` limits the number of reads and is useful for testing.

```bash
ribotxalign map ... --maximum-reads 100000
```

`--worker-count` controls the number of mapping processes. `--reads-per-worker` controls the number of reads passed to each worker in one batch.

The default multiprocessing start method is `spawn`, which is the safest general setting. On Linux, `fork` may start faster but should be tested carefully:

```bash
ribotxalign map ... --multiprocessing-start-method fork
```

## 3. Redistribute multi-mapping reads

The `redistribute` command corresponds to `TestRedistribute.py`.

```bash
ribotxalign redistribute \
  --coordinate-tsv /path/to/HEK293T_RiboTxAlign.coordinates.tsv \
  --unique-bam /path/to/output/SRR26942051.unique.bam \
  --multiplemap-bam /path/to/output/SRR26942051.multiplemap.bam \
  --mapping-parquet /path/to/output/SRR26942051.multiplemap.parquet \
  --mappability-parquet /path/to/HEK293T_RiboTxAlign.mappability.parquet \
  --output-prefix /path/to/output/SRR26942051.redistributed \
  --pseudocount 0.01 \
  --maximum-iterations 100 \
  --convergence-threshold 1e-7 \
  --maximum-allocation-error 0.10 \
  --random-seed 13 \
  --sort-threads 8
```

### Redistribution outputs

For an output prefix named `SRR26942051.redistributed`, the command creates files including:

```text
SRR26942051.redistributed.gene_groups.parquet
SRR26942051.redistributed.initial_tpm.parquet
SRR26942051.redistributed.final_tpm.parquet
SRR26942051.redistributed.em_history.tsv
SRR26942051.redistributed.em_assignments.parquet
SRR26942051.redistributed.redistributed_multi.bam
SRR26942051.redistributed.final.bam
SRR26942051.redistributed.redistribute.log
```

The final BAM contains unique reads and resolved multi-mapping reads. Sequence groups whose integer allocation differs from the EM expectation by more than `--maximum-allocation-error` remain unresolved.

## 4. Run mapping and redistribution together

The `run` command combines the previous `TestMapping.py` and `TestRedistribute.py` workflows. Mapping is completed first, and its output files are passed directly to redistribution.

```bash
ribotxalign run \
  --fastq /path/to/reads.fastq.gz \
  --genome-fasta /path/to/hg38.fa \
  --transcript-fasta /path/to/HEK293T_RiboTxAlign.fa \
  --coordinate-tsv /path/to/HEK293T_RiboTxAlign.coordinates.tsv \
  --mappability-parquet /path/to/HEK293T_RiboTxAlign.mappability.parquet \
  --output-prefix /path/to/output/SRR26942051 \
  --redistribution-output-prefix /path/to/output/SRR26942051.redistributed \
  --transcript-candidate-count 8 \
  --worker-count 6 \
  --reads-per-worker 50000 \
  --sort-threads 8 \
  --maximum-iterations 100 \
  --maximum-allocation-error 0.10
```

When `--redistribution-output-prefix` is omitted, the default is:

```text
<output-prefix>.redistributed
```

For example, `--output-prefix /data/SRR26942051` produces redistributed outputs beginning with:

```text
/data/SRR26942051.redistributed
```

## Running the command without installation

All commands can also be run through the source file:

```bash
python RiboTxAlign_Cmd.py create-reference --help
python RiboTxAlign_Cmd.py map --help
python RiboTxAlign_Cmd.py redistribute --help
python RiboTxAlign_Cmd.py run --help
```

## Example end-to-end workflow

```bash
# 1. Build the sample-specific reference.
ribotxalign create-reference \
  --gencode-gff3 /refs/gencode.v44.annotation.gff3 \
  --stringtie-gtf /data/HEK293T.transcripts.gtf \
  --genome-fasta /refs/hg38.fa \
  --output-prefix /data/reference/HEK293T_RiboTxAlign

# 2. Map and redistribute one Ribo-seq library.
ribotxalign run \
  --fastq /data/reads/SRR26942051.26_34.fastq.gz \
  --genome-fasta /refs/hg38.fa \
  --transcript-fasta /data/reference/HEK293T_RiboTxAlign.fa \
  --coordinate-tsv /data/reference/HEK293T_RiboTxAlign.coordinates.tsv \
  --mappability-parquet /data/reference/HEK293T_RiboTxAlign.mappability.parquet \
  --output-prefix /data/results/SRR26942051 \
  --worker-count 6 \
  --reads-per-worker 50000 \
  --sort-threads 8
```

## Notes for Git repositories

Do not commit generated BAM, BAI, Parquet, FASTA index, or large log files. A typical `.gitignore` may include:

```text
.venv/
__pycache__/
*.py[cod]
*.bam
*.bai
*.parquet
*.log
*.mmi
build/
dist/
*.egg-info/
```

## License

Add the intended license before public release. Common choices for scientific software include MIT, BSD-3-Clause, and GPL-3.0.
