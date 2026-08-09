# RiboTxAlign

RiboTxAlign builds a sample-specific transcript reference, maps short single-end Ribo-seq reads, estimates transcript-group abundance by EM, rescues transcript-unmapped reads by contiguous and conservative splice-aware genome mapping, creates depth tracks, and writes a final genome-coordinate BAM.

The sample-specific values are supplied as command-line options. The YAML file contains stable reference paths and analysis parameters. Output names are generated centrally from `--sample-name`, so the names remain consistent throughout the workflow.

## Workflow

The original test scripts `a` to `l` are represented by three commands:

| Command | Steps | Purpose |
|---|---:|---|
| `create-reference` | a | Build the transcript reference from GENCODE, StringTie2, and the genome FASTA |
| `deduplicate` | b | Perform mapping-free UMI deduplication |
| `run` | c–l | Transcript mapping, Parquet ordering, EM, depth tracks, genome assignment, genome rescue, BAM merge, sorting, indexing, and validation |

During `run`, one progress line is printed before each stage:

```text
[c] Mapping reads to the transcript reference
[d] Sorting mapping Parquet files by group
[e] Preparing group-level EM input
[f] Estimating group-level TPM by EM
[g] Estimating transcript probabilities within groups
[h] Creating transcript and genome depth tracks
[i] Assigning one genomic alignment to each mapped read
[j] Writing the unsorted transcript-supported genome BAM
[k] Rescuing unique genome alignments from transcript-unmapped reads
[l] Merging, coordinate-sorting, indexing, and validating the final BAM
```

## Requirements

- Linux
- Python 3.10 or later
- Genome FASTA and FASTA index (`.fai`)
- GENCODE GFF3 annotation
- StringTie2 GTF
- Single-end Ribo-seq FASTQ or FASTQ.GZ

Python dependencies:

- Click
- PyYAML
- mappy
- NumPy
- pandas
- PyArrow
- pysam

## Installation

```bash
git clone https://github.com/uedaLabR/RiboTxAlign.git
cd RiboTxAlign

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Confirm the command:

```bash
ribotxalign --help
ribotxalign --version
```

## Configuration

The YAML file contains values that normally remain unchanged between samples:

```yaml
reference:
  gencode_gff3: /path/to/gencode.annotation.gff3
  genome_fasta: /path/to/genome.fa
  genome_fai: /path/to/genome.fa.fai

pipeline:
  use_deduplicated_fastq: true
  delete_intermediate_files: true
```

The following sample-specific values are command-line options and are not stored in the YAML file:

```text
--sample-name
--output-dir
--fastq
--stringtie-gtf
```

`reference.genome_fai` may be omitted when the index is located at `<genome_fasta>.fai`. All analysis parameters have defaults. The complete `RiboTxAlign.yaml` file contains the available settings used by the current pipeline. Relative paths in the YAML file are resolved relative to the YAML file itself. Command-line paths are resolved from the current working directory.

## Commands

Use the same `--sample-name` and `--output-dir` for all three commands belonging to one sample.

### 1. Build the reference

```bash
ribotxalign create-reference \
  --config RiboTxAlign.yaml \
  --sample-name SRR26942051 \
  --output-dir /path/to/RiboTxAlign_output \
  --stringtie-gtf /path/to/stringtie.transcripts.gtf
```

Standard outputs:

```text
<output_dir>/<sample_name>/reference/<sample_name>.fa
<output_dir>/<sample_name>/reference/<sample_name>.gtf
<output_dir>/<sample_name>/reference/<sample_name>.coordinates.tsv
<output_dir>/<sample_name>/reference/<sample_name>.metadata.tsv
```

### 2. Deduplicate reads

```bash
ribotxalign deduplicate \
  --config RiboTxAlign.yaml \
  --sample-name SRR26942051 \
  --output-dir /path/to/RiboTxAlign_output \
  --fastq /path/to/reads.fastq.gz
```

Standard outputs:

```text
<output_dir>/<sample_name>/reads/<sample_name>.deduplicated.fastq.gz
<output_dir>/<sample_name>/reads/<sample_name>.deduplication.stats.tsv
```

When UMI deduplication is not required, set:

```yaml
pipeline:
  use_deduplicated_fastq: false
```

The `run` command then uses the FASTQ supplied with `--fastq` directly.

### 3. Run steps c through l

```bash
ribotxalign run \
  --config RiboTxAlign.yaml \
  --sample-name SRR26942051 \
  --output-dir /path/to/RiboTxAlign_output \
  --fastq /path/to/reads.fastq.gz
```

The command uses the standard reference paths created by `create-reference`. The same `--sample-name` and `--output-dir` must therefore be used. When `pipeline.use_deduplicated_fastq` is true, `--fastq` identifies the original input and the command reads the standard deduplicated output for that sample.

## Standard output layout

With `pipeline.delete_intermediate_files: true`, temporary mapping shards, group-EM input files, FASTQ-order assignment parts, the temporary final assignment Parquet, and unsorted BAM files are removed only after genome rescue, BAM merge, coordinate sorting, indexing, and read-count validation have completed successfully.

```text
<output_dir>/<sample_name>/
├── reference/
│   ├── <sample_name>.fa
│   ├── <sample_name>.gtf
│   ├── <sample_name>.coordinates.tsv
│   └── <sample_name>.metadata.tsv
├── reads/
│   ├── <sample_name>.deduplicated.fastq.gz
│   ├── <sample_name>.unmapped.fastq.gz
│   └── <sample_name>.deduplication.stats.tsv
├── mapping/
│   ├── <sample_name>.mapping.log
│   └── <sample_name>.mapping_stats.json
├── em/
│   ├── <sample_name>.group_tpm.tsv
│   ├── <sample_name>.group_em.log
│   ├── <sample_name>.transcript_tpm.tsv
│   ├── <sample_name>.transcript_merge.tsv
│   └── <sample_name>.reads_info.parquet
├── depth/
│   └── <sample_name>.*
└── bam/
    ├── <sample_name>.final_assignments.stats.tsv
    ├── <sample_name>.genome.sorted.bam
    └── <sample_name>.genome.sorted.bam.bai
```

The retained `<sample_name>.reads_info.parquet` is a primary output rather
than a temporary file. It contains the final read-level transcript information:

```text
group_id
transcript_id
alignment_transcript_id
read_id
start
end
read_length
hit_index
query_start
query_end
mapping_strand
cigar
prior_probability
transcript_probability
other_transcript_ids
```

Set the following only when intermediate files are required for debugging or
resuming individual internal stages:

```yaml
pipeline:
  delete_intermediate_files: false
```



`transcript_id` is the representative transcript used by transcript-level EM.
`alignment_transcript_id` is the original transcript reference against which
the read coordinates and CIGAR were produced. Genome projection always uses
`alignment_transcript_id`, so transcript merging cannot invalidate the stored
alignment coordinates.

The original transcript-alignment CIGAR is retained through transcript EM.
During genome assignment, exon-junction `N` operations are inserted without
discarding `I`, `D`, `S`, `=`, or `X` operations. For minus-strand
transcripts, the genomic CIGAR operation order is reversed. The BAM writer
reverse-complements the sequence and reverses quality scores when the final
genomic alignment strand is reverse.

## Unmapped-read genome rescue

During transcript-reference mapping, reads with no accepted transcript alignment are written to:

```text
<output_dir>/<sample_name>/reads/<sample_name>.unmapped.fastq.gz
```

Reads removed as decoy-derived are not included. Step `k` loads two mappy
aligners, each with `n_threads=1`, and processes every read in two stages:

1. contiguous genome rescue using the short-read preset;
2. splice rescue only when contiguous rescue fails.

The splice rescue is deliberately conservative. An alignment is retained only
when it has one uniquely best score, exactly one intron, at least the configured
anchor length on both sides, an intron within the configured length range, and
an allowed canonical splice motif. The default settings are:

```text
minimum anchor length: 11 nt
minimum intron length: 20 nt
maximum intron length: 2000 nt
canonical motifs: GT-AG, GC-AG, AT-AC
```

The two rescued classes are written to the same unsorted rescue BAM and are
distinguished by the `XA` tag:

```text
XA:Z:genome_rescued_unique
XA:Z:genome_splice_rescued_unique
```

Reads that remain unresolved after both stages are counted in the mapping statistics but are not written to another FASTQ file.

The machine-readable mapping summary is updated throughout the workflow:

```text
<output_dir>/<sample_name>/mapping/<sample_name>.mapping_stats.json
```

It separately records contiguous and splice-rescued counts, their sum, final
mapped/unique/multi/unmapped counts, and final BAM validation results. Rescue
parameters are configured in the YAML `rescue` section.

## Multi-mapping assignment

Transcript candidates that produce the same genome alignment are collapsed before BAM output. Genome-unique alignments are fixed directly. Remaining multi-mapping reads are assigned at transcript-group level.

The assignment objective minimizes the change in relative squared error between assigned read counts and posterior-derived group expectations. It therefore preserves the group TPM pattern rather than minimizing absolute count deficits, which would otherwise favor highly expressed groups.

## BAM output

The final assignment table is ordered according to the FASTQ read order and divided into contiguous Parquet parts. The FASTQ is then streamed once while the corresponding Parquet part is exchanged at each boundary.

Step `j` writes an unsorted BAM containing transcript-supported genome alignments. Step `k` first performs contiguous genome rescue and then applies conservative splice rescue only to unresolved reads. Rescued records are marked with one of:

```text
XA:Z:genome_rescued_unique
XA:Z:genome_splice_rescued_unique
NH:i:1
```

Step `l` merges the transcript-supported and rescued BAMs, coordinate-sorts and indexes the merged BAM, and verifies that BAM read counts agree with `<sample_name>.mapping_stats.json`. The final BAM contains all retained alignments; the `XA` tag distinguishes original genome-unique, multi-map-assigned, contiguous-rescued, and splice-rescued reads.

## Command help

```bash
ribotxalign create-reference --help
ribotxalign deduplicate --help
ribotxalign run --help
```
