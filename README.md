# RiboTxAlign

RiboTxAlign maps short, single-end ribosome profiling (Ribo-seq) reads to a sample-specific spliced transcript reference built from StringTie assemblies and GENCODE annotation. Two sequential expectation-maximization (EM) stages estimate support between transcript groups and within each group. The workflow produces fractional abundance estimates, posterior-weighted profiles, and a genome-coordinate BAM, including eligible reads recovered by genome rescue.

RNA-seq assembly is performed upstream. Either short-read- or long-read-derived StringTie models can be supplied. Estimated transcript support is not, by itself, proof that an individual isoform is actively translated.

## Entry points

| File | Role |
|---|---|
| `RiboTxAlign_Cmd.py` | User-facing command-line interface: `create-reference`, `deduplicate`, and `run` |
| `RiboTxAlign_pipeline.py` | Internal workflow implementation: paths, configuration, stages, child processes, validation, and cleanup |

Start analyses with **`python RiboTxAlign_Cmd.py`**. The pipeline module is imported by the CLI; it has no standalone command-line entry point. Keep both files and the repository's package directories together.

## Installation

Use Linux with Python 3.10 or later. From a source checkout:

```bash
git clone https://github.com/uedaLabR/RiboTxAlign.git
cd RiboTxAlign
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
# For BigWig tracks and profile plots:
python -m pip install pyBigWig matplotlib
python RiboTxAlign_Cmd.py --help
```

The examples below run from the repository root. The current packaging metadata does not include all workflow modules and dependencies, so use the source-checkout commands rather than relying on a standalone wheel installation. The CLI and package metadata also currently report different version strings; record the source commit for reproducible analyses.

## Inputs and configuration

Required inputs are:

- Genome FASTA and its `.fai` index, with chromosome names matching the annotation.
- GENCODE GFF3 annotation.
- A StringTie transcript GTF containing transcript structures and RNA abundance information.
- Single-end Ribo-seq FASTQ or FASTQ.GZ, after the desired adapter, length, and contaminant filtering.

Copy `RiboTxAlign_human.yaml` or `RiboTxAlign_mouse.yaml` and replace the reference paths. These files contain example machine-specific paths, not portable reference downloads.

```yaml
reference:
  gencode_gff3: /path/to/gencode.annotation.gff3
  genome_fasta: /path/to/genome.fa
  genome_fai: /path/to/genome.fa.fai

pipeline:
  delete_intermediate_files: false
  start_method: spawn
```

`genome_fai` defaults to `<genome_fasta>.fai`. YAML paths are resolved relative to the YAML file; command-line paths are resolved relative to the current directory. Sample name, output directory, FASTQ, and StringTie GTF are supplied on the command line.

**`run` always reads the file supplied by `--fastq`.** The legacy `pipeline.use_deduplicated_fastq` setting present in the example YAML files is not used by the current pipeline. Deduplication must be run separately, and its output passed explicitly to `run`.

## Quick start

### 1. Build the sample reference

```bash
python RiboTxAlign_Cmd.py create-reference \
  --config RiboTxAlign_human.yaml \
  --sample-name sample1 \
  --output-dir results \
  --stringtie-gtf /path/to/stringtie.transcripts.gtf
```

This creates `results/sample1/reference/sample1.{fa,gtf,coordinates.tsv,metadata.tsv}`. Use the same sample name and output directory in the mapping command.

### 2. Optional mapping-free UMI deduplication

```bash
python RiboTxAlign_Cmd.py deduplicate \
  --fastq /path/to/sample1.clean.fastq.gz \
  --output-fastq results/sample1/reads/sample1.deduplicated.fastq.gz
```

This command accepts only `--fastq` and `--output-fastq`; it does not take a YAML file or sample name. Statistics are written to `<output-fastq>.stats.tsv`, for example `sample1.deduplicated.fastq.gz.stats.tsv`. Temporary UMI buckets are removed after successful processing.

UMIs must already be extracted into the read identifiers, for example `read0001_ACGTACGTA`. The terminal run of uppercase A/C/G/T bases is interpreted as the UMI. UMI lengths must be consistent and at least five bases. The current implementation supports sequences up to 40 nt and uses encoded sequence blocks to compare equal-length sequences within a UMI, allowing up to three substitutions. This is sequence/UMI-based preprocessing before alignment, not genomic-coordinate deduplication. Mapping counts the retained representative records without expansion by family size.

### 3. Map and allocate reads

If deduplication was performed:

```bash
python RiboTxAlign_Cmd.py run \
  --config RiboTxAlign_human.yaml \
  --sample-name sample1 \
  --output-dir results \
  --fastq results/sample1/reads/sample1.deduplicated.fastq.gz
```

Otherwise pass the cleaned original FASTQ to `--fastq`. `run` neither builds the reference nor performs deduplication automatically.

## Workflow and parallel execution

1. Construct a sample-specific reference, supplementing compatible transcript structures and CDS annotation with GENCODE. Consolidate eligible CDS-equivalent models and assign groups connected by same-strand exon overlap.
2. Align footprints to the spliced reference. Among highest-score ties, retain CDS-overlapping candidates preferentially. A transcript's own CDS takes precedence over the group-derived fallback. Unique non-CDS hits and all-non-CDS candidate sets remain eligible.
3. Run group EM (EM1), followed by transcript EM (EM2) with fixed group weights. Before EM2, candidates without candidate-specific reads are merged only when their non-empty ambiguous read sets are identical. Lack of specific reads alone does not establish absence of translation.
4. Generate fractional estimates and weighted depth/P-site profiles. Project transcript alignments to genomic coordinates and allocate discrete placements for BAM output.
5. Rescue eligible transcript-unmapped reads by contiguous genome mapping, then conservative splice-aware mapping. Combine accepted rescued and transcript-supported placements, sort, index, and validate the BAM.

Stage labels `c` through `l` identify processing steps, **not a strictly serial execution order**:

```text
c: transcript mapping
├─ k: genome rescue (child process)
└─ d: sort mapping records → e: prepare group input → f: EM1 → g: EM2
                                                            ├─ h: depth/profiles (child process)
                                                            └─ i: genomic allocation → j: BAM writing
j + k → l: merge / sort / index / validate BAM
h + l → cleanup and successful completion
```

Rescue overlaps downstream EM processing, and depth generation overlaps genomic allocation/BAM writing. Account for concurrent branches when choosing worker counts and memory limits. A child-process failure fails the workflow and terminates the other active child branch.

Rescued reads contribute to the final BAM, not transcript TPM or posterior-weighted transcript profiles. Reads excluded as decoy-derived are not eligible for rescue.

## Outputs

```text
results/sample1/
├── reference/
│   ├── sample1.fa
│   ├── sample1.gtf
│   ├── sample1.coordinates.tsv
│   └── sample1.metadata.tsv
├── reads/                         # Optional deduplication output location chosen above
│   ├── sample1.deduplicated.fastq.gz
│   └── sample1.deduplicated.fastq.gz.stats.tsv
├── mapping/
│   ├── sample1.mapping.log
│   ├── sample1.mapping_stats.json
│   ├── sample1.unmapped.fastq.gz
│   └── sample1.genome_rescue_stats.json
├── em/
│   ├── sample1.group_tpm.tsv
│   ├── sample1.group_em.log
│   ├── sample1.transcript_tpm.tsv
│   ├── sample1.transcript_em.log
│   ├── sample1.transcript_merge.tsv
│   └── sample1.reads_info.parquet
├── depth/                         # Coverage, P-site profiles and configured tracks
├── bam/
│   ├── sample1.final_assignments.stats.tsv
│   ├── sample1.genome.sorted.bam
│   └── sample1.genome.sorted.bam.bai
└── logs/child_processes/           # Child-process error reports, when a stage fails
```

The transcript-unmapped FASTQ is the input to rescue, not a FASTQ of reads remaining unmapped after rescue. Final mapping and rescue counts are recorded in the JSON summaries.

`reads_info.parquet` retains read-level probabilities and original alignment coordinates. `transcript_id` identifies the EM representative; `alignment_transcript_id` identifies the reference used for the original alignment and genomic projection. Merged source IDs should not be interpreted as independently quantified isoforms.

The final BAM contains discrete placements. Its `XA` assignment tags distinguish transcript-supported and rescued classes, including `genome_rescued_unique` and `genome_splice_rescued_unique`. These assignment categories are not directly interchangeable with another mapper's definition of sequence uniqueness.

## Logs, runtime, and memory

Save both standard output and standard error to preserve workflow messages. On Linux, GNU `time` can additionally record elapsed time and maximum resident set size:

```bash
mkdir -p results/sample1/resource
/usr/bin/time -v -o results/sample1/resource/run.resource_usage.txt \
  python RiboTxAlign_Cmd.py run \
    --config RiboTxAlign_human.yaml \
    --sample-name sample1 \
    --output-dir results \
    --fastq results/sample1/reads/sample1.deduplicated.fastq.gz \
  > results/sample1/resource/run.log 2>&1
```

This measures `run`; reference construction and deduplication require separate measurements. The maximum RSS reported by GNU `time` is not a sampled sum of the simultaneous RSS of every worker. Retain the raw measurement and command when comparing resource use. Group and transcript EM logs provide their own processing details and timings.

## Rerunning and intermediate files

`run` restarts mapping from the beginning; it does not resume a checkpoint. With `mapping.overwrite: true`, old mapping shards are removed before remapping. `sorting.overwrite` and `group_input.overwrite` similarly control their dedicated intermediate outputs. Setting these to `false` makes existing outputs an error.

Use a fresh output directory for a new reference, changed settings, or benchmark run. Never run concurrent jobs with the same output directory and sample name.

With `pipeline.delete_intermediate_files: true`, temporary mapping, EM-input, assignment, and unsorted BAM files are removed after both output branches complete successfully. Keep it `false` when intermediate files are needed for diagnostics. Final read-level probabilities, abundance tables, statistics, and final BAM outputs are retained.

## Validation and paper data

Validation scripts are separate from the mapping pipeline. Their dataset paths may need editing for your environment.

- [BAM comparison](test/validation/README.md)
- [Mapping-count summaries](test/validation/README_counts.md)
- [Collecting paper results and resource logs](test/validation/README_collect_paper_results.md)
- [Simulation TPM plots](test/similation/README_plot_tpm.md)
- [Exon-level validation](test/similation/README_exon_validation.md)

The simulation directory is named `test/similation` in this repository.

## Command help

```bash
python RiboTxAlign_Cmd.py create-reference --help
python RiboTxAlign_Cmd.py deduplicate --help
python RiboTxAlign_Cmd.py run --help
```
