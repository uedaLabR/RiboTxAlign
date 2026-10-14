#!/usr/bin/env python3
"""Collect completed validation outputs and run evidence; never rerun mapping.

Edit the configuration below and run with Python 3.10+. Standard library only.
Tables are UTF-8 TSV (Excel-readable); the ZIP is the input to manuscript,
figure and supplementary-workbook revision. No BAM/FASTQ is included in it.
"""
from __future__ import annotations

import ast
import csv
import gzip
import hashlib
import json
import math
import re
import shlex
import shutil
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

# ---------------- Edit server paths here; no command-line options needed. -----
PROJECT = Path(__file__).resolve().parents[2]
VALIDATION = Path('/mnt/share/project/riboML/validation')
SIMULATION = Path('/mnt/share/project/riboML/validation_expression/similation')
MAPPING_COUNTS = VALIDATION / 'bam_comparison_four_mappers_counts'
TPM_REPORT = SIMULATION / 'publication_tpm'
PLOT_DIRS = [SIMULATION / 'plots', SIMULATION / 'plots_simple']
OUTPUT_PARENT = VALIDATION / 'paper_collection'
RUN_ROOTS = {
    'RiboTxAlign_longread_stringtie': VALIDATION / 'ribotx_validation',
    'RiboTxAlign_shortread_stringtie': VALIDATION / 'ribotx_validation_rnaseq',
}
DATASETS = {
    'HEK293T': ['SRR26942051', 'SRR26942052', 'SRR26942053', 'SRR26942054'],
    'mouse_brain': ['SRR2163083', 'SRR2163093', 'SRR2163103'],
}
MAPPERS = ['STAR_genome', 'STAR_transcript', *RUN_ROOTS]
# Current supplied shells deduplicate HEK293T, but not mouse brain.
DEDUP_EXPECTED = {'HEK293T': True, 'mouse_brain': False}
# Only needed if a completed dedup stats file is missing. Counts stream through
# FASTQ once and are cached by absolute path, size and nanosecond mtime.
COUNT_FASTQ_IF_NEEDED = True
# Add outputs/log directories if the server's *_1.py files use different paths.
# These directories should contain only publication tables/figures, not reads.
EXTRA_RESULT_DIRS: list[Path] = []
EXTRA_RESOURCE_FILES: list[Path] = []
# Optional exact session choice per run directory: '.' or 'profile.<timestamp>.<id>'.
RESOURCE_SESSION_OVERRIDES: dict[str, str] = {}
MAX_ARTIFACT_BYTES = 512 * 1024 * 1024


def tsv_rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row)) or ['source']
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)


def file_signature(path):
    path = Path(path).resolve()
    stat = path.stat()
    return {'path': str(path), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def elapsed_seconds(value):
    parts = value.strip().split(':')
    if not 1 <= len(parts) <= 3:
        raise ValueError(f'Invalid GNU time duration: {value}')
    result = 0.0
    for part in parts:
        number = float(part)
        if not math.isfinite(number) or number < 0:
            raise ValueError(f'Invalid GNU time duration: {value}')
        result = result * 60 + number
    return result


def parse_resource(path):
    """Preserve each invocation; never mix fields from concatenated runs."""
    keys = {
        'Command being timed': 'command',
        'User time (seconds)': 'user_seconds',
        'System time (seconds)': 'system_seconds',
        'Percent of CPU this job got': 'cpu_percent',
        'Elapsed (wall clock) time (h:mm:ss or m:ss)': 'elapsed_raw',
        'Maximum resident set size (kbytes)': 'max_rss_kbytes',
        'Exit status': 'exit_status',
    }
    blocks, current = [], {}
    for line in Path(path).read_text(encoding='utf-8', errors='replace').splitlines():
        line = line.strip()
        for label, key in keys.items():
            if line.startswith(label + ':'):
                if key == 'command' and current:
                    blocks.append(current)
                    current = {}
                current[key] = line[len(label)+1:].strip()
                break
    if current:
        blocks.append(current)
    if not blocks:
        raise ValueError('No GNU /usr/bin/time -v record found')
    for index, values in enumerate(blocks, 1):
        values.update(invocation=index, invocations_in_file=len(blocks))
        values['complete'] = all(k in values for k in ('command', 'elapsed_raw', 'max_rss_kbytes', 'exit_status'))
        values['command'] = values.get('command', '').strip('"')
        if 'elapsed_raw' in values:
            values['elapsed_seconds'] = elapsed_seconds(values['elapsed_raw'])
            values['elapsed_minutes'] = values['elapsed_seconds'] / 60
        for key in ('max_rss_kbytes', 'exit_status'):
            if key in values:
                values[key] = int(values[key])
        if 'max_rss_kbytes' in values:
            if values['max_rss_kbytes'] < 0:
                raise ValueError('Negative maximum RSS')
            values['max_rss_gib'] = values['max_rss_kbytes'] / 1048576
        for key in ('user_seconds', 'system_seconds'):
            if key in values:
                values[key] = float(values[key])
        values['successful'] = values.get('exit_status') == 0
    return blocks


def command_option(command, option):
    try:
        words = shlex.split(command)
    except ValueError:
        return None
    for i, word in enumerate(words):
        if word == option and i + 1 < len(words):
            return words[i+1]
        if word.startswith(option + '='):
            return word.split('=', 1)[1]
    return None


def parse_dedup(path):
    values = {row['metric']: int(row['value']) for row in tsv_rows(path)}
    before, after = values['total_reads'], values['output_reads']
    eligible = values.get('written_records')
    if not 0 <= after <= before or (eligible is not None and not after <= eligible <= before):
        raise ValueError('Dedup counts must satisfy output <= written <= total')
    skip_keys = ['skipped_short', 'skipped_long', 'skipped_non_acgt']
    if eligible is not None and all(k in values for k in skip_keys):
        if before != eligible + sum(values[k] for k in skip_keys):
            raise ValueError('Dedup total != written_records + skipped counts')
    return {
        'input_reads': before, 'output_reads': after,
        'eligible_reads': eligible, 'filtered_reads': None if eligible is None else before-eligible,
        'duplicate_reads_removed': None if eligible is None else eligible-after,
        'total_reduction_reads': before-after,
        'total_reduction_percent': (before-after)/before*100 if before else None,
        'duplicate_percent_of_eligible': (eligible-after)/eligible*100 if eligible else None,
        **{k: values.get(k) for k in skip_keys},
    }


def count_fastq(path, cache):
    """Count records, including wrapped FASTQ; validate lengths and truncation."""
    path = Path(path)
    signature = file_signature(path)
    key = str(path.resolve())
    if key in cache and cache[key].get('signature') == signature:
        return cache[key]['records']
    print(f'Counting FASTQ (streaming): {path}', flush=True)
    opener = gzip.open if path.suffix == '.gz' else open
    count = 0
    with opener(path, 'rt', encoding='ascii') as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            if not header.startswith('@'):
                raise ValueError(f'Invalid FASTQ header at record {count+1}: {path}')
            nseq = 0
            while True:
                line = handle.readline()
                if not line:
                    raise ValueError(f'Truncated FASTQ sequence: {path}')
                if line.startswith('+'):
                    break
                nseq += len(line.rstrip('\r\n'))
            if not nseq:
                raise ValueError(f'Empty FASTQ sequence: {path}')
            nqual = 0
            while nqual < nseq:
                line = handle.readline()
                if not line:
                    raise ValueError(f'Truncated FASTQ quality: {path}')
                nqual += len(line.rstrip('\r\n'))
            if nseq != nqual:
                raise ValueError(f'FASTQ sequence/quality length mismatch: {path}')
            count += 1
    if signature != file_signature(path):
        raise ValueError(f'FASTQ changed while counting: {path}')
    cache[key] = {'signature': signature, 'records': count}
    return count


def static_paths(script):
    """Read simple path constants without importing/executing analysis scripts."""
    env = {}

    def value(node):
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Dict):
            return {value(k): value(v) for k,v in zip(node.keys, node.values)}
        if isinstance(node, (ast.List, ast.Tuple)):
            return [value(v) for v in node.elts]
        if isinstance(node, ast.Name):
            return env[node.id]
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'Path' and len(node.args) == 1:
            return Path(value(node.args[0]))
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
            return value(node.left) / value(node.right)
        raise ValueError('Dynamic expression')

    for node in ast.parse(Path(script).read_text(encoding='utf-8-sig')).body:
        if isinstance(node, ast.Assign):
            try:
                result = value(node.value)
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        env[target.id] = result
            except (KeyError, ValueError, TypeError):
                pass
    paths = {}
    def flatten(prefix, item):
        if isinstance(item, Path):
            paths[prefix] = item
        elif isinstance(item, dict):
            for k, v in item.items():
                flatten(f'{prefix}.{k}', v)
        elif isinstance(item, list):
            for k, v in enumerate(item):
                flatten(f'{prefix}.{k}', v)
    for k, v in env.items():
        flatten(k, v)
    return paths


class Collector:
    def __init__(self, output):
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=False)
        self.manifest, self.issues, self.tables = [], [], []
        self.copied = {}

    def issue(self, source, message, severity='warning'):
        self.issues.append({'severity': severity, 'source': str(source), 'message': message})

    def capture(self, path, relative, purpose):
        path = Path(path)
        if not path.is_file():
            self.issue(path, 'File missing', 'error')
            return None
        key = str(path.resolve())
        if key in self.copied:
            return self.copied[key]
        sig = file_signature(path)
        if sig['size'] > MAX_ARTIFACT_BYTES:
            self.issue(path, 'Artifact exceeds MAX_ARTIFACT_BYTES; not copied', 'error')
            return None
        target = self.output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        if file_signature(path) != sig:
            target.unlink()
            self.issue(path, 'Source changed during collection; repeat after completion', 'error')
            return None
        self.manifest.append({**sig, 'collected_path': str(target.relative_to(self.output)),
                              'sha256': sha256(target), 'purpose': purpose})
        self.copied[key] = target
        return target

    def table(self, name, rows, purpose, fields=None):
        rows = list(rows)
        write_tsv(self.output/name, rows, fields)
        self.tables.append({'file': name, 'rows': len(rows), 'purpose': purpose})


def collect_mapping(c, root, datasets, mappers):
    merged, nh, rescue = [], [], []
    for dataset, samples in datasets.items():
        for sample in samples:
            directory = root/dataset/sample
            source = directory/'mapping_summary.tsv'
            snapshot = c.capture(source, f'raw/mapping/{dataset}/{sample}/{source.name}', 'Fig2 real-library counts')
            if snapshot is None:
                continue
            try:
                rows = [r for r in tsv_rows(snapshot) if r.get('scope') == 'all_mapper_reads']
                if sorted(r['mapper'] for r in rows) != sorted(mappers):
                    raise ValueError('Missing, unexpected or duplicate mapper rows')
                for r in rows:
                    for key in ('total_reads', 'mapped_reads', 'unmapped_reads', 'unique_reads', 'multi_reads'):
                        r[key] = int(r[key])
                        if r[key] < 0:
                            raise ValueError('Negative mapping count')
                    if r['mapped_reads'] + r['unmapped_reads'] != r['total_reads'] or r['unique_reads'] + r['multi_reads'] != r['mapped_reads']:
                        raise ValueError('Mapping counts do not reconcile')
                denominator = max(r['total_reads'] for r in rows)
                for r in rows:
                    r.update(dataset=dataset, sample=sample, source=str(source),
                             common_denominator=denominator, denominator_type='maximum_observed_mapper_total',
                             eligible_for_summary=True)
                    for category in ('mapped', 'unique', 'multi'):
                        r[category+'_percent_common'] = 100*r[category+'_reads']/denominator if denominator else None
                    r['residual_reads_common'] = denominator-r['mapped_reads']
                    r['residual_percent_common'] = 100*r['residual_reads_common']/denominator if denominator else None
                    # Verify completed-count cache still refers to this BAM version.
                    cache_path = directory/(r['mapper']+'.counts.json')
                    r['bam_signature_matches'] = None
                    if cache_path.is_file():
                        cp = c.capture(cache_path, f'raw/mapping/{dataset}/{sample}/{cache_path.name}', 'BAM count provenance')
                        cache = json.loads(cp.read_text()) if cp else {}
                        sig = cache.get('signature', {})
                        bam = Path(sig.get('path', ''))
                        if bam.is_file():
                            r['bam_signature_matches'] = all(file_signature(bam)[k] == sig.get(k) for k in ('path','size','mtime_ns'))
                            if not r['bam_signature_matches']:
                                r['eligible_for_summary'] = False
                                c.issue(cache_path, 'BAM changed after counts were generated; rerun Bam4Comapre_1.py', 'error')
                        else:
                            c.issue(cache_path, 'BAM not accessible for cache freshness check')
                    else:
                        c.issue(cache_path, 'No count-cache signature; freshness unverified')
                merged.extend(rows)
            except (ValueError, KeyError, TypeError) as exc:
                c.issue(source, str(exc), 'error')
            for filename, destination in [('multi_reads_by_nh.tsv', nh), ('ribotx_internal_rescue_summary.tsv', rescue)]:
                p = directory/filename
                copied = c.capture(p, f'raw/mapping/{dataset}/{sample}/{filename}', 'Mapping categories')
                if copied:
                    destination.extend({**r, 'dataset': dataset, 'sample': sample, 'source': str(p)} for r in tsv_rows(copied))
    c.table('mapping_by_library.tsv', merged, 'Fig2 C-D / main mapping table; tool-defined unique/multi')
    c.table('mapping_nh.tsv', nh, 'Supplementary reported hit multiplicity')
    c.table('genome_rescue.tsv', rescue, 'Reported rescue counters; preserve producer definitions')
    summaries = []
    groups = defaultdict(list)
    for r in merged:
        groups[(r['dataset'], r['mapper'])].append(r)
    for (dataset, mapper), rows in groups.items():
        for metric in ('mapped_percent_common','unique_percent_common','multi_percent_common','residual_percent_common'):
            values = [r[metric] for r in rows if r[metric] is not None and r['eligible_for_summary']]
            summaries.append(dict(dataset=dataset, mapper=mapper, metric=metric,
                n=len(values), expected_n=len(datasets[dataset]), complete=len(values)==len(datasets[dataset]),
                mean=statistics.mean(values) if values else None,
                sample_sd=statistics.stdev(values) if len(values)>1 else None))
    c.table('mapping_summary.tsv', summaries, 'Means and sample SD for Fig2 C-D; no imputed missing libraries')


def collect_results(c, simulation, tpm_report, plot_dirs):
    for name in ('simulation_mapping_accuracy_summary.tsv', 'ribotx_multihit_accuracy_by_nh.tsv', 'star_multi_ribotx_accuracy_by_nh.tsv'):
        c.capture(simulation/name, f'raw/simulation/{name}', 'Fig2 E-F and simulation mapping evidence')
    for dataset in DATASETS:
        directory = simulation/dataset
        matches = sorted(p for p in directory.glob('*.tsv') if not p.name.startswith('truth.'))
        if not any('source_rescue_matching' in p.name for p in matches):
            c.issue(directory, 'compare_similation output missing: source_rescue_matching.tsv', 'error')
        for p in matches:
            c.capture(p, f'raw/simulation/{dataset}/{p.name}', 'Simulation comparison or input metadata')
        # Native output and merge membership needed to interpret reported IDs.
        for method in ('RiboTxAlign_LongRead', 'RiboTxAlign_RNASeq'):
            for suffix, subdir in (('transcript_tpm_with_source.tsv','em'), ('transcript_merge.tsv','em'), ('metadata.tsv','reference')):
                found = list((directory/method).glob(f'*/{subdir}/*.{suffix}'))
                if len(found) != 1:
                    c.issue(directory/method, f'Expected one {suffix}; found {len(found)}', 'error')
                for p in found:
                    c.capture(p, 'raw/native/'+str(p.relative_to(simulation)), 'Native TPM / reference / consolidation evidence')
    for name in ('overlap_summary.tsv', 'correlation_summary.tsv', 'scatter_source_data.tsv'):
        c.capture(tpm_report/name, 'raw/tpm/'+name, 'Fig3 canonical publication table; AND/OR populations remain distinct')
    for index, root in enumerate([tpm_report, *plot_dirs, *EXTRA_RESULT_DIRS]):
        if not root.is_dir():
            c.issue(root, 'Figure/result directory missing')
            continue
        for p in sorted(root.rglob('*')):
            if p.is_file() and p.suffix.lower() in ('.tsv','.csv','.svg','.pdf','.png','.xlsx','.txt'):
                c.capture(p, f'raw/result_directory_{index}/'+str(p.relative_to(root)), 'Existing figure and source data; not recalculated')


def collect_run(c, directory, dataset, sample, workflow, cache, expected_dedup):
    context = dict(dataset=dataset, sample=sample, workflow=workflow, run_directory=str(directory))
    resources, stages = [], []
    resource_root = directory/'resource'
    profiles = sorted(p for p in resource_root.glob('profile.*') if p.is_dir())
    root_measurement = (resource_root/'ribotxalign.resource_usage.txt').is_file()
    selected = RESOURCE_SESSION_OVERRIDES.get(str(directory))
    if selected is None:
        if profiles and root_measurement:
            selected = ''
            c.issue(resource_root, 'Both direct and profile measurements exist; set RESOURCE_SESSION_OVERRIDES to select the intended session', 'error')
        elif profiles:
            dated = [(re.fullmatch(r'profile\.(\d{8}T\d{6}Z)\..+', p.name), p) for p in profiles]
            if all(match for match, p in dated):
                latest = max(match[1] for match, p in dated)
                candidates = [p for match, p in dated if match[1] == latest]
                selected = candidates[0].name if len(candidates)==1 else ''
            else:
                selected = ''
            if selected:
                c.issue(resource_root, f'Latest profile attempt selected: {selected}; earlier attempts retained separately. Verify correspondence to final BAM.')
            else:
                c.issue(resource_root, 'Ambiguous profile timestamps; set RESOURCE_SESSION_OVERRIDES', 'error')
        else:
            selected = '.'
    if selected and selected != '.' and selected not in [p.name for p in profiles]:
        c.issue(resource_root, f'Selected profile does not exist: {selected}', 'error')
    resource_paths = sorted(resource_root.rglob('*.resource_usage.txt'))
    for p in resource_paths:
        if p.is_file() and p.suffix in ('.txt', '.log'):
            relative = p.relative_to(resource_root)
            session = str(relative.parent).replace('\\', '/')
            is_selected = session == selected
            copied = c.capture(p, f'raw/runs/{workflow}/{directory.name}/resource/{relative.as_posix()}', 'Original resource measurement')
            if copied:
                try:
                    phase = 'deduplicate' if 'dedup' in p.name.lower() else ('run' if 'ribotxalign' in p.name.lower() else p.stem)
                    for r in parse_resource(copied):
                        r.update(context, phase=phase, source=str(p), resource_session=session,
                                 selected_for_summary=is_selected)
                        resources.append(r)
                except (ValueError, KeyError) as exc:
                    c.issue(p, str(exc), 'error' if is_selected else 'warning')
    run_records = [r for r in resources if r['phase']=='run' and r['selected_for_summary']]
    if len(run_records) != 1:
        c.issue(directory/'resource', f'Expected one run measurement, found {len(run_records)}', 'error')
    for r in resources:
        if not r['complete'] or not r['successful'] or r['invocations_in_file'] != 1:
            c.issue(r['source'], 'Incomplete, failed or multiple invocation record; excluded from resource summary', 'error' if r['selected_for_summary'] else 'warning')
    # Retain process telemetry as raw evidence, without conflating it with GNU
    # time's maximum RSS or summing process peaks from different timestamps.
    performance = directory/'logs'/'performance'
    for p in sorted(performance.rglob('*')):
        if p.is_file() and p.suffix.lower() in ('.tsv', '.json', '.txt', '.log'):
            c.capture(p, f'raw/runs/{workflow}/{directory.name}/logs/performance/{p.relative_to(performance).as_posix()}',
                      'Per-process profiling evidence; separate from GNU time peak RSS')
    # Preserve log time entries individually: nested/parallel durations are not additive.
    for subdir in ('em', 'logs', 'mapping'):
        for p in sorted((directory/subdir).glob('*.log')):
            if p.stat().st_size > MAX_ARTIFACT_BYTES:
                c.issue(p, 'Large log omitted', 'error')
                continue
            copied = c.capture(p, f'raw/runs/{workflow}/{directory.name}/{subdir}/{p.name}', 'Stage timing evidence')
            if copied:
                with copied.open(encoding='utf-8', errors='replace') as handle:
                    for line_number, line in enumerate(handle, 1):
                        match = re.search(r'elapsed\s*[:=]\s*([\d.]+)\s*s\b', line, re.I)
                        if match:
                            stages.append({**context, 'source': str(p), 'line': line_number,
                                           'elapsed_seconds': float(match[1]), 'log_text': line.strip()})
    for p in sorted((directory/'mapping').glob('*.mapping_stats.json')):
        c.capture(p, f'raw/runs/{workflow}/{directory.name}/mapping/{p.name}', 'Pipeline mapping counters')

    dedup_records = [r for r in resources if r['phase']=='deduplicate' and r['selected_for_summary']]
    stats_paths = sorted((directory/'reads').glob('*.stats.tsv'))
    stats_paths = [p for p in stats_paths if 'dedup' in p.name.lower()]
    row = {**context, 'status': 'missing', 'source': ''}
    try:
        if len(stats_paths) > 1:
            raise ValueError('Multiple dedup stats files; choose the intended file before collection')
        if stats_paths:
            p = stats_paths[0]
            copied = c.capture(p, f'raw/runs/{workflow}/{directory.name}/reads/{p.name}', 'Dedup before/after counts and filters')
            if copied:
                row.update(parse_dedup(copied), status='dedup_stats', source=str(p))
                output = Path(str(p).removesuffix('.stats.tsv'))
                if output.is_file() and output.stat().st_mtime_ns > p.stat().st_mtime_ns:
                    row['status'] = 'unverified_stats'
                    c.issue(p, 'Output FASTQ newer than stats; rerun may have replaced reads', 'error')
        elif expected_dedup or dedup_records:
            c.issue(directory/'reads', 'Dedup stats missing; FASTQ difference includes filtering, not just duplicates')
            if COUNT_FASTQ_IF_NEEDED and len(dedup_records) == 1 and dedup_records[0]['successful']:
                command = dedup_records[0]['command']
                before = command_option(command, '--fastq')
                after = command_option(command, '--output-fastq')
                if before and after:
                    n_before, n_after = count_fastq(before, cache), count_fastq(after, cache)
                    if n_after > n_before:
                        raise ValueError('Dedup output count exceeds input count')
                    row.update(status='fastq_counts_only', input_reads=n_before, output_reads=n_after,
                        total_reduction_reads=n_before-n_after,
                        total_reduction_percent=100*(n_before-n_after)/n_before if n_before else None,
                        source=f'{before} ; {after}')
            if row['status']=='missing':
                c.issue(directory, 'Dedup before/after counts unavailable', 'error')
        else:
            row.update(status='not_requested_by_supplied_shell', input_reads=None, output_reads=None)
        if expected_dedup and not dedup_records:
            c.issue(directory/'resource', 'Dedup time/RSS not recorded; cannot reconstruct retrospectively')
    except (ValueError, KeyError, OSError) as exc:
        row['status'] = 'invalid_or_missing'
        c.issue(directory, f'Dedup: {exc}', 'error')
    if len(run_records) == 1:
        row['mapping_input_fastq'] = command_option(run_records[0]['command'], '--fastq')
    if len(dedup_records) == 1:
        row['dedup_input_fastq'] = command_option(dedup_records[0]['command'], '--fastq')
        row['dedup_output_fastq'] = command_option(dedup_records[0]['command'], '--output-fastq')
        if row.get('mapping_input_fastq') and row.get('dedup_output_fastq'):
            row['mapping_uses_dedup_output_path'] = row['mapping_input_fastq'] == row['dedup_output_fastq']
            if not row['mapping_uses_dedup_output_path']:
                c.issue(directory, 'Mapping input path differs from dedup output path; do not substitute dedup counts as mapping denominator', 'error')
    return resources, stages, row


def resource_summary(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row['dataset'], row['workflow'], row['phase'])].append(row)
    result = []
    for (dataset, workflow, phase), records in groups.items():
        sample_counts = defaultdict(int)
        for r in records:
            if r.get('selected_for_summary', True):
                sample_counts[r['sample']] += 1
        kept = [r for r in records if r.get('selected_for_summary', True) and r['complete'] and r['successful'] and r['invocations_in_file']==1 and sample_counts[r['sample']]==1]
        values = [r['elapsed_minutes'] for r in kept]
        rss = [r['max_rss_gib'] for r in kept]
        result.append(dict(dataset=dataset, workflow=workflow, phase=phase, records=len(records),
            included_runs=len(kept), excluded_records=len(records)-len(kept),
            expected_runs=len(DATASETS[dataset]) if dataset in DATASETS else None,
            all_expected_runs_included=len(kept)==len(DATASETS[dataset]) if dataset in DATASETS else None,
            mean_minutes=statistics.mean(values) if values else None,
            sd_minutes=statistics.stdev(values) if len(values)>1 else None,
            min_peak_rss_gib=min(rss) if rss else None, max_peak_rss_gib=max(rss) if rss else None))
    return result


def collect_scripts(c):
    names = ['test/validation/Bam4Comapre_1.py', 'test/similation/compare_similation_1.py',
             'test/similation/evaluate_mapped_bam_1.py', 'test/similation/plot_graph_1.py',
             'test/similation/plotgraph_tpm_1.py']
    paths = []
    for name in names:
        source = PROJECT/name
        copied = c.capture(source, 'scripts/'+name, 'Requested validation source version at collection time')
        if copied:
            try:
                paths.extend(dict(script=name, constant=k, path=str(v)) for k,v in static_paths(copied).items())
            except (SyntaxError, ValueError) as exc:
                c.issue(source, f'Cannot inspect path constants: {exc}')
    for name in ('RiboTxAlign_human.yaml','RiboTxAlign_mouse.yaml','pyproject.toml'):
        c.capture(PROJECT/name, 'scripts/'+name, 'Collection-time settings; not proof of historical run version')
    c.capture(Path(__file__), 'scripts/collect_paper_results.py', 'Collection program')
    c.table('script_paths.tsv', paths, 'Statically readable paths from *_1.py; scripts are not executed')
    expected = {
        ('test/validation/Bam4Comapre_1.py', 'OUTPUT_DIR'): MAPPING_COUNTS,
        ('test/similation/evaluate_mapped_bam_1.py', 'BASE'): SIMULATION,
        ('test/similation/plotgraph_tpm_1.py', 'REPORT_DIR'): TPM_REPORT,
        ('test/similation/plot_graph_1.py', 'BASE'): SIMULATION,
    }
    for dataset in DATASETS:
        expected[('test/similation/compare_similation_1.py', f'DATASET_PATHS.{dataset}.output_dir')] = SIMULATION/dataset
    for row in paths:
        want = expected.get((row['script'], row['constant']))
        if want is not None and Path(row['path']) != want:
            c.issue(row['script'], f"Output path differs: {row['constant']}={row['path']}; edit collector configuration", 'error')
        if row['script']=='test/similation/plot_graph_1.py' and row['constant']=='OUTPUT_DIR':
            if Path(row['path']) not in [*PLOT_DIRS, *EXTRA_RESULT_DIRS]:
                c.issue(row['script'], f"Add plot directory to PLOT_DIRS: {row['path']}", 'error')


def main():
    if not VALIDATION.is_dir() or not SIMULATION.is_dir():
        raise SystemExit('Server paths not found. Run this script on the benchmark server, or edit its path constants.')
    OUTPUT_PARENT.mkdir(parents=True, exist_ok=True)
    output = OUTPUT_PARENT / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    c = Collector(output)
    cache_path = OUTPUT_PARENT/'fastq_count_cache.json'
    try:
        cache = json.loads(cache_path.read_text()) if cache_path.is_file() else {}
    except (ValueError, OSError):
        cache = {}
    collect_scripts(c)
    collect_mapping(c, MAPPING_COUNTS, DATASETS, MAPPERS)
    collect_results(c, SIMULATION, TPM_REPORT, PLOT_DIRS)
    resources, stages, dedup = [], [], []
    for workflow, root in RUN_ROOTS.items():
        for dataset, samples in DATASETS.items():
            for sample in samples:
                directory = root/f'{dataset}_{sample}.clean'
                print(f'Collecting {workflow}: {directory.name}', flush=True)
                r, s, d = collect_run(c, directory, dataset, sample, workflow, cache, DEDUP_EXPECTED[dataset])
                resources.extend(r); stages.extend(s); dedup.append(d)
    # Simulation timing is kept apart from real-library timing, when measured.
    for dataset in DATASETS:
        for method, suffix in [('RiboTxAlign_LongRead','LongRead'), ('RiboTxAlign_RNASeq','RNASeq')]:
            directory = SIMULATION/dataset/method/f'{dataset}_simulation_{suffix}'
            if (directory/'resource').is_dir():
                r, s, d = collect_run(c, directory, dataset+'_simulation', dataset+'_simulation',
                                      method, cache, False)
                resources.extend(r); stages.extend(s); dedup.append(d)
            else:
                c.issue(directory/'resource', 'Simulation resource logs unavailable; no time/RSS inferred')
    for i, p in enumerate(EXTRA_RESOURCE_FILES):
        copied = c.capture(p, f'raw/extra_resources/{i}_{p.name}', 'Additional unclassified GNU time record')
        if copied:
            try:
                for r in parse_resource(copied):
                    resources.append({**r, 'dataset':'unclassified', 'sample':str(p), 'workflow':'unclassified',
                                      'phase':'unclassified', 'source':str(p)})
            except ValueError as exc:
                c.issue(p, str(exc), 'error')
    cache_temp = cache_path.with_suffix('.tmp')
    cache_temp.write_text(json.dumps(cache, indent=2), encoding='utf-8')
    cache_temp.replace(cache_path)
    c.table('resource_runs.tsv', resources, 'GNU time per invocation: wall time / peak RSS; phases not summed')
    c.table('resource_summary.tsv', resource_summary(resources), 'Successful unambiguous runs only; RSS range is not simultaneous process total')
    c.table('stage_times.tsv', stages, 'Individual log timing lines; repeated/nested/parallel stages not summed')
    c.table('deduplication.tsv', dedup, 'Input/output, filtered and duplicate counts; missing values are not zero')
    c.table('collection_issues.tsv', c.issues, 'Resolve missing/stale/ambiguous inputs before publication')
    write_tsv(output/'file_manifest.tsv', c.manifest)
    write_tsv(output/'contents.tsv', c.tables)
    (output/'collection.json').write_text(json.dumps(dict(
        collected_utc=datetime.now(timezone.utc).isoformat(), files=len(c.manifest),
        complete=not any(r['severity']=='error' for r in c.issues),
        errors=sum(r['severity']=='error' for r in c.issues), warnings=sum(r['severity']=='warning' for r in c.issues),
        note='Snapshot of completed results; no mapper or evaluator was rerun.'), indent=2), encoding='utf-8')
    (output/'README.txt').write_text(
        'Read contents.tsv and collection_issues.tsv first.\n'
        'Mapping rates preserve maximum observed mapper total as denominator.\n'
        'Dedup counts do not automatically prove all mappers used the same FASTQ.\n'
        'Do not combine compare_similation metrics with plotgraph_tpm AND/OR metrics.\n'
        'Missing values remain blank. True/False fields are booleans, not checkmarks.\n'
        'RSS is GNU time maximum RSS, not the concurrent process-tree total.\n'
        'Stage timings can overlap. Reference preparation has no retrospective measurement unless separately logged.\n'
        'skipped_short in the current deduplicator includes missing/invalid UMI extraction.\n'
        'total_reduction includes both filtering and duplicate removal; duplicate removal uses written_records-output_reads.\n'
        'Copied scripts/configs describe collection-time files, not necessarily historical benchmark versions.\n', encoding='utf-8')
    archive = output.with_suffix('.zip')
    with ZipFile(archive, 'w', ZIP_DEFLATED) as z:
        for p in sorted(output.rglob('*')):
            if p.is_file():
                z.write(p, p.relative_to(output))
    print(f'\nCollected results: {output}\nTransfer this ZIP: {archive}\nIssues: {output / "collection_issues.tsv"}', flush=True)
    if any(r['severity']=='error' for r in c.issues):
        print('Collection is PARTIAL. ZIP was saved; inspect collection_issues.tsv before using its numbers.')
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())