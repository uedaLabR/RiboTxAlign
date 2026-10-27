#!/usr/bin/env python3
"""Collect HEK293T/mouse_brain reference and transcript merge evidence (stdlib only)."""
import argparse
import csv
from collections import Counter
from pathlib import Path
import re
import sys


FIELDS = ['dataset', 'sample', 'reference_kind', 'stage', 'metric', 'value',
          'unit', 'status', 'source', 'line', 'note']


def dataset_for(name):
    if name.lower().startswith(('hek293_', 'hek293t_')):
        return 'HEK293T'
    if name.lower().startswith('mouse_brain_'):
        return 'mouse_brain'
    return None


def row(context, stage, metric, value, source, line='', unit='transcripts',
        status='observed', note=''):
    return dict(zip(FIELDS, [*context, stage, metric, value, unit, status,
                             str(source), line, note]))


def parse_reference(path, context):
    """Keep executions separate; never sum repeated executions or duplicate log formats."""
    lines = path.read_text(encoding='utf-8', errors='replace').splitlines()
    starts = [i for i, s in enumerate(lines) if s.strip() == '1. Read GENCODE GFF3.']
    if not starts:
        return [row(context, 'reference', 'log_parse', '', path, status='unavailable',
                    note='No reference construction start found')]
    result = []
    for run, start in enumerate(starts):
        end = starts[run + 1] if run + 1 < len(starts) else len(lines)
        block = lines[start:end]
        ctx = list(context)
        # A shell-wide redirected log can contain several samples.
        prefixes = [s.split('Reference prefix:', 1)[1].strip() for s in block
                    if 'Reference prefix:' in s]
        if prefixes:
            sample = prefixes[-1].replace('\\', '/').rsplit('/', 1)[-1]
            dataset = dataset_for(sample)
            if not dataset:
                continue
            ctx[0:2] = [dataset, sample]
        complete = any('Reference construction flow completed.' in s for s in block)
        status = 'observed' if complete else 'incomplete'
        note = 'execution_in_log=%d; do not sum across executions' % (run + 1)
        def emit(stage, metric, value, index='', unit='transcripts', extra=''):
            result.append(row(ctx, stage, metric, value, path, index, unit,
                              status, note + ('; ' + extra if extra else '')))
        patterns = {
            'gencode_input': r'^GENCODE transcripts:\s*([\d,]+)',
            'stringtie_input': r'^StringTie2 transcripts:\s*([\d,]+)',
            'stringtie_records': r'^StringTie transcript records:\s*([\d,]+)',
            'unknown_strand_excluded': r'^Unknown-strand transcripts excluded:\s*([\d,]+)',
            'final_reference_transcripts': r'^Final reference transcripts:\s*([\d,]+)',
        }
        values = {}
        before = after = removed = 0
        collapse_count = 0
        rules = Counter()
        in_rules = False
        for index, line in enumerate(block, start + 1):
            for key, pattern in patterns.items():
                match = re.search(pattern, line)
                if match:
                    value = int(match[1].replace(',', ''))
                    values[key] = value
                    emit('reference', key, value, index)
            # Ignore the second, chromosome-labelled summary of the same collapse.
            match = re.match(r'^StringTie conservative collapse: before=\s*(\d+) '
                             r'after=\s*(\d+) removed=\s*(\d+)', line)
            if match:
                b, a, r = map(int, match.groups())
                if b - a != r or r < 0:
                    raise ValueError('Invalid collapse counts at %s:%d' % (path, index))
                before += b
                after += a
                removed += r
                collapse_count += 1
            if line.strip() == 'StringTie collapse groups by rule:':
                in_rules = True
                continue
            if in_rules:
                match = re.fullmatch(r'\s+(stringtie_\S+)\s+(\d+)\s*', line)
                if match:
                    rules[match[1]] += int(match[2])
                else:
                    in_rules = False
            if line.startswith('FASTA transcripts:'):
                for key, val in re.findall(r'(\w+):\s*(\d+)', line):
                    emit('reference_fasta', key, int(val), index,
                         'bytes' if key == 'bytes' else 'bases' if 'bases' in key else 'records')
        if collapse_count:
            if complete and before != values.get('stringtie_input'):
                raise ValueError('Collapse coverage differs from StringTie input: %s' % path)
            for metric, value in [('input', before), ('output', after), ('merged_away', removed)]:
                emit('stringtie_conservative_collapse', metric, value)
            for rule, count in sorted(rules.items()):
                emit('stringtie_conservative_collapse', rule, count, unit='merge_groups',
                     extra='Groups, not the number of transcripts removed')
            if all(k in values for k in ['gencode_input', 'final_reference_transcripts']):
                emit('reference_resolution', 'net_reduction',
                     values['gencode_input'] + after - values['final_reference_transcripts'],
                     extra='Includes filtering/selection and merging; NOT a merge count')
        else:
            result.append(row(ctx, 'stringtie_conservative_collapse', 'merged_away', '',
                              path, status='unavailable', note=note + '; no collapse summary'))
        result.append(row(ctx, 'reference_resolution', 'stage_specific_merged_away', '',
                          path, status='unavailable',
                          note=note + '; existing logs do not separate same-CDS merging, '
                          'low-TPM absorption, micro-exon filtering and reference selection'))
    return result


def read_table(path, required):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError('Missing columns in %s: %s' % (path, required))
        return list(reader)


def parse_merges(path, context):
    records = read_table(path, ['group_id', 'source_transcript_id', 'target_transcript_id'])
    links = {}
    for record in records:
        group, source, target = (record[k] for k in
                                ['group_id', 'source_transcript_id', 'target_transcript_id'])
        if not all((group, source, target)) or source == target:
            raise ValueError('Invalid merge relationship in %s' % path)
        key = (group, source)
        if key in links and links[key] != target:
            raise ValueError('Conflicting merge targets in %s' % path)
        links[key] = target
    families = {}
    for (group, source), target in links.items():
        seen = {source}
        while (group, target) in links:
            if target in seen:
                raise ValueError('Cyclic merge relationship in %s' % path)
            seen.add(target)
            target = links[group, target]
        families.setdefault((group, target), set()).update([source, target])
    metrics = [('merged_away', len(links), 'transcripts'),
               ('merge_groups', len(families), 'merge_groups'),
               ('participating_transcripts', sum(map(len, families.values())), 'transcripts'),
               ('duplicate_rows_ignored', len(records) - len(links), 'rows')]
    em_log = path.with_name(path.name.replace('.transcript_merge.tsv', '.transcript_em.log'))
    completed = (em_log.is_file() and 'Transcript EM total elapsed:' in
                 em_log.read_text(encoding='utf-8', errors='replace'))
    return [row(context, 'transcript_em', key, val, path, unit=unit,
                status='observed' if completed else 'completion_unverified',
                note='Current TSV snapshot, not linked to archived reference runs; '
                     'EM completion marker ' + ('found' if completed else 'not found') +
                     ' in ' + str(em_log))
            for key, val, unit in metrics]


def parse_metadata(path, context):
    records = read_table(path, ['reference_id', 'resolution_type', 'transcript_length',
                                'source_transcript_ids'])
    if len({r['reference_id'] for r in records}) != len(records):
        raise ValueError('Duplicate reference IDs in %s' % path)
    counts = Counter(r['resolution_type'] or 'unspecified' for r in records)
    result = [row(context, 'reference_metadata', 'final_reference_transcripts', len(records), path),
              row(context, 'reference_metadata', 'transcript_bases',
                  sum(int(r['transcript_length']) for r in records), path, unit='bases')]
    for label, count in sorted(counts.items()):
        result.append(row(context, 'reference_resolution_type', label, count, path,
                          note='Final reference categories, not independent merge events'))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=Path('/mnt/share/project/riboML/validation'))
    parser.add_argument('--output-dir', type=Path, default=Path('transcript_merge_stats'))
    parser.add_argument('--reference-log', action='append', default=[], metavar='KIND=PATH',
                        help='Additional redirected shell log; KIND is longread or shortread. '
                             'Sample is read from Reference prefix lines. Repeatable.')
    args = parser.parse_args(argv)
    rows = []
    seen = set()
    for kind, directory in [('longread', 'ribotx_validation'),
                            ('shortread', 'ribotx_validation_rnaseq')]:
        for sample_dir in sorted((args.base / directory).glob('*')):
            dataset = dataset_for(sample_dir.name)
            if not sample_dir.is_dir() or not dataset:
                continue
            context = [dataset, sample_dir.name, kind]
            logs = sorted(sample_dir.glob('resource/**/create_reference.console.log'))
            for log in logs:
                rows.extend(parse_reference(log, context))
                seen.add(log.resolve())
            if not logs:
                rows.append(row(context, 'reference', 'console_log', '', sample_dir,
                                status='missing', note='Use --reference-log for redirected shell stdout'))
            for suffix, folder, function, stage in [
                ('metadata.tsv', 'reference', parse_metadata, 'reference_metadata'),
                ('transcript_merge.tsv', 'em', parse_merges, 'transcript_em'),
            ]:
                path = sample_dir / folder / (sample_dir.name + '.' + suffix)
                if path.is_file():
                    rows.extend(function(path, context))
                else:
                    rows.append(row(context, stage, 'input_file', '', path, status='missing'))
    for spec in args.reference_log:
        kind, sep, name = spec.partition('=')
        if not sep or kind not in ('longread', 'shortread'):
            parser.error('--reference-log must be longread=PATH or shortread=PATH')
        path = Path(name)
        if path.resolve() not in seen:
            rows.extend(parse_reference(path, ['', '', kind]))
            seen.add(path.resolve())
    if not rows:
        parser.error('No HEK293T/mouse_brain sample directories or reference logs found')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / 'transcript_merge_statistics.tsv'
    with output.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print('Wrote %s (%d rows)' % (output, len(rows)))
    missing = sum(r['status'] != 'observed' for r in rows)
    print('Missing/unavailable/incomplete rows: %d (see status and note columns)' % missing)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as error:
        sys.exit(str(error))
