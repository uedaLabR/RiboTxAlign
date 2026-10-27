#!/usr/bin/env python3
"""Replay reference construction in memory; collect two-stage merge counts.

Uses the installed project's reference functions, without changing its files or
writing references. Existing per-sample metadata must match the replay before
reference-stage counts are paired with that sample's existing EM merge table.
"""
import argparse
from collections import Counter
from contextlib import ExitStack, redirect_stdout
import csv
import gzip
import hashlib
import inspect
import json
import math
from pathlib import Path
import re
import sys
from unittest.mock import patch


SUMMARY_FIELDS = ['dataset', 'sample', 'reference_type', 'stringtie_input',
                  'reference_construction_merged_away', 'reference_filtered_out',
                  'gencode_models_added', 'final_reference_transcripts',
                  'before_em2_merged_away', 'status', 'metadata_source',
                  'em_merge_source']
STEPS = ['stringtie_precollapse', 'matched_cds_merge', 'unmatched_subset_merge',
         'final_cds_merge', 'microexon_filter', 'spanning_exon_filter',
         'low_tpm_filter', 'missing_start_filter', 'gencode_additions']


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def measured_build(gencode, stringtie):
    """Temporary wrappers count actual length changes, not metadata memberships."""
    from db_construct import overlap_group as grouping
    from db_construct import overlap_group_worker as worker
    from db_construct import transcript_resolver as resolver
    counts = Counter({step: 0 for step in STEPS})
    original = {name: getattr(resolver, name) for name in [
        '_collapse_cds_duplicates', '_collapse_exon_subset_duplicates',
        '_filter_low_tpm_within_gene', '_add_gencode_representative_when_uncovered',
        '_has_no_micro_exons', '_spans_too_many_gencode_exons']}
    collapse = worker.collapse_stringtie_subsets

    def precollapse(*args, **kwargs):
        inputs = kwargs['stringtie_transcripts'] if 'stringtie_transcripts' in kwargs else args[0]
        before = len(inputs)
        result = collapse(*args, **kwargs)
        counts['stringtie_precollapse'] += before - len(result)
        counts['missing_start_filter'] += sum(t.start is None for t in result)
        return result

    def cds(*args, **kwargs):
        before = len(kwargs['candidates'] if 'candidates' in kwargs else args[0])
        final = kwargs.get('prefer_gencode', args[1] if len(args) > 1 else False)
        result = original['_collapse_cds_duplicates'](*args, **kwargs)
        counts['final_cds_merge' if final else 'matched_cds_merge'] += before - len(result)
        return result

    def sized(name, metric, parameter, additions=False):
        def wrapped(*args, **kwargs):
            before = len(kwargs[parameter] if parameter in kwargs else args[0])
            result = original[name](*args, **kwargs)
            counts[metric] += (len(result) - before) if additions else (before - len(result))
            return result
        return wrapped

    def micro(*args, **kwargs):
        result = original['_has_no_micro_exons'](*args, **kwargs)
        counts['microexon_filter'] += int(not result)
        return result

    def spanning(*args, **kwargs):
        result = original['_spans_too_many_gencode_exons'](*args, **kwargs)
        counts['spanning_exon_filter'] += int(result)
        return result

    with ExitStack() as stack:
        stack.enter_context(patch.object(worker, 'collapse_stringtie_subsets', precollapse))
        for name, function in {
            '_collapse_cds_duplicates': cds,
            '_collapse_exon_subset_duplicates': sized('_collapse_exon_subset_duplicates', 'unmatched_subset_merge', 'transcripts'),
            '_filter_low_tpm_within_gene': sized('_filter_low_tpm_within_gene', 'low_tpm_filter', 'transcripts'),
            '_add_gencode_representative_when_uncovered': sized('_add_gencode_representative_when_uncovered', 'gencode_additions', 'resolved_transcripts', True),
            '_has_no_micro_exons': micro,
            '_spans_too_many_gencode_exons': spanning,
        }.items():
            stack.enter_context(patch.object(resolver, name, function))
        result = grouping.build_reference_transcripts(gencode, stringtie)
    merged = sum(counts[s] for s in STEPS[:4])
    filtered = sum(counts[s] for s in STEPS[4:8])
    expected = len(stringtie) - merged - filtered + counts['gencode_additions']
    if expected != len(result) or any(v < 0 for v in counts.values()):
        raise ValueError('Reference count balance failed: input - merges - filters + additions '
                         f'= {expected}, observed {len(result)}; code may be incompatible')
    return result, dict(counts), merged, filtered


def metadata_equal(actual, expected):
    """Compare all saved metadata fields, allowing numeric serialization of TPM."""
    if len(actual) != len(expected):
        return False
    keyed = {r['reference_id']: r for r in actual}
    if len(keyed) != len(actual) or set(keyed) != set(expected):
        return False
    for rid, reference in expected.items():
        row = keyed[rid]
        for key, value in reference.items():
            if key not in row:
                return False
            if key == 'tpm':
                if not math.isclose(float(row[key]), float(value), rel_tol=1e-9, abs_tol=1e-12):
                    return False
            elif row[key] != str(value):
                return False
    return True


def count_em(path, reference_ids):
    records = read_tsv(path)
    # Header-only TSV is an observed zero only when the caller verifies completion.
    with path.open(encoding='utf-8-sig') as f:
        if f.readline().strip().split('\t') != ['group_id', 'source_transcript_id', 'target_transcript_id']:
            raise ValueError(f'Unexpected EM merge columns: {path}')
    links = {}
    for row in records:
        source, target = row['source_transcript_id'], row['target_transcript_id']
        if source not in reference_ids or target not in reference_ids or source == target:
            raise ValueError(f'Invalid EM merge IDs in {path}')
        pair = (row['group_id'], target)
        if source in links and links[source] != pair:
            raise ValueError(f'Conflicting EM merge source in {path}')
        links[source] = pair
    for source in links:
        seen = set()
        cursor = source
        while cursor in links:
            if cursor in seen:
                raise ValueError(f'Cyclic EM merge in {path}')
            seen.add(cursor)
            cursor = links[cursor][1]
    return len(links)


def metadata_differences(actual, expected, limit=None):
    """Report why a replay cannot be paired with the historical EM results."""
    keyed = {r['reference_id']: r for r in actual}
    differences = []
    for rid in sorted(set(keyed) | set(expected)):
        a, e = keyed.get(rid), expected.get(rid)
        if a is None or e is None:
            differences.append(dict(reference_id=rid, field='reference_id',
                                    saved='present' if a is not None else 'missing',
                                    replay='present' if e is not None else 'missing'))
        else:
            for field, value in e.items():
                equal = str(value) == a.get(field)
                if field == 'tpm' and field in a:
                    try:
                        equal = math.isclose(float(a[field]), float(value), rel_tol=1e-9, abs_tol=1e-12)
                    except ValueError:
                        equal = False
                if not equal:
                    differences.append(dict(reference_id=rid, field=field, saved=a.get(field, '<missing column>'), replay=value))
                if limit is not None and len(differences) >= limit:
                    return differences
        if limit is not None and len(differences) >= limit:
            return differences
    return differences


def default_jobs(project, base):
    import yaml
    configs = {}
    for species in ['human', 'mouse']:
        with (project / f'RiboTxAlign_{species}.yaml').open() as f:
            configs[species] = yaml.safe_load(f)['reference']['gencode_gff3']
    definitions = [
        ('HEK293T', 'longread', 'HEK293T/longread/mapping/HEK293T.longread.stringtie.gtf', 'human'),
        ('HEK293T', 'shortread', 'HEK293T/rnaseq/stringtie/HEK293T_DMSO.covered_reference.gtf', 'human'),
        ('mouse_brain', 'longread', 'mouse_brain/longread/mapping/mouse_brain.longread.stringtie.gtf', 'mouse'),
        ('mouse_brain', 'shortread', 'mouse_brain/rnaseq/mapping/mouse_brain.rnaseq.stringtie.gtf', 'mouse'),
    ]
    return [dict(dataset=d, reference_type=k, stringtie_gtf=str(base / gtf),
                 gencode_gff3=configs[species],
                 output_root=str(base / ('ribotx_validation' if k == 'longread' else 'ribotx_validation_rnaseq')))
            for d, k, gtf, species in definitions]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project', type=Path, default=Path('/share/ueda/RiboTxAlign'))
    p.add_argument('--base', type=Path, default=Path('/mnt/share/project/riboML/validation'))
    p.add_argument('--manifest', type=Path, help='Optional TSV: dataset, reference_type, stringtie_gtf, gencode_gff3, output_root')
    p.add_argument('--output-dir', type=Path, help='Default: two_stage_merge_stats beside this script')
    args = p.parse_args()
    out = (args.output_dir or Path(__file__).resolve().parent / 'two_stage_merge_stats').resolve()
    project = args.project.resolve()
    sys.path.insert(0, str(project))
    from db_construct.gtf_reader import load_gencode_gff3, load_stringtie_transcripts
    from db_construct.transcript_sort import sort_transcripts
    from db_construct.metadata_writer import build_metadata_row
    from db_construct import overlap_group
    if Path(inspect.getfile(overlap_group)).resolve().parent.parent != project:
        raise ValueError('Imported reference code is not from --project')
    jobs = read_tsv(args.manifest) if args.manifest else default_jobs(project, args.base)
    if not jobs:
        p.error('No jobs')
    keys = [(j['dataset'], j['reference_type']) for j in jobs]
    if len(set(keys)) != len(keys):
        p.error('Duplicate dataset/reference_type jobs')
    for job in jobs:
        if job['dataset'] not in ['HEK293T', 'mouse_brain'] or job['reference_type'] not in ['longread', 'shortread']:
            p.error('Only HEK293T/mouse_brain and longread/shortread are supported')
        for key in ['stringtie_gtf', 'gencode_gff3', 'output_root']:
            if not Path(job[key]).exists():
                p.error(f'Missing {key}: {job[key]}; supply corrected paths with --manifest')
    out.mkdir(parents=True, exist_ok=True)
    hashes = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((project / 'db_construct').glob('*.py'))}
    (out / 'provenance.json').write_text(json.dumps(dict(project=str(project), jobs=jobs, code_sha256=hashes), indent=2), encoding='utf-8')
    summary, details, replay_rows, difference_counts = [], [], [], []
    for job in jobs:
        dataset, kind = job['dataset'], job['reference_type']
        log = out / f'{dataset}_{kind}.reference_replay.log'
        print(f'Replaying reference only: {dataset} {kind}', flush=True)
        with log.open('w', encoding='utf-8') as f, redirect_stdout(f):
            gencode = sort_transcripts(load_gencode_gff3(job['gencode_gff3']))
            stringtie = sort_transcripts(load_stringtie_transcripts(job['stringtie_gtf']))
            models, counts, merged, filtered = measured_build(gencode, stringtie)
        expected = {r['reference_id']: r for r in map(build_metadata_row, models)}
        # Preserve the complete replay for later comparisons without rebuilding.
        with gzip.open(out / f'{dataset}_{kind}.replay_metadata.tsv.gz', 'wt', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, list(next(iter(expected.values()))), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(expected.values())
        replay_rows.append(dict(dataset=dataset, reference_type=kind,
                                stringtie_input=len(stringtie),
                                reference_construction_merged_away=merged,
                                reference_filtered_out=filtered,
                                gencode_models_added=counts['gencode_additions'],
                                final_reference_transcripts=len(models),
                                status='current_code_replay_not_historical_measurement'))
        write_tsv(out / 'reference_replay_summary.tsv', replay_rows, list(replay_rows[0]))
        if len(expected) != len(models):
            raise ValueError('Replay produced duplicate reference IDs')
        for step, count in counts.items():
            details.append(dict(dataset=dataset, reference_type=kind, step=step, count=count, replay_log=str(log)))
        directories = [s for s in sorted(Path(job['output_root']).glob(dataset + '_SRR*.clean')) if s.is_dir()]
        if not directories:
            raise ValueError(f'No current sample directories: {job}')
        for sample in directories:
            metadata = sample / 'reference' / (sample.name + '.metadata.tsv')
            em = sample / 'em' / (sample.name + '.transcript_merge.tsv')
            em_log = sample / 'em' / (sample.name + '.transcript_em.log')
            status = []
            saved = read_tsv(metadata) if metadata.is_file() else []
            match = metadata.is_file() and metadata_equal(saved, expected)
            if not match:
                status.append('reference_metadata_missing_or_mismatch')
                differences = metadata_differences(saved, expected)
                if not differences:
                    differences = [dict(reference_id='',field='row_count_or_duplicate_ids',saved=len(saved),replay=len(expected))]
                field_counts = Counter(d['field'] for d in differences)
                examples, seen = [], Counter()
                for difference in differences:
                    field = difference['field']
                    if seen[field] < 5:
                        examples.append(difference)
                    seen[field] += 1
                difference_counts.extend(dict(dataset=dataset, sample=sample.name, reference_type=kind,
                                              field=field, difference_count=n)
                                         for field, n in sorted(field_counts.items()))
                write_tsv(out / 'metadata_diff_counts.tsv', difference_counts,
                          ['dataset','sample','reference_type','field','difference_count'])
                write_tsv(out / f'{sample.name}_{kind}.metadata_diff.tsv', examples,
                          ['reference_id','field','saved','replay'])
                print(f'  Metadata mismatch: {sample.name}; saved rows={len(saved)}, replay rows={len(expected)}', flush=True)
            em_complete = em_log.is_file() and 'Transcript EM total elapsed:' in em_log.read_text(encoding='utf-8', errors='replace')
            if not em_complete:
                status.append('em_completion_unverified')
            if not em.is_file():
                status.append('em_merge_table_missing')
            saved_ids = {r['reference_id'] for r in saved}
            saved_valid = bool(saved) and len(saved_ids) == len(saved)
            em_count = count_em(em, saved_ids) if saved_valid and em_complete and em.is_file() else ''
            summary.append(dict(dataset=dataset, sample=sample.name, reference_type=kind,
                                stringtie_input=len(stringtie) if match else '',
                                reference_construction_merged_away=merged if match else '',
                                reference_filtered_out=filtered if match else '',
                                gencode_models_added=counts['gencode_additions'] if match else '',
                                final_reference_transcripts=len(saved) if saved_valid else '',
                                before_em2_merged_away=em_count,
                                status=';'.join(status) or 'metadata_matched_replay_and_completed_em',
                                metadata_source=str(metadata), em_merge_source=str(em)))
        write_tsv(out / 'two_stage_merge_summary.tsv', summary, SUMMARY_FIELDS)
        write_tsv(out / 'reference_merge_details.tsv', details, ['dataset','reference_type','step','count','replay_log'])
        print(f'  StringTie {len(stringtie):,} - merges {merged:,} - filters {filtered:,} '
              f'+ GENCODE additions {counts["gencode_additions"]:,} = reference {len(models):,}', flush=True)
    print(f'Saved: {out / "two_stage_merge_summary.tsv"}')
    bad = sum(r['status'] != 'metadata_matched_replay_and_completed_em' for r in summary)
    if bad:
        print(f'{bad} rows require inspection; unavailable counts were left blank.', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, AttributeError) as exc:
        sys.exit(str(exc))
