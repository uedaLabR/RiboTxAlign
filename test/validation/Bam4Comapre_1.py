"""Bounded-memory publication counts. Edit constants; no CLI required.

Only all_mapper_reads summaries are produced. Detailed pairwise tables and
common_to_all are intentionally left to the original Bam4Comapre.py.
Requires pysam and GNU sort on the benchmark server.
"""
from collections import Counter
from itertools import groupby
from pathlib import Path
import csv
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
import time

BASE = Path('/mnt/share/project/riboML/validation')
OUTPUT_DIR = BASE / 'bam_comparison_four_mappers_counts'
# Prefer node-local scratch, e.g. Path('/scratch/ueda/bam_counts').
# None uses TMPDIR (or the operating system temporary directory).
SCRATCH_DIR = None
SORT_MEMORY = '256M'
SORT_THREADS = 2
LOG_EVERY = 2_000_000
DATASETS = {
    name: dict(samples=samples, star_genome_dir=BASE/'star_validation',
               star_transcript_dir=BASE/f'star_transcript_validation_{name}',
               ribotx_longread_stringtie_dir=BASE/'ribotx_validation',
               ribotx_shortread_stringtie_dir=BASE/'ribotx_validation_rnaseq')
    for name, samples in [('HEK293T', ['SRR26942051', 'SRR26942052', 'SRR26942053', 'SRR26942054']),
                          ('mouse_brain', ['SRR2163083', 'SRR2163093', 'SRR2163103'])]
}
MAPPERS = [('STAR_genome', 'star_genome_bam', False),
           ('STAR_transcript', 'star_transcript_bam', False),
           ('RiboTxAlign_longread_stringtie', 'ribotx_longread_bam', True),
           ('RiboTxAlign_shortread_stringtie', 'ribotx_shortread_bam', True)]

def normalize_read_id(
    read_id: str,
) -> str:
    """
    Normalize read IDs across mapping methods.

    Examples
    --------

    SRR26942051.34957075_CGACTTAGC

    SRR26942051.34957075|umi=CCGTCGCCC|family_size=1

    SRR26942051.34957075

    become

    SRR26942051.34957075
    """

    read_id = (
        read_id
        .split()[0]
    )

    match = re.search(
        r"(SRR\d+\.\d+)",
        read_id,
    )

    if match is not None:
        return match.group(1)

    return read_id

def remove_bam_suffix(
    bam_path: Path,
) -> str:
    """
    Obtain sample name from a BAM filename.

    The longest suffixes must be tested first.
    """

    name = bam_path.name

    suffixes = [
        ".Aligned.sortedByCoord.out.dedup.bam",
        ".Aligned.sortedByCoord.out.deduplicated.bam",
        ".deduplicated.bam",
        ".dedup.bam",
        ".Aligned.sortedByCoord.out.bam",
        ".genome.sorted.bam",
        ".sorted.bam",
        ".bam",
    ]

    for suffix in suffixes:

        if name.endswith(
            suffix
        ):
            return name[
                :-len(suffix)
            ]

    return bam_path.stem

def choose_single_bam(
    candidates,
    description: str,
) -> Path:

    candidates = sorted(
        set(candidates)
    )

    if not candidates:

        raise FileNotFoundError(
            f"{description}: BAM not found"
        )

    if len(candidates) == 1:
        return candidates[0]

    # --------------------------------------------------------
    # Prefer UMI-deduplicated BAM
    # --------------------------------------------------------

    dedup_candidates = [
        path
        for path in candidates
        if (
            path.name.endswith(
                ".dedup.bam"
            )
            or path.name.endswith(
                ".deduplicated.bam"
            )
        )
    ]

    if len(dedup_candidates) == 1:
        return dedup_candidates[0]

    if len(dedup_candidates) > 1:

        text = "\n".join(
            f"    {path}"
            for path in dedup_candidates
        )

        raise RuntimeError(
            f"{description}: "
            f"multiple deduplicated BAM files found:\n"
            f"{text}"
        )

    # --------------------------------------------------------
    # Prefer STAR coordinate-sorted BAM
    # --------------------------------------------------------

    coordinate_sorted = [
        path
        for path in candidates
        if (
            "Aligned.sortedByCoord.out.bam"
            in path.name
        )
    ]

    if len(coordinate_sorted) == 1:
        return coordinate_sorted[0]

    # --------------------------------------------------------
    # Prefer RiboTxAlign genome.sorted.bam
    # --------------------------------------------------------

    genome_sorted = [
        path
        for path in candidates
        if path.name.endswith(
            ".genome.sorted.bam"
        )
    ]

    if len(genome_sorted) == 1:
        return genome_sorted[0]

    text = "\n".join(
        f"    {path}"
        for path in candidates
    )

    raise RuntimeError(
        f"{description}: "
        f"multiple BAM files found:\n"
        f"{text}"
    )

def find_star_genome_bam(
    directory: Path,
    srr_id: str,
) -> Path:
    """
    STAR genome BAMs are stored directly under star_validation.
    """

    candidates = [
        path
        for path in directory.glob(
            "*.bam"
        )
        if (
            path.is_file()
            and srr_id in path.name
        )
    ]

    return choose_single_bam(
        candidates,
        f"STAR genome {srr_id}",
    )

def find_star_transcript_bam(
    directory: Path,
    srr_id: str,
) -> Path:
    """
    STAR transcript BAMs are inside sample directories.

    The UMI-deduplicated BAM is explicitly preferred:

        SRRxxxx.Aligned.sortedByCoord.out.dedup.bam

    Example:

        star_transcript_validation_HEK293T/
            SRR26942051/
                SRR26942051.Aligned.sortedByCoord.out.bam
                SRR26942051.Aligned.sortedByCoord.out.dedup.bam
    """

    sample_dir = (
        directory
        / srr_id
    )

    # --------------------------------------------------------
    # First choice: exact deduplicated filename
    # --------------------------------------------------------

    dedup_bam = (
        sample_dir
        / (
            f"{srr_id}."
            "Aligned.sortedByCoord.out.dedup.bam"
        )
    )

    if dedup_bam.is_file():
        return dedup_bam

    # --------------------------------------------------------
    # Second choice: any dedup BAM in sample directory
    # --------------------------------------------------------

    dedup_candidates = [
        path
        for path in sample_dir.glob(
            "*.dedup.bam"
        )
        if path.is_file()
    ]

    if len(dedup_candidates) == 1:
        return dedup_candidates[0]

    if len(dedup_candidates) > 1:
        return choose_single_bam(
            dedup_candidates,
            f"STAR transcript dedup {srr_id}",
        )

    # --------------------------------------------------------
    # Fallback: original STAR transcript BAM
    # --------------------------------------------------------

    original_bam = (
        sample_dir
        / (
            f"{srr_id}."
            "Aligned.sortedByCoord.out.bam"
        )
    )

    if original_bam.is_file():

        print(
            "WARNING: deduplicated STAR transcript BAM "
            f"not found for {srr_id}."
        )

        print(
            f"Using original BAM instead: {original_bam}"
        )

        return original_bam

    # --------------------------------------------------------
    # Last fallback: search recursively
    # --------------------------------------------------------

    candidates = []

    if sample_dir.exists():

        candidates.extend(
            sample_dir.rglob(
                "*.bam"
            )
        )

    if not candidates:

        for path in directory.iterdir():

            if not path.is_dir():
                continue

            if srr_id not in path.name:
                continue

            candidates.extend(
                path.rglob(
                    "*.bam"
                )
            )

    return choose_single_bam(
        candidates,
        f"STAR transcript {srr_id}",
    )

def find_ribotx_bam(
    directory: Path,
    sample_name: str,
    srr_id: str,
    description: str,
) -> Path:
    """
    Locate RiboTxAlign genome BAM.

    Expected structure:

        directory/
            sample_name/
                bam/
                    sample_name.genome.sorted.bam
    """

    expected = (
        directory
        / sample_name
        / "bam"
        / (
            f"{sample_name}."
            "genome.sorted.bam"
        )
    )

    if expected.is_file():
        return expected

    # Some directories may use SRR ID directly.
    expected_srr = (
        directory
        / srr_id
        / "bam"
        / (
            f"{srr_id}."
            "genome.sorted.bam"
        )
    )

    if expected_srr.is_file():
        return expected_srr

    candidates = [
        path
        for path in directory.rglob(
            "*.genome.sorted.bam"
        )
        if (
            srr_id
            in str(path)
        )
    ]

    return choose_single_bam(
        candidates,
        f"{description} {srr_id}",
    )

def find_sample_bams(
    dataset_name: str,
    config: dict,
):

    results = []

    for srr_id in config[
        "samples"
    ]:

        print()
        print(
            "----------------------------------------"
        )
        print(
            f"Searching BAMs: "
            f"{dataset_name} / {srr_id}"
        )
        print(
            "----------------------------------------"
        )

        # ----------------------------------------------------
        # STAR genome
        # ----------------------------------------------------

        star_genome_bam = (
            find_star_genome_bam(
                config[
                    "star_genome_dir"
                ],
                srr_id,
            )
        )

        sample_name = (
            remove_bam_suffix(
                star_genome_bam
            )
        )

        # ----------------------------------------------------
        # STAR transcript
        # ----------------------------------------------------

        star_transcript_bam = (
            find_star_transcript_bam(
                config[
                    "star_transcript_dir"
                ],
                srr_id,
            )
        )

        # ----------------------------------------------------
        # RiboTxAlign long-read StringTie
        # ----------------------------------------------------

        ribotx_longread_bam = (
            find_ribotx_bam(
                directory=(
                    config[
                        "ribotx_longread_stringtie_dir"
                    ]
                ),
                sample_name=(
                    sample_name
                ),
                srr_id=(
                    srr_id
                ),
                description=(
                    "RiboTxAlign long-read StringTie"
                ),
            )
        )

        # ----------------------------------------------------
        # RiboTxAlign short-read StringTie
        # ----------------------------------------------------

        ribotx_shortread_bam = (
            find_ribotx_bam(
                directory=(
                    config[
                        "ribotx_shortread_stringtie_dir"
                    ]
                ),
                sample_name=(
                    sample_name
                ),
                srr_id=(
                    srr_id
                ),
                description=(
                    "RiboTxAlign short-read StringTie"
                ),
            )
        )

        print(
            f"Sample name                  : "
            f"{sample_name}"
        )

        print(
            f"STAR genome                  : "
            f"{star_genome_bam}"
        )

        print(
            f"STAR transcript              : "
            f"{star_transcript_bam}"
        )

        print(
            f"RiboTx long-read StringTie   : "
            f"{ribotx_longread_bam}"
        )

        print(
            f"RiboTx short-read StringTie  : "
            f"{ribotx_shortread_bam}"
        )

        results.append(
            {
                "dataset_name": (
                    dataset_name
                ),
                "sample_name": (
                    sample_name
                ),
                "srr_id": (
                    srr_id
                ),
                "star_genome_bam": (
                    star_genome_bam
                ),
                "star_transcript_bam": (
                    star_transcript_bam
                ),
                "ribotx_longread_bam": (
                    ribotx_longread_bam
                ),
                "ribotx_shortread_bam": (
                    ribotx_shortread_bam
                ),
            }
        )

    return results

def compact_records(records, ribotx):
    """Keep only counting fields; retain original BAM order within each ID."""
    start = time.monotonic()
    for n, record in enumerate(records, 1):
        if n % LOG_EVERY == 0:
            logging.info('BAM scanned %s records; elapsed %.1fs', f'{n:,}', time.monotonic()-start)
        if record.is_supplementary or (ribotx and record.is_secondary):
            continue
        read_id = normalize_read_id(record.query_name)
        mapped = not record.is_unmapped
        nh = int(record.get_tag('NH')) if mapped and record.has_tag('NH') else -1
        splice = int(mapped and bool(record.cigarstring) and 'N' in record.cigarstring)
        assignment = str(record.get_tag('XA')) if ribotx and mapped and record.has_tag('XA') else 'unknown'
        # Rescue categories alone affect the rescue summary.
        code = {'multimap_assigned': 1, 'genome_rescued_unique': 2,
                'genome_splice_rescued_unique': 3}.get(assignment, 0)
        yield (read_id, int(mapped), int(not record.is_secondary), nh, splice, code)

def collapse_read(rows, ribotx):
    """Exact legacy rules: STAR first mapped primary; RiboTx last primary."""
    chosen = None
    mapped_count = 0
    spliced = 0
    for row in rows:
        _, mapped, primary, nh, splice, code = row
        if ribotx:
            chosen = row
        elif mapped:
            mapped_count += 1
            spliced |= splice
            if chosen is None or (primary and not chosen[2]):
                chosen = row
    if chosen is None or not chosen[1]:
        return 'unmapped', 0, 0, 0
    nh = chosen[3]
    if nh == -1:
        nh = 1 if ribotx else mapped_count
    multi = chosen[5] == 1 if ribotx else nh != 1
    return ('multi_assigned' if ribotx else 'multi') if multi else 'unique', nh, chosen[4] if ribotx else spliced, chosen[5]

def aggregate(rows, ribotx):
    counts = Counter()
    nh_counts = Counter()
    start = time.monotonic()
    for n, (_, group) in enumerate(groupby(rows, key=lambda row: row[0]), 1):
        status, nh, splice, code = collapse_read(group, ribotx)
        counts['total_reads'] += 1
        counts[status] += 1
        counts['splice_reads'] += splice
        if status in ('multi', 'multi_assigned'):
            nh_counts[nh] += 1
        if status != 'unmapped':
            counts['contiguous_rescue'] += code == 2
            counts['splice_rescue'] += code == 3
        if n % LOG_EVERY == 0:
            logging.info('Grouped %s reads; elapsed %.1fs', f'{n:,}', time.monotonic()-start)
    return dict(counts), dict(nh_counts)

def sorted_rows(path):
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            fields = line.rstrip('\n').split('\t')
            yield (fields[0], *(int(v) for v in fields[1:]))

def scan_bam(path, ribotx):
    import pysam
    sort = shutil.which('sort')
    if not sort:
        raise RuntimeError('GNU sort is required on the server.')
    version = subprocess.run([sort, '--version'], capture_output=True, text=True, check=True)
    if 'GNU' not in version.stdout:
        raise RuntimeError('Expected GNU sort, not the Windows sort command.')
    if SCRATCH_DIR is not None:
        Path(SCRATCH_DIR).mkdir(parents=True, exist_ok=True)
    # GNU sort spills to disk instead of holding tens of millions of read IDs.
    with tempfile.TemporaryDirectory(prefix='ribotx_counts_', dir=SCRATCH_DIR) as temp:
        work = Path(temp)
        target = work/'reads.sorted.tsv'
        logging.info('Scratch: %s; free %.1f GiB; sort buffer %s', work,
                     shutil.disk_usage(work).free/2**30, SORT_MEMORY)
        command = [sort, '-s', '-t', '\t', '-k1,1', '-S', SORT_MEMORY,
                   f'--parallel={SORT_THREADS}', '-T', str(work)]
        started = time.monotonic()
        with target.open('wb') as output, (work/'sort.stderr').open('w+b') as error:
            proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=output, stderr=error,
                                    env={**os.environ, 'LC_ALL': 'C'}, text=True, encoding='utf-8')
            try:
                with pysam.AlignmentFile(path, 'rb') as bam:
                    for row in compact_records(bam.fetch(until_eof=True), ribotx):
                        proc.stdin.write('\t'.join(map(str, row))+'\n')
                proc.stdin.close()
                logging.info('BAM scan finished; waiting for external sort')
                result = proc.wait()
                if result:
                    error.seek(0)
                    raise RuntimeError(f'GNU sort failed ({result}): {error.read().decode(errors="replace")}')
            except BaseException:
                if proc.poll() is None:
                    proc.terminate()
                proc.wait()
                raise
        logging.info('BAM scan + sort: %.1fs', time.monotonic()-started)
        return aggregate(sorted_rows(target), ribotx)

def signature(path, ribotx):
    stat = path.stat()
    return dict(path=str(path.resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                ribotx=ribotx, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

def atomic_json(path, data):
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(data, indent=2), encoding='utf-8')
    temp.replace(path)

def cached_counts(path, ribotx, cache):
    sig = signature(path, ribotx)
    if cache.exists():
        try:
            saved = json.loads(cache.read_text(encoding='utf-8'))
            if saved['signature'] == sig:
                logging.info('Reusing completed BAM count: %s', path)
                return saved['counts'], {int(k): v for k, v in saved['nh'].items()}
        except (ValueError, KeyError, TypeError):
            logging.warning('Ignoring invalid cache %s', cache)
    started = time.monotonic()
    logging.info('START %s', path)
    counts, nh = scan_bam(path, ribotx)
    if signature(path, ribotx) != sig:
        raise RuntimeError(f'BAM changed during counting: {path}')
    atomic_json(cache, dict(signature=sig, counts=counts, nh=nh))
    logging.info('DONE %s: %s reads, %.1fs', path.name, counts.get('total_reads', 0), time.monotonic()-started)
    return counts, nh

def percent(a, b):
    return a/b*100 if b else 0.0

def tables(mapper, c, nh):
    total = c.get('total_reads', 0)
    unique = c.get('unique', 0)
    multi = c.get('multi', 0)+c.get('multi_assigned', 0)
    unmapped = c.get('unmapped', 0)
    mapped = unique+multi
    assert total == mapped+unmapped
    splice = c.get('splice_reads', 0) if mapper != 'STAR_transcript' else None
    row = dict(scope='all_mapper_reads', mapper=mapper, total_reads=total,
               mapped_reads=mapped, unmapped_reads=unmapped, unique_reads=unique, multi_reads=multi)
    for label, count in [('mapped', mapped), ('unmapped', unmapped), ('unique', unique), ('multi', multi)]:
        row[label+'_percent_of_total'] = percent(count, total)
    row.update(unique_percent_of_mapped=percent(unique, mapped), multi_percent_of_mapped=percent(multi, mapped),
               splice_reads=splice, splice_percent_of_mapped=None if splice is None else percent(splice, mapped),
               splice_percent_of_total=None if splice is None else percent(splice, total))
    hits = [dict(scope='all_mapper_reads', mapper=mapper, nh=n, read_count=v,
                 percent_of_multi=percent(v, multi), percent_of_total=percent(v, total)) for n,v in sorted(nh.items())]
    rescue = None
    if mapper.startswith('RiboTxAlign'):
        contiguous, spliced = c.get('contiguous_rescue', 0), c.get('splice_rescue', 0)
        rescued = contiguous+spliced
        input_count = rescued+unmapped
        rescue = dict(mapper=mapper, total_ribotx_reads=total, rescue_input_reads=input_count,
                      contiguous_rescued_unique_reads=contiguous, splice_rescued_unique_reads=spliced,
                      rescued_unique_reads=rescued, remaining_unmapped_reads=unmapped,
                      rescue_percent_of_rescue_input=percent(rescued,input_count),
                      contiguous_rescue_percent_of_rescue_input=percent(contiguous,input_count),
                      splice_rescue_percent_of_rescue_input=percent(spliced,input_count),
                      rescue_percent_of_all_ribotx_reads=percent(rescued,total))
    return row, hits, rescue

def write_tsv(path, rows, fields):
    temp = path.with_suffix('.tsv.tmp')
    with temp.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.StreamHandler(), logging.FileHandler(OUTPUT_DIR/'counts.log')])
    started = time.monotonic()
    failures = []
    for dataset, config in DATASETS.items():
        for sample in config['samples']:
            try:
                paths = find_sample_bams(dataset, {**config, 'samples': [sample]})[0]
                directory = OUTPUT_DIR/dataset/sample
                directory.mkdir(parents=True, exist_ok=True)
                summary, hits, rescues = [], [], []
                for mapper, key, ribotx in MAPPERS:
                    counts, nh = cached_counts(paths[key], ribotx, directory/(mapper+'.counts.json'))
                    row, rows, rescue = tables(mapper, counts, nh)
                    summary.append(row); hits.extend(rows)
                    if rescue is not None: rescues.append(rescue)
                write_tsv(directory/'mapping_summary.tsv', summary, list(summary[0]))
                write_tsv(directory/'multi_reads_by_nh.tsv', hits,
                          ['scope','mapper','nh','read_count','percent_of_multi','percent_of_total'])
                write_tsv(directory/'ribotx_internal_rescue_summary.tsv', rescues, list(rescues[0]))
                logging.info('SAMPLE COMPLETE %s/%s', dataset, sample)
            except Exception:
                failures.append(f'{dataset}/{sample}')
                logging.exception('Sample failed: %s/%s; continuing with other samples', dataset, sample)
    logging.info('Total elapsed %.1fs; failed samples: %s', time.monotonic()-started, failures)
    if failures:
        raise SystemExit(1)

if __name__ == '__main__':
    main()
