import csv
from collections import Counter, defaultdict
from pathlib import Path


def load_reference_metadata(metadata_tsv):
    """Load reference transcript definitions from metadata TSV."""

    reference_definitions = {}

    with Path(metadata_tsv).open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            reference_id = row["reference_id"]

            reference_definitions[reference_id] = {
                "reference_id": reference_id,
                "group_id": row["overlap_group_id"],
                "transcript_length": int(row["transcript_length"]),
                "cds_length": int(row["cds_length"]),
                "strand": row.get("strand", ""),
                "chromosome": row.get("chromosome", ""),
                **{key: (int(row[key]) if row.get(key) not in (None, "", ".") else None)
                   for key in ("cds_start", "cds_end", "cds_reference_start", "cds_reference_end")},
                "exons": [],
            }

    return reference_definitions


def load_reference_exons(reference_gtf, reference_definitions):
    """Add exon intervals from the existing reference GTF."""

    with Path(reference_gtf).open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            columns = line.rstrip("\n").split("\t")

            if len(columns) != 9 or columns[2] != "exon":
                continue

            attributes = {}

            for item in columns[8].split(";"):
                item = item.strip()

                if not item or " " not in item:
                    continue

                key, value = item.split(" ", 1)
                attributes[key] = value.strip().strip('"')

            transcript_id = attributes.get("transcript_id")

            if transcript_id not in reference_definitions:
                continue

            reference_definitions[transcript_id]["strand"] = columns[6]
            reference_definitions[transcript_id]["chromosome"] = columns[0]

            reference_definitions[transcript_id]["exons"].append(
                (
                    int(columns[3]) - 1,
                    int(columns[4]),
                )
            )

    for definition in reference_definitions.values():
        definition["exons"].sort()

    return reference_definitions


def load_reference_definitions(
    metadata_tsv,
    reference_gtf,
):
    """Load reference metadata and exon definitions."""

    definitions = load_reference_metadata(metadata_tsv)

    return load_reference_exons(
        reference_gtf=reference_gtf,
        reference_definitions=definitions,
    )


def merge_intervals(intervals):
    """Merge overlapping half-open intervals."""

    if not intervals:
        return []

    merged = [list(interval) for interval in sorted(intervals)]
    output = [merged[0]]

    for start, end in merged[1:]:
        if start <= output[-1][1]:
            output[-1][1] = max(output[-1][1], end)
        else:
            output.append([start, end])

    return output


def calculate_unique_exon_length(own_exons, other_exons):
    """Return exon length not covered by other transcripts."""

    other_exons = merge_intervals(other_exons)
    unique_length = 0

    for own_start, own_end in own_exons:
        overlap_length = 0

        for other_start, other_end in other_exons:
            if other_end <= own_start:
                continue
            if other_start >= own_end:
                break

            overlap_length += max(
                0,
                min(own_end, other_end)
                - max(own_start, other_start),
            )

        unique_length += own_end - own_start - overlap_length

    return unique_length


def get_junctions(exons):
    """Return splice junctions from ordered exon intervals."""

    exons = sorted(exons)

    return {
        (exons[index][1], exons[index + 1][0])
        for index in range(len(exons) - 1)
    }


def determine_unique_transcript_lengths(
    transcript_ids,
    reference_definitions,
    unique_junction_length=50,
):
    """Calculate approximate transcript-unique effective lengths."""

    transcript_ids = list(dict.fromkeys(transcript_ids))
    junctions = {
        transcript_id: get_junctions(
            reference_definitions[transcript_id]["exons"]
        )
        for transcript_id in transcript_ids
    }

    # Count each junction once per transcript instead of rebuilding the union
    # of every other transcript's junctions for each candidate.
    junction_counts = Counter(junction for values in junctions.values() for junction in values)
    unique_exon_lengths = dict.fromkeys(transcript_ids, 0)

    # Sweep exon boundaries once for the whole group. A segment is unique
    # exactly when one transcript covers it. Track exon multiplicity as well:
    # the former per-exon calculation counted duplicate/overlapping own exons
    # separately, while coverage by ANY other transcript excluded the segment.
    events = defaultdict(list)
    for tid in transcript_ids:
        for start, end in reference_definitions[tid]["exons"]:
            if start < end:
                events[start].append((tid, 1))
                events[end].append((tid, -1))

    active = {}
    previous_position = 0
    for position in sorted(events):
        if len(active) == 1:
            tid, multiplicity = next(iter(active.items()))
            unique_exon_lengths[tid] += (position - previous_position) * multiplicity
        # Process every event at a shared boundary before the next segment;
        # touching half-open exons do not overlap.
        for tid, change in events[position]:
            count = active.get(tid, 0) + change
            if count:
                active[tid] = count
            else:
                active.pop(tid, None)
        previous_position = position

    return {
        tid: unique_exon_lengths[tid]
        + sum(junction_counts[junction] == 1 for junction in junctions[tid])
        * unique_junction_length
        for tid in transcript_ids
    }


def determine_transcript_lengths(
    transcript_ids,
    group_records,
    reference_definitions,
    minimum_mapped_reads=10,
    minimum_mapped_span=100,
    minimum_transcript_fraction=1.0 / 3.0,
):
    """Determine CDS-like lengths used after EM initialization."""

    records_by_transcript = defaultdict(list)

    for record in group_records:
        records_by_transcript[record["transcript_id"]].append(record)

    transcript_lengths = {}

    for transcript_id in transcript_ids:
        definition = reference_definitions[transcript_id]
        if "em_region_length" in definition:
            transcript_lengths[transcript_id] = definition["em_region_length"]
            continue
        annotated_cds_length = definition["cds_length"]

        if annotated_cds_length > 0:
            transcript_lengths[transcript_id] = annotated_cds_length
            continue

        transcript_records = records_by_transcript.get(transcript_id, [])
        distinct_read_count = len(
            {record["read_id"] for record in transcript_records}
        )
        transcript_length = max(
            int(definition["transcript_length"]),
            1,
        )

        if transcript_records:
            mapped_start = min(
                int(record["reference_start"])
                for record in transcript_records
            )
            mapped_end = max(
                int(record["reference_end"])
                for record in transcript_records
            )
            mapped_span = max(mapped_end - mapped_start, 1)
        else:
            mapped_span = 0

        use_mapped_span = (
            distinct_read_count > minimum_mapped_reads
            and mapped_span >= minimum_mapped_span
            and mapped_span >= transcript_length * minimum_transcript_fraction
        )

        transcript_lengths[transcript_id] = (
            mapped_span
            if use_mapped_span
            else transcript_length
        )

    return transcript_lengths
