#!/usr/bin/env python3

from collections import defaultdict
import pysam


BAM_FILE = (
    "/mnt/share/project/riboML/validation_expression/"
    "riboMap/HEK293T/mapping/sorted.bam"
)

GFF3_FILE = (
    "/mnt/share/project/riboML/validation/reference/"
    "gencode.v44.annotation.gff3"
)

CDS_FILE = (
    "/mnt/share/project/riboML/validation_expression/"
    "riboMap/HEK293T/cds_range.txt"
)

OUTPUT_FILE = (
    "/mnt/share/project/riboML/validation_expression/"
    "deepShape/HEK293T/HEK293T.absolute"
)

# BAM_FILE = (
#     "/mnt/share/project/riboML/validation_expression/"
#     "riboMap/mouse_brain/mapping/sorted.bam"
# )
#
# GFF3_FILE = (
#     "/mnt/share/project/riboML/validation/reference/"
#     "gencode.vM25.annotation.gff3"
# )
#
# CDS_FILE = (
#     "/mnt/share/project/riboML/validation_expression/"
#     "riboMap/mouse_brain/cds_range.txt"
# )
#
# OUTPUT_FILE = (
#     "/mnt/share/project/riboML/validation_expression/"
#     "deepShape/mouse_brain/mouse_brain.absolute"
# )

OFFSETS = {
    25: 12,
    26: 12,
    27: 12,
    28: 12,
    29: 12,
    30: 12,
    31: 13,
    32: 13,
    33: 13,
    34: 14,
    35: 14,
    36: 14,
}


def clean_id(value):
    if value is None:
        return None

    if ":" in value:
        value = value.split(":", 1)[1]

    return value.split(".")[0]


def parse_gff3_attributes(text):
    attributes = {}

    for item in text.strip().split(";"):
        if "=" not in item:
            continue

        key, value = item.split("=", 1)
        attributes[key] = value

    return attributes


def load_cds_ranges(path):
    cds_ranges = {}

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")

            transcript_id = clean_id(fields[0])
            cds_start = int(fields[1])
            cds_stop = int(fields[2])

            cds_ranges[transcript_id] = (
                cds_start,
                cds_stop,
            )

    return cds_ranges


def load_cds_exons(gff3_file):
    cds_exons = defaultdict(list)
    transcript_gene = {}

    with open(gff3_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if len(fields) != 9:
                continue

            chrom = fields[0]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            attrs = parse_gff3_attributes(
                fields[8]
            )

            if feature in ("transcript", "mRNA"):

                transcript_id = clean_id(
                    attrs.get("ID")
                )

                gene_id = (
                    attrs.get("gene_id")
                    or attrs.get("geneID")
                    or attrs.get("Parent")
                )

                if transcript_id is not None:
                    transcript_gene[
                        transcript_id
                    ] = clean_id(gene_id)

            elif feature == "CDS":

                parents = attrs.get("Parent")

                if parents is None:
                    continue

                for parent in parents.split(","):

                    transcript_id = clean_id(
                        parent
                    )

                    cds_exons[
                        transcript_id
                    ].append(
                        (
                            chrom,
                            strand,
                            start,
                            end,
                        )
                    )

    for transcript_id, blocks in cds_exons.items():

        strand = blocks[0][1]

        if strand == "+":
            blocks.sort(
                key=lambda x: x[2]
            )
        else:
            blocks.sort(
                key=lambda x: x[2],
                reverse=True,
            )

    return cds_exons, transcript_gene


def cds_position_to_genome(
    blocks,
    cds_position,
):
    remaining = cds_position

    for chrom, strand, start, end in blocks:

        block_length = end - start + 1

        if remaining > block_length:
            remaining -= block_length
            continue

        if strand == "+":
            genome_position = (
                start + remaining - 1
            )
        else:
            genome_position = (
                end - remaining + 1
            )

        return (
            chrom,
            strand,
            genome_position,
        )

    return None


def main():

    print("Loading CDS ranges...")
    cds_ranges = load_cds_ranges(
        CDS_FILE
    )

    print(
        f"CDS transcripts: "
        f"{len(cds_ranges):,}"
    )

    print("Loading GFF3...")
    cds_exons, transcript_gene = (
        load_cds_exons(
            GFF3_FILE
        )
    )

    print(
        f"GFF3 CDS transcripts: "
        f"{len(cds_exons):,}"
    )

    print("Reading BAM...")

    bam = pysam.AlignmentFile(
        BAM_FILE,
        "rb",
    )

    total = 0
    written = 0
    unmapped = 0
    skipped_length = 0
    skipped_cds = 0
    skipped_position = 0

    with open(
        OUTPUT_FILE,
        "w",
    ) as out:

        for read in bam.fetch(
            until_eof=True
        ):

            total += 1

            if read.is_unmapped:
                unmapped += 1
                continue

            read_length = read.query_length

            if read_length not in OFFSETS:
                skipped_length += 1
                continue

            transcript_id = clean_id(
                bam.get_reference_name(
                    read.reference_id
                )
            )

            if (
                transcript_id
                not in cds_ranges
            ):
                skipped_cds += 1
                continue

            if (
                transcript_id
                not in cds_exons
            ):
                skipped_cds += 1
                continue

            cds_start, cds_stop = (
                cds_ranges[
                    transcript_id
                ]
            )

            # BAM reference_start is 0-based.
            # DeepShape original script uses
            # 1-based transcript coordinates.
            mapping_start = (
                read.reference_start + 1
            )

            ribosome_position = (
                mapping_start
                + OFFSETS[read_length]
            )

            if (
                ribosome_position
                < cds_start
            ):
                if (
                    ribosome_position
                    < cds_start - 3
                ):
                    continue

                ribosome_position = cds_start

            if (
                ribosome_position
                > cds_stop
            ):
                if (
                    ribosome_position
                    > cds_stop + 3
                ):
                    continue

                ribosome_position = cds_stop

            ribosome_position_in_cds = (
                ribosome_position
                - cds_start
                + 1
            )

            genome_result = (
                cds_position_to_genome(
                    cds_exons[
                        transcript_id
                    ],
                    ribosome_position_in_cds,
                )
            )

            if genome_result is None:
                skipped_position += 1
                continue

            (
                chrom,
                strand,
                genome_position,
            ) = genome_result

            gene_id = transcript_gene.get(
                transcript_id,
                transcript_id,
            )

            out.write(
                f"{read.query_name}\t"
                f"{chrom}\t"
                f"{strand}\t"
                f"{genome_position}\t"
                f"{ribosome_position_in_cds}\t"
                f"{transcript_id}\t"
                f"{gene_id}\n"
            )

            written += 1

            if total % 1_000_000 == 0:
                print(
                    f"Processed: {total:,} "
                    f"Written: {written:,}"
                )

    bam.close()

    print()
    print("Finished")
    print(f"Total alignments : {total:,}")
    print(f"Written          : {written:,}")
    print(f"Unmapped         : {unmapped:,}")
    print(
        f"Length filtered  : "
        f"{skipped_length:,}"
    )
    print(
        f"No CDS annotation: "
        f"{skipped_cds:,}"
    )
    print(
        f"Invalid position : "
        f"{skipped_position:,}"
    )
    print(f"Output           : {OUTPUT_FILE}")


if __name__ == "__main__":
    main()