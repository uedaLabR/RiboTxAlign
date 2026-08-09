from depth_output.transcript_depth import (
    create_transcript_depth_tracks,
)


READS_INFO_PARQUET = (
    "/mnt/share/project/riboML/sample_ref/"
    "transcript_reads_info.parquet"
)

REFERENCE_METADATA_TSV = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.metadata.tsv"
)

REFERENCE_GTF = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.gtf"
)

TRANSCRIPT_FASTA = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.fa"
)

COORDINATE_TABLE = (
    "/mnt/share/project/riboML/sample_ref/"
    "HEK293T_RiboTxAlign.coordinates.tsv"
)

GENOME_FAI = (
    "/mnt/ssdnas07/pipeline/rna_modscore_v2/source/"
    "genome/hg38.fa.fai"
)

OUTPUT_PREFIX = (
    "/mnt/share/project/riboML/sample_ref/depth/"
    "HEK293T_RiboTxAlign"
)


def main() -> None:
    """Create transcript and genome Ribo-seq depth tracks."""

    output_paths = create_transcript_depth_tracks(
        reads_info_parquet=READS_INFO_PARQUET,
        reference_metadata_tsv=REFERENCE_METADATA_TSV,
        reference_gtf=REFERENCE_GTF,
        transcript_fasta=TRANSCRIPT_FASTA,
        coordinate_table=COORDINATE_TABLE,
        genome_fai=GENOME_FAI,
        output_prefix=OUTPUT_PREFIX,

        profile_minimum=-40,
        profile_maximum=20,

        minimum_offset=8,
        maximum_offset=20,
        fallback_offset=12,
        minimum_profile_weight=10.0,

        metagene_minimum=-50,
        metagene_maximum=100,

        batch_size=100000,

        create_bigwig=True,
        create_genome_bigwig=True,
    )

    print("Finished.")

    for output_name, output_path in output_paths.items():
        print(f"{output_name}: {output_path}")


if __name__ == "__main__":
    main()