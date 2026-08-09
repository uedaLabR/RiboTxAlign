from __future__ import annotations

from pathlib import Path
from typing import TextIO

import pysam

from db_construct.fasta_utils import (
    write_fasta_record,
)


DECOY_FASTA_SUFFIXES = (
    ".fa",
    ".fasta",
    ".fna",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
)


def is_decoy_fasta(
    fasta_path: Path,
) -> bool:
    """Return whether a file has a supported FASTA suffix."""

    file_name = fasta_path.name.lower()

    return any(
        file_name.endswith(suffix)
        for suffix in DECOY_FASTA_SUFFIXES
    )


def get_decoy_fasta_paths(
    decoy_directory: str | Path | None,
) -> list[Path]:
    """Return sorted FASTA files found in the decoy directory."""

    if decoy_directory is None:
        return []

    directory = Path(
        decoy_directory
    )

    if not directory.exists():
        raise FileNotFoundError(
            "Decoy directory was not found: "
            f"{directory}"
        )

    if not directory.is_dir():
        raise NotADirectoryError(
            "Decoy path is not a directory: "
            f"{directory}"
        )

    return sorted(
        path
        for path in directory.iterdir()
        if (
            path.is_file()
            and is_decoy_fasta(
                path
            )
        )
    )


def create_decoy_reference_id(
    original_id: str,
) -> str:
    """Create a normalized reference ID with a DECOY_ prefix."""

    original_id = original_id.strip()

    if not original_id:
        raise ValueError(
            "A decoy FASTA record has an empty ID."
        )

    normalized_id = (
        original_id
        .replace(" ", "_")
        .replace("\t", "_")
    )

    if normalized_id.startswith(
        "DECOY_"
    ):
        return normalized_id

    return f"DECOY_{normalized_id}"


def write_decoy_fasta_records(
    output_file: TextIO,
    decoy_directory: str | Path | None,
    seen_reference_ids: set[str],
) -> tuple[int, int]:
    """
    Append decoy FASTA records to an open reference FASTA file.

    Decoy sequences are written in the supplied orientation only.
    Reverse-complement records are not generated.

    Returns
    -------
    tuple[int, int]
        Number of written decoy records and total decoy bases.
    """

    written_decoy_count = 0
    total_decoy_bases = 0

    for decoy_fasta_path in get_decoy_fasta_paths(
        decoy_directory
    ):
        with pysam.FastxFile(
            str(decoy_fasta_path)
        ) as decoy_fasta:

            for record in decoy_fasta:
                reference_id = (
                    create_decoy_reference_id(
                        record.name
                    )
                )

                if reference_id in seen_reference_ids:
                    raise ValueError(
                        "Duplicate reference ID in FASTA output: "
                        f"{reference_id}; "
                        f"source={decoy_fasta_path}"
                    )

                sequence = record.sequence.upper()

                if not sequence:
                    raise ValueError(
                        "Empty decoy FASTA sequence: "
                        f"{reference_id}; "
                        f"source={decoy_fasta_path}"
                    )

                write_fasta_record(
                    output_file=output_file,
                    reference_id=reference_id,
                    sequence=sequence,
                )

                seen_reference_ids.add(
                    reference_id
                )

                written_decoy_count += 1
                total_decoy_bases += len(
                    sequence
                )

    return (
        written_decoy_count,
        total_decoy_bases,
    )
