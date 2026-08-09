from __future__ import annotations

from typing import TextIO


def write_fasta_record(
    output_file: TextIO,
    reference_id: str,
    sequence: str,
    line_width: int = 80,
) -> None:
    """Write one FASTA record."""

    output_file.write(
        f">{reference_id}\n"
    )

    for start in range(
        0,
        len(sequence),
        line_width,
    ):
        output_file.write(
            sequence[
                start:start + line_width
            ]
            + "\n"
        )
