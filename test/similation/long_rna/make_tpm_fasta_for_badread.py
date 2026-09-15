from pathlib import Path
import pandas as pd
import re

# BASE = Path(
#     "/mnt/share/project/riboML/"
#     "validation_expression/similation/mouse_brain"
# )
#
# FASTA = BASE / "mouse_brain.simulation_transcripts.fa"
# ABUNDANCE = BASE / "mouse_brain.simulation_abundance.sf"
# OUTPUT = BASE / "mouse_brain.badread_transcripts.fa"


BASE = Path(
    "/mnt/share/project/riboML/"
    "validation_expression/similation/HEK293T"
)

STRINGTIE_GTF = Path(
    "/mnt/share/project/riboML/"
    "validation/HEK293T/longread/mapping/"
    "HEK293T.longread.stringtie.gtf"
)

SIM_FASTA = BASE / "HEK293T.simulation_transcripts.fa"
INFO_TSV = BASE / "HEK293T.simulation_transcript_info.tsv"
OUTPUT = BASE / "HEK293T.badread_transcripts.fa"

BASE = Path(
    "/mnt/share/project/riboML/"
    "validation_expression/similation/HEK293T"
)

STRINGTIE_GTF = Path(
    "/mnt/share/project/riboML/"
    "validation/mouse_brain/longread/mapping/"
    "mouse_brain.longread.stringtie.gtf"
)

BASE = Path(
    "/mnt/share/project/riboML/"
    "validation_expression/similation/mouse_brain"
)
#
# FASTA = BASE / "mouse_brain.simulation_transcripts.fa"
# ABUNDANCE = BASE / "mouse_brain.simulation_abundance.sf"
# OUTPUT = BASE / "mouse_brain.badread_transcripts.fa"

SIM_FASTA = BASE / "mouse_brain.simulation_transcripts.fa"
INFO_TSV = BASE / "mouse_brain.simulation_transcript_info.tsv"
OUTPUT = BASE / "mouse_brain.badread_transcripts.fa"

MIN_DEPTH = 3.0


def read_stringtie_cov(gtf_path):
    cov = {}

    with open(gtf_path) as f:
        for line in f:

            if "\ttranscript\t" not in line:
                continue

            tid_match = re.search(
                r'transcript_id "([^"]+)"',
                line,
            )

            cov_match = re.search(
                r'(?:^|;\s)cov "([0-9.eE+-]+)"',
                line,
            )

            if not tid_match or not cov_match:
                continue

            tid = tid_match.group(1)
            value = float(
                cov_match.group(1)
            )

            cov[tid] = value

    return cov


def main():

    info = pd.read_csv(
        INFO_TSV,
        sep="\t",
    )

    stringtie_cov = read_stringtie_cov(
        STRINGTIE_GTF
    )

    # Map simulation ENST ID to original StringTie transcript ID.
    sim_to_stringtie = dict(
        zip(
            info["simulation_transcript_id"],
            info["stringtie_transcript_id"],
        )
    )

    depths = {}

    for sim_id, stringtie_id in (
        sim_to_stringtie.items()
    ):

        cov = stringtie_cov.get(
            stringtie_id
        )

        if cov is None:
            continue

        depths[sim_id] = max(
            cov,
            MIN_DEPTH,
        )

    with open(SIM_FASTA) as src, open(OUTPUT, "w") as out:

        for line in src:

            if line.startswith(">"):

                tid = (
                    line[1:]
                    .strip()
                    .split()[0]
                    .split(".")[0]
                )

                if tid not in depths:
                    raise ValueError(
                        f"No StringTie coverage for {tid}"
                    )

                out.write(
                    f">{tid} depth={depths[tid]:.8f}\n"
                )

            else:
                out.write(line)

    values = pd.Series(
        list(
            depths.values()
        )
    )

    print("Transcripts:", len(values))
    print("Min depth:", values.min())
    print("Median depth:", values.median())
    print("Mean depth:", values.mean())
    print("q90:", values.quantile(0.90))
    print("q95:", values.quantile(0.95))
    print("q99:", values.quantile(0.99))
    print("Max depth:", values.max())
    print()
    print("Output:", OUTPUT)


if __name__ == "__main__":
    main()

