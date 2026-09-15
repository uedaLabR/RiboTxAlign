from pathlib import Path
import pandas as pd

BASE = Path(
    "/mnt/share/project/riboML/validation_expression/similation/HEK293T"
)

INPUT = BASE / "HEK293T.simulation_abundance.sf"
OUTPUT = BASE / "HEK293T.art_modern_coverage.tsv"


# BASE = Path(
#     "/mnt/share/project/riboML/validation_expression/similation/mouse_brain"
# )
#
# INPUT = BASE / "mouse_brain.simulation_abundance.sf"
# OUTPUT = BASE / "mouse_brain.art_modern_coverage.tsv"

TARGET_PAIRS = 50_000_000
READ_LENGTH = 150


df = pd.read_csv(INPUT, sep="\t")

scale = (
    TARGET_PAIRS
    * 2
    * READ_LENGTH
    / (df["TPM"] * df["Length"]).sum()
)

df["coverage"] = df["TPM"] * scale

df[
    ["Name", "coverage"]
].to_csv(
    OUTPUT,
    sep="\t",
    index=False,
    header=False,
)

print("Scale:", scale)
print("Transcripts:", len(df))
print("Mean coverage:", df["coverage"].mean())
print("Median coverage:", df["coverage"].median())
print("Max coverage:", df["coverage"].max())
print("Expected read pairs:", TARGET_PAIRS)
print("Output:", OUTPUT)