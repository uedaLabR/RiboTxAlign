import re
import pandas as pd

GTF = (
    "/mnt/share/project/riboML/validation/"
    "HEK293T/longread/mapping/"
    "HEK293T.longread.stringtie.gtf"
)

INFO = (
    "/mnt/share/project/riboML/"
    "validation_expression/similation/HEK293T/"
    "HEK293T.simulation_transcript_info.tsv"
)

info = pd.read_csv(INFO, sep="\t")

selected = set(
    info["stringtie_transcript_id"]
)

records = []

with open(GTF) as f:
    for line in f:
        if "\ttranscript\t" not in line:
            continue

        tid = re.search(
            r'transcript_id "([^"]+)"',
            line,
        )

        cov = re.search(
            r'(?:^|;\s)cov "([0-9.eE+-]+)"',
            line,
        )

        if not tid or not cov:
            continue

        tid = tid.group(1)

        if tid in selected:
            records.append(
                float(cov.group(1))
            )

s = pd.Series(records)

print("n =", len(s))
print("min =", s.min())
print("median =", s.median())
print("mean =", s.mean())
print("max =", s.max())
print("q90 =", s.quantile(0.90))
print("q95 =", s.quantile(0.95))
print("q99 =", s.quantile(0.99))