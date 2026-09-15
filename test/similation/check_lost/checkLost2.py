import pandas as pd
import numpy as np

df = pd.read_csv(
    "/mnt/share/project/riboML/validation_expression/similation/HEK293T/plots_transcript/transcript_pairwise_values.RiboTxAlign_LongRead.tsv",
    sep="\t",
)

df = df[
    (df["Truth_tpm"] > 0)
    & (df["RiboTxAlign_LongRead_tpm"] > 0)
].copy()

df["error"] = (
    np.log10(df["RiboTxAlign_LongRead_tpm"] + 1)
    - np.log10(df["Truth_tpm"] + 1)
)

print(
    df.sort_values(
        "error",
        key=abs,
        ascending=False,
    ).head(20).to_string(index=False)
)