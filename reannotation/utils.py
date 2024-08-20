import re

import pandas as pd

ACC_COL = 11
DESC_COL = 12

def extract_accessions_from_transcript(tran):
    try:
        info = tran.attributes.get("info")[0]
    except TypeError:
        print(f"No InterPro accessions found for transcript {tran.id}")
        return
    print(info)
    for acc in info.split("\n"):
        yield re.search(r'IPR\d+', acc).group(), acc.split("description:")[1]


def extract_accessions_from_tsv(tsv_file):
    data = pd.read_csv(tsv_file, sep="\t")
    for acc, desc in data[data.iloc[:, ACC_COL].str.startswith("IPR")].iloc[:, [ACC_COL, DESC_COL]].to_records(index=False):
        yield acc, desc
