import re

import pandas as pd
import requests

from utils.generic import flatten_list_to_list

ACC_COL = 11
DESC_COL = 12
INTERPRO_API = "https://www.ebi.ac.uk/interpro/api/entry/interpro/"


def extract_accessions_from_transcript(tran, lookup=False, acc_product=None):
    try:
        info = tran.attributes.get("info")[0]
    except TypeError:
        return
    for attr in info.split("\n"):
        acc = re.search(r'IPR\d+', attr).group()
        desc = attr.split("description:")[1].strip()
        if acc_product is not None and acc in flatten_list_to_list(acc_product.values()):
            continue
        if lookup:
            resp = requests.get(INTERPRO_API + acc, timeout=30).json()
            try:
                yield acc, desc, resp["metadata"]["type"]
            except KeyError:
                pass
        else:
            yield acc, desc


def extract_accessions_from_tsv(tsv_file, lookup=False, acc_product=None):
    data = pd.read_csv(tsv_file, sep="\t")
    for acc, desc in data[data.iloc[:, ACC_COL].str.startswith("IPR")].iloc[:, [ACC_COL, DESC_COL]].to_records(index=False):
        if acc_product is not None and acc in flatten_list_to_list(acc_product.values()):
            continue
        if lookup:
            resp = requests.get(INTERPRO_API + acc, timeout=30).json()
            try:
                yield acc, desc, resp["metadata"]["type"]
            except KeyError:
                pass
        else:
            yield acc, desc

