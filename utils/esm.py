import re

def extract_esm_means(esm_stdouts_path):
    out = {}
    with open(esm_stdouts_path, "r") as f:
        for l in f.read().splitlines():
            match = re.search(r'structure for (?P<acc>[^\s]+).*pLDDT\s(?P<mean>\d+\.\d+)\,', l)
            if match:
                out[match["acc"]] = match["mean"]
    return out
