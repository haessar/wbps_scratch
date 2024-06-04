#! /usr/bin/env python3
import argparse
from collections import Counter
import json

from natsort import natsorted
import pandas as pd

from utils.generic import flatten_list_to_set


with open("data/transposable_elements/transposable_element_pfams.json") as f:
    TE_PFAMS = json.loads(f.read())

# pfamout column indices
TACC_IDX = 1
QNAM_IDX = 3
QLEN_IDX = 5


def find_biologically_interesting_genes(pfamout_path, unique_prots):
    encountered_te_pfams = []
    biologically_interesting_genes = set()
    with open(pfamout_path, newline="\n") as tsv_file:
        for line in tsv_file:
            for line in tsv_file:
                if not line.startswith("#"):
                    row = line.split()
                    # Is it greater than 150 amino acids?
                    if int(row[QLEN_IDX]) >= 150:
                        pfam = row[TACC_IDX].split('.')[0]
                        prot = row[QNAM_IDX]
                        if prot.startswith("g11954"):
                            print()
                        # Is it associated with transposable elements?
                        if pfam in TE_PFAMS:
                            encountered_te_pfams.append(pfam)
                        else:
                            if prot in unique_prots:
                                biologically_interesting_genes.add(prot.split(".")[0])
    print("Encountered transposable elements:")
    for k, v in Counter(encountered_te_pfams).most_common():
        print(k, v)
    print(f"TOTAL:\t{len(encountered_te_pfams)}")
    return biologically_interesting_genes


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('pfamout_path')
    parser.add_argument('--output', '-o', type=str, default=None)
    parser.add_argument("--combined", "-c", action="store_true")
    parser.add_argument("--odb", choices=['exists', 'missing'])
    args = parser.parse_args()
    
    orthogroups_path = "data/from_MARS/OrthoFinder/Results_Jun03_1/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
    cols = ["hog", "og", "clade", "wbps", "odb", "braker", "helixer"]

    df = pd.read_csv(orthogroups_path, delimiter="\t")
    df.columns = cols
    if args.odb == "exists":
        df = df[~df["odb"].isna()]
    elif args.odb == "missing":
        df = df[df["odb"].isna()]
    # Orthologues from Helixer (and Nematoda ODB10?) only
    helixer_only = flatten_list_to_set(df[~(df["helixer"].isna()) & (df["braker"].isna()) & (df["wbps"].isna())]["helixer"].str.split(",").tolist())
    # Orthologues from BRAKER (and Nematoda ODB10?) only
    braker_only = flatten_list_to_set(df[(df["helixer"].isna()) & ~(df["braker"].isna()) & (df["wbps"].isna())]["braker"].str.split(",").tolist())
    # Orthologues from BRAKER, Helixer (and Nematoda ODB10?) but missing from WBPS.
    novel_df = df[~(df["helixer"].isna()) & ~(df["braker"].isna()) & (df["wbps"].isna())]
    novel = set()
    novel.update(flatten_list_to_set(novel_df["braker"].str.split(",").tolist()))
    novel.update(flatten_list_to_set(novel_df["helixer"].str.split(",").tolist()))

    if args.combined:
        unique_prots = novel
    elif "helixer" in args.pfamout_path:
        unique_prots = helixer_only
    elif "braker" in args.pfamout_path:
        unique_prots = braker_only

    genes = find_biologically_interesting_genes(args.pfamout_path, unique_prots)

    print(f"\nBiologically interesting gene candidates ({len(genes)} total):")
    for g in natsorted(genes):
        print(g)

    if args.output:
        with open(args.output, "w") as f:
            for line in natsorted(genes):
                f.write(f"{line};\n")
