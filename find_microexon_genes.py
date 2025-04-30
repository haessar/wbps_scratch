#! /usr/bin/env python3
import os.path
import sys

import gffutils

helixer = ("data/from_MARS/Schistosoma_mansoni_helixer_full.gff3", "db/Sman_helixer_full.db")
wbps = ("data/from_WBPS/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3", "db/Sman_wbps.db")
braker3 = ("data/from_MARS/Schistosoma_mansoni_braker3.gff3", "db/Sman_braker3.db")

INPUTS = helixer
MAX_LENGTH = 36
REQ_CONSECUTIVE = 2


def is_microexon(exon_length):
    if 0 < exon_length <= MAX_LENGTH:
        if exon_length % 3 == 0:
            return True
    return False


def format_microexon_length(exon_length):
    if is_microexon(exon_length):
        return "<ins>" + str(exon_length) + "</ins>"
    return str(exon_length)


if __name__ == "__main__":
    output_dir = sys.argv[1] if len(sys.argv) > 1 else None
    if not os.path.exists(INPUTS[1]):
        db = gffutils.create_db(INPUTS[0], INPUTS[1], merge_strategy="create_unique")
    else:
        db = gffutils.FeatureDB(INPUTS[1])

    if output_dir:
        output_path = os.path.join(output_dir, os.path.splitext(INPUTS[1])[0] + ".md")
        f = open(output_path, "w")

    ME_GENE_COUNT = 0
    for gene in db.all_features(featuretype="gene"):
        IS_ME_GENE = False
        exons = list(db.children(gene.id, featuretype="exon"))
        if len(exons) > 1:
            PREVIOUS_COMPLIED = False
            CONSECUTIVE = 1
            for exon in exons:
                if is_microexon(len(exon)):
                    if PREVIOUS_COMPLIED:
                        CONSECUTIVE += 1
                        if CONSECUTIVE == REQ_CONSECUTIVE:
                            IS_ME_GENE = True
                            break
                    PREVIOUS_COMPLIED = True
                    continue
                PREVIOUS_COMPLIED = False
        if IS_ME_GENE:
            print(gene)
            if output_dir:
                f.write(gene.id + ":&emsp;")
                f.write("-".join([format_microexon_length(len(e)) for e in exons]))
                f.write("\\\n")
            ME_GENE_COUNT += 1
    print(ME_GENE_COUNT)

    if output_dir:
        f.close()
