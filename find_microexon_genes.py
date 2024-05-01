import os.path
import sys

import gffutils

helixer = ("data/from_MARS/Schistosoma_mansoni_helixer_full.gff3", "Sm_helixer_full.db")
wbps = ("data/from_WBP/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3", "Sm_wbps.db")
braker3 = ("data/from_MARS/Schistosoma_mansoni_braker3.gff3", "Sm_braker3.db")

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
    
    ME_gene_count = 0

    if output_dir:
        output_path = os.path.join(output_dir, os.path.splitext(INPUTS[1])[0] + ".md")
        f = open(output_path, "w")

    for gene in db.all_features(featuretype="gene"):
        is_ME_gene = False
        exons = list(db.children(gene.id, featuretype="exon"))
        if len(exons) > 1:
            previous_complied = False
            consecutive = 1
            for exon in exons:
                if is_microexon(len(exon)):
                    if previous_complied:
                        consecutive += 1
                        if consecutive == REQ_CONSECUTIVE:
                            is_ME_gene = True
                            break
                    previous_complied = True
                    continue
                previous_complied = False
        if is_ME_gene:
            print(gene)
            if output_dir:
                f.write(gene.id + ":&emsp;")
                f.write("-".join([format_microexon_length(len(e)) for e in exons]))
                f.write("\\\n")
            ME_gene_count += 1
    print(ME_gene_count)

    if output_dir:
        f.close()
