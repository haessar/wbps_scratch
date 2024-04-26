import os.path

import gffutils

helixer_gff3 = "data/from_MARS/Schistosoma_mansoni_helixer_full.gff3"
helixer_db_path = "Sm_helixer_full.db"
wbps_gff3 = "/home/will/tests/Schistosoma/mansoni/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
wbps_db_path = "Sm_wbps.db"

gff_path = helixer_gff3
db_path = helixer_db_path

MAX_LENGTH = 36
REQ_CONSECUTIVE = 2

if __name__ == "__main__":
    if not os.path.exists(db_path):
        db = gffutils.create_db(gff_path, db_path, merge_strategy="create_unique")

    else:
        db = gffutils.FeatureDB(db_path)
    
    ME_gene_count = 0
    for gene in db.all_features(featuretype="gene"):
        is_ME_gene = False
        exons = list(db.children(gene.id, featuretype="exon"))
        if len(exons) > 1:
            previous_complied = False
            consecutive = 1
            for exon in exons:
                if 0 < len(exon) <= MAX_LENGTH:
                    if len(exon) % 3 == 0:
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
            ME_gene_count += 1
    print(ME_gene_count)

