from orthologue_analysis.orthogroups import init_orthogroup_df
from orthologue_analysis.utils import orthofinder_paths
from reannotation.pipelines import interpro_accession_pipeline
from utils.gffutils import init_db


if __name__ == "__main__":
    results_label = "Results_Jul31"
    wbps_ann_path = "data/from_WBPS/haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3"
    braker_path = "data/from_MARS/Haemonchus_contortus_braker3_full.gff3"
    db = init_db(wbps_ann_path, "db/Hcon_wbps.db")
    of =  orthofinder_paths(results_label)

    wbps_col = "Hcon_LT"
    braker_col = "Hcon_braker3_LT"
    hog_df = init_orthogroup_df(of["orthogroups"])

    interpro_accession_pipeline(db, hog_df, braker_col, wbps_col)
