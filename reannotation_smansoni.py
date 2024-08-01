from orthologue_analysis.orthogroups import init_orthogroup_df
from orthologue_analysis.utils import SequenceIDMapping, orthofinder_paths
from reannotation.pipelines import interpro_accession_pipeline
from utils.gffutils import init_db


if __name__ == "__main__":
    results_label = "Results_Aug01"
    wbps_ann_path = "data/from_WBPS/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
    braker_path = "data/from_MARS/Schistosoma_mansoni_braker3_full.gff3"
    db = init_db(wbps_ann_path, "db/Sman_wbps.db")
    of = orthofinder_paths(results_label)

    wbps_col = "Sman_LT"
    braker_col = "Sman_braker3_LT"
    hog_df = init_orthogroup_df(of["orthogroups"])
    seq_id_map = SequenceIDMapping(of["wd"])

    interpro_accession_pipeline(db, hog_df, braker_col, wbps_col)
