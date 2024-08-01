from collections import Counter
import contextlib

import numpy as np

from orthologue_analysis.orthogroups import init_orthogroup_df
from orthologue_analysis.utils import SequenceIDMapping, orthofinder_paths
from reannotation.utils import extract_accessions_from_transcript
from reannotation.statistics import fisher_exact_for_two_lists_of_accessions
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
    seq_id_map = SequenceIDMapping(of["wd"])

    # InterPro accessions from all mRNA features in WBPS H. contortus annotation.
    acc_product = {}
    acc_list1 = []
    for tran in db.all_features(featuretype="mRNA"):
        for acc, prod in extract_accessions_from_transcript(tran):
            acc_product[acc] = prod
            acc_list1.append(acc)

    braker_only = set()
    wbps_only = set()
    shared_orth = {}
    acc_list3 = []
    for _, row in hog_df.iterrows():
        if all(row[[braker_col, wbps_col]].isna()):
            continue
        if row[braker_col] is np.nan and not row[wbps_col] is np.nan:
            for tid in (p.split("transcript_")[1].strip() for p in row[wbps_col].split(",")):
                wbps_only.add(tid)
                tran = db["transcript:" + tid]
                acc_list3.extend(acc for acc, _ in extract_accessions_from_transcript(tran))
        elif row[wbps_col] is np.nan and not row[braker_col] is np.nan:
            braker_only.update(p.strip() for p in row[braker_col].split(","))
        else:
            shared_orth[row[braker_col]] = row[wbps_col]

    # Full WBPS accession list with acc_list3 removed
    acc_list1a = list((Counter(acc_list1) - Counter(acc_list3)).elements())

    # Most common InterPro accessions in H contortus gene missed by BRAKER3
    l3_more_expressed, l3_evenly_expressed, l3_less_expressed = fisher_exact_for_two_lists_of_accessions(acc_list1a, acc_list3)
    for acc, freq in Counter(acc_list3).most_common():
        try:
            if acc in l3_more_expressed:
                print(f"{acc}: {freq} - significantly more expressed ({l3_more_expressed[acc]}) - {acc_product[acc]}")
            elif acc in l3_less_expressed:
                print(f"{acc}: {freq} - significantly less expressed ({l3_less_expressed[acc]}) - {acc_product[acc]}")
            else:
                print(f"{acc}: {freq} - expressed as expected ({l3_evenly_expressed[acc]}) - {acc_product[acc]}")
        except KeyError:
            pass

    # Filter accessions with high frequency and 
    filt = [k for k, v in Counter(acc_list3).items() if v >= 10]
    filt = [k for k, v in l3_more_expressed.items() if v < 0.1 and k in filt]

    for tran in db.all_features(featuretype="mRNA"):
        if tran.id.strip("transcript:") in wbps_only:
            with contextlib.redirect_stdout(None):
                tran_filt_accs = set(acc for acc, _ in extract_accessions_from_transcript(tran) if acc in filt)
            if tran_filt_accs:
                print(f"{tran.seqid} - {tran.id.strip('transcript:')} - {tran_filt_accs}")
