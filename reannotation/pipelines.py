from collections import Counter
import contextlib
import re

from gffutils.exceptions import FeatureNotFoundError
import numpy as np
from tqdm import tqdm

from reannotation.utils import extract_accessions_from_transcript
from reannotation.statistics import fisher_exact_for_two_lists_of_accessions


def interpro_accession_pipeline(db, hog_df, wbps_col, tool_col):
    # InterPro accessions from all mRNA features in WBPS annotation.
    acc_product = {}
    acc_list1 = []
    for tran in db.all_features(featuretype="mRNA"):
        for acc, prod in extract_accessions_from_transcript(tran):
            acc_product[acc] = prod
            acc_list1.append(acc)

    tool_only = set()
    wbps_only = set()
    shared_orth = {}
    acc_list3 = []
    for _, row in hog_df.iterrows():
        if all(row[[tool_col, wbps_col]].isna()):
            continue
        if row[tool_col] is np.nan and not row[wbps_col] is np.nan:
            for tid in (p.split("transcript_")[1].strip() for p in row[wbps_col].split(",")):
                wbps_only.add(tid)
                tran = db["transcript:" + tid]
                acc_list3.extend(acc for acc, _ in extract_accessions_from_transcript(tran))
        elif row[wbps_col] is np.nan and not row[tool_col] is np.nan:
            tool_only.update(p.strip() for p in row[tool_col].split(","))
        else:
            shared_orth[row[tool_col]] = row[wbps_col]

    # Full WBPS accession list with acc_list3 removed
    acc_list1a = list((Counter(acc_list1) - Counter(acc_list3)).elements())

    # Most common InterPro accessions in WBPS gene missed by annotation tool
    l3_more_expressed, l3_evenly_expressed, l3_less_expressed = fisher_exact_for_two_lists_of_accessions(acc_list1a, acc_list3)
    for acc, freq in Counter(acc_list3).most_common():
        try:
            if acc in l3_more_expressed:
                print(f"{acc}: {freq} - significantly more frequent ({l3_more_expressed[acc]}) - {acc_product[acc]}")
            elif acc in l3_less_expressed:
                print(f"{acc}: {freq} - significantly less frequent ({l3_less_expressed[acc]}) - {acc_product[acc]}")
            else:
                print(f"{acc}: {freq} - occurring as expected ({l3_evenly_expressed[acc]}) - {acc_product[acc]}")
        except KeyError:
            pass

    # Filter accessions with high frequency and low test statistic.
    filt = [k for k, v in Counter(acc_list3).items() if v >= 10]
    filt = [k for k, v in l3_more_expressed.items() if v < 0.1 and k in filt]

    for tran in db.all_features(featuretype="mRNA"):
        if tran.id.strip("transcript:") in wbps_only:
            with contextlib.redirect_stdout(None):
                tran_filt_accs = set(acc for acc, _ in extract_accessions_from_transcript(tran) if acc in filt)
            if tran_filt_accs:
                print(f"{tran.seqid} - {tran.id.strip('transcript:')} - {tran_filt_accs}")

    return {
        "l3_more_expressed": l3_more_expressed,
        "l3_evenly_expressed": l3_evenly_expressed,
        "l3_less_expressed": l3_less_expressed,
        "acc_list3": acc_list3,
        "acc_product": acc_product,
    }


def count_prod_word_occurrence_for_signif_accs(acc_sig_diff_iter, acc_all_occ_iter, acc_product_dict):
    """
    For an iterable of significant accessions "acc_sig_diff_iter", count words encountered
    in their product descriptions (using dict "acc_product_dict") from the iterable of all
    encountered accessions "acc_all_occ_iter". Prints most commonly occurring words.
    """
    words = []
    for acc in acc_sig_diff_iter:
        words.extend(re.split(r'\W', acc_product_dict[acc].lower()) * Counter(acc_all_occ_iter)[acc])
    print(Counter(words).most_common())

