from collections import Counter
import os
import os.path

from scipy.stats import fisher_exact
from .utils import extract_accessions_from_tsv


def fisher_exact_for_two_lists_of_accessions(l1, l2):
    more_frequent = {}
    less_frequent = {}
    no_significant_difference = {}
    not_occurring = []
    for acc, freq1 in Counter(l1).items():
        try:
            freq2 = dict(Counter(l2))[acc]
        except KeyError:
            not_occurring.append(acc)
            continue
        other1 = len(l1) - freq1
        other2 = len(l2) - freq2
        res = fisher_exact([[freq1, freq2], [other1, other2]])
        if res.pvalue < 0.05:
            if res.statistic > 1:
                more_frequent[acc] = res.statistic
            elif res.statistic < 1:
                less_frequent[acc] = res.statistic
        else:
            no_significant_difference[acc] = res.statistic
    return {
        "more_frequent": more_frequent, 
        "as_expected": no_significant_difference,
        "less_frequent": less_frequent,
        "not_occurring": not_occurring
    }


def count_transcripts_with_accession(accession, prefix, interproscan_dir):
    count = 0
    unique = set()
    for tsv in os.listdir(interproscan_dir):
        tsv_path = os.path.join(interproscan_dir, tsv)
        if os.stat(tsv_path).st_size != 0 and tsv.startswith(prefix):
            for acc, _ in extract_accessions_from_tsv(tsv_path):
                if acc == accession:
                    count += 1
                    unique.add(tsv)
                    break
    return count, unique
