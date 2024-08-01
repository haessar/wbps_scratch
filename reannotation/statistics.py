from collections import Counter

from scipy.stats import fisher_exact


def fisher_exact_for_two_lists_of_accessions(l1, l2):
    more_expressed = {}
    less_expressed = {}
    evenly_expressed = {}
    for acc, freq2 in Counter(l2).items():
        try:
            freq1 = dict(Counter(l1))[acc]
        except KeyError:
            print(f"{acc} not present in any l1 gene")
            continue
        other1 = len(l1) - freq1
        other2 = len(l2) - freq2
        res = fisher_exact([[freq1, freq2], [other1, other2]])
        if res.pvalue < 0.05:
            if res.statistic > 1:
                less_expressed[acc] = res.statistic
            elif res.statistic < 1:
                more_expressed[acc] = res.statistic
        else:
            evenly_expressed[acc] = res.statistic
    return more_expressed, evenly_expressed, less_expressed
