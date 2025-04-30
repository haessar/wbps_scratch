from collections import OrderedDict
from itertools import islice
import os.path


import matplotlib.pyplot as plt

TE_DIR = "data/transposable_elements/"

braker_tes = os.path.join(TE_DIR, "counts_braker.txt")
helixer_tes = os.path.join(TE_DIR, "counts_helixer.txt")

def read_counts(path):
    counts = {}
    with open(path, "r") as f:
        for line in f:
            pfam, count = line.split()
            counts[pfam] = int(count)
    return counts

braker_counts = read_counts(braker_tes)
helixer_counts = read_counts(helixer_tes)


def most_different_counts(c1, c2, num=5):
    differences = {}
    for pfam in c1.keys():
        differences[pfam] = c1[pfam] - c2.get(pfam, 0)
    for idx, (k, v) in enumerate(sorted(differences.items(), key=lambda item: item[1], reverse=True)):
        print(f"{k}: {v} ({c1[k]}-{c2.get(k, 0)})")
        yield k, v
        if idx == num-1:
            break

NUM = 10
print(f"HELIXER - Total TEs: {sum(helixer_counts.values())}")
print(f"Top {NUM} most differentiated from BRAKER:")
mdch = list(most_different_counts(helixer_counts, braker_counts, NUM))
print()
print(f"BRAKER - Total TEs: {sum(braker_counts.values())}")
print(f"Top {NUM} most differentiated from HELIXER:")
mdcb = list(most_different_counts(braker_counts, helixer_counts, NUM))
print()

# plt.barh(y=[k for k, v in mdch], width=[v for k, v in mdch])
# plt.barh(y=[k for k, v in mdcb], width=[v for k, v in mdcb])

toph = OrderedDict(islice(OrderedDict(helixer_counts).items(), 10))
topb = OrderedDict(islice(OrderedDict(braker_counts).items(), 10))

plt.barh(y=[k for k, v in toph.items()], width=[v for k, v in toph.items()])
plt.barh(y=[k for k, v in topb.items()], width=[v for k, v in topb.items()], height=0.5)
