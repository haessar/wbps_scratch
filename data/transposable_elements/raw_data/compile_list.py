import csv
import json


with open("data/transposable_elements/raw_data/1.json") as f:
    prev_te_pfams = json.loads(f.read())

te_pfams = {}
with open("data/transposable_elements/raw_data/2.csv") as f:
    for idx, row in enumerate(csv.reader(f)):
        if idx == 0:
            continue
        pfam = row[1]
        desc = row[0]
        if pfam in te_pfams:
            assert desc == te_pfams[pfam]
        te_pfams[pfam] = desc

for pfam in set(te_pfams.keys()).intersection(prev_te_pfams.keys()):
    if te_pfams[pfam] != prev_te_pfams[pfam]:
        print(f"{pfam}: {te_pfams[pfam]} != {prev_te_pfams[pfam]}")

te_pfams.update(prev_te_pfams)

with open("data/transposable_elements/transposable_element_pfams.json", "w") as write_file:
    json.dump(te_pfams, write_file, sort_keys=True, indent=4)
