from collections import defaultdict
import json
import os.path

from orthologue_analysis.species import HaemonchusFromTool
from reannotation.utils import populate_accession_product_dict

if __name__ == "__main__":
    filepath = os.path.join("data", "acc_product_hcontortus.json")
    mars_data_dir = os.path.join("data", "from_MARS", "")

    species = HaemonchusFromTool("contortus", data_dir=mars_data_dir, data_label="Hcon_LT", prot_filename_suffix=".fa")
    acc_product = defaultdict(dict)
    if os.path.exists(filepath):
        with open(filepath, "r") as f:
            acc_product = json.loads(f.read())
    acc_product = populate_accession_product_dict(species, acc_product)
    with open(filepath, "w") as f:
        json.dump(acc_product, f, sort_keys=True, indent=4)
