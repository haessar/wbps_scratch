from collections import defaultdict
import json
import os.path

from orthologue_analysis.species import HaemonchusFromTool
from reannotation.utils import populate_accession_product_dict

if __name__ == "__main__":
    mars_data_dir = os.path.join("data", "from_MARS", "")

    species = HaemonchusFromTool("contortus", data_dir=mars_data_dir, data_label="Hcon_LT", prot_filename_suffix=".fa")
    acc_product = defaultdict(dict)
    if os.path.exists("acc_product_hcontortus.json"):
        with open("acc_product_hcontortus.json", "r") as f:
            acc_product = json.loads(f.read())
    acc_product = populate_accession_product_dict(species, acc_product)
    with open("acc_product_hcontortus.json", "w") as f:
        json.dump(acc_product, f, sort_keys=True, indent=4)
