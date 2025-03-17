#! /usr/bin/env python3
import os.path

from orthologue_analysis.utils import SequenceIDMapping, orthofinder_paths
from orthologue_analysis.orthogroups import init_orthogroup_df
from orthologue_analysis.species import SpeciesList, PristionchusFromTool

from ppac_merged_split_run_utils import pickle_cache_suspicious_orthologue_pipeline

wbps_col = "Ppac_LT"
helixer_col = "Ppac_helixer_LT"
results_label = "Results_Aug21"
of = orthofinder_paths(results_label, subdir="Orthogroups")
seq_id_map = SequenceIDMapping(of["wd"])
og_df = init_orthogroup_df(of["orthogroups"])
mars_data_dir = os.path.join("data", "from_MARS", "")
ebi_data_dir = os.path.join("data", "from_EBI", "")


species_list = SpeciesList([
    PristionchusFromTool("pacificus", data_dir=mars_data_dir, data_label="Ppac_LT", prot_filename_suffix=".fa"),
    PristionchusFromTool("pacificus_braker3_reann", data_dir=mars_data_dir, data_label="Ppac_braker3_LT", prot_filename_suffix=".fa"),
    PristionchusFromTool("pacificus_helixer_reann", data_dir=mars_data_dir, data_label="Ppac_helixer_LT", prot_filename_suffix=".fa"),
    PristionchusFromTool("pacificus_anno_reann", data_dir=ebi_data_dir, data_label="Ppac_anno_LT", prot_filename_suffix=".fa")],
    wd_path=of["wd"],
    load_blast=True
)

helixer_merged, helixer_split = pickle_cache_suspicious_orthologue_pipeline("helixer", og_df, wbps_col, helixer_col, species_list, seq_id_map, wbps_prefix="Transcript")

num_genes = len(list(species_list.get_species_with_data_label("Ppac_helixer_LT").db.all_features(featuretype="gene")))
print(f"Helixer: merged={len(helixer_merged)}, split={len(helixer_split)}, total={round(100*(len(helixer_split) + len(helixer_merged)*2)/num_genes, 2)}")
