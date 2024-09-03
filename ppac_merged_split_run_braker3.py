#! /usr/bin/env python3
import os.path

from orthologue_analysis.utils import SequenceIDMapping, orthofinder_paths
from orthologue_analysis.orthogroups import init_orthogroup_df
from orthologue_analysis.species import SpeciesList

from ppac_merged_split_run_utils import PristionchusFromTool, pickle_cache_suspicious_orthologue_pipeline

wbps_col = "Ppac_LT"
braker_col = "Ppac_braker3_LT"
results_label = "Results_Aug21"
of = orthofinder_paths(results_label)
seq_id_map = SequenceIDMapping(of["wd"])
hog_df = init_orthogroup_df(of["orthogroups"])
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

braker_merged, braker_split = pickle_cache_suspicious_orthologue_pipeline("braker", hog_df, wbps_col, braker_col, species_list, seq_id_map, wbps_prefix="Transcript")
