import os.path

from utils.generic import get_project_root


class SequenceIDMapping:
    def __init__(self, wd_path):
        self.map = {}
        with open(os.path.join(wd_path, "SequenceIDs.txt"), "r") as f:
            for l in f:
                sid, info = l.strip().split(": ")
                tid = info.split(" ")[0]
                self.map[tid] = sid
        self.inv_map = {v: k for k, v in self.map.items()}

    def __getitem__(self, item):
        if str(item)[0].isnumeric():
            return self.inv_map[item]
        prefixes = ["", "transcript:", "transcript_"]
        while True:
            try:
                return self.map[prefixes[0] + item]
            except KeyError:
                prefixes.pop(0)

    def get(self, item):
        return self[item]


class SpeciesIDMapping:
    def __init__(self, wd_path, ext):
        self.map = {}
        with open(os.path.join(wd_path, "SpeciesIDs.txt"), "r") as f:
            for l in f:
                sid, prot_path = l.strip().split(": ")
                data_label = prot_path.split(ext)[0]
                self.map[data_label] = int(sid)

    def __getitem__(self, item):
        return self.map[item]


def orthofinder_paths(label, subdir = "Phylogenetic_Hierarchical_Orthogroups"):
    orthofinder_dir = get_project_root() / "data" / "from_MARS" / "OrthoFinder"
    paths = {}
    paths["wd"] = orthofinder_dir / "WorkingDirectory" / label
    if subdir == "Phylogenetic_Hierarchical_Orthogroups":
        paths["orthogroups"] = orthofinder_dir / subdir / label / "N0.tsv"
    elif subdir == "Orthogroups":
        paths["orthogroups"] = orthofinder_dir / subdir / label / "Orthogroups.tsv"
        paths["orthogroups_unassigned_genes"] = orthofinder_dir / subdir / label / "Orthogroups_UnassignedGenes.tsv"
    return paths
