from collections import defaultdict
import configparser
import csv
import os.path
import statistics

import matplotlib 
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import pandas as pd
import pygenometracks.tracksClass
from pygenometracks.tracks.BedTrack import DEFAULT_BED_COLOR
from tqdm import tqdm

from utils.generic import flatten_list_to_list, flatten_list_to_set
from utils.gffutils import init_db



hog_path = "data/from_MARS/OrthoFinder/Results_May10/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"

df = pd.read_csv(hog_path, delimiter="\t")
df = df.iloc[df.isnull().sum(1).sort_values(ascending=True).index]

DO_PLOT = False
OVERWRITE = True
CHROM_LABEL = "X"


class Species:
    data_dir = os.path.join("data", "from_WBPS")
    tmp_dir = os.path.join("data", "tmp")
    db_dir = "db"

    def __init__(self, name, gff, prots_label):
        self.name = name
        self.prefix = self.abbr + self.name.lower()[:2]
        self.prots_label = self._replace_wildcard_prefix(prots_label)
        self.last_tid = ""
        self.tmp_bed_path = os.path.join(self.tmp_dir, "tmp_{}.bed".format(self.prefix.lower()))
        db_path = os.path.join(self.db_dir, "{}_wbps.db".format(self.prefix))
        gff_path = os.path.join(self.data_dir, self._replace_wildcard_prefix(gff))
        self.db = init_db(gff_path, db_path)
    
    def _replace_wildcard_prefix(self, s):
        if s.startswith("*"):
            return s.replace("*", f"{self.genus}_{self.name}")
        return s

    def _parse_first_transcript_id(self, prot_ids):
        transcript_ids = set(p.strip() for p in prot_ids.split(","))
        for tid in transcript_ids:
            self.last_tid = tid
            return tid

    def write_bed_for_transcript(self, prot_ids):
        with open(self.tmp_bed_path, "w") as f:
            tid = "transcript:" + self._parse_first_transcript_id(prot_ids)
            tran = self.db[tid]
            exons = sorted(self.db.children(tid, featuretype="CDS"), key=lambda x: x.start, reverse=tran.strand=="-")
            for jdx, cds in enumerate(exons):
                if jdx == 0:
                    offset = cds.end if cds.strand == "-" else cds.start
                    writer = csv.writer(f, delimiter="\t")
                if cds.strand == "+":
                    tsvdata = [CHROM_LABEL, cds.start - offset, cds.end - offset]
                elif cds.strand == "-":
                    tsvdata = [CHROM_LABEL, abs(cds.end - offset), abs(cds.start - offset)]
                # Avoid 0 length exons
                if tsvdata[2] - tsvdata[1] == 0:
                    continue
                writer.writerow(tsvdata)
        return tran, exons

    def bed_config(self, strand):
        return {
            "file": self.tmp_bed_path,
            "type": "bed",
            "height": "3",
            "title": f"{self.abbr} {self.name}: \n {self.last_tid}",
            "fontsize": 10,
            "file_type": "bed",
            "line_width": 1,
            "display": "collapsed",
            "color": DEFAULT_BED_COLOR if strand == "-" else "red" 
        }


class Schistosoma(Species):
    abbr = "S"
    genus = "schistosoma"


class HaematobiumClade(Schistosoma):
    clade = 1


class MansoniClade(Schistosoma):
    clade = 2


class JaponicumClade(Schistosoma):
    clade = 3


species = [
    HaematobiumClade("bovis", "*.TD2_PRJEB44434.WBPS19.annotations.gff3", "*.TD2_PRJEB44434.WBPS19.protein"),
    HaematobiumClade("curassoni", "*.PRJEB44434.WBPS19.annotations.gff3", "*.PRJEB44434.WBPS19.protein"),
    HaematobiumClade("guineensis", "*.PRJEB44434.WBPS19.annotations.gff3", "*.PRJEB44434.WBPS19.protein"),
    HaematobiumClade("haematobium", "*.TD2_PRJEB44434.WBPS19.annotations.gff3", "*.TD2_PRJEB44434.WBPS19.protein"),
    HaematobiumClade("intercalatum", "*.TD2_PRJEB44434.WBPS19.annotations.gff3", "*.TD2_PRJEB44434.WBPS19.protein"),
    MansoniClade("rodhaini", "*.TD2_PRJEB44434.WBPS19.annotations.gff3", "*.TD2_PRJEB44434.WBPS19.protein"),
    MansoniClade("mansoni", "*.PRJEA36577.WBPS19.annotations.gff3", "*.PRJEA36577.WBPS19.protein"),
    JaponicumClade("japonicum", "*.PRJNA520774.WBPS19.annotations.gff3", "*.PRJNA520774.WBPS19.protein"),
]

true_orths = set()
one_exon_diff_orths = set()
multi_exon_diff_orths = set()

for idx, row in tqdm(df.iterrows(), total=len(df.dropna())):
    if row["HOG"] == "N0.HOG0004956":
        print()
    if any(row.isna()):
        break
    conf_path = f"data/configs/{row['HOG']}.ini"
    plot_path = f"plots/{row['HOG']}.png"
    if DO_PLOT:
        if os.path.exists(conf_path) and os.path.exists(plot_path) and not OVERWRITE:
            continue
        fig, ax = plt.subplots(1)
    parser = configparser.ConfigParser()
    max_end = 0
    clade_exons = defaultdict(list)
    for sp in species:
        prot_ids = row[sp.prots_label]
        transcript, exons = sp.write_bed_for_transcript(prot_ids)
        max_end = max(max_end, transcript.end - transcript.start)
        clade_exons[sp.clade].append(len(exons))
        parser["test bed {} {}".format(idx, sp.prefix)] = sp.bed_config(transcript.strand)
    uniq_exon_nums = flatten_list_to_set(clade_exons.values())
    if len(uniq_exon_nums) == 1:
        true_orths.add(row['HOG'])
    else:
        # med_num_exons = statistics.median(flatten_list_to_list(clade_exons.values()))
        max_diff = max(uniq_exon_nums) - min(uniq_exon_nums)
        if max_diff == 1:
            one_exon_diff_orths.add(row["HOG"])
        elif max_diff > 1:
            multi_exon_diff_orths.add(row["HOG"])
        # for clade, exon_nums in clade_exons.items():
        #     if len(set(exon_nums)) > 1:
        #         print()
    if DO_PLOT:
        with open(conf_path, "w") as ff:
            parser.write(ff)
        ax.set_xlim(0, max_end)
        trp = pygenometracks.tracksClass.PlotTracks(conf_path, dpi=150, plot_regions=[(CHROM_LABEL, 0, max_end)])
        current_fig = trp.plot(plot_path, CHROM_LABEL, 0, max_end)
