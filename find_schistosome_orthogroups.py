from collections import defaultdict
import configparser
import csv
import os.path
import statistics

from matplotlib import pyplot as plt
import pandas as pd
import pygenometracks.tracksClass
from pygenometracks.tracks.BedTrack import DEFAULT_BED_COLOR
import pysam
from tqdm import tqdm

from utils.generic import flatten_list_to_list, flatten_list_to_set
from utils.gffutils import init_db


hog_path = "data/from_MARS/OrthoFinder/Results_May10/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"

df = pd.read_csv(hog_path, delimiter="\t")
df = df.iloc[df.isnull().sum(1).sort_values(ascending=True).index]

DO_PLOT = True
OVERWRITE = True
CHROM_LABEL = "X"
WBPS_RELEASE = "WBPS19"


class Species:
    data_dir = os.path.join("data", "from_WBPS")
    tmp_dir = os.path.join("data", "tmp")
    db_dir = "db"
    release = WBPS_RELEASE

    def __init__(self, name, acc):
        self.name = name
        self.acc = acc
        self.prefix = self.abbr + self.name.lower()[:2]
        self.prots_label = self.data_label + ".protein"
        self.tmp_bed_path = os.path.join(self.tmp_dir, "tmp_{}.bed".format(self.prefix.lower()))
        db_path = os.path.join(self.db_dir, "{}_wbps.db".format(self.prefix))
        gff_path = os.path.join(self.data_dir, ".".join((self.data_label, "annotations", "gff3")))
        self.prots_path = os.path.join(self.data_dir, ".".join((self.data_label, "protein", "fa")))
        self.db = init_db(gff_path, db_path)
    
    @property
    def data_label(self):
        return f"{self.genus}_{self.name}.{self.acc}.{self.release}"
        
    def get_transcript_with_most_exons(self, prot_ids):
        transcript_ids = set(p.strip() for p in prot_ids.split(","))
        exon_counts = {}
        for tid in transcript_ids:
            cds = list(self.db.children("transcript:" + tid, featuretype="CDS"))
            exon_counts[len(cds)] = {self.db["transcript:" + tid]: cds}
        for tran, exons in exon_counts[max(exon_counts)].items():
            return tran, sorted(exons, key=lambda x: x.start, reverse=tran.strand=="-")

    def write_bed_for_single_transcript_exons(self, exons):
        with open(self.tmp_bed_path, "w") as f:
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

    def bed_track_config(self, transcript):
        return {
            "file": self.tmp_bed_path,
            "type": "bed",
            "height": "3",
            "title": f"{self.abbr} {self.name}: \n {transcript.id.split('transcript:')[1]}",
            "fontsize": 10,
            "file_type": "bed",
            "line_width": 1,
            "display": "collapsed",
            "color": DEFAULT_BED_COLOR if transcript.strand == "-" else "red"
        }
    
    def get_protein_sequence(self, tid):
        # Slow method to obtain amino acids using pysam...
        seq = pysam.faidx(self.prots_path, tid)
        return "".join(seq.split("\n")[1:])
    
    def get_amino_acid_count(self, cds_exons):
        # ...or faster method by inferring from CDS lengths
        return int(sum(len(cds) for cds in cds_exons) / 3)


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
    HaematobiumClade("bovis", "TD2_PRJEB44434"),
    HaematobiumClade("curassoni", "PRJEB44434"),
    HaematobiumClade("guineensis", "PRJEB44434"),
    HaematobiumClade("haematobium", "TD2_PRJEB44434"),
    HaematobiumClade("intercalatum", "TD2_PRJEB44434"),
    MansoniClade("rodhaini", "TD2_PRJEB44434"),
    MansoniClade("mansoni", "PRJEA36577"),
    JaponicumClade("japonicum", "PRJNA520774"),
]

true_orths = set()
one_exon_diff_orths = set()
multi_exon_diff_orths = set()
prot_lengths = defaultdict(list)

for idx, row in tqdm(df.iterrows(), total=len(df.dropna())):
    if any(row.isna()):
        break
    if DO_PLOT:
        conf_path = f"data/configs/{row['HOG']}.ini"
        plot_path = f"plots/{row['HOG']}.png"
        if os.path.exists(conf_path) and os.path.exists(plot_path) and not OVERWRITE:
            continue
        parser = configparser.ConfigParser()
        max_end = 0
    clade_exons = defaultdict(list)
    for sp in species:
        prot_ids = row[sp.prots_label]
        transcript, exons = sp.get_transcript_with_most_exons(prot_ids)
        if DO_PLOT:
            sp.write_bed_for_single_transcript_exons(exons)
            max_end = max(max_end, transcript.end - transcript.start)
            parser["test bed {} {}".format(idx, sp.prefix)] = sp.bed_track_config(transcript)
        clade_exons[sp.clade].append(len(exons))
        # prot_lengths[row['HOG']].append(len(sp.get_protein_sequence(transcript.id.split("transcript:")[1])))
        prot_lengths[row['HOG']].append(sp.get_amino_acid_count(exons))
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
        trp = pygenometracks.tracksClass.PlotTracks(conf_path, dpi=150, plot_regions=[(CHROM_LABEL, 0, max_end)])
        fig = trp.plot(plot_path, CHROM_LABEL, 0, max_end)
        fig.clear()
        plt.close(fig)
