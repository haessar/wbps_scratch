from abc import ABC, abstractmethod
import argparse
from collections import Counter, UserList, defaultdict
import configparser
import csv
from datetime import datetime
import os.path

from matplotlib import pyplot as plt
import pandas as pd
import pygenometracks.tracksClass
from pygenometracks.tracks.BedTrack import DEFAULT_BED_COLOR
import pysam
from tqdm import tqdm

from utils.generic import flatten_list_to_list, makedirs
from utils.gffutils import init_db

CHROM_LABEL = "X"
WBPS_RELEASE = "WBPS19"


def main(args, species_list):
    df = init_orthogroup_df(args.hog_path)
    
    do_plot = args.do_plot or bool(args.hog)
    table_path = format_output_table_path(args.results_label, args.hog, args.clade)
    makedirs(table_path)
    conf_dir = os.path.join("data", "configs", args.results_label, "")
    plot_dir = os.path.join("plots", args.results_label, "")
    if do_plot:
        makedirs([conf_dir, plot_dir])
    if args.hog:
        df = df[df["HOG"] == args.hog]
    for _, row in tqdm(df.iterrows(), total=len(df.dropna())):
        if any(row.isna()) and not bool(args.hog):
            break
        with OrthoGroup(
            row=row,
            species_list=species_list,
            do_plot=do_plot,
            overwrite=args.overwrite,
            table_path=table_path,
            conf_dir=conf_dir,
            plot_dir=plot_dir,
        ) as og:
            og.parse_exons()
            if not og.skip_plot:
                og.plot_tracks()
            if args.load_blast:
                og.find_transcripts_with_lowest_pident(args.clade)


def parse_args():
    parser = argparse.ArgumentParser(
        prog="analyse_schistosome_orthogroups.py",
        description="visualise and gather stats for orthogroups from OrthoFinder output",
    )
    parser.add_argument('of_out_dir')
    parser.add_argument('--hog', type=str, default=None)
    parser.add_argument('--do-plot', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--load-blast', action='store_true')
    parser.add_argument('--clade', type=int, default=None)
    args = parser.parse_args()
    args.hog_path = os.path.join(args.of_out_dir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv")
    args.wd_path = os.path.join(args.of_out_dir, "WorkingDirectory", "")
    args.results_label = os.path.basename(os.path.normpath(args.of_out_dir))
    return args


def init_orthogroup_df(hog_path):
    df = pd.read_csv(hog_path, delimiter="\t")
    return df.iloc[df.isnull().sum(1).sort_values(ascending=True).index]


def format_output_table_path(results_label, hog=None, clade=None):
    return f"""data/schistosome_orthogroups/{results_label}/table_{
        "_".join(filter(None, (
            datetime.now().strftime("%Y_%m_%d_%H_%M_%S"),
            str(hog) if hog else None,
            "clade" + str(clade) if clade else None
    )))}.tsv"""


def init_BLASTOUT():
    global BLASTOUT
    BLASTOUT = \
    pd.DataFrame(
        columns=("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                 "sstart", "send", "evalue", "bitscore")
    )


def init_SEQUENCE_ID_map(wd_path):
    global SEQUENCE_ID
    SEQUENCE_ID = {}
    with open(os.path.join(wd_path, "SequenceIDs.txt"), "r") as f:
        for l in f:
            sid, info = l.strip().split(": ")
            tid = info.split(" ")[0]
            SEQUENCE_ID[tid] = sid


def init_SPECIES_ID_map(wd_path):
    global SPECIES_ID
    SPECIES_ID = {}
    with open(os.path.join(wd_path, "SpeciesIDs.txt"), "r") as f:
        for l in f:
            sid, prot_path = l.strip().split(": ")
            data_label = prot_path.split(".protein.fa")[0]
            SPECIES_ID[data_label] = int(sid)


class SpeciesList(UserList):
    def __init__(self, initlist, args):
        wd_path = args.wd_path if args.load_blast else None
        tmp_dir = os.path.join("data", "tmp", args.results_label, "")
        makedirs(tmp_dir)
        if wd_path:
            init_SEQUENCE_ID_map(wd_path)
            init_SPECIES_ID_map(wd_path)
            init_BLASTOUT()
        self.data = []
        for sp in initlist:
            sp.species_id = SPECIES_ID[sp.data_label] if wd_path else None
            sp.wd_path = wd_path
            sp.tmp_bed_path = os.path.join(tmp_dir, "tmp_{}.bed".format(sp.prefix.lower()))
            if args.load_blast:
                sp.load_blastout()
            self.data.append(sp)

    def get_species_ids_for_clade(self, clade):
        species_ids = set()
        for sp in self:
            if sp.clade == clade:
                species_ids.add(str(sp.species_id))
        return tuple(species_ids)


class Species(ABC):
    data_dir = os.path.join("data", "from_WBPS", "")
    db_dir = "db"
    release = WBPS_RELEASE
    transcript_selection_method = "_get_longest_transcript"

    def __init__(self, name, acc):
        self.name = name
        self.acc = acc
        self.prefix = self.abbr + self.name.lower()[:3]
        self.data_label = f"{self.genus}_{self.name}.{self.acc}.{self.release}"
        self.prots_label = self.data_label + ".protein"
        db_path = os.path.join(self.db_dir, "{}_wbps.db".format(self.prefix))
        gff_path = os.path.join(self.data_dir, ".".join((self.data_label, "annotations", "gff3")))
        self.prots_path = os.path.join(self.data_dir, ".".join((self.data_label, "protein", "fa")))
        self.db = init_db(gff_path, db_path)

    @property
    @abstractmethod
    def abbr(self):
        raise NotImplementedError()

    @property
    @abstractmethod
    def genus(self):
        raise NotImplementedError()

    def load_blastout(self):
        global BLASTOUT
        for sp_id in sorted(SPECIES_ID.values()):
            if sp_id > self.species_id:
                df = pd.read_csv(
                        os.path.join(self.wd_path, f"Blast{self.species_id}_{sp_id}.txt"),
                        names=BLASTOUT.columns,
                        delimiter="\t"
                )
                BLASTOUT = pd.concat([BLASTOUT, df])
        print("Loaded blastout for species {}".format(self.species_id))

    def select_transcript(self, prot_ids):
        return getattr(self, self.transcript_selection_method)(prot_ids)

    def _get_transcript_with_most_exons(self, transcript_ids):
        exon_counts = {}
        for tid in transcript_ids:
            cds = list(self.db.children("transcript:" + tid, featuretype="CDS"))
            exon_counts[len(cds)] = {self.db["transcript:" + tid]: cds}
        for tran, exons in exon_counts[max(exon_counts)].items():
            return tran, sorted(exons, key=lambda x: x.start, reverse=tran.strand=="-")
    
    def _get_longest_transcript(self, transcript_ids):
        prot_lengths = {}
        for tid in transcript_ids:
            cds = list(self.db.children("transcript:" + tid, featuretype="CDS"))
            prot_lengths[self.get_amino_acid_count(cds)] = {self.db["transcript:" + tid]: cds}
        for tran, exons in prot_lengths[max(prot_lengths)].items():
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
    
    @staticmethod
    def get_amino_acid_count(cds_exons):
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


class IndicumClade(Schistosoma):
    clade = 4


class NewSpeciesClade(Schistosoma):
    clade = 5


class OrthoGroup:
    def __init__(self, row, species_list, do_plot, overwrite, table_path, conf_dir, plot_dir):
        self.label = row["HOG"]
        self.row = row
        self.species_list = species_list
        self.do_plot = do_plot
        self.skip_plot = False
        self.overwrite = overwrite
        self.table_path = table_path
        self.conf_path = os.path.join(conf_dir, self.label + ".ini")
        self.plot_path = os.path.join(plot_dir, self.label + ".png")
        self.parser = configparser.ConfigParser() if self.do_plot else None
        self.max_end = 0
        self.selected_transcripts = {}
        self.clade_exons = defaultdict(list)
        self.prot_lengths = []
        self.exon_counts = []
        self.gene_counts = []
        self.transcript_counts = []
        self.min_pident = 0
        self.worst_pair = {}
        self.worst_transcript = {}
    
    def parse_exons(self):
        for sp in self.species_list:
            prot_ids = self.row[sp.prots_label]
            if type(prot_ids) == float:
                continue
            transcript_ids = set(p.strip() for p in prot_ids.split(","))
            gene_ids = set(p.strip().split(".")[0] for p in prot_ids.split(","))
            self.transcript_counts.append(len(transcript_ids))
            self.gene_counts.append(len(gene_ids))
            transcript, exons = sp.select_transcript(transcript_ids)
            self.selected_transcripts[sp] = transcript.id.split(":")[1]
            self.max_end = max(self.max_end, transcript.end - transcript.start)
            self.clade_exons[sp.clade].append(len(exons))
            self.prot_lengths.append(sp.get_amino_acid_count(exons))
            self.exon_counts.append(len(exons))
            if not self.skip_plot:
                sp.write_bed_for_single_transcript_exons(exons)
                self.parser["test bed {}".format(sp.prefix)] = sp.bed_track_config(transcript)
    
    def plot_tracks(self):
        with open(self.conf_path, "w") as ff:
            self.parser.write(ff)
        trp = pygenometracks.tracksClass.PlotTracks(self.conf_path, dpi=150, plot_regions=[(CHROM_LABEL, 0, self.max_end)])
        fig = trp.plot(self.plot_path, CHROM_LABEL, 0, self.max_end)
        fig.clear()
        plt.close(fig)
    
    def find_transcripts_with_lowest_pident(self, clade=None):
        if not BLASTOUT.empty:
            orthologue_seq_ids = [SEQUENCE_ID[tid] for tid in self.selected_transcripts.values()]
            if clade:
                species_ids = self.species_list.get_species_ids_for_clade(clade)
                orthologue_seq_ids = [i for i in orthologue_seq_ids if i.startswith(species_ids)]
            worst_batch = \
                BLASTOUT[(BLASTOUT["qseqid"].isin(orthologue_seq_ids)) & (BLASTOUT["sseqid"].isin(orthologue_seq_ids))] \
                    .sort_values("pident").head(7)
            culprit = Counter(flatten_list_to_list(worst_batch[["qseqid", "sseqid"]].values)).most_common(1)[0][0]
            worst_pair = worst_batch[(worst_batch["qseqid"] == culprit) | (worst_batch["sseqid"] == culprit)].head(1)[["qseqid", "sseqid", "pident"]].values[0]
            seq_id_map = {k: SEQUENCE_ID[k] for k in self.selected_transcripts.values()}
            self.worst_transcript = {sp.data_label: seqid for sp, seqid in self.selected_transcripts.items() if seqid in [k for k, v in seq_id_map.items() if v==culprit]}
            self.worst_pair = {sp.data_label: seqid for sp, seqid in self.selected_transcripts.items() if seqid in (k for k, v in seq_id_map.items() if v in worst_pair[:2])}
            self.min_pident = worst_pair[2]
        else:
            raise Exception("This method requires BLASTOUT to be loaded for each species.")

    def __enter__(self):
        if self.do_plot and os.path.exists(self.conf_path) and os.path.exists(self.plot_path) and not self.overwrite:
            print(f"{self.plot_path} already exists.")
            self.skip_plot = True
        else:
            self.skip_plot = False
        return self

    def __exit__(self, type, value, traceback):
        if self.selected_transcripts:
            data = {
                "HOG": self.label,
                "selected_transcripts": ", ".join(self.selected_transcripts.values()),
                "exon_counts": ", ".join(map(str, self.exon_counts)),
                "protein_lengths": ", ".join(map(str, self.prot_lengths)),
                "gene_counts": ", ".join(map(str, self.gene_counts)),
                "transcript_counts": ", ".join(map(str, self.transcript_counts)),
            }
            if self.worst_transcript:
                data.update({
                    "worst_transcript": ", ".join(self.worst_transcript.values()),
                    "worst_pair": ", ".join(self.worst_pair.values()),
                    "min_pident": self.min_pident,
                })
            df = pd.DataFrame([data])
            df.to_csv(self.table_path, mode="a", index=False, header=False, sep="\t")
