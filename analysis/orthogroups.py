from abc import ABC, abstractmethod
import argparse
from collections import Counter, UserList, defaultdict
import configparser
import csv
from datetime import datetime
from glob import glob
from itertools import combinations
import os
import os.path
import re
import sys

import dask.dataframe as dd
from matplotlib import pyplot as plt
from minineedle import needle
from natsort import natsorted
import pandas as pd
import pygenometracks.tracksClass
from pygenometracks.tracks.BedTrack import DEFAULT_BED_COLOR
import pysam
from tqdm import tqdm
from tqdm.dask import TqdmCallback

from utils.generic import flatten_list_to_list, flatten_list_to_set, makedirs
from utils.gffutils import init_db

CHROM_LABEL = "X"
WBPS_RELEASE = "WBPS19"
MAX_NEEDLE_SEQ_LEN = 200


def main(args, species_list):
    global SEQ_ID_MAP
    SEQ_ID_MAP = SequenceIDMapping(args.wd_path)
    df = init_orthogroup_df(args.hog_path)
    
    do_plot = args.do_plot or bool(args.hog)
    table_path = format_output_table_path(args.results_label, args.hog, args.clade)
    makedirs(table_path)
    table_cols = ["HOG", "selected_transcripts", "exon_counts", "protein_lengths", "gene_counts", "transcript_counts"]
    if args.global_ident or args.load_blast:
        table_cols += ["worst_transcript", "worst_pair", "blast_pident"]
    if args.global_ident:
        table_cols += ["align_pident"]
    conf_dir = os.path.join("data", "configs", args.results_label, "")
    plot_dir = os.path.join("plots", args.results_label, "")
    if do_plot:
        makedirs([conf_dir, plot_dir])
    if args.resume is not None:
        output_dir = os.path.dirname(table_path)
        # If path to resume given explicitly...
        if args.resume:
            if not os.path.exists(args.resume):
                raise Exception("--resume path {} does not exist.".format(args.resume))
            else:
                file_to_resume = args.resume
                if os.path.basename(os.path.dirname(file_to_resume)) != os.path.basename(output_dir):
                    raise Exception("Ensure you are resuming from the same set of OrthoFinder results.")
        # ... else use latest output file.
        else:
            output_files = os.listdir(output_dir)
            dtpat = re.compile("\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2}")
            file_to_resume = sorted(output_files, key=lambda x: datetime.strptime(dtpat.search(x).group(0), '%Y_%m_%d_%H_%M_%S'))[-1]
        output_df = pd.read_csv(os.path.join(output_dir, file_to_resume), sep="\t")
        if output_df.columns.to_list() != table_cols:
            raise Exception("Ensure arguments match the latest run you want to --resume from.")
        df = df[~df["HOG"].isin(output_df["HOG"])]
        table_path = file_to_resume
    else:
        if args.hog:
            df = df[df["HOG"] == args.hog]
        pd.DataFrame(columns=table_cols).to_csv(table_path, mode="a", index=False, header=True, sep="\t")
    for _, row in tqdm(df.iterrows(), total=len(df.dropna())):
        if any(row.isna()) and not bool(args.hog):
            break
        with OrthoGroup(
            row=row,
            species_list=species_list,
            do_plot=do_plot,
            overwrite=args.overwrite_plot,
            table_path=table_path,
            table_cols=table_cols,
            conf_dir=conf_dir,
            plot_dir=plot_dir,
            args=args
        ) as og:
            og.process()


def parse_args():
    parser = argparse.ArgumentParser(
        prog="analyse_schistosome_orthogroups.py",
        description="visualise and gather stats for orthogroups from OrthoFinder output",
    )
    parser.add_argument('of_out_dir')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--hog', type=str, default=None)
    group.add_argument('--resume', '-r', nargs='?', const='')
    parser.add_argument('--do-plot', '-p', action='store_true')
    parser.add_argument('--overwrite-plot', action='store_true')
    parser.add_argument('--load-blast', '-l', action='store_true',
                        required="--global-ident" in sys.argv and sys.argv[sys.argv.index("--global-ident") + 1] == "infer")
    parser.add_argument('--global-ident', '-g', choices=[None, 'needle', 'infer'], default=None)
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
        else:
            return self.map[item]
    
    def get(self, item):
        return self.__getitem__(item)


class SpeciesIDMapping:
    def __init__(self, wd_path):
        self.map = {}
        with open(os.path.join(wd_path, "SpeciesIDs.txt"), "r") as f:
            for l in f:
                sid, prot_path = l.strip().split(": ")
                data_label = prot_path.split(".protein.fa")[0]
                self.map[data_label] = int(sid)
    
    def __getitem__(self, item):
        return self.map[item]


class SpeciesList(UserList):
    def __init__(self, initlist, args):
        self.data = []
        tmp_dir = os.path.join("data", "tmp", args.results_label, "")
        makedirs(tmp_dir)
        blast_columns = ("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                "sstart", "send", "evalue", "bitscore")        
        for sp in initlist:
            sp.species_id = SpeciesIDMapping(args.wd_path)[sp.data_label]
            seen_pairs = []
            blast_paths = []
            for bp in natsorted(glob(os.path.join(args.wd_path, "Blast*"))):
                pair = sorted(map(int, re.findall(r'\d+', os.path.basename(bp))))
                if int(sp.species_id) in pair and pair not in seen_pairs and len(set(pair)) == 2:
                    if args.load_blast:
                        print(f"loading {bp}...")
                    blast_paths.append(bp)
                    seen_pairs.append(pair)
            blastout = dd.read_csv(blast_paths, names=blast_columns, delimiter="\t", dtype={'bitscore': 'float64'})
            sp.wd_path = args.wd_path
            sp.tmp_bed_path = os.path.join(tmp_dir, "tmp_{}.bed".format(sp.prefix.lower()))
            sp.slice_blast_data(blastout)
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
    default_transcript_selection_method = "_get_longest_transcript"

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
        self.blast_slice = None

    @property
    @abstractmethod
    def abbr(self):
        raise NotImplementedError()

    @property
    @abstractmethod
    def genus(self):
        raise NotImplementedError()
    
    def slice_blast_data(self, blastout):
        left = blastout[(blastout["qseqid"].str.startswith(str(self.species_id))) & ~(blastout["sseqid"].str.startswith(str(self.species_id)))]
        left = left.rename(columns={"qseqid": "transcript_id", "sseqid": "other_transcript_id"})
        right = blastout[(blastout["sseqid"].str.startswith(str(self.species_id))) & ~(blastout["qseqid"].str.startswith(str(self.species_id)))]
        right = right.rename(columns={"qseqid": "other_transcript_id", "sseqid": "transcript_id"})
        df = dd.concat([left, right], axis=0)
        self.blast_slice = df

    def select_transcript(self, prot_ids, other_prot_ids=set(), load_blast=False):
        if load_blast and self.blast_slice is not None and other_prot_ids:
            if hasattr(self.blast_slice, "compute"):
                with TqdmCallback(desc=f"compute blast slice for species {self.species_id}"):
                    self.blast_slice = self.blast_slice.compute()
            return self._get_transcript_with_best_blast(prot_ids, other_prot_ids)
        return getattr(self, self.default_transcript_selection_method)(prot_ids)

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
    
    def _get_transcript_with_best_blast(self, transcript_ids, other_transcript_ids):
        filt = self.blast_slice[(self.blast_slice["transcript_id"].isin(map(SEQ_ID_MAP.get, transcript_ids))) &
                                (self.blast_slice["other_transcript_id"].isin(map(SEQ_ID_MAP.get, other_transcript_ids)))]
        tid = SEQ_ID_MAP[filt.groupby("transcript_id").median(["pident", "bitscore"]).sort_values(["pident", "bitscore"], ascending=False).head(1).index[0]]
        tran = self.db["transcript:" + tid]
        exons = list(self.db.children("transcript:" + tid, featuretype="CDS"))
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
    def __init__(self, row, species_list, do_plot, overwrite, table_path, table_cols, conf_dir, plot_dir, args):
        self.label = row["HOG"]
        self.row = row
        self.species_list = species_list
        self.do_plot = do_plot
        self.skip_plot = True
        self.overwrite = overwrite
        self.table_path = table_path
        self.table_cols = table_cols
        self.conf_path = os.path.join(conf_dir, self.label + ".ini")
        self.plot_path = os.path.join(plot_dir, self.label + ".png")
        self.parser = configparser.ConfigParser() if self.do_plot else None
        self.max_end = 0
        self.selected_transcripts = {}
        self.clade_exons = defaultdict(list)
        self.prot_lengths = {}
        self.exon_counts = {}
        self.gene_counts = []
        self.transcript_counts = []
        self.worst_pair = tuple()
        self.worst_transcript = ""
        self.blast_pident = None
        self.align_pident = None
        self.filtered_blast = pd.DataFrame()
        self.load_blast = args.load_blast
        self.global_ident = args.global_ident
        self.clade = args.clade
    
    def process(self):
        self.ingest_species_data()
        if not self.skip_plot:
            self.plot_tracks()
        if self.global_ident:
            if self.global_ident == "needle":
                self.find_transcript_with_worst_global_alignment(self.clade)
            elif self.global_ident == "infer":
                self.find_transcript_with_worst_infered_global_alignment(self.clade)
        if self.load_blast:
            self.find_transcripts_with_worst_blast(self.clade)
    
    def ingest_species_data(self):
        for sp in self.species_list:
            prot_ids = self.row[sp.prots_label]
            if type(prot_ids) == float:
                continue
            transcript_ids = set(p.strip() for p in prot_ids.split(","))
            gene_ids = set(p.strip().split(".")[0] for p in prot_ids.split(","))
            self.transcript_counts.append(len(transcript_ids))
            self.gene_counts.append(len(gene_ids))
            other_prots_labels = [c for c in self.row.index if c.startswith(sp.genus) and c != sp.prots_label]
            other_transcript_ids = flatten_list_to_set(list(map(str.strip, p.split(", "))) for p in self.row.filter(items=other_prots_labels).to_list())
            transcript, exons = sp.select_transcript(transcript_ids, other_transcript_ids, self.load_blast)
            transcript_id = transcript.id.split(":")[1]
            self.selected_transcripts[sp] = transcript_id
            self.max_end = max(self.max_end, transcript.end - transcript.start)
            self.clade_exons[sp.clade].append(len(exons))
            self.prot_lengths[transcript_id] = sp.get_amino_acid_count(exons)
            self.exon_counts[transcript_id] = len(exons)
            if self.load_blast:
                self.filtered_blast = pd.concat([self.filtered_blast, sp.blast_slice[sp.blast_slice["transcript_id"]==SEQ_ID_MAP[transcript_id]]])
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

    def find_transcripts_with_worst_blast(self, clade=None):
        if self.worst_pair:
            for sp, tid in self.selected_transcripts.items():
                if tid in self.worst_pair:
                    worst_batch = sp.blast_slice[(sp.blast_slice["transcript_id"] == SEQ_ID_MAP[tid]) &
                                                 (sp.blast_slice["other_transcript_id"].isin(map(SEQ_ID_MAP.get, self.worst_pair)))]
                    break
        else:
            if not clade:
                selected_transcripts = self.selected_transcripts
            else:
                species_ids = self.species_list.get_species_ids_for_clade(clade)
                selected_transcripts = {k: v for k, v in self.selected_transcripts.items() if SEQ_ID_MAP[v].startswith(species_ids)}
            batch = self.filtered_blast[
                (self.filtered_blast["transcript_id"].isin(map(SEQ_ID_MAP.get, selected_transcripts.values()))) & 
                (self.filtered_blast["other_transcript_id"].isin(map(SEQ_ID_MAP.get, selected_transcripts.values())))]
            culprit = Counter(flatten_list_to_list(batch.sort_values("pident").head(len(selected_transcripts) - 1)[["transcript_id"]].values)).most_common(1)[0][0]
            worst_batch = batch[(batch["transcript_id"] == culprit) | (batch["other_transcript_id"] == culprit)].sort_values(["pident", "bitscore"]).head(1)
            self.worst_pair = tuple(worst_batch[["transcript_id", "other_transcript_id"]].map(SEQ_ID_MAP.get).values[0])
            self.worst_transcript = SEQ_ID_MAP[culprit]
        try:
            self.blast_pident = float(worst_batch["pident"])
        except TypeError:
            print(f"No BLAST result for {self.worst_pair}.")

    def find_transcript_with_worst_global_alignment(self, clade=None):
        orthologue_seq_ids = self.selected_transcripts.items()
        if clade:
            orthologue_seq_ids = dict([i for i in orthologue_seq_ids if i[0].clade == clade])
        alignments = {}
        safe_tid = set()
        for pair in combinations(orthologue_seq_ids, 2):
            pair_tids = tuple(tid for _, tid in pair)
            if safe_tid.issuperset(pair_tids):
                continue
            seqs = [sp.get_protein_sequence(tid) for sp, tid in pair]
            # Takes a real performance hit if protein sequences are too long
            if min([len(s) for s in seqs]) > MAX_NEEDLE_SEQ_LEN:
                return
            alg = needle.NeedlemanWunsch(*seqs)
            alg.align()
            pident = alg.get_identity()
            # If a given pair of sequences have good global alignment, chances are neither is the culprit.
            if pident > 90:
                safe_tid.update(pair_tids)
            alignments[pair_tids] = pident
        batch = sorted(alignments.items(), key=lambda x: x[1])[:len(orthologue_seq_ids) - 1]
        culprit = Counter(flatten_list_to_list([p for p, _ in batch])).most_common(1)[0][0]
        self.worst_pair, min_pident = [p for p in batch if culprit in p[0]][0]
        self.worst_transcript = culprit
        self.align_pident = min_pident

    def find_transcript_with_worst_infered_global_alignment(self, clade=None):
        if clade:
            orthologue_seq_ids = [tid for sp, tid in self.selected_transcripts.items() if sp.clade == clade]
        else:
            orthologue_seq_ids = self.selected_transcripts.values()
        alignments = {}
        for pair in combinations(orthologue_seq_ids, 2):
            tr1, tr2 = pair[0], pair[1]
            blast = self.filtered_blast[(self.filtered_blast["transcript_id"]==SEQ_ID_MAP[tr1]) &
                                    (self.filtered_blast["other_transcript_id"]==SEQ_ID_MAP[tr2])]
            align_length = blast["length"].astype(int)
            if align_length.empty:
                continue
            lens = [self.prot_lengths[tr1], self.prot_lengths[tr2]]
            pident = float(align_length / max(lens) * 100)
            alignments[pair] = pident
        batch = sorted(alignments.items(), key=lambda x: x[1])[:len(orthologue_seq_ids) - 1]
        culprit = Counter(flatten_list_to_list([p for p, _ in batch])).most_common(1)[0][0]
        self.worst_pair, min_pident = [p for p in batch if culprit in p[0]][0]
        self.worst_transcript = culprit
        self.align_pident = min_pident

    def __enter__(self):
        if self.do_plot:
            if os.path.exists(self.conf_path) and os.path.exists(self.plot_path) and not self.overwrite:
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
                "exon_counts": ", ".join(map(str, self.exon_counts.values())),
                "protein_lengths": ", ".join(map(str, self.prot_lengths.values())),
                "gene_counts": ", ".join(map(str, self.gene_counts)),
                "transcript_counts": ", ".join(map(str, self.transcript_counts)),
            }
            if self.worst_transcript:
                data.update({
                    "worst_transcript": self.worst_transcript,
                    "worst_pair": ", ".join(self.worst_pair),
                    "blast_pident": self.blast_pident,
                })
            if self.align_pident:
                data.update({
                    "align_pident": self.align_pident,
                })
            df = pd.DataFrame([data], columns=self.table_cols)
            df.to_csv(self.table_path, mode="a", index=False, header=False, sep="\t")
