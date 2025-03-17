from collections import Counter, defaultdict
import configparser
import csv
import dataclasses
from datetime import datetime
from itertools import combinations
import os
import os.path
import re

from matplotlib import pyplot as plt
from minineedle import needle
import pandas as pd
from pygenometracks.tracks.BedTrack import DEFAULT_BED_COLOR
import pygenometracks.tracksClass

from utils.generic import flatten_list_to_list, flatten_list_to_set

CHROM_LABEL = "X"
MAX_NEEDLE_SEQ_LEN = 200


def init_orthogroup_df(og_path):
    df = pd.read_csv(og_path, delimiter="\t")
    return df.iloc[df.isnull().sum(1).sort_values(ascending=True).index]


def format_output_table_path(results_label, hog=None, clade=None):
    return f"""data/schistosome_orthogroups/{results_label}/table_{
        "_".join(filter(None, (
            datetime.now().strftime("%Y_%m_%d_%H_%M_%S"),
            str(hog) if hog else None,
            "clade" + str(clade) if clade else None
    )))}.tsv"""


@dataclasses.dataclass
class Counts:
    exon: dict
    gene: list
    transcript: list
    clade_exons: defaultdict
    prot_amino_acids: dict


class Plotter:
    def __init__(self, do_plot, overwrite, path, conf_path, tmp_dir):
        self.do_plot = do_plot
        self.overwrite = overwrite
        self.path = path
        self.conf_path = conf_path
        self.tmp_dir = tmp_dir
        self.parser = configparser.ConfigParser() if do_plot else None
        self.max_end = 0
        self.skip = False

    def plot_tracks(self):
        with open(self.conf_path, "w") as ff:
            self.parser.write(ff)
        trp = pygenometracks.tracksClass.PlotTracks(self.conf_path, dpi=150, plot_regions=[(CHROM_LABEL, 0, self.max_end)])
        fig = trp.plot(self.path, CHROM_LABEL, 0, self.max_end)
        fig.clear()
        plt.close(fig)

    def tmp_bed_file_path(self, prefix):
        return os.path.join(self.tmp_dir, f"tmp_{prefix}.bed")

    def write_bed_for_single_transcript_exons(self, species, exons):
        with open(self.tmp_bed_file_path(species.prefix.lower()), "w") as f:
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

    def bed_track_config(self, species, transcript):
        tid = transcript.id.split('transcript:')[1] if transcript.id.startswith("transcript:") else transcript.id
        return {
            "file": self.tmp_bed_file_path(species.prefix.lower()),
            "type": "bed",
            "height": "3",
            "title": f"{species.abbr} {species.name}: \n {tid}",
            "fontsize": 10,
            "file_type": "bed",
            "line_width": 1,
            "display": "collapsed",
            "color": DEFAULT_BED_COLOR if transcript.strand == "-" else "red"
        }

    def __enter__(self):
        if self.do_plot:
            if os.path.exists(self.conf_path) and os.path.exists(self.path) and not self.overwrite:
                print(f"{self.path} already exists.")
                self.skip = True
            else:
                self.skip = False
        else:
            self.skip = True
        return self

    def __exit__(self, *args, **kwargs):
        if not self.skip:
            self.plot_tracks()


class OrthoGroup:
    def __init__(self, label, species_list, table_path, table_cols, seq_id_map):
        self.label = label
        self.species_list = species_list
        self.table_path = table_path
        self.table_cols = table_cols
        self.seq_id_map = seq_id_map
        self.selected_transcripts = {}
        self.worst_pair = tuple()
        self.worst_transcript = ""
        self.blast_pident = None
        self.align_pident = None
        self.filtered_blast = pd.DataFrame()
        self.counts = Counts(
            exon={},
            gene=[],
            transcript=[],
            clade_exons=defaultdict(list),
            prot_amino_acids={}
        )

    def process(self, row, prefix_cut="", plotter=None, **kwargs):
        load_blast = kwargs.get("load_blast", False)
        global_ident = kwargs.get("global_ident")
        clade = kwargs.get("clade", None)
        self.ingest_species_data(row, prefix_cut, load_blast, plotter)
        if global_ident:
            if global_ident == "needle":
                self.find_transcript_with_worst_global_alignment(clade)
            elif global_ident == "infer":
                self.find_transcript_with_worst_infered_global_alignment(clade)
        if load_blast:
            self.find_transcripts_with_worst_blast(clade)

    def ingest_species_data(self, row, prefix_cut="", load_blast=False, plotter=None):
        for sp in self.species_list:
            prot_ids_string = row[sp.prot_meta.label]
            if isinstance(prot_ids_string, float):
                continue
            prot_ids = [p.strip() for p in prot_ids_string.split(",")]
            if prefix_cut:
                to_cut = [tid for tid in prot_ids if tid.startswith(prefix_cut)]
                if to_cut:
                    to_leave = list(set(prot_ids).difference(to_cut))
                    prot_ids = list(map(lambda p: p.split(prefix_cut)[1], prot_ids)) + to_leave
            transcript_ids = set(prot_ids)
            gene_ids = set(p.split(".")[0] for p in prot_ids)
            self.counts.transcript.append(len(transcript_ids))
            self.counts.gene.append(len(gene_ids))
            other_prots_labels = row.index.intersection([s.prot_meta.label for s in self.species_list if s.prot_meta.label != sp.prot_meta.label]).to_list()
            other_transcript_ids = flatten_list_to_set(list(map(str.strip, p.split(", "))) for p in row.filter(items=other_prots_labels).to_list())
            if prefix_cut:
                to_cut = [tid for tid in other_transcript_ids if tid.startswith(prefix_cut)]
                if to_cut:
                    to_leave = list(other_transcript_ids.difference(to_cut))
                    other_transcript_ids = list(map(lambda p: p.split(prefix_cut)[1], to_cut)) + to_leave
                else:
                    other_transcript_ids = list(other_transcript_ids)
            transcript, exons = sp.select_transcript(transcript_ids, other_transcript_ids, self.seq_id_map, load_blast)
            transcript_id = transcript.id.split(":")[1] if ":" in transcript.id else transcript.id
            self.selected_transcripts[sp] = transcript_id
            self.counts.clade_exons[sp.clade].append(len(exons))
            self.counts.prot_amino_acids[transcript_id] = sp.get_amino_acid_count(exons)
            self.counts.exon[transcript_id] = len(exons)
            if load_blast:
                self.filtered_blast = pd.concat([self.filtered_blast, sp.blast_slice[sp.blast_slice["transcript_id"]==self.seq_id_map[transcript_id]]])
            if plotter and not plotter.skip and not sp.skip_plot:
                plotter.max_end = max(plotter.max_end, transcript.end - transcript.start)
                plotter.write_bed_for_single_transcript_exons(sp, exons)
                plotter.parser[f"test bed {sp.prefix}"] = plotter.bed_track_config(sp, transcript)

    def find_transcripts_with_worst_blast(self, clade=None):
        if self.worst_pair:
            for sp, tid in self.selected_transcripts.items():
                if tid in self.worst_pair:
                    worst_batch = sp.blast_slice[(sp.blast_slice["transcript_id"] == self.seq_id_map[tid]) &
                                                 (sp.blast_slice["other_transcript_id"].isin(map(self.seq_id_map.get, self.worst_pair)))]
                    break
        else:
            if not clade:
                selected_transcripts = self.selected_transcripts
            else:
                species_ids = self.species_list.get_species_ids_for_clade(clade)
                selected_transcripts = {k: v for k, v in self.selected_transcripts.items() if self.seq_id_map[v].startswith(species_ids)}
            batch = self.filtered_blast[
                (self.filtered_blast["transcript_id"].isin(map(self.seq_id_map.get, selected_transcripts.values()))) &
                (self.filtered_blast["other_transcript_id"].isin(map(self.seq_id_map.get, selected_transcripts.values())))]
            culprit = Counter(flatten_list_to_list(batch.sort_values("pident").head(len(selected_transcripts) - 1)[["transcript_id"]].values)).most_common(1)[0][0]
            worst_batch = batch[(batch["transcript_id"] == culprit) | (batch["other_transcript_id"] == culprit)].sort_values(["pident", "bitscore"]).head(1)
            self.worst_pair = tuple(worst_batch[["transcript_id", "other_transcript_id"]].map(self.seq_id_map.get).values[0])
            self.worst_transcript = self.seq_id_map[culprit]
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
            if min(len(s) for s in seqs) > MAX_NEEDLE_SEQ_LEN:
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
            blast = self.filtered_blast[(self.filtered_blast["transcript_id"]==self.seq_id_map[tr1]) &
                                    (self.filtered_blast["other_transcript_id"]==self.seq_id_map[tr2])]
            align_length = blast["length"].astype(int)
            if align_length.empty:
                continue
            lens = [self.counts.prot_amino_acids[tr1], self.counts.prot_amino_acids[tr2]]
            pident = float(align_length / max(lens) * 100)
            alignments[pair] = pident
        batch = sorted(alignments.items(), key=lambda x: x[1])[:len(orthologue_seq_ids) - 1]
        culprit = Counter(flatten_list_to_list([p for p, _ in batch])).most_common(1)[0][0]
        self.worst_pair, min_pident = [p for p in batch if culprit in p[0]][0]
        self.worst_transcript = culprit
        self.align_pident = min_pident

    def extract_product(self, acc_product=None):
        acc_list = []
        for sp, tid in self.selected_transcripts.items():
            info = sp.db["transcript:" + tid].attributes.get("info")
            if info:
                print(info[0])
                acc_list.extend(re.findall(r'IPR\d+', info[0]))
                if acc_product is not None:
                    for acc in info[0].split("\n"):
                        acc_product[re.search(r'IPR\d+', acc).group()] = acc.split("description:")[1]
        return acc_list

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        if self.selected_transcripts:
            data = {
                "HOG": self.label,
                "selected_transcripts": ", ".join(self.selected_transcripts.values()),
                "exon_counts": ", ".join(map(str, self.counts.exon.values())),
                "protein_lengths": ", ".join(map(str, self.counts.prot_amino_acids.values())),
                "gene_counts": ", ".join(map(str, self.counts.gene)),
                "transcript_counts": ", ".join(map(str, self.counts.transcript)),
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
