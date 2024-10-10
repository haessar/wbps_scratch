import argparse
from datetime import datetime
import os.path
import re
import sys

import pandas as pd
from tqdm import tqdm

from utils.generic import makedirs
from .orthogroups import format_output_table_path, init_orthogroup_df, OrthoGroup, Plotter
from .utils import SequenceIDMapping


class OrthologueAnalysisResumeError(Exception):
    pass


def main(args, species_list):
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
    tmp_dir = os.path.join("data", "tmp", args.results_label, "")
    if do_plot:
        makedirs([conf_dir, plot_dir, tmp_dir])
    if args.resume is not None:
        output_dir = os.path.dirname(table_path)
        # If path to resume given explicitly...
        if args.resume:
            if not os.path.exists(args.resume):
                raise OrthologueAnalysisResumeError(f"--resume path {args.resume} does not exist.")
            file_to_resume = args.resume
            if os.path.basename(os.path.dirname(file_to_resume)) != os.path.basename(output_dir):
                raise OrthologueAnalysisResumeError("Ensure you are resuming from the same set of OrthoFinder results.")
        # ... else use latest output file.
        else:
            output_files = os.listdir(output_dir)
            dtpat = re.compile(r"\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2}")
            file_to_resume = sorted(output_files, key=lambda x: datetime.strptime(dtpat.search(x).group(0), '%Y_%m_%d_%H_%M_%S'))[-1]
        output_df = pd.read_csv(os.path.join(output_dir, file_to_resume), sep="\t")
        if output_df.columns.to_list() != table_cols:
            raise OrthologueAnalysisResumeError("Ensure arguments match the latest run you want to --resume from.")
        df = df[~df["HOG"].isin(output_df["HOG"])]
        table_path = file_to_resume
    else:
        if args.hog:
            df = df[df["HOG"] == args.hog]
        pd.DataFrame(columns=table_cols).to_csv(table_path, mode="a", index=False, header=True, sep="\t")
    for _, row in tqdm(df.iterrows(), total=len(df.dropna())):
        if any(row.isna()) and not bool(args.hog):
            break
        label = row["HOG"]
        with Plotter(
            do_plot=do_plot,
            overwrite=args.overwrite_plot,
            path=os.path.join(plot_dir, label + ".png"),
            conf_path=os.path.join(conf_dir, label + ".ini"),
            tmp_dir=tmp_dir,
        ) as plotter:
            with OrthoGroup(
                label=label,
                species_list=species_list,
                table_path=table_path,
                table_cols=table_cols,
                seq_id_map=SequenceIDMapping(args.wd_path)
            ) as og:
                og.process(row, plotter=plotter, **vars(args))


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
