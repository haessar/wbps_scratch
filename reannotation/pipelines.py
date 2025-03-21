from collections import Counter
import os
import os.path
import pickle
import re

from gffutils.exceptions import FeatureNotFoundError
import numpy as np
from tqdm import tqdm

from reannotation.utils import extract_accessions_from_transcript, extract_accessions_from_tsv
from utils.generic import makedirs

MAX_GENE_DISTANCE = 2000
MIN_PIDENT = 95


def interpro_accession_pipeline_all_tools(wbps_species, og_df, wbps_col, tool_cols, interproscan_dir):
    acc_tally_no_tool = []
    acc_tally_one_plus_tool_shared = []
    acc_tally_one_plus_tool_novel = []
    for _, row in og_df.iterrows():
        # WBPS transcripts without an orthologue with any tool
        if not row[wbps_col] is np.nan and all(row[tool_col] is np.nan for tool_col in tool_cols):
            for tid in (p.split("transcript_")[1].strip() for p in row[wbps_col].split(",")):
                tran = wbps_species.db["transcript:" + tid]
                acc_tally_no_tool.extend(acc for acc, _ in extract_accessions_from_transcript(tran))
        # WBPS transcripts sharing at least one orthologue with an automated tool
        elif not row[wbps_col] is np.nan and any(not row[tool_col] is np.nan for tool_col in tool_cols):
            for tid in (p.split("transcript_")[1].strip() for p in row[wbps_col].split(",")):
                tran = wbps_species.db["transcript:" + tid]
                acc_tally_one_plus_tool_shared.extend(acc for acc, _ in extract_accessions_from_transcript(tran))
        # Novel transcripts from any automated tool
        elif row[wbps_col] is np.nan and any(not row[tool_col] is np.nan for tool_col in tool_cols):
            for tool_col in tool_cols:
                row_novel_accs = set()
                if row[tool_col] is not np.nan:
                    for tid in (p.strip() for p in row[tool_col].split(",")):
                        tsv_path = os.path.join(interproscan_dir, tid + ".fa.tsv")
                        if os.path.isfile(tsv_path) and os.stat(tsv_path).st_size != 0:
                            for acc, _ in extract_accessions_from_tsv(tsv_path):
                                row_novel_accs.add(acc)
                acc_tally_one_plus_tool_novel.extend(row_novel_accs)
    return acc_tally_no_tool, acc_tally_one_plus_tool_shared, acc_tally_one_plus_tool_novel


def interpro_accession_pipeline(wbps_db, og_df, wbps_col, tool_col, interproscan_dir, prefix="transcript"):
    novel_transcripts = set()
    missed_transcripts = set()
    shared_orth = {}
    acc_tally_novel = []
    acc_tally_missed = []
    acc_tally_shared = []
    for _, row in og_df.iterrows():
        if all(row[[tool_col, wbps_col]].isna()):
            continue
        # WBPS transcripts without a shared orthologue with tool
        if row[tool_col] is np.nan and not row[wbps_col] is np.nan:
            for tid in (p.split(prefix + "_")[1].strip() for p in row[wbps_col].split(",")):
                missed_transcripts.add(tid)
                tran = wbps_db[prefix + ":" + tid]
                acc_tally_missed.extend(acc for acc, _ in extract_accessions_from_transcript(tran))
        # Novel transcripts from tool
        elif row[wbps_col] is np.nan and not row[tool_col] is np.nan:
            novel_transcripts.update(p.strip() for p in row[tool_col].split(","))
            for t in row[tool_col].split(","):
                tsv_path = os.path.join(interproscan_dir, t + ".fa.tsv")
                if os.path.isfile(tsv_path) and os.stat(tsv_path).st_size != 0:
                    for acc, _ in set(extract_accessions_from_tsv(tsv_path)):
                        acc_tally_novel.append(acc)
        # WBPS transcripts sharing at least one orthologue with automated
        else:
            shared_orth[row[tool_col]] = row[wbps_col]
            for tid in (p.split(prefix + "_")[1].strip() for p in row[wbps_col].split(",")):
                tran = wbps_db[prefix + ":" + tid]
                acc_tally_shared.extend(acc for acc, _ in extract_accessions_from_transcript(tran))

    return acc_tally_shared, acc_tally_missed, acc_tally_novel, missed_transcripts


def count_prod_word_occurrence_for_signif_accs(acc_sig_diff_iter, acc_all_occ_iter, acc_product_dict):
    """
    For an iterable of significant accessions "acc_sig_diff_iter", count words encountered
    in their product descriptions (using dict "acc_product_dict") from the iterable of all
    encountered accessions "acc_all_occ_iter". Prints most commonly occurring words.
    """
    words = []
    for acc in acc_sig_diff_iter:
        words.extend(re.split(r'\W', acc_product_dict[acc].lower()) * Counter(acc_all_occ_iter)[acc])
    print(Counter(words).most_common())


def find_adjacent_merged_genes(seq_id_map, species1, transcript1, blast1, transcript2, prefix):
    left_gene = right_gene = None
    left_genes = sorted(
        species1.db.region(
            seqid=transcript1.seqid,
            start=transcript1.start - MAX_GENE_DISTANCE,
            end=transcript1.end,
            strand=transcript1.strand, featuretype="gene"),
        key=lambda x: x.start, reverse=False)
    left_genes = [g for g in left_genes if g not in list(species1.db.parents(transcript1, featuretype="gene"))]
    if len(left_genes) > 0:
        left_genes = [g for g in left_genes if any(c.featuretype == "mRNA" for c in species1.db.children(g))]
        if left_genes:
            left_gene = left_genes[-1]
    right_genes = sorted(
        species1.db.region(
            seqid=transcript1.seqid,
            start=transcript1.start,
            end=transcript1.end + MAX_GENE_DISTANCE,
            strand=transcript1.strand, featuretype="gene"),
        key=lambda x: x.end, reverse=False)
    right_genes = [g for g in right_genes if g not in list(species1.db.parents(transcript1, featuretype="gene"))]
    if len(right_genes) > 0:
        right_genes = [g for g in right_genes if any(c.featuretype == "mRNA" for c in species1.db.children(g))]
        if right_genes:
            right_gene = right_genes[0]
    if left_gene or right_gene:
        tool_col_slice = blast1[blast1["other_transcript_id"]==seq_id_map[transcript2.id]]
        filter_ids = []
        if left_gene:
            tid = list(species1.db.children(left_gene, featuretype="mRNA"))[0].id.split("gene:")[-1].split(".")[0]
            filter_ids.extend([v for k, v in seq_id_map.map.items() if k.startswith(tid)])
        if right_gene:
            tid = list(species1.db.children(right_gene, featuretype="mRNA"))[0].id.split("gene:")[-1].split(".")[0]
            filter_ids.extend([v for k, v in seq_id_map.map.items() if k.startswith(tid)])
        tool_col_slice = tool_col_slice[tool_col_slice["transcript_id"].isin(filter_ids)]
        if not tool_col_slice.empty:
            tool_col_slice = tool_col_slice[tool_col_slice["pident"] > MIN_PIDENT]
            if not tool_col_slice.empty:
                good_left = good_right = False
                if left_gene:
                    if left_gene.end - transcript2.start >= 0:
                        good_left = True
                if right_gene:
                    if not right_gene.start - transcript2.end <= 0:
                        good_right = True
                if any([good_left, good_right]):
                    return list(map(lambda x: x.split(prefix + ":")[-1], list(map(seq_id_map.get, tool_col_slice.transcript_id.to_list())) + [transcript1.id]))
    return []


def suspicious_orthologue_pipeline(og_df, wbps_col, tool_col, species_list, seq_id_map, wbps_prefix="transcript", tool_prefix="transcript"):
    split = {}
    merged = {}
    tool_species = species_list.get_species_with_data_label(tool_col)
    wbps_species = species_list.get_species_with_data_label(wbps_col)
    wbps_blast = wbps_species.blast_slice.compute()
    tool_blast = tool_species.blast_slice.compute()
    for _, row in tqdm(og_df.iterrows(), total=len(og_df)):
        if not row[tool_col] is np.nan and not row[wbps_col] is np.nan:
            # Selecting just the first orthologue for simplicity
            wbps_id = list(map(str.strip, row[wbps_col].split(",")))[0].split(wbps_prefix + "_")[1]
            tool_id = list(map(str.strip, row[tool_col].split(",")))[0]
            wbps_transcript = wbps_species.db[wbps_prefix + ":" + wbps_id]
            wbps_cds_exons = list(wbps_species.db.children(wbps_transcript, featuretype="CDS"))
            wbps_prot_len = wbps_species.get_amino_acid_count(wbps_cds_exons)
            try:
                tool_transcript = tool_species.db[tool_id.split(tool_prefix + "_")[-1]]
            except FeatureNotFoundError:
                tool_transcript = tool_species.db[tool_prefix + ":" + tool_id.split(tool_prefix + "_")[-1]]
            tool_cds_exons = list(tool_species.db.children(tool_transcript, featuretype="CDS"))
            tool_prot_len = tool_species.get_amino_acid_count(tool_cds_exons)
            if len(wbps_cds_exons) != len(tool_cds_exons) or abs(1 - wbps_prot_len / tool_prot_len) > 0.1:
                split_transcripts = find_adjacent_merged_genes(seq_id_map, wbps_species, wbps_transcript, wbps_blast, tool_transcript, wbps_prefix)
                if split_transcripts:
                    split[tool_id] = split_transcripts
                merged_transcripts = find_adjacent_merged_genes(seq_id_map, tool_species, tool_transcript, tool_blast, wbps_transcript, tool_prefix)
                if merged_transcripts:
                    merged[wbps_id] = merged_transcripts

    # Some of the "merged" and "split" genes seem to be subsets.
    subset_count = 0
    genuine_merged = {}
    for k, v in merged.items():
        try:
            t1_range = range(tool_species.db[v[0]].start, tool_species.db[v[0]].end)
        except FeatureNotFoundError:
            t1_range = range(tool_species.db[tool_prefix + ":" + v[0]].start, tool_species.db[tool_prefix + ":" + v[0]].end)
        try:
            t2_range = range(tool_species.db[v[1]].start, tool_species.db[v[1]].end)
        except FeatureNotFoundError:
            t2_range = range(tool_species.db[tool_prefix + ":" + v[1]].start, tool_species.db[tool_prefix + ":" + v[1]].end)
        if set(t1_range).issubset(t2_range) or set(t2_range).issubset(t1_range):
            subset_count += 1
        else:
            genuine_merged[k] = v
    subset_count = 0
    genuine_split = {}
    for k, v in split.items():
        try:
            t1_range = range(wbps_species.db[wbps_prefix + ":" + v[0]].start, wbps_species.db[wbps_prefix + ":" + v[0]].end)
            t2_range = range(wbps_species.db[wbps_prefix + ":" + v[1]].start, wbps_species.db[wbps_prefix + ":" + v[1]].end)
        except FeatureNotFoundError:
            t1_range = range(wbps_species.db[v[0]].start, wbps_species.db[v[0]].end)
            t2_range = range(wbps_species.db[v[1]].start, wbps_species.db[v[1]].end)
        if set(t1_range).issubset(t2_range) or set(t2_range).issubset(t1_range):
            subset_count += 1
        else:
            genuine_split[k] = v
    return genuine_merged, genuine_split


def pickle_cache_suspicious_orthologue_pipeline(tool, sp_prefix, *args, **kwargs):
    merged_path = os.path.join("data", "tmp", f"{sp_prefix}_{tool}_merged.pickle")
    split_path = os.path.join("data", "tmp", f"{sp_prefix}_{tool}_split.pickle")
    if os.path.isfile(merged_path) and os.path.isfile(split_path):
        with open(merged_path, "rb") as f:
            merged = pickle.load(f)
        with open(split_path, "rb") as f:
            split = pickle.load(f)
    else:
        merged, split = suspicious_orthologue_pipeline(*args, **kwargs)
        with open(merged_path, 'wb') as f:
            pickle.dump(merged, f, protocol=pickle.HIGHEST_PROTOCOL)
        with open(split_path, 'wb') as f:
            pickle.dump(split, f, protocol=pickle.HIGHEST_PROTOCOL)
    return merged, split


def novel_orthologue_pipeline(og_df, wbps_col, tool_col, species_list, out_dir="data/novel_orthologue_sequences/"):
    makedirs(out_dir)
    count = 0
    tool_species = species_list.get_species_with_data_label(tool_col)
    for _, row in tqdm(og_df.iterrows(), total=len(og_df)):
        if not row[tool_col] is np.nan and row[wbps_col] is np.nan:
            count += 1
            # Selecting just the first orthologue for simplicity
            tool_id = list(map(str.strip, row[tool_col].split(",")))[0]
            try:
                tool_transcript = tool_species.db[tool_id.split("transcript_")[-1]]
            except FeatureNotFoundError:
                tool_transcript = tool_species.db["transcript:" + tool_id.split("transcript_")[-1]]
            with open(os.path.join(out_dir, tool_id + ".fa"), 'a') as f:
                f.write(">" + tool_transcript.id + "\n")
                f.write(tool_species.get_protein_sequence(tool_transcript.id).strip("*"))
