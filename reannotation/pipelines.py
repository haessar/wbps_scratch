from collections import Counter
import contextlib
import os.path
import re

from gffutils.exceptions import FeatureNotFoundError
import numpy as np
from tqdm import tqdm

from reannotation.utils import extract_accessions_from_transcript
from reannotation.statistics import fisher_exact_for_two_lists_of_accessions
from utils.generic import makedirs

MAX_GENE_DISTANCE = 2000
MIN_PIDENT = 95


def interpro_accession_pipeline(db, hog_df, wbps_col, tool_col, min_freq=5):
    # InterPro accessions from all mRNA features in WBPS annotation.
    acc_product = {}
    acc_list1 = []
    for tran in db.all_features(featuretype="mRNA"):
        with contextlib.redirect_stdout(None):
            for acc, prod in extract_accessions_from_transcript(tran):
                acc_product[acc] = prod
                acc_list1.append(acc)

    tool_only = set()
    wbps_only = set()
    shared_orth = {}
    acc_list3 = []
    for _, row in hog_df.iterrows():
        if all(row[[tool_col, wbps_col]].isna()):
            continue
        if row[tool_col] is np.nan and not row[wbps_col] is np.nan:
            for tid in (p.split("transcript_")[1].strip() for p in row[wbps_col].split(",")):
                wbps_only.add(tid)
                tran = db["transcript:" + tid]
                with contextlib.redirect_stdout(None):
                    acc_list3.extend(acc for acc, _ in extract_accessions_from_transcript(tran))
        elif row[wbps_col] is np.nan and not row[tool_col] is np.nan:
            tool_only.update(p.strip() for p in row[tool_col].split(","))
        else:
            shared_orth[row[tool_col]] = row[wbps_col]

    # Full WBPS accession list with acc_list3 removed
    acc_list1a = list((Counter(acc_list1) - Counter(acc_list3)).elements())

    # Find InterPro accessions occurring with significantly different frequency than in control (acc_list1a)
    l3_more_frequent, l3_as_expected, l3_less_frequent, l3_only = fisher_exact_for_two_lists_of_accessions(acc_list3, acc_list1a)

    print("InterPro accessions that are completely missing from control, with high frequency in test:")
    for acc, freq in Counter(acc_list3).most_common():
        if acc in l3_only:
            if freq >= min_freq:
                print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences)")
    print()

    # InterPro accessions missed by automated tool, sorted by greatest odds ratio
    print("InterPro accessions occurring with significantly higher frequency than in control:")
    for acc, stat in sorted(l3_more_frequent.items(), key=lambda x: x[1], reverse=True):
        freq = Counter(acc_list3)[acc]
        if freq >= min_freq:
            print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences, {round(freq/stat)} expected)")
    print()

    print("InterPro accessions occurring as expected with high frequency:")
    for acc, stat in l3_as_expected.items():
        if stat < 1.5 and stat >= 0.5:
            freq = Counter(acc_list3)[acc]
            if freq >= min_freq:
                print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences)")
    print()

    print("InterPro accessions occurring less frequently than expected:")
    for acc, stat in l3_less_frequent.items():
        freq = Counter(acc_list3)[acc]
        print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences, {round(freq/stat)} expected)")
    print()


    # Filter accessions with high frequency and high test statistic [print out is easy to copy/paste for Artemis].
    filt = [acc for acc, freq in Counter(acc_list3).items() if freq >= min_freq]
    filt = [acc for acc, stat in l3_more_frequent.items() if stat > 5 and acc in filt]

    for tran in db.all_features(featuretype="mRNA"):
        if tran.id.strip("transcript:") in wbps_only:
            with contextlib.redirect_stdout(None):
                tran_filt_accs = set(acc for acc, _ in extract_accessions_from_transcript(tran) if acc in filt)
            if tran_filt_accs:
                print(f"{tran.seqid} - {tran.id.strip('transcript:')} - {tran_filt_accs}")

    return {
        "l3_more_frequent": l3_more_frequent,
        "l3_as_expected": l3_as_expected,
        "l3_less_frequent": l3_less_frequent,
        "acc_list3": acc_list3,
        "acc_product": acc_product,
    }


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


def find_adjacent_merged_genes(seq_id_map, species1, transcript1, blast1, id1, transcript2, id2):
    left_gene = right_gene = None
    left_genes = sorted(
        species1.db.region(
            seqid=transcript1.seqid,
            start=transcript1.start - MAX_GENE_DISTANCE,
            end=transcript1.end,
            strand=transcript1.strand, featuretype="gene"),
        key=lambda x: x.start, reverse=False)
    left_genes = [g for g in left_genes if not id1.startswith(g.id.split("gene:")[-1] + ".")]
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
    right_genes = [g for g in right_genes if not id1.startswith(g.id.split("gene:")[-1] + ".")]
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
                    return list(map(lambda x: x.split("transcript:")[-1], list(map(seq_id_map.get, tool_col_slice.transcript_id.to_list())) + [transcript1.id]))
    return []


def suspicious_orthologue_pipeline(hog_df, wbps_col, tool_col, species_list, seq_id_map):
    split = {}
    merged = {}
    tool_species = species_list.get_species_with_data_label(tool_col)
    wbps_species = species_list.get_species_with_data_label(wbps_col)
    wbps_blast = wbps_species.blast_slice.compute()
    tool_blast = tool_species.blast_slice.compute()
    for _, row in tqdm(hog_df.iterrows(), total=len(hog_df)):
        if not row[tool_col] is np.nan and not row[wbps_col] is np.nan:
            # Selecting just the first orthologue for simplicity
            wbps_id = list(map(str.strip, row[wbps_col].split(",")))[0].split("transcript_")[1]
            tool_id = list(map(str.strip, row[tool_col].split(",")))[0]
            wbps_transcript = wbps_species.db["transcript:" + wbps_id]
            wbps_cds_exons = list(wbps_species.db.children(wbps_transcript, featuretype="CDS"))
            wbps_prot_len = wbps_species.get_amino_acid_count(wbps_cds_exons)
            try:
                tool_transcript = tool_species.db[tool_id.split("transcript_")[-1]]
            except FeatureNotFoundError:
                tool_transcript = tool_species.db["transcript:" + tool_id.split("transcript_")[-1]]
            tool_cds_exons = list(tool_species.db.children(tool_transcript, featuretype="CDS"))
            tool_prot_len = tool_species.get_amino_acid_count(tool_cds_exons)
            if len(wbps_cds_exons) != len(tool_cds_exons) or abs(1 - wbps_prot_len / tool_prot_len) > 0.1:
                split_transcripts = find_adjacent_merged_genes(seq_id_map, wbps_species, wbps_transcript, wbps_blast, wbps_id, tool_transcript, tool_id)
                if split_transcripts:
                    split[tool_id] = split_transcripts
                merged_transcripts = find_adjacent_merged_genes(seq_id_map, tool_species, tool_transcript, tool_blast, tool_id, wbps_transcript, wbps_id)
                if merged_transcripts:
                    merged[wbps_id] = merged_transcripts

    # Some of the "merged" and "split" genes seem to be subsets.
    subset_count = 0
    genuine_merged = {}
    for k, v in merged.items():
        try:
            t1_range = range(tool_species.db[v[0]].start, tool_species.db[v[0]].end)
        except FeatureNotFoundError:
            t1_range = range(tool_species.db["transcript:" + v[0]].start, tool_species.db["transcript:" + v[0]].end)
        try:
            t2_range = range(tool_species.db[v[1]].start, tool_species.db[v[1]].end)
        except FeatureNotFoundError:
            t2_range = range(tool_species.db["transcript:" + v[1]].start, tool_species.db["transcript:" + v[1]].end)
        if set(t1_range).issubset(t2_range) or set(t2_range).issubset(t1_range):
            subset_count += 1
        else:
            genuine_merged[k] = v
    subset_count = 0
    genuine_split = {}
    for k, v in split.items():
        t1_range = range(wbps_species.db["transcript:" + v[0]].start, wbps_species.db["transcript:" + v[0]].end)
        t2_range = range(wbps_species.db["transcript:" + v[1]].start, wbps_species.db["transcript:" + v[1]].end)
        if set(t1_range).issubset(t2_range) or set(t2_range).issubset(t1_range):
            subset_count += 1
        else:
            genuine_split[k] = v
    return genuine_merged, genuine_split


def novel_orthologue_pipeline(hog_df, wbps_col, tool_col, species_list):
    out_dir = "data/novel_orthologue_sequences/"
    makedirs(out_dir)
    count = 0
    tool_species = species_list.get_species_with_data_label(tool_col)
    for _, row in tqdm(hog_df.iterrows(), total=len(hog_df)):
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
                f.write(tool_species.get_protein_sequence(tool_transcript.id))
