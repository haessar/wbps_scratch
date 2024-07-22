#! /usr/bin/env python3
import argparse
import os.path

from analysis.orthogroups import Species
from utils.gffutils import init_db


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_gff3')
    args = parser.parse_args()

    db = init_db(args.input_gff3, os.path.splitext(os.path.basename(args.input_gff3))[0] + ".db")
    for gene in db.all_features(featuretype="gene"):
        transcripts = list(db.children(gene, featuretype="mRNA"))
        if len(transcripts) > 1:
            longest_prot_length = 0
            for transcript in transcripts:
                cds_exons = db.children(transcript, featuretype="CDS")
                prot_length = Species.get_amino_acid_count(cds_exons)
                if prot_length > longest_prot_length:
                    longest_transcript = transcript
                    longest_prot_length = prot_length
        elif len(transcripts) == 1:
            longest_transcript = transcripts[0]
        else:
            # Non-coding RNA
            continue
        print(gene)
        print(longest_transcript)
        for c in db.children(longest_transcript):
            print(c)
