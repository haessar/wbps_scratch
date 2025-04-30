#! /usr/bin/env python3
from collections import defaultdict
from glob import glob
import os
import os.path
import re
import shutil
import sys


if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    pattern = r"^(?P<species>[a-z0-9\_]+)\.(?P<acc>[^\.]+)\.(?P<rel>[^\.]+)\..*"
    filetree = defaultdict(set)

    for file in os.listdir(input_dir):
        match = re.match(pattern, file)
        try:
            if match.group("species"):
                filetree[match["species"]].add(match["acc"])
        except AttributeError:
            print(f"{file} does not match input file pattern", file=sys.stderr)

    for species, accs in filetree.items():
        for acc in accs:
            if len(accs) > 1:
                output_path = os.path.join(output_dir, "_".join((species, acc)), "input")
            else:
                output_path = os.path.join(output_dir, species, "input")
            os.makedirs(output_path, exist_ok=True)
            for f in glob(os.path.join(input_dir, f"{species}*{acc}*")):
                shutil.move(f, output_path)
            if len(os.listdir(output_path)) != 2:
                print(f"There are {len(os.listdir(output_path))}!=2 files in {output_path}.", file=sys.stderr)
            if not glob(os.path.join(output_path, "*genomic.fa*")):
                print(f"Missing genomic.fa file in {output_path}", file=sys.stderr)
            if not glob(os.path.join(output_path, "*annotations.gff3*")):
                print(f"Missing annotations.gff3 file in {output_path}", file=sys.stderr)
