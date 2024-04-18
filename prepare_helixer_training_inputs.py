#! /usr/bin/env python3
from glob import glob
import os
import os.path
import sys


def extract_species(s):
    return "_".join(s.split("_")[:2])


def extract_accession(s):
    return "".join(s.split("_")[2:]).lower()


def make_symlinks_for_species_file(species_file, type="training", start=1):
    with open(species_file, "r") as f:
        species_list = f.read().splitlines()
    
    for idx, species in enumerate(species_list, start=start):
        sp = extract_species(species)
        src_paths = glob(os.path.join(h5s_dir, "{}*".format(sp), "*.h5"))
        if len(src_paths) == 0:
            print("Couldn't determine src path for {}".format(species), file=sys.stderr)
            continue
        elif len(src_paths) > 1:
            acc = extract_accession(species)
            src_paths = [f for f in glob if "_".join((sp, acc)) in f]
            assert len(src_paths) == 1
        src_path = src_paths[0]
        os.symlink(
            os.path.relpath(
                src_path,
                train_dir
            ),
            os.path.join(train_dir, "{}_data.species_{:02d}.h5".format(type, idx))
        )
    return idx


if __name__ == "__main__":
    train_set_file = sys.argv[1]
    valid_set_file = sys.argv[2]
    base_dir = sys.argv[3]

    h5s_dir = os.path.join(base_dir, "h5s")
    train_dir = os.path.join(base_dir, "train")

    os.makedirs(train_dir, exist_ok=True)

    last_num = make_symlinks_for_species_file(train_set_file)
    make_symlinks_for_species_file(valid_set_file, type="validation", start=last_num + 1)
