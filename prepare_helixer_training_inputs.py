#! /usr/bin/env python3
from glob import glob
import os
import os.path
import sys


def extract_species(s):
    return "_".join(s.split("_")[:2])


def extract_accession(s):
    return "".join(s.split("_")[2:]).lower()


def make_symlinks_for_species_file(species_file, h5s_dir, train_dir, typ="training", start=1):
    with open(species_file, "r") as f:
        species_list = f.read().splitlines()

    for idx, species in enumerate(species_list, start=start):
        try:
            sp = extract_species(species)
            src_paths = glob(os.path.join(h5s_dir, f"{sp}*", "*.h5"))
            if len(src_paths) == 0:
                raise RuntimeError(f"Couldn't determine src path for {species}")
            elif len(src_paths) > 1:
                src_matching_accession_paths = []
                for fp in src_paths:
                    acc = extract_accession(os.path.basename(os.path.dirname(fp)))
                    if acc == extract_accession(species):
                        src_matching_accession_paths.append(fp)
                src_paths = src_matching_accession_paths
                assert len(src_paths) == 1
            src_path = src_paths[0]
            os.symlink(
                os.path.relpath(
                    src_path,
                    train_dir
                ),
                # pylint: disable=C0209
                os.path.join(train_dir, "{}_data.species_{:02d}.h5".format(typ, idx))
            )
        except FileExistsError:
            pass
        except Exception as e:
            print(f"{species}: {e}", file=sys.stderr)
    # pylint: disable=W0631
    return idx


if __name__ == "__main__":
    train_set_file = sys.argv[1]
    valid_set_file = sys.argv[2]
    train_dir = sys.argv[3]

    h5s_dir = os.path.join(os.path.dirname(train_dir), "h5s")

    os.makedirs(train_dir, exist_ok=True)

    last_num = make_symlinks_for_species_file(train_set_file, h5s_dir, train_dir)
    make_symlinks_for_species_file(valid_set_file, h5s_dir, train_dir, typ="validation", start=last_num + 1)
