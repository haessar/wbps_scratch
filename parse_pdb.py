import os
import os.path
from statistics import mean, median
import sys

import Bio
import Bio.PDB


def read_pdb(pdbcode, pdbfilenm):
    """
    Read a PDB structure from a file.
    :param pdbcode: A PDB ID string
    :param pdbfilenm: The PDB file
    :return: a Bio.PDB.Structure object or None if something went wrong
    """
    try:
        pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
        struct = pdbparser.get_structure(pdbcode, pdbfilenm)
        return struct
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None


def parse_af(afid):
    return afid.split('-')[1]


if __name__ == "__main__":
    batchdir = sys.argv[1]
    files = os.listdir(batchdir)
    for pdb in files:
        all_plDDT = []
        try:
            struct = read_pdb(parse_af(pdb), os.path.join(batchdir, pdb))
            for res in struct.get_residues():
                pLDDT = set(a.bfactor for a in res.child_list)
                assert len(pLDDT) == 1
                all_plDDT.append(pLDDT.pop())
            print(",".join((pdb, str(int(mean(all_plDDT))), str(int(median(all_plDDT))))))
        except Exception as e:
            print(pdb, file=sys.stderr)
