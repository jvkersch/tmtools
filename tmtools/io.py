""" Simple tools to read chain data from PDB files.
"""
import os
import warnings

from Bio.PDB.PDBParser import PDBParser
try:
    from Bio.PDB.Polypeptide import protein_letters_3to1
except ImportError:
    # pre 1.80
    from Bio.PDB import protein_letters_3to1

import numpy as np


def get_structure(fname):
    fname = str(fname)
    b = os.path.basename(fname)
    sname, _ = os.path.splitext(b)

    parser = PDBParser(PERMISSIVE=1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure(sname, fname)

    return structure


def get_residue_data(chain):
    """Extract residue coordinates and sequence from PDB chain.

    Uses the coordinates of the Cα atom as the center of the residue.

    """
    coords = []
    seq = []
    for residue in chain.get_residues():
        if "CA" in residue.child_dict:
            coords.append(residue.child_dict["CA"].coord)
            seq.append(protein_letters_3to1[residue.resname])

    return np.vstack(coords), "".join(seq)
