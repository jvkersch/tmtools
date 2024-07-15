"""Simple tools to read chain data from PDB files."""

import os
import warnings

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser

try:
    from Bio.PDB.Polypeptide import protein_letters_3to1
except ImportError:
    # pre 1.80
    from Bio.PDB import protein_letters_3to1

import numpy as np


def get_structure(fname, format="pdb"):
    """
    Load and parse a structure file.

    Parameters
    ----------
    fname : str
        Path to the structure file.
    format : {'pdb', 'mmcif'}, optional
        Format of the structure file. Default is 'pdb'.

    Returns
    -------
    Bio.PDB.Structure.Structure
        Parsed structure object.

    Raises
    ------
    ValueError
        If an invalid format is provided.

    Notes
    -----
    This function uses Biopython's structure parsers to load the file.
    Warnings are suppressed during parsing.

    Examples
    --------
    >>> structure = get_structure('1abc.pdb')
    >>> structure = get_structure('1xyz.cif', format='mmcif')
    """
    # Check if format provided is valid
    if format not in ["pdb", "mmcif"]:
        raise ValueError(f"Invalid structure format: {format!r}")

    # Select parser based on format
    if format == "pdb":
        parser = PDBParser(PERMISSIVE=True)
    elif format == "mmcif":
        parser = MMCIFParser()

    # Load structure
    fname = str(fname)
    b = os.path.basename(fname)
    sname, _ = os.path.splitext(b)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure(sname, fname)

    return structure


def get_residue_data(chain, ignore_hetero=True):
    """
    Extract residue coordinates and sequence from a PDB chain.

    This function extracts the coordinates of the Cα atoms and the amino acid sequence
    from a given protein chain. It uses the coordinates of the Cα atom as the center
    of each residue.

    Parameters
    ----------
    chain : Bio.PDB.Chain.Chain
        Protein chain object from which to extract data.
    ignore_hetero : bool, optional
        Whether to ignore heteroatoms from the chain. Default is True.

    Returns
    -------
    coords : numpy.ndarray
        Array of shape (n, 3) containing the Cα atom coordinates for each residue,
        where n is the number of residues.
    seq : str
        String representing the amino acid sequence of the chain.

    Notes
    -----
    This function assumes that the `protein_letters_3to1` dictionary is available
    for converting 3-letter amino acid codes to 1-letter codes.

    The function skips residues that do not have a Cα atom.

    Examples
    --------
    >>> from Bio.PDB import PDBParser
    >>> structure = PDBParser().get_structure('1a1q', '1a1q.pdb')
    >>> chain = structure[0]['A']
    >>> coords, seq = get_residue_data(chain)
    >>> print(coords.shape)
    (141, 3)
    >>> print(len(seq))
    141
    """
    coords = []
    seq = []
    for residue in chain.get_residues():
        if residue.id[0] == " " and ignore_hetero:
            if "CA" in residue.child_dict:
                coords.append(residue.child_dict["CA"].coord)
                seq.append(protein_letters_3to1[residue.resname])

    return np.vstack(coords), "".join(seq)
