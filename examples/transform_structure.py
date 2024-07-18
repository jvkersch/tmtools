from tmtools.io import get_residue_data, get_structure
from tmtools.helpers import transform_structure
from tmtools.testing import get_pdb_path
from tmtools import tm_align

from Bio.PDB.PDBIO import PDBIO

# Load structures
guide_struct = get_structure(get_pdb_path("2gtl"))
mobile_struct = get_structure(get_pdb_path("7ok9"))

# Get relevant data from first chain
guide_coords, guide_seq = get_residue_data(next(guide_struct.get_chains()))
mobile_coords, mobile_seq = get_residue_data(next(mobile_struct.get_chains()))

# Run TMalign
res = tm_align(mobile_coords, guide_coords, mobile_seq, guide_seq)
aligned_mobile_struct = transform_structure(mobile_struct, res)

# Save aligned structure
io = PDBIO()
io.set_structure(aligned_mobile_struct)
io.save("7ok9_tmalign.pdb")
