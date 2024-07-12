from Bio.PDB.Structure import Structure
import numpy as np


def apply_tmalign_transformation(original_structure, tmalign_result):
    # Check if inputs are correct
    assert isinstance(
        original_structure, Structure
    ), "Input must be a Bio.PDB.Structure.Structure object"
    assert hasattr(tmalign_result, "u") and hasattr(
        tmalign_result, "t"
    ), "Alignment data must have 'u' and 't' attributes"
    assert len(tmalign_result.u) == 3 and all(
        len(row) == 3 for row in tmalign_result.u
    ), "Rotation matrix must be 3x3"
    assert len(tmalign_result.t) == 3, "Translation vector must have 3 elements"

    # Create a copy of the structure to transform
    aligned_structure = original_structure.copy()
    aligned_structure.id = aligned_structure.id + "_tmalign"

    # Create rotation matrix and translation vector as numpy arrays
    rotation_matrix = np.array(tmalign_result.u, dtype=float)
    translation_vector = np.array(tmalign_result.t, dtype=float)

    # Iterate through all atoms in the structure applying transformation
    for model in aligned_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    original_coord = atom.coord
                    transformed_coord = (
                        np.dot(rotation_matrix, original_coord) + translation_vector
                    )
                    atom.set_coord(transformed_coord)

    return aligned_structure
