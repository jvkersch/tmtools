import numpy as np


def apply_tmalign_transformation(structure, tmalign_result):
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
