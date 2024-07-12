import numpy as np


def apply_tmalign_transformation(structure, tmalign_result):
    # Create a copy of the structure to transform
    transformed_structure = structure.copy()
    transformed_structure.id = transformed_structure.id + "_tmalign"

    # Create rotation matrix and translation vector as numpy arrays
    rotation = np.array(tmalign_result.u, dtype=float)
    translation = np.array(tmalign_result.t, dtype=float)

    # Iterate through all atoms in the structure applying transformation
    for model in transformed_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    old_coord = atom.coord
                    new_coord = np.dot(rotation, old_coord) + translation
                    atom.set_coord(new_coord)

    return transformed_structure
