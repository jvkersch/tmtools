"""Helper tools to work with TM-Align."""

import numpy as np


def transform_structure(original_structure, tmalign_result):
    """
    Transform structure based on TM-Align result.

    This function creates a copy of the input structure and applies the rotation
    matrix and translation vector from a TM-align result to all atoms in the copied
    structure.

    Parameters
    ----------
    original_structure : Bio.PDB.Structure.Structure
        The original structure to transform.
    tmalign_result : TM_result
        An object containing the TM-align transformation data, with attributes:
        - u : array_like
            3x3 rotation matrix.
        - t : array_like
            Translation vector with 3 elements.

    Returns
    -------
    Bio.PDB.Structure.Structure
        A new structure object with the TM-align transformation applied. The ID of
        this structure is the original ID with "_tmalign" appended.

    Note
    ----
    To align atomic coordinates based on the TM-Align result, one can use the following code:
    ```python
    >>> res = tm_align(coords1, coords2, seq1, seq2)
    >>> aligned_coords1 = coords1 @ res.u.T + res.t
    ```
    """
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
