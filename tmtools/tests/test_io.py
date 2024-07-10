import unittest

import numpy.testing as nptest

from ..io import get_residue_data, get_structure
from ..testing import get_pdb_path, get_mmcif_path


class TestIO(unittest.TestCase):
    def test_get_structure(self):
        # Test valid PDB file loading
        pdb = get_pdb_path("2gtl")
        pdb_structure = get_structure(pdb)
        self.assertEqual(pdb_structure.id, "2gtl")

        # Test valid MMCIF file loading
        mmcif = get_mmcif_path("2gtl")
        mmcif_structure = get_structure(mmcif, format="mmcif")
        self.assertEqual(mmcif_structure.id, "2gtl")

        # Test if PDB and MMCIF structures are equal
        self.assertEqual(pdb_structure, mmcif_structure)

    def test_get_residue_coordinates(self):
        # Given
        pdb = get_pdb_path("2gtl")
        structure = get_structure(pdb)
        chain0 = next(structure.get_chains())

        # When
        coords, seq = get_residue_data(chain0)

        # Then
        self.assertEqual(coords.shape, (147, 3))
        nptest.assert_array_almost_equal(coords[0], [14.58, 114.133, 44.707], decimal=5)
        nptest.assert_array_almost_equal(
            coords[-1], [18.708, 134.427, 44.33], decimal=5
        )

        self.assertEqual(len(seq), coords.shape[0])
        self.assertEqual(seq[:5], "DCCSY")
        self.assertEqual(seq[-5:], "AKDLP")
