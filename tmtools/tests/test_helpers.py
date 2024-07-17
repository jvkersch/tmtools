import unittest

import numpy.testing as nptest

from ..io import get_residue_data, get_structure
from ..helpers import transform_structure
from ..testing import get_pdb_path

from tmtools import tm_align


def _coords_from_pdb(sname):
    pdb = get_pdb_path(sname)
    s = get_structure(pdb)
    c = next(s.get_chains())
    return get_residue_data(c)


def _structure_from_pdb(sname):
    pdb = get_pdb_path(sname)
    s = get_structure(pdb)
    return s


class TestHelpers(unittest.TestCase):
    def test_call_different(self):
        # Given
        coords1, seq1 = _coords_from_pdb("2gtl")
        coords2, seq2 = _coords_from_pdb("7ok9")
        struct2 = _structure_from_pdb("7ok9")

        # When
        res = tm_align(coords1, coords2, seq1, seq2)
        aligned_coords2 = coords2 @ res.u.T + res.t
        aligned_struct2 = transform_structure(struct2, res)

        # Then
        nptest.assert_array_almost_equal(
            aligned_coords2, get_residue_data(aligned_struct2[0]["A"])[0]
        )
