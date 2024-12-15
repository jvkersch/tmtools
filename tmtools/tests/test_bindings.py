import contextlib
import os
import sys
import threading
import unittest

import numpy.testing as nptest
import numpy as np

from ..io import get_residue_data, get_structure
from ..testing import get_pdb_path

from tmtools import tm_align, print_version


def _coords_from_pdb(sname):
    pdb = get_pdb_path(sname)
    s = get_structure(pdb)
    c = next(s.get_chains())
    return get_residue_data(c)


@contextlib.contextmanager
def _capture_stdout():
    # Based on https://stackoverflow.com/a/24277852
    stdout_fileno = sys.stdout.fileno()
    stdout_save = os.dup(stdout_fileno)
    stdout_pipe = os.pipe()
    os.dup2(stdout_pipe[1], stdout_fileno)
    os.close(stdout_pipe[1])

    captured_stdout = ''
    def drain_pipe():
        nonlocal captured_stdout
        while True:
            data = os.read(stdout_pipe[0], 1024).decode()
            if not data:
                break
            captured_stdout += data

    t = threading.Thread(target=drain_pipe)
    t.start()

    result = _StdoutData()

    try:
        yield result
    finally:
        os.close(stdout_fileno)
        t.join()

        os.close(stdout_pipe[0])
        os.dup2(stdout_save, stdout_fileno)
        os.close(stdout_save)

        result.stdout = captured_stdout


class _StdoutData:
    pass


class TestBindings(unittest.TestCase):
    def test_call_identical(self):
        # Given
        coords, seq = _coords_from_pdb("2gtl")

        # When
        res = tm_align(coords, coords, seq, seq)

        # Then
        nptest.assert_array_almost_equal(res.t, np.zeros(3))
        nptest.assert_array_almost_equal(res.u, np.eye(3))
        nptest.assert_almost_equal(res.rmsd, 0.0, decimal=4)
        self.assertEqual(res.seqxA, seq)
        self.assertEqual(res.seqyA, seq)
        self.assertEqual(res.seqM, ":" * len(seq))
        self.assertEqual(res.seqxA, res.seqyA)

    def test_call_segments(self):
        # Given
        coords, seq = _coords_from_pdb("2gtl")
        coords1 = coords[0:40]
        seq1 = seq[0:40]
        coords2 = coords[10:50]
        seq2 = seq[10:50]

        # When
        res = tm_align(coords1, coords2, seq1, seq2)

        # Then
        nptest.assert_array_almost_equal(res.t, np.zeros(3))
        nptest.assert_array_almost_equal(res.u, np.eye(3))
        nptest.assert_almost_equal(res.rmsd, 0.0, decimal=4)
        self.assertEqual(res.seqxA, f"{seq1}----------")
        self.assertEqual(res.seqyA, f"----------{seq2}")
        self.assertEqual(res.seqM, " " * 10 + ":" * 30 + " " * 10)


    def test_call_different(self):
        # Given
        coords1, seq1 = _coords_from_pdb("2gtl")
        coords2, seq2 = _coords_from_pdb("7ok9")

        # Verified by running TMalign manually
        t_expected = np.array([10.42888594, 42.96954856, 74.43889102])
        u_expected = np.array(
            [
                [-0.93809129, -0.30032785, 0.17259179],
                [-0.21695917, 0.12102283, -0.96864968],
                [0.27002492, -0.94612719, -0.17868934],
            ]
        )
        tm_norm2 = 0.15158
        tm_norm1 = 0.38759
        rmsd = 5.37

        # When
        res = tm_align(coords1, coords2, seq1, seq2)

        # Then
        nptest.assert_array_almost_equal(res.t, t_expected)
        nptest.assert_array_almost_equal(res.u, u_expected)
        self.assertAlmostEqual(res.tm_norm_chain1, tm_norm1, places=4)
        self.assertAlmostEqual(res.tm_norm_chain2, tm_norm2, places=4)
        self.assertAlmostEqual(res.rmsd, rmsd, places=1)

    def test_call_error(self):
        # Given
        coords1 = np.array([[0, 0, 0],
                            [1, 1, 1],
                            [2, 2, 2],
                            [3, 3, 3]])
        coords2 = np.array(
            [[0, 0, 0],
             [1, 1, 1]])
        seq1 = "AAAA"
        seq2 = "CC"

        # When/then
        with self.assertRaises(RuntimeError):
            tm_align(coords1, coords2, seq1, seq2)

    def test_print_version(self):
        # When
        with _capture_stdout() as captured:
            print_version()

        # Then
        self.assertIn("Version 20210224", captured.stdout)

