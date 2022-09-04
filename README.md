TM-Tools
========

Python bindings for the TM-align algorithm and code [developed by Zhang et
al](https://zhanggroup.org/TM-align/) for protein structure comparison.


Installation
------------

You can install the released version of the package directly from PyPI by
running
```console
    pip install tmtools
```
Pre-built wheels are available for Linux, macOS, and Windows, for Python 3.6
and up.

To build the package from scratch, e.g. because you want to contribute to it,
clone this repository, and then from the root of the repository, run
```console
    pip install -e . -v
```
This requires a C++ compiler to be installed with support for C++ 14.

Usage
-----

The function `tmtools.tm_align` takes two NumPy arrays with coordinates for the
residues (with shape `(N, 3)`) and two sequences of peptide codes, performs the
alignment, and returns the optimal rotation matrix and translation, along with
the TM score:
```python
>>> import numpy as np
>>> from tmtools import tm_align
>>>
>>> coords1 = np.array(
...     [[1.2, 3.4, 1.5],
...      [4.0, 2.8, 3.7],
...      [1.2, 4.2, 4.3],
...      [0.0, 1.0, 2.0]])
>>> coords2 = np.array(
...     [[2.3, 7.4, 1.5],
...      [4.0, 2.9, -1.7],
...      [1.2, 4.2, 4.3]])
>>>
>>> seq1 = "AYLP"
>>> seq2 = "ARN"
>>>
>>> res = tm_align(coords1, coords2, seq1, seq2)
>>> res.t
array([ 2.94676159,  5.55265245, -1.75151383])
>>> res.u
array([[ 0.40393231,  0.04161396, -0.91384187],
       [-0.59535733,  0.77040999, -0.22807475],
       [ 0.69454181,  0.63618922,  0.33596866]])
>>> res.tm_norm_chain1
0.3105833326322145
>>> res.tm_norm_chain2
0.414111110176286
```

If you already have some PDB files, you can use the functions from `tmalign.io`
to retrieve the coordinate and sequence data. These functions rely on
`BioPython`, which is not installed by default to keep dependencies
lightweight. To use them, you have to install `BioPython` first (`pip install
biopython`). Then run:

```python
>>> from tmtools.io import get_structure, get_residue_data
>>> from tmtools.testing import get_pdb_path
>>> s = get_structure(get_pdb_path("2gtl"))
>>> s
<Structure id=2gtl>
>>> chain = next(s.get_chains())
>>> coords, seq = get_residue_data(chain)
>>> seq
'DCCSYEDRREIRHIWDDVWSSSFTDRRVAIVRAVFDDLFKHYPTSKALFERVKIDEPESGEFKSHLVRVANGLKLLINLLDDTLVLQSHLGHLADQHIQRKGVTKEYFRGIGEAFARVLPQVLSCFNVDAWNRCFHRLVARIAKDLP'
>>> coords.shape
(147, 3)
```

Credits
-------

This package arose out of a personal desire to better understand both the
TM-score algorithm and the
[pybind11](https://pybind11.readthedocs.io/en/stable/index.html) library to
interface with C++ code. At this point in time it contains no original research
code.

If you use the package for research, you should cite the [original TM-score
papers](https://zhanggroup.org/TM-score/):

- Y. Zhang, J. Skolnick, _Scoring function for automated assessment of protein
  structure template quality_, Proteins, 57: 702-710 (2004).
- J. Xu, Y. Zhang, How significant is a protein structure similarity with
  TM-score=0.5? Bioinformatics, 26, 889-895 (2010).

License
-------

The original TM-align software (version 20210224, released under the MIT
license) is bundled with this repository (`src/extern/TMalign.cpp`). Some small
tweaks had to be made to compile the code on macOS and to embed it as a
library. This modifications are also released under the MIT license.

The rest of the codebase is released under the GPL v3 license.
