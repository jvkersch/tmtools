import numpy as np
from tmtools import tm_align

coords1 = np.array(
    [[1.2, 3.4, 1.5],
     [4.0, 2.8, 3.7],
     [1.2, 4.2, 4.3],
     [0.0, 1.0, 2.0]])
coords2 = np.array(
    [[2.3, 7.4, 1.5],
     [4.0, 2.9, -1.7],
     [1.2, 4.2, 4.3]])

seq1 = "AYLP"
seq2 = "ARN"

res = tm_align(coords1, coords2, seq1, seq2)
print(res.t)
print(res.u)
print(res.tm_norm_chain1)
print(res.tm_norm_chain2)
                   
