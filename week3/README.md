1) distutils removed in Python 3.12+
Easiest: pin CI to Python 3.11 (see above), or


The handout suggested loading distances.txt via the Python bridge:

```
from python import numpy as pnp
import numpy.pybridge
distances: np.ndarray[int,2] = pnp.loadtxt("tests/sequence/data/distances.txt", dtype=pnp.int64)
# Now you can use distances
tree = upgma(distances)
```

This relies on CODON_PYTHON and numpy.pybridge. On my setup, that destabilized the environment (Python shared-lib conflicts). I tried to fight it but I just messed everything up so it was easier pure-Python file loader. 

```
def get_distances() -> ndarray:
    rows: List[List[int]] = []
    with open("../../data/sequence/distances.txt", "r") as fh:
        for ln in fh:
            s = ln.strip()
            if not s:
                continue
            rows.append([int(tok) for tok in s.split()])

    # Build an int32 numpy array with the parsed shape
    n = len(rows)
    m = len(rows[0]) if n else 0
    mat = zeros((n, m), dtype=int32)

    for r in range(n):
        vals = rows[r]
        for c in range(m):
            mat[r, c] = int32(vals[c])

    return mat
```
