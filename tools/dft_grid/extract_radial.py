"""Extract radial-grid tables from PySCF and a reference mesh to compare against."""
import json, sys
import numpy as np
from pyscf.dft import radi

xi = list(radi._treutler_ahlrichs_xi)
bragg = list(radi.BRAGG_RADII)

ref = {}
for z in (1, 6, 8, 26, 79):
    for n in (25, 75):
        r, dr = radi.treutler_ahlrichs(n, z)
        ref[f"{z}_{n}"] = {"r": r.tolist(), "dr": dr.tolist()}

json.dump({"xi": xi, "bragg": bragg, "reference": ref}, sys.stdout)
