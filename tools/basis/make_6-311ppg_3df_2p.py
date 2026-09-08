#!/usr/bin/env python3
"""6-311++G(3df,2p), the EFP-recommended basis, from the BSE sets we ship.

The Basis Set Exchange has no 6-311++G(3df,2p) entry. It is 6-311++G(3df,3pd)
on every atom heavier than hydrogen, and on hydrogen the same s functions
(three contracted, two valence, one diffuse) with the Pople split-polarization
pair of p exponents, 1.5 and 0.375 (the standard 0.75 times two and half),
and no d. That is what GAMESS builds for NDFUNC=3 NFFUNC=1 NPFUNC=2 with
DIFFSP and DIFFS on, and what PySCF ships under the same name.

    python tools/basis/make_6-311ppg_3df_2p.py

writes basis_sets/pople/6-311++g(3df,2p).json from
basis_sets/6-311++g(3df,3pd).json (the BSE bundle, unpacked by CMake).
"""
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "basis_sets" / "6-311++g(3df,3pd).json"
DST = REPO / "basis_sets" / "pople" / "6-311++g(3df,2p).json"
H_P_EXPONENTS = ["0.1500000000E+01", "0.3750000000E+00"]

src = json.loads(SRC.read_text())
out = dict(src)
out["name"] = "6-311++G(3df,2p)"
out["description"] = ("6-311++G(3df,2p): 6-311++G(3df,3pd) on heavy atoms, and on hydrogen "
                      "6-311++G with the (2p) polarization pair 1.5 and 0.375 and no d. "
                      "Assembled from the BSE 6-311++G(3df,3pd) entry by "
                      "tools/basis/make_6-311ppg_3df_2p.py; the EFP-recommended basis.")
h = json.loads(json.dumps(src["elements"]["1"]))
s_shells = [sh for sh in h["electron_shells"] if sh["angular_momentum"] == [0]]
p_template = next(sh for sh in h["electron_shells"] if sh["angular_momentum"] == [1])
p_shells = []
for e in H_P_EXPONENTS:
    sh = json.loads(json.dumps(p_template))
    sh["exponents"] = [e]
    sh["coefficients"] = [["0.1000000000E+01"]]
    p_shells.append(sh)
# The BSE order for H: contracted s, two valence s, polarization, then the
# diffuse s last. Keep it: the diffuse s stays where the source has it.
diffuse = [sh for sh in s_shells if float(sh["exponents"][0]) < 0.05]
core = [sh for sh in s_shells if sh not in diffuse]
h["electron_shells"] = core + p_shells + diffuse
out["elements"] = dict(src["elements"])
out["elements"]["1"] = h
DST.parent.mkdir(parents=True, exist_ok=True)
DST.write_text(json.dumps(out, indent=1) + "\n")
print(DST, "hydrogen shells:", [(sh["angular_momentum"], sh["exponents"]) for sh in h["electron_shells"]])
