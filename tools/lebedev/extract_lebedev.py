"""Record the Lebedev-Laikov orbit parameters PySCF uses, by instrumenting it.

Parsing 5000 lines of source would be fragile. Instead we monkeypatch
SphGenOh, call every MakeAngularGrid_N, and record the (code, a, b, v) tuples
it is actually handed. The output is the published Lebedev-Laikov parameter
data -- numbers, not anyone's code -- which is what the Fortran side needs.

Emits JSON on stdout: {order: [[code, a, b, v], ...]}
"""
import json
import sys

import pyscf.dft.LebedevGrid as L

calls = []
_real = L.SphGenOh


def spy(code, a, b, v):
    calls.append((int(code), float(a), float(b), float(v)))
    return _real(code, a, b, v)


L.SphGenOh = spy

tables = {}
reference = {}
for order in L.LEBEDEV_NGRID:
    order = int(order)  # LEBEDEV_NGRID is a numpy array; np.int32 is not a JSON key
    fn = getattr(L, f"MakeAngularGrid_{order}", None)
    if fn is None:
        print(f"# no generator for order {order}", file=sys.stderr)
        continue
    calls.clear()
    grid = fn()
    tables[order] = list(calls)
    # Keep the points too, as the comparison target for the Fortran gen_oh.
    reference[order] = grid.tolist()

json.dump({"params": tables, "reference": reference}, sys.stdout)
