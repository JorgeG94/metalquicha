"""Na+(H2O)n and F-(H2O)n come apart into the ion and n waters, against a real libmqc.

    python3 python/tests/test_ionic_monomers.py

Needs a build (`build/libmqc.so` or MQC_LIBRARY); skipped when there is none.
"""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

try:
    import mqc
except ImportError:
    mqc = None

WATERS = [
    ("O", (2.30, 0.00, 0.00)), ("H", (2.90, 0.76, 0.00)), ("H", (2.90, -0.76, 0.00)),
    ("O", (-2.30, 0.00, 0.00)), ("H", (-2.90, 0.76, 0.00)), ("H", (-2.90, -0.76, 0.00)),
    ("O", (0.00, 2.30, 0.00)), ("H", (0.76, 2.90, 0.00)), ("H", (-0.76, 2.90, 0.00)),
]


def ion_water(ion, n):
    atoms = [(ion, (0.0, 0.0, 0.0))] + WATERS[: 3 * n]
    lines = [str(len(atoms)), "%s (H2O)%d" % (ion, n)]
    lines += ["%s %.4f %.4f %.4f" % (el, *xyz) for el, xyz in atoms]
    return "\n".join(lines) + "\n"


def check(ion, charge):
    if mqc is None:
        print("skipped: libmqc not available")
        return
    for n in (1, 2, 3):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ion.xyz")
            with open(path, "w") as handle:
                handle.write(ion_water(ion, n))
            s = mqc.System.from_xyz(path, charge=charge)
            s.perceive_monomers(atom_charges={0: charge})
            assert len(s.monomers) == n + 1, s.monomers
            assert s.monomers[0] == [0], s.monomers
            assert s.missing_bonds() == 0


def test_sodium_alone():
    check("Na", 1)


def test_fluoride_alone():
    check("F", -1)


if __name__ == "__main__":
    test_sodium_alone()
    test_fluoride_alone()
    print("ok")
