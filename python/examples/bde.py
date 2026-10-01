"""Bond dissociation energies, end to end: GFN2-xTB in the gas phase.

    python bde.py            parent conformer search, then optimization, Hessians, single points
    python bde.py --quick    skip the conformer search

The O-H and C-H bonds of methanol (one `scan` over every bond to hydrogen) and
the O-H of phenol. The geometries are hand-built and deliberately rough -- the
optimization stage is what makes them minima.

**What this asserts is identities and sanity, not agreement with experiment.**
GFN2 radicals are xTB's weak point -- mqc runs a doublet as an unpolarized
calculation with the unpaired electron set (see `mqc.bde`) -- so errors of
several kcal/mol are expected, and the experimental values are printed beside the
computed ones to show how many. What is asserted:

  * every energy is finite and the three methanol C-H bonds are one calculation
    (four bonds, three distinct fragments);
  * D0 is below D_e (zero-point energy is lost), every bond is bound, and
    dG < BDE(298) (a dissociation gains entropy);
  * the BDEs are ordered, and the table, the JSON and the kJ/mol column agree;
  * phenol's O-H is weaker than methanol's. Experimentally by about 18
    kcal/mol, because the phenoxyl radical delocalizes the spin over the ring --
    an effect of pi conjugation any semiempirical method carries, and the one
    ordering here that an unpolarized xTB doublet should get right.

Two of the four stages run the ``mqc`` executable (the library refuses the
``optimize`` and ``conformers`` drivers through a session). A build without
tblite or without the executable is reported and skipped rather than failed;
one without CREST runs without a conformer search and says so.

Run from a scratch directory: the Hessian and energy runs leave
``output_bde_*.json`` in the working directory and the executable stages work in
``bde_example_work/``.
"""

import math
import sys

import mqc
from mqc import bde

# -- hand-built geometries, Angstrom ----------------------------------------

#: C0 O1 H2 H3 H4 (on C) H5 (on O)
METHANOL = (
    ["C", "O", "H", "H", "H", "H"],
    [[-0.047, 0.664, 0.000], [-0.047, -0.758, 0.000], [-1.092, 0.969, 0.000],
     [0.439, 1.073, 0.890], [0.439, 1.073, -0.890], [0.858, -1.070, 0.000]],
)


def phenol():
    """C0-C5 ring (C0 bears the oxygen), O6, H7 on O, then H8-H12 on C1-C5."""
    symbols = ["C"] * 6 + ["O", "H"] + ["H"] * 5
    coords = [[1.395 * math.cos(math.radians(60 * k)), 1.395 * math.sin(math.radians(60 * k)), 0.0]
              for k in range(6)]
    coords.append([1.395 + 1.36, 0.0, 0.0])
    coords.append([1.395 + 1.36 + 0.96 * math.cos(math.radians(71)), 0.96 * math.sin(math.radians(71)), 0.0])
    coords += [[2.48 * math.cos(math.radians(60 * k)), 2.48 * math.sin(math.radians(60 * k)), 0.0]
               for k in range(1, 6)]
    return symbols, coords


#: kcal/mol, gas phase, 298 K (Luo, Comprehensive Handbook of Chemical Bond Energies).
EXPERIMENT = {"methanol O-H": 105.0, "methanol C-H": 96.0, "phenol O-H": 87.0}

NO_TBLITE = ("tblite", "not built", "build with")
NO_EXECUTABLE = ("none was found",)
NO_CREST = ("cannot sample conformers", "MQC_ENABLE_CREST")


def tblite_available():
    """One xTB single point that either returns or names the missing option."""
    system = mqc.System(symbols=["O", "H", "H"],
                        coords=[[0.0, 0.0, 0.1], [0.0, 0.77, -0.47], [0.0, -0.77, -0.47]])
    system.set_monomers([[0, 1, 2]])
    try:
        mqc.MBE(system, level=0, method="gfn2", verbosity="error").run(
            label="probe_bde_tblite", write_to_file=False)
    except mqc.MQCError as exc:
        return not any(hint in str(exc) for hint in NO_TBLITE)
    return True


def finite(x):
    return isinstance(x, float) and math.isfinite(x)


def run_scan(parent, bonds, protocol):
    """`bde.scan`, retried without the conformer search if the build has no CREST."""
    try:
        return bde.scan(parent, bonds, protocol, prefix="bde")
    except bde.BDEError as exc:
        if protocol.conformers is None or not any(hint in str(exc) for hint in NO_CREST):
            raise
        print("  this build has no CREST: running without a conformer search")
        protocol.conformers = None
        return bde.scan(parent, bonds, protocol, prefix="bde")


def check_rows(result):
    for row in result.bonds:
        for value in (row.d_e_kcal, row.d0_kcal, row.bde_kcal, row.dg_kcal):
            assert finite(value), f"{row.label}: {value}"
        assert row.d0_kcal < row.d_e_kcal, f"{row.label}: D0 {row.d0_kcal:.2f} not below D_e {row.d_e_kcal:.2f}"
        assert row.bde_kcal > 0.0, f"{row.label}: not bound, BDE {row.bde_kcal:.2f}"
        assert row.dg_kcal < row.bde_kcal, f"{row.label}: dG {row.dg_kcal:.2f} >= BDE {row.bde_kcal:.2f}"
        assert abs(row.bde_kj - 4.184 * row.bde_kcal) < 1e-9
    bdes = [r.bde_kcal for r in result.bonds]
    assert bdes == sorted(bdes), "rows are not sorted by BDE"
    for warning in result.warnings:
        print("  note:", warning)


def methanol(protocol):
    print("\nmethanol, every bond to hydrogen")
    result = run_scan(METHANOL, "X-H", protocol)
    check_rows(result)
    print(result.table())

    assert len(result.bonds) == 4, "methanol has four bonds to hydrogen"
    assert len(result.species) == 3, "CH2OH, CH3O and H: the three methyl C-H bonds are one fragment"
    if result.parent.n_imaginary:  # a soft methyl rotor can come back imaginary; reported, not fatal
        print(f"  note: the parent has {result.parent.n_imaginary} imaginary frequency(ies)")
    for sp in result.species.values():
        if sp.n_atoms > 1:
            assert sp.relaxation_kcal is not None and sp.relaxation_kcal > -0.5, sp.name

    oh = next(r for r in result.bonds if r.i == 1)
    ch = [r for r in result.bonds if r.i == 0]
    assert max(r.bde_kcal for r in ch) - min(r.bde_kcal for r in ch) < 1e-9, "equivalent C-H bonds differ"
    print(f"\n  O-H  {oh.bde_kcal:7.1f}   experiment {EXPERIMENT['methanol O-H']:.0f}")
    print(f"  C-H  {ch[0].bde_kcal:7.1f}   experiment {EXPERIMENT['methanol C-H']:.0f}")
    assert '"bonds"' in result.to_json()
    return oh


def phenol_oh(protocol):
    print("\nphenol O-H")
    result = run_scan(phenol(), [(6, 7)], protocol)
    check_rows(result)
    (row,) = result.bonds
    print(f"  O-H  {row.bde_kcal:7.1f}   experiment {EXPERIMENT['phenol O-H']:.0f}")
    print(f"  D_e {row.d_e_kcal:.1f}   D0 {row.d0_kcal:.1f}   BDE {row.bde_kcal:.1f}   dG {row.dg_kcal:.1f}  kcal/mol")
    return row


def main(argv):
    quick = "--quick" in argv
    with mqc.session():
        if not tblite_available():
            print("skipped: this build has no tblite")
            return 0
        try:
            bde.find_executable()
        except bde.BDEError as exc:
            if any(hint in str(exc) for hint in NO_EXECUTABLE):
                print("skipped: the optimize and conformers stages need the mqc executable "
                      "(set MQC_EXECUTABLE)")
                return 0
            raise
        protocol = bde.Protocol(workdir="bde_example_work", max_conformers=3,
                                conformers=None if quick else bde.Protocol().conformers)
        methanol_oh = methanol(protocol)
        phenol = phenol_oh(protocol)
        assert phenol.bde_kcal < methanol_oh.bde_kcal, (
            f"phenol O-H ({phenol.bde_kcal:.1f}) is not weaker than methanol O-H ({methanol_oh.bde_kcal:.1f})")
        print(f"\n  phenol O-H is {methanol_oh.bde_kcal - phenol.bde_kcal:.1f} kcal/mol weaker than "
              f"methanol's (experiment: {EXPERIMENT['methanol O-H'] - EXPERIMENT['phenol O-H']:.0f})")
    print("\nBDE workflow: identities hold")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
