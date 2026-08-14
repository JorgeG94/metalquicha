"""Atomic partial charges, and what they are for.

    python charges.py

Two schemes are exposed, and the interesting part is where they disagree.

**Mulliken** divides the density by which basis function carries it. Cheap --
one trace over `D S` -- and it inherits every arbitrariness of the basis,
because a function's atom label is what assigns its density and a diffuse
function centred on hydrogen reaches well over the oxygen while still counting
as hydrogen's.

**CHELPG** asks a physical question instead: what point charges reproduce the
molecule's own electrostatic potential on a shell of points outside its van der
Waals surface. Costs an SCF plus an ESP evaluation on a few thousand points,
and the answer barely moves when the basis changes.

The last case here is the one that decides which to use, and it is a
measurement rather than an assertion of taste: the same water in two bases.

Run as a test rather than a demonstration: each case asserts, and the script
exits non-zero if any of them is wrong.
"""

import sys

import mqc

#: Water, Angstrom. The C API converts; internals are Bohr.
WATER = dict(
    symbols=["O", "H", "H"],
    coords=[[0.0, 0.0, 0.0], [0.0, -0.7572, 0.5865], [0.0, 0.7572, 0.5865]],
)

def system(spec, charge=0):
    s = mqc.System()
    s.set_geometry(spec["symbols"], spec["coords"], charge=charge)
    return s


def show(label, symbols, q):
    body = "  ".join(f"{s}{i} {v:+.4f}" for i, (s, v) in enumerate(zip(symbols, q)))
    print(f"  {label:10s} {body}")


def case_both_schemes():
    """Both schemes on water, and the properties any scheme has to have."""
    print("water, 6-31G")
    for scheme in ("mulliken", "chelpg"):
        s = system(WATER).compute_charges(scheme=scheme, basis="6-31g")
        q = s.charges()
        show(scheme, WATER["symbols"], q)

        # They sum to the molecular charge. For CHELPG that is the constraint
        # the fit is solved under; for Mulliken it is an identity of the trace.
        assert abs(sum(q)) < 1e-8, f"{scheme} charges sum to {sum(q)}"
        # Oxygen is the electronegative one and both schemes must say so.
        assert q[0] < 0 < q[1], f"{scheme} got the polarity backwards"
        # The two hydrogens are related by symmetry. CHELPG's grid is a cubic
        # lattice that does not share the molecule's symmetry, so it gets a
        # loose tolerance; Mulliken has no such excuse.
        tol = 5e-3 if scheme == "chelpg" else 1e-9
        assert abs(q[1] - q[2]) < tol, f"{scheme} split equivalent hydrogens"

        # The single-atom accessor agrees with the vector.
        assert abs(s.charge_on(0) - q[0]) < 1e-14
        # And the handle remembers which question it answered.
        assert s.charge_scheme == scheme, f"handle says {s.charge_scheme!r}"


def case_total_charge():
    """A charged system: the fit is constrained to the charge, not to zero."""
    print("hydroxide, 6-31G")
    s = mqc.System()
    s.set_geometry(["O", "H"], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.97]], charge=-1)
    s.compute_charges(scheme="chelpg", basis="6-31g")
    q = s.charges()
    show("chelpg", ["O", "H"], q)
    assert abs(sum(q) + 1.0) < 1e-8, f"charges sum to {sum(q)}, not -1"


def case_basis_sensitivity():
    """The measurement that decides which scheme to feed downstream."""
    print("water: how much does the basis move the answer?")
    moved = {}
    for scheme in ("mulliken", "chelpg"):
        small = system(WATER).compute_charges(scheme=scheme, basis="6-31g").charges()
        big = system(WATER).compute_charges(scheme=scheme, basis="aug-cc-pvdz").charges()
        moved[scheme] = max(abs(b - a) for a, b in zip(small, big))
        print(
            f"  {scheme:10s} O {small[0]:+.4f} -> {big[0]:+.4f}"
            f"    largest move {moved[scheme]:.4f}"
        )

    # The whole reason CHELPG is here. Not a close call in practice -- Mulliken
    # moves by about half an electron on the oxygen -- so the assertion is
    # generous and would still catch the two being swapped.
    assert moved["chelpg"] < moved["mulliken"], (
        "CHELPG moved at least as much as Mulliken between bases, which is the "
        "opposite of why it is exposed"
    )
    print(f"  -> CHELPG is {moved['mulliken'] / moved['chelpg']:.1f}x the more stable")


def case_errors():
    """The failure paths report rather than returning something plausible."""
    print("failure paths")
    s = system(WATER)
    assert not s.has_charges
    assert s.charge_scheme == "", "an uncomputed handle claims a scheme"

    for fn, what in (
        (lambda: s.charges(), "reading before computing"),
        (lambda: s.charge_on(0), "one atom before computing"),
    ):
        try:
            fn()
        except mqc.MQCError:
            print(f"  ok  {what} raises")
        else:
            raise AssertionError(f"{what} did not raise")

    s.compute_charges(basis="6-31g")
    assert s.has_charges

    try:
        s.charge_on(99)
    except mqc.MQCError:
        print("  ok  an index past the end raises")
    else:
        raise AssertionError("an out-of-range index did not raise")

    try:
        s.compute_charges(scheme="lowdin")
    except mqc.MQCError:
        print("  ok  an unimplemented scheme raises rather than falling back")
    else:
        raise AssertionError("an unknown scheme did not raise")

    # Closed-shell RHF only. A radical has to be refused rather than quietly
    # paired up, because a caller fragmenting one needs to know.
    r = mqc.System()
    r.set_geometry(["H"], [[0.0, 0.0, 0.0]])
    try:
        r.compute_charges(basis="6-31g")
    except mqc.MQCError:
        print("  ok  an odd electron count raises")
    else:
        raise AssertionError("an open-shell system did not raise")


CASES = {
    "schemes": case_both_schemes,
    "total-charge": case_total_charge,
    "basis": case_basis_sensitivity,
    "errors": case_errors,
}


def main(argv):
    wanted = argv[1] if len(argv) > 1 else ""
    ran = 0
    for name, fn in CASES.items():
        if wanted and wanted not in name:
            continue
        fn()
        ran += 1
    if ran == 0:
        print(f"no case matching {wanted!r}; have: {', '.join(CASES)}")
        return 1
    print(f"\n{ran} case(s) passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
