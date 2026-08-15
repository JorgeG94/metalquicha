"""Every backend the Python interface can reach, on the same two molecules.

    python backends.py            run them all
    python backends.py xtb        run the ones whose name contains "xtb"

Four cases, because they are four different paths through the program and
three of them can break without the others noticing:

    standalone xTB      one molecule, tblite
    standalone HF       one molecule, libcint on the CPU
    fragmented xTB      MBE(2), tblite per fragment
    fragmented HF       MBE(2), libcint per fragment

The reference energies for the Hartree-Fock cases come from PySCF on the same
geometry and basis. The xTB ones come from this program, which is weaker --
they are there to catch a change, not to prove the physics -- and they are
marked as such below.

Run as a test rather than a demonstration: each case asserts, and the script
exits non-zero if any of them is wrong. An example nobody executes stops being
true without anyone finding out, which is the whole reason this is in CI.
"""

import sys

import mqc

#: Water, and two waters far enough apart to be separate monomers.
WATER = dict(
    symbols=["O", "H", "H"],
    coords=[[0.0, 0.0, 0.10077199], [0.0, 0.77250895, -0.46780200],
            [0.0, -0.77250895, -0.46780200]],
)
DIMER = dict(
    symbols=["O", "H", "H", "O", "H", "H"],
    coords=[[0.0, 0.0, 0.10077199], [0.0, 0.77250895, -0.46780200],
            [0.0, -0.77250895, -0.46780200],
            [4.0, 0.0, 0.10077199], [4.0, 0.77250895, -0.46780200],
            [4.0, -0.77250895, -0.46780200]],
)

#: PySCF RHF/sto-3g on the WATER geometry above. An external reference: if this
#: disagrees, one of the two codes is wrong, which is the point.
HF_WATER_STO3G = -74.962005687948

#: PySCF density-fitted RHF/sto-3g, def2-universal-jkfit, on WATER. Fitting is
#: an approximation, so this differs from the exact number by 8.6e-5 -- which
#: is the point: a wired-up fitting path must reproduce the fitted reference,
#: not the exact one.
HF_WATER_STO3G_DF = -74.962092123400

#: PySCF RHF/sto-3g on the whole DIMER. Not a monomer sum -- see
#: fragmented_hf() for why that is the reference an MBE(2) must reproduce.
HF_DIMER_STO3G = -149.922856255009

#: This program's own answer for one water. It checks that nothing changed, not
#: that anything is right: tblite has no second implementation here to check
#: against, and a tighter tolerance would be false precision.
XTB_WATER = -5.070544

#: PySCF RHF/sto-3g nuclear gradient norm on the WATER geometry, in
#: Hartree/Bohr, from the same basis JSON this program reads. The CPU backend
#: refused gradients when this file was written and now computes them, so what
#: was a check that it refused cleanly is now a check that it is right --
#: the property the refusal was standing in for.
HF_WATER_STO3G_GRAD_NORM = 0.086151389909

TOL_HF = 1.0e-6      # against PySCF
TOL_GRAD = 1.0e-6    # against PySCF, on the gradient norm
TOL_XTB = 1.0e-4     # against ourselves; loose on purpose
TOL_CLOSURE = 1.0e-8 # MBE(2) against the supermolecule; exact, so tight


def standalone_xtb():
    """One molecule through tblite.

    The whole molecule is declared as one monomer. The validator wants a
    partition whether or not the expansion will use one, and auto_monomers()
    refuses a single covalently connected molecule -- rightly, since where to
    cut it is a chemical choice and not something to guess. At level 0 nothing
    is cut, so saying "all of it, once" is the honest way to ask.
    """
    system = mqc.System(**WATER)
    system.set_monomers([[0, 1, 2]])
    return mqc.MBE(system, level=0, method="gfn2").run(
        label="py_standalone_xtb", write_to_file=False).energy


def standalone_hf():
    """One molecule through libcint on the CPU."""
    system = mqc.System(**WATER)
    system.set_monomers([[0, 1, 2]])
    return mqc.MBE(system, level=0, method="hf", basis="sto-3g").run(
        label="py_standalone_hf", write_to_file=False).energy


def standalone_hf_df():
    """One molecule, density-fitted J and K on the CPU."""
    system = mqc.System(**WATER)
    system.set_monomers([[0, 1, 2]])
    return mqc.MBE(system, level=0, method="hf", basis="sto-3g",
                   aux_basis="def2-universal-jkfit", density_fitting=True).run(
        label="py_standalone_hf_df", write_to_file=False).energy


def gradient_support():
    """Both backends produce a gradient, and the CPU one produces the right one.

    This used to assert that libcint *refused*, which was correct when no CPU
    derivative existed: a clean refusal beats a wrong number. Now that the
    gradients are implemented the refusal is gone, and checking merely that
    something came back would drop the guard without replacing it. So the
    number itself is checked, against PySCF on the same geometry and the same
    basis data -- which is what the refusal was protecting all along.

    Returns (cpu_norm, cpu_ok, xtb_ok). `cpu_norm` is None if the run produced
    no gradient at all, which is a failure now rather than the expected result.
    """
    system = mqc.System(**WATER)
    system.set_monomers([[0, 1, 2]])

    # write_to_file=True because the norm travels in the output document, the
    # same route `gap_ev` takes.
    cpu_norm = None
    try:
        result = mqc.MBE(system, level=0, method="hf", basis="sto-3g",
                         driver="gradient").run(label="py_grad_hf", write_to_file=True)
        cpu_norm = result.gradient_norm
    except mqc.MQCError:
        cpu_norm = None
    cpu_ok = (cpu_norm is not None
              and abs(cpu_norm - HF_WATER_STO3G_GRAD_NORM) <= TOL_GRAD)

    xtb_ok = False
    try:
        mqc.MBE(system, level=0, method="gfn2", driver="gradient").run(
            label="py_grad_xtb", write_to_file=False)
        xtb_ok = True
    except mqc.MQCError:
        xtb_ok = False

    return cpu_norm, cpu_ok, xtb_ok


def fragmented_xtb():
    """MBE(2) over two monomers against the same dimer computed whole.

    With exactly two monomers the expansion closes:

        E = E_1 + E_2 + (E_12 - E_1 - E_2) = E_12

    so MBE(2) is not an approximation here, it is the supermolecular energy by
    a longer route. That makes this a real test of the fragmentation machinery
    -- the partition, the n-body deltas, the assembly -- rather than a recorded
    constant, and it needs no external reference because both sides come from
    the same method. Anything but agreement to rounding is a bug in the
    expansion.
    """
    system = mqc.System(**DIMER)
    system.auto_monomers()
    fragmented = mqc.MBE(system, level=2, method="gfn2").run(
        label="py_frag_xtb", write_to_file=False).energy

    whole = mqc.System(**DIMER)
    whole.set_monomers([[0, 1, 2, 3, 4, 5]])
    supermolecule = mqc.MBE(whole, level=0, method="gfn2").run(
        label="py_whole_xtb", write_to_file=False).energy
    return fragmented, supermolecule


def fragmented_hf():
    """MBE(2) over two monomers, libcint per fragment, against PySCF.

    Same closure as fragmented_xtb: two monomers means the expansion reproduces
    the dimer exactly. So the reference is PySCF on the whole six-atom dimer --
    external, and it checks the fragmentation and the CPU SCF together.

    A monomer sum would be the wrong reference and is what this originally
    used: it omits the interaction, which here is 1.2 millihartree, and the
    test failed until the reference was corrected rather than the code.
    """
    system = mqc.System(**DIMER)
    system.auto_monomers()
    return mqc.MBE(system, level=2, method="hf", basis="sto-3g").run(
        label="py_frag_hf", write_to_file=False).energy


def main(argv):
    pattern = argv[1] if len(argv) > 1 else ""
    failures = []

    with mqc.session():
        cases = [
            ("standalone xtb", standalone_xtb, XTB_WATER, TOL_XTB, "ours"),
            ("standalone hf", standalone_hf, HF_WATER_STO3G, TOL_HF, "PySCF"),
            ("standalone hf/df", standalone_hf_df, HF_WATER_STO3G_DF, TOL_HF, "PySCF fitted"),
            ("fragmented xtb", fragmented_xtb, None, TOL_CLOSURE, "same dimer, whole"),
            ("fragmented hf", fragmented_hf, HF_DIMER_STO3G, TOL_HF, "PySCF dimer"),
        ]
        for name, fn, expected, tol, source in cases:
            if pattern and pattern not in name:
                continue
            energy = fn()
            if expected is None:
                # The case returned (fragmented, supermolecule): its reference
                # is the second, computed in this same run.
                energy, expected = energy
            delta = abs(energy - expected)
            ok = delta <= tol
            print(f"  {name:<18} {energy:>18.9f}   ref {expected:>18.9f} "
                  f"({source})   diff {delta:.2e}   {'ok' if ok else 'FAIL'}")
            if not ok:
                failures.append(f"{name}: {energy} vs {expected}, diff {delta:.2e} > {tol:.1e}")

        cpu_norm, cpu_ok, xtb_ok = gradient_support()
        shown = f"{cpu_norm:.9f}" if cpu_norm is not None else "none"
        delta = (abs(cpu_norm - HF_WATER_STO3G_GRAD_NORM)
                 if cpu_norm is not None else float("inf"))
        print(f"  {'gradient hf':<18} {shown:>18}   ref "
              f"{HF_WATER_STO3G_GRAD_NORM:>18.9f} (PySCF)   diff {delta:.2e}   "
              f"{'ok' if cpu_ok else 'FAIL'}")
        print(f"  {'gradient xtb':<18} {'provided' if xtb_ok else 'MISSING':>18}")
        if cpu_norm is None:
            failures.append("the CPU backend returned no gradient at all")
        elif not cpu_ok:
            failures.append(f"CPU gradient norm {cpu_norm} vs "
                            f"{HF_WATER_STO3G_GRAD_NORM}, diff {delta:.2e} > {TOL_GRAD:.1e}")
        if not xtb_ok:
            failures.append("tblite must still provide gradients")

    if failures:
        print("\nFAILED:")
        for f in failures:
            print("  " + f)
        return 1
    print("\nall backends agree with their references")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
