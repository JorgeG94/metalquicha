"""Every method the backends implement, run through the Python interface.

    python methods.py              run them all
    python methods.py cc           run the ones whose name contains "cc"

`backends.py` covers the four *paths* through the program -- standalone and
fragmented, tblite and libcint. This covers the *methods*: what a deck can ask
for, asked for from Python instead. The two are separate scripts because they
fail for different reasons. A break in backends.py means the plumbing moved; a
break here means a method a deck can run is no longer reachable, or is reachable
and returns the wrong number.

The references are PySCF on the same geometry, fed this repository's own basis
JSON through `tools/cpu_validation/gen_cpu_validation.py`. That last part
matters: PySCF's internal basis tables differ from the Basis Set Exchange data
shipped here in the eighth decimal on some sets, which looks exactly like a bug
in whichever code you are checking.

Where a reference would be less informative than an identity, an identity is
used instead -- RI against conventional, spin-adapted against spin-orbital,
CASSCF against CASCI on the same space. Those hold exactly or to a stated
approximation, and they check the two routes against each other rather than
both against a third code.

Runs as a test: every case asserts, and the script exits non-zero if any is
wrong. A method the build cannot provide -- no libxc, no libcint -- is reported
as a skip rather than a failure, because a reduced build is a build, not a
regression.
"""

import os
import sys

import mqc

#: Water. The geometry backends.py uses, so the two scripts' references are
#: comparable and the Hartree-Fock number below is literally the same one.
WATER = dict(
    symbols=["O", "H", "H"],
    coords=[[0.0, 0.0, 0.10077199], [0.0, 0.77250895, -0.46780200],
            [0.0, -0.77250895, -0.46780200]],
)

#: The hydroxyl radical: a doublet, which is how the unrestricted path is
#: reached. Nothing asks for UHF -- the multiplicity selects it.
HYDROXYL = dict(symbols=["O", "H"], coords=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.96966]])

#: Two waters 3 A apart, for the methods that return an interaction rather than
#: an energy. SAPT0 and EFP2 both read the monomers as the interacting pieces,
#: so the partition is the physics here and not bookkeeping.
DIMER = dict(
    symbols=["O", "H", "H", "O", "H", "H"],
    coords=[[0.0, 0.0, 0.10077199], [0.0, 0.77250895, -0.46780200],
            [0.0, -0.77250895, -0.46780200],
            [3.0, 0.0, 0.10077199], [3.0, 0.77250895, -0.46780200],
            [3.0, -0.77250895, -0.46780200]],
)

#: The same dimer pulled to 30 A, where any interaction energy must have
#: decayed to nothing. Cheaper than a reference and harder to satisfy by
#: accident.
FAR_DIMER = dict(
    symbols=DIMER["symbols"],
    coords=[list(xyz) for xyz in DIMER["coords"][:3]]
    + [[30.0 + xyz[0] - 3.0, xyz[1], xyz[2]] for xyz in DIMER["coords"][3:]],
)

# -- references, all PySCF on the geometries above --------------------------
#
# The correlated ones are frozen-core because this program freezes the core by
# default (keywords.correlation.freeze_core), so an all-electron reference
# would disagree by 2.3 mHa and read as a broken MP2.
HF_STO3G = -74.962005687948
UHF_OH_STO3G = -74.362633303788
B3LYP_STO3G = -75.311397471446          # grid level 3, the default here
MP2_CCPVDZ = -76.227939363772
MP2_CCPVDZ_ALL_ELECTRON = -76.230274593250
SCS_MP2_CCPVDZ = -76.224332224785
SOS_MP2_CCPVDZ = -76.222528655292
RI_MP2_CCPVDZ = -76.227922438649        # fitted with cc-pvdz-rifit
CCSD_STO3G = -75.011101318288
CCSD_T_STO3G = -75.011174129989
RI_CCSD_T_STO3G = -75.011173699760      # fitted with cc-pvdz-rifit
CASCI_44_STO3G = -74.969266882206       # CAS(4,4) on the RHF orbitals

# CASSCF has no reference here on purpose. On this molecule and basis the
# orbital optimisation settles on a different stationary point from PySCF's --
# a symmetry-trapped one, which a first-order method cannot leave without the
# second-order augmented Hessian step that would tell it the valley it is in is
# not a minimum. Not a wrong energy: it is a converged CASSCF of a different
# solution, and CASCI on the same space agrees to 1e-10. So a number pinned in
# here would be a number about the optimiser rather than about the interface,
# and would need changing the moment that step lands. The ordering below is
# what the interface can assert either way.

TOL = 1.0e-6         # against PySCF
TOL_TIGHT = 1.0e-9   # between two routes to the same quantity
TOL_FIT = 1.0e-3     # fitted against conventional: the fitting error itself


def water(**kwargs):
    system = mqc.System(**WATER, **kwargs)
    system.set_monomers([[0, 1, 2]], multiplicities=[kwargs.get("multiplicity", 1)])
    return system


def hydroxyl():
    system = mqc.System(**HYDROXYL, multiplicity=2)
    system.set_monomers([[0, 1]], multiplicities=[2])
    return system


def dimer(far=False):
    system = mqc.System(**(FAR_DIMER if far else DIMER))
    system.auto_monomers()
    return system


#: Quiet, because thirty calculations at the default level bury the table this
#: script exists to print. A failure still says why: refusals come back as
#: exceptions, not as log lines.
VERBOSITY = "warning"


def energy(label, **kwargs):
    """One unfragmented calculation. `level=0` is the whole system at once."""
    system = kwargs.pop("system", None) or water()
    kwargs.setdefault("verbosity", VERBOSITY)
    return mqc.MBE(system, level=0, **kwargs).run(
        label=label, write_to_file=False).energy


# ---------------------------------------------------------------------------
#  The cases
# ---------------------------------------------------------------------------
#
# Each returns (value, reference, tolerance, what the reference is). A case
# whose reference is None returns a bool instead: an identity that holds or
# does not.


def hartree_fock():
    return energy("m_hf", method="hf", basis="sto-3g"), HF_STO3G, TOL, "PySCF"


def hartree_fock_fitted():
    """Fitted J and K. A different energy from exact HF, not an approximation
    to be checked against it -- so its reference is PySCF's fitted number."""
    got = energy("m_hf_df", method="hf", basis="sto-3g",
                 aux_basis="def2-universal-jkfit", density_fitting=True)
    return got, -74.962092123400, TOL, "PySCF fitted"


def unrestricted():
    return (energy("m_uhf", system=hydroxyl(), method="hf", basis="sto-3g"),
            UHF_OH_STO3G, TOL, "PySCF UHF")


def kohn_sham():
    return (energy("m_b3lyp", method="dft", functional="b3lyp", basis="sto-3g"),
            B3LYP_STO3G, TOL, "PySCF, grid level 3")


def kohn_sham_grid():
    """The grid is a keyword block, and a finer one moves the energy.

    Not a reference case: what is being checked is that `dft=` reaches the
    quadrature at all. A grid level that changed nothing would mean the block
    was read, validated and dropped -- which is what happens to a setting that
    has no way through the interface.
    """
    coarse = energy("m_grid3", method="dft", functional="pbe", basis="sto-3g",
                    dft={"grid_level": 3})
    fine = energy("m_grid5", method="dft", functional="pbe", basis="sto-3g",
                  dft={"grid_level": 5})
    return abs(coarse - fine) > 1.0e-9, None, None, "grid level changes the energy"


def double_hybrid():
    """B2PLYP: a functional whose perturbative part is part of the functional.

    Reachable exactly like any other, which is the point of the case -- the
    only thing it needs beyond a hybrid is an auxiliary basis.
    """
    got = energy("m_b2plyp", method="dft", functional="b2plyp", basis="cc-pvdz")
    return got < 0.0, None, None, "runs and returns an energy"


def mp2():
    return energy("m_mp2", method="mp2", basis="cc-pvdz"), MP2_CCPVDZ, TOL, "PySCF"


def mp2_all_electron():
    """`correlation=` reaches the frozen core, and unfreezing it moves the energy."""
    got = energy("m_mp2_ae", method="mp2", basis="cc-pvdz",
                 correlation={"freeze_core": False})
    return got, MP2_CCPVDZ_ALL_ELECTRON, TOL, "PySCF all-electron"


def scs_mp2():
    return energy("m_scs", method="scs-mp2", basis="cc-pvdz"), SCS_MP2_CCPVDZ, TOL, "PySCF"


def sos_mp2():
    return energy("m_sos", method="sos-mp2", basis="cc-pvdz"), SOS_MP2_CCPVDZ, TOL, "PySCF"


def ri_mp2():
    got = energy("m_rimp2", method="ri-mp2", basis="cc-pvdz", aux_basis="cc-pvdz-rifit")
    return got, RI_MP2_CCPVDZ, TOL, "PySCF fitted"


def ccsd():
    return energy("m_ccsd", method="ccsd", basis="sto-3g"), CCSD_STO3G, TOL, "PySCF"


def ccsd_t():
    return energy("m_ccsd_t", method="ccsd(t)", basis="sto-3g"), CCSD_T_STO3G, TOL, "PySCF"


def ri_ccsd_t():
    got = energy("m_ri_ccsd_t", method="ri-ccsd(t)", basis="sto-3g",
                 aux_basis="cc-pvdz-rifit")
    return got, RI_CCSD_T_STO3G, TOL, "PySCF fitted"


def cc_spin_adapted_identity():
    """Spatial and spin-orbital coupled cluster are the same wavefunction.

    For a closed shell they are exact restatements of each other, so this is an
    identity rather than a comparison: anything above rounding is a bug in one
    of the two, and `cc={"spin_adapted": False}` is the only way to ask for the
    second from here.
    """
    spatial = energy("m_cc_spatial", method="ccsd", basis="sto-3g")
    spin_orbital = energy("m_cc_spinorb", method="ccsd", basis="sto-3g",
                          cc={"spin_adapted": False})
    return (abs(spatial - spin_orbital) <= TOL_TIGHT, None, None,
            f"spatial vs spin-orbital, diff {abs(spatial - spin_orbital):.2e}")


def cc_tolerance_reaches_the_solver():
    """`cc=` reaches the amplitude solver: a loose one stops somewhere else."""
    tight = energy("m_cc_tight", method="ccsd", basis="sto-3g",
                   cc={"tolerance": 1.0e-10})
    loose = energy("m_cc_loose", method="ccsd", basis="sto-3g",
                   cc={"tolerance": 1.0e-3})
    return abs(tight - loose) > 0.0, None, None, "amplitude tolerance is honoured"


def casci():
    got = energy("m_casci", method="casci", basis="sto-3g",
                 mcscf={"n_active_electrons": 4, "n_active_orbitals": 4})
    return got, CASCI_44_STO3G, TOL, "PySCF CASCI(4,4)"


def casscf_below_casci():
    """Optimising the orbitals cannot raise the energy of the same CI space.

    An ordering rather than a reference -- see the note beside the constants.
    """
    fixed = energy("m_casci2", method="casci", basis="sto-3g",
                   mcscf={"n_active_electrons": 4, "n_active_orbitals": 4})
    optimised = energy("m_casscf", method="casscf", basis="sto-3g",
                       mcscf={"n_active_electrons": 4, "n_active_orbitals": 4})
    return (optimised < fixed, None, None,
            f"CASSCF {optimised:.6f} < CASCI {fixed:.6f}")


def ormas():
    """Occupation windows over subspaces of the active space.

    A partition that allows every occupation the complete space allows is the
    complete space, so this asks for a restricted one and checks it lands above
    the CASCI it is a subset of.
    """
    full = energy("m_casci3", method="casci", basis="sto-3g",
                  mcscf={"n_active_electrons": 4, "n_active_orbitals": 4})
    restricted = energy("m_ormas", method="casci", basis="sto-3g",
                        mcscf={"n_active_electrons": 4, "n_active_orbitals": 4,
                               "ormas": {"subspaces": [1, 3],
                                         "min_electrons": [2, 0],
                                         "max_electrons": [4, 2]}})
    return (restricted >= full - TOL_TIGHT, None, None,
            f"ORMAS {restricted:.6f} >= CAS {full:.6f}")


def avas():
    """The active space chosen from atomic orbital character rather than counts."""
    got = energy("m_avas", method="casscf", basis="cc-pvdz",
                 mcscf={"avas": {"orbitals": ["O 2p"], "threshold": 0.2}})
    return got < 0.0, None, None, "AVAS selects a space and it runs"


def sapt0():
    """The interaction energy of two monomers, and its decomposition.

    The terms travel in the output document rather than through the energy, so
    this is the one case that must write files. Checking that they sum to the
    total is checking that the decomposition returned is the one the total was
    built from.
    """
    system = dimer()
    result = mqc.MBE(system, level=0, method="sapt0", basis="sto-3g",
                     verbosity=VERBOSITY).run(label="m_sapt0", write_to_file=True)
    terms = result.sapt or {}
    parts = ("elst10", "exch10", "ind20_r", "exch_ind20_r", "disp20",
             "exch_disp20", "delta_hf")
    if not all(p in terms for p in parts):
        return False, None, None, f"terms missing: {sorted(terms)}"
    total = sum(terms[p] for p in parts)
    return (abs(total - result.energy) <= TOL_TIGHT, None, None,
            f"terms sum to the total, diff {abs(total - result.energy):.2e}")


def makefp_and_efp2():
    """The pair that only makes sense together.

    MAKEFP solves one wavefunction and writes it out as a potential; EFP2
    solves none and evaluates the interaction between potentials already
    written. So the first is what makes the second possible, and neither was
    reachable from Python until the potentials could be named on the system.

    The check is the long-range limit: two waters 30 A apart interact by
    essentially nothing, whatever the terms do at 3 A.
    """
    monomer = water()
    mqc.MBE(monomer, level=0, method="hf", basis="6-31g", driver="makefp",
            verbosity=VERBOSITY).run(label="m_water", write_to_file=False)
    if not os.path.exists("m_water.efp"):
        return False, None, None, "makefp wrote no potential"

    close = dimer()
    close.set_fragment_potentials(["m_water.efp", "m_water.efp"])
    near = mqc.MBE(close, level=0, method="efp2", verbosity=VERBOSITY).run(
        label="m_efp_near", write_to_file=False).energy

    apart = dimer(far=True)
    apart.set_fragment_potentials(["m_water.efp", "m_water.efp"])
    far = mqc.MBE(apart, level=0, method="efp2", verbosity=VERBOSITY).run(
        label="m_efp_far", write_to_file=False).energy

    return (abs(far) < 1.0e-5 and abs(near) > abs(far), None, None,
            f"3 A {near:+.6f}, 30 A {far:+.2e} Hartree")


def gradient():
    """Any driver a deck accepts, including the ones that return no energy."""
    result = mqc.MBE(water(), level=0, method="hf", basis="sto-3g",
                     driver="gradient", verbosity=VERBOSITY).run(
        label="m_grad", write_to_file=True)
    norm = result.gradient_norm
    return (norm is not None and abs(norm - 0.086151389909) <= TOL,
            None, None, f"gradient norm {norm}")


def bonding_analysis():
    """A property, not a keyword: it asks for something further to be done with
    the wavefunction and changes nothing about the energy."""
    plain = energy("m_plain", method="hf", basis="6-31g")
    analysed = energy("m_quao", method="hf", basis="6-31g",
                      properties={"bonding_analysis": {"type": "gms_quao"}})
    return (abs(plain - analysed) <= TOL_TIGHT, None, None,
            "the analysis leaves the energy alone")


# -- what must be refused ---------------------------------------------------
#
# Each of these returned a plausible number or killed the interpreter before
# it was fixed, which is why they are cases rather than documentation.


def refuses_unknown_method():
    return _raises(lambda: energy("m_bogus", method="mp3", basis="sto-3g"),
                   "unknown model.method")


def refuses_unimplemented_method():
    return _raises(lambda: energy("m_f12", method="mp2-f12", basis="sto-3g"),
                   "not implemented")


def refuses_unknown_driver():
    return _raises(lambda: energy("m_baddriver", method="hf", basis="sto-3g",
                                  driver="banana"), "unknown driver")


def refuses_optimize():
    """The optimizer drives run_calculation rather than being driven by it, so
    a session cannot host one. It used to abort the communicator, which took
    the interpreter with it."""
    return _raises(lambda: energy("m_opt", method="hf", basis="sto-3g",
                                  driver="optimize"), "not available")


def refuses_makefp_that_wrote_nothing():
    """MAKEFP returns no energy, so a failed one is invisible without this.

    A driver that produces a file and not a number has only two outcomes a
    caller can see: an error, or a file. Without the error it had one.
    """
    return _raises(lambda: energy("m_bad_fp", method="hf", basis="not-a-basis",
                                  driver="makefp"), "MAKEFP failed")


def refuses_efp_without_potentials():
    """This returned 0.0 as a success, which is a plausible interaction energy."""
    return _raises(lambda: energy("m_efp_bare", system=dimer(), method="efp2"),
                   "no fragment carries a potential")


def refuses_sapt_on_three_monomers():
    """SAPT partitions H into H_A + H_B + V. There is no slot for a third.

    Three whole waters rather than a dimer cut into three pieces: cutting a
    molecule trips the uncapped-valence check first, which is a different
    refusal and would pass this case for the wrong reason.
    """
    coords = [list(xyz) for xyz in FAR_DIMER["coords"]]
    coords += [[x + 60.0, y, z] for x, y, z in FAR_DIMER["coords"][:3]]
    system = mqc.System(symbols=FAR_DIMER["symbols"] + ["O", "H", "H"], coords=coords)
    system.auto_monomers()
    return _raises(lambda: energy("m_sapt3", system=system, method="sapt0",
                                  basis="sto-3g"), "exactly two fragments")


def _raises(call, expected):
    try:
        call()
    except mqc.MQCError as exc:
        message = str(exc)
        return (expected in message, None, None,
                f"refused: {message.split('.')[0][:70]}")
    return False, None, None, "was not refused"


CASES = [
    ("hf", hartree_fock),
    ("hf/df", hartree_fock_fitted),
    ("uhf", unrestricted),
    ("dft b3lyp", kohn_sham),
    ("dft grid", kohn_sham_grid),
    ("dft b2plyp", double_hybrid),
    ("mp2", mp2),
    ("mp2 all-electron", mp2_all_electron),
    ("scs-mp2", scs_mp2),
    ("sos-mp2", sos_mp2),
    ("ri-mp2", ri_mp2),
    ("ccsd", ccsd),
    ("ccsd(t)", ccsd_t),
    ("ri-ccsd(t)", ri_ccsd_t),
    ("cc spin-adapted", cc_spin_adapted_identity),
    ("cc tolerance", cc_tolerance_reaches_the_solver),
    ("casci", casci),
    ("casscf", casscf_below_casci),
    ("ormas", ormas),
    ("avas", avas),
    ("sapt0", sapt0),
    ("makefp + efp2", makefp_and_efp2),
    ("gradient", gradient),
    ("bonding analysis", bonding_analysis),
    ("refuse method", refuses_unknown_method),
    ("refuse f12", refuses_unimplemented_method),
    ("refuse driver", refuses_unknown_driver),
    ("refuse optimize", refuses_optimize),
    ("refuse makefp", refuses_makefp_that_wrote_nothing),
    ("refuse bare efp", refuses_efp_without_potentials),
    ("refuse sapt(3)", refuses_sapt_on_three_monomers),
]

#: What a build without a piece says when asked for it. A reduced build is not
#: a regression, so these are skipped and named rather than failed -- but only
#: these: any other refusal is a failure.
MISSING_BUILD = ("libxc", "libcint", "cuEST", "not built", "build with")


def main(argv):
    pattern = argv[1] if len(argv) > 1 else ""
    failures, skipped = [], []

    with mqc.session():
        for name, case in CASES:
            if pattern and pattern not in name:
                continue
            try:
                value, reference, tolerance, source = case()
            except mqc.MQCError as exc:
                message = str(exc)
                if any(hint in message for hint in MISSING_BUILD):
                    skipped.append((name, message.split(".")[0][:70]))
                    print(f"  {name:<18} {'skipped':>18}   {message.split('.')[0][:60]}")
                    continue
                failures.append(f"{name}: {message}")
                print(f"  {name:<18} {'RAISED':>18}   {message[:70]}")
                continue

            if reference is None:
                ok = bool(value)
                print(f"  {name:<18} {('ok' if ok else 'FAIL'):>18}   {source}")
                if not ok:
                    failures.append(f"{name}: {source}")
            else:
                delta = abs(value - reference)
                ok = delta <= tolerance
                print(f"  {name:<18} {value:>18.9f}   ref {reference:>18.9f} "
                      f"({source})   diff {delta:.2e}   {'ok' if ok else 'FAIL'}")
                if not ok:
                    failures.append(
                        f"{name}: {value} vs {reference}, diff {delta:.2e} > {tolerance:.1e}")
            sys.stdout.flush()

    if skipped:
        print(f"\n{len(skipped)} case(s) skipped, this build cannot run them:")
        for name, why in skipped:
            print(f"  {name}: {why}")
    if failures:
        print("\nFAILED:")
        for failure in failures:
            print("  " + failure)
        return 1
    print("\nevery method reachable from Python agrees with its reference")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
