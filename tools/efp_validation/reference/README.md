# GAMESS MAKEFP reference potentials

Generated locally, so they can be regenerated. These are the reference for M4
(distributed multipoles) and M5 (distributed polarizabilities), which are the
first milestones PySCF cannot check.

## Regenerating

There is a working GAMESS build at `../mgga/gamess`:

```bash
cd ../mgga/gamess
cp <this dir>/water_6-31gs_boys.inp water_makefp.inp
./misc/automation/rungms water_makefp 00 1 1 > water.log 2>&1
# the potential lands in restart/water_makefp/<timestamp>/water_makefp.efp
```

Use **`misc/automation/rungms`**, not the top-level `rungms` — the latter dies
with `Illegal variable name` in this checkout, and the test harness
(`tests/runtest.py`) uses the automation copy for the same reason.

The potential is written to a `.efp` file *beside* the punch file, in
`restart/<job>/<timestamp>/`, not to the working directory.

## What is here

| file | |
|---|---|
| `water_6-31gs_boys.efp` | complete. 5 expansion points, 4 polarizable points, 12 dynamic frequencies |
| `water_6-31gs_boys.inp` | `RUNTYP=MAKEFP LOCAL=BOYS`, 6-31G\*, `ISPHER` left at its Cartesian default |
| `hcn_6-31gs_boys.inp` | |
| `hcn_6-31gs_boys.partial.efp` | **incomplete — see below** |

## HCN does not finish, and it is not our fault

The HCN run aborts with

```
*** TOO MANY ITERATIONS IN AOCPCG *** MAX CPHF=   50
AS LONG AS EQUATIONS REMAIN CONVERGENT, GRACE LIMIT= 299
MOST OFTEN THIS IS CAUSED BY USE OF AN INAPPROPRIATE WAVEFUNCTION.
CHANGE SCFTYP, OR CHECK HOMO/LUMO FILLING ORDER.
```

so `hcn_6-31gs_boys.partial.efp` stops after `OCTUPOLES` — the multipoles are
there and usable, the polarizabilities are not.

**The diagnostic is misleading.** The reference is a closed shell with a 0.685
Hartree HOMO-LUMO gap, and the orbital Hessian for this molecule and basis is
positive definite with eigenvalues in [0.219, 19.7] — condition number 90, and
**6.3 after Jacobi preconditioning**. Our own CPHF converges it in 14 iterations
(`validation/check_cphf.f90` case 4, agreeing with a dense exact solve to 1e-9
relative). The residual in the GAMESS log oscillates between 2.8e-3 and 5.8e-3
without trending, which is a stalled AO-basis iteration and not an ill-posed
problem.

Consequence for the plan: **the M5 reference has to come from a molecule GAMESS
can finish.** Water does. Do not assume a missing section means the physics is
hard.

## What these files already established

Checked by `../check_sum_rules.py`, which needs no metalquicha build:

- The distributed multipoles satisfy the **sum rule against PySCF** — total dipole
  to 4.8e-8 and the full second moment to 2.8e-7, at the file's printed precision.
  That confirms, all at once: `MONOPOLES` columns are (electronic, nuclear);
  dipoles and quadrupoles are electronic only with the nucleus in the monopole
  alone; the quadrupole packing is `XX YY ZZ XY XZ YZ`; the moments are raw
  Cartesian, not traceless (trace −12.55); and bond midpoints are plain
  arithmetic means.
- Per-LMO polarizability tensors are **asymmetric** (0.17 for the bonds, 0.06 for
  the lone pairs) while their **sum is symmetric to 4.6e-6**. M5 must not
  symmetrize them.
- **The core orbital is excluded.** Water has 5 occupied MOs and 4 polarizable
  points, and the 4 tensors sum to an isotropic 4.8658 against 4.8666 for an
  exact all-occupied solve. The missing 7.7e-4 is diagonal-only, which is what a
  spherical 1s should contribute.
- Polarizable points sit exactly at the LMO centroids — identical coordinates in
  both sections.
- A default MAKEFP run on **linear HCN adds no dummy point**: 3 atoms and 2 bond
  midpoints, nothing else. The shipped `HCN.efp` has an `A04D4`, so that comes
  from its own input rather than from MAKEFP's handling of linear molecules. One
  less unknown for M4 than the plan feared.


## Does GAMESS accept a potential *we* wrote

`../dimer_energy.py` is the end-to-end test, and the only one here that judges the
file rather than the numbers in it:

```bash
./build/check_makefp                                # writes /tmp/mqc_water.efp
python3 tools/efp_validation/dimer_energy.py
python3 tools/efp_validation/dimer_energy.py --screen-from-gamess
```

It builds the same water dimer twice -- once from our potential, once from
GAMESS's -- as an all-EFP job and compares the energy terms GAMESS reports.

| term | ours vs GAMESS | |
|---|---|---|
| exchange repulsion | **0** | exact |
| charge transfer | **0** | exact, against `water_6-31gs_cmo.efp` |
| electrostatics | 1.0e-10 | with `--screen-from-gamess` |
| polarization | 1.2e-09 | |
| dispersion E6 | 1.7e-09 | |
| dispersion E7 | 4.0e-10 | |
| dispersion E8 | 2.0e-07 | |
| dispersion total | 2.0e-07 | |
| electrostatics | 4.4e-05 | with our own screening fit |

All seventeen sections GAMESS's reader recognises are written, so every EFP energy
component above comes from parameters this codebase computed.

**Two `CTVEC` forms, and they are different bases.** By default MAKEFP writes five
occupied orbitals plus two quasi-atomic valence virtuals built by `VVOS`; with
`$MAKEFP CTVVO=.FALSE.` it writes the whole canonical MO matrix instead, which is
the form we emit and the one GAMESS itself recommends when the valence virtuals
cannot be formed. GAMESS reads either. Compared against `water_6-31gs_boys.efp`
charge transfer therefore differs by 2.6e-05 -- 38 EFP MOs against 14, a larger
virtual space -- and against `water_6-31gs_cmo.efp` it agrees exactly. That second
file is here for exactly this reason.

**Three things about the deck, each of which cost a run to find.** An all-EFP
system needs `COORD=FRAGONLY` and *no* `$DATA` group -- with `$DATA` GAMESS asks
for basis functions it has no atoms to put a basis on. The fragment must not be
named `WATER`, because GAMESS ships a built-in fragment by that name, so a deck
asking for `WATER` silently gets the internal EFP1 potential and reads neither
file, while running to completion and printing plausible energies; both are
renamed to `EFPTEST`. And `CTFOK` must not be written unless `CTVEC` is: it is a
subsection of `CTVEC` rather than a section, and GAMESS aborts on a standalone one
(its reader rejects it). That last one is why this test exists -- every parameter
in the file already matched GAMESS's own and the file was still unreadable.

**Why the `--screen-from-gamess` flag exists.** Electrostatics is a multipole
expansion plus a charge-penetration correction, and both are ours, so a
disagreement cannot be attributed by looking at it. Substituting only GAMESS's
screening moves the term from 4.4e-05 to 1.0e-10, which places the entire
difference in the damping fit and none of it in the moments. The fit difference is
itself understood: our optimizer finds a real minimum on the bond midpoint where
GAMESS's search reaches its `alpha = 10` "off" bound, an absorbing state its
objective is flat in.

**Dispersion is not zero, which was worth learning.** Our potential carries the
dipole-dipole dynamic polarizabilities and GAMESS builds `E6` from them, so it
reports a dispersion energy from our file -- that section is live, not merely
present. GAMESS's own value additionally contains `E7` and `E8` from the two
blocks we cannot write, so the two numbers are not comparable term by term.


## What sets the floor on the agreement

Not the integrals and not the units, both of which are checkable:

* GAMESS reports `FINAL RHF ENERGY IS -76.0105013381` for the reference geometry
  and this code gets `-76.0105013381` -- identical to all ten printed digits, so
  libcint and GAMESS's Rys quadrature agree far below anything visible in the
  parameters;
* everything is computed and written in Bohr, from a geometry converted with
  GAMESS's own old CODATA value `0.52917724924`, so no unit conversion is left to
  disagree about.

**It is the basis table.** Comparing the `PROJECTION BASIS SET` block we emit
against GAMESS's own, 5 of 42 numbers differ and all five are *exponents* -- the
contraction coefficients are exact:

```
ours   15.5396162500      GAMESS 15.5396162453
ours    2.8253943650      GAMESS  2.8253943648
```

Worst relative difference 3.0e-10. Basis Set Exchange ships `15.53961625` to ten
significant figures and GAMESS's internal `N31` table carries more, so our JSON is
the truncated one. That is the floor under the 1e-10 to 1e-09 agreements above, and
it is the same class of problem noted for the PySCF ab initio references: feeding a
code a basis table that differs in the tenth digit looks exactly like a bug in
whichever code you are checking.

**The one term this does not explain is `E8`, at 2.0e-07.** That is three orders
above this floor. It is the per-orbital pair-antisymmetric difference in `LMOQQPOL`
described in `backends/libcint/mqc_efp_potential.f90`, which is structural rather
than numerical -- worth keeping distinct, because "it is just the basis" is a
comfortable and wrong explanation for it.


## Does it describe the same *surface*

`../hexamer_optimize.py` relaxes a water hexamer -- the prism from
`validation/inputs/sample_inputs/` -- as six rigid EFP fragments, once from our
potential and once from GAMESS's, and compares where they settle:

```bash
./build/check_makefp
python3 tools/efp_validation/hexamer_optimize.py --screen-from-gamess
```

A single-point comparison tests the parameters; it cannot tell whether they describe
a surface, because a potential can have the right energy at one geometry and the
wrong gradient. Thirty-six degrees of freedom over about thirty-four steps does test
that, and the answer is the structure.

| | with the reference's screening | with our own fit |
|---|---|---|
| final energy | 1.3e-06 Hartree | 1.2e-03 Hartree |
| 15 oxygen-oxygen distances | 4.7e-05 Angstrom | 1.2e-02 Angstrom |
| optimization steps | 34 against 34 | 37 against 34 |

So the surface is right, and the charge-penetration screening is the one thing that
differs -- the same conclusion the dimer single point reaches, now with a physical
size attached. **The screening fit is worth 0.74 kcal/mol and 0.012 Angstrom on a
water hexamer**, which is not negligible for fragment work, and it is worth
remembering that ours finds a genuine minimum where GAMESS's search reaches its
`alpha = 10` "off" bound on one bond midpoint. Neither is obviously the better fit;
they are different fits, and only one of them reproduces GAMESS.

**Basis.** This runs in 6-31G*, not the 6-311++G(3df,2p) usually recommended for
EFP2, for two reasons. That set is not in the Basis Set Exchange bundle at all -- it
carries `6-311++g`, `(2d,2p)`, `(3df,3pd)` and `**` -- and any basis with f functions
cannot be emitted yet, because the AO-ordering map in `mqc_efp_potential` handles
Cartesian d and stops there. `6-311++G**` is d-only and would work once added to
`MQC_BASIS_SETS`.
