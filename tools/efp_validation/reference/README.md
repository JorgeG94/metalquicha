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
