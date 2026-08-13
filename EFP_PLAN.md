# Toward MAKEFP: generating effective fragment potentials

Scratch plan. The destination is GAMESS's `RUNTYP=MAKEFP`: take one molecule,
emit the parameters that let it act as a polarizable classical fragment in
somebody else's calculation. EFP2 — the general *ab initio* flavour — not EFP1,
which is a water-specific fitted potential and needs none of this.

## The target format, read off a real file

GAMESS is built locally at `../mgga/gamess` and ships eight reference potentials
in `auxdata/EFP/` (HCN, MEOH, BENZEN, URACIL, INDOLE, PEPTID, FORMID, PYRIDO).
`HCN.efp` is 520 lines and lays out exactly what has to be produced:

```
$HCN
Hydrogen cyanide  Elec.: 6-31+G*  Rest: 6-311++G(3df,2p)
 COORDINATES (BOHR)     label  x y z  mass  Z
 MONOPOLES              label  electronic  nuclear
 DIPOLES                label  x y z
 QUADRUPOLES            label  6 components
 OCTUPOLES              label  10 components
 POLARIZABLE POINTS     CT n   centroid + 9 polarizability components
 DYNAMIC POLARIZABLE POINTS     the same, 12 imaginary frequencies
 PROJECTION BASIS SET   shells and exponents, per atom
 MULTIPLICITY
 PROJECTION WAVEFUNCTION   n_occ x n_basis MO coefficients
 FOCK MATRIX ELEMENTS
 LMO CENTROIDS
 <per-atom screening>   label  coefficient  exponent
```

Six things in there are worth noticing before planning anything:

1. **Two basis sets.** Electrostatics and polarization in one (`6-31+G*`),
   exchange-repulsion in a much larger one (`6-311++G(3df,2p)`). MAKEFP is
   effectively two wavefunctions, not one.
2. **Expansion points are atoms *and* bond midpoints *and* dummy points.**
   `A01N1 A02C2 A03H3` are atoms, `BO21 BO32` are bond midpoints, and `A04D4` is
   a dummy carrying no mass or charge, placed beyond the terminal H. Linear and
   small molecules need extra points to represent the density adequately, and
   where they go is a GAMESS convention we will have to match.
3. **Monopoles are split** into an electronic part and the nuclear charge, not
   summed.
4. **Polarizable points sit at LMO centroids** (`CT n`), one tensor each, and the
   nine components are the full 3×3 — not a symmetrized six.
5. **Twelve imaginary frequencies**, confirmed by counting: 60 dynamic entries
   over 5 centroids.
6. **The exchange-repulsion data is dumped, not reduced.** MO coefficients, the
   Fock matrix, the basis. The energy is formed at *use* time from inter-fragment
   LMO overlap, so MAKEFP only has to emit ingredients — which makes that
   milestone serialization rather than physics.

## Our format is JSON; GAMESS's format is a view onto it

**We emit JSON. A Python tool renders it into GAMESS's `.efp` text, and parses
GAMESS's `.efp` back, purely so the two can be compared.** The `.efp` layout above
is the *specification we must be able to express*, not the file our Fortran
writes.

This is the better split for four reasons:

- **Comparison becomes structured.** Parse both sides into the same canonical
  dictionary and compare parameter by parameter with per-quantity tolerances, the
  way `validation_tests_cpu.json` already holds PySCF numbers. Text diffing
  fixed-width columns would fail on trailing zeros and tell us nothing about
  which parameter moved.
- **The `>` continuations stay out of Fortran.** GAMESS wraps quadrupoles,
  octopoles and polarizability tensors across lines with a continuation marker
  and specific column widths. That is unpleasant to emit from Fortran and trivial
  in Python.
- **The durable artifact is ours.** If GAMESS's format shifts we fix a Python
  script, and the parameters themselves never move.
- **It matches the program.** Decks are JSON in, `output_<basename>.json` out,
  through `mqc_json_output_types`. An EFP potential is just another output
  document.

**The Python tool is test scaffolding, not a product.** It lives in `tools/`
beside `cpu_validation/`, never in `python/`, and nothing ships it. So it carries
no API to keep stable, needs no packaging or docs beyond a docstring, and may
depend on whatever is convenient — numpy, whatever GAMESS parsing wants. The only
thing it owes anyone is that a failure points at the right parameter.

That also sets how good the renderer has to be: good enough for GAMESS to read,
not byte-perfect. Byte fidelity in the round-trip test below is a *means* of
proving the parser understood every field, not a requirement of its own.

The tool has two jobs, and the *parser* is the one that does the scientific work:

1. **parse** a GAMESS `.efp` into the canonical dictionary — this is what makes
   parameter-by-parameter comparison possible at every milestone;
2. **render** our JSON into a GAMESS `.efp` — this is what makes the end-to-end
   check possible, where GAMESS consumes our potential.

Neither needs to be written up front. The parser is wanted at M4, the first
milestone that produces comparable numbers; the renderer not until M7.

## The reference strategy

**Generate our own reference with the local GAMESS, in a basis we already
ship.** The shipped files are the authority on *format*, but they were made with
`6-31+G*` and `6-311++G(3df,2p)` and we ship neither — `basis_sets/` has 42
sets including `6-31G`, `6-31G*` and `6-31G**`, none of them diffuse. Chasing
those two bases from BSE just to diff against a shipped file would put a
basis-acquisition task on the critical path for no physics.

So: run MAKEFP ourselves on small molecules in `6-31G*`, and diff parameter by
parameter. That is the whole reason this project is tractable — every milestone
below has an external check available the day it produces its first number, and a
disagreement localises to one parameter set instead of to "the interaction energy
is wrong".

Two conditions, and both go in the comparison script rather than in anyone's head:

- **Same basis**, since EFP2 parameters are basis-dependent quantities, not
  converged properties.
- **Same localization.** `LOCAL=BOYS` and `RUEDNBRG` give different LMOs and
  therefore different distributed parameters. Both are legitimate; only agreement
  is at stake, so ours must match the reference run.

Check the generated `.efp` files into the repo as reference data, the way the CPU
manifest holds PySCF numbers. They are small text files, and a reference that
needs GAMESS installed to reproduce is a reference that rots on every machine
without it.

## What we already have

- A validated RHF/UHF wavefunction on the CPU, with MO coefficients and orbital
  energies in hand.
- `mqc_bond_perception` — bond midpoints come free and are already tested.
- `pic_gemm`, DIIS, a working direct and in-core Fock build.
- On the GPU side `cuestMultipoleCompute` is already bound and uncalled, exactly
  the state PCM was in before it was wired.

## What is missing

- **Multipole integrals.** libcint has `int1e_r` and friends; none are exposed in
  the Fortran interface — the same gap `int1e_grids` had before PCM. The `env`
  slot constants are libcint's own 0-based offsets and are *not* converted by the
  interface; that cost a day on PCM, so `+ 1` and a test.
- **Orbital localization.** Nothing anywhere.
- **CPHF.** Nothing, and it is the shared bottleneck: analytical Hessians,
  MP2/CC gradients and SAPT induction all want it too, so it is the one piece
  here that pays for itself more than once.
- **Frequency-dependent CPHF** at imaginary frequency.
- **A JSON potential document**, and the Python tool that parses and renders
  GAMESS's `.efp` around it.
- **A dipole on the CPU path at all** — only cuEST computes one today. M1 fixes
  that as a byproduct.

## Who validates what

Checked against the installed PySCF, not assumed:

| Step | Reference | Available? |
|---|---|---|
| Multipole integrals | `mol.intor('int1e_r'/'int1e_rr'/'int1e_rrr')` | **yes**, through `int1e_rrrr` |
| Localization | `pyscf.lo.Boys`, `pyscf.lo.ER`, `pyscf.lo.PM` | **yes** — and `ER` is GAMESS's `RUEDNBRG`, so either localization choice is covered |
| CPHF, total polarizability | `pyscf.scf.cphf` + finite field | **yes** (`pyscf.prop.polarizability` is *not* installed, so use `cphf.solve` and finite difference) |
| **Distributed** multipoles | — | **no.** PySCF has no Stone DMA. GAMESS only |
| **Distributed** polarizabilities | — | **no.** GAMESS only |
| Dynamic α at imaginary ω | — | GAMESS only |
| Exchange-repulsion data | trivially, it is MOs and a Fock matrix | yes |

**That splits the project cleanly, and it is the reason to reorder.** Everything
that builds *machinery* — integrals, localization, CPHF — has a fast local
reference. Only the *distribution* steps need GAMESS, and they come after. So the
order below front-loads every milestone with a seconds-long feedback loop and
leaves the slow, format-archaeology work for when the machinery is already known
good.

## Settled from the GAMESS source

Facts established by reading `source/` (citations are GAMESS's own files), each
either corroborated against the 31 reference potentials or flagged where it was
not. **Two claims that came back from the source reading were checked against the
files and found wrong** — recorded here so nobody acts on them: it is *not* true
that `ADENI2.efp` and `FORMIC.efp` lack a screening block (all 31 carry exactly
one `SCREEN2`), and it is *not* true that a default run emits both `SCREEN2` and
`SCREEN` (no shipped file has both). Treat source-derived claims as hypotheses
until a reference file agrees.

### Multipoles — `prppop.src` (`STONE`/`STNRD`/`STNMOM`/`STNXYZ`)

- **Partitioning is nearest-expansion-point, winner-take-all**, on the *Gaussian
  product centre* of each primitive pair. No overlap weighting, no Becke weights.
  Ties within 1e-6 bohr² split **equally**, each tied point taking `1/n`. The
  Becke-grid alternative is gated behind `BIGEXP`, which defaults to zero, so a
  default deck never reaches it. This was the plan's stated top risk and it is
  much simpler than feared.
- **Moments are raw primitive Cartesian, not traceless, no prefactor.** The
  Buckingham conversion happens in the *consumer* (`efelec.src`), not the
  generator. Confirmed empirically: HCN's quadrupole trace is -12.40, not 0. So we
  emit raw moments — less work than planned.
- **Packing.** Quadrupole `XX YY ZZ XY XZ YZ`; octopole
  `XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ`. Cross-checked from GAMESS's *reader*
  as well as its writer.
- Each moment is about **its own** expansion point; there is no common origin.
  Dipole and higher are **electronic only** — the nucleus enters the monopole
  alone.
- `MONOPOLES` columns are (electronic, nuclear). Verified on all 31 files:
  Sum(electronic) = -Sum(Z) to the printed precision, and the dipole assembled as
  Sum(q·R + mu) is origin-independent to 1e-10.
- Bonds: `|Ri-Rj| <= rad_i + rad_j` in Angstrom on an Emsley covalent-radius
  table; midpoint is the plain arithmetic mean; label `BO<hi><lo>`.

### Screening — `chgpen.src` (`CGP`/`CGPFIT`/`CGPPRT`), grid in `prpel.src`

- The two numbers are **(beta, alpha)**: linear coefficient and damping exponent,
  alpha in bohr^-1.
- **They are fitted, not computed.** Powell direction-set optimization of alpha
  against the *electronic-only* QM electrostatic potential, on a geodesic grid
  over fused scaled van der Waals spheres centred on every expansion point
  (`VDWSCL=0.7`, `VDWINC=0.1`, 25 layers), objective a layer-weighted RMS in
  kcal/mol, `FTOL=1e-6`. Initial alpha 2.0 on atoms, 4.0 on bond midpoints.
- **beta is structurally 1.0** for the exponential fits — forced whenever the fit
  type is exponential — which is why the first column is always `1.000000000`.
- **alpha is bounded to [0.5, 10.0], and 10.0 means "off" for that centre.**
  Corroborated: exactly `10.000000000` appears 134 times across the 31 files.
  One value is `10.000078788`, so the bound constrains the search rather than the
  returned number.
- **The fit form and the energy form are different functions.** The fit damps only
  the monopole term with a one-centre `1 - exp(-alpha r)`; the EFP-EFP energy uses
  a two-centre Freitag/Slipchenko expression built from pairs of alphas. Anyone
  reimplementing must not assume one is the other.
- The block is optional in principle, and absent means simply no charge-penetration
  correction.

**Consequence for the plan: screening is its own milestone, not part of M4.** It
needs a geodesic ESP grid, layer weighting and a Powell optimizer — none of which
the multipoles require. It is also independently checkable, since the fitted
alphas are in every reference file.

## Milestones

Each states its own verification, and none depends on the next being right.

### M1 — Multipole integrals  *(PySCF)*

Expose `int1e_r`, `int1e_rr`, `int1e_rrr` and add them to `one_electron`'s
dispatch.

*The detail that will bite.* libcint returns **full Cartesian tensors** — 3, 9 and
27 components — while GAMESS's `.efp` stores **6 quadrupole and 10 octopole**
unique components in its own order. So there is a packing-and-ordering map between
them, and it is exactly the kind of thing that silently transposes a tensor. Get
it wrong and the sum rule at M4 still passes, because a transposed quadrupole sums
to a transposed total.

*Verify.* Elementwise against `mol.intor`, on the full tensors, before any
packing. Then `Tr(D·r)` against the molecular dipole — which the CPU path cannot
currently produce at all, so this is the first time that number exists there.
Fix an origin and say which: a neutral molecule's dipole is origin-independent and
its quadrupole is not.

*Size.* Small. Bindings are mechanical; the ordering is the thinking.

### M2 — Localization  *(PySCF)*

Boys: maximize Σ|⟨i|**r**|i⟩|² over occupied–occupied rotations. Jacobi sweeps are
plenty at fragment sizes.

*Verify.* Internal first, needing nothing: the SCF energy is invariant under the
rotation (exact, and what catches a non-unitary transform), the functional
increases monotonically, water gives two bonds and two lone pairs. Then against
`pyscf.lo.Boys` — orbital by orbital up to sign and degenerate-subspace mixing,
and centroid by centroid, which is the comparison that actually matters since the
centroids are what M5 places polarizabilities on.

`pyscf.lo.ER` is there too, so if the reference `MAKEFP` runs use `RUEDNBRG` we can
follow without losing a reference.

*Size.* Small-to-medium, self-contained.

### M3 — CPHF and total static polarizability  *(PySCF)*  ✅

`mqc_libcint_cphf.f90`. Preconditioned conjugate gradients on the occupied-virtual
block; the two-electron term of the response operator is one ordinary Fock build
against a symmetrized response density, so the full Hessian is never formed and
nothing new had to be written to contract it. RHF only — a functional would need
the XC kernel, and EFP2 is defined over a Hartree-Fock reference anyway.

Converges in 8 (STO-3G) to 15 (aug-cc-pVDZ) iterations. `run_libcint_rhf` gained
an optional `h_extra`, the hook a finite field needs.

*Verified.* Three references, and their disagreements were the interesting part:

| reference | agreement | notes |
|---|---|---|
| finite field through our own SCF | 2e-6 | in the unit tests; needs nothing external, runs every commit |
| dense exact solve of the same Hessian, from PySCF's integrals | **4e-9** | the primary reference: LAPACK, no solver tolerance in it at all |
| `pyscf.scf.cphf.solve` | 5.7e-7 | *PySCF's* residual, not ours — see below |
| Richardson finite field through PySCF | 3.5e-8 | independent of both response solvers |

Two things worth remembering. **PySCF's Krylov CPHF is the least accurate of the
three** on aug-cc-pVDZ, and its answer does not move between `tol=1e-9` and
`tol=1e-15` — so it cannot be tightened from outside, and had it been the primary
reference the choice would have been to accept 1e-6 agreement or to "fix" code
already correct to 4e-9. The dense solve settled it. **A finite field is limited by
SCF convergence amplified by 1/F**, not by truncation: PySCF defaults
`conv_tol_grad` to sqrt(conv_tol) ≈ 3e-7, and the giveaway was that Richardson
extrapolation made agreement *worse* — cancelling a term that was not dominant
while amplifying the noise that was.

Also worth carrying forward: coupling *raises* water's polarizability here (20% in
x, 90% out of plane), because the occupied-virtual Coulomb repulsion `(aa|ii)`
enters the response diagonal negatively and dominates it. Uncoupled Hartree-Fock
underestimates polarizabilities. The direction is basis dependent, so the unit test
treats it as a regression guard rather than a law.

α_iso: 2.46 (STO-3G) → 5.02 (cc-pVDZ) → 8.14 (aug-cc-pVDZ) Bohr³, against ~9.6
experimental. The diffuse basis is the only one where the number means anything —
STO-3G gives water an out-of-plane α of 0.04, having no out-of-plane flexibility at
all. Worth remembering for M5: distributed polarizabilities will need a diffuse
basis to be physical, whatever GAMESS's default happens to be.

*Leverage.* `cphf_solve` takes arbitrary AO perturbations, so analytical Hessians,
MP2/CC gradients and SAPT induction can reuse it unchanged. M5 needs the same
solutions resolved per localized orbital rather than summed.

**Everything above is PySCF-checkable and independent of GAMESS. Everything below
needs it.**

### M4 — Distributed multipoles  *(GAMESS)*  ✅

Partition the density among expansion points — atoms and bond midpoints — and
compute charge through octopole at each.

**There is a working GAMESS at `../mgga/gamess`, and the reference now lives in
`tools/efp_validation/reference/`.** Run it with `misc/automation/rungms`, not the
top-level `rungms`, which dies with `Illegal variable name`; the potential is
written to a `.efp` beside the punch file under `restart/<job>/<timestamp>/`.
Water/6-31G\*/`LOCAL=BOYS` takes two seconds.

*Verify.* **The sum rule, an exact identity:** translate every distributed
multipole to a common origin, add, and recover M1's total molecular multipoles.
Note its limit — it constrains the *total*, not the *distribution*, and not fully
the component ordering. Only the `.efp` comparison constrains those.

`tools/efp_validation/check_sum_rules.py` runs it against the reference *before we
can generate one*, which turned the following from hypotheses into facts (total
dipole to 4.9e-8, full second moment to 3.0e-7, at the file's printed precision):

- `MONOPOLES` columns are (electronic, nuclear); electron count recovers 10.0000000
- dipoles and quadrupoles are **electronic only**, nucleus in the monopole alone
- quadrupole packing `XX YY ZZ XY XZ YZ`. Water's XZ is −0.477, so the sum rule
  genuinely constrains it — but XY and YZ vanish by symmetry and stay unpinned.
  **An asymmetric molecule is needed to close that gap.**
- raw Cartesian, not traceless (trace −12.55)
- bond midpoints are plain arithmetic means

*Two risks the plan listed are now retired.* A default MAKEFP run on **linear HCN
adds no dummy point** — 3 atoms, 2 bond midpoints, nothing else. The shipped
`HCN.efp`'s `A04D4` comes from its own input, not from MAKEFP. And the partition
being nearest-Gaussian-product-centre winner-take-all is corroborated by the
electron count closing exactly.

*Done.* `mqc_libcint_dma.f90`, plus `validation/check_dma.{f90,py}`.

Stone's partition, at the primitive-pair level as it has to be: the product centre
depends on the two exponents, so different primitives of the *same* contracted
shell pair can land on different expansion points, and assigning at the shell-pair
level would give a different answer rather than a rounder one. So `uncontract`
builds a one-primitive-per-shell basis with unit stored coefficients, the density
is transformed into it, and the assignment and accumulation happen there. That
transform is the only real machinery in the module.

*Agreement with GAMESS, water/6-31G\*, all five points:*

| | |
|---|---|
| expansion points | 4.7e-11 Bohr |
| monopoles | 5.0e-9 |
| dipoles | 2.1e-9 |
| quadrupoles | 5.4e-9 |
| octopoles | 8.7e-9 |

which is the printed precision of the reference. The sum rule holds separately in
the Fortran, needing no reference: electrons recovered to 1e-12, and the total
dipole and second moment to 4e-16 and 6e-15.

**GAMESS uses the old CODATA Bohr, and it mattered.** 0.52917724924 Angstrom per
Bohr against the current 0.529177210903 -- 7e-8 relative, invisible in an energy.
Converting water's geometry with our constant put every coordinate 1.3e-7 Bohr from
GAMESS's, and that discrepancy would have been charged against the partition. The
tell was that the ratio was a *uniform* 1.00000007 across every coordinate, and the
implied constant matched 1/0.52917724924 to 1.2e-11. GAMESS-facing validation
programs now convert the way GAMESS does; PySCF-facing ones must not, since PySCF
is on the current value.

That fix also tightened M5's centroids from 4.6e-8 to **1.7e-9** while leaving its
tensors at 2e-5 -- which is the right way round, and confirms that the tensor
residual really is GAMESS's own CPHF precision rather than geometry.

*Size.* Medium, as estimated. The transposition risk did not materialise: the
component packing was already pinned by `check_sum_rules.py` before any code was
written, which is exactly what that script was for.

*One live gap.* HCN cannot be the M5 reference: **GAMESS's own CPHF fails to
converge it** (`TOO MANY ITERATIONS IN AOCPCG`, residual stalling near 3e-3
through its 299-iteration grace limit) and MAKEFP aborts before the
polarizabilities. Its diagnostic blames the wavefunction, but the gap is 0.685
Hartree and the orbital Hessian is positive definite with a preconditioned
condition number of 6.3 — our CG does it in 14 iterations. So pick M5 reference
molecules GAMESS can finish, and do not read a missing section as hard physics.

### M5 — Distributed (LMO) polarizabilities  *(GAMESS)*  ✅

Decompose M3's response onto M2's localized orbitals: one full 3×3 tensor per LMO
at its centroid.

*Settled from the water reference, all three checkable by `check_sum_rules.py`:*

- **The polarizable points are exactly the LMO centroids** — identical coordinates
  in both sections, to 0.
- **Per-LMO tensors are asymmetric** (0.17 for the O–H bonds, 0.06 for the lone
  pairs) and **their sum is symmetric to 4.6e-6**. So M5 must emit all nine
  components and must *not* symmetrize; only the sum rule may assume symmetry.
- **The core orbital is excluded.** Water has 5 occupied MOs and 4 polarizable
  points, and the four tensors sum to α_iso 4.865822 against 4.866587 for an exact
  all-occupied solve. The missing 7.7e-4 is diagonal-only — what a spherical 1s
  should contribute. So the sum rule is Σα_i = α(valence), *not* α(all occupied),
  and we need GAMESS's core-counting convention.

*Verify.* Σ_i α_i against the valence-only total from M3, exact by construction.
Then tensor by tensor against `POLARIZABLE POINTS`.

*Done.* `distributed_polarizability` in `mqc_libcint_cphf.f90`, plus
`validation/check_distributed_polarizability.{f90,py}` against the GAMESS
reference. The decomposition needed no new physics: the occupied index of
`alpha_kl = -4 sum_ai h^k_ai U^l_ai` is summed over and a rotation `W` among the
occupied orbitals is orthogonal, so inserting `W W^T` is free and each localized
orbital's share falls out. M2 did the real work.

*Agreement with GAMESS, water/6-31G\*:*

| | ours | GAMESS |
|---|---|---|
| centroids | — | **4.6e-8 Bohr apart** |
| per-tensor components | — | 2.0e-5 |
| alpha_iso of the sum | 4.865815 | 4.865822 |
| per-LMO asymmetry | 0.1721 / 0.1720 / 0.0609 | 0.1721 / 0.1720 / 0.0609 |

Our Boys optimizer lands on GAMESS's *exact* local maximum, which was not
guaranteed and is why the comparison pairs orbitals by centroid before comparing
tensors at all. The residual 2e-5 is GAMESS's own CPHF precision, not ours — their
valence sum carries 4.6e-6 of asymmetry where ours carries 3.7e-7.

**The transpose the plan warned about was real, and it was in the reference
convention.** The source-derived component order `xx yy zz xy xz yz yx zx zy` put
every tensor exactly one transpose from ours: the worst component gap came out at
1.72e-01, which *is* the tensor's own asymmetry — the signature of a transpose and
nothing else. Reading the two off-diagonal triples the other way round drops it to
2.0e-05. So GAMESS's labels are (field, dipole) where ours are (dipole, field), and
the correct order in our convention is `xx yy zz | yx zx zy | xy xz yz`. Nothing
else could have caught this: the sum over orbitals is symmetric, so M4's sum rule is
blind to it, as is every isotropic average, as is the asymmetry magnitude. **It
matters for M7 and M8** — an induced dipole computed from a transposed alpha is
wrong.

*Also settled:* the core is excluded from the **distribution** but kept in the
**response**. Freezing it out of the response too gives an exactly symmetric
valence sum (1e-15), which GAMESS's 4.6e-6 rules out; keeping it gives 3.7e-7, and
the isotropic value then matches GAMESS to 7e-6. That is the discriminating
measurement.

*Size.* Small once M2 and M3 exist — both did. **Reference molecules must be ones
GAMESS can converge**, which rules out HCN; see M4.

**Stop and ship here.** Electrostatics plus polarization is most of the
interaction energy for polar systems. Everything below is a second phase.

### M6 — Dynamic polarizabilities at twelve imaginary frequencies  *(GAMESS)*  ◐

**`DYNAMIC POLARIZABLE POINTS` done and validated**; the two quadrupole blocks are
not, and the obstacle is a definition rather than missing code.

`dynamic_polarizability` and `distributed_dynamic_cross` in `mqc_libcint_cphf.f90`,
`validation/check_dynamic.{f90,py}`. Imaginary frequency is what makes this cheap:
eliminating the antisymmetric part gives `[(A+B) + ν²(A−B)⁻¹] S = −2h`, a sum of
positive definite pieces, so there are no poles on the axis and CG still applies
with an inner CG against `A−B`. No explicit Hessian, nothing on disk — where GAMESS
needs GMRES against a stored Hessian because it solves the general real-ω case.

*Verified.* Quadrature points 4.8e-7 against the frequencies GAMESS stamps on its
own blocks; worst tensor element **8.2e-5** over 12 frequencies × 4 orbitals × 9
components, which is GAMESS's own precision. Two exact identities need no
reference: ω→0 reproduces `static_polarizability` to 1.1e-15, and the distributed
set at ω→0 reproduces M5's distributed tensors.

#### The dipole-quadrupole and quad-quad blocks: what is known and what is not

**Known.** The file stores 27 and 81 values per point, which are `3×9` and `9×9`
expansions of DAF records holding `18 = 3×6` and `36 = 6×6` — six *unique*
quadrupole components, duplicated into nine slots. Component order is
`XX YY ZZ XY XZ YZ`, from `prpel.src:5643`. And the operator is the **traceless**
Buckingham quadrupole, read off `prpel.src:5625` rather than guessed:

```
QXX = 0.5*(2*QMXX - QMYY - QMZZ)      QXY = 1.5*QMXY
```

i.e. `Θ_ab = (3 r_a r_b − r² δ_ab)/2`, with `QM` the raw second moments.

**Ruled out**, each by a full run over four orbitals and twelve frequencies and a
constant-ratio search over all 27×27 component pairings:

1. raw second moment `r_l r_m` as the driving operator;
2. the same referred to each orbital's own centroid — the physically motivated
   guess, since the dipole's occupied-virtual block is origin independent (`<a|i>=0`)
   while the quadrupole's is not. Testable exactly as
   `S_quad − R_l S_dip,m − R_m S_dip,l`;
3. the traceless operator above, from GAMESS's source;
4. the traceless operator referred to each orbital's own centroid — the combination
   of 2 and 3, which neither test alone covered. Computed exactly as
   `1.5(raw_ab − R_a dip_b − R_b dip_a) − 0.5 δ_ab Σ_c(raw_cc − 2R_c dip_c)`.

*Also verified, so it need not be re-checked:* both quadrupole blocks are
**frequency-major** — `CT1..CT4` then the next frequency, like the dip-dip block —
so the sample ordering used in the searches above was right. That mattered because a
mismatched sample order would produce exactly the symptom observed (no component
proportional to any other) even with a correct operator, since the constant-ratio
test compares one component's 48 samples against another's.

*And the tensor's shape is settled:* `NDIM=3, NRES=6, NPOL=18` at `locpol.src:1444`,
so GAMESS drives the response with the **dipole** and measures it with the six
traceless quadrupole components — the transpose of what our search generated, though
the response function is symmetric in its two operators so that alone cannot explain
the mismatch.

None of the 27 reference components is proportional to any of ours in any of the
three, so the difference is not a permutation and scale — which is what the
nine-component static polarizability and the Cartesian d shell both turned out to
be, and why those were eliminated first.

**Next.** Read the assembly rather than the operator: `POLDB` and `SOLVCPDYN` in
`locpol.src` are what build the tensor that reaches DAF 816, and the factor is
likely in there rather than in the integrals. Two constraints any candidate must
satisfy before an energy is computed: water is planar in `xz`, so every component
with an odd number of `y` indices must vanish, and `xy`/`yx` must hold identical
values in the nine-slot storage.


The same CPHF at ω on the imaginary axis, at the Casimir–Polder quadrature points.
Imaginary frequency is the friendly case: no poles, better conditioned than real ω.

*Verify.* The ω → 0 limit must reproduce M5 exactly — free and sharp. Then
`DYNAMIC POLARIZABLE POINTS`, then a water-dimer C6 against literature.

*Size.* Medium. The physics is a short step from M3; the bookkeeping of 12
frequencies × 3 perturbations × n_LMO is where mistakes will live.

### M7 — The JSON potential, and the renderer  *(GAMESS end-to-end)*

**M7a, Fortran:** emit the whole potential as JSON — projection basis, MO
coefficients, Fock matrix in the LMO basis, LMO centroids, plus M4 and M5. The
second, larger basis becomes real here, and it touches how a fragment calculation
is set up.

**M7b, Python:** the renderer, JSON to `.efp`.

*Verify.* Round-trip our JSON through our own reader. Then render and **feed it to
GAMESS**: if GAMESS accepts our potential and gets the energy it gets from its own,
the chain is validated by a program with no reason to be kind to us.

*Size.* M7a medium; M7b small, and testable against the eight shipped files with
no quantum chemistry at all — parse each, re-render, require a match modulo
whitespace. That task has zero prerequisites and can be done today.

### M8 — Stretch: use the potential

EFP–EFP, then EFP/*ab initio*. Only once M7 round-trips.

## Decisions before M1

0. **Nothing blocks the parser.** Parsing the eight shipped `.efp` files and
   re-rendering them needs no wavefunction, no integrals and no decisions at all.
   If you want to start today, start there — it pins the format and gives every
   later milestone something to compare against.
1. **Reference molecule set, basis and localization**, pinned in the comparison
   script. Suggest HCN and water in `6-31G*` with `LOCAL=BOYS` — small, fast, and
   one of them is linear so the dummy-point convention shows up immediately
   rather than at M4.
2. **Where the reference lives.** A generated manifest like
   `tools/cpu_validation/`, but holding `.efp` files produced by the local GAMESS
   instead of PySCF energies.
3. **CPU first**, for the same reason the DFT ladder went that way: the checks
   are cheap and the loop is seconds. cuEST's bound multipole entry point can
   follow.

## What is actually hard

Honest, worst first.

- **The DMA partition and the dummy-point convention (M4).** Not conceptually
  hard, but matching GAMESS's specific choices is easy to get subtly wrong in a
  way that satisfies the sum rule and still disagrees with the reference.
- **Two basis sets (M7).** A structural wrinkle rather than a difficulty, but it
  touches how a fragment calculation is set up, so better designed for early
  than retrofitted.
- **Twelve-frequency bookkeeping (M6).** Nothing deep, many indices.
- **Nothing about CPHF.** A preconditioned linear solve with a well-understood
  failure mode, verifiable by finite difference from the first iteration.

## Explicitly out of scope

MCSCF (EFP2 is HF-based), SAPT (wants M3's CPHF, so strictly cheaper afterwards),
UMP2 (unrelated; do it when someone hits the refusal), and gradients — which
share nothing with this chain, since SCF gradients need no response equations.
A parallel track, not a step toward MAKEFP.
