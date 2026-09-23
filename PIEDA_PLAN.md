# Pair interaction energy decomposition for a fragmented protein: a plan

Per-fragment-pair interaction energies, decomposed into electrostatics,
exchange repulsion, charge transfer and dispersion, over a protein fragmented
by residue. Built from three surveys run on 2026-09-22: our FMO and EFMO
paths, our EFP2 terms, and GAMESS's PIEDA.

**Baseline: `origin/main` at `31e4e0bd73`.** The FMO/EFMO and GAMESS surveys
were run against main and are sound. The EFP survey was run against a working
tree 308 commits behind; its physics is expression-level and stands, but
**every file:line it gives must be re-checked before use.** Do not plan from
`QUAO_PROTEIN_PLAN.md`, written the same day from the same stale tree.

Not committed on purpose, like the other plan files.

---

## The finding that sets the shape

**EFMO already computes four of PIEDA's five buckets, per pair, under those
exact names.** `efmo_pair_t` on main carries `electrostatics`, `dispersion`,
`exchange_repulsion` and `charge_transfer` for every far pair, filled from
`efp_pair_terms`, and prints them in a per-pair table at info level. And
**FMO already computes the undecomposed pair interaction** — after
`subtract_subsets`, a two-member term's correction is exactly
`E'_IJ − E'_I − E'_J + Tr(ΔD_IJ u_IJ)`, the quantity PIEDA decomposes. It is
summed into one scalar and the locals go out of scope.

So the first deliverable is not a calculation. It is keeping numbers that are
already correct and already computed.

Three layers each throw data away: the backend result type has no pairs array,
the bridge narrows FMO to a single real, and the reporter writes a total.
`fmo_report` declares MBE output mode with nothing allocated behind it, so an
FMO user gets one JSON number and an empty levels array.

---

## Decisions

1. **EFMO is the host, not FMO.** It already has the pair type, the four
   terms and the table. FMO gets the undecomposed pair energy exported
   (cheap, useful on its own) but the decomposition lands in EFMO.

2. **Report unconnected pairs by default and refuse to print a bare connected
   one — an FMO-side concern only.** EFMO has no bond-breaking option at all,
   so no EFMO fragment is half a severed bond and no row of its table can be
   the covalently-joined case. The same reason means the monopole problem of
   Layer 0 does not reach an EFMO export today. Every adjacent-residue pair is covalently joined, and an uncorrected
   connected pair reads about −14 Hartree in a table of kcal/mol. GAMESS
   prints a separate unconnected-only block for exactly this reason. For
   ligand binding, which is what people actually do with this, nothing of
   interest is connected and the default costs nothing.

3. **The cut-bond reference is a repartition, not a correction, and the
   invariant is the test.** Subtracting a per-bond reference from each
   component and moving it into the one-body sum must leave the total energy
   **bit-identical** to the same run with the analysis off. That is the one
   cheap regression the whole feature can be held to, and it should be a test
   before it is a feature.

4. **Charge transfer is the residual on the FMO side, and must be named so
   there — but not in EFMO.** EFMO's charge transfer is an explicit
   perturbative sum over occupied-on-A to virtual-on-B amplitudes, computed
   independently of the other three terms, so it carries none of the
   basis-set superposition error a residual absorbs and may be exported
   under its own name. What follows applies to the FMO-side term only. GAMESS defines it
   as what is left after the other three, which is why it is called
   `Ect+mix`. Basis-set superposition error and any unaccounted coupling land
   there. Report it under a name that says so; never as "charge transfer"
   unqualified.

5. **Dispersion is the correlation contribution to the pair interaction**, and
   is zero at Hartree-Fock. EFMO already correlates fragments and already
   reports damped dispersion per far pair, so this comes free there and is a
   bounded addition to FMO only if someone wants it.

6. **There is no cap in the potential generation, and no double counting.**
   Two earlier drafts of this decision were wrong in opposite directions.
   GAMESS adds no atom for a cut fragment's effective potential: it splits
   the **nucleus**. The detached atom keeps Z−1 in its own fragment and loses
   an electron; an extra centre at its exact coordinates carrying **+1 and
   the detached atom's full heavy-atom basis** joins the neighbour, which
   gains an electron. Both fragments stay neutral closed shells, so an
   ordinary parameterisation run on them is well defined, and summing the two
   fragments' multipoles reconstructs the full nucleus exactly once. The +1
   centre is a cap in every respect except that it was not moved to 1.09 Å
   and not given a hydrogen basis. Real hydrogen caps do exist in GAMESS, but
   in the auxiliary model system used to manufacture the hybrid, in the free
   monomers of the reference state, and in a QM/MM path — never in the
   fragment itself. The widespread reading that this method needs no caps is
   true of the fragment and false of the model system; the documentation does
   not actually claim otherwise.

   The bond region is therefore described once, not twice. The frozen orbital
   is a constraint inside one fragment's SCF and removes no charge; the
   potential is the expansion of what that SCF produced. And a covalently
   connected pair never uses the classical term at all — it is computed
   quantum mechanically, where the nuclear charges add back and the bond is
   whole.

7. **Polarization is not reported per pair.** The induction is solved
   self-consistently over the whole system by construction, and a per-pair
   number would be a new definition with two defensible answers. EFMO already
   subtracts a pair polarization for its own energy; that is a different
   quantity and should not be relabelled. Leave it out of the table and say
   why.

8. **Energy runs only.** GAMESS refuses gradients for ab initio PIEDA. Follow
   that rather than discovering it.

---

## The layers

### Layer 0: make a cut fragment neutral

**For EFMO alone. It does nothing for FMO, and an earlier draft of this
section claimed otherwise.** Implemented and measured on
`feat/fmo-neutral-cut`.

Our fragments carried net charges of plus and minus one at a cut, because we
kept whole nuclei and shifted only electrons. GAMESS splits the nucleus so
both sides stay neutral, and we now do too.

**The FMO accuracy argument was wrong, and is disproved rather than merely
doubted.** With point-charge embedding the total is *exactly invariant* to
the convention: the field a group feels is the whole-system charge minus its
own share, so a unit of charge taken out of a group's own nucleus reappears
in the field it sees at the same point, and the summed populations do not
change either. Measured across the convention change: 0.48908757203956554
against 0.48908757203680864 on propane, agreeing to convergence noise. There
was never an accuracy gain here to collect.

**Field-free runs get worse, so the split is gated on there being a field.**
With nothing to supply the other half of the nucleus, the owning fragment is
solved around a nucleus short by a proton, which is a worse model of a methyl
group than the cation it replaces: 0.180 to 0.304 on propane, and 0.125 to
1.549 on a numbering where the middle carbon is detached twice. GAMESS cannot
arbitrate this because it never runs the construction without a field; its
field-free reference state uses methyl caps instead. So whole nuclei are kept
when `esp = "none"` and the nucleus is split otherwise, with both conventions
stated in the log rather than switched silently.

**Why it matters for EFMO.** The far-pair terms are EFP2 expressions between
two fragments' own potentials, and there is no subtraction of a group's own
share in a pair term, so a fragment's net charge simply sits in its multipole
expansion. A charge on each of two adjacent residues gives a monopole term of
roughly 96 kcal/mol at 3 to 4 Angstrom, against interactions of single-digit
kcal/mol, plus monopole-dipole and monopole-induced-dipole terms — and it
lands in the electrostatics column of the table this plan exists to produce.
The neutral convention removes all of it by construction.

*Gate:* not the telescoping identity, which is blind to the field, and not
the two-body number, which is invariant. The gates are per-fragment
neutrality, and an explicit test that the embedded total agrees under both
conventions, which turns the invariance from a derivation into something CI
enforces.

*Also established:* there is no three-electron lone-pair case to handle. A
bond carrying more than one localized orbital is refused by name, and a
hydrogen detached atom likewise, so every cut moves exactly one electron and
one unit of charge. Several separate cuts at one atom compose correctly and
were exercised.

### Layer 1: covalent fragments

The prerequisite. EFMO declares covalent fragments to be Phase 5 and not
started. A protein cannot be fragmented by residue without cutting the
backbone, so nothing downstream matters until this exists.

**It is much smaller than this plan first assumed, for energies.** Measured
in GAMESS: the Fock-zeroing core is twelve lines, the charge split and the
merge-back-inside-an-n-mer rule are about thirty each, the side-detection and
frozen-vector assembly around eighty, and **EFMO's entire covalent-specific
contribution is about sixty lines** — it reuses FMO's machinery unchanged and
only forces the generalized variant on. What is genuinely large is the
analytic gradient, about 1800 lines, because the model system's caps move
with the real atoms and the constraint needs its own response solve. Energies
first; budget separately for derivatives.

The zeroing is per iteration, in a **fixed** orthonormal basis stored once
rather than the current orbitals, applied after the Fock build and before
diagonalisation, with the frozen-occupied block left intact, every cross
block zeroed, and the frozen-virtual diagonal pinned to a large shift. The
energy is then repaired by half the change in the one-electron trace. Zeroing
inside the frozen-occupied block, or leaving the frozen-virtual couplings
alive, are both different and worse methods; the second is fine for a
Hartree-Fock energy and silently corrupts a correlated one.

For FMO the neighbouring problem is bond breaking under embedding, refused
because the detached atom's share of the field is not subtracted. That refusal
is precisely scoped: clean for point charges, undefined for an exact density.
The point-charge case is a mask change rather than a subtraction — mark the
ghosted detached atom as inside the group so its term is never built — plus
the awkward part, widening the monomer density layout to keep the ghost block,
which is the surface that makes a parallel answer rank-independent.

*Gates:* at expansion level equal to the fragment count the expansion
telescopes to the supermolecular energy whatever the partition did, now with
point charges and frozen orbitals together. **Done, on `feat/fmo-afo-embedding`.**

**But that identity is a weaker gate than this plan first claimed.** At full
expansion order the top group is the whole system and has nothing outside
it, so there is no field at all, and every intermediate group's embedding
cancels out of the sum. It catches assembly faults, ghost and size
mismatches and crashes — it did catch a live size mismatch — but it is blind
to the field, to the charge convention, and to any embedding error. Pair it
with a charge-sum assertion, which is a physics check, and with reported
two-body numbers, which are not bounded by anything the identity says.

A ring cannot be used as a second gate molecule: every partition of a
three-ring severs two bonds between the same pair of fragments, which is
refused by name as a ring cut. Methylcyclopropane, cut at the exocyclic
bond, is the reachable case.

*Hazards recorded by the survey:* a live size mismatch where fragment charges
are computed against a ghosted molecule with a truncated density, unreachable
only because of the current refusal and reachable the moment it lifts; and the
detached atom's population arriving from two fragments with only one writer.

### Layer 2: export the pair energies that already exist

FMO keeps a pairs array on its result and the bridge and reporter carry it
through to JSON. EFMO's per-pair table reaches JSON the same way, with the
four terms per far pair. Purely additive, low risk, and on its own it gives a
per-pair interaction map of a protein with distances, which is most of what a
consumer of this actually reads.

*Gate, corrected against measurement:* the sum of the exported pairs
reproduces the **interaction** energy to round-off (1.4e-17 on the prism,
exactly zero on the cage), **not** the total — the monomer sum and the
polarization total are not pair quantities. And above level 2 the pair map
is a pair map while the near expansion is not: pairs are filled only for
two-member terms, so a level-3 deck falls short by exactly the difference
between its third-level vacuum and induction contributions, measured at
1.26e-4. Assert that shortfall rather than working around it. The second
half of the gate does hold unchanged: the total must not move.

### Layer 3: decompose the QM pairs

Far pairs already have four terms; close pairs have only the undecomposed
interaction. Electrostatics between the two unperturbed monomer densities uses
machinery that exists on the FMO path. Exchange repulsion needs the
Heitler-London reference, which means the dimer SCF must start from
orthonormalised monomer orbitals and its first-iteration energy captured
before the embedding trace is added back. Charge transfer is then the
residual.

*The structural blocker:* fragments deliberately keep only a density, because
orbitals must not cross a wire, and the exchange term needs orbitals. Either
the carrier grows, with an MPI cost on the layout that guarantees
rank-independence, or fragments are re-solved for the analysis pass. Decide
this before writing any of Layer 3.

*Gate:* the four terms sum to the pair interaction exactly, by construction,
and the identity is asserted.

### Layer 4: the cut-bond reference

Generate a per-bond reference from a two-fragment model molecule at the real
bond length, subtract it component by component, move it into the one-body
sum. Ethane split across its carbon-carbon bond is the model for a backbone
cut.

*Gate:* decision 3's invariant, and the model molecule and bond lengths named
in the output, because the numbers are model-dependent and a reader must be
able to reproduce them.

*Known limit to exceed:* GAMESS gives up when an atom carries two cut bonds,
and says so in a warning. A backbone cut per residue puts two cuts on many
atoms immediately.

### Layer 5: reporting

A per-pair table sorted by interaction strength, which nothing in the tree
does today, and a JSON section. The intrinsic decomposition's pair block and
SAPT's named terms are the two patterns to copy; the EFP path is the one not
to, since it flattens its breakdown to an unnamed six-vector and exports only
the total.

---

## What this costs

GAMESS's single-run PIEDA is 1.0 to 1.2 times a plain FMO2 run: one extra
electrostatic-only evaluation per close pair, the first-iteration energy
captured for free, and arithmetic. No second SCF, no second integral
transformation, no change in scaling. The four-job variant that also
decomposes polarization is twice the dimer work plus manual plumbing between
runs, and buys interpretation that is ill-defined once a covalent bond is cut.
Do the cheap one.

---

## Explicit non-goals

- The free-state decomposition and its polarization terms. Four coupled runs,
  and the free state is ill-defined the moment a bond is cut.
- Gradients.
- A per-pair polarization number, per decision 7.
- Solvation terms. Continuum solvation exists but has no coupling to any
  fragment decomposition.

---

## Hazards

- **This repository's working directory is far behind main.** Re-baseline
  before reading any file reference in an older plan, including this one's
  EFP citations.
- The intrinsic decomposition and the quasi-atomic analysis are a different
  question and a different quantity; do not let them merge into this.
- Combining the bonding analysis with fragmentation is still unguarded on
  main: it runs per fragment, on capped geometries, and writes nothing.
- No output anywhere is sorted by interaction strength, and no screening
  report names a pair.
- The many-body deltas are copied to JSON only when not computing a Hessian,
  so a frequency run discards every one of them.

## References

Fedorov and Kitaura's PIEDA, as implemented in GAMESS. `mqc_docs/source/fmo.rst`
and `efmo.rst` on main. The EFP2 term implementations, whose expressions
transfer to ab initio fragment pairs without a potential file.
