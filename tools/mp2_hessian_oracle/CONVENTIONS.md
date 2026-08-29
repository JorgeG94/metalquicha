# Reading pycc's numbers in metalquicha's conventions

Unit 0.3 of the MP2 Hessian ladder (`mp2-hessian-phased-plan.md`). Every gate on
that ladder compares one of our arrays against one of pycc's, and the two codes
do not spell the same tensor the same way. Nothing here is deep; all of it is
the kind of thing that costs a night when it is discovered by a failing gate
instead of read beforehand.

Each claim below cites where it was checked. Anything not yet measured says so.

## 1. Index order: chemist and physicist both appear, in the same file

* **pycc's `Gam` (the 2-PDM) is physicist**, `<ab|cd>`. `_effective_2pdm_ao`
  converts it with `.swapaxes(1, 2)` before doing anything else
  (`correlatedderivs.py:1015`).
* **pycc's `erix`/`ao_eri2` are chemist**, `(μν|λσ)`.
* **Our `gamma_ao` is chemist**, the same slot layout as the AO ERI. Fixed by
  its consumer: `contract_gamma_eri` computes
  `Imat(p,q) = Σ_{u,r,s} (u p|r s) γ(u,q,r,s)`
  (`mqc_libcint_mp2_gradient.f90:1068`) — γ's four slots line up with the
  integral's, so it is `(μν|λσ)`-ordered, not `<μν|λσ>`.

So our γ and pycc's `Gam` differ by the `1↔2` swap before any other difference
is considered.

## 2. The two ERI spellings

`(μν|λσ)` (chemist, charge-density pairs) is what the spatial/closed-shell path
uses on both sides. `<pq||rs>` (physicist, antisymmetrised) appears only in
pycc's spin-orbital branch, which we do not implement — if a formula you are
transcribing has a `||` in it, you are reading the wrong branch. The
closed-shell companion of `<pq||rs>` is
`L = 2(μν|λσ) − (μν|σλ)`, which pycc build per perturbation as `wx`.

## 3. pycc's naming rule

**A trailing `X` or `x` marks a skeleton *integral* derivative — never a
density derivative.** `fX` is `f^(X)`, `SX` is `S^(X)`, `erix` is
`(μν|λσ)^(X)`. Density derivatives are spelled `d…`: `dDrel`, `dGam`, `dI`.
Stated in `_hessian_blocks`' docstring; worth re-reading because `Xx` and `I2x`
break the eye — they are the *augmented skeleton carriers*, still integral-side.

## 4. The `_effective_2pdm_ao` transpose chain

Four operations, in this order (`correlatedderivs.py:1014-1021`):

```python
GamD = 2.0 * einsum('ac,bd->abcd', Drel, Pocc) - einsum('ad,bc->abcd', Drel, Pocc)
G = (Gam + GamD).swapaxes(1, 2)          # 1. physicist -> chemist
G = 0.5 * (G + G.transpose(2, 3, 0, 1))  # 2. bra<->ket symmetrise (chemist)
G = ...back-transform all four indices with C...   # 3. MO -> AO
return G.transpose(0, 1, 3, 2)           # 4. invert mo_eri_helper's ket reorder
```

Steps 2 and 4 exist **to match psi4's raw `ao_tei_deriv2` layout**, not because
the physics needs them — their docstring says so outright ("so Psi4's raw,
bra↔ket-doubled `ao_tei_deriv2` is exact without `_complete_deriv2`", "invert
`mo_eri_helper`'s internal reorder").

**Consequence for Gate 1.2a — settled 2026-08-27.** The dumped `GeffAO` is a
psi4-convention object, and the comparison is against
`GeffAO.transpose(0, 1, 3, 2)`: step 4 undone. Measured, our
`build_effective_2pdm_ao` against the dump — raw **8.4e-2**, ket-swapped
**2.09e-11**. libcint hands the ket pair in the natural order, so our builder
must not carry that step, and it does not. Step 2, the bra↔ket symmetrisation, is *not* optional on our side
either: it is what makes the density contract correctly against raw derivative
integrals with no completion step.

## 4a. The three reconciliations, measured

Unit 1.1. Each of the gradient's three quantities is the same physics as pycc's
counterpart and none of them is the same array. Measured 2026-08-27 on
water/6-31G all-electron, at one thread.

**Relaxed density** — direct, after aligning the MO phase (a sign per orbital,
§7): `relaxed_density_mo` against `Drel`, **1.75e-11**.

**Energy-weighted density.** `imat` is *not* the Lagrangian, and comparing it to
pycc's `I` is comparing a part to a whole — its own assembly comment says the
two-particle density's share is there and the unrelaxed one-particle correction
is in `veff`. What both codes mean is the matrix multiplying `dS/dR`, exposed as
`energy_weighted_ao`:

```
W_ours(AO) = 2 C I_pycc C^T  -  4 sum_{i in occ} eps_i C_i C_i^T     (1.92e-10)
```

The second term because our gradient returns the *total*: our W carries the
reference's energy-weighted density too, which is most of its norm (81.98
against 0.33 for the correlation part). Compared without it, the disagreement is
a factor of 245, which reads like a broken routine rather than a different
decomposition.

**Two-particle density.** The two differ in index *placement*, not just
convention. pycc's `_mp2_tpdm` is `Gamma_ijab = 2 t2_ijab - t2_ijba` in the
`oovv`/`vvoo` blocks — order `(o,o,v,v)`. Ours carries the energy's coefficient
`4 t2_ijab - 2 t2_ijba` = `2 u` and back-transforms with
`C_occ, C_vir, C_vir, C_occ` — order `(o,v,v,o)`. So in chemist pairs ours is
`(ia|bj)` where pycc's is `(ia|jb)`: the **ket pair is reversed**. Then ours is
symmetrised `1<->2` (summed, not averaged) where pycc's effective density is
symmetrised bra<->ket. With `u = 2 t2 - t2.transpose(0,1,3,2)`:

```
T[m,n,l,k] = sum_ijab 2 u_ijab C_o[m,i] C_v[n,a] C_v[l,b] C_o[k,j]
gamma_ours = T + T.transpose(1,0,2,3)                                (9.27e-12)
```

**Why this matters beyond the gate.** The `1<->2` symmetrisation is there because
the first-order contraction exploits `(grad i j|k l) = (grad i j|l k)`; pycc's
bra<->ket one is there so the density contracts against a bra<->ket-doubled
*second*-derivative integral. They are not interchangeable. `Gamma_eff` for the
Hessian has to be built with the bra<->ket symmetrisation rather than inherited
from `gamma_ao`, which is the plan's Unit 1.2 and is more work than reusing what
the gradient already builds.

## 5. Storage and architecture

pycc hold `nmo⁴` MO-basis tensors and stream the per-perturbation ones from a
disk store. We are AO-blocked and direct. The one place this matters for a gate:
our gradient takes a **dense** `nao⁴` path whenever `2·nao⁴·8 ≤ IN_CORE_LIMIT`
(4e9, about 118 basis functions; `mqc_libcint_mp2_gradient.f90:223`) and blocks
only above it. Water/6-31G is far below, so every Phase 1 gate runs on the dense
path and the blocked path is gated separately, later.

## 6. AO index map, psi4 ↔ libcint

**Our conventions** (`mqc_libcint_ao.f90:10-30`, which documents them as three
ways to be wrong quietly):

* normalisation is folded into `env` by `molecule_build`; nothing may normalise
  again, or the overlap diagonal is not 1;
* **Cartesian component order is libcint's**: `x^(l-i) y^(i-j) z^j` over
  `i = 0..l`, `j = 0..i` — for d that is `xx, xy, xz, yy, yz, zz`;
* the spherical transform is libcint's own table
  (`mqc_libcint_ao_data.f90`), whose `l = 1` block is the px/py/pz identity.
  That file's header records that libcint's C table carries **two p orderings**
  behind `#ifdef PYPZPX`, so anything scraped from it must say which it took.

For the pinned case this reduces to shell ordering plus p-component order:
water/6-31G has s and p shells only, and no d, so the Cartesian-vs-spherical
question does not arise at all.

**Measured 2026-08-27: the map is the identity.** Comparing the AO overlap
matrix at the pinned geometry — `mol%overlap(s)` on our side against
`MintsHelper.ao_overlap()` on psi4's, both fed `basis_sets/6-31g.json` —

```
nao 13 both;   max |S_ours - S_psi4| = 5.55e-16     (identity ordering)
```

so no permutation and no sign. AO-basis dumps can be compared element-wise as
they stand.

**Why, and where it stops being true.** Four things happen to agree here, and
only the first is a general fact:

* both expand an `SP` shell into an S then a P in file order;
* **no d functions**, so Cartesian-versus-spherical never arises — 6-31G on
  water is s and p only;
* Cartesian `p` is `x, y, z` on both sides;
* normalisation matches, which the unit overlap diagonal confirms on both.

### 6-31G\*: the ordering survives d functions, the normalisation does not

Measured the same way, 2026-08-27. 6-31G\* puts a Cartesian d shell on oxygen
(`function_type: gto_cartesian`, so 6d), nao = 19 on both sides.

* **Shell layout is identical** — psi4 builds S,S,P,S,P,D on O and S,S on each
  H, which is our layout exactly; the SP shells expand in file order on both.
* **Cartesian component order is identical.** Both overlap diagonals carry the
  `1, 1/3, 1/3, 1, 1/3, 1` pattern position for position, which is
  `xx, xy, xz, yy, yz, zz` in both.
* **The d block is scaled.** Our d functions are `sqrt(4*pi/5)` = 1.58533092
  times psi4's — libcint does not normalise individual Cartesian components,
  psi4 does. Raw, `max |S_ours - S_psi4| = 1.51`; after a diagonal rescale
  `N = diag(1..1, sqrt(4pi/5) x6, 1..1)`, `max |S_ours - N S_psi4 N| = 6.66e-16`.

**This is a convention, not a bug.** A uniform scaling of a subset of basis
functions is a change of basis and the SCF is invariant under it, which is why
our 6-31G\* energies validate against PySCF — PySCF is libcint-based and shares
our convention. psi4 is the outlier.

**The rule for AO-basis gates.** Comparing an AO-basis dump element-wise
(`GeffAO`, `ao_eri2`) for any basis with d or higher needs the diagonal rescale
applied first — a permutation is not the correction, a similarity transform by
`N` is, and for a rank-4 object that is one factor of `N` per index. The pattern
looks like `sqrt(4*pi/(2l+1))` per shell, but only `l = 2` has been measured;
f has not been checked (`6-31g(2df,p)` is extracted and would settle it).

None of this touches Phases 1 and 2, which are pinned to 6-31G and have no d
functions at all. It is written down because the first basis upgrade will hit
it, and a factor of 1.585 on part of a tensor looks far more like a bug in a new
Hessian contraction than like a basis convention.

The measurement was made with a throwaway dump, deliberately not committed as a
test: `test_mqc_hess_ints`' own header records why — a stored comparison against
an external library pins our layout to that library's conventions rather than to
anything true. The recipe is two lines (`build_libcint_molecule` then
`mol%overlap`), so re-running it costs less than maintaining it.

**This does not settle Gate 1.2a.** The AO *ordering* being the identity is a
different question from the ket-swap in §4, which is a tensor-layout convention
inside `_effective_2pdm_ao` and still has to be decided there.

## 7. MO phase alignment

An SCF eigenvector is defined up to a sign, and the two codes may land on
opposite ones. Before any **MO-basis** element-wise gate, fix a phase per
orbital — e.g. force the largest-magnitude coefficient of each column positive
in both codes, or align column by column against the dump.

C1 water/6-31G has no degeneracies, so a sign per column is the whole job. A
degenerate pair would need a rotation rather than a sign, which is a reason to
keep the pinned case non-degenerate rather than to solve the general problem.

## 8. Gauge

`perturbed_mo_gauge` is `'non-canonical'` for MP2 (`'canonical'` only for
CCSD(T)). The dumps record it in the manifest. Perturbed intermediates compared
across a gauge difference disagree *while both are right*, which is
indistinguishable from a bug by inspection — this is also why perturbed
quantities are not finite-difference-checkable (design plan, §Validation).

## 9. Frozen core

`P` in the `Γ_D` fold is the projector onto the **full** occupied space, core
plus active. pycc build it as `arange(nmo) < wfn.o.stop`
(`correlatedderivs.py:1012-1013`) — `o.stop`, not `no`.
