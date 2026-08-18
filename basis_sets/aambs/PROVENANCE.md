# AAMBS — accurate atomic minimal basis set

`aambs.json` is the free-atom minimal basis that the quasi-atomic orbital (QUAO)
analysis projects onto. It is **tracked**, unlike everything else under
`basis_sets/`, which is extracted from the Basis Set Exchange bundle at configure
time. That is also why it lives in a subdirectory: `mqc_extract_basis_sets()`
globs `basis_sets/*.json` and deletes anything `MQC_BASIS_SETS` does not name, and
`.gitignore` ignores the same pattern. Neither reaches into here.

## Where it came from

Transcribed from GAMESS, `source/modf77_aambs.src`, `module aambsnorel`
(lines 10943–15353). The data is tabulated there as Fortran `DATA` statements —
nothing is computed at runtime — so this is transcription rather than
regeneration. Element and shell bookkeeping comes from `AAMBS_ATOM`
(`aambs.src`), `ETSHELL` (`vvos.src`), `NVVOS_NUMCOR` (`vvos.src:1166`) and
`LOCAL_NUMVAL` (`locsvd.src:1148`).

Upstream provenance, from the comment block at `modf77_aambs.src:11501`:

- **Exponents, H–Ar**: even-tempered, Schmidt & Ruedenberg, *J. Chem. Phys.* **71**,
  3951 (1979). Verified here — carbon's s-exponent ratios are constant at
  2.578990 to 1e-9.
- **Exponents, K onwards**: well-tempered, Huzinaga & Klobukowski, *Chem. Phys.
  Lett.* **212**, 260 (1993); Huzinaga & Miguel, *CPL* **175**, 289 (1990).
- **Coefficients**: spherically-averaged free-atom calculations, method varying by
  element — ROHF (H, Li, N, Na, P, K, Rb), RHF (He, Be, Ne, Mg, Ar, Ca, Kr, Sr,
  Xe), GVB (B, C, O, F, Al, Si, S, Cl), ORMAS (Ga–Br, In–I), and state-averaged
  MCSCF for the transition metals.

The exponents were generated on the fly until 2022, when they were frozen into
`DATA` arrays; some coefficients were then recomputed against the fixed-precision
exponents. The superseded values survive as 14 commented-out `DATA` blocks, which
the extraction skips — taking them would silently select the wrong Cl/Ar 3p data.

## Format

Basis Set Exchange schema, with two additions per element that the algorithm
cannot run without:

- `n` on each shell — the principal quantum number. BSE shells are unlabelled,
  but the core/valence split downstream is index arithmetic over shells in a
  specific order, so the label has to survive.
- `n_core_orbitals` / `n_valence_orbitals` — *chemical* core (semicore counted as
  core, unlike an ordinary frozen-core count) and the valence remainder.

This is a **general contraction**: per element and angular momentum there is one
exponent set, and every shell of that L is a separate contraction over it. Carbon
is `1s`, `2s` (14 shared s primitives) and `2p` (7 p primitives).

**Shell order is not simply by (n, l), and it matters.** Sc–Zn are ordered
`1s 2s 2p 3s 3p 4s 3d` — 4s *before* 3d — while Ga–Kr are `... 3p 3d 4s 4p`. The
swap exists so that each atom's chemical core is always a contiguous prefix of
that atom's block, which is the invariant the core/valence split relies on.

## What was checked

All 54 elements (Z = 1–54), 376 shells:

- every contracted shell normalizes to 1 — worst deviation 2.5e-6, on Ru 3p,
  which GAMESS itself special-cases with a relaxed tolerance;
- `n_core + n_valence` equals the total spherical orbital count;
- the derived valence count agrees with `LOCAL_NUMVAL` for every element;
- `n_core` falls on a shell boundary, i.e. the core really is a whole number of
  leading shells.

## Two things to know before trusting it

**Zn and Xe are not internally orthogonal.** These are free-atom SCF orbitals and
same-L shells on one atom should be orthonormal. Almost all are, to ~1e-6. But
the outermost valence shell of two elements is not:

| pair | overlap |
|---|---|
| Zn 1s–4s | 0.638 |
| Zn 3s–4s | 0.534 |
| Xe 2p–5p | 0.627 |
| Xe 4p–5p | 0.449 |

Both are the last element of their block, which smells like a regeneration slip
upstream. It is transcribed verbatim rather than quietly orthogonalized, because
a "fixed" file would silently disagree with the code it came from. It matters
because Paper I's construction assumes the AAMBS are orthonormal within an atom
(the SVD alignment theorem needs both bases orthonormal). GAMESS's VVO path is
immune — it Löwdin-orthogonalizes the whole AAMBS overlap including intra-atomic
blocks — but its QUAO path reads the un-rotated overlap and is not. **Anything
touching Zn or Xe should be treated as suspect until this is resolved upstream.**

**The relativistic set is not here yet.** `module aambsrel` covers Z = 1–86 and
extracts cleanly (886 shells, all norms within 8e-7), but the lanthanide shell
order is not `(n, l)`: for Z = 57–70 the core count of 27 only lands on a shell
boundary if 5s and 5p precede 4f, mirroring the 4s-before-3d rule. That branch of
`ETSHELL` has not been read carefully enough to commit the file. It is deferred
regardless — it needs a relativistic Hamiltonian, and this code has none.

## Coverage

Z = 1–54. Beyond that, GAMESS silently abandons VVO generation (as it does when
ECPs are present); we should refuse loudly instead.
