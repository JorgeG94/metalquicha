# Properties and analyses

Excited states, NMR shieldings and couplings, Raman activities, and the
analyses run on a converged wave function. Read from `AGENTS.md`, which says when.

## Excited states

Excited states are not a method name and not a driver. `keywords.excited_states`
on an ordinary energy deck gives linear-response excitation energies -- `tda`
or `rpa`, singlets, triplets or `both` out of a closed shell, and unrestricted
roots out of an open one -- over Hartree-Fock or Kohn-Sham up to a
range-separated hybrid, with transition dipoles and oscillator strengths in
both gauges and natural transition orbital weights. Amplitudes come back at
`sum X^2 - sum Y^2 = 1/2` restricted and `sum_sigma (...) = 1` unrestricted,
and the factor in a transition moment moves with that convention rather than
separately. meta-GGA, VV10, density fitting, PCM, capped fragments, a
correlated reference and cuEST are refused by name rather than approximated.
See `mqc_docs/source/excited_states.rst`.

## NMR shieldings

NMR shieldings are `properties.nmr` on an energy deck: the shielding tensor
of every nucleus over gauge-including (London) orbitals, closed-shell
Hartree-Fock or Kohn-Sham (LDA, GGA, global hybrids) with exact integrals,
everything else refused by name. Two
modules. `mqc_czt_giao` is the integrals, through the field-derivative
entry points `libcint_fortran` exports from libfint v0.9.0 on
(`libcint_1e_igovlp_*` through `libcint_2e_ig1_*`). `mqc_czt_nmr` is the physics, term
for term `pyscf.prop.nmr.rhf`: the field is an imaginary antisymmetric
perturbation, so the coupled-perturbed operator is exchange only
(`response_mean_field` with `minus`), and the occupied-occupied block of the
first-order coefficients is fixed at `-s1/2`. **The tensor is not
symmetric**, first index the field and second the nuclear moment, and a
transposed result passes any symmetric check. The reference is
`pyscf.prop.nmr.RHF` from the separate `pyscf-properties` package, whose
packaged response callback assumes three vectors a call and fails under
current PySCF; pass it one that reshapes to any count.
**The two-electron field derivative is folded and screened**
(`giao_coulomb_exchange`): the gauge factor of a pair is `(1/2)(R_i - R_j) x r`,
so the integral vanishes on one centre, changes sign under `i <-> j` and keeps
the ordinary symmetry in `k, l` -- a quarter of the quartets at most, each put
in the four places the full sum would. The screen's estimate of that factor
grows with the distance from the origin, so **test it with the molecule moved
off the origin**; `test_folding` holds the routine to the full sum written
out, to 1e-11, and `test_screen` does the same on three waters seven
angstroms apart, where the screen skips quartets and says how many -- on one
water it skips none, and a test there guards nothing. **The field response builds its coupling matrix once**
(`exchange_coupling`, through `build_hessian_mo`): the operator's two-electron
part is `(A-B)` less its diagonal, scaled by the exact-exchange fraction, and
`mo_transform_fits` decides whether it fits, per process, against a
budget that `set_memory_sharers` divides among the ranks on the node; otherwise every pass recomputes
the integrals (`integral_direct` forces that, and a test holds the two routes
together). A profile of a pass is all libfint -- the digestion is a tenth of
it -- so nothing in the direct route is left to tune but the pass count. Two
things that did not work: one Krylov subspace shared by the three directions
took the same fourteen passes as three separate ones, and applying the
operator to each new vector at its own shrinking length so the build could
screen on it saved a third and then **failed to converge** on triglycine,
its noise sitting at the 1e-10 tolerance. Measured on four threads, the
shielding step of a water hexamer in cc-pVDZ went from 27.7 s to 8.8 s and
triglycine's in 6-31G* from 78.8 s to 39.5 s direct, 28.9 s stored.
**Under a functional the potential does not respond to the field**, so a
pure functional has no coupled equations at all and `exchange_response`
returns zero without touching an integral; a hybrid scales exchange by
`xc%exx_fraction` in the first-order Fock matrix and the response alike.
What the functional adds is `giao_xc_potential`: the potential between a
function and the gauge derivative of another, `(1/2)(R x r) chi`, less its
transpose, on the grid, with that derivative's gradient for a GGA.
`pyscf.prop.nmr.rks`'s own `get_vxc_giao` fails under current PySCF (a block
size that is not a multiple of the grid's), so the reference is the same
formulas in one pass, held to the two numbers in that module's self-test --
and those were made with the old B3LYP, `b3lyp5` now, which is 0.02 ppm from
today's. **In a continuum the surface charges do not respond either**, to
the field or to a nucleus's contact field, so the couplings need nothing but
the solvated orbitals and the shieldings need one term: the gauge derivative
of the charges' potential, `giao_smeared_charges`. The charges are Gaussians,
and `int1e_ignuc` sums integer charges over every atom it is told of, so each
charge is put first in the atom list as a unit Gaussian nucleus and the
integral told there is one atom -- it still finds each shell's centre further
down the list. One pass over the shell pairs per surface point. No other
program here does a shielding in a continuum, so **the check is translation**:
6e-8 ppm with the term, 87 ppm without.
**In a fragmented run the shieldings are assembled like a gradient**:
`compute_mbe` notices that every term carries `nmr_shielding`, and
`accumulate_fragment_derivatives` folds each term's tensors onto the system
atoms with the expansion's collapsed weights, dropping caps. Nothing asks for
it -- the caller allocates nothing, unlike the gradient. The tensors travel
over MPI flattened to `(9, n_atoms)` after the dipole derivatives, in both
`send` routines and through `recv_shielding` in both receives; **the four have
to stay in step**. On four ranks the result is the serial one to 1e-11 ppm. The check that
needs no reference: at the level that reaches the whole system the result is
the unfragmented one, 1e-11 ppm on a water trimer. **Overlapping fragments
(GMBE) assemble them too, and are the way to cut through covalent bonds**:
`assemble_pie_shieldings` in `mqc_gmbe_fragment_distribution_scheme`, called
by the serial processor and the MPI coordinator alike, sums each atom over
the inclusion-exclusion terms that contain it. An intersection carries its
caps where its two fragments carry theirs, so a cap's disturbance enters
with both signs; measured on decane, six-carbon windows advancing two leave
carbon 0.41 ppm off, where disjoint two-carbon fragments at MBE(2) leave
18. **Reading each atom from the one fragment where it is furthest from a
cap was tried and is worse** -- 5.5 ppm on the same windows -- because a cap
two bonds away is still worth ppm to a carbon and nothing cancels it. A
window that overlaps its neighbour by one atom leaves that atom beside a cap
in every term. Couplings, counterpoise and embedding are not assembled,
and **a fragmented Hessian deck with `properties.nmr` is refused**: its terms
are displaced gradients, which lost the tensors and computed them 6N+1
times over. The validation manifest can pin shieldings
(`expected_nmr_isotropic`), because the energy does not notice whether they
were assembled. See
`mqc_docs/source/nmr_shielding.rst`.

## Spin-spin couplings

Spin-spin couplings are `properties.spin_spin_coupling`, a separate request
from the shieldings: a deck may want either without the other.
`mqc_czt_spin_coupling` computes **the Fermi-contact term only**, term for
term `pyscf.prop.ssc`'s `make_fc`: the right-hand side for nucleus `A` is
`(4 pi / 3) phi_a(R_A) phi_i(R_A)`, from `eval_ao_block` at the nuclei, and
the response is a *triplet* one, through `cphf_solve(..., triplet=.true.)`.
**That flag has to reach `xc_kernel_cache_fill` as well as the operator**: a
cache filled for singlets is refused by the triplet kernel, which is how it
was found. A triplet-unstable reference makes the operator indefinite and the
solve says so; at Hartree-Fock that is common and is the method's failure,
not a bug. What is stored is the reduced coupling in atomic units;
`mqc_nmr_isotopes` holds the observed isotope and g factor per element (13C,
not 12C) and turns it into hertz in the writer. **No basis set made for
couplings is in `basis_sets/`**, and the contact term is the property most
sensitive to that. See `mqc_docs/source/spin_spin_coupling.rst`.

## Raman activities

Raman activities are `properties.raman` on a Hessian deck, closed-shell HF or
KS with exact integrals. `mqc_czt_raman` has three routes to the polarizability's
Cartesian derivatives. The default, `polarizability_derivatives_analytic`, is
Hartree-Fock only and differences nothing: the gradient formula holds in a
field, so its second field derivative is the same contraction at `D^kl` and
`W^kl`, plus `Tr(D^k G^x[D^l])`. `D^kl` is six more `cphf_solve`s through
`mo_rhs`, with its diagonal blocks fixed by idempotency; `W^kl` is the product
rule on `(1/2) D F D`. **The two-electron terms pair two different
densities**, `D^kl` with `D` and `D^k` with `D^l`, and come from one
`two_electron_deriv_many` pass over a stack of all ten: each pairing is then
two traces over an atom's functions, *both orders*, because the sum over an
atom is on the first index alone (`add_pairing`). It used to be eighteen
`two_electron_gradient` calls, a product of two densities taken as the
difference of the gradient at their sum and at each -- ten times the time,
and a subtraction of nearly equal numbers. **The nine response solves share
one `(A+B)`**: `cphf_solve` iterates one right-hand side at a time with an
integral pass per iteration, about a hundred passes for nine, so where
`mo_transform_fits` allows the matrix is built once by `build_hessian_mo`,
factorized, and each solve is a back-substitution (`integral_direct` forces
the iteration; a test holds the two to 1e-9, where they land 2e-11 apart,
and asserts each call took the route it was meant to). The factorized route
refuses what the iterated one does: a reference with no gap, and one where
`b . (A+B)^-1 b` is not positive, which proves a saddle. The LU goes through
`factorize_response` in `mqc_czt_cphf`, which is where `MQC_SEQUENTIAL_BLAS`
is visible. The derivative pass over the ten-density stack holds two
potentials per density per thread, so it is taken in pieces against
`memory_budget`. Measured on four threads,
the Raman step of a Hessian run went from 24.8 s to 2.2 s on propane in
cc-pVDZ and from 55.2 s to 5.5 s on a water hexamer in 6-31G*. Its
check is the translational sum rule, 1e-13, and PySCF's differences
*extrapolated to zero step*; the plain ones are 3e-5 off, their own error. **A
Kohn-Sham reference takes the same route**, LDA, GGA and global hybrids
(`xc=` to `polarizability_derivatives_analytic`); a meta-GGA, a
range-separated hybrid and VV10 take the next one. Three additions, each a
derivative of the functional further: the kernel in every response solve,
`xc_kernel2_apply(D; D^k, D^l)` in the second-order Fock matrix -- there
whatever `D^kl` is, because the potential is not linear in the density -- and
in the gradient `xc_potential_gradient(D, D^kl)` plus `xc_kernel_gradient(D;
D^k, D^l)`, the derivative of `Tr(P_a f_xc P_b)` with all three densities
held and the grid moving. Both are threaded over grid blocks with a private
gradient per thread; **an internal procedure called from the parallel region
reaches the host's variable, not the thread's private copy**, so the block's
error is passed to `report_failure` as an argument. The exchange term is
`two_electron_deriv_many(..., exx_fraction=)`. **A context that was never
created is present, inactive, and says its exchange fraction is zero**: it is
what the bridge holds for Hartree-Fock, and passing it down made the mean
field one with no exchange, activities of 33.8 for 45.8 on water with every
test green. The routine forwards the context only when it is active, and
`has_functional` in `mqc_czt_response_product` is the test everywhere there
-- never `present(xc)` alone. No other program gives a
Kohn-Sham reference number here, so the check is the finite-field route,
which shares no formula with it, *extrapolated over two field strengths*:
1e-7 for an LDA and a GGA, 2e-7 for B3LYP at fields of 4e-3 and 2e-3.
**Unextrapolated it is 1e-4 off and hides a missing term; at smaller fields
the hybrid's reference is the noisy one** (7e-6 at 2e-3 and 1e-3), the
exchange build's noise divided by the field squared. The KS solves are
iterated, since the factorized matrix has no kernel, and **`cphf_solve` now
advances its right-hand sides in lockstep**, `LOCKSTEP_GROUP` at a time, one
`response_product` per iteration for every member still unconverged: each
member's iterates are the ones it would have had alone, and nine solves cost
two solves' passes. Every caller of `cphf_solve` gets that. On four
threads at B3LYP the Raman step is 22 s against the field route's 56 on
propane in cc-pVDZ, and 50 against 157 on a water hexamer in 6-31G*. `polarizability_derivatives_field`
differences the analytic gradient twice over twelve electric fields --
thirteen SCFs and gradients whatever the molecule's size, the field entering
through `run_czt_rhf`'s `h_extra` and its own gradient term added by hand,
since `czt_scf_gradient` does not know about it. The other,
`polarizability_derivatives`, differences `static_polarizability`
over `6N` displaced geometries, each a fresh molecule, SCF and response solve.
**The mixed second difference does not vanish where symmetry says the
derivative does**: its error carries `F_k^2 F_l^2` terms, which are even, so
a small nonzero `xy` component is truncation and not noise -- it falls by
four when the field is halved. **A fragmented Hessian run assembles them
too**: `compute_mbe` sets `compute_raman` when the dipole derivatives are
being assembled and every term carries `polarizability_derivatives`, and they
go through `redistribute_cap_dipole_derivatives`, which moves whole columns
and so takes six rows as it takes three. Over MPI they follow the dipole
derivatives in both sends and both receives. **The MPI coordinator ordinarily
flattens a fragmented Hessian into displaced gradients, and a queue of
gradients never computes a polarizability derivative**: the activities were
silently absent on more than one rank until `expand_displacements` learned to
stand down for `properties%raman`, leaving each fragment whole on one rank.
Found by running it under `mpirun`, not by reading it. The full expansion of a water
trimer reproduces the unfragmented activities to 8e-10; at the pair level a
*band's* total is within a per cent while single modes are off by half,
because the assembled Hessian mixes near-degenerate stretches differently --
compare bands. The bridge calls any route on the analytic-Hessian branch,
and the analytic and finite-field ones on a declined Hessian too
(`raman_derivatives`, contained in `run_czt_hf`): a capped fragment has no
analytic Hessian, and `hf_calc_hessian` carries the derivatives across the
`result%destroy()` that precedes the differenced one. **That capped path runs
and is unvalidated** -- caps change the fragments, so no expansion order
reproduces the whole molecule to check against;
`compute_raman_activities` in `mqc_vibrational_analysis` projects the result
on the modes beside the IR intensities. Two things cost time to find. **The
displaced SCFs need an orbital gradient of 1e-10**: at 1e-8 the derivatives
were 3e-6 off, the noise divided by the step. **A reference built from
another program's Hessian has to project out rotations too**, or a mode that
mixes with one differs by 0.1 cm-1 and 0.1 in its activity at a geometry
that is not a stationary point, which every sample geometry is not. The
reference is PySCF's Hessian with PySCF's own coupled-perturbed
polarizability differenced, through `cphf.solve`.

## Analyses on a converged wave function

Under `properties`, beside `keywords`: `bonding_analysis` (the quasi-atomic
analysis), `charges`, `bond_orders`, `fukui` and `natural_orbitals`. None
changes an energy, and a failure in one is a warning, not an error.

`properties.natural_orbitals` diagonalises the most complete one-particle
density the method has (`mqc_czt_natural_orbitals`, density in the AO basis,
overlap orthogonalised canonically): the spin-summed density of an
unrestricted reference, or the relaxed density of closed-shell `ri-mp2` --
the same `relaxed_density` PTD uses, fetched once when both ask. The name of
the density is stored beside the occupations because occupations of different
densities are not comparable. **A relaxed density's occupations are not bound
to [0, 2]**; the summary clamps them for the unpaired-electron counts and the
report says when one strays. Conventional, spin-scaled and unrestricted MP2,
coupled cluster, and MP2 under PCM `pte` are refused by name. CASSCF reports
its own active-space occupations and does not go through this.

`properties.intrinsic_orbitals` is Knizia's IAO/IBO (`mqc_czt_iao`): IAO
charges and the occupied orbitals localised by Jacobi rotations on
`sum_i sum_A (Q_A^i)^p`, `p` 4 or 2, HF or KS only. Its formulas
are `pyscf.lo.iao` and `pyscf.lo.ibo`'s and its tests are pinned to them.
**An unrestricted reference is two separate analyses**, one per spin
(`intrinsic_one_spin` in the bridge): the IAOs depend on the occupied space,
so alpha's and beta's differ, and there is no one set to take a spin density
into. Charges are the two populations added, spin populations their
difference.
**Its minimal basis is `basis_sets/minao/iao.json` (Knizia's MINAO, cc-pVTZ
contractions), not `minao/minao.json`** (ANO-RCC, the initial guess's): two
bases share the name upstream and give intrinsic charges 0.03 apart on water.
**Two orbitals on one atom are indistinguishable to the functional**, so the
sweep skips a pair whose `A^2 + B^2` vanishes and converges on the gradient;
converging on the rotation angle never terminates.

Charge schemes are `mulliken`, `lowdin`, `chelpg`, `hirshfeld`, `mbis`; bond
order schemes are `mayer` and `wiberg_lowdin`. **Lowdin and Wiberg are Mulliken and
Mayer in the symmetrically orthogonalised basis**: `lowdin_density` in
`mqc_population_analysis` returns `S^1/2 D S^1/2`, and the existing routines
are called with it and `identity_overlap`. **The bare `wiberg` is refused by
the schema, by name**: an NBO analysis means the natural-atomic-orbital index
by it, which is a different number and is not implemented, so the scheme
here carries its basis in its name and the bare word is left for that one.

`properties.density: "relaxed"` points `charges` and `bond_orders` at that
same relaxed density; the default, `"scf"`, is the reference's whatever runs
on top, and `result%analysis_density` (the `density` field of both output
sections) says which. **The reference's dipole is withheld after a correlated
method on purpose**; the relaxed density's is published whenever that density
is built, through `publish_dipole` in the bridge.
