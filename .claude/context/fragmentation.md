# Fragmentation, the MBE and how a run is distributed

What an expansion is, how a covalent bond is cut, the types a fragment and a
result are held in, and how the work goes over MPI. Read from `AGENTS.md`, which says when.

## Many-body expansion (MBE)
- Breaks system into fragments (monomers, dimers, trimers...)
- E_MBE = ΣE_I - ΣΔE_IJ + ΣΔE_IJK - ...
- Truncate at level N for O(N^k) scaling

## Generalized MBE (GMBE)
- Handles overlapping fragments
- Uses Principle of Inclusion-Exclusion (PIE)
- Automatically computes intersection terms

## Cutting a covalent bond

`keywords.fragmentation.bond_breaking` is `"none"` by default and refuses a
partition that severs one. `"afo"` detaches it with an adjusted frozen orbital
instead: a model system around the bond is solved and localized, the orbital on
the bond is reduced to the detached atom's own basis functions, and that hybrid
is frozen -- empty in the fragment that gets nothing of the bond, occupied in
the one that gets all of it. Runs with `embedding = "none"` and `"ptc"`, and is
refused with `"exact"`.

**Under a field the detached atom's nucleus is split along with the bond
pair** -- `Z-1` with the fragment that owns the atom, `+1` on the ghost centre
its neighbour already carries there -- so both fragments are neutral closed
shells. That is what a fragment *potential* needs: a net charge per fragment is
a monopole of order 1/R on every adjacent-residue pair. It leaves an embedded
FMO total *exactly* unchanged, because the field a group feels is
`q_all - own_q` and the charge reappears in it at the same point; that
invariance is asserted by a test rather than believed. With `embedding =
"none"` the nucleus stays whole, because nothing would supply the other half.
See `mqc_docs/source/fmo.rst`.

The bookkeeping rule that matters: **a cut belongs to a group, not a fragment.**
A bond severed between two monomers is whole again inside the dimer holding both
ends, so that dimer carries no ghost, no frozen orbital and no electron shift
there. Deciding it per fragment and inheriting it is what made an earlier capped
version 11 Hartree wrong. See `mqc_czt_afo.f90` and `mqc_docs/source/fmo.rst`.

## Hydrogen capping
- Broken covalent bonds get H-cap atoms
- Forces/Hessian elements redistributed to original heavy atoms
- Implemented in `mqc_physical_fragment.f90`

## Distance screening
- Filter fragments by inter-monomer distance
- Per-level cutoffs (dimer, trimer, etc.)
- 5-10x speedup for large systems

## Key Types

```fortran
! Main result container (src/core/mqc_result_types.f90)
type :: calculation_result_t
  type(energy_t) :: energy           ! SCF, MP2, CC energies
  real(dp), allocatable :: gradient(:,:)
  real(dp), allocatable :: hessian(:,:)
  real(dp) :: dipole(3)
  logical :: has_energy, has_gradient, has_hessian, has_error
end type

! Fragment representation (src/fragmentation/mqc_physical_fragment.f90)
type :: physical_fragment_t
  integer :: n_atoms
  integer, allocatable :: element_numbers(:)    ! Atomic numbers (Z)
  real(dp), allocatable :: coordinates(:,:)     ! (3, n_atoms) in Bohr
  integer :: charge, multiplicity, nelec
  integer :: n_caps                             ! Number of H-caps
  integer, allocatable :: cap_replaces_atom(:)  ! Original atom each cap replaces
  integer, allocatable :: local_to_global(:)    ! Map to system indices
end type

! System geometry (src/fragmentation/mqc_physical_fragment.f90)
type :: system_geometry_t
  integer :: n_monomers, total_atoms
  integer, allocatable :: element_numbers(:)
  real(dp), allocatable :: coordinates(:,:)     ! (3, total_atoms) in Bohr
  integer :: charge, multiplicity
end type

! QC method base class (src/methods/dispatch/mqc_method_base.f90)
type, abstract :: qc_method_t
contains
  procedure(calc_energy_interface), deferred :: calc_energy
  procedure(calc_gradient_interface), deferred :: calc_gradient
  procedure(calc_hessian_interface), deferred :: calc_hessian
end type
```

## Calculation Types

- `CALC_TYPE_ENERGY` - Single-point energy
- `CALC_TYPE_GRADIENT` - First derivatives
- `CALC_TYPE_HESSIAN` - Second derivatives + vibrational analysis
- `CALC_TYPE_INTERACTION_ENERGY` - `driver: "InteractionEnergy"`: the MBE
  terms containing `keywords.fragmentation.reference_fragment` (0-based), from
  a term list reduced to their closure under subsets. Fragments run as
  energies; the output has an `interaction_energy` block and **no**
  `total_energy`. Under `method: "fmo"` at level 2 (FMO2, PIEDA allowed) every
  monomer is still solved and only the reference's pairs follow; FMO3, EE-MBE
  and EFMO are refused. See `mqc_docs/source/interaction_energy.rst`

## MPI Architecture

- **Global coordinator**: Distributes fragments across nodes
- **Node coordinator**: Manages workers within a node
- **Workers**: Compute individual fragments
- Tags defined in `src/parallel/mqc_mpi_tags.f90`
- **Hessian parallelization**: Finite difference displacements are distributed across MPI ranks, even for unfragmented systems

## Output Format

JSON output (`output_<basename>.json`) contains:
- Total energy and per-fragment breakdown
- Gradients (if requested)
- Hessian (if requested)
- Vibrational frequencies and thermochemistry
- GMBE coefficients and intersection terms
