!! Mulliken population analysis, over plain arrays
module mqc_population_analysis
   !! The Mulliken partition needs a density, an overlap and a statement of
   !! which atom each basis function belongs to. None of that is particular to
   !! an integrals backend: this module holds the arithmetic and each backend
   !! supplies the three arrays its own way.
   !!
   !! Backends differ only in how the AO-to-atom map is arrived at -- libcint
   !! reads it off its shell table, cuEST counts functions per atom -- and both
   !! end up as the same `owner` array.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: ao_owner_from_counts
   public :: mulliken_atomic_charges
   public :: mulliken_atomic_spin_populations

contains

   subroutine ao_owner_from_counts(counts, owner)
      !! Expand a per-atom count of basis functions into a per-function owner
      !!
      !! Valid only where the AO index runs over atoms in order, with each
      !! atom's functions contiguous. That is a property of how the basis was
      !! built, so the caller is asserting it rather than this checking it.
      integer, intent(in) :: counts(:)              !! (natm) functions on each atom
      integer, allocatable, intent(out) :: owner(:)  !! (nao) atom index, 1-based

      integer :: iatom, k, mu

      allocate (owner(sum(counts)))
      mu = 0
      do iatom = 1, size(counts)
         do k = 1, counts(iatom)
            mu = mu + 1
            owner(mu) = iatom
         end do
      end do
   end subroutine ao_owner_from_counts

   subroutine mulliken_atomic_charges(owner, nuclear_charges, density, overlap, &
                                      charges, error)
      !! q_A = Z_A - sum_{mu in A} (D S)_mu,mu
      !!
      !! The diagonal of `D S` is the gross population of each basis function,
      !! and summing it over an atom's functions charges that atom with every
      !! overlap it takes part in, half of which belongs to its neighbour.
      !!
      !! `density` is the *total* density: an unrestricted caller adds its two
      !! spin blocks before calling, and a closed-shell one already has it,
      !! since that build carries the factor of two.
      integer, intent(in) :: owner(:)               !! (nao) atom of each function
      real(dp), intent(in) :: nuclear_charges(:)    !! (natm), the Z the SCF saw
      real(dp), intent(in) :: density(:, :), overlap(:, :)
      real(dp), allocatable, intent(out) :: charges(:)
      type(error_t), intent(inout) :: error

      real(dp) :: population
      integer :: mu, nu, nao

      nao = size(owner)
      if (size(density, 1) /= nao .or. size(overlap, 1) /= nao) then
         call error%set(ERROR_VALIDATION, "mulliken charges: density and overlap must "// &
                        "be the size of the basis")
         return
      end if

      allocate (charges(size(nuclear_charges)))
      charges = nuclear_charges     ! nuclear charge to start from

      do mu = 1, nao
         population = 0.0_dp
         do nu = 1, nao
            population = population + density(mu, nu)*overlap(nu, mu)
         end do
         charges(owner(mu)) = charges(owner(mu)) - population
      end do
   end subroutine mulliken_atomic_charges

   subroutine mulliken_atomic_spin_populations(owner, natm, spin_density, overlap, &
                                               populations, error)
      !! The Mulliken partition applied to P_alpha - P_beta
      !!
      !! Same trace, different matrix: where `mulliken_atomic_charges` asks how
      !! many electrons sit on an atom, this asks how many *unpaired* ones do.
      !! The nuclear charge does not enter and nothing is subtracted from it --
      !! a nucleus carries no spin -- so this is a population rather than a
      !! charge, and it sums to `n_alpha - n_beta` and not to the molecular
      !! charge.
      !!
      !! `natm` is passed rather than inferred: `maxval(owner)` would silently
      !! shrink the result for a molecule whose last atom carries no basis
      !! functions.
      integer, intent(in) :: owner(:)               !! (nao) atom of each function
      integer, intent(in) :: natm                   !! Number of atoms
      real(dp), intent(in) :: spin_density(:, :), overlap(:, :)
      real(dp), allocatable, intent(out) :: populations(:)
      type(error_t), intent(inout) :: error

      real(dp) :: population
      integer :: mu, nu, nao

      nao = size(owner)
      if (size(spin_density, 1) /= nao .or. size(overlap, 1) /= nao) then
         call error%set(ERROR_VALIDATION, "mulliken spin populations: density and "// &
                        "overlap must be the size of the basis")
         return
      end if

      allocate (populations(natm), source=0.0_dp)

      do mu = 1, nao
         population = 0.0_dp
         do nu = 1, nao
            population = population + spin_density(mu, nu)*overlap(nu, mu)
         end do
         populations(owner(mu)) = populations(owner(mu)) + population
      end do
   end subroutine mulliken_atomic_spin_populations

end module mqc_population_analysis
