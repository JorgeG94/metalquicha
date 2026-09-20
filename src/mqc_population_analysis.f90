!! Population analysis and bond orders, over plain arrays
module mqc_population_analysis
   !! Everything here needs the same three things: a density, an overlap and a
   !! statement of which atom each basis function belongs to. None of that is
   !! particular to an integrals backend: this module holds the arithmetic and
   !! each backend supplies the three arrays its own way.
   !!
   !! Backends differ only in how the AO-to-atom map is arrived at -- libcint
   !! reads it off its shell table, cuEST counts functions per atom -- and both
   !! end up as the same `owner` array.
   !!
   !! Two quantities live here, and they answer different questions from the
   !! same matrices. The Mulliken partition asks how much of the density sits
   !! *on* an atom; Mayer's bond order asks how much of it is shared *between*
   !! two. Neither is the other's diagonal.
   !!
   !! **What is not here.** The Wiberg-Mayer orders `mqc_method_xtb` reports
   !! come out of a semi-empirical Hamiltonian in its own minimal basis, and
   !! the kinetic bond orders in `mqc_czt_quao` are Ruedenberg's definition in
   !! the quasi-atomic basis. All three are called bond orders, none of them
   !! agree to more than a trend, and each is only interpretable alongside the
   !! name of the scheme that produced it.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: ao_owner_from_counts
   public :: mulliken_atomic_charges
   public :: mulliken_atomic_spin_populations
   public :: mayer_atomic_bond_orders
   public :: mayer_atomic_bond_orders_open_shell
   public :: mayer_atomic_valences

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

   subroutine mayer_atomic_bond_orders(owner, natm, density, overlap, orders, error)
      !! B_AB = sum_{mu in A} sum_{nu in B} (D S)_mu,nu (D S)_nu,mu
      !!
      !! Mayer's bond order over a *total* closed-shell density: the part of
      !! the charge-density matrix `D S` that spans the two atoms, squared. It
      !! counts electron pairs shared between A and B, so a covalent single
      !! bond lands near one and two atoms that share nothing land near zero,
      !! with no distance ever entering.
      !!
      !! **Only valid for a density with D_alpha = D_beta.** The open-shell
      !! expression is not this one with the total density -- see
      !! `mayer_atomic_bond_orders_open_shell`, which reduces to this when the
      !! two spin densities are equal. Handing an unrestricted total density to
      !! this routine gives numbers that look reasonable and are wrong, which
      !! is why the two entries are separate rather than one with an optional
      !! argument.
      !!
      !! The diagonal is zeroed: `B_AA` is a valence-like quantity of a
      !! different kind, not a bond of an atom with itself, and leaving it in
      !! makes every sum over a row wrong.
      integer, intent(in) :: owner(:)               !! (nao) atom of each function
      integer, intent(in) :: natm                   !! Number of atoms
      real(dp), intent(in) :: density(:, :), overlap(:, :)
      real(dp), allocatable, intent(out) :: orders(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ds(:, :)

      call charge_density_matrix(owner, natm, density, overlap, "mayer bond orders", &
                                 ds, error)
      if (error%has_error()) return
      call accumulate_pair_blocks(owner, natm, ds, 1.0_dp, orders)
   end subroutine mayer_atomic_bond_orders

   subroutine mayer_atomic_bond_orders_open_shell(owner, natm, density_alpha, &
                                                  density_beta, overlap, orders, error)
      !! B_AB = 2 sum_{mu in A, nu in B} [ (Da S)_mu,nu (Da S)_nu,mu
      !!                                 + (Db S)_mu,nu (Db S)_nu,mu ]
      !!
      !! The open-shell Mayer order, per spin and then doubled. With
      !! `Da = Db = D/2` the two terms are each a quarter of the closed-shell
      !! sum and the factor of two restores it, so this *is* the closed-shell
      !! formula for a closed shell -- and differs from it everywhere else. A
      !! doublet radical is where the difference shows: the singly occupied
      !! orbital contributes to one spin only.
      integer, intent(in) :: owner(:)               !! (nao) atom of each function
      integer, intent(in) :: natm                   !! Number of atoms
      real(dp), intent(in) :: density_alpha(:, :), density_beta(:, :)
      real(dp), intent(in) :: overlap(:, :)
      real(dp), allocatable, intent(out) :: orders(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ds(:, :)

      call charge_density_matrix(owner, natm, density_alpha, overlap, &
                                 "mayer bond orders", ds, error)
      if (error%has_error()) return
      call accumulate_pair_blocks(owner, natm, ds, 2.0_dp, orders)

      deallocate (ds)
      call charge_density_matrix(owner, natm, density_beta, overlap, &
                                 "mayer bond orders", ds, error)
      if (error%has_error()) return
      call accumulate_pair_blocks(owner, natm, ds, 2.0_dp, orders, add=.true.)
   end subroutine mayer_atomic_bond_orders_open_shell

   subroutine mayer_atomic_valences(orders, valences)
      !! V_A = sum_{B /= A} B_AB, the Mayer valence
      !!
      !! Free once the matrix exists, and the number a reader looks at first:
      !! a carbon in a hydrocarbon comes out near four, a nitrogen near three,
      !! and a large shortfall says the atom has a lone electron or the
      !! reference is not describing it.
      !!
      !! The diagonal is not summed because the matrix carries none -- see
      !! `mayer_atomic_bond_orders`.
      real(dp), intent(in) :: orders(:, :)
      real(dp), allocatable, intent(out) :: valences(:)

      integer :: iatom

      allocate (valences(size(orders, 1)))
      do iatom = 1, size(orders, 1)
         valences(iatom) = sum(orders(iatom, :))
      end do
   end subroutine mayer_atomic_valences

   subroutine charge_density_matrix(owner, natm, density, overlap, what, ds, error)
      !! `D S`, with the shapes and the owner map checked first
      !!
      !! Through `pic_gemm` rather than a triple loop: this is the only O(N^3)
      !! step in a bond-order analysis and everything after it is O(N^2).
      integer, intent(in) :: owner(:)
      integer, intent(in) :: natm
      real(dp), intent(in) :: density(:, :), overlap(:, :)
      character(len=*), intent(in) :: what
      real(dp), allocatable, intent(out) :: ds(:, :)
      type(error_t), intent(inout) :: error

      integer :: nao

      nao = size(owner)
      if (size(density, 1) /= nao .or. size(density, 2) /= nao .or. &
          size(overlap, 1) /= nao .or. size(overlap, 2) /= nao) then
         call error%set(ERROR_VALIDATION, what//": density and overlap must be "// &
                        "square and the size of the basis")
         return
      end if
      ! Checked rather than assumed: the accumulation below indexes the result
      ! by `owner`, so a map naming an atom the matrix does not have would
      ! write past it. Cheap next to the GEMM that follows.
      if (nao > 0) then
         if (minval(owner) < 1 .or. maxval(owner) > natm) then
            call error%set(ERROR_VALIDATION, what//": a basis function is owned by an "// &
                           "atom outside the molecule")
            return
         end if
      end if

      allocate (ds(nao, nao))
      call pic_gemm(density, overlap, ds)
   end subroutine charge_density_matrix

   subroutine accumulate_pair_blocks(owner, natm, ds, factor, orders, add)
      !! sum_{mu in A, nu in B} factor * (D S)_mu,nu (D S)_nu,mu, per atom pair
      !!
      !! One pass over the basis-function pairs, each one landing in the block
      !! its two owners name. O(N_ao^2) and no atom-by-atom search.
      integer, intent(in) :: owner(:)
      integer, intent(in) :: natm
      real(dp), intent(in) :: ds(:, :)
      real(dp), intent(in) :: factor
      real(dp), allocatable, intent(inout) :: orders(:, :)
      logical, intent(in), optional :: add   !! Add into `orders` rather than start it

      integer :: mu, nu, iatom, jatom
      logical :: accumulate

      accumulate = .false.
      if (present(add)) accumulate = add
      if (.not. accumulate) then
         if (allocated(orders)) deallocate (orders)
         allocate (orders(natm, natm), source=0.0_dp)
      end if

      do nu = 1, size(owner)
         jatom = owner(nu)
         do mu = 1, size(owner)
            iatom = owner(mu)
            orders(iatom, jatom) = orders(iatom, jatom) + &
                                   factor*ds(mu, nu)*ds(nu, mu)
         end do
      end do

      ! An atom is not bonded to itself. What sits there is a different
      ! quantity -- roughly twice the atom's own population -- and summing a
      ! row with it left in would report it as bonding.
      do iatom = 1, natm
         orders(iatom, iatom) = 0.0_dp
      end do
   end subroutine accumulate_pair_blocks

end module mqc_population_analysis
