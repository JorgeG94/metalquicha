!! The effective-density contraction against the SCF gradient
module test_mqc_density_contraction
   !! An SCF gradient is the one derivative whose effective densities are
   !! known in closed form: P is the density, W the energy-weighted density,
   !! and the two-particle part is all separable pairs of the density with
   !! itself. Handing those to `contract_effective_density` through a fixed
   !! provider and adding the nuclear repulsion must give what
   !! `czt_scf_gradient` gives, which builds the same derivative by its own
   !! single-pass route. That fixes the factors of every term the contraction
   !! forms: the Pulay factor two, the core-Hamiltonian trace, and the 1/2 and
   !! 1/4 that make a pair `D . (J^x - K^x/2)(D)/2`.
   !!
   !! Restricted checks Coulomb and exchange in their fixed RHF ratio.
   !! Unrestricted takes them apart: Coulomb on the total density alone,
   !! exchange on each spin alone, so a pair weight read the wrong way round,
   !! or one term dropped, cannot cancel. Splitting the density in two across
   !! the pool checks a cross pair counts both of its orderings.
   !!
   !! The two routes screen differently -- `czt_scf_gradient` at 1e-12 on the
   !! gradient, the many-density sweep at 1e-10 on the potential -- so in
   !! general they agree to the screening rather than bit for bit. On water in
   !! 6-31G nothing that matters is screened and they agree to 5e-15; a 10 per
   !! cent error in the exchange weight lands at 2e-2.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf
   use mqc_czt_gradient, only: czt_scf_gradient, nuclear_repulsion_gradient
   use mqc_czt_effective_density, only: effective_density_t, separable_pair_t, column_t, &
                                        COLUMN_KIND_GRADIENT
   use mqc_czt_density_contraction, only: contract_effective_density
   implicit none
   private

   public :: collect_mqc_density_contraction_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.190459_dp, &
                           0.0_dp, 1.459853_dp, -0.884001_dp, &
                           0.0_dp, -1.459853_dp, -0.884001_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: AGREEMENT = 1.0e-10_dp
      !! Hartree/Bohr between the two routes; see the module header

   type, extends(effective_density_t) :: fixed_density_t
      !! One gradient column, every matrix given
      real(dp), allocatable :: p(:, :), w(:, :), pool(:, :, :)
      type(separable_pair_t), allocatable :: pairs(:)
   contains
      procedure :: n_columns => fixed_n_columns
      procedure :: column => fixed_column
      procedure :: energy_weighted => fixed_w
      procedure :: one_particle => fixed_p
      procedure :: density_pool => fixed_pool
      procedure :: separable_pairs => fixed_pairs
   end type fixed_density_t

contains

   subroutine collect_mqc_density_contraction_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("restricted_scf_through_the_contraction", test_restricted), &
                  new_unittest("a_density_split_across_the_pool", test_split_pool), &
                  new_unittest("unrestricted_coulomb_and_exchange_apart", test_unrestricted), &
                  new_unittest("a_pair_outside_the_pool_is_refused", test_bad_pair) &
                  ]
   end subroutine collect_mqc_density_contraction_tests

   ! ---- the fixed provider -------------------------------------------------

   function fixed_n_columns(self) result(n)
      class(fixed_density_t), intent(in) :: self
      integer :: n
      associate (unused_self => self)
      end associate
      n = 1
   end function fixed_n_columns

   function fixed_column(self, ic) result(col)
      class(fixed_density_t), intent(in) :: self
      integer, intent(in) :: ic
      type(column_t) :: col
      associate (unused_self => self, unused_ic => ic)
      end associate
      col = column_t(kind=COLUMN_KIND_GRADIENT, bra=1, ket=1)
   end function fixed_column

   subroutine fixed_w(self, ic, matrix, error)
      class(fixed_density_t), intent(in) :: self
      integer, intent(in) :: ic
      real(dp), allocatable, intent(out) :: matrix(:, :)
      type(error_t), intent(inout) :: error
      associate (unused_ic => ic, unused_error => error)
      end associate
      matrix = self%w
   end subroutine fixed_w

   subroutine fixed_p(self, ic, matrix, error)
      class(fixed_density_t), intent(in) :: self
      integer, intent(in) :: ic
      real(dp), allocatable, intent(out) :: matrix(:, :)
      type(error_t), intent(inout) :: error
      associate (unused_ic => ic, unused_error => error)
      end associate
      matrix = self%p
   end subroutine fixed_p

   subroutine fixed_pool(self, pool, error)
      class(fixed_density_t), intent(in) :: self
      real(dp), allocatable, intent(out) :: pool(:, :, :)
      type(error_t), intent(inout) :: error
      associate (unused_error => error)
      end associate
      pool = self%pool
   end subroutine fixed_pool

   subroutine fixed_pairs(self, ic, pairs, error)
      class(fixed_density_t), intent(in) :: self
      integer, intent(in) :: ic
      type(separable_pair_t), allocatable, intent(out) :: pairs(:)
      type(error_t), intent(inout) :: error
      associate (unused_ic => ic, unused_error => error)
      end associate
      pairs = self%pairs
   end subroutine fixed_pairs

   ! ---- helpers ------------------------------------------------------------

   function weighted_density(orbitals, energies, n_occ, occupation) result(w)
      !! `sum_i occupation eps_i C_i C_i^T`
      real(dp), intent(in) :: orbitals(:, :), energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: occupation
      real(dp), allocatable :: w(:, :)
      integer :: i, mu
      allocate (w(size(orbitals, 1), size(orbitals, 1)))
      w = 0.0_dp
      do i = 1, n_occ
         do mu = 1, size(orbitals, 1)
            w(:, mu) = w(:, mu) + occupation*energies(i)*orbitals(:, i)*orbitals(mu, i)
         end do
      end do
   end function weighted_density

   subroutine through_contraction(mol, provider, gradient, err)
      !! The provider's one column, plus the nuclear repulsion
      type(czt_molecule_t), intent(in) :: mol
      type(fixed_density_t), intent(in) :: provider
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: err
      real(dp), allocatable :: columns(:, :, :)
      call contract_effective_density(mol, provider, columns, err)
      if (err%has_error()) return
      allocate (gradient(3, mol%natm))
      gradient = 0.0_dp
      call nuclear_repulsion_gradient(mol, gradient)
      gradient = gradient + columns(:, :, 1)
   end subroutine through_contraction

   subroutine restricted_water(mol, scf, error)
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "6-31g", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return
      call run_czt_rhf(mol, 10, 300, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error, scf%converged, "the SCF should converge")
   end subroutine restricted_water

   ! ---- tests --------------------------------------------------------------

   subroutine test_restricted(error)
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(fixed_density_t) :: provider
      type(error_t) :: err
      real(dp), allocatable :: g_ref(:, :), g(:, :)

      call restricted_water(mol, scf, error)
      if (allocated(error)) return
      call czt_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                            orbital_energies=scf%orbital_energies, &
                            n_occupied=scf%n_occupied, gradient=g_ref, error=err)
      call check(error,.not. err%has_error(), "the SCF gradient should build")
      if (allocated(error)) return

      provider%p = scf%density
      provider%w = weighted_density(scf%orbitals, scf%orbital_energies, scf%n_occupied, 2.0_dp)
      allocate (provider%pool(mol%nao, mol%nao, 1))
      provider%pool(:, :, 1) = scf%density
      provider%pairs = [separable_pair_t(left=1, right=1, coulomb=0.5_dp, exchange=0.25_dp)]

      call through_contraction(mol, provider, g, err)
      call check(error,.not. err%has_error(), "the contraction should run")
      if (allocated(error)) return
      write (*, "(a, es10.2)") "    max |contracted - czt_scf_gradient| =", maxval(abs(g - g_ref))
      call check(error, maxval(abs(g - g_ref)) < AGREEMENT, &
                 "the contracted RHF gradient matches czt_scf_gradient")
   end subroutine test_restricted

   subroutine test_split_pool(error)
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(fixed_density_t) :: provider
      type(error_t) :: err
      real(dp), allocatable :: g_ref(:, :), g(:, :)
      real(dp) :: share

      call restricted_water(mol, scf, error)
      if (allocated(error)) return
      call czt_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                            orbital_energies=scf%orbital_energies, &
                            n_occupied=scf%n_occupied, gradient=g_ref, error=err)
      if (err%has_error()) return

      ! D = D_1 + D_2 with unequal shares, so (D_1 + D_2)(J - K/2)(D_1 + D_2)/2
      ! needs the cross pair at full weight and both of its orderings.
      share = 0.3_dp
      provider%p = scf%density
      provider%w = weighted_density(scf%orbitals, scf%orbital_energies, scf%n_occupied, 2.0_dp)
      allocate (provider%pool(mol%nao, mol%nao, 2))
      provider%pool(:, :, 1) = share*scf%density
      provider%pool(:, :, 2) = (1.0_dp - share)*scf%density
      provider%pairs = [separable_pair_t(left=1, right=1, coulomb=0.5_dp, exchange=0.25_dp), &
                        separable_pair_t(left=1, right=2, coulomb=1.0_dp, exchange=0.5_dp), &
                        separable_pair_t(left=2, right=2, coulomb=0.5_dp, exchange=0.25_dp)]

      call through_contraction(mol, provider, g, err)
      call check(error,.not. err%has_error(), "the contraction should run")
      if (allocated(error)) return
      write (*, "(a, es10.2)") "    max |contracted - czt_scf_gradient| =", maxval(abs(g - g_ref))
      call check(error, maxval(abs(g - g_ref)) < AGREEMENT, &
                 "a density split across the pool gives the same gradient")
   end subroutine test_split_pool

   subroutine test_unrestricted(error)
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(fixed_density_t) :: provider
      type(error_t) :: err
      real(dp), allocatable :: g_ref(:, :), g(:, :)
      integer :: nao

      ! The water cation: a doublet, so the two spin densities differ.
      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "6-31g", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return
      call run_czt_uhf(mol, 9, 2, 300, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                       grad_tol=1.0e-9_dp)
      call check(error, scf%converged, "the UHF should converge")
      if (allocated(error)) return

      call czt_scf_gradient(mol, scf%density, scf%density_beta, scf%orbitals, &
                            scf%orbitals_beta, scf%orbital_energies, &
                            scf%orbital_energies_beta, scf%n_occupied, &
                            scf%n_occupied_beta, g_ref, err)
      call check(error,.not. err%has_error(), "the UHF gradient should build")
      if (allocated(error)) return

      nao = mol%nao
      provider%p = scf%density + scf%density_beta
      provider%w = weighted_density(scf%orbitals, scf%orbital_energies, scf%n_occupied, 1.0_dp) &
                   + weighted_density(scf%orbitals_beta, scf%orbital_energies_beta, &
                                      scf%n_occupied_beta, 1.0_dp)
      allocate (provider%pool(nao, nao, 3))
      provider%pool(:, :, 1) = scf%density + scf%density_beta
      provider%pool(:, :, 2) = scf%density
      provider%pool(:, :, 3) = scf%density_beta
      ! E2 = D J(D)/2 - (D_a K(D_a) + D_b K(D_b))/2: Coulomb and exchange never
      ! meet in one pair.
      provider%pairs = [separable_pair_t(left=1, right=1, coulomb=0.5_dp, exchange=0.0_dp), &
                        separable_pair_t(left=2, right=2, coulomb=0.0_dp, exchange=0.5_dp), &
                        separable_pair_t(left=3, right=3, coulomb=0.0_dp, exchange=0.5_dp)]

      call through_contraction(mol, provider, g, err)
      call check(error,.not. err%has_error(), "the contraction should run")
      if (allocated(error)) return
      write (*, "(a, es10.2)") "    max |contracted - czt_scf_gradient| =", maxval(abs(g - g_ref))
      call check(error, maxval(abs(g - g_ref)) < AGREEMENT, &
                 "the contracted UHF gradient matches czt_scf_gradient")
   end subroutine test_unrestricted

   subroutine test_bad_pair(error)
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(fixed_density_t) :: provider
      type(error_t) :: err
      real(dp), allocatable :: columns(:, :, :)

      call restricted_water(mol, scf, error)
      if (allocated(error)) return
      provider%p = scf%density
      provider%w = scf%density
      allocate (provider%pool(mol%nao, mol%nao, 1))
      provider%pool(:, :, 1) = scf%density
      provider%pairs = [separable_pair_t(left=1, right=2, coulomb=1.0_dp, exchange=0.5_dp)]

      call contract_effective_density(mol, provider, columns, err)
      call check(error, err%has_error(), "a pair indexing past the pool is an error")
      if (allocated(error)) return
      call check(error, err%get_code() == ERROR_VALIDATION, "and a validation error")
   end subroutine test_bad_pair

end module test_mqc_density_contraction

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_density_contraction, only: collect_mqc_density_contraction_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_density_contraction", collect_mqc_density_contraction_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
