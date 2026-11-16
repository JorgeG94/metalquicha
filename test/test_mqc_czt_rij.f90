!! The fitted Coulomb matrix without the fitted tensor
module test_mqc_czt_rij
   !! `rij_coulomb` builds `J` from the raw three-centre integrals and a solve
   !! against the metric. The fitted tensor `B = (mu nu|P) V^-1/2` gives the
   !! same matrix as `B (B^T D)`, by a different route through different
   !! storage, so the two agreeing to rounding is the whole test: both are the
   !! same fit, and there is no approximation between them.
   !!
   !! Both ways of holding the integrals are checked against that reference
   !! and against each other: stored, and recomputed in every build.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, build_df_tensor
   use mqc_czt_rij, only: rij_t, build_rij, rij_coulomb
   implicit none
   private

   public :: collect_mqc_czt_rij_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

contains

   subroutine collect_mqc_czt_rij_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rij_matches_the_fitted_tensor_spherical", test_spherical), &
                  new_unittest("rij_matches_the_fitted_tensor_cartesian", test_cartesian), &
                  new_unittest("rij_keeps_only_the_symmetric_part", test_symmetric_part), &
                  new_unittest("rij_refuses_what_it_cannot_build", test_refusals) &
                  ]
   end subroutine collect_mqc_czt_rij_tests

   subroutine model_density(n, density)
      !! A symmetric matrix with elements of every size and sign
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: density(:, :)

      integer :: mu, nu

      allocate (density(n, n))
      do nu = 1, n
         do mu = 1, n
            density(mu, nu) = cos(0.7_dp*real(mu + nu, dp))/(1.0_dp + abs(real(mu - nu, dp))) &
                              + merge(1.0_dp, 0.0_dp, mu == nu)
         end do
      end do
   end subroutine model_density

   subroutine reference_coulomb(mol, aux, density, j, error)
      !! `J = B (B^T D)` out of the fitted tensor
      type(czt_molecule_t), intent(in) :: mol, aux
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: j(:, :)
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: b(:, :), gamma(:)
      type(error_t) :: err

      call build_df_tensor(mol, aux, b, err)
      call check(error,.not. err%has_error(), "the fitted tensor failed: "//err%get_full_trace())
      if (allocated(error)) return
      gamma = matmul(transpose(b), reshape(density, [mol%nao*mol%nao]))
      j = reshape(matmul(b, gamma), [mol%nao, mol%nao])
   end subroutine reference_coulomb

   subroutine compare(basis, aux_basis, expect_cartesian, error)
      !! Stored and recomputed integrals, both against the fitted tensor
      character(len=*), intent(in) :: basis, aux_basis
      logical, intent(in) :: expect_cartesian
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rij_t) :: stored, direct
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), j_ref(:, :), j_stored(:, :), j_direct(:, :)

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, basis, mol, err)
      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, aux_basis, aux, err)
      call check(error,.not. err%has_error(), "the molecules should build: "//err%get_full_trace())
      if (allocated(error)) return
      ! Guards the guard: the case is named for an angular form.
      call check(error, mol%cartesian .eqv. expect_cartesian, &
                 "the basis is not in the angular form this case is for")
      if (allocated(error)) return

      call model_density(mol%nao, density)
      call reference_coulomb(mol, aux, density, j_ref, error)
      if (allocated(error)) return

      call build_rij(mol, aux, stored, err, in_core=.true.)
      call build_rij(mol, aux, direct, err, in_core=.false.)
      call check(error,.not. err%has_error(), "the Coulomb build should set up: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, allocated(stored%three) .and. .not. allocated(direct%three), &
                 "the integrals are not held the way each build was asked to hold them")
      if (allocated(error)) return

      call rij_coulomb(mol, stored, density, j_stored, err)
      call rij_coulomb(mol, direct, density, j_direct, err)
      call check(error,.not. err%has_error(), "the Coulomb build failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, maxval(abs(j_stored - j_ref)) < 1.0e-9_dp, &
                 "the Coulomb matrix from stored integrals is not the fitted tensor's")
      if (allocated(error)) return
      call check(error, maxval(abs(j_direct - j_stored)) < 1.0e-11_dp, &
                 "recomputing the integrals gives a different Coulomb matrix from storing them")
      if (allocated(error)) return
      call check(error, maxval(abs(j_ref)) > 1.0_dp, "the reference Coulomb matrix is implausibly small")
      if (allocated(error)) return
      call check(error, maxval(abs(j_stored - transpose(j_stored))) < 1.0e-14_dp, &
                 "the Coulomb matrix is not symmetric")
      call mol%destroy()
      call aux%destroy()
   end subroutine compare

   subroutine test_spherical(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("def2-svp", "def2-universal-jkfit", .false., error)
   end subroutine test_spherical

   subroutine test_cartesian(error)
      !! 6-31G* fitted with itself: the one Cartesian pairing in the tree
      type(error_type), allocatable, intent(out) :: error
      call compare("6-31g*", "6-31g*", .true., error)
   end subroutine test_cartesian

   subroutine test_symmetric_part(error)
      !! An antisymmetric addition to the density changes nothing
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rij_t) :: rij
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), skewed(:, :), j_plain(:, :), j_skewed(:, :)
      integer :: mu, nu

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "def2-svp", mol, err)
      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "def2-universal-jkfit", aux, err)
      call build_rij(mol, aux, rij, err)
      call check(error,.not. err%has_error(), "setup failed: "//err%get_full_trace())
      if (allocated(error)) return
      call model_density(mol%nao, density)
      skewed = density
      do nu = 1, mol%nao
         do mu = 1, mol%nao
            skewed(mu, nu) = skewed(mu, nu) + 0.4_dp*sin(real(mu - nu, dp))
         end do
      end do
      call rij_coulomb(mol, rij, density, j_plain, err)
      call rij_coulomb(mol, rij, skewed, j_skewed, err)
      call check(error,.not. err%has_error(), "the Coulomb build failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(j_skewed - j_plain)) < 1.0e-11_dp, &
                 "an antisymmetric part of the density reached the Coulomb matrix")
      if (allocated(error)) return
      call check(error, maxval(abs(skewed - transpose(skewed))) > 0.1_dp, &
                 "the test density was not unsymmetric")
      call mol%destroy()
      call aux%destroy()
   end subroutine test_symmetric_part

   subroutine test_refusals(error)
      !! Mixed angular forms, a build that was never set up, a density of the wrong size
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux, cart
      type(rij_t) :: rij, empty
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), j(:, :)

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "def2-svp", mol, err)
      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "def2-universal-jkfit", aux, err)
      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "6-31g*", cart, err)
      call check(error,.not. err%has_error(), "the molecules should build")
      if (allocated(error)) return

      call build_rij(cart, aux, rij, err)
      call check(error, err%has_error(), "a Cartesian basis fitted with a spherical one was accepted")
      if (allocated(error)) return
      call err%clear()

      call model_density(mol%nao, density)
      call rij_coulomb(mol, empty, density, j, err)
      call check(error, err%has_error(), "a Coulomb build that was never set up was accepted")
      if (allocated(error)) return
      call err%clear()

      call build_rij(mol, aux, rij, err)
      call rij_coulomb(mol, rij, density(1:3, 1:3), j, err)
      call check(error, err%has_error(), "a density of the wrong size was accepted")
      call mol%destroy()
      call aux%destroy()
      call cart%destroy()
   end subroutine test_refusals

end module test_mqc_czt_rij

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_rij, only: collect_mqc_czt_rij_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_rij", collect_mqc_czt_rij_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
