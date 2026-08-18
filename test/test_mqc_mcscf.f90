!! The orbital gradient, and the optimisation that drives it to zero
module test_mqc_mcscf
   !! CASSCF is CASCI with the orbitals allowed to move. What has to be right is
   !! the derivative of the energy with respect to moving them, and there is
   !! exactly one way to establish that without circularity: differentiate the
   !! energy numerically and compare.
   !!
   !! **The finite-difference test is the load-bearing one here.** Everything
   !! else -- that the optimiser converges, that it reaches PySCF's energy --
   !! can be true of a subtly wrong gradient, because an optimiser will happily
   !! find a stationary point of whatever function its gradient actually
   !! describes. And a gradient with the wrong *sign* does not look like a sign
   !! error; it looks like a convergence problem, and gets diagnosed as one for
   !! a long time. That sign was in fact wrong when this was first written, and
   !! this is what caught it.
   !!
   !! The reference SCFs run to 1e-12 in the energy and 1e-10 in the density,
   !! matching the rest of the suite. An earlier version asked for 1e-13, which
   !! on a -75 hartree energy is 1.3e-15 relative -- the last bit of a double --
   !! and the threaded Fock builds are not bit-reproducible, so it reached it
   !! only sometimes. The result was a test that passed four, five or six of its
   !! six cases depending on the run. Nothing downstream needs those digits:
   !! the energies here are compared at 1e-9.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_casci, only: run_libcint_casci, casci_result_t
   use mqc_determinants, only: link_table_t, build_link_table
   use mqc_rdm, only: active_space_rdms
   use mqc_libcint_mcscf, only: mcscf_fock_t, generalized_fock, orbital_gradient, &
                                rotation_matrix, is_redundant, run_libcint_casscf, &
                                casscf_result_t
   implicit none
   private

   public :: collect_mqc_mcscf_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

   real(dp), parameter :: NITROGEN(3, 2) = reshape( &
                          [0.0_dp, 0.0_dp, -1.0371_dp, &
                           0.0_dp, 0.0_dp, 1.0371_dp], [3, 2])
   integer, parameter :: NITROGEN_Z(2) = [7, 7]
   character(len=2), parameter :: NITROGEN_SYM(2) = ["N ", "N "]

contains

   subroutine collect_mqc_mcscf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rotation_is_orthogonal", test_rotation), &
                  new_unittest("redundant_blocks_are_excluded", test_redundant), &
                  new_unittest("gradient_against_finite_differences", test_gradient), &
                  new_unittest("nitrogen_casscf_against_pyscf", test_nitrogen), &
                  new_unittest("water_casscf_against_pyscf", test_water), &
                  new_unittest("casscf_improves_on_casci", test_variational) &
                  ]
   end subroutine collect_mqc_mcscf_tests

   subroutine test_rotation(error)
      !! `exp(kappa)` is orthogonal, for small and for large steps
      !!
      !! This is the reason orbital changes are parametrised as an exponential
      !! rather than added directly: an orthogonal transformation of orthonormal
      !! orbitals is orthonormal whatever the step size, so a bad step is never
      !! an invalid one and nothing has to be reorthogonalised afterwards. The
      !! large-angle case is the one that matters -- a truncated series would
      !! pass the small one and quietly fail here.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: kappa(5, 5)
      real(dp), allocatable :: rotation(:, :), product(:, :)
      integer :: i, j, case_index
      real(dp), parameter :: SIZES(3) = [1.0e-6_dp, 0.2_dp, 3.0_dp]

      do case_index = 1, 3
         kappa = 0.0_dp
         do j = 1, 5
            do i = j + 1, 5
               kappa(i, j) = SIZES(case_index)/real(i + j, dp)
               kappa(j, i) = -kappa(i, j)
            end do
         end do

         call rotation_matrix(kappa, rotation)
         allocate (product(5, 5))
         call pic_gemm(rotation, rotation, product, transa="T")
         do i = 1, 5
            product(i, i) = product(i, i) - 1.0_dp
         end do
         call check(error, maxval(abs(product)) < 1.0e-13_dp, &
                    "the rotation should be orthogonal whatever the step size")
         deallocate (product, rotation)
         if (allocated(error)) return
      end do
   end subroutine test_rotation

   subroutine test_redundant(error)
      !! Rotations that change nothing are recognised as such
      type(error_type), allocatable, intent(out) :: error

      ! Two inactive, two active, two virtual: orbitals 1-2, 3-4, 5-6.
      call check(error, is_redundant(1, 2, 2, 2), "inactive with inactive")
      if (allocated(error)) return
      call check(error, is_redundant(3, 4, 2, 2), "active with active, in a "// &
                 "complete active space")
      if (allocated(error)) return
      call check(error, is_redundant(5, 6, 2, 2), "virtual with virtual")
      if (allocated(error)) return
      call check(error,.not. is_redundant(3, 1, 2, 2), "active with inactive")
      if (allocated(error)) return
      call check(error,.not. is_redundant(5, 1, 2, 2), "virtual with inactive")
      if (allocated(error)) return
      call check(error,.not. is_redundant(5, 3, 2, 2), "virtual with active")
   end subroutine test_redundant

   subroutine test_gradient(error)
      !! The analytic gradient against central differences of the CASCI energy
      !!
      !! The orbitals are rotated away from the SCF solution first, on purpose.
      !! At the Hartree-Fock orbitals several blocks of the gradient are already
      !! stationary and come out at 1e-15, which agrees with a finite difference
      !! of zero no matter what the sign convention is -- so testing there would
      !! establish nothing about the thing most likely to be wrong.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casci_result_t) :: ci
      type(mcscf_fock_t) :: fock
      type(link_table_t) :: alpha, beta
      type(error_t) :: err
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :), gradient(:, :)
      real(dp), allocatable :: start(:, :), moved(:, :), kappa(:, :), rotation(:, :)
      real(dp) :: plus, minus, difference
      integer :: n_ao, n_mo, i, j, pair
      integer, parameter :: N_INACTIVE = 3, N_ACTIVE = 4
      integer, parameter :: LEFT(4) = [4, 5, 6, 7], RIGHT(4) = [1, 2, 3, 1]
      real(dp), parameter :: STEP = 1.0e-4_dp

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call run_libcint_rhf(mol, 10, 300, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return
      n_ao = size(scf%orbitals, 1)
      n_mo = size(scf%orbitals, 2)

      allocate (start(n_ao, n_mo), moved(n_ao, n_mo), kappa(n_mo, n_mo))
      kappa = 0.0_dp
      do j = 1, n_mo
         do i = j + 1, n_mo
            kappa(i, j) = 0.03_dp/real(i + j, dp)
            kappa(j, i) = -kappa(i, j)
         end do
      end do
      call rotation_matrix(kappa, rotation)
      call pic_gemm(scf%orbitals, rotation, start)
      deallocate (rotation)

      call run_libcint_casci(mol, start, N_INACTIVE, N_ACTIVE, 2, 2, ci, err, &
                             tolerance=1.0e-12_dp)
      call build_link_table(N_ACTIVE, 2, alpha, err)
      call build_link_table(N_ACTIVE, 2, beta, err)
      call active_space_rdms(ci%ci_vector, alpha, beta, dm1, dm2, err)
      call generalized_fock(mol, start, N_INACTIVE, N_ACTIVE, dm1, dm2, fock, err)
      call orbital_gradient(fock, N_INACTIVE, N_ACTIVE, gradient)
      call check(error,.not. err%has_error(), "the gradient should build")
      if (allocated(error)) return

      ! Something must actually be non-zero, or the comparison below is vacuous.
      call check(error, maxval(abs(gradient)) > 1.0e-3_dp, &
                 "the rotated orbitals should have a gradient worth differencing")
      if (allocated(error)) return

      do pair = 1, 4
         kappa = 0.0_dp
         kappa(LEFT(pair), RIGHT(pair)) = STEP
         kappa(RIGHT(pair), LEFT(pair)) = -STEP
         call rotation_matrix(kappa, rotation)
         call pic_gemm(start, rotation, moved)
         call run_libcint_casci(mol, moved, N_INACTIVE, N_ACTIVE, 2, 2, ci, err, &
                                tolerance=1.0e-12_dp)
         plus = ci%energy
         deallocate (rotation)

         kappa = -kappa
         call rotation_matrix(kappa, rotation)
         call pic_gemm(start, rotation, moved)
         call run_libcint_casci(mol, moved, N_INACTIVE, N_ACTIVE, 2, 2, ci, err, &
                                tolerance=1.0e-12_dp)
         minus = ci%energy
         deallocate (rotation)

         difference = (plus - minus)/(2.0_dp*STEP)
         ! A central difference with this step carries a truncation error of
         ! order STEP squared, so 1e-7 is the tightest this can be asked for.
         call check(error, abs(gradient(LEFT(pair), RIGHT(pair)) - difference) &
                    < 1.0e-7_dp, "the analytic gradient should equal the numerical "// &
                    "derivative of the CASCI energy")
         if (allocated(error)) return
      end do

      ! And the redundant directions really are flat.
      call check(error, abs(gradient(6, 5)) < 1.0e-14_dp, &
                 "a rotation between two active orbitals should have no gradient")
      if (allocated(error)) return
      call check(error, maxval(abs(gradient + transpose(gradient))) < 1.0e-12_dp, &
                 "the gradient should be antisymmetric")

      call alpha%destroy()
      call beta%destroy()
      call mol%destroy()
   end subroutine test_gradient

   subroutine run_casscf(atomic_numbers, symbols, coordinates, basis, n_electrons, &
                         n_inactive, n_active, n_alpha, n_beta, cycles, &
                         casci_energy, result, err, ok)
      !! An SCF, a CASCI on its orbitals, then a CASSCF from the same start
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: n_electrons, n_inactive, n_active, n_alpha, n_beta, cycles
      real(dp), intent(out) :: casci_energy
      type(casscf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casci_result_t) :: ci

      ok = .false.
      casci_energy = 0.0_dp
      call build_libcint_molecule(atomic_numbers, symbols, coordinates, basis, mol, err)
      call run_libcint_rhf(mol, n_electrons, 300, 1.0e-12_dp, 1.0e-10_dp, .false., &
                           scf, err)
      if (err%has_error() .or. .not. scf%converged) return
      call run_libcint_casci(mol, scf%orbitals, n_inactive, n_active, n_alpha, n_beta, &
                             ci, err, tolerance=1.0e-12_dp)
      if (err%has_error()) return
      casci_energy = ci%energy

      call run_libcint_casscf(mol, scf%orbitals, n_inactive, n_active, n_alpha, n_beta, &
                              result, err, max_iterations=cycles, &
                              gradient_tol=1.0e-6_dp)
      ok = .not. err%has_error()
      call mol%destroy()
   end subroutine run_casscf

   subroutine test_nitrogen(error)
      !! N2, cc-pVDZ, CAS(6,6)
      !!
      !! The textbook active space: six electrons in the three bonding and three
      !! antibonding orbitals of the triple bond. Twenty-one macro-iterations
      !! with DIIS, thirty-three without, which is what a well-conditioned case
      !! looks like for a first-order method.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(casscf_result_t) :: result
      real(dp) :: casci_energy
      logical :: ok

      real(dp), parameter :: REFERENCE = -109.090018071393_dp

      call run_casscf(NITROGEN_Z, NITROGEN_SYM, NITROGEN, "cc-pvdz", 14, 4, 6, 3, 3, &
                      120, casci_energy, result, err, ok)
      call check(error, ok, "the calculation should succeed")
      if (allocated(error)) return
      call check(error, result%converged, "the orbital gradient should reach the "// &
                 "threshold")
      if (allocated(error)) return
      call check(error, result%energy, REFERENCE, "the CASSCF energy against PySCF", &
                 thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, result%gradient_norm < 1.0e-6_dp, &
                 "and the gradient it stopped on should really be small")
   end subroutine test_nitrogen

   subroutine test_water(error)
      !! Water, cc-pVDZ, CAS(4,4)
      !!
      !! Slower than nitrogen by a factor of three, and the reason is worth
      !! recording rather than hiding behind a larger iteration cap: the energy
      !! surface here has a long flat valley, and a first-order optimiser crawls
      !! along it -- iterations eleven to twenty-five gain a few tens of
      !! microhartree between them, taking steps well inside the trust radius.
      !! DIIS brings it from about 102 macro-iterations to about 60; a
      !! second-order method would not care at all. This is the case that says
      !! what the two-step algorithm costs.
      !!
      !! The count varies by ten or so between runs at identical settings,
      !! because the threaded Fock builds are not bit-reproducible and the
      !! valley is flat enough to amplify that. The iteration cap is set well
      !! above it for that reason; it is not a convergence criterion.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(casscf_result_t) :: result
      real(dp) :: casci_energy
      logical :: ok

      real(dp), parameter :: REFERENCE = -76.077846973211_dp

      call run_casscf(WATER_Z, WATER_SYM, WATER, "cc-pvdz", 10, 3, 4, 2, 2, &
                      250, casci_energy, result, err, ok)
      call check(error, ok, "the calculation should succeed")
      if (allocated(error)) return
      call check(error, result%converged, "it should converge, given the iterations")
      if (allocated(error)) return
      call check(error, result%energy, REFERENCE, "the CASSCF energy against PySCF", &
                 thr=1.0e-9_dp)
   end subroutine test_water

   subroutine test_variational(error)
      !! Optimising the orbitals lowers the energy, and the density is sane
      !!
      !! CASSCF minimises the same energy CASCI evaluates, over a strictly
      !! larger set of variables, so it cannot come out higher. That is worth
      !! asserting because it is the one property that needs no reference at all
      !! and it is exactly what a sign error in the gradient breaks -- an
      !! optimiser given the negative of its gradient climbs, and reports
      !! convergence at a maximum with every appearance of success.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(casscf_result_t) :: result
      real(dp) :: casci_energy, trace
      logical :: ok
      integer :: t

      call run_casscf(WATER_Z, WATER_SYM, WATER, "sto-3g", 10, 2, 5, 3, 3, &
                      200, casci_energy, result, err, ok)
      call check(error, ok, "the calculation should succeed")
      if (allocated(error)) return
      call check(error, result%energy < casci_energy, &
                 "CASSCF minimises over more variables than CASCI and cannot come "// &
                 "out above it")
      if (allocated(error)) return

      ! -75.009717745820 from pyscf.mcscf.CASSCF on the same active space.
      call check(error, result%energy, -75.009717745820_dp, &
                 "water STO-3G CAS(6,5) against PySCF", thr=1.0e-9_dp)
      if (allocated(error)) return

      ! The active occupations have to be occupations: between zero and two,
      ! summing to the active electron count.
      trace = 0.0_dp
      do t = 1, size(result%dm1, 1)
         trace = trace + result%dm1(t, t)
         call check(error, result%dm1(t, t) > -1.0e-10_dp .and. &
                    result%dm1(t, t) < 2.0_dp + 1.0e-10_dp, &
                    "every natural occupation should lie between zero and two")
         if (allocated(error)) return
      end do
      call check(error, trace, 6.0_dp, "and they should sum to the active electrons", &
                 thr=1.0e-9_dp)
   end subroutine test_variational

end module test_mqc_mcscf

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mcscf, only: collect_mqc_mcscf_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_mcscf", collect_mqc_mcscf_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
