!! The MP2 correlation Hessian by difference of our own correlation gradient
module test_mqc_mp2_hessian_fd
   !! Unit 0.4 of the MP2 Hessian ladder (`mp2-hessian-phased-plan.md`): the
   !! finite-difference reference the analytic correlation Hessian will be
   !! checked against, and the first thing on that ladder to touch our Fortran.
   !!
   !! **Difference the gradient, not the energy.** A `1/h` stencil over an
   !! analytic first derivative is about three orders tighter than a `1/h^2`
   !! second difference of the energy, and the gradient being differenced here
   !! is the one `check_mp2_gradient` verifies -- so this rests on a derivative
   !! with an independent reference behind it rather than on its own internal
   !! consistency.
   !!
   !! **The correlation block alone.** `libcint_mp2_gradient` returns the total,
   !! so the correlation gradient is it minus `libcint_scf_gradient` over the
   !! same converged reference. That is what makes the column comparable to
   !! pycc's correlation Hessian, and it is what every later unit of the ladder
   !! works on. pycc difference exactly the same quantity for the same reason.
   !!
   !! **One column, not all nine.** pycc gives the cross-code reference in
   !! seconds; what a finite difference gives that pycc cannot is a check in our
   !! own conventions on our own gradient, and one column proves that as well as
   !! nine while costing six SCF+gradient evaluations instead of fifty-four.
   !!
   !! The reference literals are pycc's, regenerated on the basis *we* read
   !! (`basis_sets/6-31g.json`) rather than on psi4's internal tables, which
   !! carry two fewer significant figures and move this column by ~3e-9 -- over
   !! the tolerance here. Recipe and provenance:
   !! `tools/mp2_hessian_oracle/`, and the columns in `hess_col_bse_631g.json`.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_gradient, only: libcint_scf_gradient
   use mqc_libcint_mp2_gradient, only: libcint_mp2_gradient
   use mqc_libcint_mp2_hessian, only: mp2_correlation_hessian
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none
   private

   public :: collect_mqc_mp2_hessian_fd_tests

   !! pycc's pinned case: bohr, frame pinned. A displacement must not move the
   !! centre of mass or re-orient, or it stops matching analytic (bohr)
   !! integral derivatives.
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.000000_dp, 0.000000_dp, &
                           0.0_dp, 0.000000_dp, 1.814137_dp, &
                           0.0_dp, 1.756000_dp, -0.454300_dp], [3, 3])
   integer, parameter :: WATER_NELEC = 10

   !! The perturbed coordinate: atom 1 (O), z. Column index 3 in the flattened
   !! `(3*natom)` ordering this file uses, which is Fortran's -- pycc's column 2.
   integer, parameter :: PERT_ATOM = 1
   integer, parameter :: PERT_CART = 3

   !! Bohr. The 7-point O(h^6) stencil's step; pycc use the same.
   real(dp), parameter :: STEP = 2.0e-3_dp

   !! pycc's analytic correlation Hessian column `H[:, (0,2)]` for this case,
   !! 6-31G all-electron, on the BSE basis. Hartree/bohr^2, flattened
   !! `(x,y,z)` per atom.
   real(dp), parameter :: PYCC_COL(9) = [ &
                          -5.858801618662E-18_dp, -1.863531593298E-02_dp, -1.236323597338E-02_dp, &
                          5.244853948085E-18_dp, 4.451383450513E-03_dp, -3.435962860391E-03_dp, &
                          6.139476705757E-19_dp, 1.418393248247E-02_dp, 1.579919883377E-02_dp]

   !! The gate. Loose against the stencil's own reach, not against the
   !! derivative's: the step error here is far below this, and what the number
   !! has to absorb is the two codes' SCF convergence and integral rounding.
   real(dp), parameter :: TOL = 1.0e-8_dp

contains

   subroutine collect_mqc_mp2_hessian_fd_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("correlation_hessian_column_against_pycc", column_against_pycc), &
                  new_unittest("correlation_hessian_column_sums_to_nothing", column_translates), &
                  new_unittest("analytic_column_against_the_same_literals", analytic_column), &
                  new_unittest("frozen_core_column_against_finite_difference", frozen_column) &
                  ]
   end subroutine collect_mqc_mp2_hessian_fd_tests

   !! The correlation gradient at one geometry: the MP2 total less the
   !! reference's own, over the same converged orbitals.
   subroutine correlation_gradient_at(coords, gradient, err, n_frozen)
      real(dp), intent(in) :: coords(:, :)
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: err
      integer, intent(in), optional :: n_frozen

      integer :: use_frozen

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: total(:, :), reference(:, :)

      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "6-31g", mol, err)
      if (err%has_error()) return

      ! Tighter than the default on purpose: a loose SCF is the usual cause of
      ! a noisy finite-difference column, and it is the cheapest thing to rule
      ! out before doubting the derivative.
      call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if

      use_frozen = 0
      if (present(n_frozen)) use_frozen = n_frozen
      call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, &
                                WATER_NELEC/2, total, err, n_frozen=use_frozen)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if

      call libcint_scf_gradient(mol, density=scf%density, orbitals=scf%orbitals, &
                                orbital_energies=scf%orbital_energies, &
                                n_occupied=WATER_NELEC/2, gradient=reference, error=err)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if

      gradient = total - reference
      call mol%destroy()
   end subroutine correlation_gradient_at

   !! `H[:, X] = d(correlation gradient)/dX` by the 7-point O(h^6) stencil
   !!
   !!   (-g(-3h) + 9g(-2h) - 45g(-h) + 45g(h) - 9g(2h) + g(3h)) / (60h)
   !!
   !! No `g(0)` term, which is the formula and not an omission.
   subroutine fd_column(column, err, ok, n_frozen)
      real(dp), intent(out) :: column(9)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok
      integer, intent(in), optional :: n_frozen

      real(dp), parameter :: WEIGHT(6) = [-1.0_dp, 9.0_dp, -45.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]
      integer, parameter :: OFFSET(6) = [-3, -2, -1, 1, 2, 3]
      real(dp), allocatable :: gradient(:, :)
      real(dp) :: coords(3, 3), accumulated(3, 3)
      integer :: point, threads

      ! One thread, restored afterwards. The stencil divides by 60h = 0.12, so
      ! it multiplies gradient-level noise by about eight: measured here, the
      ! residual is 2.778e-10 and bit-reproducible at one thread but scatters
      ! between 8e-10 and 3.1e-9 over repeat runs at four, because the
      ! unordered OpenMP critical merges do not sum in a fixed order. A gate
      ! this tight has to be handed a deterministic number.
      threads = omp_get_max_threads()
      call omp_set_num_threads(1)

      ok = .false.
      accumulated = 0.0_dp
      do point = 1, 6
         coords = WATER
         coords(PERT_CART, PERT_ATOM) = coords(PERT_CART, PERT_ATOM) &
                                        + real(OFFSET(point), dp)*STEP
         call correlation_gradient_at(coords, gradient, err, n_frozen=n_frozen)
         if (err%has_error()) then
            call omp_set_num_threads(threads)
            return
         end if
         accumulated = accumulated + WEIGHT(point)*gradient
      end do
      column = reshape(accumulated/(60.0_dp*STEP), [9])
      call omp_set_num_threads(threads)
      ok = .true.
   end subroutine fd_column

   subroutine column_against_pycc(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: column(9)
      logical :: ok

      call fd_column(column, err, ok)
      call check(error, ok, "the correlation gradient did not evaluate: "//err%get_message())
      if (allocated(error)) return

      ! Printed rather than only asserted: a tolerance that passes says
      ! nothing about how much room was left, and this column is the reference
      ! every later unit of the ladder is measured against.
      write (*, "(a, es10.3)") "        max |ours - pycc| = ", &
         maxval(abs(column - PYCC_COL))
      call check(error, maxval(abs(column - PYCC_COL)) < TOL, &
                 "correlation Hessian column disagrees with pycc")
   end subroutine column_against_pycc

   !! Translational invariance of the correlation energy, read off the column:
   !! summing the rows that share a Cartesian direction must vanish. Free, and
   !! weak -- a gradient in this tree was once wrong by 14 Hartree/bohr while
   !! leaving the net force at 8e-14 -- so it stands beside the comparison
   !! above and never in place of it.
   subroutine column_translates(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: column(9)
      logical :: ok
      integer :: cart

      call fd_column(column, err, ok)
      call check(error, ok, "the correlation gradient did not evaluate: "//err%get_message())
      if (allocated(error)) return

      write (*, "(a, 3es10.3)") "        sum over atoms   = ", &
         (sum(column(cart::3)), cart=1, 3)
      do cart = 1, 3
         call check(error, abs(sum(column(cart::3))) < 1.0e-8_dp, &
                    "the correlation Hessian column does not sum to zero")
         if (allocated(error)) return
      end do
   end subroutine column_translates

   !! Gate 1.9b, closed transitively through the literals both sides of it
   !! are measured against: the finite-difference column above sits 2.778e-10
   !! from `PYCC_COL` (bit-reproducible at one thread), so the analytic
   !! column agreeing with the same literals ties analytic to FD without
   !! re-running the stencil -- which is the difference between a seconds
   !! test and a minutes one. Measured 1.096e-12 when the assembly landed;
   !! the shared `TOL` is the stencil comparison's and leaves three orders.
   subroutine analytic_column(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: hess_corr(:, :, :, :), hess_ref(:, :, :, :)
      real(dp) :: column(9)
      integer :: threads, ia, a

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "6-31g", mol, err)
      if (.not. err%has_error()) then
         call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                              .false., scf, err)
      end if
      if (.not. err%has_error()) then
         call mp2_correlation_hessian(mol, scf%orbitals, scf%orbital_energies, &
                                      scf%density, WATER_NELEC/2, 0, &
                                      hess_corr, hess_ref, err)
      end if
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the analytic Hessian did not evaluate: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      do ia = 1, 3
         do a = 1, 3
            column(3*(ia - 1) + a) = hess_corr(a, PERT_CART, ia, PERT_ATOM)
         end do
      end do
      write (*, "(a, es10.3)") "        max |analytic - pycc| = ", &
         maxval(abs(column - PYCC_COL))
      call check(error, maxval(abs(column - PYCC_COL)) < TOL, &
                 "the analytic correlation Hessian column disagrees with the "// &
                 "literals the FD column is measured against")
      call mol%destroy()
   end subroutine analytic_column

   !! The frozen-core tie-out, direct rather than through literals: the
   !! same 7-point stencil over our own frozen-core correlation gradient
   !! (n_frozen = 1) against the same column of the analytic frozen-core
   !! assembly. Entirely in our own conventions -- the cross-code element-wise
   !! gates ran when Phase 2 landed and live in the commit messages; what
   !! this pins is that the analytic core rotations keep differentiating the
   !! gradient actually shipped, and it is exactly the check that catches a
   !! quotient-rule Sylvester block (a ~7e-7 signature against a 1e-8 gate).
   subroutine frozen_column(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: hess_corr(:, :, :, :), hess_ref(:, :, :, :)
      real(dp) :: fd(9), analytic(9)
      logical :: ok
      integer :: threads, ia, a

      call fd_column(fd, err, ok, n_frozen=1)
      call check(error, ok, "the frozen-core correlation gradient did not "// &
                 "evaluate: "//err%get_message())
      if (allocated(error)) return

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "6-31g", mol, err)
      if (.not. err%has_error()) then
         call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                              .false., scf, err)
      end if
      if (.not. err%has_error()) then
         call mp2_correlation_hessian(mol, scf%orbitals, scf%orbital_energies, &
                                      scf%density, WATER_NELEC/2, 1, &
                                      hess_corr, hess_ref, err)
      end if
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the frozen-core analytic Hessian did not evaluate: "// &
                 err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      do ia = 1, 3
         do a = 1, 3
            analytic(3*(ia - 1) + a) = hess_corr(a, PERT_CART, ia, PERT_ATOM)
         end do
      end do
      write (*, "(a, es10.3)") "        max |analytic - FD|   = ", &
         maxval(abs(analytic - fd))
      call check(error, maxval(abs(analytic - fd)) < TOL, &
                 "the frozen-core analytic column disagrees with the finite "// &
                 "difference of the frozen-core gradient")
      call mol%destroy()
   end subroutine frozen_column

end module test_mqc_mp2_hessian_fd

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mp2_hessian_fd, only: collect_mqc_mp2_hessian_fd_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_mp2_hessian_fd", collect_mqc_mp2_hessian_fd_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
