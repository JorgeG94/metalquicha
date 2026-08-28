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
   implicit none
   private

   public :: collect_mqc_mp2_hessian_fd_tests

   !> pycc's pinned case: bohr, frame pinned. A displacement must not move the
   !> centre of mass or re-orient, or it stops matching analytic (bohr)
   !> integral derivatives.
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.000000_dp, 0.000000_dp, &
                           0.0_dp, 0.000000_dp, 1.814137_dp, &
                           0.0_dp, 1.756000_dp, -0.454300_dp], [3, 3])
   integer, parameter :: WATER_NELEC = 10

   !> The perturbed coordinate: atom 1 (O), z. Column index 3 in the flattened
   !> `(3*natom)` ordering this file uses, which is Fortran's -- pycc's column 2.
   integer, parameter :: PERT_ATOM = 1
   integer, parameter :: PERT_CART = 3

   !> Bohr. The 7-point O(h^6) stencil's step; pycc use the same.
   real(dp), parameter :: STEP = 2.0e-3_dp

   !> pycc's analytic correlation Hessian column `H[:, (0,2)]` for this case,
   !> 6-31G all-electron, on the BSE basis. Hartree/bohr^2, flattened
   !> `(x,y,z)` per atom.
   real(dp), parameter :: PYCC_COL(9) = [ &
                          -5.858801618662E-18_dp, -1.863531593298E-02_dp, -1.236323597338E-02_dp, &
                          5.244853948085E-18_dp, 4.451383450513E-03_dp, -3.435962860391E-03_dp, &
                          6.139476705757E-19_dp, 1.418393248247E-02_dp, 1.579919883377E-02_dp]

   !> The gate. Loose against the stencil's own reach, not against the
   !> derivative's: the step error here is far below this, and what the number
   !> has to absorb is the two codes' SCF convergence and integral rounding.
   real(dp), parameter :: TOL = 1.0e-8_dp

contains

   subroutine collect_mqc_mp2_hessian_fd_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("correlation_hessian_column_against_pycc", column_against_pycc), &
                  new_unittest("correlation_hessian_column_sums_to_nothing", column_translates) &
                  ]
   end subroutine collect_mqc_mp2_hessian_fd_tests

   !> The correlation gradient at one geometry: the MP2 total less the
   !> reference's own, over the same converged orbitals.
   subroutine correlation_gradient_at(coords, gradient, err)
      real(dp), intent(in) :: coords(:, :)
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: err

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

      call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, &
                                WATER_NELEC/2, total, err)
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

   !> `H[:, X] = d(correlation gradient)/dX` by the 7-point O(h^6) stencil
   !>
   !>   (-g(-3h) + 9g(-2h) - 45g(-h) + 45g(h) - 9g(2h) + g(3h)) / (60h)
   !>
   !> No `g(0)` term, which is the formula and not an omission.
   subroutine fd_column(column, err, ok)
      real(dp), intent(out) :: column(9)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      real(dp), parameter :: WEIGHT(6) = [-1.0_dp, 9.0_dp, -45.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]
      integer, parameter :: OFFSET(6) = [-3, -2, -1, 1, 2, 3]
      real(dp), allocatable :: gradient(:, :)
      real(dp) :: coords(3, 3), accumulated(3, 3)
      integer :: point

      ok = .false.
      accumulated = 0.0_dp
      do point = 1, 6
         coords = WATER
         coords(PERT_CART, PERT_ATOM) = coords(PERT_CART, PERT_ATOM) &
                                        + real(OFFSET(point), dp)*STEP
         call correlation_gradient_at(coords, gradient, err)
         if (err%has_error()) return
         accumulated = accumulated + WEIGHT(point)*gradient
      end do
      column = reshape(accumulated/(60.0_dp*STEP), [9])
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

   !> Translational invariance of the correlation energy, read off the column:
   !> summing the rows that share a Cartesian direction must vanish. Free, and
   !> weak -- a gradient in this tree was once wrong by 14 Hartree/bohr while
   !> leaving the net force at 8e-14 -- so it stands beside the comparison
   !> above and never in place of it.
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
