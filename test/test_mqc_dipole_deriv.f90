!! The nuclear derivative of the dipole integrals, against a difference of them
module test_mqc_dipole_deriv
   !! What an analytic infrared intensity rests on, one rung down
   !!
   !! `dipole_integral_derivatives` reads `int1e_irp` -- nine components that
   !! are `<u| r_a d_b |v>` in one of two orders. Both orders return a matrix of
   !! the right shape and the right magnitude, and picking the wrong one gives a
   !! dipole derivative transposed in its two Cartesian indices. That error
   !! survives the translational sum rule, because the sum runs over atoms and
   !! not over the component pair, so nothing cheap catches it.
   !!
   !! So the layout is measured here rather than taken from a convention: the
   !! same integrals are differenced numerically, and the analytic array is
   !! compared against that both as it stands and with its two Cartesian indices
   !! swapped. The first has to match and the second has to fail, or the test
   !! cannot tell the two readings apart and is not doing its job.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_multipole, only: multipole_matrices, dipole_integral_derivatives
   implicit none
   private

   public :: collect_mqc_dipole_deriv_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.000000_dp, 0.000000_dp, &
                           0.0_dp, 0.000000_dp, 1.814137_dp, &
                           0.0_dp, 1.756000_dp, -0.454300_dp], [3, 3])

   !> Deliberately not the centre of mass and not zero: an origin that happens
   !> to be special can hide an index that is being read from the wrong place.
   real(dp), parameter :: ORIGIN(3) = [0.13_dp, -0.07_dp, 0.21_dp]

   real(dp), parameter :: STEP = 1.0e-4_dp

contains

   subroutine collect_mqc_dipole_deriv_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("derivatives_difference_the_dipole_integrals", &
                               test_against_difference), &
                  new_unittest("d_functions_difference_the_dipole_integrals", &
                               test_against_difference_d), &
                  new_unittest("the_transposed_reading_is_distinguishable", &
                               test_transpose_fails) &
                  ]
   end subroutine collect_mqc_dipole_deriv_tests

   !> Central difference of `multipole_matrices` for one atom and direction
   subroutine fd_block(atom, cart, basis, block, err, ok)
      integer, intent(in) :: atom, cart
      character(len=*), intent(in) :: basis
      real(dp), allocatable, intent(out) :: block(:, :, :)   !! (n_ao, n_ao, 3)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      real(dp), allocatable :: plus(:, :, :), minus(:, :, :)
      real(dp) :: coords(3, 3)

      ok = .false.
      coords = WATER
      coords(cart, atom) = coords(cart, atom) + STEP
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, basis, mol, err)
      if (err%has_error()) return
      call multipole_matrices(mol, ORIGIN, 1, plus, err)
      call mol%destroy()
      if (err%has_error()) return

      coords = WATER
      coords(cart, atom) = coords(cart, atom) - STEP
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, basis, mol, err)
      if (err%has_error()) return
      call multipole_matrices(mol, ORIGIN, 1, minus, err)
      call mol%destroy()
      if (err%has_error()) return

      block = (plus - minus)/(2.0_dp*STEP)
      ok = .true.
   end subroutine fd_block

   subroutine worst_errors(basis, direct, transposed, err, ok)
      !! `max |analytic - fd|` for the array as read, and with `a` and `b` swapped
      character(len=*), intent(in) :: basis
      real(dp), intent(out) :: direct, transposed
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      real(dp), allocatable :: ddip(:, :, :, :, :), fd(:, :, :)
      integer :: atom, cart, a

      ok = .false.
      direct = 0.0_dp
      transposed = 0.0_dp

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, basis, mol, err)
      if (err%has_error()) return
      call dipole_integral_derivatives(mol, ORIGIN, ddip, err)
      call mol%destroy()
      if (err%has_error()) return

      do atom = 1, 3
         do cart = 1, 3
            call fd_block(atom, cart, basis, fd, err, ok)
            if (.not. ok) return
            do a = 1, 3
               direct = max(direct, maxval(abs(ddip(:, :, a, cart, atom) - fd(:, :, a))))
               ! The other reading of the nine components: dipole direction and
               ! gradient direction exchanged.
               transposed = max(transposed, &
                                maxval(abs(ddip(:, :, cart, a, atom) - fd(:, :, a))))
            end do
         end do
      end do
      ok = .true.
   end subroutine worst_errors

   subroutine test_against_difference(error)
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: direct, transposed
      type(error_t) :: err
      logical :: ok

      call worst_errors("6-31g", direct, transposed, err, ok)
      call check(error, ok, "the dipole integral derivatives did not evaluate: "// &
                 err%get_message())
      if (allocated(error)) return

      write (*, "(a, es10.3)") "        max |analytic - fd| = ", direct
      ! The difference's own floor, bounded for a toolchain rather than a box:
      ! measured 6.9e-10 here, and the same comparison one rung up varies by a
      ! factor of twenty between this machine and CI. A bound at ten times one
      ! measurement does not survive that; this one does, and still sits eight
      ! orders below the transposed-layout error the next case pins.
      call check(error, direct < 5.0e-8_dp, &
                 "the dipole integral derivatives do not difference the dipole "// &
                 "integrals")
   end subroutine test_against_difference

   subroutine test_against_difference_d(error)
      !! The same on a basis with d functions
      !!
      !! `int1e_irp` at `l = 2` is a different code path in the integral
      !! library, and the layout and the assembly were both wrong at s and p
      !! before they were measured, so neither is assumed to carry over.
      !! 6-31G* is Cartesian in our convention, making this 6d rather than 5d.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: direct, transposed
      type(error_t) :: err
      logical :: ok

      call worst_errors("6-31g*", direct, transposed, err, ok)
      call check(error, ok, "the dipole integral derivatives did not evaluate: "// &
                 err%get_message())
      if (allocated(error)) return

      write (*, "(a, es10.3)") "        6-31g*  max |analytic - fd| = ", direct
      ! Measured 5.476e-09, larger than the s-and-p case's 1.543e-09 for the
      ! stencil's reason rather than the integrals': d functions are sharper, so
      ! the same step resolves them less well. Its own bound rather than the one
      ! above, and both now carry room for the twenty-fold spread between
      ! machines that this ladder's tighter cases turned out to have.
      call check(error, direct < 2.0e-7_dp, &
                 "the dipole integral derivatives do not difference the dipole "// &
                 "integrals on a basis with d functions")
   end subroutine test_against_difference_d

   subroutine test_transpose_fails(error)
      !! The check above can tell the two component orders apart
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: direct, transposed
      type(error_t) :: err
      logical :: ok

      call worst_errors("6-31g", direct, transposed, err, ok)
      call check(error, ok, "the dipole integral derivatives did not evaluate: "// &
                 err%get_message())
      if (allocated(error)) return

      write (*, "(a, es10.3)") "        max |transposed - fd| = ", transposed
      call check(error, transposed > 1.0e-3_dp, &
                 "reading int1e_irp with its two Cartesian indices swapped is "// &
                 "indistinguishable here, so the check above proves nothing")
   end subroutine test_transpose_fails

end module test_mqc_dipole_deriv

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_dipole_deriv, only: collect_mqc_dipole_deriv_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_dipole_deriv", collect_mqc_dipole_deriv_tests)]
   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
