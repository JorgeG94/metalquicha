!! The range-separated gradient, against a difference of its own energy
module test_mqc_rsh_gradient
   !! A range-separated functional splits its exchange over an error-function
   !! kernel, so its gradient needs two exchange derivatives rather than one:
   !! the full-range pass every hybrid takes, and a second at the screened
   !! omega. Only the density-fitted path built the second one for a while;
   !! the exact-ERI path refused outright rather than return a gradient
   !! missing wB97X's dominant exchange term.
   !!
   !! Differencing the total energy is the right check for that second pass.
   !! It needs no reference implementation, it is sensitive to the coefficient
   !! and the sign as well as to the algebra, and -- unlike the Hessian tests
   !! next door -- nothing here is held fixed, so the comparison is exact up
   !! to the step error and whatever the quadrature contributes.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_libcint_gradient, only: libcint_scf_gradient
   implicit none
   private

   public :: collect_mqc_rsh_gradient

   ! Water, bent, in Bohr -- the same geometry the Hessian tests use, so a
   ! number from one is comparable with a number from the other.
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape([ &
                                                0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                0.0000_dp, 1.4300_dp, 1.1075_dp, &
                                                0.0000_dp, -1.4300_dp, 1.1075_dp], [3, 3])

contains

   subroutine collect_mqc_rsh_gradient(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("wb97x-gradient", test_wb97x), &
                  new_unittest("cam-b3lyp-gradient", test_cam_b3lyp), &
                  new_unittest("global-hybrid-unchanged", test_b3lyp_unchanged) &
                  ]
   end subroutine collect_mqc_rsh_gradient

   subroutine test_wb97x(error)
      !! The long-range coefficient is 1.0 here, so the second pass is not a
      !! correction to the first -- it is most of the exchange gradient.
      type(error_type), allocatable, intent(out) :: error
      call gradient_against_energy("wb97x", error)
   end subroutine test_wb97x

   subroutine test_cam_b3lyp(error)
      !! 0.19 short range and 0.65 long, so a pass dropped here is a smaller
      !! error than in wB97X and a coefficient confused between the two is a
      !! larger one.
      type(error_type), allocatable, intent(out) :: error
      call gradient_against_energy("cam-b3lyp", error)
   end subroutine test_cam_b3lyp

   subroutine test_b3lyp_unchanged(error)
      !! A global hybrid takes no second pass at all. Here to catch an omega
      !! left set on a shared `env`, which would attenuate a functional that
      !! never asked for it and is exactly the failure the local copy exists
      !! to prevent.
      type(error_type), allocatable, intent(out) :: error
      call gradient_against_energy("b3lyp", error)
   end subroutine test_b3lyp_unchanged

   subroutine gradient_against_energy(functional, error)
      character(len=*), intent(in) :: functional
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-3_dp
      real(dp) :: shifted(3, 3), fd, worst, ep, em, net(3)
      real(dp), allocatable :: gradient(:, :)
      type(error_t) :: err
      logical :: ok
      integer :: ia, a

      if (.not. xc_available()) return

      call gradient_at(WATER, functional, gradient, err, ok)
      call check(error, ok, "the gradient failed for "//functional)
      if (allocated(error)) return

      ! Every basis function moves with an atom, so the forces sum to zero
      ! whatever the functional. Cheap, and it fails loudly if the second pass
      ! is added to the wrong atom's block.
      do a = 1, 3
         net(a) = sum(gradient(a, :))
      end do
      call check(error, maxval(abs(net)) < 1.0e-8_dp, &
                 "the gradient has a net force for "//functional)
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            call energy_at(shifted, functional, ep, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced energy failed")
               return
            end if
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            call energy_at(shifted, functional, em, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced energy failed")
               return
            end if
            fd = (ep - em)/(2.0_dp*H)
            worst = max(worst, abs(gradient(a, ia) - fd))
         end do
      end do

      ! Bounded by the difference, not by the gradient. All three land
      ! 1.47e-7 from it, to three figures -- that sameness is the point: it is
      ! the h^2 step error, so the range-separated gradients are as good as the
      ! global hybrid's rather than merely close enough. Dropping the
      ! long-range pass entirely misses by 1e-1 on wB97X and a coefficient
      ! confused between alpha and alpha+beta by 1e-2, so the tolerance
      ! separates right from wrong by several decades either way.
      !
      ! Against PySCF rather than against itself, at grid level 6: cam-B3LYP
      ! agrees to 4.1e-8 and B3LYP to 6.3e-8. wB97X is 3.9e-6, and that is the
      ! functional and not this code -- its semilocal part is stiff enough that
      ! PySCF's own wB97X energy is still moving in the seventh decimal between
      ! grid levels 6 and 8. Its disagreement tracks the grid, 1.2e-4 at level
      ! 3 and 9.0e-7 at level 8, while cam-B3LYP -- same second pass, long-range
      ! coefficient 0.65 rather than 1.0 -- is grid-converged at 1e-8.
      call check(error, worst < 1.0e-5_dp, &
                 "the gradient disagrees with a difference of the energy for "//functional)
   end subroutine gradient_against_energy

   subroutine gradient_at(coords, functional, gradient, err, ok)
      real(dp), intent(in) :: coords(3, 3)
      character(len=*), intent(in) :: functional
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, functional, ctx, err, level=3)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if
      call run_libcint_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      if (.not. err%has_error()) then
         call libcint_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                                   orbital_energies=scf%orbital_energies, &
                                   n_occupied=5, gradient=gradient, error=err, xc=ctx)
      end if
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine gradient_at

   subroutine energy_at(coords, functional, energy, err, ok)
      real(dp), intent(in) :: coords(3, 3)
      character(len=*), intent(in) :: functional
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf

      ok = .false.
      energy = 0.0_dp
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, functional, ctx, err, level=3)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if
      call run_libcint_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      if (.not. err%has_error()) energy = scf%energy
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine energy_at

end module test_mqc_rsh_gradient

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_rsh_gradient, only: collect_mqc_rsh_gradient
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_rsh_gradient", collect_mqc_rsh_gradient)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
