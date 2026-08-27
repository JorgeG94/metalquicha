!! Range-separated and non-local gradients, against a difference of the energy
module test_mqc_rsh_gradient
   !! Two families of functional need more than the standard gradient, and
   !! both are checked here the same way.
   !!
   !! A range-separated functional splits its exchange over an error-function
   !! kernel, so its gradient needs two exchange derivatives rather than one:
   !! the full-range pass every hybrid takes, and a second at the screened
   !! omega. Only the density-fitted path built the second one for a while;
   !! the exact-ERI path refused outright rather than return a gradient
   !! missing wB97X's dominant exchange term.
   !!
   !! A functional carrying VV10 has a non-local correlation term whose energy
   !! is a double integral over the grid. Its gradient turns out not to need
   !! new machinery -- the pair sum is spent producing `vrho` and `vsigma`, and
   !! what follows is the ordinary GGA contraction -- but it does need two
   !! terms no semilocal functional has: the kernel depends on where the grid
   !! points are, and a point's weight enters the energy twice rather than
   !! once.
   !!
   !! Differencing the total energy is the right check for both.
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
   use mqc_libcint_hessian, only: ks_hessian
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
                  new_unittest("global-hybrid-unchanged", test_b3lyp_unchanged), &
                  new_unittest("rsh-hessian", test_rsh_hessian), &
                  new_unittest("vv10-gradient", test_vv10) &
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

   subroutine test_rsh_hessian(error)
      !! The assembled Hessian, for both range-separated functionals at once
      type(error_type), allocatable, intent(out) :: error
      call hess_against_gradient("cam-b3lyp", error)
      if (allocated(error)) return
      call hess_against_gradient("wb97x", error)
   end subroutine test_rsh_hessian

   subroutine hess_against_gradient(functional, error)
      !! `ks_hessian` against a difference of `libcint_scf_gradient`
      !!
      !! The gradient this differences is the one checked against the energy
      !! above, so a failure here is in the second derivative rather than
      !! somewhere below it. Two bugs of exactly one kind were found this way,
      !! both in code that built an attenuated `env` and then handed the
      !! integral the unattenuated one: `h1_contract` passed `tab%env`, and
      !! `build_fock_direct_many` passed `tab%env`. Neither is detectable from
      !! the result's shape -- the Hessian stays symmetric and translationally
      !! invariant, because full-range integrals are perfectly good integrals
      !! of the wrong operator -- which is why this differences a validated
      !! gradient rather than checking invariants.
      character(len=*), intent(in) :: functional
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-3_dp
      real(dp) :: shifted(3, 3), fd, worst, wtrans, rowsum
      real(dp), allocatable :: hess(:, :, :, :), plus(:, :), minus(:, :)
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, ja, b

      if (.not. xc_available()) return

      call hess_at(WATER, functional, hess, err, ok)
      call check(error, ok, "the Kohn-Sham Hessian failed for "//functional)
      if (allocated(error)) return

      wtrans = 0.0_dp
      do a = 1, 3
         do b = 1, 3
            do ia = 1, 3
               rowsum = sum(hess(a, b, ia, :))
               wtrans = max(wtrans, abs(rowsum))
            end do
         end do
      end do

      worst = 0.0_dp
      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            call gradient_at(shifted, functional, plus, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced gradient failed")
               return
            end if
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            call gradient_at(shifted, functional, minus, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced gradient failed")
               return
            end if
            do ja = 1, 3
               do b = 1, 3
                  fd = (plus(b, ja) - minus(b, ja))/(2.0_dp*H)
                  worst = max(worst, abs(hess(a, b, ia, ja) - fd))
               end do
            end do
         end do
      end do

      ! Loose for the reason the neighbouring Hessian tests are loose: the
      ! omitted grid response breaks this by an amount that depends on the
      ! functional and the grid, and it is the same omission PySCF makes by
      ! default. What matters is that a range-separated functional is now no
      ! worse than a global hybrid -- cam-B3LYP lands at 4.5e-4 against B3LYP's
      ! 4.5e-4, and wB97X at 5.8e-3 against its own 5.9e-3 translational
      ! invariance, so in both cases the residual *is* the grid response.
      !
      ! Before the two `env` fixes these read 0.53 and 2.15, so the tolerance
      ! separates right from wrong by two decades. The tight comparison is
      ! against PySCF's analytic Hessian, where cam-B3LYP agrees to 1.7e-8 and
      ! wB97X to 1.3e-8 -- the same 1e-8 every non-separated functional gets.
      ! That one lives in the validation suite because it needs PySCF.
      call check(error, wtrans < 2.0e-2_dp, &
                 "the Kohn-Sham Hessian is not translationally invariant for "//functional)
      if (allocated(error)) return
      call check(error, worst < 2.0e-2_dp, &
                 "the Kohn-Sham Hessian disagrees with differences of the gradient for "// &
                 functional)
   end subroutine hess_against_gradient

   subroutine hess_at(coords, functional, hess, err, ok)
      real(dp), intent(in) :: coords(3, 3)
      character(len=*), intent(in) :: functional
      real(dp), allocatable, intent(out) :: hess(:, :, :, :)
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
         call ks_hessian(mol, WATER_Z, scf%density, scf%orbitals, scf%orbital_energies, &
                         5, ctx, ctx%exx_fraction, hess, err)
      end if
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine hess_at

   subroutine test_vv10(error)
      !! The non-local term's gradient, alone and on top of everything else
      !!
      !! `b97m-v` is VV10 over a meta-GGA and nothing else; `wb97x-v` adds
      !! range separation; `wb97m-v` has all three at once. All three land on
      !! the same 1.47e-7 as a plain hybrid, which is the point -- the
      !! non-local term is not approximated here, it is differentiated.
      type(error_type), allocatable, intent(out) :: error
      call gradient_against_energy("b97m-v", error)
      if (allocated(error)) return
      call gradient_against_energy("wb97x-v", error)
      if (allocated(error)) return
      call gradient_against_energy("wb97m-v", error)
   end subroutine test_vv10

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

      ! Bounded by the difference, not by the gradient. Every functional here
      ! lands between 1.46e-7 and 1.48e-7 -- that sameness is the point: it is
      ! the h^2 step error, so the range-separated and non-local gradients are
      ! as good as the global hybrid's rather than merely close enough.
      !
      ! Dropping the long-range exchange pass misses by 1e-1 on wB97X and a
      ! coefficient confused between alpha and alpha+beta by 1e-2. Omitting
      ! VV10's explicit coordinate derivative, or taking its weight derivative
      ! as `rho*exc` instead of `dE/dw`, misses by 4.5e-4 -- and that one does
      ! not shrink when the grid is refined, because it is not a quadrature
      ! error. The tolerance separates right from wrong by decades in each
      ! case.
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
