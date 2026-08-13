module test_mqc_libcint_rsh
   !! Pins the range-separation parameter to the kernel it is supposed to change.
   !!
   !! The validation manifest checks that ωB97X and CAM-B3LYP land on PySCF's
   !! energies, which is the answer. This checks the mechanism, because the way it
   !! failed was invisible from the answer's side: the parameter is passed to
   !! libcint by writing a slot of the `env` array, the slot constants are libcint's
   !! own 0-based offsets, and writing the neighbouring one instead put omega into
   !! `PTR_RINV_ZETA` -- which a plain two-electron integral ignores. So the
   !! long-range pass returned *full-range* exchange, the short- and long-range
   !! coefficients summed back to unscaled exchange, and the SCF converged happily
   !! several Hartree out of place. Nothing in that run looked wrong.
   !!
   !! What is checked here is the two limits of the erf kernel, which no wrong slot
   !! can satisfy at once:
   !!
   !!   * as omega goes to zero, erf(omega r)/r vanishes and so must the long-range
   !!     exchange. This is the limit the dead-slot bug fails: an ignored parameter
   !!     leaves the full-range matrix, which does not vanish;
   !!   * as omega grows, erf(omega r)/r becomes 1/r and the long-range exchange
   !!     approaches the full-range matrix as omega^-2. That rate, and not merely
   !!     the approach, is what is asserted -- it is the limit a *mis-scaled* omega
   !!     fails, and the pair together leave no room for one that is only close;
   !!   * at the physical omega of an actual functional the two must differ
   !!     substantially, which is what makes the second pass worth doing at all;
   !!   * and the magnitude has to be monotone in between, so that omega attenuates
   !!     rather than merely perturbs.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, direct_stats_t
   use mqc_libcint_rhf, only: SCF_GUESS_SAD
   use mqc_libcint_atomic_guess, only: build_atomic_guess, clear_atomic_cache
   implicit none
   private
   public :: collect_mqc_libcint_rsh_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_libcint_rsh_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("small_omega_leaves_no_long_range_exchange", test_small_omega), &
                  new_unittest("large_omega_recovers_full_range_exchange", test_large_omega), &
                  new_unittest("physical_omega_changes_the_exchange", test_physical_omega), &
                  new_unittest("attenuation_is_monotone_in_omega", test_monotone) &
                  ]
   end subroutine collect_mqc_libcint_rsh_tests

   subroutine exchange_only(omega, k, peak, error)
      !! The long-range exchange matrix alone, at one omega
      !!
      !! A zero core Hamiltonian and no Coulomb term, which is exactly how the
      !! Kohn-Sham build asks for this: what comes back is the exchange
      !! contribution and nothing that would have to be subtracted out again. The
      !! density is the SAD guess rather than a converged one -- any symmetric
      !! density exercises the same quartet loop, and not converging an SCF keeps
      !! this a fast test.
      real(dp), intent(in), optional :: omega
      real(dp), allocatable, intent(out) :: k(:, :)
      real(dp), intent(out) :: peak
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(direct_stats_t) :: stats
      type(error_t) :: err
      real(dp), allocatable :: bounds(:, :), d_a(:, :), d_b(:, :), density(:, :)
      real(dp), allocatable :: h_zero(:, :)
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call clear_atomic_cache()
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "water must build")
      if (allocated(error)) return

      call build_atomic_guess(mol, SCF_GUESS_SAD, d_a, d_b, err)
      call check(error,.not. err%has_error(), "the SAD guess must build")
      if (allocated(error)) return

      call schwarz_bounds(mol, bounds, err)
      call check(error,.not. err%has_error(), "the Schwarz bounds must build")
      if (allocated(error)) return

      allocate (density(mol%nao, mol%nao), h_zero(mol%nao, mol%nao), k(mol%nao, mol%nao))
      density = d_a + d_b
      h_zero = 0.0_dp

      call build_fock_direct(mol, h_zero, density, bounds, k, stats, err, &
                             k_scale=1.0_dp, j_scale=0.0_dp, omega=omega)
      call check(error,.not. err%has_error(), "the exchange build must succeed")
      if (allocated(error)) return

      peak = maxval(abs(k))
      call mol%destroy()
   end subroutine exchange_only

   subroutine test_small_omega(error)
      !! omega -> 0: no long-range exchange left
      !!
      !! The test the wrong `env` slot fails. erf(omega r)/r tends to the constant
      !! 2 omega / sqrt(pi), so the matrix is linear in omega near zero and at
      !! 1e-6 must be six orders below the full-range one. An omega that never
      !! arrived would leave the full-range matrix here, off by a factor of 1e6.
      !!
      !! Note that omega is not passed as exactly zero, which libcint reads as "no
      !! attenuation" and would make this test assert the opposite of what it
      !! means.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: k_small(:, :), k_full(:, :)
      real(dp) :: peak_small, peak_full

      call exchange_only(1.0e-6_dp, k_small, peak_small, error)
      if (allocated(error)) return
      call exchange_only(k=k_full, peak=peak_full, error=error)
      if (allocated(error)) return

      call check(error, peak_full > 0.1_dp, &
                 "the full-range exchange must be of order one, or this test is "// &
                 "comparing two small numbers and asserting nothing")
      if (allocated(error)) return
      call check(error, peak_small < 1.0e-5_dp*peak_full, &
                 "a vanishing range parameter must leave essentially no long-range "// &
                 "exchange; a full-range matrix here means omega never reached libcint")
   end subroutine test_small_omega

   subroutine test_large_omega(error)
      !! omega -> infinity: the full Coulomb kernel, at the right rate
      !!
      !! What is left over at finite omega is the short-range complement,
      !! erfc(omega r)/r, whose integral against a smooth density is pi/omega^2 --
      !! so the long-range matrix approaches the full-range one from below as
      !! omega^-2, and measured here it does exactly that: 5.79e-3, 5.81e-5,
      !! 5.81e-7 at omega = 1e2, 1e3, 1e4.
      !!
      !! Asserting the power law rather than one threshold is what makes this a
      !! test of omega's *scale* and not only of its arrival. A parameter that
      !! reached libcint multiplied by some constant c would still fall off as
      !! omega^-2 -- but shifted by c^2, which the absolute band below catches;
      !! while a residual that is not the short-range complement at all would
      !! miss the factor of a hundred per decade.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: k_full(:, :), k_hundred(:, :), k_thousand(:, :)
      real(dp) :: peak, d_hundred, d_thousand

      call exchange_only(k=k_full, peak=peak, error=error)
      if (allocated(error)) return
      call exchange_only(1.0e2_dp, k_hundred, peak, error)
      if (allocated(error)) return
      call exchange_only(1.0e3_dp, k_thousand, peak, error)
      if (allocated(error)) return

      d_hundred = maxval(abs(k_hundred - k_full))
      d_thousand = maxval(abs(k_thousand - k_full))

      ! The band: a factor of two either way on omega moves this by four.
      call check(error, d_thousand > 1.0e-5_dp .and. d_thousand < 3.0e-4_dp, &
                 "at omega = 1e3 the long-range exchange must sit ~6e-5 below the "// &
                 "full-range one; outside that band omega arrives mis-scaled")
      if (allocated(error)) return
      ! And the rate, which says the remainder really is the short-range kernel.
      call check(error, d_hundred/d_thousand > 80.0_dp .and. &
                 d_hundred/d_thousand < 120.0_dp, &
                 "the residual must fall as omega^-2, a hundredfold per decade")
   end subroutine test_large_omega

   subroutine test_physical_omega(error)
      !! At ωB97X's omega the two kernels really are different
      !!
      !! Without this the pair of limits above could both pass on a build that
      !! attenuates far too sharply or far too weakly to matter at the omega any
      !! real functional uses -- and the whole point of the second pass is that at
      !! 0.3 the long-range matrix is a substantially different object.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: k_mid(:, :), k_full(:, :)
      real(dp) :: peak_mid, peak_full

      call exchange_only(0.3_dp, k_mid, peak_mid, error)
      if (allocated(error)) return
      call exchange_only(k=k_full, peak=peak_full, error=error)
      if (allocated(error)) return

      call check(error, maxval(abs(k_mid - k_full)) > 0.01_dp*peak_full, &
                 "at omega = 0.3 the long-range exchange must differ substantially "// &
                 "from the full-range one")
      if (allocated(error)) return
      call check(error, peak_mid < peak_full, &
                 "attenuation cannot make the exchange larger")
   end subroutine test_physical_omega

   subroutine test_monotone(error)
      !! Larger omega, more exchange -- all the way up
      !!
      !! An erf-attenuated kernel is a fraction of the full one that grows with
      !! omega, so the magnitudes have to be ordered. Cheap, and it rules out a
      !! parameter that arrives with the wrong sign: negative omega is libcint's
      !! *short-range* complement, which runs the other way.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: k(:, :)
      real(dp) :: peak(4)
      integer :: i
      real(dp), parameter :: OMEGAS(4) = [0.05_dp, 0.3_dp, 1.0_dp, 5.0_dp]

      do i = 1, 4
         call exchange_only(OMEGAS(i), k, peak(i), error)
         if (allocated(error)) return
      end do

      do i = 2, 4
         call check(error, peak(i) > peak(i - 1), &
                    "long-range exchange must grow with the range parameter")
         if (allocated(error)) return
      end do
   end subroutine test_monotone

end module test_mqc_libcint_rsh

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_rsh, only: collect_mqc_libcint_rsh_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_rsh", collect_mqc_libcint_rsh_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
