!! The RI-MP2 analytic gradient, against finite differences of its own energy
program check_ri_mp2_gradient
   !! No external analytic reference exists for this one: PySCF implements no
   !! RI-MP2 gradient at all. What plays that role is
   !! `tools/cpu_validation/ri_mp2_gradient.py`, a numpy reference written
   !! against the same formulation and itself checked against finite
   !! differences -- so the printed components here are meant to be compared
   !! with it component by component, at machine precision rather than at the
   !! 1e-6 a difference formula can offer.
   !!
   !! Finite difference of this program's own RI-MP2 energy is what this file
   !! asserts, and it is the check that cannot be fooled by a misunderstanding
   !! shared between the two implementations.
   !!
   !! Translational invariance is printed and is worth what it is worth
   !! elsewhere: during development of the numpy reference the one-particle
   !! terms were wrong through three revisions while it held at 1e-15.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_ri_mp2
   use mqc_libcint_ri_mp2_gradient, only: libcint_ri_mp2_gradient
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3
   integer, parameter :: BENCH_THREADS = 4
   real(dp), parameter :: STEP = 2.5e-3_dp

   integer :: n_bad
   character(len=128) :: filter
   integer :: filter_len, arg_status

   call omp_set_num_threads(min(BENCH_THREADS, omp_get_max_threads()))
   filter = ""
   call get_command_argument(1, filter, length=filter_len, status=arg_status)
   if (arg_status /= 0) filter = ""

   n_bad = 0
   write (*, "(a)") "== RI-MP2 gradients"

   call check_case("H2 / sto-3g", [1, 1], ["H", "H"], &
                   reshape([0.0_dp, 0.0_dp, -0.7_dp, &
                            0.0_dp, 0.0_dp, 0.7_dp], [N_DIM, 2]), &
                   "sto-3g", "cc-pvdz-rifit", 2, n_bad)

   call check_case("H2O / sto-3g", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", "cc-pvdz-rifit", 10, n_bad)

   ! Asymmetric, so a term wrong in a way symmetry cancels shows up.
   call check_case("HCN / sto-3g", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", "cc-pvdz-rifit", 14, n_bad)

   ! With a core frozen, the finite difference *must* freeze the same core:
   ! the relaxed density's occupied-frozen block comes from the Lagrangian
   ! rather than the amplitudes, and differencing the frozen energy is the one
   ! check that sees it as a derivative rather than a number.
   call check_case("H2O / 6-31g, frozen 1", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "6-31g", "cc-pvtz-rifit", 10, n_bad, n_frozen=1)

   ! More than one frozen orbital, on an asymmetric molecule: two cores on two
   ! different atoms, so an occupied-frozen block indexed off by one lands in
   ! the wrong atom's terms rather than cancelling by symmetry.
   call check_case("HCN / sto-3g, frozen 2", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", "cc-pvdz-rifit", 14, n_bad, n_frozen=2)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all RI-MP2 gradient checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " case(s)"
      stop 1
   end if

contains

   subroutine check_case(label, numbers, symbols, coords, basis, aux_basis, nelec, &
                         n_bad, n_frozen)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec
      integer, intent(inout) :: n_bad
      integer, intent(in), optional :: n_frozen
         !! Core orbitals frozen on *both* sides of the comparison -- the
         !! analytic gradient and the differenced energy have to be
         !! derivatives of the same thing.

      real(dp), allocatable :: analytic(:, :), direct(:, :), numeric(:, :)
      real(dp) :: translation(3)
      real(dp) :: worst, energy
      integer :: natm, ia, ic, frozen
      type(error_t) :: error

      natm = size(numbers)
      frozen = 0
      if (present(n_frozen)) frozen = n_frozen
      if (len_trim(filter) > 0) then
         if (index(label, trim(filter)) == 0) return
      end if

      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label
      flush (output_unit)

      call gradient_at(numbers, symbols, coords, basis, aux_basis, nelec, frozen, &
                       analytic, .false., error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      ! The same gradient with the reference integrals recomputed instead of
      ! stored. Past the in-core limit that is the only path there is, and the
      ! limit is far above anything a test case reaches -- so it is forced here
      ! rather than waited for. The two differ only in where the integrals came
      ! from, so they have to agree to rounding, not to a tolerance.
      call gradient_at(numbers, symbols, coords, basis, aux_basis, nelec, frozen, &
                       direct, .true., error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL (direct): ", error%get_message()
         n_bad = n_bad + 1
         return
      end if
      write (*, "(a,es14.4)") "  stored against recomputed integrals:       ", &
         maxval(abs(analytic - direct))
      if (maxval(abs(analytic - direct)) > 1.0e-10_dp) then
         write (*, "(a)") "  FAIL: the two integral paths disagree"
         n_bad = n_bad + 1
      end if

      call numeric_gradient(numbers, symbols, coords, basis, aux_basis, nelec, &
                            frozen, numeric, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL (finite difference): ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      ! Printed so the components can be read against
      ! `tools/cpu_validation/ri_mp2_gradient.py`, which is the only
      ! machine-precision reference there is. The energy goes with them: the
      ! two codes take their basis data from different places, and a gradient
      ! of an energy that differs in the ninth decimal differs in the eighth.
      call energy_at(numbers, symbols, coords, basis, aux_basis, nelec, frozen, &
                     energy, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL (energy): ", error%get_message()
         n_bad = n_bad + 1
         return
      end if
      write (*, "(a,f18.12)") "  E(RI-MP2) = ", energy

      write (*, "(a)") "  atom        analytic (x, y, z)                            finite difference"
      do ia = 1, natm
         write (*, "(i6,3f16.11,3x,3f16.11)") ia, analytic(:, ia), numeric(:, ia)
      end do

      worst = maxval(abs(analytic - numeric))
      write (*, "(a,es14.4)") "  largest deviation from finite difference: ", worst
      do ic = 1, 3
         translation(ic) = sum(analytic(ic, :))
      end do
      write (*, "(a,es14.4)") "  |sum over atoms| (should be zero):        ", &
         maxval(abs(translation))
      flush (output_unit)

      if (worst > 1.0e-5_dp) then
         write (*, "(a)") "  FAIL: analytic and numeric disagree"
         n_bad = n_bad + 1
      end if
      if (maxval(abs(translation)) > 1.0e-8_dp) then
         write (*, "(a)") "  FAIL: the gradient does not sum to zero"
         n_bad = n_bad + 1
      end if
   end subroutine check_case

   subroutine gradient_at(numbers, symbols, coords, basis, aux_basis, nelec, &
                          frozen, gradient, force_direct, error)
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec, frozen
      real(dp), allocatable, intent(out) :: gradient(:, :)
      logical, intent(in) :: force_direct
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf

      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      call build_libcint_molecule(numbers, symbols, coords, aux_basis, aux, error)
      if (error%has_error()) return
      call run_libcint_rhf(mol, nelec, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, error)
      if (error%has_error()) return
      call libcint_ri_mp2_gradient(mol, aux, scf%orbitals, scf%orbital_energies, &
                                   nelec/2, gradient, error, n_frozen=frozen, &
                                   force_direct=force_direct)
      call mol%destroy()
      call aux%destroy()
   end subroutine gradient_at

   subroutine energy_at(numbers, symbols, coords, basis, aux_basis, nelec, frozen, &
                        energy, error)
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec, frozen
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2

      energy = 0.0_dp
      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      call build_libcint_molecule(numbers, symbols, coords, aux_basis, aux, error)
      if (error%has_error()) return
      call run_libcint_rhf(mol, nelec, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, error)
      if (error%has_error()) return
      call run_libcint_ri_mp2(mol, aux, scf%orbitals, scf%orbital_energies, nelec/2, &
                              scf%energy, mp2, error, n_frozen=frozen)
      if (error%has_error()) return
      energy = mp2%total
      call mol%destroy()
      call aux%destroy()
   end subroutine energy_at

   subroutine numeric_gradient(numbers, symbols, coords, basis, aux_basis, nelec, &
                               frozen, gradient, error)
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec, frozen
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: shifted(:, :)
      real(dp) :: plus, minus
      integer :: natm, ia, ic

      natm = size(numbers)
      allocate (gradient(3, natm), shifted(3, natm))
      gradient = 0.0_dp

      do ia = 1, natm
         do ic = 1, 3
            shifted = coords
            shifted(ic, ia) = coords(ic, ia) + STEP
            call energy_at(numbers, symbols, shifted, basis, aux_basis, nelec, &
                           frozen, plus, error)
            if (error%has_error()) return
            shifted(ic, ia) = coords(ic, ia) - STEP
            call energy_at(numbers, symbols, shifted, basis, aux_basis, nelec, &
                           frozen, minus, error)
            if (error%has_error()) return
            gradient(ic, ia) = (plus - minus)/(2.0_dp*STEP)
         end do
      end do
   end subroutine numeric_gradient

end program check_ri_mp2_gradient
