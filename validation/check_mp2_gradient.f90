!! The MP2 analytic gradient, against finite differences of the MP2 energy
program check_mp2_gradient
   !! **Finite difference is the only check that can fail here.** Unlike the
   !! variational gradients, an MP2 gradient has a term -- the orbital
   !! relaxation -- with no counterpart in the energy expression, so agreeing
   !! with a total energy says nothing at all about it. Differentiating this
   !! program's own MP2 energy is what pins it.
   !!
   !! Translational invariance is printed because it is free, and is worth
   !! exactly as much as it is elsewhere: it catches a missing Hellmann-Feynman
   !! term and is blind to almost everything else. The relaxation term is
   !! translationally invariant whether or not it is right.
   !!
   !! The printed values are for comparison against `pyscf.grad.mp2` fed this
   !! repository's own basis JSON.
   !!
   !! Not a ctest case: one converged SCF and one MP2 per displaced geometry.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_mp2
   use mqc_libcint_mp2_gradient, only: libcint_mp2_gradient
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3
   integer, parameter :: BENCH_THREADS = 4
   real(dp), parameter :: STEP = 2.5e-3_dp
      !! Bohr. Chosen by measurement rather than taste: at 5e-3 the HCN case
      !! sat 2.17e-5 from the analytic gradient and at 2.5e-3 it sits 5.42e-6,
      !! a factor of 4.0 for a halved step. That is the central difference's own
      !! O(h^2) error being seen, which is what says the disagreement is the
      !! difference formula and not the derivative.

   integer :: n_bad
   character(len=128) :: filter
   integer :: filter_len, arg_status

   call omp_set_num_threads(min(BENCH_THREADS, omp_get_max_threads()))

   filter = ""
   call get_command_argument(1, filter, length=filter_len, status=arg_status)
   if (arg_status /= 0) filter = ""

   n_bad = 0
   write (*, "(a)") "== MP2 gradients"

   call check_case("H2 / sto-3g", [1, 1], ["H", "H"], &
                   reshape([0.0_dp, 0.0_dp, -0.7_dp, &
                            0.0_dp, 0.0_dp, 0.7_dp], [N_DIM, 2]), &
                   "sto-3g", 2, n_bad)

   call check_case("H2O / sto-3g", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, n_bad)

   ! Asymmetric, so a term that is wrong in a way symmetry cancels shows up.
   call check_case("HCN / sto-3g", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", 14, n_bad)

   ! A general contraction with d functions, where the amplitudes stop being
   ! a handful of numbers and the two-particle density has real structure.
   call check_case("H2O / cc-pvdz", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "cc-pvdz", 10, n_bad)

   call check_case("H2O / 6-31g", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "6-31g", 10, n_bad)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all MP2 gradient checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " case(s)"
      stop 1
   end if

contains

   subroutine check_case(label, numbers, symbols, coords, basis, nelec, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec
      integer, intent(inout) :: n_bad

      real(dp), allocatable :: analytic(:, :), numeric(:, :)
      real(dp) :: translation(3)
      real(dp) :: worst
      integer :: natm, ia, ic
      type(error_t) :: error

      natm = size(numbers)

      if (len_trim(filter) > 0) then
         if (index(label, trim(filter)) == 0) return
      end if

      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label
      flush (output_unit)

      call gradient_at(numbers, symbols, coords, basis, nelec, analytic, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      call numeric_gradient(numbers, symbols, coords, basis, nelec, numeric, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL (finite difference): ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      write (*, "(a)") "  atom        analytic (x, y, z)                      finite difference"
      do ia = 1, natm
         write (*, "(i6,3f14.9,3x,3f14.9)") ia, analytic(:, ia), numeric(:, ia)
      end do

      worst = maxval(abs(analytic - numeric))
      write (*, "(a,es14.4)") "  largest deviation from finite difference: ", worst
      do ic = 1, 3
         translation(ic) = sum(analytic(ic, :))
      end do
      write (*, "(a,es14.4)") "  |sum over atoms| (should be zero):        ", &
         maxval(abs(translation))
      flush (output_unit)

      ! Loose against the analytic gradients' own accuracy, because it is the
      ! difference formula being tolerated rather than the derivative: a
      ! central difference at this step carries its own error of about this
      ! size, which is why the printed columns matter more than the bound.
      if (worst > 1.0e-5_dp) then
         write (*, "(a)") "  FAIL: analytic and numeric disagree"
         n_bad = n_bad + 1
      end if
      if (maxval(abs(translation)) > 1.0e-8_dp) then
         write (*, "(a)") "  FAIL: the gradient does not sum to zero"
         n_bad = n_bad + 1
      end if
   end subroutine check_case

   subroutine gradient_at(numbers, symbols, coords, basis, nelec, gradient, error)
      !! Converge an SCF, then the MP2 gradient over it
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      call run_libcint_rhf(mol, nelec, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, error)
      if (error%has_error()) return
      call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                                gradient, error)
      call mol%destroy()
   end subroutine gradient_at

   subroutine energy_at(numbers, symbols, coords, basis, nelec, energy, error)
      !! One converged MP2 total energy, for the finite difference
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2

      energy = 0.0_dp
      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error)
      if (error%has_error()) return
      call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                           scf%energy, mp2, error)
      if (error%has_error()) return
      energy = mp2%total
      call mol%destroy()
   end subroutine energy_at

   subroutine numeric_gradient(numbers, symbols, coords, basis, nelec, gradient, error)
      !! Central differences of the MP2 total energy
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec
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
            call energy_at(numbers, symbols, shifted, basis, nelec, plus, error)
            if (error%has_error()) return
            shifted(ic, ia) = coords(ic, ia) - STEP
            call energy_at(numbers, symbols, shifted, basis, nelec, minus, error)
            if (error%has_error()) return
            gradient(ic, ia) = (plus - minus)/(2.0_dp*STEP)
         end do
      end do
   end subroutine numeric_gradient

end program check_mp2_gradient
