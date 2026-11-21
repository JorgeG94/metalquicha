!! Natural orbitals of a density, held to what a density must satisfy
module test_mqc_czt_natural_orbitals
   !! The occupations themselves are compared with PySCF end to end, on
   !! unrestricted references and on RI-MP2 relaxed densities, and agree to
   !! 1e-10. What is pinned here is what needs no other program:
   !!
   !!   * a closed-shell determinant has occupations of exactly two and zero,
   !!     and both unpaired-electron counts vanish;
   !!   * the natural orbitals are orthonormal under the overlap and rebuild
   !!     the density they came from, `D = C n C^T`, which fails for a routine
   !!     that returns the right eigenvalues with vectors in the wrong basis;
   !!   * a high-spin determinant has one orbital of occupation one per
   !!     unpaired electron, and the linear count says so exactly;
   !!   * an unrestricted singlet that has not broken symmetry is a closed
   !!     shell, whatever driver produced it;
   !!   * the summary is arithmetic, and is checked on occupations written
   !!     down by hand, including ones a relaxed density leaves a hair outside
   !!     [0, 2].
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf
   use mqc_czt_natural_orbitals, only: natural_occupation_summary_t, &
                                       density_natural_orbitals, &
                                       summarize_natural_occupations
   implicit none
   private

   public :: collect_mqc_czt_natural_orbitals_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_czt_natural_orbitals_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("closed_shell_occupations_are_two_and_zero", &
                               test_closed_shell), &
                  new_unittest("natural_orbitals_rebuild_their_density", test_rebuild), &
                  new_unittest("high_spin_has_one_electron_per_open_orbital", &
                               test_high_spin), &
                  new_unittest("summary_counts_what_it_says", test_summary), &
                  new_unittest("mismatched_shapes_are_refused", test_shapes) &
                  ]
   end subroutine collect_mqc_czt_natural_orbitals_tests

   subroutine water(mol, err)
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
   end subroutine water

   subroutine test_closed_shell(error)
      !! Water, restricted and as an unrestricted singlet: five twos, the rest zero
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: r, u
      type(error_t) :: err
      type(natural_occupation_summary_t) :: summary
      real(dp), allocatable :: s(:, :), natural(:, :), occ(:)

      call water(mol, err)
      call check(error,.not. err%has_error(), "water must build")
      if (allocated(error)) return
      call mol%overlap(s)

      call run_czt_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., r, err, in_core=.true.)
      call density_natural_orbitals(r%density, s, natural, occ, err)
      call check(error,.not. err%has_error(), "natural orbitals: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(occ(1:5) - 2.0_dp)) < 1.0e-10_dp, &
                 "the five occupied orbitals must hold two electrons each")
      if (allocated(error)) return
      call check(error, maxval(abs(occ(6:))) < 1.0e-10_dp, "the virtuals must be empty")
      if (allocated(error)) return
      call summarize_natural_occupations(occ, summary)
      call check(error, summary%unpaired_linear < 1.0e-9_dp .and. &
                 summary%unpaired_quadratic < 1.0e-9_dp .and. summary%n_fractional == 0, &
                 "a closed-shell determinant has no unpaired electrons")
      if (allocated(error)) return

      ! The same molecule through the unrestricted driver: the two spin
      ! densities are equal, so their sum is the restricted density.
      call run_czt_uhf(mol, 10, 1, 200, 1.0e-12_dp, 1.0e-10_dp, .false., u, err, in_core=.true.)
      call density_natural_orbitals(u%density + u%density_beta, s, natural, occ, err)
      call check(error,.not. err%has_error(), "natural orbitals: "//err%get_full_trace())
      if (allocated(error)) return
      call summarize_natural_occupations(occ, summary)
      call check(error, summary%unpaired_linear < 1.0e-7_dp, &
                 "an unrestricted singlet that kept its symmetry is a closed shell")
      call mol%destroy()
   end subroutine test_closed_shell

   subroutine test_rebuild(error)
      !! Orthonormal under S, and `C n C^T` is the density that went in
      !!
      !! On the cation, where the density is not idempotent and so the
      !! eigenvectors are not determined up to a rotation within a degenerate
      !! block of twos.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: u
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), natural(:, :), occ(:), density(:, :)
      real(dp), allocatable :: metric(:, :), rebuilt(:, :)
      integer :: i, n

      call water(mol, err)
      call check(error,.not. err%has_error(), "water must build")
      if (allocated(error)) return
      call mol%overlap(s)
      call run_czt_uhf(mol, 9, 2, 300, 1.0e-12_dp, 1.0e-10_dp, .false., u, err, in_core=.true.)
      call check(error, u%converged, "the cation must converge")
      if (allocated(error)) return
      density = u%density + u%density_beta
      call density_natural_orbitals(density, s, natural, occ, err)
      call check(error,.not. err%has_error(), "natural orbitals: "//err%get_full_trace())
      if (allocated(error)) return

      n = size(occ)
      metric = matmul(transpose(natural), matmul(s, natural))
      do i = 1, n
         metric(i, i) = metric(i, i) - 1.0_dp
      end do
      call check(error, maxval(abs(metric)) < 1.0e-10_dp, &
                 "the natural orbitals are not orthonormal under the overlap")
      if (allocated(error)) return

      allocate (rebuilt(size(density, 1), size(density, 2)))
      rebuilt = 0.0_dp
      do i = 1, n
         rebuilt = rebuilt + occ(i)*matmul(natural(:, i:i), transpose(natural(:, i:i)))
      end do
      call check(error, maxval(abs(rebuilt - density)) < 1.0e-10_dp, &
                 "the natural orbitals and occupations do not rebuild the density")
      if (allocated(error)) return
      call check(error, all(occ(1:n - 1) >= occ(2:n)), "the occupations must descend")
      if (allocated(error)) return
      call check(error, abs(sum(occ) - 9.0_dp) < 1.0e-9_dp, &
                 "the occupations must sum to the electron count")
      call mol%destroy()
   end subroutine test_rebuild

   subroutine test_high_spin(error)
      !! Triplet H2: both electrons unpaired, two occupations of exactly one
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: u
      type(error_t) :: err
      type(natural_occupation_summary_t) :: summary
      real(dp), allocatable :: s(:, :), natural(:, :), occ(:)
      real(dp) :: c(3, 2)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.7414_dp*ANG], [3, 2])
      call build_czt_molecule([1, 1], ["H ", "H "], c, "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "H2 must build")
      if (allocated(error)) return
      call mol%overlap(s)
      call run_czt_uhf(mol, 2, 3, 100, 1.0e-12_dp, 1.0e-10_dp, .false., u, err)
      call check(error, u%converged, "triplet H2 must converge")
      if (allocated(error)) return
      call density_natural_orbitals(u%density + u%density_beta, s, natural, occ, err)
      call check(error,.not. err%has_error(), "natural orbitals: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(occ(1:2) - 1.0_dp)) < 1.0e-10_dp, &
                 "the two open orbitals must hold one electron each")
      if (allocated(error)) return
      call summarize_natural_occupations(occ, summary)
      call check(error, abs(summary%unpaired_linear - 2.0_dp) < 1.0e-9_dp, &
                 "the linear count must be the number of unpaired electrons")
      if (allocated(error)) return
      call check(error, abs(summary%unpaired_quadratic - 2.0_dp) < 1.0e-9_dp, &
                 "the quadratic count must agree at occupations of one")
      if (allocated(error)) return
      call check(error, summary%n_fractional, 2)
      call mol%destroy()
   end subroutine test_high_spin

   subroutine test_summary(error)
      !! The counts, on occupations written by hand
      type(error_type), allocatable, intent(out) :: error

      type(natural_occupation_summary_t) :: summary
      real(dp) :: occ(6)

      ! Two outside [0, 2] by rounding or by a relaxed density's licence.
      occ = [2.0000004_dp, 1.9_dp, 1.0_dp, 0.1_dp, 0.015_dp, -0.0000004_dp]
      call summarize_natural_occupations(occ, summary)
      call check(error, abs(summary%n_electrons - sum(occ)) < 1.0e-14_dp, &
                 "the electron count is the plain sum, clamping nothing")
      if (allocated(error)) return
      ! min(n, 2-n): 0, 0.1, 1, 0.1, 0.015, 0 with the two strays clamped.
      call check(error, abs(summary%unpaired_linear - 1.215_dp) < 1.0e-12_dp, &
                 "the linear unpaired count is wrong")
      if (allocated(error)) return
      call check(error, abs(summary%unpaired_quadratic - &
                            ((1.9_dp*0.1_dp)**2*2.0_dp + 1.0_dp &
                             + (0.015_dp*1.985_dp)**2)) < 1.0e-12_dp, &
                 "the quadratic unpaired count is wrong")
      if (allocated(error)) return
      ! 1.9, 1.0 and 0.1 are inside the window; 0.015 is not.
      call check(error, summary%n_fractional, 3)
      if (allocated(error)) return
      call check(error, abs(summary%largest_deviation - 1.0_dp) < 1.0e-14_dp)
   end subroutine test_summary

   subroutine test_shapes(error)
      !! A density and an overlap of different sizes are not silently truncated
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: d(3, 3), s(4, 4)
      real(dp), allocatable :: natural(:, :), occ(:)

      d = 0.0_dp
      s = 0.0_dp
      call density_natural_orbitals(d, s, natural, occ, err)
      call check(error, err%has_error(), "mismatched shapes must be refused")
   end subroutine test_shapes

end module test_mqc_czt_natural_orbitals

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_natural_orbitals, only: collect_mqc_czt_natural_orbitals_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_natural_orbitals", &
                               collect_mqc_czt_natural_orbitals_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
