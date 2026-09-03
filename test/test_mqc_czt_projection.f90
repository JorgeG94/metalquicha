!! Projecting a density between basis sets, against what it must satisfy
module test_mqc_czt_projection
   !! The basis-set-projection guess solves a cheap SCF in a small basis and
   !! carries its density into the target basis. What comes out is a starting
   !! point, so there is no reference value to compare against -- but "a
   !! starting point" is not the same as "any matrix", and the properties it
   !! owes are exactly the ones that make it usable:
   !!
   !!   * **Idempotency against the target overlap**: `D S D = 2 D` for a
   !!     closed-shell density in the target basis. A projection that lost the
   !!     orthogonalisation produces something density-shaped that the first SCF
   !!     iteration then has to undo.
   !!   * **The electron count survives the crossing**: `Tr(D S) = 2 n_occ`
   !!     exactly. This is the one a wrong `S_BA` breaks, and it breaks it by an
   !!     amount that looks like a converging SCF rather than a bug.
   !!   * **Projecting a basis into itself changes nothing.** With target and
   !!     source the same basis, the projected density must equal the density
   !!     built directly from those orbitals. That is the identity that fixes
   !!     the orientation of the cross overlap -- transposed, it still produces
   !!     a symmetric, plausible matrix.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: run_czt_rhf, rhf_result_t
   use mqc_czt_projection, only: cross_overlap, project_occupied
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_projection_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_czt_projection_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("cross_overlap_of_a_basis_with_itself", test_self_overlap), &
                  new_unittest("projection_keeps_the_electron_count", test_electron_count), &
                  new_unittest("projected_density_is_idempotent", test_idempotent), &
                  new_unittest("projecting_into_the_same_basis_is_identity", test_identity) &
                  ]
   end subroutine collect_mqc_czt_projection_tests

   subroutine water(basis, mol, err)
      character(len=*), intent(in) :: basis
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule(z, symbols, c, basis, mol, err)
   end subroutine water

   subroutine small_scf(mol, scf, err)
      !! A converged STO-3G water, which is what the ladder starts from.
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      call water("sto-3g", mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 10, 60, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
   end subroutine small_scf

   subroutine test_self_overlap(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: s_target(:, :), s_cross(:, :), plain(:, :)

      call water("sto-3g", mol, err)
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      ! Target and source the same basis: both blocks of the merged overlap are
      ! then that basis's own overlap, which the molecule can produce directly.
      call cross_overlap(mol, mol, s_target, s_cross, err)
      call check(error,.not. err%has_error(), "the cross overlap")
      if (allocated(error)) return

      call mol%overlap(plain)
      call check(error, maxval(abs(s_target - plain)) < 1.0e-12_dp, &
                 "S_BB from the merged molecule is not the target basis overlap")
      if (allocated(error)) return
      call check(error, maxval(abs(s_cross - plain)) < 1.0e-12_dp, &
                 "S_BA between a basis and itself is not that basis's overlap")
   end subroutine test_self_overlap

   subroutine test_electron_count(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: small, target_mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), s_target(:, :)
      real(dp) :: trace
      integer :: i, j

      call small_scf(small, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the small-basis SCF")
      if (allocated(error)) return
      call water("6-31g", target_mol, err)
      call check(error,.not. err%has_error(), "the target basis")
      if (allocated(error)) return

      call project_occupied(target_mol, small, scf%orbitals, scf%n_occupied, density, err)
      call check(error,.not. err%has_error(), "projecting")
      if (allocated(error)) return
      call check(error, size(density, 1) == target_mol%nao, &
                 "the projected density is not in the target basis")
      if (allocated(error)) return

      call target_mol%overlap(s_target)
      trace = 0.0_dp
      do j = 1, target_mol%nao
         do i = 1, target_mol%nao
            trace = trace + density(i, j)*s_target(j, i)
         end do
      end do
      call check(error, abs(trace - 2.0_dp*scf%n_occupied) < 1.0e-9_dp, &
                 "the projected density does not hold the electrons it started with")
   end subroutine test_electron_count

   subroutine test_idempotent(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: small, target_mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), s_target(:, :), dsd(:, :)

      call small_scf(small, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the small-basis SCF")
      if (allocated(error)) return
      call water("6-31g", target_mol, err)
      call check(error,.not. err%has_error(), "the target basis")
      if (allocated(error)) return

      call project_occupied(target_mol, small, scf%orbitals, scf%n_occupied, density, err)
      call check(error,.not. err%has_error(), "projecting")
      if (allocated(error)) return

      call target_mol%overlap(s_target)
      dsd = matmul(density, matmul(s_target, density))
      call check(error, maxval(abs(dsd - 2.0_dp*density)) < 1.0e-8_dp, &
                 "D S D /= 2 D: the projected density is not a projector")
   end subroutine test_idempotent

   subroutine test_identity(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: small
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: projected(:, :), direct(:, :)

      call small_scf(small, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the small-basis SCF")
      if (allocated(error)) return

      ! Projecting a basis into itself is the identity, so this compares the
      ! projection against the density the SCF already built.
      call project_occupied(small, small, scf%orbitals, scf%n_occupied, projected, err)
      call check(error,.not. err%has_error(), "projecting into the same basis")
      if (allocated(error)) return

      direct = 2.0_dp*matmul(scf%orbitals(:, 1:scf%n_occupied), &
                             transpose(scf%orbitals(:, 1:scf%n_occupied)))
      call check(error, maxval(abs(projected - direct)) < 1.0e-9_dp, &
                 "projecting a basis into itself changed the density")
   end subroutine test_identity

end module test_mqc_czt_projection

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_projection, only: collect_mqc_czt_projection_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_projection", collect_mqc_czt_projection_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
