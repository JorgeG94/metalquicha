!! Quasi-atomic orbitals and the population-bond-order matrix
module test_mqc_quao
   !! Orbitals spanning the same space as the occupied valence plus
   !! valence-virtual orbitals, but each belonging to one atom -- deformed
   !! free-atom orbitals rather than delocalized molecular ones.
   !!
   !! The properties worth asserting are the ones that hold whatever the
   !! molecule, because the published numbers to compare against are all for
   !! *oriented* quasi-atomic orbitals, which is the next step. What is checkable
   !! here is that the orbitals are orthonormal and span the right space, that
   !! each really is atomic, that the density in their basis accounts for every
   !! valence electron, and that water comes out with the bonding pattern a
   !! first-year course would draw.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t
   use mqc_aambs, only: aambs_dimensions, aambs_dimensions_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule, &
                                    mixed_basis_overlap
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_quao, only: build_aambs_molecule, mo_aambs_overlap, &
                               valence_virtual_orbitals, vvo_result_t, &
                               aambs_atom_ranges, quasi_atomic_orbitals, quao_result_t
   implicit none
   private

   public :: collect_mqc_quao_tests

   integer, parameter :: N_DIM = 3
   real(dp), parameter :: WATER(N_DIM, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

contains

   subroutine collect_mqc_quao_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("orbitals_are_orthonormal", test_orthonormal), &
                  new_unittest("population_counts_the_electrons", test_populations), &
                  new_unittest("orbitals_are_actually_atomic", test_atomic_character), &
                  new_unittest("water_bonding_pattern", test_water_pattern) &
                  ]
   end subroutine collect_mqc_quao_tests

   subroutine water_quaos(basis_name, quao, overlap, dims, err, ok)
      !! Everything from an SCF through to the quasi-atomic orbitals
      character(len=*), intent(in) :: basis_name
      type(quao_result_t), intent(out) :: quao
      real(dp), allocatable, intent(out) :: overlap(:, :)
      type(aambs_dimensions_t), intent(out) :: dims
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol, aambs
      type(rhf_result_t) :: scf
      type(vvo_result_t) :: vvo
      real(dp), allocatable :: mixed(:, :), s_mbs(:, :), projection(:, :)
      real(dp), allocatable :: valence_internal(:, :)
      integer, allocatable :: core_off(:), core_n(:), val_off(:), val_n(:)
      integer :: n_ao, n_valocc

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, basis_name, mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (err%has_error() .or. .not. scf%converged) return

      call build_aambs_molecule(WATER_Z, WATER_SYM, WATER, aambs, err)
      if (err%has_error()) return
      call mol%overlap(overlap)
      call aambs%overlap(s_mbs)
      call mixed_basis_overlap(mol, aambs, mixed, err)
      if (err%has_error()) return

      call aambs_dimensions(WATER_Z, 10, dims, err)
      if (err%has_error()) return
      call mo_aambs_overlap(scf%orbitals, mixed, s_mbs, projection, err)
      if (err%has_error()) return
      call valence_virtual_orbitals(scf%orbitals, projection, dims, vvo, err)
      if (err%has_error()) return

      call aambs_atom_ranges(WATER_Z, aambs, core_off, core_n, val_off, val_n, err)
      if (err%has_error()) return

      ! The valence-internal space: occupied orbitals that are not core,
      ! followed by the valence-virtual ones.
      n_ao = size(scf%orbitals, 1)
      n_valocc = dims%n_valocc
      allocate (valence_internal(n_ao, n_valocc + vvo%n_vvo))
      valence_internal(:, 1:n_valocc) = &
         scf%orbitals(:, dims%n_core + 1:dims%n_occupied)
      valence_internal(:, n_valocc + 1:) = vvo%orbitals

      call quasi_atomic_orbitals(WATER_Z, valence_internal, n_valocc, mixed, &
                                 val_off, val_n, quao, err)
      ok = .not. err%has_error()

      call mol%destroy()
      call aambs%destroy()
   end subroutine water_quaos

   subroutine test_orthonormal(error)
      !! Orthonormal, and there is one per valence minimal-basis orbital
      type(error_type), allocatable, intent(out) :: error
      type(quao_result_t) :: quao
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), work(:, :), gram(:, :)
      logical :: ok
      integer :: n, i

      call water_quaos("cc-pvdz", quao, overlap, dims, err, ok)
      call check(error, ok, "the water quasi-atomic construction should succeed")
      if (allocated(error)) return

      ! Oxygen contributes four valence orbitals and each hydrogen one.
      call check(error, quao%n_quao, 6, "water has six valence quasi-atomic orbitals")
      if (allocated(error)) return
      call check(error, quao%n_quao, dims%n_valence, &
                 "there is one quasi-atomic orbital per valence minimal-basis orbital")
      if (allocated(error)) return

      n = quao%n_quao
      allocate (work(size(overlap, 1), n), gram(n, n))
      call pic_gemm(overlap, quao%orbitals, work)
      call pic_gemm(quao%orbitals, work, gram, transa="T")
      do i = 1, n
         gram(i, i) = gram(i, i) - 1.0_dp
      end do
      call check(error, maxval(abs(gram)) < 1.0e-10_dp, &
                 "the quasi-atomic orbitals must be orthonormal")
   end subroutine test_orthonormal

   subroutine test_populations(error)
      !! The density in the quasi-atomic basis accounts for every valence electron
      !!
      !! Water has ten electrons, two of them in the oxygen 1s core, so the
      !! valence-internal space holds eight. The trace of the
      !! population-bond-order matrix is that count -- it is the same density
      !! written in a different orthonormal basis, and a trace is invariant to
      !! that. A failure means the transformation is not orthogonal.
      type(error_type), allocatable, intent(out) :: error
      type(quao_result_t) :: quao
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :)
      logical :: ok
      real(dp) :: trace
      integer :: i

      call water_quaos("cc-pvdz", quao, overlap, dims, err, ok)
      call check(error, ok, "the construction should succeed")
      if (allocated(error)) return

      trace = 0.0_dp
      do i = 1, quao%n_quao
         trace = trace + quao%population_bond_order(i, i)
      end do
      call check(error, abs(trace - 8.0_dp) < 1.0e-9_dp, &
                 "the valence populations must sum to the eight valence electrons")
      if (allocated(error)) return

      ! No orbital may hold more than a pair or fewer than none. Paper III
      ! quotes the intrinsic bound on a bond order as one, and the same
      ! reasoning bounds a population by two.
      do i = 1, quao%n_quao
         call check(error, quao%population_bond_order(i, i) > -1.0e-9_dp .and. &
                    quao%population_bond_order(i, i) < 2.0_dp + 1.0e-9_dp, &
                    "every population must lie between zero and two")
         if (allocated(error)) return
      end do

      call check(error, maxval(abs(quao%population_bond_order - &
                                   transpose(quao%population_bond_order))) < 1.0e-12_dp, &
                 "the population-bond-order matrix must be symmetric")
   end subroutine test_populations

   subroutine test_atomic_character(error)
      !! Each orbital really does sit on its own atom
      !!
      !! This is what the Paper II refinement maximizes, and the number is
      !! interpretable: one would mean every orbital lies entirely inside its own
      !! atom's free-atom space. It cannot reach one -- the orbitals are
      !! orthogonal to their neighbours, and that orthogonalization tail is
      !! exactly the deformation the molecule imposes -- but it should be close,
      !! and Paper I's tables of atomic contributions run 0.75 to 0.97.
      type(error_type), allocatable, intent(out) :: error
      type(quao_result_t) :: quao
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :)
      logical :: ok

      call water_quaos("cc-pvdz", quao, overlap, dims, err, ok)
      call check(error, ok, "the construction should succeed")
      if (allocated(error)) return

      call check(error, quao%atomic_character > 0.75_dp, &
                 "the quasi-atomic orbitals should be mostly on their own atoms")
      if (allocated(error)) return
      call check(error, quao%atomic_character < 1.0_dp + 1.0e-9_dp, &
                 "the atomic character is a sum of squared projections onto an "// &
                 "orthonormal set and cannot exceed one")
      if (allocated(error)) return

      ! It settles. GAMESS allows ten thousand sweeps; anything near that many
      ! would mean the functional is not being reduced.
      call check(error, quao%sweeps > 0 .and. quao%sweeps < 100, &
                 "the refinement should converge in a handful of sweeps")
      if (allocated(error)) return

      ! Two atoms' worth of oxygen orbitals and one per hydrogen, in that order.
      call check(error, quao%atom_of(1), 1, "the first orbitals belong to oxygen")
      if (allocated(error)) return
      call check(error, quao%atom_of(5), 2, "then one for each hydrogen")
      if (allocated(error)) return
      call check(error, quao%atom_of(6), 3, "then one for each hydrogen")
   end subroutine test_atomic_character

   subroutine test_water_pattern(error)
      !! Water comes out looking like water
      !!
      !! Before orientation the bonding is *there* but not yet concentrated.
      !! Paper I eq (5.4) is the reason: the sum of squares of each interatomic
      !! block of the population-bond-order matrix is invariant under rotations
      !! within an atom. So an O-H bond of 0.93 can sit as one element of 0.93
      !! or as two of 0.66 and 0.63, and which one it is has no physical content
      !! until the orientation step chooses. Concentrating it is precisely what
      !! orientation does, by maximizing the fourth powers -- which is only a
      !! meaningful thing to maximize because the squares cannot change.
      !!
      !! So the check here is the block norm, which is the invariant, and the
      !! per-element concentration belongs to the next stage. Magnitudes only:
      !! nothing has fixed a phase, and the papers are explicit that a negative
      !! bond order routinely describes a bonding interaction.
      type(error_type), allocatable, intent(out) :: error
      type(quao_result_t) :: quao
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :)
      logical :: ok
      real(dp) :: oh1, oh2, hh
      integer :: i, n_lone

      call water_quaos("cc-pvdz", quao, overlap, dims, err, ok)
      call check(error, ok, "the construction should succeed")
      if (allocated(error)) return

      ! Orbitals 1-4 are oxygen's, 5 and 6 are the hydrogens'.
      oh1 = sqrt(sum(quao%population_bond_order(5, 1:4)**2))
      oh2 = sqrt(sum(quao%population_bond_order(6, 1:4)**2))
      hh = abs(quao%population_bond_order(5, 6))

      ! Paper I reports O-H bond orders of 0.97 in quinone and Paper IV
      ! 0.96-0.98 along a whole reaction path, so a bonded pair should land
      ! near one however the block is distributed.
      call check(error, oh1 > 0.85_dp .and. oh1 < 1.05_dp, &
                 "the first O-H block norm should be a single bond")
      if (allocated(error)) return
      call check(error, oh2 > 0.85_dp .and. oh2 < 1.05_dp, &
                 "the second O-H block norm should be a single bond")
      if (allocated(error)) return

      ! Symmetric molecule, symmetric bonds.
      call check(error, abs(oh1 - oh2) < 1.0e-6_dp, &
                 "the two O-H bonds are related by symmetry and should match")
      if (allocated(error)) return

      ! The hydrogens are bonded to oxygen, not to each other.
      call check(error, hh < 0.1_dp, &
                 "the two hydrogens should not be bonded to each other")
      if (allocated(error)) return

      ! Two of oxygen's four valence orbitals are lone pairs.
      n_lone = 0
      do i = 1, 4
         if (quao%population_bond_order(i, i) > 1.5_dp) n_lone = n_lone + 1
      end do
      call check(error, n_lone >= 2, &
                 "oxygen should carry two nearly doubly occupied lone pairs")
   end subroutine test_water_pattern

end module test_mqc_quao

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_quao, only: collect_mqc_quao_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_quao", collect_mqc_quao_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
