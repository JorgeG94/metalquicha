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
                               aambs_atom_ranges, quasi_atomic_orbitals, quao_result_t, &
                               orient_quasi_atomic_orbitals, kinetic_bond_orders, &
                               split_localize
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
                  new_unittest("water_bonding_pattern", test_water_pattern), &
                  new_unittest("orientation_concentrates_bonding", test_orientation), &
                  new_unittest("kinetic_bond_orders", test_kbo), &
                  new_unittest("split_localization_pairs_bonds", test_split) &
                  ]
   end subroutine collect_mqc_quao_tests

   subroutine water_quaos(basis_name, quao, overlap, dims, err, ok, vi_out)
      !! Everything from an SCF through to the quasi-atomic orbitals
      character(len=*), intent(in) :: basis_name
      type(quao_result_t), intent(out) :: quao
      real(dp), allocatable, intent(out) :: overlap(:, :)
      type(aambs_dimensions_t), intent(out) :: dims
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok
      real(dp), allocatable, intent(out), optional :: vi_out(:, :)

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
      if (present(vi_out)) vi_out = valence_internal

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

   subroutine test_orientation(error)
      !! Orientation concentrates the bonding without creating any
      !!
      !! Two claims, and the pair of them is the whole content of Paper I
      !! eqs (5.4) and (5.5). The interatomic block norms are *invariant* --
      !! rotating inside an atom cannot change how much bonding there is between
      !! two atoms -- while the fourth-power sum *increases*, because the same
      !! total gets concentrated into fewer, larger elements.
      !!
      !! Before orientation each water hydrogen couples to two oxygen orbitals
      !! at roughly 0.66 and 0.63; after it, one of them should carry almost the
      !! whole bond.
      type(error_type), allocatable, intent(out) :: error
      type(quao_result_t) :: quao
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :)
      logical :: ok
      real(dp) :: before_norm, after_norm, before_sum, after_sum, best
      integer :: i

      call water_quaos("cc-pvdz", quao, overlap, dims, err, ok)
      call check(error, ok, "the construction should succeed")
      if (allocated(error)) return

      before_norm = sqrt(sum(quao%population_bond_order(5, 1:4)**2))
      before_sum = 0.0_dp
      do i = 1, 4
         before_sum = before_sum + quao%population_bond_order(5, i)**4
      end do

      call orient_quasi_atomic_orbitals(quao, err)
      call check(error,.not. err%has_error(), "orientation should succeed")
      if (allocated(error)) return
      call check(error, quao%oriented, "the result should be marked oriented")
      if (allocated(error)) return

      after_norm = sqrt(sum(quao%population_bond_order(5, 1:4)**2))
      after_sum = 0.0_dp
      do i = 1, 4
         after_sum = after_sum + quao%population_bond_order(5, i)**4
      end do

      ! eq (5.4): the block norm cannot move.
      call check(error, abs(after_norm - before_norm) < 1.0e-9_dp, &
                 "orientation must not change how much bonding there is between "// &
                 "two atoms, only how it is distributed")
      if (allocated(error)) return

      ! eq (5.5): the fourth-power sum is what it maximizes.
      call check(error, after_sum > before_sum - 1.0e-12_dp, &
                 "the fourth-power sum should not decrease")
      if (allocated(error)) return

      ! And the point of it all: one oxygen orbital now carries the O-H bond.
      best = 0.0_dp
      do i = 1, 4
         best = max(best, abs(quao%population_bond_order(5, i)))
      end do
      call check(error, best > 0.85_dp, &
                 "after orientation one oxygen orbital should carry the O-H bond, "// &
                 "rather than the bond being spread over several")
      if (allocated(error)) return

      ! The trace is a property of the density, not of the basis it is written in.
      before_sum = 0.0_dp
      do i = 1, quao%n_quao
         before_sum = before_sum + quao%population_bond_order(i, i)
      end do
      call check(error, abs(before_sum - 8.0_dp) < 1.0e-9_dp, &
                 "orientation must not change the electron count")
   end subroutine test_orientation

   subroutine test_kbo(error)
      !! Kinetic bond orders: negative for bonds, and on the published scale
      type(error_type), allocatable, intent(out) :: error
      type(quao_result_t) :: quao
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      type(libcint_molecule_t) :: mol
      real(dp), allocatable :: overlap(:, :), kinetic(:, :), kbo(:, :)
      logical :: ok
      real(dp) :: strongest
      integer :: i, best

      call water_quaos("cc-pvdz", quao, overlap, dims, err, ok)
      call check(error, ok, "the construction should succeed")
      if (allocated(error)) return

      ! Refused before orientation, because the individual couplings are
      ! arbitrary until then.
      call kinetic_bond_orders(quao, overlap, kbo, err)
      call check(error, err%has_error(), &
                 "kinetic bond orders should be refused before orientation")
      if (allocated(error)) return
      call err%clear()

      call orient_quasi_atomic_orbitals(quao, err)
      call check(error,.not. err%has_error(), "orientation should succeed")
      if (allocated(error)) return

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "the molecule should rebuild")
      if (allocated(error)) return
      call mol%kinetic(kinetic)

      call kinetic_bond_orders(quao, kinetic, kbo, err)
      call check(error,.not. err%has_error(), "kinetic bond orders should build")
      if (allocated(error)) return

      ! The O-H bond. Paper IV reports C-H kinetic bond orders of -35 to -37
      ! kcal/mol across a whole reaction path and Paper III -42.9 for N-H, so an
      ! O-H bond should land in the same tens-of-kcal/mol range and be negative.
      best = 1
      strongest = 0.0_dp
      do i = 1, 4
         if (abs(quao%population_bond_order(5, i)) > strongest) then
            strongest = abs(quao%population_bond_order(5, i))
            best = i
         end if
      end do
      call check(error, kbo(5, best) < -10.0_dp .and. kbo(5, best) > -200.0_dp, &
                 "the O-H kinetic bond order should be bonding and on the scale of "// &
                 "a chemical bond energy")
      if (allocated(error)) return

      ! Symmetric, like the matrix it is built from.
      call check(error, abs(kbo(5, best) - kbo(best, 5)) < 1.0e-9_dp, &
                 "the kinetic bond order matrix should be symmetric")
      if (allocated(error)) return

      ! The two O-H bonds are equivalent by symmetry. Not to machine precision:
      ! the Jacobi sweep visits pairs in index order, so two orbitals related by
      ! a symmetry the algorithm knows nothing about converge to the same answer
      ! from different directions and stop a rotation apart.
      strongest = 0.0_dp
      do i = 1, 4
         strongest = min(strongest, kbo(6, i))
      end do
      call check(error, abs(strongest - kbo(5, best)) < 1.0e-4_dp, &
                 "the two O-H kinetic bond orders should match by symmetry")
      if (allocated(error)) return

      ! Water's shape, in one assertion: two lone pairs holding a pair each and
      ! coupling to nothing, two bonds of 0.93, and the H->O charge transfer
      ! that makes each bonding pair sum to exactly two electrons.
      call check(error, abs(quao%population_bond_order(5, 5) + &
                            quao%population_bond_order(best, best) - 2.0_dp) < 0.01_dp, &
                 "the two orbitals forming an O-H bond should hold two electrons "// &
                 "between them, however unevenly they share it")
      call mol%destroy()
   end subroutine test_kbo

   subroutine test_split(error)
      !! Localizing the two blocks separately pairs bonds with antibonds
      !!
      !! Water's valence-internal space is eight orbitals: four occupied and two
      !! empty. Localizing them separately should give two O-H bonds and two
      !! lone pairs among the occupied, and two O-H antibonds among the empty --
      !! and each antibond should sit on the same two atoms as a bond.
      !!
      !! The assertion is that every localized orbital is concentrated: the
      !! criterion counts how many quasi-atomic orbitals each one covers, so a
      !! successful localization leaves each with one or two large components
      !! and the rest near zero. That is checkable without knowing which orbital
      !! turned out to be which.
      type(error_type), allocatable, intent(out) :: error
      type(quao_result_t) :: quao
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), vi(:, :), localized(:, :)
      real(dp), allocatable :: weights(:), work(:, :)
      logical :: ok
      real(dp) :: spread_before, spread_after, top_two
      integer :: n, i, j

      call water_quaos("cc-pvdz", quao, overlap, dims, err, ok, vi)
      call check(error, ok, "the construction should succeed")
      if (allocated(error)) return

      ! Refused before orientation: the criterion is stated against the oriented
      ! quasi-atomic orbitals.
      call split_localize(quao, dims%n_valocc, vi, localized, err)
      call check(error, err%has_error(), &
                 "split localization should be refused before orientation")
      if (allocated(error)) return
      call err%clear()

      call orient_quasi_atomic_orbitals(quao, err)
      call check(error,.not. err%has_error(), "orientation should succeed")
      if (allocated(error)) return

      ! How concentrated the *unlocalized* orbitals are, for comparison.
      n = quao%n_quao
      spread_before = 0.0_dp
      do i = 1, dims%n_valocc
         spread_before = spread_before + sum(quao%to_valence_internal(i, :)**4)
      end do

      call split_localize(quao, dims%n_valocc, vi, localized, err)
      call check(error,.not. err%has_error(), "split localization should succeed")
      if (allocated(error)) return
      call check(error, size(localized, 2), n, &
                 "localization produces one orbital per valence-internal orbital")
      if (allocated(error)) return

      ! Still an orthonormal set: the transformation is orthogonal, so the
      ! space is unchanged and only the description of it moved.
      allocate (work(size(overlap, 1), n), weights(n))
      call pic_gemm(overlap, localized, work)
      do i = 1, n
         weights(i) = dot_product(localized(:, i), work(:, i))
      end do
      call check(error, maxval(abs(weights - 1.0_dp)) < 1.0e-9_dp, &
                 "the localized orbitals must stay normalized")
      if (allocated(error)) return

      ! Each localized orbital should be carried by one or two quasi-atomic
      ! orbitals. A bond covers two; a lone pair covers one.
      do i = 1, n
         weights = 0.0_dp
         do j = 1, n
            weights(j) = dot_product(quao%orbitals(:, j), work(:, i))**2
         end do
         call sort_descending(weights)
         top_two = weights(1) + weights(2)
         call check(error, top_two > 0.85_dp, &
                    "every localized orbital should be carried by at most two "// &
                    "quasi-atomic orbitals -- a bond covers two, a lone pair one")
         if (allocated(error)) return
      end do

      ! And the criterion it maximizes did go up.
      spread_after = 0.0_dp
      do i = 1, dims%n_valocc
         do j = 1, n
            spread_after = spread_after + &
                           dot_product(quao%orbitals(:, j), work(:, i))**4
         end do
      end do
      call check(error, spread_after > spread_before - 1.0e-10_dp, &
                 "the fourth-power sum should not decrease")
   end subroutine test_split

   subroutine sort_descending(a)
      !! Small insertion sort, largest first
      real(dp), intent(inout) :: a(:)
      real(dp) :: t
      integer :: i, j

      do i = 2, size(a)
         t = a(i)
         j = i - 1
         do while (j >= 1)
            if (a(j) >= t) exit
            a(j + 1) = a(j)
            j = j - 1
         end do
         a(j + 1) = t
      end do
   end subroutine sort_descending

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
