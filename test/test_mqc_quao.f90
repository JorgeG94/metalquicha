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
   use mqc_aambs, only: aambs_dimensions, aambs_dimensions_t, aambs_element_counts
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

   !> Formyl chloride, from the GAMESS test deck that ships reference QUAO
   !> populations. Angstrom in the deck, Bohr here.
   integer, parameter :: FC_Z(4) = [6, 8, 17, 1]
   character(len=2), parameter :: FC_SYM(4) = ["C ", "O ", "Cl", "H "]
   real(dp), parameter :: FC(3, 4) = reshape([ &
                                             -1.3271546579_dp, -0.8536648799_dp, 0.0000000000_dp, &
                                             -3.0122423412_dp, 0.5424647816_dp, 0.0000000000_dp, &
                                             1.8500985686_dp, 0.1267628285_dp, 0.0000000000_dp, &
                                             -1.4456215887_dp, -2.8971202257_dp, 0.0000000000_dp], [3, 4])

   integer, parameter :: BQ_Z(12) = [6, 6, 6, 6, 6, 6, 8, 8, 1, 1, 1, 1]
   character(len=2), parameter :: BQ_SYM(12) = ["C ", "C ", "C ", "C ", "C ", "C ", "O ", "O ", "H ", "H ", "H ", "H "]
   real(dp), parameter :: BQ(3, 12) = reshape([ &
                                              1.6803690372_dp, 2.1213574156_dp, -0.1919168058_dp, &
                                              2.6490048546_dp, -0.4712787984_dp, 0.1259370182_dp, &
                                              1.0841264295_dp, -2.4470479338_dp, 0.3038963069_dp, &
                                              -1.6803709269_dp, -2.1213479669_dp, 0.1919394826_dp, &
                                              -2.6490086340_dp, 0.4712844676_dp, -0.1259067826_dp, &
                                              -1.0841320987_dp, 2.4470536030_dp, -0.3038641815_dp, &
                                              -3.1157709870_dp, -3.9333250904_dp, 0.3550870979_dp, &
                                              3.1157728768_dp, 3.9333175315_dp, -0.3551985917_dp, &
                                              4.6825807122_dp, -0.6674682750_dp, 0.2038069626_dp, &
                                              1.7542478801_dp, -4.3645190183_dp, 0.5366595429_dp, &
                                              -1.7542592184_dp, 4.3645227977_dp, -0.5366519840_dp, &
                                              -4.6825807122_dp, 0.6674701648_dp, -0.2038031832_dp], [3, 12])

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
                  new_unittest("split_localization_pairs_bonds", test_split), &
                  new_unittest("benzoquinone_against_paper_one", test_benzoquinone), &
                  new_unittest("formyl_chloride_against_gamess", test_formyl_chloride) &
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
      real(dp), allocatable :: interference(:, :), vi(:, :), work(:, :)
      logical :: ok
      real(dp) :: strongest, trace_mo
      integer :: i, best

      call water_quaos("cc-pvdz", quao, overlap, dims, err, ok, vi)
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

   subroutine test_benzoquinone(error)
      !! Paper I Figure 5 and Table IV, on 1,4-benzoquinone
      !!
      !! The first check against published chemistry rather than against an
      !! invariant. Twelve atoms, 44 minimal-basis orbitals, 16 valence-virtual
      !! orbitals to recover -- enough that the per-atom machinery is actually
      !! exercised, where water needed almost none of it.
      !!
      !! The geometry is not Paper I's: it gives none, saying only
      !! "HF-optimized", so this uses an ordinary molecular-mechanics structure.
      !! Bond orders and populations are therefore expected to agree to about
      !! the second decimal rather than exactly -- which is still a real test,
      !! because every one of these numbers is a specific published value and
      !! the pattern of six distinct bond types is not something a broken
      !! implementation reproduces by accident.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol, aambs
      type(rhf_result_t) :: scf
      type(vvo_result_t) :: vvo
      type(quao_result_t) :: quao
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: mixed(:, :), s_mbs(:, :), proj(:, :), vi(:, :)
      real(dp), allocatable :: pop(:)
      integer, allocatable :: co(:), cn(:), vo(:), vn(:)
      integer :: i, core, valence

      call build_libcint_molecule(BQ_Z, BQ_SYM, BQ, "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "benzoquinone should build")
      if (allocated(error)) return
      call run_libcint_rhf(mol, 56, 300, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err)
      call check(error, scf%converged, "the benzoquinone SCF should converge")
      if (allocated(error)) return

      call build_aambs_molecule(BQ_Z, BQ_SYM, BQ, aambs, err)
      call aambs%overlap(s_mbs)
      call mixed_basis_overlap(mol, aambs, mixed, err)
      call aambs_dimensions(BQ_Z, 56, dims, err)
      call check(error,.not. err%has_error(), "the setup should succeed")
      if (allocated(error)) return

      ! Six carbons at five minimal-basis orbitals each, two oxygens at five,
      ! four hydrogens at one.
      call check(error, dims%n_mbs, 44, "benzoquinone has 44 minimal-basis orbitals")
      if (allocated(error)) return
      call check(error, dims%n_vvo, 16, "sixteen valence-virtual orbitals to recover")
      if (allocated(error)) return

      call mo_aambs_overlap(scf%orbitals, mixed, s_mbs, proj, err)
      call valence_virtual_orbitals(scf%orbitals, proj, dims, vvo, err)
      call check(error,.not. err%has_error(), "the valence-virtual step should succeed")
      if (allocated(error)) return

      ! Paper I Table I quotes 0.99999 retained against 0.105-0.272 rejected
      ! across its eight molecules; this should land in the same band.
      call check(error, vvo%smallest_retained > 0.99_dp, &
                 "the retained singular values should be near one")
      if (allocated(error)) return
      call check(error, vvo%largest_rejected < 0.3_dp, &
                 "the rejected singular values should be well below")
      if (allocated(error)) return

      call aambs_atom_ranges(BQ_Z, aambs, co, cn, vo, vn, err)
      allocate (vi(size(scf%orbitals, 1), dims%n_valocc + vvo%n_vvo))
      vi(:, 1:dims%n_valocc) = scf%orbitals(:, dims%n_core + 1:dims%n_occupied)
      vi(:, dims%n_valocc + 1:) = vvo%orbitals
      call quasi_atomic_orbitals(BQ_Z, vi, dims%n_valocc, mixed, vo, vn, quao, err)
      call orient_quasi_atomic_orbitals(quao, err)
      call check(error,.not. err%has_error(), "the quasi-atomic steps should succeed")
      if (allocated(error)) return

      ! Table IV, cc-pVDZ column: O 8.442, carbonyl C 5.568, CH carbon 6.163,
      ! H 0.832. Core electrons included, which Paper I does not say outright but
      ! its own sum to 56 electrons requires.
      allocate (pop(size(BQ_Z)))
      pop = 0.0_dp
      do i = 1, quao%n_quao
         pop(quao%atom_of(i)) = pop(quao%atom_of(i)) + quao%population_bond_order(i, i)
      end do
      do i = 1, size(BQ_Z)
         call aambs_element_counts(BQ_Z(i), core, valence, err)
         pop(i) = pop(i) + 2.0_dp*real(core, dp)
      end do

      call check(error, abs(pop(1) - 5.568_dp) < 0.05_dp, &
                 "carbonyl carbon population, Paper I Table IV")
      if (allocated(error)) return
      call check(error, abs(pop(2) - 6.163_dp) < 0.05_dp, &
                 "CH carbon population, Paper I Table IV")
      if (allocated(error)) return
      call check(error, abs(pop(7) - 8.442_dp) < 0.05_dp, &
                 "oxygen population, Paper I Table IV")
      if (allocated(error)) return
      call check(error, abs(pop(9) - 0.832_dp) < 0.05_dp, &
                 "hydrogen population, Paper I Table IV")
      if (allocated(error)) return

      ! Every electron accounted for.
      call check(error, abs(sum(pop) - 56.0_dp) < 1.0e-8_dp, &
                 "the populations must sum to benzoquinone's 56 electrons")
      if (allocated(error)) return

      ! Figure 5: the six bond types, by magnitude. Assert |p| -- phases are a
      ! convention, and the papers note negative bond orders describe bonding.
      call check(error, best_bond(quao, 1, 8) > 0.95_dp, &
                 "C-O sigma bond order, Paper I Figure 5 gives 0.97")
      if (allocated(error)) return
      call check(error, best_bond(quao, 2, 3) > 0.97_dp, &
                 "C-C sigma bond order, Paper I Figure 5 gives 0.99")
      if (allocated(error)) return
      call check(error, best_bond(quao, 2, 9) > 0.95_dp, &
                 "C-H bond order, Paper I Figure 5 gives 0.97")
      if (allocated(error)) return

      ! The two oxygen lone pairs, 1.99 and 1.90 in Figure 5.
      call check(error, count_above(quao, 7, 1.98_dp) >= 1, &
                 "oxygen should carry a nearly full sigma lone pair")
      if (allocated(error)) return
      call check(error, count_above(quao, 7, 1.85_dp) >= 2, &
                 "oxygen should carry two lone pairs, Paper I Figure 5 gives "// &
                 "1.99 and 1.90")

      call mol%destroy()
      call aambs%destroy()
   end subroutine test_benzoquinone

   function best_bond(quao, atom_a, atom_b) result(largest)
      !! The largest bond-order magnitude between two atoms
      type(quao_result_t), intent(in) :: quao
      integer, intent(in) :: atom_a, atom_b
      real(dp) :: largest
      integer :: i, j

      largest = 0.0_dp
      do i = 1, quao%n_quao
         if (quao%atom_of(i) /= atom_a) cycle
         do j = 1, quao%n_quao
            if (quao%atom_of(j) /= atom_b) cycle
            largest = max(largest, abs(quao%population_bond_order(i, j)))
         end do
      end do
   end function best_bond

   function count_above(quao, atom, threshold) result(n)
      !! How many of an atom's orbitals hold more than `threshold` electrons
      type(quao_result_t), intent(in) :: quao
      integer, intent(in) :: atom
      real(dp), intent(in) :: threshold
      integer :: n, i

      n = 0
      do i = 1, quao%n_quao
         if (quao%atom_of(i) /= atom) cycle
         if (quao%population_bond_order(i, i) > threshold) n = n + 1
      end do
   end function count_above

   subroutine test_formyl_chloride(error)
      !! Against GAMESS's own reference QUAO populations
      !!
      !! From the GAMESS test deck for formyl chloride, which ships the answer
      !! in its header:
      !!
      !!     ATOM     QUAO POP.      CHARGE
      !!     1 C       5.639027     0.360973
      !!     2 O       8.384177    -0.384177
      !!     3 CL     17.116362    -0.116362
      !!     4 H       0.860434     0.139566
      !!
      !! Unlike the benzoquinone comparison this fixes the geometry, the basis
      !! and the method exactly -- RHF, cc-pVDZ, spherical -- so it is a
      !! genuine head-to-head against the reference implementation rather than
      !! a comparison against a figure whose structure is unavailable. It is
      !! also the first element past the second row: chlorine's chemical core is
      !! five orbitals, so ten of its electrons never enter the valence
      !! construction and have to be added back.
      !!
      !! The deck notes the populations are the same for `local=svd` and
      !! `local=ruednbrg`, which is worth knowing: they are a property of the
      !! quasi-atomic space, not of how it was localized.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol, aambs
      type(rhf_result_t) :: scf
      type(vvo_result_t) :: vvo
      type(quao_result_t) :: quao
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: mixed(:, :), s_mbs(:, :), proj(:, :), vi(:, :)
      real(dp), allocatable :: pop(:)
      integer, allocatable :: co(:), cn(:), vo(:), vn(:)
      integer :: i, core, valence
      real(dp), parameter :: REFERENCE(4) = &
                             [5.639027_dp, 8.384177_dp, 17.116362_dp, 0.860434_dp]

      call build_libcint_molecule(FC_Z, FC_SYM, FC, "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "formyl chloride should build")
      if (allocated(error)) return
      call run_libcint_rhf(mol, 32, 300, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err)
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return

      call build_aambs_molecule(FC_Z, FC_SYM, FC, aambs, err)
      call aambs%overlap(s_mbs)
      call mixed_basis_overlap(mol, aambs, mixed, err)
      call aambs_dimensions(FC_Z, 32, dims, err)
      call check(error,.not. err%has_error(), "the setup should succeed")
      if (allocated(error)) return

      ! C and O give five each, chlorine nine, hydrogen one.
      call check(error, dims%n_mbs, 20, "formyl chloride has 20 minimal-basis orbitals")
      if (allocated(error)) return
      call check(error, dims%n_core, 7, "seven chemical-core orbitals, five of them "// &
                 "on chlorine")
      if (allocated(error)) return

      call mo_aambs_overlap(scf%orbitals, mixed, s_mbs, proj, err)
      call valence_virtual_orbitals(scf%orbitals, proj, dims, vvo, err)
      call aambs_atom_ranges(FC_Z, aambs, co, cn, vo, vn, err)
      call check(error,.not. err%has_error(), "the valence-virtual step should succeed")
      if (allocated(error)) return

      allocate (vi(size(scf%orbitals, 1), dims%n_valocc + vvo%n_vvo))
      vi(:, 1:dims%n_valocc) = scf%orbitals(:, dims%n_core + 1:dims%n_occupied)
      vi(:, dims%n_valocc + 1:) = vvo%orbitals
      call quasi_atomic_orbitals(FC_Z, vi, dims%n_valocc, mixed, vo, vn, quao, err)
      call orient_quasi_atomic_orbitals(quao, err)
      call check(error,.not. err%has_error(), "the quasi-atomic steps should succeed")
      if (allocated(error)) return

      allocate (pop(4))
      pop = 0.0_dp
      do i = 1, quao%n_quao
         pop(quao%atom_of(i)) = pop(quao%atom_of(i)) + quao%population_bond_order(i, i)
      end do
      do i = 1, 4
         call aambs_element_counts(FC_Z(i), core, valence, err)
         pop(i) = pop(i) + 2.0_dp*real(core, dp)
      end do

      ! Every electron, first: a population analysis that does not conserve
      ! charge is wrong before any individual number is looked at.
      call check(error, abs(sum(pop) - 32.0_dp) < 1.0e-8_dp, &
                 "the populations must sum to the 32 electrons of HCOCl")
      if (allocated(error)) return

      ! Agreement is 0.004 at worst, on chlorine. That is not machine
      ! precision and should not be: GAMESS Loewdin-orthogonalizes the free-atom
      ! basis before the valence-virtual projection -- a step added in 2013 and
      ! absent from the published equations -- and its valence and quasi-atomic
      ! halves disagree about whether to use the rotated set or the raw one.
      ! Those choices are recorded in the plan; the residual here is their size.
      do i = 1, 4
         call check(error, abs(pop(i) - REFERENCE(i)) < 0.01_dp, &
                    "atom "//trim(FC_SYM(i))//" population against GAMESS")
         if (allocated(error)) return
      end do

      call mol%destroy()
      call aambs%destroy()
   end subroutine test_formyl_chloride

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
