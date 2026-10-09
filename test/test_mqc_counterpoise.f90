!! Ghost centres, and the basis-set superposition error they exist to expose
!!
!! A counterpoise-corrected many-body term computes each monomer in the *pair's*
!! basis rather than its own. These are the properties that makes rest on, all
!! checkable before any expansion is assembled:
!!
!!   * a ghost centre carries basis functions and no nucleus
!!   * ghosting changes nothing about the AO space -- same count, same ordering
!!   * a monomer in the pair's basis is *lower* than in its own, and the gap is
!!     the BSSE, which is the whole quantity counterpoise removes
!!
!! The last one is the point. Without it a "counterpoise correction" could be
!! wired up, run, and produce a number that is simply the uncorrected one.
!!
!! The second half checks the full-cluster-basis scheme (SSFC) against the
!! identity it rests on. At L = N the many-body sum telescopes to
!!
!!     E_SSFC(N) = E_cluster + sum_i [ E_i(i) - E_i(FB) ]
!!
!! the cluster plus the superposition error of each monomer, where `E_i(i)` is
!! monomer i in its own basis and `E_i(FB)` the same monomer in the basis of the
!! whole cluster. Every side of that is assembled twice, once by the expansion
!! and once by hand from single calculations built with explicit ghost masks.
module test_mqc_counterpoise
   use pic_types, only: dp, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, &
                                    build_fragment_from_indices
   use mqc_combinatorics, only: vmfc_subset_key, vmfc_row_subset_key, is_auxiliary_row, &
                                real_count_of, counterpoise_scheme_of, &
                                COUNTERPOISE_VMFC, COUNTERPOISE_SSFC
   use mqc_error, only: error_t
   use mqc_mbe, only: compute_mbe
   use mqc_result_types, only: calculation_result_t, mbe_result_t
   use mqc_frag_utils, only: generate_mbe_term_list
   use mqc_config_adapter, only: driver_config_t
   use mqc_calc_types, only: CALC_TYPE_ENERGY
   use mqc_method_types, only: METHOD_TYPE_HF
   use mqc_json_output_types, only: json_output_data_t
   use mqc_mbe_fragment_distribution_scheme, only: serial_fragment_processor
   use pic_logger, only: logger => global_logger, error_level
   use mqc_elements, only: core_orbital_count
   use mqc_physical_constants, only: PI
   use pic_io, only: to_char
   implicit none
   private

   public :: collect_counterpoise

   real(dp), parameter :: ANG = 1.8897261254578281_dp

   !! Two waters, 3 Angstrom apart along x. Far enough that the interaction is
   !! small and close enough that each borrows the other's functions.
   real(dp), parameter :: SEP = 3.0_dp

   integer, parameter :: N_ATOMS = 6
   integer, parameter :: N_MONOMER = 3

   !! SAPT's counterpoise-corrected supermolecular Hartree-Fock interaction
   !! energy for this dimer in 6-31G, reached through the dimer-centred basis
   !! and none of this code. An outside number, not one this program printed.
   real(dp), parameter :: SAPT_E_INT_HF = 0.009315543671_dp

contains

   subroutine collect_counterpoise(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("a_ghost_carries_no_nucleus", test_no_nucleus), &
                  new_unittest("ghosting_preserves_the_ao_space", test_ao_space), &
                  new_unittest("the_pair_basis_lowers_a_monomer", test_bsse_is_real), &
                  new_unittest("an_isolated_ghost_changes_nothing", test_far_ghost), &
                  new_unittest("a_signed_index_ghosts_its_monomer", test_signed_indices), &
                  new_unittest("a_ghost_has_no_core_to_freeze", test_ghost_core), &
                  new_unittest("vmfc_reproduces_the_supermolecule", test_vmfc_identity), &
                  new_unittest("the_subset_key_ghosts_the_complement", test_subset_key), &
                  new_unittest("a_ghosted_row_keeps_its_ghosts", test_row_subset_key), &
                  new_unittest("an_auxiliary_row_is_never_summed", test_auxiliary), &
                  new_unittest("bsse_shrinks_as_the_monomers_separate", test_bsse_decays), &
                  new_unittest("ssfc_of_a_generated_list_is_the_boys_bernardi_total", &
                               test_ssfc_synthetic_identity), &
                  new_unittest("ssfc_of_a_water_trimer_is_the_boys_bernardi_total", &
                               test_ssfc_trimer_identity), &
                  new_unittest("ssfc_of_a_water_dimer_is_vmfc", test_ssfc_dimer_is_vmfc), &
                  new_unittest("ssfc_two_body_is_not_vmfc_on_a_trimer", test_ssfc_is_not_vmfc), &
                  new_unittest("ssfc_interaction_energy_is_the_references_full_basis_terms", &
                               test_ssfc_interaction_energy) &
                  ]
   end subroutine collect_counterpoise

   subroutine dimer_geometry(z, sym, c)
      !! Two waters, monomer A first then monomer B
      integer, intent(out) :: z(N_ATOMS)
      character(len=2), intent(out) :: sym(N_ATOMS)
      real(dp), intent(out) :: c(3, N_ATOMS)

      integer :: i

      z = [8, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "O ", "H ", "H "]
      c(:, 1) = [0.0_dp, 0.0_dp, 0.10077199_dp]
      c(:, 2) = [0.0_dp, 0.77250895_dp, -0.46780200_dp]
      c(:, 3) = [0.0_dp, -0.77250895_dp, -0.46780200_dp]
      do i = 1, 3
         c(:, i + 3) = c(:, i)
         c(1, i + 3) = c(1, i) + SEP
      end do
      c = c*ANG
   end subroutine dimer_geometry

   subroutine test_no_nucleus(error)
      !! Ghosting monomer B removes its charge and leaves its functions
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: dimer, ghosted
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6)
      logical :: ghost(6)

      call dimer_geometry(z, sym, c)
      ghost = [.false., .false., .false., .true., .true., .true.]

      call build_czt_molecule(z, sym, c, "sto-3g", dimer, err)
      call check(error,.not. err%has_error(), "dimer: "//err%get_full_trace())
      if (allocated(error)) return

      call build_czt_molecule(z, sym, c, "sto-3g", ghosted, err, ghost=ghost)
      call check(error,.not. err%has_error(), "ghosted: "//err%get_full_trace())
      if (allocated(error)) return

      ! Ten electrons' worth of nucleus gone -- one water -- and nothing else.
      call check(error, nint(sum(dimer%charges)), 20, &
                 "the dimer should carry twenty protons")
      if (allocated(error)) return
      call check(error, nint(sum(ghosted%charges)), 10, &
                 "ghosting one water should leave ten")

      call dimer%destroy()
      call ghosted%destroy()
   end subroutine test_no_nucleus

   subroutine test_ao_space(error)
      !! Same number of basis functions, in the same order
      !!
      !! This is the invariant the whole correction rests on. If ghosting moved
      !! or dropped a function, a monomer-in-pair-basis energy would not be
      !! comparable with the pair's and the difference would be meaningless.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: dimer, ghost_a, ghost_b
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6)
      real(dp), allocatable :: s_dimer(:, :), s_ghost(:, :)

      call dimer_geometry(z, sym, c)
      call build_czt_molecule(z, sym, c, "sto-3g", dimer, err)
      call build_czt_molecule(z, sym, c, "sto-3g", ghost_a, err, &
                              ghost=[.true., .true., .true., .false., .false., .false.])
      call build_czt_molecule(z, sym, c, "sto-3g", ghost_b, err, &
                              ghost=[.false., .false., .false., .true., .true., .true.])
      call check(error,.not. err%has_error(), "build: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, ghost_a%nao, dimer%nao, &
                 "ghosting A changed the number of basis functions")
      if (allocated(error)) return
      call check(error, ghost_b%nao, dimer%nao, &
                 "ghosting B changed the number of basis functions")
      if (allocated(error)) return

      ! Ordering too, not just the count: the overlap is built from the same
      ! functions in the same places, so it must come back identical.
      call dimer%overlap(s_dimer)
      call ghost_b%overlap(s_ghost)
      call check(error, maxval(abs(s_dimer - s_ghost)), 0.0_dp, &
                 "ghosting changed the overlap matrix, so the AO ordering moved", &
                 thr=1.0e-14_dp)

      call dimer%destroy()
      call ghost_a%destroy()
      call ghost_b%destroy()
   end subroutine test_ao_space

   subroutine test_bsse_is_real(error)
      !! A monomer in the pair's basis lies below the same monomer alone
      !!
      !! The gap is the basis-set superposition error: monomer A, having
      !! borrowed B's functions, describes itself better than its own basis
      !! allows. In a plain expansion that borrowing lands in the pair term and
      !! nowhere else, which is why the dimer looks more bound than it is.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: alone, in_pair
      type(rhf_result_t) :: scf_alone, scf_pair
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6), bsse

      call dimer_geometry(z, sym, c)

      ! Monomer A in its own basis: three atoms, nothing else.
      call build_czt_molecule(z(1:3), sym(1:3), c(:, 1:3), "sto-3g", alone, err)
      call check(error,.not. err%has_error(), "alone: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_rhf(alone, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_alone, err)
      call check(error,.not. err%has_error(), "alone SCF: "//err%get_full_trace())
      if (allocated(error)) return

      ! The same monomer, same ten electrons, in the pair's basis.
      call build_czt_molecule(z, sym, c, "sto-3g", in_pair, err, &
                              ghost=[.false., .false., .false., .true., .true., .true.])
      call check(error,.not. err%has_error(), "in pair: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_rhf(in_pair, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_pair, err)
      call check(error,.not. err%has_error(), "pair SCF: "//err%get_full_trace())
      if (allocated(error)) return

      bsse = scf_alone%energy - scf_pair%energy

      ! Strictly lower. A variational method handed more functions cannot do
      ! worse, and at 3 Angstrom in a minimal basis it does measurably better.
      call check(error, bsse > 0.0_dp, &
                 "the pair basis did not lower the monomer, so the ghost "// &
                 "functions are not reaching the SCF")
      if (allocated(error)) return

      ! And it is a real effect rather than convergence noise: sto-3g at this
      ! separation is worth more than a microhartree and less than a hartree.
      call check(error, bsse > 1.0e-6_dp, &
                 "the lowering is too small to be superposition error")
      if (allocated(error)) return
      call check(error, bsse < 1.0_dp, &
                 "the lowering is far too large to be superposition error")

      call alone%destroy()
      call in_pair%destroy()
   end subroutine test_bsse_is_real

   subroutine test_far_ghost(error)
      !! Ghost functions a long way off change nothing
      !!
      !! The counterpart to the test above. Superposition error comes from
      !! functions near enough to be borrowed, so at two hundred Angstrom the
      !! ghosted monomer must return to its isolated energy. Without this, a
      !! ghost that was silently ignored and a ghost that worked would look the
      !! same from the sign of one difference.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: alone, far
      type(rhf_result_t) :: scf_alone, scf_far
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6)
      integer :: i

      call dimer_geometry(z, sym, c)
      call build_czt_molecule(z(1:3), sym(1:3), c(:, 1:3), "sto-3g", alone, err)
      call run_czt_rhf(alone, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_alone, err)
      call check(error,.not. err%has_error(), "alone: "//err%get_full_trace())
      if (allocated(error)) return

      do i = 4, 6
         c(1, i) = c(1, i) + 200.0_dp*ANG
      end do
      call build_czt_molecule(z, sym, c, "sto-3g", far, err, &
                              ghost=[.false., .false., .false., .true., .true., .true.])
      call run_czt_rhf(far, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_far, err)
      call check(error,.not. err%has_error(), "far: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, scf_far%energy, scf_alone%energy, &
                 "distant ghost functions changed the energy, which they cannot do", &
                 thr=1.0e-8_dp)

      call alone%destroy()
      call far%destroy()
   end subroutine test_far_ghost

   subroutine two_water_system(sys_geom)
      !! The dimer as a two-monomer system, three atoms each
      type(system_geometry_t), intent(out) :: sys_geom

      integer :: z(N_ATOMS)
      character(len=2) :: sym(N_ATOMS)
      real(dp) :: c(3, N_ATOMS)

      call dimer_geometry(z, sym, c)
      sys_geom%total_atoms = N_ATOMS
      sys_geom%n_monomers = 2
      sys_geom%atoms_per_monomer = N_MONOMER
      sys_geom%charge = 0
      sys_geom%multiplicity = 1
      allocate (sys_geom%element_numbers(N_ATOMS), source=z)
      allocate (sys_geom%coordinates(3, N_ATOMS), source=c)
   end subroutine two_water_system

   subroutine test_signed_indices(error)
      !! A negative monomer index contributes atoms and no electrons
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: pair, a_in_pair
      type(error_t) :: err

      call two_water_system(sys_geom)

      call build_fragment_from_indices(sys_geom, [1, 2], pair, err)
      call check(error,.not. err%has_error(), "pair: "//err%get_full_trace())
      if (allocated(error)) return

      ! Monomer A real, monomer B as ghosts: same atoms, half the electrons.
      call build_fragment_from_indices(sys_geom, [1, -2], a_in_pair, err)
      call check(error,.not. err%has_error(), "A in pair: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, a_in_pair%n_atoms, pair%n_atoms, &
                 "ghosting a monomer changed the atom count")
      if (allocated(error)) return
      call check(error, pair%nelec, 20, "the pair should have twenty electrons")
      if (allocated(error)) return
      call check(error, a_in_pair%nelec, 10, &
                 "one water ghosted should leave ten electrons")
      if (allocated(error)) return

      call check(error, allocated(a_in_pair%is_ghost), &
                 "the ghost mask was never set")
      if (allocated(error)) return
      call check(error, count(a_in_pair%is_ghost), N_MONOMER, &
                 "the wrong number of atoms was ghosted")
      if (allocated(error)) return
      ! Second monomer, so the last three atoms and not the first three.
      call check(error, all(a_in_pair%is_ghost(N_MONOMER + 1:)), &
                 "the ghosts landed on the wrong monomer")
      if (allocated(error)) return

      ! All-positive is the ordinary path, untouched.
      call check(error,.not. allocated(pair%is_ghost), &
                 "an unghosted fragment should carry no mask at all")
   end subroutine test_signed_indices

   subroutine test_ghost_core(error)
      !! A frozen core is counted over the real atoms, never over ghost centres
      !!
      !! A ghost keeps its element in `element_numbers`, for its basis. Counted
      !! there, water in the pair basis froze two orbitals rather than one --
      !! its own oxygen 1s and then its lowest valence orbital in the ghost
      !! oxygen's place -- which put the counterpoise monomer 52 mHartree above
      !! the isolated one at MP2/cc-pVDZ, and the VMFC dimer interaction at -75
      !! kcal/mol.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: pair, a_in_pair
      type(error_t) :: err

      call two_water_system(sys_geom)
      call build_fragment_from_indices(sys_geom, [1, 2], pair, err)
      call build_fragment_from_indices(sys_geom, [1, -2], a_in_pair, err)
      call check(error,.not. err%has_error(), "building: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, size(a_in_pair%real_element_numbers()), N_MONOMER, &
                 "the ghosted fragment should report three real atoms")
      if (allocated(error)) return
      call check(error, all(a_in_pair%real_element_numbers() == pair%element_numbers(1:N_MONOMER)), &
                 "the real atoms should be the first monomer's, in order")
      if (allocated(error)) return
      call check(error, core_orbital_count(a_in_pair%real_element_numbers()), 1, &
                 "water in the pair basis has one core orbital, not one per oxygen centre")
      if (allocated(error)) return
      call check(error, core_orbital_count(pair%real_element_numbers()), 2, &
                 "an unghosted fragment should count every atom")
   end subroutine test_ghost_core

   subroutine test_vmfc_identity(error)
      !! VMFC(2) on two fragments is the supermolecule, exactly
      !!
      !! At full expansion level there is nothing left to truncate, so
      !!
      !!     E_AB + (E_AB - E_A(b) - E_B(a)) - (E_AB - E_A(b) - E_B(a)) = E_AB
      !!
      !! trivially. The content is in the pieces: E_A(b) and E_B(a) must come
      !! from the *pair's* basis, and the interaction energy they give must be
      !! the counterpoise-corrected one rather than the raw one. Checked against
      !! the same quantity assembled by hand from explicit ghost masks, so a
      !! wrapper that silently dropped the ghosts would fail here rather than
      !! quietly returning the uncorrected number.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: pair, a_ghosted, b_ghosted
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf_pair, scf_a, scf_b
      type(error_t) :: err
      real(dp) :: e_int_cp

      call two_water_system(sys_geom)
      call build_fragment_from_indices(sys_geom, [1, 2], pair, err)
      call build_fragment_from_indices(sys_geom, [1, -2], a_ghosted, err)
      call build_fragment_from_indices(sys_geom, [-1, 2], b_ghosted, err)
      call check(error,.not. err%has_error(), "build: "//err%get_full_trace())
      if (allocated(error)) return

      call scf_of(pair, "6-31g", scf_pair, err)
      if (bail(error, err)) return
      call scf_of(a_ghosted, "6-31g", scf_a, err)
      if (bail(error, err)) return
      call scf_of(b_ghosted, "6-31g", scf_b, err)
      if (bail(error, err)) return

      e_int_cp = scf_pair%energy - scf_a%energy - scf_b%energy

      ! The dimer is symmetric, so the two monomer-in-pair-basis energies are
      ! the same calculation in two orientations and must agree.
      call check(error, scf_a%energy, scf_b%energy, &
                 "a symmetric dimer gave its two monomers different energies "// &
                 "in the pair basis", thr=1.0e-9_dp)
      if (allocated(error)) return

      ! And this is the counterpoise-corrected interaction energy, which SAPT
      ! reaches independently through the dimer-centred basis.
      call check(error, e_int_cp, SAPT_E_INT_HF, &
                 "the counterpoise-corrected interaction energy does not match "// &
                 "SAPT's counterpoise-corrected supermolecular HF", thr=1.0e-8_dp)

      call mol%destroy()
   end subroutine test_vmfc_identity

   subroutine scf_of(fragment, basis, scf, err)
      !! RHF on a fragment, ghosts and all, in the named basis
      type(physical_fragment_t), intent(in) :: fragment
      character(len=*), intent(in) :: basis
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
      character(len=2), allocatable :: sym(:)
      integer :: i

      allocate (sym(fragment%n_atoms))
      do i = 1, fragment%n_atoms
         if (fragment%element_numbers(i) == 8) then
            sym(i) = "O "
         else
            sym(i) = "H "
         end if
      end do

      call build_czt_molecule(fragment%element_numbers, sym, &
                              fragment%coordinates, basis, mol, err, &
                              ghost=ghost_mask(fragment))
      if (err%has_error()) return
      call run_czt_rhf(mol, fragment%nelec, 200, 1.0e-11_dp, 1.0e-9_dp, &
                       .false., scf, err)
      call mol%destroy()
   end subroutine scf_of

   function ghost_mask(fragment) result(mask)
      !! A fragment's ghost mask, all-false when it has none
      type(physical_fragment_t), intent(in) :: fragment
      logical :: mask(fragment%n_atoms)

      if (allocated(fragment%is_ghost)) then
         mask = fragment%is_ghost
      else
         mask = .false.
      end if
   end function ghost_mask

   function bail(error, err) result(failed)
      !! Turn an error_t into a test failure carrying its trace
      type(error_type), allocatable, intent(out) :: error
      type(error_t), intent(inout) :: err
      logical :: failed

      failed = err%has_error()
      if (failed) call check(error, .false., "SCF: "//err%get_full_trace())
   end function bail

   subroutine test_subset_key(error)
      !! The key names the chosen monomers real and the rest as ghosts
      !!
      !! This is the whole difference between MBE and VMFC in one function: the
      !! subset {A} of the pair {A,B} becomes {A in AB's basis}. If the
      !! complement were dropped rather than ghosted, the lookup would find the
      !! ordinary monomer and the expansion would be uncorrected -- which is a
      !! number, not a failure, so it is worth pinning here.
      type(error_type), allocatable, intent(out) :: error

      integer :: key(3)

      ! Pair [1,2], choosing the first: 1 real, 2 ghosted.
      call vmfc_subset_key([1, 2], 2, [1], 1, key(1:2))
      call check(error, key(1), 1, "the chosen monomer should stay real")
      if (allocated(error)) return
      call check(error, key(2), -2, "the complement should be ghosted")
      if (allocated(error)) return

      ! And the other way round.
      call vmfc_subset_key([1, 2], 2, [2], 1, key(1:2))
      call check(error, key(1), 2, "the chosen monomer should stay real")
      if (allocated(error)) return
      call check(error, key(2), -1, "the complement should be ghosted")
      if (allocated(error)) return

      ! Trimer [1,2,3] choosing two: both real, the third ghosted -- so a
      ! dimer-in-trimer-basis, which is what VMFC(3) subtracts.
      call vmfc_subset_key([1, 2, 3], 3, [1, 3], 2, key)
      call check(error, key(1), 1, "first chosen")
      if (allocated(error)) return
      call check(error, key(2), 3, "second chosen")
      if (allocated(error)) return
      call check(error, key(3), -2, "the unchosen monomer should be ghosted")
      if (allocated(error)) return

      ! Every monomer chosen is the parent itself, with nothing ghosted.
      call vmfc_subset_key([1, 2], 2, [1, 2], 2, key(1:2))
      call check(error, all(key(1:2) > 0), &
                 "choosing everything should ghost nothing")
   end subroutine test_subset_key

   subroutine test_row_subset_key(error)
      !! A row that is itself ghosted passes its ghosts on to its subsets
      !!
      !! `[1,2,-3]` is the pair 12 in the basis of 123. Its subsets are what
      !! Valiron-Mayer subtracts inside the trimer's correction, and they are in
      !! the trimer's basis too: `[1,-2,-3]`, not `[1,-2]`. Dropping the `-3` is
      !! the defect this key exists to prevent -- it looks up monomer 1 in the
      !! pair's basis, a real row with a real energy, so nothing fails and
      !! VMFC(3) comes out wrong.
      type(error_type), allocatable, intent(out) :: error

      integer :: key(3)
      integer :: key_len

      call vmfc_row_subset_key([1, 2, -3], [1], 1, key, key_len)
      call check(error, key_len, 3, "the key should span the whole trimer")
      if (allocated(error)) return
      call check(error, key(1), 1, "the chosen monomer should stay real")
      if (allocated(error)) return
      call check(error, all(key(2:3) == [-2, -3]), &
                 "the other real monomer should be ghosted and the row's ghost kept")
      if (allocated(error)) return

      ! Ghost first, and zero-padded: position counts real monomers only.
      call vmfc_row_subset_key([-1, 3, 2], [2], 1, key, key_len)
      call check(error, key_len, 3, "a leading ghost is still part of the row")
      if (allocated(error)) return
      call check(error, key(1), 2, "the second real monomer is 2, wherever the ghost sits")
      if (allocated(error)) return
      call check(error, all(key(2:3) == [-3, -1]), "3 ghosted, -1 kept")
      if (allocated(error)) return

      ! An unghosted row gives what vmfc_subset_key gives.
      call vmfc_row_subset_key([1, 2, 0], [2], 1, key, key_len)
      call check(error, key_len, 2, "padding is not a monomer")
      if (allocated(error)) return
      call check(error, all(key(1:2) == [2, -1]), "the plain pair rule")
   end subroutine test_row_subset_key

   subroutine test_auxiliary(error)
      !! A ghosted row is auxiliary, and its size is its real monomers
      !!
      !! Both rules exist to stop the same double count. VMFC's one-body term
      !! is each monomer in its *own* basis, so the ghosted rows belong inside
      !! the pair correction and nowhere else -- summing their deltas as well
      !! would count them twice. And [1,-2] has one real monomer, so it has no
      !! proper subsets; using the row width instead would send the recursion
      !! hunting for subsets that were never generated.
      type(error_type), allocatable, intent(out) :: error

      call check(error,.not. is_auxiliary_row([1, 2], COUNTERPOISE_VMFC), &
                 "an ordinary pair is not auxiliary")
      if (allocated(error)) return
      call check(error,.not. is_auxiliary_row([1, 0], COUNTERPOISE_VMFC), &
                 "a padded monomer row is not auxiliary")
      if (allocated(error)) return
      call check(error, is_auxiliary_row([1, -2], COUNTERPOISE_VMFC), &
                 "a ghosted row is auxiliary")
      if (allocated(error)) return
      call check(error, is_auxiliary_row([-1, 2], COUNTERPOISE_VMFC), &
                 "a ghosted row is auxiliary either way round")
      if (allocated(error)) return

      call check(error, int(real_count_of([1, 2])), 2, "a pair has two real")
      if (allocated(error)) return
      call check(error, int(real_count_of([1, -2])), 1, &
                 "a monomer in the pair basis has one real monomer")
      if (allocated(error)) return
      call check(error, int(real_count_of([1, 3, -2])), 2, &
                 "a dimer in the trimer basis has two real monomers")
      if (allocated(error)) return
      call check(error, int(real_count_of([1, 0])), 1, &
                 "padding is not a monomer")
   end subroutine test_auxiliary

   subroutine test_bsse_decays(error)
      !! Superposition error falls away as the monomers separate
      !!
      !! It exists because monomer A can reach B's basis functions, so it must
      !! shrink with the overlap that lets it. That makes it monotone in the
      !! separation, which is a stronger statement than "positive": a ghost
      !! implementation that leaked some constant error would still be positive
      !! everywhere and would not decay.
      !!
      !! At the far end the corrected and uncorrected interaction energies have
      !! to meet, which is the check that the correction invents nothing where
      !! there is no overlap left to correct.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N_R = 3
      real(dp), parameter :: R(N_R) = [3.0_dp, 4.0_dp, 6.0_dp]
      real(dp) :: bsse(N_R), raw, corrected
      integer :: k

      do k = 1, N_R
         call bsse_at(R(k), bsse(k), raw, corrected, error)
         if (allocated(error)) return
      end do

      call check(error, all(bsse > 0.0_dp), &
                 "the pair basis must lower a monomer at every separation")
      if (allocated(error)) return

      do k = 2, N_R
         call check(error, bsse(k) < bsse(k - 1), &
                    "superposition error must fall as the monomers separate")
         if (allocated(error)) return
      end do

      ! Two orders of magnitude across three Angstrom -- overlap, not a constant.
      call check(error, bsse(N_R) < bsse(1)/100.0_dp, &
                 "the decay is too slow to be an overlap effect")
      if (allocated(error)) return

      ! `raw` and `corrected` are the last row's, at the widest separation.
      call check(error, abs(raw - corrected) < 1.0e-7_dp, &
                 "with no overlap left, correcting must change nothing")
   end subroutine test_bsse_decays

   subroutine bsse_at(sep, bsse, raw, corrected, error)
      !! One separation: the BSSE on a monomer, and both interaction energies
      real(dp), intent(in) :: sep
      real(dp), intent(out) :: bsse, raw, corrected
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: alone, in_pair, dimer
      type(rhf_result_t) :: scf_alone, scf_pair, scf_dimer
      type(error_t) :: err
      integer :: z(N_ATOMS), i
      character(len=2) :: sym(N_ATOMS)
      real(dp) :: c(3, N_ATOMS)

      call dimer_geometry(z, sym, c)
      c = c/ANG
      do i = 1, N_MONOMER
         c(1, i + N_MONOMER) = c(1, i) + sep
      end do
      c = c*ANG

      call build_czt_molecule(z(1:N_MONOMER), sym(1:N_MONOMER), &
                              c(:, 1:N_MONOMER), "6-31g", alone, err)
      call run_czt_rhf(alone, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_alone, err)
      call build_czt_molecule(z, sym, c, "6-31g", in_pair, err, &
                              ghost=[.false., .false., .false., .true., .true., .true.])
      call run_czt_rhf(in_pair, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_pair, err)
      call build_czt_molecule(z, sym, c, "6-31g", dimer, err)
      call run_czt_rhf(dimer, 20, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_dimer, err)
      if (err%has_error()) then
         call check(error, .false., "SCF at "//to_char(sep)//": "//err%get_full_trace())
         return
      end if

      ! Symmetric dimer, so one monomer's BSSE is the other's.
      bsse = scf_alone%energy - scf_pair%energy
      raw = scf_dimer%energy - 2.0_dp*scf_alone%energy
      corrected = scf_dimer%energy - 2.0_dp*scf_pair%energy

      call alone%destroy()
      call in_pair%destroy()
      call dimer%destroy()
   end subroutine bsse_at

   !---------------------------------------------------------------------------
   ! Full-cluster-basis counterpoise
   !---------------------------------------------------------------------------

   subroutine water_cluster(n, radius, sys_geom)
      !! `n` waters on a circle of `radius` Angstrom, each tilted and turned
      !! differently so that no two monomers are equivalent
      !!
      !! Nothing here is symmetric under exchanging monomers, so a term looked
      !! up in the wrong basis, or a monomer mixed up with another, moves the
      !! energy rather than cancelling.
      integer, intent(in) :: n
      real(dp), intent(in) :: radius
      type(system_geometry_t), intent(out) :: sys_geom

      real(dp) :: template(3, N_MONOMER), rz(3, 3), rx(3, 3), centre(3)
      real(dp) :: theta, turn, tilt
      integer :: i, k, first

      template(:, 1) = [0.0_dp, 0.0_dp, 0.10077199_dp]
      template(:, 2) = [0.0_dp, 0.77250895_dp, -0.46780200_dp]
      template(:, 3) = [0.0_dp, -0.77250895_dp, -0.46780200_dp]

      sys_geom%total_atoms = N_MONOMER*n
      sys_geom%n_monomers = n
      sys_geom%atoms_per_monomer = N_MONOMER
      sys_geom%charge = 0
      sys_geom%multiplicity = 1
      allocate (sys_geom%element_numbers(N_MONOMER*n))
      allocate (sys_geom%coordinates(3, N_MONOMER*n))

      do i = 1, n
         theta = 2.0_dp*PI*real(i - 1, dp)/real(n, dp)
         turn = 0.9_dp*real(i, dp)
         tilt = 0.6_dp*real(i, dp)
         centre = [radius*cos(theta), radius*sin(theta), 0.35_dp*real(mod(i, 2), dp)]
         rz = reshape([cos(turn), sin(turn), 0.0_dp, &
                       -sin(turn), cos(turn), 0.0_dp, &
                       0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
         rx = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
                       0.0_dp, cos(tilt), sin(tilt), &
                       0.0_dp, -sin(tilt), cos(tilt)], [3, 3])
         first = N_MONOMER*(i - 1)
         do k = 1, N_MONOMER
            sys_geom%coordinates(:, first + k) = (centre + matmul(rz, matmul(rx, template(:, k))))*ANG
            sys_geom%element_numbers(first + k) = merge(8, 1, k == 1)
         end do
      end do
   end subroutine water_cluster

   subroutine chain_of_atoms(n, sys_geom)
      !! `n` one-atom monomers in a line, far apart enough that nothing is screened
      integer, intent(in) :: n
      type(system_geometry_t), intent(out) :: sys_geom

      integer :: i

      sys_geom%n_monomers = n
      sys_geom%atoms_per_monomer = 1
      sys_geom%total_atoms = n
      sys_geom%charge = 0
      sys_geom%multiplicity = 1
      allocate (sys_geom%element_numbers(n), source=10)
      allocate (sys_geom%coordinates(3, n), source=0.0_dp)
      do i = 1, n
         sys_geom%coordinates(1, i) = real(i - 1, dp)*4.0_dp
      end do
   end subroutine chain_of_atoms

   function synthetic_energy(row) result(e)
      !! An energy for a term-list row that depends on which monomers are real
      !! and which are ghosts, and on nothing else
      !!
      !! One-, two-, three- and four-body parts, none of them symmetric under
      !! exchanging monomers, and a different pair term for a real-ghost pair
      !! than for a real-real one. A row looked up in the wrong basis therefore
      !! returns a different number, and the order of the entries does not
      !! matter, as for a real energy.
      integer, intent(in) :: row(:)
      real(dp) :: e

      integer :: a, b, c, d
      real(dp) :: u, v, w, x

      e = 0.0_dp
      do a = 1, size(row)
         if (row(a) == 0) cycle
         u = real(abs(row(a)), dp)
         if (row(a) > 0) then
            e = e - 20.0_dp - 0.731_dp*u - 0.0113_dp*u*u
         else
            e = e - 1.9_dp - 0.0617_dp*u
         end if
         do b = a + 1, size(row)
            if (row(b) == 0) cycle
            v = real(abs(row(b)), dp)
            if (row(a) > 0 .and. row(b) > 0) then
               e = e - 0.0173_dp*u*v
            else if (row(a) < 0 .and. row(b) < 0) then
               e = e - 0.0009_dp*(u + v)
            else
               ! One real and one ghost: weigh them by role, not by position
               e = e - 0.0041_dp*(merge(u, v, row(a) > 0) + 2.0_dp*merge(v, u, row(a) > 0))
            end if
            do c = b + 1, size(row)
               if (row(c) == 0) cycle
               w = real(abs(row(c)), dp)
               e = e + 0.0007_dp*u*v*w*real(count([row(a), row(b), row(c)] > 0), dp)
               do d = c + 1, size(row)
                  if (row(d) == 0) cycle
                  x = real(abs(row(d)), dp)
                  e = e - 0.00031_dp*u*v*w*x*real(count([row(a), row(b), row(c), row(d)] > 0), dp)
               end do
            end do
         end do
      end do
   end function synthetic_energy

   function same_signed_set(x, y) result(equal)
      !! Whether two zero-padded rows name the same signed monomers
      integer, intent(in) :: x(:), y(:)
      logical :: equal

      integer :: i

      equal = count(x /= 0) == count(y /= 0)
      if (.not. equal) return
      do i = 1, size(x)
         if (x(i) == 0) cycle
         if (.not. any(y == x(i))) then
            equal = .false.
            return
         end if
      end do
   end function same_signed_set

   function find_row(polymers, count, wanted) result(index)
      !! The row of the list that names exactly the signed monomers of `wanted`,
      !! or 0 when there is none
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: count
      integer, intent(in) :: wanted(:)
      integer :: index

      integer(int64) :: i

      index = 0
      do i = 1_int64, count
         if (same_signed_set(polymers(i, :), wanted)) then
            index = int(i)
            return
         end if
      end do
   end function find_row

   function full_basis_row(real_monomers, n_monomers) result(row)
      !! The monomers in `real_monomers` real, every other monomer of the
      !! system ghosted; the whole system when nothing is left over
      integer, intent(in) :: real_monomers(:)
      integer, intent(in) :: n_monomers
      integer :: row(n_monomers)

      integer :: m, k

      row = 0
      k = size(real_monomers)
      row(1:k) = real_monomers
      do m = 1, n_monomers
         if (any(real_monomers == m)) cycle
         k = k + 1
         row(k) = -m
      end do
   end function full_basis_row

   subroutine test_ssfc_synthetic_identity(error)
      !! SSFC(4) over four monomers, from a generated list, is the Boys-Bernardi total
      !!
      !! The list is the one `generate_mbe_term_list` builds for
      !! `counterpoise = "ssfc"`, not a hand-built one, and the energies are
      !! made up (`synthetic_energy`), a different number for every real/ghost
      !! pattern. They are run through `compute_mbe`, which is the public face
      !! of the by-level recursion, with the scheme the deck would name. At
      !! L = N the sum must telescope to
      !!
      !!     E_cluster + sum_i [ E_i(i) - E_i(FB) ]
      !!
      !! with the three kinds of row read straight out of the generated list:
      !! the all-real whole system, the own-basis monomer `[i]`, and the
      !! full-basis monomer `[i,-others]`. A list that lacked any of them, or a
      !! recursion that looked a subset up in the wrong basis, misses it.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 4
      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      type(mbe_result_t) :: mbe_result
      type(calculation_result_t), allocatable :: results(:)
      integer, allocatable :: polymers(:, :)
      integer(int64) :: count, i
      integer :: m, row, console_level
      real(dp) :: expected

      call chain_of_atoms(N, sys_geom)
      config%nlevel = N
      config%counterpoise = "ssfc"
      call logger%configuration(level=console_level)
      call logger%configure(level=error_level)
      call generate_mbe_term_list(sys_geom, config, N, polymers, count)
      call logger%configure(level=console_level)

      ! Four own-basis monomers, four full-basis monomers, six pairs, four
      ! triples, and the whole system as an all-real row.
      call check(error, count, 19_int64, "SSFC(4) over four monomers should have nineteen rows")
      if (allocated(error)) return

      allocate (results(count))
      do i = 1_int64, count
         results(i)%has_energy = .true.
         results(i)%energy%scf = synthetic_energy(polymers(i, :))
      end do

      row = find_row(polymers, count, [(m, m=1, N)])
      call check(error, row > 0, "the list has no all-real whole-system row")
      if (allocated(error)) return
      expected = results(row)%energy%scf

      do m = 1, N
         row = find_row(polymers, count, [m])
         call check(error, row > 0, "the list has no own-basis row for monomer "//to_char(m))
         if (allocated(error)) return
         expected = expected + results(row)%energy%scf

         row = find_row(polymers, count, full_basis_row([m], N))
         call check(error, row > 0, "the list has no full-basis row for monomer "//to_char(m))
         if (allocated(error)) return
         expected = expected - results(row)%energy%scf
      end do

      call logger%configure(level=error_level)
      call compute_mbe(polymers, count, N, results, mbe_result, sys_geom, &
                       counterpoise_scheme=counterpoise_scheme_of("ssfc"))
      call logger%configure(level=console_level)

      call check(error, mbe_result%total_energy, expected, thr=1.0e-12_dp, &
                 message="SSFC(4) of a generated list is not E_cluster + sum_i [E_i(i) - E_i(FB)]")
      if (allocated(error)) return

      ! The made-up energies must actually separate the schemes' ingredients,
      ! or an identity between equal numbers would hold for any wiring.
      call check(error, abs(expected - synthetic_energy([(m, m=1, N)])) > 1.0e-3_dp, &
                 "the made-up energies carry no superposition error, so the identity proves nothing")
      call mbe_result%destroy()
   end subroutine test_ssfc_synthetic_identity

   subroutine rhf_of_row(sys_geom, row, energy, err)
      !! RHF/STO-3G of one term-list row, built from explicit ghost masks and
      !! run through the SCF directly, none of the expansion's machinery
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: row(:)
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err

      type(physical_fragment_t) :: fragment
      type(rhf_result_t) :: scf

      energy = 0.0_dp
      call build_fragment_from_indices(sys_geom, row, fragment, err)
      if (err%has_error()) return

      call scf_of(fragment, "sto-3g", scf, err)
      if (.not. err%has_error()) energy = scf%energy
      call fragment%destroy()
   end subroutine rhf_of_row

   subroutine expansion_of(sys_geom, level, scheme, json_data, reference)
      !! The real code path for a water cluster: the term list the driver would
      !! build, every fragment solved by RHF/STO-3G, and the expansion over them
      !!
      !! This is `generate_mbe_term_list` followed by `serial_fragment_processor`,
      !! which is what the driver runs on one rank, with the scheme named by the
      !! deck's spelling. Console logging is silenced for the run.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: level
      character(len=*), intent(in) :: scheme
         !! The deck's spelling: "none", "vmfc" or "ssfc"
      type(json_output_data_t), intent(out) :: json_data
      integer, intent(in), optional :: reference
         !! Monomer whose interaction energy is wanted, 1-based; absent for a total

      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: count
      integer :: console_level

      config%nlevel = level
      config%counterpoise = scheme
      config%method_config%method_type = METHOD_TYPE_HF
      config%method_config%basis_set = "sto-3g"
      config%method_config%scf%energy_convergence = 1.0e-11_dp
      config%method_config%scf%density_convergence = 1.0e-9_dp
      if (present(reference)) config%reference_fragment = reference

      call logger%configuration(level=console_level)
      call logger%configure(level=error_level)
      call generate_mbe_term_list(sys_geom, config, level, polymers, count)
      call serial_fragment_processor(count, polymers, level, sys_geom, config%method_config, &
                                     CALC_TYPE_ENERGY, json_data, &
                                     reference_fragment=config%reference_fragment, &
                                     counterpoise_scheme=counterpoise_scheme_of(scheme))
      call logger%configure(level=console_level)
   end subroutine expansion_of

   subroutine test_ssfc_trimer_identity(error)
      !! SSFC(3) of a water trimer, HF/STO-3G, is the Boys-Bernardi total
      !!
      !! One side is the whole expansion: the generated term list, ten RHF
      !! calculations, and the many-body recursion. The other is assembled by
      !! hand from seven single calculations -- the supermolecule, the three
      !! monomers in their own basis and the three monomers in the trimer's
      !! basis, the latter built with explicit ghost masks:
      !!
      !!     E_cluster + sum_i [ E_i(i) - E_i(FB) ]
      !!
      !! They must agree to 1e-10 Eh. The superposition error summed over the
      !! monomers is also required to be larger than a microhartree, so that the
      !! corrected total is visibly not the uncorrected one and a scheme that
      !! silently ghosted nothing could not satisfy the identity.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 3
      integer, parameter :: MONOMERS(N) = [1, 2, 3]
      type(system_geometry_t) :: sys_geom
      type(json_output_data_t) :: json_data
      type(error_t) :: err
      real(dp) :: e_cluster, e_own, e_full, bsse, by_hand
      integer :: m
      integer :: others(2)

      call water_cluster(N, 1.75_dp, sys_geom)

      call expansion_of(sys_geom, N, "ssfc", json_data)
      call check(error, json_data%has_energy, "the SSFC(3) run reported no total energy")
      if (allocated(error)) return
      ! Three own-basis monomers, three full-basis monomers, three full-basis
      ! pairs, and the whole system as an all-real row.
      call check(error, json_data%fragment_count, 10_int64, &
                 "SSFC(3) over three monomers should be ten rows")
      if (allocated(error)) return

      call rhf_of_row(sys_geom, [1, 2, 3], e_cluster, err)
      if (bail(error, err)) return
      bsse = 0.0_dp
      do m = 1, N
         others = pack(MONOMERS, MONOMERS /= m)
         call rhf_of_row(sys_geom, [m], e_own, err)
         if (bail(error, err)) return
         call rhf_of_row(sys_geom, [m, -others(1), -others(2)], e_full, err)
         if (bail(error, err)) return
         bsse = bsse + (e_own - e_full)
      end do
      by_hand = e_cluster + bsse

      call check(error, json_data%total_energy, by_hand, thr=1.0e-10_dp, &
                 message="SSFC(3) of a water trimer is not E_cluster + sum_i [E_i(i) - E_i(FB)]")
      if (allocated(error)) return

      call check(error, bsse > 1.0e-6_dp, &
                 "the monomers gain nothing from the trimer's basis, so the identity is vacuous")
   end subroutine test_ssfc_trimer_identity

   subroutine test_ssfc_dimer_is_vmfc(error)
      !! With two monomers the two schemes are one: the pair is the cluster
      !!
      !! Both lists are `[1] [2] [1,-2] [-1,2] [1,2]`, so a real water dimer
      !! under SSFC(2) and VMFC(2) must give the same energy to 1e-12 Eh, and
      !! that energy is the counterpoise-corrected dimer assembled by hand.
      !! Asserted so that the coincidence is not "fixed" later by making SSFC
      !! different where it should not be.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(json_output_data_t) :: by_ssfc, by_vmfc
      type(error_t) :: err
      real(dp) :: e_a, e_b, e_ab, e_a_in_ab, e_b_in_ab, by_hand

      call water_cluster(2, 1.5_dp, sys_geom)

      call expansion_of(sys_geom, 2, "ssfc", by_ssfc)
      call expansion_of(sys_geom, 2, "vmfc", by_vmfc)

      call check(error, by_ssfc%fragment_count, 5_int64, "SSFC(2) over two monomers has five rows")
      if (allocated(error)) return
      call check(error, by_vmfc%fragment_count, 5_int64, "VMFC(2) over two monomers has five rows")
      if (allocated(error)) return

      call check(error, by_ssfc%total_energy, by_vmfc%total_energy, thr=1.0e-12_dp, &
                 message="SSFC(2) and VMFC(2) differ on a water dimer")
      if (allocated(error)) return

      call rhf_of_row(sys_geom, [1], e_a, err)
      if (bail(error, err)) return
      call rhf_of_row(sys_geom, [2], e_b, err)
      if (bail(error, err)) return
      call rhf_of_row(sys_geom, [1, 2], e_ab, err)
      if (bail(error, err)) return
      call rhf_of_row(sys_geom, [1, -2], e_a_in_ab, err)
      if (bail(error, err)) return
      call rhf_of_row(sys_geom, [-1, 2], e_b_in_ab, err)
      if (bail(error, err)) return
      by_hand = e_a + e_b + (e_ab - e_a_in_ab - e_b_in_ab)

      call check(error, by_ssfc%total_energy, by_hand, thr=1.0e-10_dp, &
                 message="SSFC(2) of a water dimer is not the counterpoise-corrected dimer")
   end subroutine test_ssfc_dimer_is_vmfc

   subroutine test_ssfc_is_not_vmfc(error)
      !! On three monomers SSFC(2) and VMFC(2) are different numbers
      !!
      !! VMFC puts a pair in the pair's basis and SSFC in the trimer's, so the
      !! third water's functions reach the pair term under SSFC and not under
      !! VMFC. A run in which SSFC fell back to VMFC rows would agree with VMFC
      !! to the last digit; they must differ by well above numerical noise. The
      !! SSFC(2) total is also checked against the same expansion assembled by
      !! hand, so the difference is shown to come from SSFC being right rather
      !! than from either being wrong:
      !!
      !!     sum_i E_i(i) + sum_{i<j} [ E_ij(FB) - E_i(FB) - E_j(FB) ]
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 3
      type(system_geometry_t) :: sys_geom
      type(json_output_data_t) :: by_ssfc, by_vmfc
      type(error_t) :: err
      real(dp) :: by_hand, e, e_pair_fb
      real(dp) :: e_fb(N)
      integer :: i, j

      call water_cluster(N, 1.75_dp, sys_geom)

      call expansion_of(sys_geom, 2, "ssfc", by_ssfc)
      call expansion_of(sys_geom, 2, "vmfc", by_vmfc)

      by_hand = 0.0_dp
      do i = 1, N
         call rhf_of_row(sys_geom, [i], e, err)
         if (bail(error, err)) return
         by_hand = by_hand + e
         call rhf_of_row(sys_geom, full_basis_row([i], N), e_fb(i), err)
         if (bail(error, err)) return
      end do
      do i = 1, N - 1
         do j = i + 1, N
            call rhf_of_row(sys_geom, full_basis_row([i, j], N), e_pair_fb, err)
            if (bail(error, err)) return
            by_hand = by_hand + e_pair_fb - e_fb(i) - e_fb(j)
         end do
      end do

      call check(error, by_ssfc%total_energy, by_hand, thr=1.0e-10_dp, &
                 message="SSFC(2) of a water trimer is not the full-basis inclusion-exclusion total")
      if (allocated(error)) return

      ! 0.11 mEh at this geometry; the SCF noise is 1e-12, so anything above ten
      ! microhartree is the schemes differing.
      call check(error, abs(by_ssfc%total_energy - by_vmfc%total_energy) > 1.0e-5_dp, &
                 "SSFC(2) and VMFC(2) agree on a water trimer, so SSFC is running as VMFC: "// &
                 to_char(by_ssfc%total_energy - by_vmfc%total_energy))
   end subroutine test_ssfc_is_not_vmfc

   subroutine test_ssfc_interaction_energy(error)
      !! An InteractionEnergy run under SSFC reports the reference's full-basis terms
      !!
      !! A water tetramer at level 3 is run twice: as an ordinary expansion, and
      !! as an interaction energy with monomer 2 as the reference, which
      !! reduces the term list to the terms holding the reference and their
      !! subsets. The reference's share, level by level, must equal the sum of
      !! the deltas of the same rows in the ordinary run, read from its fragment
      !! table; its one-body energy must be the own-basis monomer's, and not the
      !! full-basis one the same list also holds. The level-two share is
      !! assembled once more from the raw full-basis energies of that table, so
      !! the deltas are not taken on trust either.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 4, LEVEL = 3, REFERENCE = 2
      type(system_geometry_t) :: sys_geom
      type(json_output_data_t) :: ordinary, reduced
      integer :: row, n_real, level_of, j
      integer(int64) :: i
      integer(int64) :: terms(LEVEL)
      real(dp) :: pair_share, e_own, e_fb_ref
      real(dp) :: share(LEVEL)

      call water_cluster(N, 2.1_dp, sys_geom)

      call expansion_of(sys_geom, LEVEL, "ssfc", ordinary)
      call expansion_of(sys_geom, LEVEL, "ssfc", reduced, reference=REFERENCE)

      call check(error, reduced%has_interaction, "a reference run reported no interaction energy")
      if (allocated(error)) return
      call check(error,.not. reduced%has_energy, &
                 "an interaction energy run must not also report a total")
      if (allocated(error)) return
      call check(error, reduced%fragment_count < ordinary%fragment_count, &
                 "the reference run kept every term, so nothing was reduced")
      if (allocated(error)) return

      ! The reference's terms in the ordinary run: the rows SSFC sums, with the
      ! reference among their real monomers.
      share = 0.0_dp
      terms = 0_int64
      do i = 1_int64, ordinary%fragment_count
         if (is_auxiliary_row(ordinary%polymers(i, :), COUNTERPOISE_SSFC)) cycle
         if (.not. any(ordinary%polymers(i, :) == REFERENCE)) cycle
         n_real = int(real_count_of(ordinary%polymers(i, :)))
         if (n_real < 2) cycle
         share(n_real) = share(n_real) + ordinary%delta_energies(i)
         terms(n_real) = terms(n_real) + 1_int64
      end do

      ! With monomer 2 of four at level 3, three pairs and three triples hold it.
      ! Fixed here rather than counted by the predicate the expansion also uses,
      ! which would zero both sides together if it were wrong.
      call check(error, terms(2), 3_int64, "three pairs of four monomers hold the reference")
      if (allocated(error)) return
      call check(error, terms(3), 3_int64, "three triples of four monomers hold the reference")
      if (allocated(error)) return

      do level_of = 2, LEVEL
         call check(error, reduced%interaction_count_by_level(level_of), terms(level_of), &
                    "wrong number of terms holding the reference at level "//to_char(level_of))
         if (allocated(error)) return
         call check(error, reduced%interaction_by_level(level_of), share(level_of), thr=1.0e-10_dp, &
                    message="the reference's share at level "//to_char(level_of)// &
                    " is not the sum of its terms' full-basis deltas")
         if (allocated(error)) return
      end do
      call check(error, reduced%interaction_energy, sum(share), thr=1.0e-10_dp, &
                 message="the interaction energy is not the sum of the reference's terms")
      if (allocated(error)) return

      ! One-body energy: own basis.
      row = find_row(ordinary%polymers, ordinary%fragment_count, [REFERENCE])
      call check(error, row > 0, "the ordinary list has no own-basis row for the reference")
      if (allocated(error)) return
      e_own = ordinary%fragment_energies(row)
      row = find_row(ordinary%polymers, ordinary%fragment_count, full_basis_row([REFERENCE], N))
      call check(error, row > 0, "the ordinary list has no full-basis row for the reference")
      if (allocated(error)) return
      e_fb_ref = ordinary%fragment_energies(row)
      call check(error, reduced%reference_energy, e_own, thr=1.0e-10_dp, &
                 message="the reference energy is not the own-basis monomer's")
      if (allocated(error)) return
      call check(error, abs(e_own - e_fb_ref) > 1.0e-6_dp, &
                 "own-basis and full-basis reference energies coincide, so the choice is untested")
      if (allocated(error)) return

      ! Level two from raw energies: E_ij(FB) - E_i(FB) - E_j(FB) for each partner.
      pair_share = 0.0_dp
      do j = 1, N
         if (j == REFERENCE) cycle
         row = find_row(ordinary%polymers, ordinary%fragment_count, &
                        full_basis_row([min(REFERENCE, j), max(REFERENCE, j)], N))
         call check(error, row > 0, "the ordinary list lacks a full-basis pair row")
         if (allocated(error)) return
         pair_share = pair_share + ordinary%fragment_energies(row) - e_fb_ref
         row = find_row(ordinary%polymers, ordinary%fragment_count, full_basis_row([j], N))
         call check(error, row > 0, "the ordinary list lacks a full-basis monomer row")
         if (allocated(error)) return
         pair_share = pair_share - ordinary%fragment_energies(row)
      end do
      call check(error, reduced%interaction_by_level(2), pair_share, thr=1.0e-10_dp, &
                 message="the level-two share is not built from full-basis energies")
   end subroutine test_ssfc_interaction_energy

end module test_mqc_counterpoise

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_counterpoise, only: collect_counterpoise
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_counterpoise", collect_counterpoise)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
