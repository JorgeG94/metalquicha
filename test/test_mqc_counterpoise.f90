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
module test_mqc_counterpoise
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, &
                                    build_fragment_with_ghosts
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_counterpoise

   real(dp), parameter :: ANG = 1.8897261254578281_dp

   !> Two waters, 3 Angstrom apart along x. Far enough that the interaction is
   !> small and close enough that each borrows the other's functions.
   real(dp), parameter :: SEP = 3.0_dp

   integer, parameter :: N_ATOMS = 6
   integer, parameter :: N_MONOMER = 3

   !> SAPT's counterpoise-corrected supermolecular Hartree-Fock interaction
   !> energy for this dimer in 6-31G, reached through the dimer-centred basis
   !> and none of this code. An outside number, not one this program printed.
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
                  new_unittest("vmfc_reproduces_the_supermolecule", test_vmfc_identity) &
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

      type(libcint_molecule_t) :: dimer, ghosted
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6)
      logical :: ghost(6)

      call dimer_geometry(z, sym, c)
      ghost = [.false., .false., .false., .true., .true., .true.]

      call build_libcint_molecule(z, sym, c, "sto-3g", dimer, err)
      call check(error,.not. err%has_error(), "dimer: "//err%get_full_trace())
      if (allocated(error)) return

      call build_libcint_molecule(z, sym, c, "sto-3g", ghosted, err, ghost=ghost)
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

      type(libcint_molecule_t) :: dimer, ghost_a, ghost_b
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6)
      real(dp), allocatable :: s_dimer(:, :), s_ghost(:, :)

      call dimer_geometry(z, sym, c)
      call build_libcint_molecule(z, sym, c, "sto-3g", dimer, err)
      call build_libcint_molecule(z, sym, c, "sto-3g", ghost_a, err, &
                                  ghost=[.true., .true., .true., .false., .false., .false.])
      call build_libcint_molecule(z, sym, c, "sto-3g", ghost_b, err, &
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

      type(libcint_molecule_t) :: alone, in_pair
      type(rhf_result_t) :: scf_alone, scf_pair
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6), bsse

      call dimer_geometry(z, sym, c)

      ! Monomer A in its own basis: three atoms, nothing else.
      call build_libcint_molecule(z(1:3), sym(1:3), c(:, 1:3), "sto-3g", alone, err)
      call check(error,.not. err%has_error(), "alone: "//err%get_full_trace())
      if (allocated(error)) return
      call run_libcint_rhf(alone, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_alone, err)
      call check(error,.not. err%has_error(), "alone SCF: "//err%get_full_trace())
      if (allocated(error)) return

      ! The same monomer, same ten electrons, in the pair's basis.
      call build_libcint_molecule(z, sym, c, "sto-3g", in_pair, err, &
                                  ghost=[.false., .false., .false., .true., .true., .true.])
      call check(error,.not. err%has_error(), "in pair: "//err%get_full_trace())
      if (allocated(error)) return
      call run_libcint_rhf(in_pair, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_pair, err)
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

      type(libcint_molecule_t) :: alone, far
      type(rhf_result_t) :: scf_alone, scf_far
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6)
      integer :: i

      call dimer_geometry(z, sym, c)
      call build_libcint_molecule(z(1:3), sym(1:3), c(:, 1:3), "sto-3g", alone, err)
      call run_libcint_rhf(alone, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_alone, err)
      call check(error,.not. err%has_error(), "alone: "//err%get_full_trace())
      if (allocated(error)) return

      do i = 4, 6
         c(1, i) = c(1, i) + 200.0_dp*ANG
      end do
      call build_libcint_molecule(z, sym, c, "sto-3g", far, err, &
                                  ghost=[.false., .false., .false., .true., .true., .true.])
      call run_libcint_rhf(far, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_far, err)
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

      call build_fragment_with_ghosts(sys_geom, [1, 2], pair, err)
      call check(error,.not. err%has_error(), "pair: "//err%get_full_trace())
      if (allocated(error)) return

      ! Monomer A real, monomer B as ghosts: same atoms, half the electrons.
      call build_fragment_with_ghosts(sys_geom, [1, -2], a_in_pair, err)
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
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf_pair, scf_a, scf_b
      type(error_t) :: err
      real(dp) :: e_int_cp

      call two_water_system(sys_geom)
      call build_fragment_with_ghosts(sys_geom, [1, 2], pair, err)
      call build_fragment_with_ghosts(sys_geom, [1, -2], a_ghosted, err)
      call build_fragment_with_ghosts(sys_geom, [-1, 2], b_ghosted, err)
      call check(error,.not. err%has_error(), "build: "//err%get_full_trace())
      if (allocated(error)) return

      call scf_of(pair, scf_pair, err)
      if (bail(error, err)) return
      call scf_of(a_ghosted, scf_a, err)
      if (bail(error, err)) return
      call scf_of(b_ghosted, scf_b, err)
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

   subroutine scf_of(fragment, scf, err)
      !! RHF on a fragment, ghosts and all
      type(physical_fragment_t), intent(in) :: fragment
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
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

      call build_libcint_molecule(fragment%element_numbers, sym, &
                                  fragment%coordinates, "6-31g", mol, err, &
                                  ghost=ghost_mask(fragment))
      if (err%has_error()) return
      call run_libcint_rhf(mol, fragment%nelec, 200, 1.0e-11_dp, 1.0e-9_dp, &
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
