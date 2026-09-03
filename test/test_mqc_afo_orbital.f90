!! Unit tests for the frozen orbital taken off a model system
module test_mqc_afo_orbital
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: system_geometry_t, to_bohr, to_angstrom
   use mqc_bond_perception, only: find_severed_bonds, severed_bond_t
   use mqc_czt_afo, only: afo_model_t, afo_options_t, build_afo_model, &
                          bond_hybrid, BOND_ORBITAL_REACH, &
                          afo_hybrid_t, build_group_frozen
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, atom_ao_blocks
   use mqc_czt_rhf, only: run_czt_rhf, rhf_result_t
   use mqc_fock_projector, only: fock_projector_t, build_frozen_basis
   implicit none

   !! The model SCF converges to 1e-10 and the localization to its own sweep
   !! tolerance, so an orbital coefficient is good to somewhere around here.
   real(dp), parameter :: TOL = 1.0e-6_dp

   !! STO-3G on carbon is 1s, 2s, then three 2p. The tests naming these indices
   !! name the basis too, so the layout is a statement about that basis rather
   !! than an assumption about all of them.
   integer, parameter :: P_FIRST = 3, P_LAST = 5

   private
   public :: collect_mqc_afo_orbital

contains

   subroutine collect_mqc_afo_orbital(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("a_single_sigma_bond_carries_one_orbital", test_one_orbital), &
                  new_unittest("cores_sit_at_half_the_bond_and_are_excluded", test_cores), &
                  new_unittest("hybrid_is_normalised_on_its_own_atom", test_normalised), &
                  new_unittest("hybrid_points_along_the_bond", test_points), &
                  new_unittest("hybrid_rotates_with_the_system", test_rotates), &
                  new_unittest("frozen_columns_land_on_their_own_atom", test_place), &
                  new_unittest("frozen_puts_the_occupied_columns_first", test_order), &
                  new_unittest("frozen_refuses_a_hybrid_from_another_basis", test_wrong_basis), &
                  new_unittest("a_frozen_hybrid_comes_out_of_the_scf_empty", test_end_to_end) &
                  ]
   end subroutine collect_mqc_afo_orbital

   subroutine hybrid_of(spin, hybrid, distance, model, error, err)
      !! Ethane, optionally rotated, cut at the C-C, through to the hybrid
      real(dp), intent(in) :: spin(3, 3)
      real(dp), allocatable, intent(out) :: hybrid(:), distance(:)
      type(afo_model_t), intent(out) :: model
      type(error_type), allocatable, intent(out) :: error
      type(error_t), intent(inout) :: err

      type(system_geometry_t) :: sys
      type(severed_bond_t), allocatable :: cuts(:)
      type(afo_options_t) :: opts
      integer :: n_cuts, n_on_bond

      call ethane(sys)
      sys%coordinates = matmul(spin, sys%coordinates)
      call find_severed_bonds(sys, [1, 1, 1, 1, 2, 2, 2, 2], cuts, n_cuts)
      call build_afo_model(sys%element_numbers, sys%coordinates, cuts(1), model, err, &
                           radius=5.0_dp)
      opts%basis = "sto-3g"
      call bond_hybrid(model, opts, hybrid, n_on_bond, err, centroid_distance=distance)
      call check(error,.not. err%has_error(), "taking the hybrid off the model failed")
   end subroutine hybrid_of

   pure function identity() result(eye)
      real(dp) :: eye(3, 3)
      eye = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                     0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
   end function identity

   subroutine test_one_orbital(error)
      !! A C-C single bond has exactly one localized orbital sitting on it
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(afo_model_t) :: model
      real(dp), allocatable :: hybrid(:), distance(:)
      real(dp) :: bond_length
      integer :: n_on

      call hybrid_of(identity(), hybrid, distance, model, error, err)
      if (allocated(error)) return

      bond_length = sqrt(sum((model%xyz(:, model%bda_local) &
                              - model%xyz(:, model%baa_local))**2))
      n_on = count(distance < BOND_ORBITAL_REACH*bond_length)
      call check(error, n_on, 1, &
                 "a single sigma bond did not report exactly one orbital on it")
      if (allocated(error)) return
      call check(error, minval(distance) < 1.0e-6_dp, &
                 "the bond orbital's centroid is not at the midpoint of a symmetric bond")
   end subroutine test_one_orbital

   subroutine test_cores(error)
      !! Why the reach must be well under a half
      !!
      !! A core orbital's centroid is on its nucleus, which is exactly half a
      !! bond length from the midpoint -- so a reach of a half admits both
      !! cores and a single bond reports three orbitals, which reads as a triple
      !! bond. That is structural and not a property of ethane, so it is worth a
      !! test rather than a comment: the cores are found where the geometry says
      !! they must be, and the reach is on the right side of them.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(afo_model_t) :: model
      real(dp), allocatable :: hybrid(:), distance(:)
      real(dp) :: bond_length, half
      integer :: n_at_half

      call hybrid_of(identity(), hybrid, distance, model, error, err)
      if (allocated(error)) return

      bond_length = sqrt(sum((model%xyz(:, model%bda_local) &
                              - model%xyz(:, model%baa_local))**2))
      half = 0.5_dp*bond_length

      n_at_half = count(abs(distance - half) < 1.0e-4_dp)
      call check(error, n_at_half, 2, &
                 "the two carbon cores are not sitting on their nuclei")
      if (allocated(error)) return
      call check(error, BOND_ORBITAL_REACH < 0.45_dp, &
                 "the reach is close enough to a half to start counting cores as bonds")
   end subroutine test_cores

   subroutine test_normalised(error)
      !! `h^T S_AA h = 1` in the bond-detached atom's own block
      !!
      !! Checked against an overlap this test builds for itself, so it is the
      !! property that is asserted and not the line of code that produced it.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(afo_model_t) :: model
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: hybrid(:), distance(:), s(:, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: first, last
      real(dp) :: norm

      call hybrid_of(identity(), hybrid, distance, model, error, err)
      if (allocated(error)) return

      call build_czt_molecule(model%z, model%sym, model%xyz, "sto-3g", mol, err)
      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      first = offsets(model%bda_local) + 1
      last = first + counts(model%bda_local) - 1
      call mol%overlap(s)

      call check(error, size(hybrid), counts(model%bda_local), &
                 "the hybrid is not the size of the atom's basis")
      if (allocated(error)) return
      norm = dot_product(hybrid, matmul(s(first:last, first:last), hybrid))
      call check(error, abs(norm - 1.0_dp) < TOL, &
                 "the hybrid is not normalised on the atom it belongs to")
   end subroutine test_normalised

   subroutine test_points(error)
      !! The hybrid points at the atom on the other end of the cut bond
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(afo_model_t) :: model
      real(dp), allocatable :: hybrid(:), distance(:)
      real(dp) :: p(3), axis(3), along

      call hybrid_of(identity(), hybrid, distance, model, error, err)
      if (allocated(error)) return
      call check(error, size(hybrid), 5, "sto-3g carbon should have five functions")
      if (allocated(error)) return

      p = hybrid(P_FIRST:P_LAST)
      p = p/sqrt(sum(p**2))
      axis = model%xyz(:, model%baa_local) - model%xyz(:, model%bda_local)
      axis = axis/sqrt(sum(axis**2))

      ! Up to sign: an orbital's overall phase is arbitrary, its axis is not.
      along = abs(dot_product(p, axis))
      call check(error, along > 1.0_dp - 1.0e-4_dp, &
                 "the hybrid's p component is not along the bond it stands on")
   end subroutine test_points

   subroutine test_rotates(error)
      !! Rotate the system and the hybrid rotates with it
      !!
      !! The test the whole construction has to pass and the one that needs no
      !! external reference. An `s` function is invariant under rotation and a
      !! `p` shell transforms as a vector, so the same orbital seen from a
      !! turned frame has the same `s` coefficients and `p` coefficients turned
      !! by the same rotation. An orbital's overall phase is arbitrary, so the
      !! comparison is up to one global sign, fixed here from the `p` block.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(afo_model_t) :: still, spun
      real(dp), allocatable :: h0(:), h1(:), d0(:), d1(:)
      real(dp) :: rot(3, 3), turn(3, 3), c, s, expected(3), phase

      c = cos(0.7_dp)
      s = sin(0.7_dp)
      rot = reshape([c, -s, 0.0_dp, s, c, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      c = cos(0.4_dp)
      s = sin(0.4_dp)
      turn = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, c, -s, 0.0_dp, s, c], [3, 3])
      rot = matmul(turn, rot)

      call hybrid_of(identity(), h0, d0, still, error, err)
      if (allocated(error)) return
      call hybrid_of(rot, h1, d1, spun, error, err)
      if (allocated(error)) return

      expected = matmul(rot, h0(P_FIRST:P_LAST))
      phase = 1.0_dp
      if (dot_product(expected, h1(P_FIRST:P_LAST)) < 0.0_dp) phase = -1.0_dp

      call check(error, maxval(abs(h1(P_FIRST:P_LAST) - phase*expected)) < TOL, &
                 "the hybrid's p block did not rotate with the molecule")
      if (allocated(error)) return
      call check(error, maxval(abs(h1(:P_FIRST - 1) - phase*h0(:P_FIRST - 1))) < TOL, &
                 "the hybrid's s coefficients changed under a rotation")
   end subroutine test_rotates

   subroutine test_end_to_end(error)
      !! Model system to hybrid to constrained SCF, and the orbital is empty
      !!
      !! Everything built for AFO, composed: solve a model, take the orbital on
      !! the cut bond, place it in another molecule's basis, orthonormalise it
      !! into a frozen basis, constrain a Fock matrix with it and solve. The
      !! assertion is physical rather than structural -- an orbital frozen as
      !! virtual has to come back with no electrons in it.
      !!
      !! Ethane stands in for a fragment here, so the hybrid is frozen in the
      !! molecule it came from. That makes the check sharp: the C-C hybrid is a
      !! large part of an occupied bond, so it is well populated unless the
      !! constraint actually removed it. The unconstrained population is
      !! measured in the same test rather than assumed, so the comparison is
      !! against this molecule and not against a remembered number.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(afo_model_t) :: model
      type(czt_molecule_t) :: mol
      type(afo_hybrid_t) :: hybrids(1)
      type(fock_projector_t) :: proj
      type(rhf_result_t) :: bare, held
      real(dp), allocatable :: h0(:), d0(:), frozen(:, :), basis(:, :), s(:, :), sh(:)
      integer :: n_occ, n_mo

      call hybrid_of(identity(), h0, d0, model, error, err)
      if (allocated(error)) return

      call ethane_molecule(mol, err)
      hybrids(1)%coeff = h0
      call build_group_frozen(mol, [model%bda_local], [.false.], hybrids, frozen, n_occ, err)
      call check(error,.not. err%has_error(), "placing the hybrid failed")
      if (allocated(error)) return

      call mol%overlap(s)
      call build_frozen_basis(frozen, 0, s, basis, n_mo, err)
      call check(error,.not. err%has_error(), "building the frozen basis failed")
      if (allocated(error)) return

      call proj%init(basis, s, 0, 1, 1.0e3_dp, err)
      call check(error,.not. err%has_error(), "the projector rejected the basis")
      if (allocated(error)) return

      call run_czt_rhf(mol, 18, 100, 1.0e-10_dp, 1.0e-8_dp, .false., bare, err)
      call run_czt_rhf(mol, 18, 100, 1.0e-10_dp, 1.0e-8_dp, .false., held, err, &
                       projector=proj)
      call check(error,.not. err%has_error(), "the constrained SCF failed")
      if (allocated(error)) return
      call check(error, held%converged, "the constrained SCF did not converge")
      if (allocated(error)) return

      ! `n = h^T S D S h` is the number of electrons in `h`, which is two for an
      ! occupied orbital and zero for an empty one.
      sh = matmul(s, frozen(:, 1))
      call check(error, dot_product(sh, matmul(bare%density, sh)) > 1.0_dp, &
                 "the hybrid is not populated even without the constraint, so this "// &
                 "test would pass for the wrong reason")
      if (allocated(error)) return
      call check(error, dot_product(sh, matmul(held%density, sh)) < 1.0e-8_dp, &
                 "the frozen virtual still holds electrons")
   end subroutine test_end_to_end

   subroutine ethane_molecule(mol, err)
      !! Ethane in STO-3G: five functions on each carbon, one on each hydrogen
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      type(system_geometry_t) :: sys
      character(len=2) :: sym(8)
      integer :: i

      call ethane(sys)
      do i = 1, 8
         if (sys%element_numbers(i) == 6) then
            sym(i) = "C "
         else
            sym(i) = "H "
         end if
      end do
      call build_czt_molecule(sys%element_numbers, sym, sys%coordinates, "sto-3g", &
                              mol, err)
   end subroutine ethane_molecule

   subroutine test_place(error)
      !! A hybrid writes into its own atom's block and nowhere else
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(afo_hybrid_t) :: hybrids(1)
      real(dp), allocatable :: frozen(:, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_occ, first, last

      call ethane_molecule(mol, err)
      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      ! The second carbon, so a wrong offset shows up as a wrong block rather
      ! than as the first one by luck.
      hybrids(1)%coeff = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp]
      call build_group_frozen(mol, [5], [.false.], hybrids, frozen, n_occ, err)
      call check(error,.not. err%has_error(), "placing the hybrid failed")
      if (allocated(error)) return

      first = offsets(5) + 1
      last = first + counts(5) - 1
      call check(error, maxval(abs(frozen(first:last, 1) - hybrids(1)%coeff)) < TOL, &
                 "the hybrid did not land on the atom it belongs to")
      if (allocated(error)) return
      call check(error, sum(abs(frozen(:first - 1, 1))) + sum(abs(frozen(last + 1:, 1))) < TOL, &
                 "the hybrid left weight outside its own atom's block")
      if (allocated(error)) return
      call check(error, n_occ, 0, "a virtual boundary was counted as occupied")
   end subroutine test_place

   subroutine test_order(error)
      !! Occupied boundaries come first, whatever order they arrived in
      !!
      !! The constraint names its blocks by index range, so an occupied orbital
      !! placed after a virtual one would be held at the level shift -- the bond
      !! pair pushed out of the fragment that is supposed to hold it.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(afo_hybrid_t) :: hybrids(2)
      real(dp), allocatable :: frozen(:, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_occ

      call ethane_molecule(mol, err)
      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      ! Virtual on carbon 1 given first, occupied on carbon 5 given second.
      hybrids(1)%coeff = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      hybrids(2)%coeff = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      call build_group_frozen(mol, [1, 5], [.false., .true.], hybrids, frozen, n_occ, err)
      call check(error,.not. err%has_error(), "placing the hybrids failed")
      if (allocated(error)) return
      call check(error, n_occ, 1, "the occupied boundary was not counted")
      if (allocated(error)) return

      ! Column 1 must be the occupied one -- carbon 5's second function.
      call check(error, abs(frozen(offsets(5) + 2, 1) - 1.0_dp) < TOL, &
                 "the occupied boundary is not in the leading column")
      if (allocated(error)) return
      call check(error, abs(frozen(offsets(1) + 1, 2) - 1.0_dp) < TOL, &
                 "the virtual boundary is not after the occupied one")
   end subroutine test_order

   subroutine test_wrong_basis(error)
      !! A hybrid of the wrong length was built against a different basis
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(afo_hybrid_t) :: hybrids(1)
      real(dp), allocatable :: frozen(:, :)
      integer :: n_occ

      call ethane_molecule(mol, err)
      hybrids(1)%coeff = [0.1_dp, 0.2_dp, 0.3_dp]
      call build_group_frozen(mol, [1], [.true.], hybrids, frozen, n_occ, err)
      call check(error, err%has_error(), &
                 "a hybrid with the wrong number of coefficients was accepted")
   end subroutine test_wrong_basis

   subroutine ethane(sys)

      type(system_geometry_t), intent(out) :: sys
      real(dp) :: xyz(3, 8)

      xyz = reshape([0.000_dp, 0.000_dp, 0.768_dp, &
                     -1.019_dp, 0.000_dp, 1.157_dp, &
                     0.510_dp, 0.883_dp, 1.157_dp, &
                     0.510_dp, -0.883_dp, 1.157_dp, &
                     0.000_dp, 0.000_dp, -0.768_dp, &
                     1.019_dp, 0.000_dp, -1.157_dp, &
                     -0.510_dp, -0.883_dp, -1.157_dp, &
                     -0.510_dp, 0.883_dp, -1.157_dp], [3, 8])

      sys%total_atoms = 8
      sys%n_monomers = 0
      allocate (sys%element_numbers(8))
      sys%element_numbers = [6, 1, 1, 1, 6, 1, 1, 1]
      allocate (sys%coordinates(3, 8))
      sys%coordinates = to_bohr(xyz)
   end subroutine ethane

end module test_mqc_afo_orbital

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_afo_orbital, only: collect_mqc_afo_orbital
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_afo_orbital", collect_mqc_afo_orbital)]
   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
