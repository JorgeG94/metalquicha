!! Unit tests for the frozen orbital taken off a model system
module test_mqc_afo_orbital
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: system_geometry_t, to_bohr, to_angstrom
   use mqc_bond_perception, only: find_severed_bonds, severed_bond_t
   use mqc_czt_afo, only: afo_model_t, afo_options_t, build_afo_model, &
                          bond_hybrid, BOND_ORBITAL_REACH, &
                          afo_hybrid_t, build_group_frozen, afo_lmo_set_t, &
                          bond_lmo_set, build_bonded_model, orient_cut
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, atom_ao_blocks
   use mqc_czt_rhf, only: run_czt_rhf, rhf_result_t
   use mqc_czt_localize, only: er_localize
   use mqc_czt_xc, only: xc_context_t, xc_available
   use mqc_czt_fragment_solver, only: fragment_xc_context
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_scf_types, only: scf_numerics_t
   use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use mqc_fock_projector, only: fock_projector_t, build_frozen_basis
   implicit none

   !! The model SCF converges to 1e-10 and the localization to its own sweep
   !! tolerance, so an orbital coefficient is good to somewhere around here.
   real(dp), parameter :: TOL = 1.0e-6_dp

   !! STO-3G on carbon is 1s, 2s, then three 2p. The tests naming these indices
   !! name the basis too, so the layout is a statement about that basis rather
   !! than an assumption about all of them.
   integer, parameter :: P_FIRST = 3, P_LAST = 5

   real(dp), parameter :: NUDGE = 0.03_dp
      !! Bohr; see `ethane_lmo_set`

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
                  new_unittest("a_frozen_hybrid_comes_out_of_the_scf_empty", test_end_to_end), &
                  new_unittest("model_system_under_pbe_is_solved_at_pbe", test_model_kohn_sham), &
                  new_unittest("model_system_without_a_method_is_hartree_fock", test_model_hartree_fock) &
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

      ! Within 1e-3 Bohr rather than exactly: an Edmiston-Ruedenberg core
      ! mixes in a little of its atom's valence and sits a few 1e-4 Bohr off
      ! the nucleus -- 1.6e-4 for ethane's in STO-3G -- where a Boys core is
      ! well inside 1e-4. The reach is a third of a bond away either way.
      n_at_half = count(abs(distance - half) < 1.0e-3_dp)
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

   subroutine ethane_lmo_set(opts, model, set, error)
      !! Ethane cut at the C-C bond, through to its frozen-orbital set
      !!
      !! Every atom is nudged off the D3d geometry by a fixed few hundredths of
      !! a Bohr. Ethane's six C-H bonds are equivalent, so the localization of
      !! the symmetric molecule is decided by rounding, and two runs of it
      !! need not pick the same set of orbitals; the nudge makes it decided by
      !! the geometry.
      type(afo_options_t), intent(in) :: opts
      type(afo_model_t), intent(out) :: model
      type(afo_lmo_set_t), intent(out) :: set
      type(error_t), intent(inout) :: error

      type(system_geometry_t) :: sys
      type(severed_bond_t), allocatable :: cuts(:)
      integer :: n_cuts, n_on_bond, k

      call ethane(sys)
      do k = 1, sys%total_atoms
         sys%coordinates(:, k) = sys%coordinates(:, k) + NUDGE*[sin(1.0_dp*k), &
                                                                cos(2.0_dp*k), sin(3.0_dp*k)]
      end do
      call find_severed_bonds(sys, [1, 1, 1, 1, 2, 2, 2, 2], cuts, n_cuts)
      call orient_cut(sys%element_numbers, sys%coordinates, cuts(1), error)
      call build_bonded_model(sys%element_numbers, sys%coordinates, cuts(1), model, error)
      call bond_lmo_set(model, opts, set, n_on_bond, error)
   end subroutine ethane_lmo_set

   subroutine independent_localized(model, functional, localized, error, mol)
      !! The model's occupied orbitals, ER-localized, from a direct `run_czt_rhf`
      !! at the model's own convergence, with `functional` or Hartree-Fock
      type(afo_model_t), intent(in) :: model
      character(len=*), intent(in) :: functional
      real(dp), allocatable, intent(out) :: localized(:, :)
      type(error_t), intent(inout) :: error
      type(czt_molecule_t), intent(out) :: mol

      type(afo_options_t) :: defaults
      type(cuest_scf_settings_t) :: method
      type(scf_numerics_t) :: numerics
      type(xc_context_t) :: xc
      type(rhf_result_t) :: scf
      real(dp), allocatable :: centroids(:, :)

      call build_czt_molecule(model%z, model%sym, model%xyz, "sto-3g", mol, error)
      if (error%has_error()) return
      numerics = defaults%scf
      numerics%incremental_fock = .false.
      if (len_trim(functional) > 0) then
         method%functional = functional
         call fragment_xc_context(method, mol, xc, error)
         if (error%has_error()) return
         call run_czt_rhf(mol, model%nelec, defaults%scf_max_iter, defaults%scf_energy_tol, &
                          defaults%scf_density_tol, .false., scf, error, scf=numerics, &
                          xc=xc, grad_tol=defaults%scf_grad_tol)
         call xc%destroy()
      else
         call run_czt_rhf(mol, model%nelec, defaults%scf_max_iter, defaults%scf_energy_tol, &
                          defaults%scf_density_tol, .false., scf, error, scf=numerics, &
                          grad_tol=defaults%scf_grad_tol)
      end if
      if (error%has_error()) return
      call er_localize(mol, scf%orbitals, scf%n_occupied, localized, centroids, error)
   end subroutine independent_localized

   function set_rows(model, mol, set) result(rows)
      !! The rows of the model's AO axis a set's coefficients are kept on
      type(afo_model_t), intent(in) :: model
      type(czt_molecule_t), intent(in) :: mol
      type(afo_lmo_set_t), intent(in) :: set
      integer, allocatable :: rows(:)

      integer, allocatable :: offsets(:), counts(:)
      integer :: k, l, n

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (rows(0))
      do k = 1, set%n_at
         do l = 1, model%n_atoms
            if (model%from_system(l) /= set%atoms(k)) cycle
            rows = [rows, [(offsets(l) + n, n=1, counts(l))]]
         end do
      end do
   end function set_rows

   pure function column_defect(a, pool) result(worst)
      !! How far the worst column of `a` is from its nearest column of `pool`,
      !! up to sign -- an orbital's sign is not a property of the orbital
      real(dp), intent(in) :: a(:, :), pool(:, :)
      real(dp) :: worst

      real(dp) :: nearest
      integer :: i, j

      worst = 0.0_dp
      do i = 1, size(a, 2)
         nearest = huge(1.0_dp)
         do j = 1, size(pool, 2)
            nearest = min(nearest, maxval(abs(a(:, i) - pool(:, j))), &
                          maxval(abs(a(:, i) + pool(:, j))))
         end do
         worst = max(worst, nearest)
      end do
   end function column_defect

   subroutine test_model_kohn_sham(error)
      !! With `opts%method` at PBE the model is the Kohn-Sham solution
      !!
      !! The reference is built here from `run_czt_rhf` and an xc context and
      !! localized with `er_localize`, so it shares no code with the model's own
      !! path beyond those two. Every column of the returned set must be one of
      !! its localized orbitals cut down to the kept atoms. The Hartree-Fock set
      !! must not be: were `opts%method` ignored the two sets would be the same
      !! and the second check is the one that would fail.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(afo_options_t) :: opts, opts_hf
      type(afo_model_t) :: model
      type(afo_lmo_set_t) :: set_pbe, set_hf
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: localized(:, :)
      integer, allocatable :: rows(:)
      real(dp) :: matched, moved

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      opts%basis = "sto-3g"
      opts_hf = opts
      allocate (opts%method)
      opts%method%functional = "pbe"

      call ethane_lmo_set(opts, model, set_pbe, err)
      call check(error,.not. err%has_error(), "the PBE model's orbital set failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if
      call ethane_lmo_set(opts_hf, model, set_hf, err)
      call check(error,.not. err%has_error(), "the Hartree-Fock model's orbital set failed")
      if (allocated(error)) return

      call independent_localized(model, "pbe", localized, err, mol)
      call check(error,.not. err%has_error(), "the reference Kohn-Sham model failed")
      if (allocated(error)) return

      rows = set_rows(model, mol, set_pbe)
      call check(error, size(rows), size(set_pbe%coeff, 1), &
                 "the kept rows do not add up to the set's coefficients")
      if (allocated(error)) return
      matched = column_defect(set_pbe%coeff, localized(rows, :))
      call check(error, matched < 1.0e-8_dp, &
                 "the PBE model's orbitals are not the Kohn-Sham ones")
      if (allocated(error)) then
         write (*, *) "   worst column defect against the reference =", matched
         return
      end if

      moved = column_defect(set_pbe%coeff, set_hf%coeff)
      write (*, *) "   PBE set against Hartree-Fock set, worst column defect =", moved
      call check(error, moved > 1.0e-4_dp, &
                 "the PBE model's orbitals are the Hartree-Fock ones: the method was ignored")
   end subroutine test_model_kohn_sham

   subroutine test_model_hartree_fock(error)
      !! An unallocated method is Hartree-Fock, bit for bit what a settings
      !! object with no functional gives, and the orbitals a direct
      !! Hartree-Fock `run_czt_rhf` localizes to
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(afo_options_t) :: opts, opts_named
      type(afo_model_t) :: model
      type(afo_lmo_set_t) :: set_bare, set_named
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: localized(:, :)
      integer, allocatable :: rows(:)
      real(dp) :: matched
      integer :: threads

      opts%basis = "sto-3g"
      opts_named = opts
      allocate (opts_named%method)   ! no functional, no correlation: Hartree-Fock

      ! One thread: a threaded Fock build sums in an order that varies from run
      ! to run, so two identical calls agree to rounding and not to the bit.
      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call ethane_lmo_set(opts, model, set_bare, err)
      call ethane_lmo_set(opts_named, model, set_named, err)
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), "the model's orbital set failed")
      if (allocated(error)) return

      call check(error, all(set_bare%coeff == set_named%coeff), &
                 "an unallocated method is not bit-identical to Hartree-Fock by name")
      if (allocated(error)) return

      call independent_localized(model, "", localized, err, mol)
      call check(error,.not. err%has_error(), "the reference Hartree-Fock model failed")
      if (allocated(error)) return
      rows = set_rows(model, mol, set_bare)
      matched = column_defect(set_bare%coeff, localized(rows, :))
      call check(error, matched < 1.0e-8_dp, &
                 "the Hartree-Fock model's orbitals are not the direct Hartree-Fock ones")
      if (allocated(error)) write (*, *) "   worst column defect =", matched
   end subroutine test_model_hartree_fock

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
