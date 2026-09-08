!! That a fragment potential with f functions survives a round trip
module test_mqc_czt_efp_f_functions
   !! **The f block of the format had never been written to a file and read back.**
   !!
   !! `test_mqc_czt_efp_read` round trips a 6-31G* potential, whose highest shell
   !! is d. An f basis is a wider file -- ten coefficients a shell instead of six,
   !! more points in the projection basis, more values a record -- and the format
   !! is fixed-column with continuation markers, so "it works for d" says nothing
   !! about whether the last value of an f record survives. That is what this
   !! pins, block by block against the in-memory object.
   !!
   !! **What it deliberately does not pin, and cannot.** Whether
   !! `F_FROM_LIBCINT`, `F_CLASS` and `F_NORMALIZATION` describe *GAMESS's*
   !! ordering is not testable from inside this program. `to_gamess_ao_order` and
   !! `from_gamess_ao_order` read the same three constants, so any error in them
   !! cancels exactly on a round trip -- including in `C^T S C = I`, since the
   !! coefficients that reach the overlap are the ones that went in. Measured:
   !! swapping two entries of `F_FROM_LIBCINT` leaves every case here passing.
   !! Only handing the file to GAMESS distinguishes a right map from a
   !! self-consistent wrong one, which is `tools/efp_validation/dimer_energy.py`.
   !!
   !! What the two orbital cases *do* pin is that the two routines remain inverses
   !! of each other and that the recovered orbitals are orthonormal against our
   !! own integrals -- so an edit to one of them alone, which is how a matched
   !! pair actually breaks, fails here rather than in a validation run.
   !!
   !! **The basis.** 6-311++G(3df,3pd), the recommended EFMO basis with hydrogen's
   !! polarisation set left as the Basis Set Exchange ships it: three d shells and
   !! an f on oxygen, three p and a d on each hydrogen, 83 Cartesian functions.
   !! The whole file runs in about seven seconds at four threads, so there was no
   !! reason to settle for a smaller set -- and the diffuse functions are what
   !! make the CTVEC block large enough to expose the format's *relative*
   !! precision, which is the one tolerance here that had to be measured rather
   !! than assumed.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_efp_potential, only: efp_potential_t, make_efp_potential, &
                                    write_efp_potential, from_gamess_ao_order, &
                                    to_gamess_ao_order
   use mqc_czt_efp_read, only: efp_fragment_t, read_efp_potential
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_efp_f_function_tests

   !! GAMESS's Bohr, matching the emitter.
   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp

   !! The file carries ten decimals on a parameter, so this is what a round trip
   !! of one can preserve.
   real(dp), parameter :: FORMAT_TOL = 1.0e-9_dp

   !! Orbital coefficients are written `ES15.8E2`, eight decimals of mantissa, so
   !! they come back an order of magnitude coarser than the parameters do.
   real(dp), parameter :: ORBITAL_TOL = 1.0e-8_dp

   !! `C^T S C - I` over coefficients carrying `ORBITAL_TOL`, summed over 83 basis
   !! functions. Loose enough for that accumulation and far tighter than any
   !! plausible mapping error, which moves an element by a factor of sqrt(5) or
   !! sqrt(15) rather than by 1e-7.
   real(dp), parameter :: ORTHONORMAL_TOL = 1.0e-6_dp

   character(len=*), parameter :: BASIS = "6-311++g(3df,3pd)"

contains

   subroutine collect_efp_f_function_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("efp_f_basis_has_an_f_shell", test_has_f_shell), &
                  new_unittest("efp_f_round_trip_blocks", test_blocks), &
                  new_unittest("efp_f_orbitals_are_orthonormal", test_orthonormal), &
                  new_unittest("efp_f_permutation_inverts", test_permutation) &
                  ]
   end subroutine collect_efp_f_function_tests

   subroutine geometry(z, symbols, c)
      !! One water, in a frame with no zero coordinate
      !!
      !! **Not the planar frame the other EFP tests use.** An f function with an
      !! odd power of y is exactly zero on every atom of a molecule lying in the
      !! xz plane, so a coefficient sent to the wrong slot could be zero on both
      !! sides and the comparison would pass. Tilting the molecule out of the
      !! plane is what makes every one of the ten f slots carry a distinct number
      !! -- the same reason `F_NORMALIZATION` was only solvable in a tilted frame.
      integer, intent(out) :: z(3)
      character(len=2), intent(out) :: symbols(3)
      real(dp), intent(out) :: c(3, 3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.1731_dp*ANG, 0.4327_dp*ANG, 0.8385_dp*ANG, &
                   0.9268_dp*ANG, 0.1102_dp*ANG, -0.2400_dp*ANG], [3, 3])
   end subroutine geometry

   subroutine water(pot, err)
      !! One water potential in a basis with an f shell
      type(efp_potential_t), intent(out) :: pot
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)

      call geometry(z, symbols, c)
      call make_efp_potential(z, symbols, c, BASIS, "WATERF", pot, err)
   end subroutine water

   subroutine molecule(mol, err)
      !! The same molecule the potential was built from, Cartesian
      !!
      !! `force_cartesian` because the potential is written in GAMESS's Cartesian
      !! ordering and the map is a Cartesian one; a spherical build would have a
      !! different `nao` and the comparison would not even be shaped right.
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)

      call geometry(z, symbols, c)
      call build_czt_molecule(z, symbols, c, BASIS, mol, err, force_cartesian=.true.)
   end subroutine molecule

   subroutine round_trip(pot, frag, path, err)
      !! Build, write, read back
      type(efp_potential_t), intent(out) :: pot
      type(efp_fragment_t), intent(out) :: frag
      character(len=*), intent(in) :: path
      type(error_t), intent(inout) :: err

      call water(pot, err)
      if (err%has_error()) return
      call write_efp_potential(pot, path, err)
      if (err%has_error()) return
      call read_efp_potential(path, frag, err)
   end subroutine round_trip

   subroutine test_has_f_shell(error)
      !! The guard that keeps the rest of this file from being vacuous
      !!
      !! If the basis on the search path ever loses its f shell -- renamed,
      !! replaced, resolved to something else by `find_basis_file` -- every other
      !! case here would still pass while testing the d map a second time.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: frag
      type(error_t) :: err
      character(len=*), parameter :: path = "test_efp_f_guard.efp"

      call molecule(mol, err)
      call check(error,.not. err%has_error(), &
                 "building the molecule failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, mol%cartesian, "the molecule did not build Cartesian")
      if (allocated(error)) return
      call check(error, mol%nao == 83, "expected 83 Cartesian functions")
      if (allocated(error)) return
      call mol%destroy()

      ! Read back from the file rather than asked of the molecule, because the
      ! projection basis is what a reader of this potential sees: an f shell that
      ! reached the file is an f shell the map had to place.
      call round_trip(pot, frag, path, err)
      call check(error,.not. err%has_error(), &
                 "the round trip failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, count(frag%shell_l == 3) > 0, &
                 "this basis has no f shell, so the f map is untested")
      if (allocated(error)) return
      call check(error, count(frag%shell_l == 2) > 0, &
                 "this basis has no d shell")

      call frag%destroy()
      call pot%destroy()
      call delete(path)
   end subroutine test_has_f_shell

   subroutine test_blocks(error)
      !! Every parameter block, in against out
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: frag
      type(error_t) :: err
      character(len=*), parameter :: path = "test_efp_f_blocks.efp"
      integer :: i, k, f, a, b, c, e, slot
      real(dp) :: scale

      call round_trip(pot, frag, path, err)
      call check(error,.not. err%has_error(), &
                 "the round trip failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, frag%n_points == pot%n_points, "point count changed")
      if (allocated(error)) return
      call check(error, frag%n_atoms == pot%n_atoms, "atom count changed")
      if (allocated(error)) return
      call check(error, frag%multiplicity == pot%multiplicity, "multiplicity changed")
      if (allocated(error)) return

      ! --- multipoles and the geometry they sit on ------------------------
      do i = 1, frag%n_points
         call worst(error, frag%points(:, i), pot%points(:, i), FORMAT_TOL, &
                    "an expansion point moved")
         if (allocated(error)) return
         call check(error, frag%mass(i), pot%mass(i), thr=FORMAT_TOL, &
                    message="a mass changed")
         if (allocated(error)) return
         call check(error, frag%charge(i), pot%charge(i), thr=FORMAT_TOL, &
                    message="a nuclear charge changed")
         if (allocated(error)) return
         call check(error, frag%q_elec(i), pot%q_elec(i), thr=FORMAT_TOL, &
                    message="an electronic monopole changed")
         if (allocated(error)) return
         call check(error, frag%q_nuc(i), pot%q_nuc(i), thr=FORMAT_TOL, &
                    message="a nuclear monopole changed")
         if (allocated(error)) return
         call worst(error, frag%dipole(:, i), pot%dipole(:, i), FORMAT_TOL, &
                    "a dipole component changed")
         if (allocated(error)) return
         call worst(error, frag%quadrupole(:, i), pot%quadrupole(:, i), FORMAT_TOL, &
                    "a quadrupole component changed")
         if (allocated(error)) return
         call worst(error, frag%octopole(:, i), pot%octopole(:, i), FORMAT_TOL, &
                    "an octupole component changed")
         if (allocated(error)) return
      end do

      ! --- screening ------------------------------------------------------
      call check(error, frag%has_screen .and. frag%has_screen2, &
                 "a screening block was not read")
      if (allocated(error)) return
      call worst(error, frag%screen, pot%screen, FORMAT_TOL, &
                 "a Gaussian screening exponent changed")
      if (allocated(error)) return
      call worst(error, frag%screen2, pot%screen2, FORMAT_TOL, &
                 "an exponential screening exponent changed")
      if (allocated(error)) return

      ! --- polarizabilities, static and dynamic ---------------------------
      call check(error, frag%has_static_pol .and. frag%has_dynamic, &
                 "a polarizability block was not read")
      if (allocated(error)) return
      call check(error, frag%n_pol == pot%n_lmo, "the static point count changed")
      if (allocated(error)) return
      call check(error, frag%n_lmo == pot%n_lmo, "the dynamic point count changed")
      if (allocated(error)) return
      call check(error, frag%n_freq == size(pot%frequencies), "the frequency count changed")
      if (allocated(error)) return
      call worst(error, frag%frequencies, pot%frequencies, 1.0e-6_dp, &
                 "a frequency changed")
      if (allocated(error)) return
      do k = 1, pot%n_lmo
         call worst(error, frag%centroids(:, k), pot%centroids(:, k), FORMAT_TOL, &
                    "a centroid moved")
         if (allocated(error)) return
         call worst(error, reshape(frag%static_pol(:, :, k), [9]), &
                    reshape(pot%static_pol(:, :, k), [9]), FORMAT_TOL, &
                    "a static polarizability changed")
         if (allocated(error)) return
         do f = 1, frag%n_freq
            call worst(error, reshape(frag%dyn_pol(:, :, k, f), [9]), &
                       reshape(pot%dynamic_pol(:, :, k, f), [9]), FORMAT_TOL, &
                       "a dynamic polarizability changed")
            if (allocated(error)) return
         end do
      end do

      ! --- the two higher dispersion blocks, flat in the fragment ---------
      !
      ! Slot for slot against the writer's own indexing, which is the only place
      ! the flat layout is defined: the dipole-quadrupole block runs the first
      ! quadrupole index fastest, the quadrupole-quadrupole one the last of four.
      call check(error, frag%has_dipquad .and. frag%has_quadquad, &
                 "a higher dispersion block was not read")
      if (allocated(error)) return
      do f = 1, frag%n_freq
         do k = 1, pot%n_lmo
            do a = 1, 3
               do b = 1, 3
                  do c = 1, 3
                     slot = (a - 1)*9 + (c - 1)*3 + b
                     call check(error, frag%dipquad(slot, k, f), &
                                pot%dipquad(a, b, c, k, f), thr=FORMAT_TOL, &
                                message="a dipole-quadrupole component changed")
                     if (allocated(error)) return
                     do e = 1, 3
                        slot = ((a - 1)*9 + (b - 1)*3 + (c - 1))*3 + e
                        call check(error, frag%quadquad(slot, k, f), &
                                   pot%quadquad(a, b, c, e, k, f), thr=FORMAT_TOL, &
                                   message="a quadrupole-quadrupole component changed")
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! --- exchange repulsion and charge transfer -------------------------
      call check(error, frag%has_fock, "the LMO Fock matrix was not read")
      if (allocated(error)) return
      call check(error, frag%n_lmo_proj == pot%n_lmo, "the LMO count changed")
      if (allocated(error)) return
      call check(error, frag%nao_proj == pot%nao, "the basis function count changed")
      if (allocated(error)) return
      do k = 1, pot%n_lmo
         call worst(error, frag%fock_lmo(:, k), pot%fock_lmo(:, k), FORMAT_TOL, &
                    "an LMO Fock element changed")
         if (allocated(error)) return
      end do

      call check(error, frag%has_lmo, "the projection wavefunction was not read")
      if (allocated(error)) return
      do k = 1, pot%n_lmo
         call worst(error, frag%lmo_gamess(:, k), pot%orbitals(:, k), ORBITAL_TOL, &
                    "a localized orbital coefficient changed")
         if (allocated(error)) return
      end do

      call check(error, frag%has_ctvec, "CTVEC was not read")
      if (allocated(error)) return
      call check(error, frag%has_ctfok, "CTFOK was not read")
      if (allocated(error)) return
      call check(error, frag%n_mo_ct == pot%nao, "the CTVEC orbital count changed")
      if (allocated(error)) return
      call check(error, frag%n_occ_ct == pot%n_occ, "the CTVEC occupied count changed")
      if (allocated(error)) return
      ! `ORBITAL_TOL` is scaled by the largest coefficient in the block, because
      ! `ES15.8E2` carries *nine significant figures* rather than eight decimals:
      ! an occupied coefficient of order one comes back to 1e-9, while a diffuse
      ! virtual of order a hundred -- which is what `++` on both atoms produces,
      ! and this block is every virtual -- can only come back to 1e-7. An absolute
      ! 1e-8 here fails on the format rather than on the mapping.
      scale = max(1.0_dp, maxval(abs(pot%canonical)))
      do k = 1, frag%n_mo_ct
         call worst(error, frag%ctvec_gamess(:, k), pot%canonical(:, k), &
                    ORBITAL_TOL*scale, "a CTVEC coefficient changed")
         if (allocated(error)) return
      end do
      call worst(error, frag%eps_occ, pot%eps_occ(1:frag%n_occ_ct), FORMAT_TOL, &
                 "an orbital energy changed")
      if (allocated(error)) return

      call frag%destroy()
      call pot%destroy()
      call delete(path)
   end subroutine test_blocks

   subroutine test_orthonormal(error)
      !! `C^T S C = I` for both orbital sets, after the GAMESS order is undone
      !!
      !! Note what this is and is not. The overlap comes from the integral code,
      !! which knows nothing about GAMESS -- but the coefficients reaching it went
      !! out through `to_gamess_ao_order` and came back through
      !! `from_gamess_ao_order`, so a shared error in the constants cancels and
      !! this passes. It fails when the two routines stop agreeing with each
      !! other, and when the file has mangled a coefficient rather than a slot.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: frag
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      character(len=*), parameter :: path = "test_efp_f_orthonormal.efp"
      real(dp), allocatable :: s(:, :), lmo(:, :), ct(:, :)

      call round_trip(pot, frag, path, err)
      call check(error,.not. err%has_error(), &
                 "the round trip failed: "//err%get_full_trace())
      if (allocated(error)) return
      call molecule(mol, err)
      call check(error,.not. err%has_error(), &
                 "building the molecule failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, mol%nao == frag%nao_proj, &
                 "the molecule and the potential disagree on the basis size")
      if (allocated(error)) return

      call mol%overlap(s)

      call from_gamess_ao_order(mol, frag%lmo_gamess, lmo, err)
      call check(error,.not. err%has_error(), &
                 "undoing the LMO ordering failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check_orthonormal(error, s, lmo, "the localized orbitals")
      if (allocated(error)) return

      call from_gamess_ao_order(mol, frag%ctvec_gamess, ct, err)
      call check(error,.not. err%has_error(), &
                 "undoing the CTVEC ordering failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check_orthonormal(error, s, ct, "the charge-transfer orbitals")
      if (allocated(error)) return

      call mol%destroy()
      call frag%destroy()
      call pot%destroy()
      call delete(path)
   end subroutine test_orthonormal

   subroutine check_orthonormal(error, s, c, what)
      !! The largest element of `C^T S C - I`
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: s(:, :), c(:, :)
      character(len=*), intent(in) :: what

      real(dp), allocatable :: metric(:, :)
      integer :: i

      metric = matmul(transpose(c), matmul(s, c))
      do i = 1, size(metric, 1)
         metric(i, i) = metric(i, i) - 1.0_dp
      end do
      call check(error, maxval(abs(metric)) < ORTHONORMAL_TOL, &
                 what//" are not orthonormal against our own overlap, so the "// &
                 "GAMESS ordering was undone wrongly")
   end subroutine check_orthonormal

   subroutine test_permutation(error)
      !! `to_gamess_ao_order` of the recovered coefficients is what was written
      !!
      !! The two routines are used at opposite ends of the format -- one on the
      !! way out of the SCF, one on the way back in from a file -- so nothing else
      !! requires them to be inverses. They are here.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: frag
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      character(len=*), parameter :: path = "test_efp_f_permutation.efp"
      real(dp), allocatable :: lmo(:, :), back(:, :)

      call round_trip(pot, frag, path, err)
      call check(error,.not. err%has_error(), &
                 "the round trip failed: "//err%get_full_trace())
      if (allocated(error)) return
      call molecule(mol, err)
      call check(error,.not. err%has_error(), &
                 "building the molecule failed: "//err%get_full_trace())
      if (allocated(error)) return

      call from_gamess_ao_order(mol, frag%lmo_gamess, lmo, err)
      call check(error,.not. err%has_error(), &
                 "undoing the LMO ordering failed: "//err%get_full_trace())
      if (allocated(error)) return
      call to_gamess_ao_order(mol, lmo, back, err)
      call check(error,.not. err%has_error(), &
                 "redoing the LMO ordering failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(back - frag%lmo_gamess)) < &
                 ORBITAL_TOL*max(1.0_dp, maxval(abs(frag%lmo_gamess))), &
                 "the two ordering maps are not inverses of each other")
      if (allocated(error)) return

      ! And the recovered coefficients are genuinely a different arrangement:
      ! were the f block silently skipped, `lmo` would equal `lmo_gamess` and the
      ! inverse test above would pass without the map having done anything.
      call check(error, maxval(abs(lmo - frag%lmo_gamess)) > 1.0e-3_dp, &
                 "the ordering map left the coefficients untouched")
      if (allocated(error)) return

      call mol%destroy()
      call frag%destroy()
      call pot%destroy()
      call delete(path)
   end subroutine test_permutation

   subroutine worst(error, got, want, tol, message)
      !! The largest disagreement in a vector, as one check
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: got(:), want(:)
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: message

      call check(error, maxval(abs(got - want)) < tol, message)
   end subroutine worst

   subroutine delete(path)
      character(len=*), intent(in) :: path

      integer :: unit, stat

      open (newunit=unit, file=path, status="old", action="read", iostat=stat)
      if (stat == 0) close (unit, status="delete")
   end subroutine delete

end module test_mqc_czt_efp_f_functions

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_efp_f_functions, only: collect_efp_f_function_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_efp_f_functions", collect_efp_f_function_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
