module test_mqc_czt_guess
   !! Pins the free-atom initial guesses against what a reference energy cannot see.
   !!
   !! The validation manifest already checks that the CPU SCF lands on the right
   !! number. None of those cases would notice most of the ways a guess can be
   !! wrong, because a guess that is merely bad still converges to the same
   !! answer -- slower, and with nothing to say so. What is checked here instead:
   !!
   !!   * every guess reaches the same energy. This is the one invariant a guess
   !!     must not break, and the only one whose failure is a correctness bug
   !!     rather than a performance one;
   !!   * SAD is spherically symmetric and SAC is not. That is the entire
   !!     difference between the two, and without it they would be one guess
   !!     under two names -- superposed atomic coefficients are row-disjoint, so
   !!     the density they imply is exactly the sum of the atomic densities;
   !!   * averaging conserves electrons, and is idempotent;
   !!   * the pseudo-orbital factorisation reproduces its density exactly, which
   !!     is what makes the fitted exchange exact on the guess iteration;
   !!   * a closed-shell free atom is already spherical, so SAC and SAD must
   !!     agree on it -- the two guesses can only differ where a shell is open;
   !!   * the guess densities are block-diagonal over atoms, since nothing in a
   !!     superposition of free atoms couples them;
   !!   * an unknown guess name is refused rather than quietly becoming another.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, &
                                subshell_layout
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, &
                          density_pseudo_orbitals, &
                          SCF_GUESS_CORE, SCF_GUESS_GWH, SCF_GUESS_SAC, SCF_GUESS_SAD
   use mqc_czt_atomic_guess, only: build_atomic_guess, parse_guess_name, &
                                   hund_multiplicity, spherical_average, &
                                   clear_atomic_cache
   implicit none
   private
   public :: collect_mqc_czt_guess_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: E_TOL = 1.0e-10_dp
   real(dp), parameter :: D_TOL = 1.0e-8_dp

contains

   subroutine collect_mqc_czt_guess_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("guess_choice_does_not_change_the_energy", test_same_energy), &
                  new_unittest("sad_is_spherical_and_sac_is_not", test_sad_symmetry), &
                  new_unittest("averaging_conserves_electrons", test_average_trace), &
                  new_unittest("cartesian_d_conserves_electrons", test_cartesian_trace), &
                  new_unittest("averaging_is_idempotent", test_average_idempotent), &
                  new_unittest("pseudo_orbitals_rebuild_their_density", test_pseudo), &
                  new_unittest("closed_shell_atom_gives_one_guess", test_closed_atom), &
                  new_unittest("guess_is_block_diagonal_over_atoms", test_block_diagonal), &
                  new_unittest("hund_multiplicities_are_right", test_hund), &
                  new_unittest("unknown_guess_name_is_refused", test_bad_name) &
                  ]
   end subroutine collect_mqc_czt_guess_tests

   subroutine water(mol, err, basis)
      !! Water at a geometry with no symmetry to flatter a guess
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
   end subroutine water

   subroutine single_atom(z, symbol, basis, mol, err, cartesian)
      !! One free atom, so its guess block is the whole matrix
      integer, intent(in) :: z
      character(len=*), intent(in) :: symbol, basis
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      logical, intent(in), optional :: cartesian
      real(dp) :: c(3, 1)
      character(len=2) :: sym(1)

      sym(1) = symbol
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1])
      if (present(cartesian)) then
         call build_czt_molecule([z], sym, c, basis, mol, err, force_cartesian=cartesian)
      else
         call build_czt_molecule([z], sym, c, basis, mol, err)
      end if
   end subroutine single_atom

   subroutine test_same_energy(error)
      !! Four guesses, one answer
      !!
      !! A guess decides how the SCF gets there, not where. The exception is an
      !! open-shell system that can converge onto a different stationary point,
      !! which is why this uses a closed-shell molecule at its minimum: here
      !! there is one solution, and any disagreement between guesses is a defect
      !! in whichever one moved.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: r(4)
      real(dp), allocatable :: g_a(:, :), g_b(:, :), total(:, :)
      integer :: i
      integer, parameter :: KINDS(4) = [SCF_GUESS_CORE, SCF_GUESS_GWH, &
                                        SCF_GUESS_SAC, SCF_GUESS_SAD]

      call clear_atomic_cache()
      call water(mol, err, "cc-pvdz")
      call check(error,.not. err%has_error(), "water must build: "//err%get_full_trace())
      if (allocated(error)) return

      do i = 1, 4
         if (KINDS(i) == SCF_GUESS_SAC .or. KINDS(i) == SCF_GUESS_SAD) then
            call build_atomic_guess(mol, KINDS(i), g_a, g_b, err)
            call check(error,.not. err%has_error(), "the atomic guess must build: "//err%get_full_trace())
            if (allocated(error)) return
            total = g_a + g_b
            call run_czt_rhf(mol, 10, 100, E_TOL, D_TOL, .false., r(i), err, &
                             in_core=.true., guess=KINDS(i), guess_density=total)
            deallocate (g_a, g_b, total)
         else
            call run_czt_rhf(mol, 10, 100, E_TOL, D_TOL, .false., r(i), err, &
                             in_core=.true., guess=KINDS(i))
         end if
         call check(error,.not. err%has_error(), "the SCF must run: "//err%get_full_trace())
         if (allocated(error)) return
         call check(error, r(i)%converged, "every guess must converge")
         if (allocated(error)) return
      end do

      do i = 2, 4
         call check(error, abs(r(i)%energy - r(1)%energy) < 1.0e-9_dp, &
                    "a guess changed the converged energy")
         if (allocated(error)) return
      end do

      call mol%destroy()
   end subroutine test_same_energy

   subroutine test_sad_symmetry(error)
      !! Carbon's 2p block: round under SAD, not round under SAC
      !!
      !! The load-bearing test of the whole module. A free carbon at Hund
      !! multiplicity puts its two 2p electrons in two of three degenerate
      !! orbitals, so the converged density singles out a direction that the
      !! molecule it is a guess for knows nothing about. SAD averages that away
      !! and SAC keeps it, and if SAD did not, the two guesses would be the same
      !! guess. Checked on the density itself rather than on an iteration count,
      !! because an iteration count would not distinguish "averaged" from
      !! "averaged wrongly".
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: sad_a(:, :), sad_b(:, :), sac_a(:, :), sac_b(:, :)
      integer, allocatable :: ang(:), first(:), ncomp(:)
      integer :: n_sub, k, p, i, j
      real(dp) :: spread_sad, spread_sac, off_sad

      call clear_atomic_cache()
      call single_atom(6, "C ", "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "carbon must build: "//err%get_full_trace())
      if (allocated(error)) return

      call build_atomic_guess(mol, SCF_GUESS_SAD, sad_a, sad_b, err)
      call check(error,.not. err%has_error(), "the SAD guess must build: "//err%get_full_trace())
      if (allocated(error)) return
      call build_atomic_guess(mol, SCF_GUESS_SAC, sac_a, sac_b, err)
      call check(error,.not. err%has_error(), "the SAC guess must build: "//err%get_full_trace())
      if (allocated(error)) return

      ! SAD is spin-restricted by construction, so the two spins must be equal.
      call check(error, maxval(abs(sad_a - sad_b)) < 1.0e-12_dp, &
                 "SAD must give both spins the same density")
      if (allocated(error)) return

      call subshell_layout(mol, ang, first, ncomp, n_sub)
      p = 0
      do k = 1, n_sub
         if (ang(k) == 1) p = k
      end do
      call check(error, p > 0, "carbon in sto-3g must have a p shell")
      if (allocated(error)) return

      ! Within the p shell: how far apart are the three diagonal entries, and how
      ! large are the entries that a spherical density has to vanish on?
      spread_sad = 0.0_dp
      spread_sac = 0.0_dp
      off_sad = 0.0_dp
      do i = 1, ncomp(p)
         do j = 1, ncomp(p)
            if (i == j) then
               spread_sad = max(spread_sad, abs(sad_a(first(p) + i, first(p) + i) &
                                                - sad_a(first(p) + 1, first(p) + 1)))
               spread_sac = max(spread_sac, abs(sac_a(first(p) + i, first(p) + i) &
                                                - sac_a(first(p) + 1, first(p) + 1)))
            else
               off_sad = max(off_sad, abs(sad_a(first(p) + i, first(p) + j)))
            end if
         end do
      end do

      call check(error, spread_sad < 1.0e-12_dp, &
                 "SAD's p diagonal must be equal in all three components")
      if (allocated(error)) return
      call check(error, off_sad < 1.0e-12_dp, &
                 "SAD must have no coupling between p components")
      if (allocated(error)) return
      ! The point of the pair: SAC really does single out a direction. Carbon's
      ! two 2p electrons make this a difference of order one electron, not a
      ! rounding difference, so a loose threshold is honest here.
      call check(error, spread_sac > 0.1_dp, &
                 "SAC on carbon should be spatially polarised; if it is not, the "// &
                 "free atom converged somewhere symmetric and SAD is testing nothing")

      call mol%destroy()
   end subroutine test_sad_symmetry

   subroutine test_average_trace(error)
      !! Averaging moves electrons around, never destroys them
      !!
      !! Tr(D S) is the electron count. A free atom's overlap is itself
      !! spherically symmetric, so projecting D onto its symmetric part cannot
      !! change that contraction -- which makes this an exact identity rather
      !! than an approximation, and a sharp check on the component bookkeeping:
      !! averaging over the wrong number of components would show up here as a
      !! fractional electron.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d_a(:, :), d_b(:, :), s(:, :), total(:, :), avg(:, :)
      real(dp) :: n_before, n_after
      logical :: exact

      call clear_atomic_cache()
      call single_atom(8, "O ", "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "oxygen must build: "//err%get_full_trace())
      if (allocated(error)) return

      call build_atomic_guess(mol, SCF_GUESS_SAC, d_a, d_b, err)
      call check(error,.not. err%has_error(), "the SAC guess must build: "//err%get_full_trace())
      if (allocated(error)) return

      call mol%overlap(s)
      allocate (total(mol%nao, mol%nao), avg(mol%nao, mol%nao))
      total = d_a + d_b
      call spherical_average(mol, total, avg, exact)

      n_before = sum(total*s)
      n_after = sum(avg*s)
      call check(error, abs(n_before - 8.0_dp) < 1.0e-6_dp, &
                 "the free oxygen density must carry eight electrons")
      if (allocated(error)) return
      call check(error, abs(n_after - n_before) < 1.0e-10_dp, &
                 "spherical averaging must conserve the electron count")
      if (allocated(error)) return
      ! cc-pVDZ oxygen has a d shell and is spherical, so the average is exact.
      call check(error, exact, "a spherical basis must average exactly")

      call mol%destroy()
   end subroutine test_average_trace

   subroutine test_cartesian_trace(error)
      !! The same electron count, in the convention where l is not a good label
      !!
      !! `test_average_trace` next door checks this on a SPHERICAL basis, where
      !! it held while this did not. A Cartesian shell of l>=2 is not pure l --
      !! (xx+yy+zz) is s-like -- so a spherically symmetric density genuinely
      !! carries s-d entries, and the averaging used to zero them along with
      !! every other cross-l block. Oxygen came back with Tr(D S) = 8.014 in
      !! cc-pVDZ and 7.933 in cc-pVTZ, which adds f.
      !!
      !! Nothing noticed, because a guess that has lost a hundredth of an
      !! electron still converges to the right energy from a slightly worse
      !! start. It took an outside consumer reading the densities out and
      !! integrating them to see it, which is the argument for checking a guess
      !! against something other than the energy it eventually reaches.
      !!
      !! cc-pVTZ is here as well as cc-pVDZ because the error grew with angular
      !! momentum: a fix that only handled d would pass on one and fail here.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d_a(:, :), d_b(:, :), s(:, :)
      real(dp) :: n_electrons
      integer :: icase
      character(len=8), parameter :: bases(2) = [character(len=8) :: "cc-pvdz", "cc-pvtz"]

      do icase = 1, size(bases)
         call clear_atomic_cache()
         call single_atom(8, "O ", trim(bases(icase)), mol, err, cartesian=.true.)
         call check(error,.not. err%has_error(), "oxygen must build: "//err%get_full_trace())
         if (allocated(error)) return

         call build_atomic_guess(mol, SCF_GUESS_SAD, d_a, d_b, err)
         call check(error,.not. err%has_error(), "the SAD guess must build: "// &
                    err%get_full_trace())
         if (allocated(error)) return

         call mol%overlap(s)
         n_electrons = sum((d_a + d_b)*s)
         call check(error, abs(n_electrons - 8.0_dp) < 1.0e-8_dp, &
                    "a Cartesian "//trim(bases(icase))//" oxygen guess must carry 8 "// &
                    "electrons, not "//real_text(n_electrons))
         if (allocated(error)) return
         deallocate (s, d_a, d_b)
         call mol%destroy()
      end do
   end subroutine test_cartesian_trace

   function real_text(x) result(text)
      !! A number in a message, so a failure says how far off it was
      real(dp), intent(in) :: x
      character(len=:), allocatable :: text
      character(len=32) :: buffer
      write (buffer, "(f16.9)") x
      text = trim(adjustl(buffer))
   end function real_text

   subroutine test_average_idempotent(error)
      !! Averaging an averaged density changes nothing
      !!
      !! A projection applied twice is the projection. This catches the
      !! plausible-looking mistake of averaging over the wrong set -- summing
      !! across shells rather than within one, say -- which conserves electrons
      !! and looks symmetric but is not a projection and so moves on the second
      !! application.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d_a(:, :), d_b(:, :), once(:, :), twice(:, :), total(:, :)
      logical :: exact

      call clear_atomic_cache()
      call single_atom(7, "N ", "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "nitrogen must build: "//err%get_full_trace())
      if (allocated(error)) return

      call build_atomic_guess(mol, SCF_GUESS_SAC, d_a, d_b, err)
      call check(error,.not. err%has_error(), "the SAC guess must build: "//err%get_full_trace())
      if (allocated(error)) return

      allocate (total(mol%nao, mol%nao), once(mol%nao, mol%nao), twice(mol%nao, mol%nao))
      total = d_a + d_b
      call spherical_average(mol, total, once, exact)
      call spherical_average(mol, once, twice, exact)

      call check(error, maxval(abs(twice - once)) < 1.0e-13_dp, &
                 "spherical averaging must be idempotent")

      call mol%destroy()
   end subroutine test_average_idempotent

   subroutine test_pseudo(error)
      !! D = 2 sum_i c_i c_i^T, exactly
      !!
      !! What makes the fitted exchange exact on the guess iteration. The DF Fock
      !! build takes occupied orbitals, a guess density has none, and this
      !! factorisation is the substitute -- so if the identity does not hold, the
      !! density-fitted path starts from an exchange matrix belonging to a
      !! different density and nothing downstream can tell.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d_a(:, :), d_b(:, :), total(:, :), c(:, :), rebuilt(:, :)
      integer :: n_modes, i, j, k

      call clear_atomic_cache()
      call water(mol, err, "cc-pvdz")
      call check(error,.not. err%has_error(), "water must build: "//err%get_full_trace())
      if (allocated(error)) return

      call build_atomic_guess(mol, SCF_GUESS_SAD, d_a, d_b, err)
      call check(error,.not. err%has_error(), "the SAD guess must build: "//err%get_full_trace())
      if (allocated(error)) return

      allocate (total(mol%nao, mol%nao), rebuilt(mol%nao, mol%nao))
      total = d_a + d_b
      call density_pseudo_orbitals(total, c, n_modes, err)
      call check(error,.not. err%has_error(), "the factorisation must succeed: "//err%get_full_trace())
      if (allocated(error)) return

      rebuilt = 0.0_dp
      do j = 1, mol%nao
         do i = 1, mol%nao
            do k = 1, n_modes
               rebuilt(i, j) = rebuilt(i, j) + 2.0_dp*c(i, k)*c(j, k)
            end do
         end do
      end do

      call check(error, maxval(abs(rebuilt - total)) < 1.0e-11_dp, &
                 "the pseudo-orbitals must rebuild their density exactly")

      call mol%destroy()
   end subroutine test_pseudo

   subroutine test_closed_atom(error)
      !! On neon, SAC and SAD are the same guess
      !!
      !! A closed-shell atom's converged density is already spherical and already
      !! spin-symmetric, so averaging has nothing to do and halving the total
      !! returns each spin's own block. The two guesses can only differ where a
      !! shell is open, and this is the case that says so.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: sad_a(:, :), sad_b(:, :), sac_a(:, :), sac_b(:, :)

      call clear_atomic_cache()
      call single_atom(10, "Ne", "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "neon must build: "//err%get_full_trace())
      if (allocated(error)) return

      call build_atomic_guess(mol, SCF_GUESS_SAD, sad_a, sad_b, err)
      call check(error,.not. err%has_error(), "the SAD guess must build: "//err%get_full_trace())
      if (allocated(error)) return
      call build_atomic_guess(mol, SCF_GUESS_SAC, sac_a, sac_b, err)
      call check(error,.not. err%has_error(), "the SAC guess must build: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, maxval(abs(sad_a - sac_a)) < 1.0e-9_dp, &
                 "a closed-shell atom must give the same alpha density either way")
      if (allocated(error)) return
      call check(error, maxval(abs(sad_b - sac_b)) < 1.0e-9_dp, &
                 "a closed-shell atom must give the same beta density either way")

      call mol%destroy()
   end subroutine test_closed_atom

   subroutine test_block_diagonal(error)
      !! Nothing in a superposition of free atoms couples two atoms
      !!
      !! Also the check that the atomic blocks land where they belong: a block
      !! written at the wrong offset would leave a nonzero entry off the diagonal
      !! blocks, and the SCF would still converge from it.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d_a(:, :), d_b(:, :)
      integer :: i, j
      real(dp) :: worst

      call clear_atomic_cache()
      call water(mol, err, "sto-3g")
      call check(error,.not. err%has_error(), "water must build: "//err%get_full_trace())
      if (allocated(error)) return

      call build_atomic_guess(mol, SCF_GUESS_SAD, d_a, d_b, err)
      call check(error,.not. err%has_error(), "the SAD guess must build: "//err%get_full_trace())
      if (allocated(error)) return

      ! sto-3g water: oxygen takes functions 1-5, the two hydrogens 6 and 7.
      worst = 0.0_dp
      do j = 1, mol%nao
         do i = 1, mol%nao
            if ((i <= 5) .neqv. (j <= 5)) worst = max(worst, abs(d_a(i, j)))
            if (i == 6 .and. j == 7) worst = max(worst, abs(d_a(i, j)))
            if (i == 7 .and. j == 6) worst = max(worst, abs(d_a(i, j)))
         end do
      end do

      call check(error, worst < 1.0e-14_dp, &
                 "the guess density must not couple two atoms")

      call mol%destroy()
   end subroutine test_block_diagonal

   subroutine test_hund(error)
      !! Ground-state multiplicities across the range the CPU path is used on
      type(error_type), allocatable, intent(out) :: error

      call check(error, hund_multiplicity(1) == 2, "H is a doublet")
      if (allocated(error)) return
      call check(error, hund_multiplicity(2) == 1, "He is a singlet")
      if (allocated(error)) return
      call check(error, hund_multiplicity(6) == 3, "C is a triplet")
      if (allocated(error)) return
      call check(error, hund_multiplicity(7) == 4, "N is a quartet")
      if (allocated(error)) return
      call check(error, hund_multiplicity(8) == 3, "O is a triplet")
      if (allocated(error)) return
      call check(error, hund_multiplicity(9) == 2, "F is a doublet")
      if (allocated(error)) return
      call check(error, hund_multiplicity(10) == 1, "Ne is a singlet")
      if (allocated(error)) return
      call check(error, hund_multiplicity(16) == 3, "S is a triplet")
      if (allocated(error)) return
      call check(error, hund_multiplicity(18) == 1, "Ar is a singlet")
      if (allocated(error)) return
      ! Every neutral atom's unpaired count has the parity of its electron count,
      ! which is what lets the atomic UHF pair them at all.
      call check(error, mod(hund_multiplicity(17) - 1, 2) == mod(17, 2), &
                 "Cl's unpaired count must have the parity of its electron count")
   end subroutine test_hund

   subroutine test_bad_name(error)
      !! An unrecognised guess is refused, not silently substituted
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer :: kind

      call parse_guess_name("sap", kind, err)
      call check(error, err%has_error(), "an unknown guess name must be refused")
      if (allocated(error)) return

      call err%clear()
      call parse_guess_name("auto", kind, err)
      call check(error,.not. err%has_error(), "'auto' must be accepted: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, kind == SCF_GUESS_SAD, "'auto' must resolve to SAD on the CPU path")
      if (allocated(error)) return

      call err%clear()
      call parse_guess_name("gwh", kind, err)
      call check(error, kind == SCF_GUESS_GWH, "'gwh' must resolve to GWH")
   end subroutine test_bad_name

end module test_mqc_czt_guess

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_guess, only: collect_mqc_czt_guess_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_guess", collect_mqc_czt_guess_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
