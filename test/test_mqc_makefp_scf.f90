!! That the deck's SCF settings reach the one SCF a MakeFP run performs
module test_mqc_makefp_scf
   !! **These exist because six settings were dropped and nothing said so.**
   !!
   !! `make_efp_potential` used to call its SCF with seven arguments and the
   !! iteration cap written into the call, so `incremental_fock`,
   !! `level_shift`, `linear_dependence`, `accelerator`, `diis_size` and
   !! `maxiter` were read from the deck, validated by the schema, carried
   !! through the config, and then never reached the SCF. A deck setting one of
   !! them ran exactly as though it had not, and printed a fragment potential
   !! that looked fine.
   !!
   !! A deck-level test cannot see that. It can only check the numbers a run
   !! produced, and the whole failure mode is that the numbers are the ones a
   !! *default* run produces. What is needed is a setting whose arrival changes
   !! something observable, which is what these do: each case picks a value
   !! extreme enough that the run cannot complete normally if it arrived, and
   !! asserts the run did not complete normally.
   !!
   !! The last case is the other half, and the more important one for physics:
   !! the settings that *should* be invisible must actually be invisible. An
   !! accelerator, a DIIS subspace and a level shift change how a density is
   !! reached and not which density it is, so the multipoles fitted to it have
   !! to come out the same. That catches a setting arriving and being wired to
   !! the wrong thing, which is the failure the arrival cases cannot see.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_scf_types, only: scf_numerics_t
   use mqc_czt_efp_potential, only: efp_potential_t, make_efp_potential
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_makefp_scf

   !! Angstrom to Bohr, so the geometry below can be written the way it is
   !! measured. The same number `mqc_physical_constants` carries; spelled here
   !! only to keep this file free of a use-association for one literal.
   real(dp), parameter :: ANG = 1.8897261254578281_dp

   !! A dipole this far apart is a different density, not a different route to
   !! the same one. Chosen well above the OpenMP reduction scatter, which on
   !! this size of problem sits around 1e-12.
   real(dp), parameter :: DIPOLE_TOL = 1.0e-8_dp

contains

   subroutine collect_makefp_scf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("makefp_baseline_converges", test_baseline), &
                  new_unittest("makefp_honours_maxiter", test_maxiter), &
                  new_unittest("makefp_honours_diis_off", test_diis_off), &
                  new_unittest("makefp_refuses_a_bad_accelerator", test_bad_accelerator), &
                  new_unittest("makefp_honours_level_shift", test_level_shift), &
                  new_unittest("makefp_honours_diis_size", test_diis_size), &
                  new_unittest("makefp_honours_linear_dependence", test_lindep), &
                  new_unittest("makefp_honours_gradient_tolerance", test_grad_tol), &
                  new_unittest("makefp_route_does_not_move_the_multipoles", test_invariance) &
                  ]
   end subroutine collect_makefp_scf

   subroutine water(pot, err, scf, max_iter, grad_tol, energy_tol)
      !! One water potential in a small basis, optionally given SCF settings
      !!
      !! STO-3G rather than the 6-31G* the other EFP tests use: nothing here
      !! looks at the quality of the potential, only at whether a setting
      !! arrived, so the basis only has to be big enough to have an SCF worth
      !! converging.
      type(efp_potential_t), intent(out) :: pot
      type(error_t), intent(inout) :: err
      type(scf_numerics_t), intent(in), optional :: scf
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: grad_tol, energy_tol

      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      ! Both optionals pass straight through. An absent actual argument arrives
      ! absent at the optional dummy, so one call serves every case here.
      call make_efp_potential(z, symbols, c, "sto-3g", "WATER", pot, err, &
                              scf_in=scf, max_iter_in=max_iter, grad_tol_in=grad_tol, &
                              energy_tol=energy_tol)
   end subroutine water

   subroutine test_baseline(error)
      !! The control. Every case below asserts a *failure*, so this is what
      !! says the failures mean something rather than the geometry being wrong.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(error_t) :: err

      call water(pot, err)
      call check(error,.not. err%has_error(), &
                 "the baseline potential must build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, pot%n_atoms == 3, "three atoms")
   end subroutine test_baseline

   subroutine test_maxiter(error)
      !! `maxiter` arrives
      !!
      !! Two iterations cannot converge this SCF, so a run that completes is a
      !! run where the cap never reached the loop -- which is what the
      !! hardcoded 200 in the call site did to every deck.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(error_t) :: err

      call water(pot, err, max_iter=2)
      call check(error, err%has_error(), &
                 "maxiter=2 must not converge; if it did, the cap never arrived")
      if (allocated(error)) return
      call err%clear()
   end subroutine test_maxiter

   subroutine test_diis_off(error)
      !! `diis` arrives, spelled as a subspace of zero
      !!
      !! The deck says `"diis": false` and the SCF says `diis_vectors = 0`; the
      !! translation happens on the way in. Unaccelerated, this SCF needs far
      !! more than twelve iterations, so completing within twelve would mean
      !! DIIS was still on.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(error_t) :: err
      type(scf_numerics_t) :: scf

      ! The control first: twelve iterations WITH DIIS is comfortable.
      scf%use_diis = .true.
      call water(pot, err, scf=scf, max_iter=12)
      call check(error,.not. err%has_error(), &
                 "twelve iterations with DIIS should converge: "//err%get_full_trace())
      if (allocated(error)) return

      scf%use_diis = .false.
      call water(pot, err, scf=scf, max_iter=12)
      call check(error, err%has_error(), &
                 "with DIIS off the same budget must not converge; if it did, "// &
                 "use_diis never arrived")
      if (allocated(error)) return
      call err%clear()
   end subroutine test_diis_off

   subroutine test_bad_accelerator(error)
      !! `accelerator` arrives, and a misspelling is refused rather than ignored
      !!
      !! Refused before any integral is computed, and the message names the key
      !! the user wrote. Silently falling back to DIIS would be the same class
      !! of bug as dropping the setting.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(error_t) :: err
      type(scf_numerics_t) :: scf

      scf%accelerator = "wibble"
      call water(pot, err, scf=scf)
      call check(error, err%has_error(), "a bad accelerator name must be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "wibble") > 0, &
                 "the refusal must name the spelling the deck used, got: "// &
                 err%get_message())
      if (allocated(error)) return
      call err%clear()
   end subroutine test_bad_accelerator

   subroutine test_level_shift(error)
      !! `level_shift` arrives
      !!
      !! The invariance case below cannot guard this one: a shift changes the
      !! route and not the answer, so a shift that never arrived produces the
      !! same multipoles. What does distinguish them is cost -- a shift of two
      !! hartree slows this SCF well past a budget it clears comfortably
      !! unshifted, so the run must fail if the shift arrived.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(error_t) :: err
      type(scf_numerics_t) :: scf

      call water(pot, err, max_iter=10)
      call check(error,.not. err%has_error(), &
                 "ten iterations unshifted should converge: "//err%get_full_trace())
      if (allocated(error)) return

      scf%level_shift = 2.0_dp
      call water(pot, err, scf=scf, max_iter=10)
      call check(error, err%has_error(), &
                 "a two-hartree shift must not converge in ten iterations; if it "// &
                 "did, level_shift never arrived")
      if (allocated(error)) return
      call err%clear()
   end subroutine test_level_shift

   subroutine test_diis_size(error)
      !! `diis_size` arrives
      !!
      !! A subspace of one is DIIS with nothing to extrapolate from, which
      !! converges like no acceleration at all. Same argument as the case
      !! above: the converged answer is identical either way, so only the cost
      !! separates arrival from silence.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(error_t) :: err
      type(scf_numerics_t) :: scf

      scf%diis_size = 1
      call water(pot, err, scf=scf, max_iter=12)
      call check(error, err%has_error(), &
                 "a DIIS subspace of one must not converge in twelve iterations; "// &
                 "if it did, diis_size never arrived")
      if (allocated(error)) return
      call err%clear()
   end subroutine test_diis_size

   subroutine test_lindep(error)
      !! `linear_dependence` arrives
      !!
      !! This one IS visible in the answer, which makes it the strongest of the
      !! set. A threshold above the overlap's smallest eigenvalue drops basis
      !! functions, and a potential fitted in a smaller space has different
      !! multipoles. Measured rather than assumed: STO-3G on this geometry has
      !! its smallest overlap eigenvalue at 3.44e-1, so 4e-1 removes something
      !! and the default 0 removes nothing. An earlier version of this case
      !! guessed 7e-2 and passed for the wrong reason -- it dropped nothing.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: base, dropped
      type(error_t) :: err
      type(scf_numerics_t) :: scf

      call water(base, err)
      call check(error,.not. err%has_error(), &
                 "the reference run must build: "//err%get_full_trace())
      if (allocated(error)) return

      scf%linear_dependence = 4.0e-1_dp
      call water(dropped, err, scf=scf)
      ! Either it dropped functions and moved the multipoles, or it refused --
      ! both are evidence the threshold arrived. Silence would be identical
      ! dipoles and no error.
      if (err%has_error()) then
         call err%clear()
         return
      end if
      call check(error, maxval(abs(dropped%dipole - base%dipole)) > DIPOLE_TOL, &
                 "a threshold above the smallest overlap eigenvalue must change "// &
                 "the fitted multipoles; if it did not, linear_dependence never arrived")
   end subroutine test_lindep

   subroutine test_grad_tol(error)
      !! The commutator threshold, and which of the three rules picked it
      !!
      !! `make_efp_potential` chooses between three: a stated
      !! `gradient_tolerance` wins outright, a stated `density_tolerance` sets
      !! it, and a lone `tolerance` derives it.
      !!
      !! **The energy tolerance has to be loosened for this to be visible.**
      !! Both gates must be met, so with MAKEFP's own 1e-10 energy default that
      !! one binds and the commutator never decides anything -- a first version
      !! of this case loosened only the commutator and saw seven iterations
      !! either way. Loosening the energy to 1e-4 hands the decision to the
      !! commutator, and then 1e-2 against 1e-8 stops in different places and
      !! fits different multipoles.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: tight, loose
      type(error_t) :: err

      call water(tight, err, energy_tol=1.0e-4_dp, grad_tol=1.0e-8_dp)
      call check(error,.not. err%has_error(), &
                 "the tight-commutator run must build: "//err%get_full_trace())
      if (allocated(error)) return

      call water(loose, err, energy_tol=1.0e-4_dp, grad_tol=1.0e-2_dp)
      call check(error,.not. err%has_error(), &
                 "the loose run must build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(loose%dipole - tight%dipole)) > DIPOLE_TOL, &
                 "a commutator gate of 1e-2 must not fit the same multipoles as "// &
                 "1e-8; if it did, gradient_tolerance never arrived")
   end subroutine test_grad_tol

   subroutine test_invariance(error)
      !! The settings that only change the route must not move the multipoles
      !!
      !! An accelerator, a DIIS subspace and a level shift change how the
      !! density is reached, not which density it is. The dipoles are fitted to
      !! that density, so they are the place a mis-wiring would show: a setting
      !! that arrived but was connected to the wrong argument would converge to
      !! a different answer rather than the same one by a different path.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: base, other
      type(error_t) :: err
      type(scf_numerics_t) :: scf
      real(dp) :: worst

      call water(base, err)
      call check(error,.not. err%has_error(), &
                 "the reference run must build: "//err%get_full_trace())
      if (allocated(error)) return

      scf%accelerator = "ediis"
      scf%diis_size = 12
      scf%level_shift = 0.2_dp
      scf%incremental_fock = .false.
      call water(other, err, scf=scf)
      call check(error,.not. err%has_error(), &
                 "the re-routed run must build: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, other%n_points == base%n_points, &
                 "the expansion points must not depend on how the SCF converged")
      if (allocated(error)) return

      worst = maxval(abs(other%dipole - base%dipole))
      call check(error, worst < DIPOLE_TOL, &
                 "a different route to the same density moved the dipoles")
   end subroutine test_invariance

end module test_mqc_makefp_scf

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_makefp_scf, only: collect_makefp_scf
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_makefp_scf", collect_makefp_scf)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
