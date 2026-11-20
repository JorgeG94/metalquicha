!! Geometric direct minimisation as the SCF's second phase
module test_mqc_czt_gdm
   !! `run_gdm_phase` end to end. The numerics it is built from -- the L-BFGS
   !! history, the preconditioner, the accept rule -- are held to dense linear
   !! algebra in `test_mqc_scf_gdm`; what is left to check here is that the
   !! phase wired around them minimises the right function.
   !!
   !! **The energy DIIS reaches is the reference.** Both are the same
   !! stationary point of the same functional on water, so the two energies
   !! agree to what the convergence threshold leaves, and a gradient with the
   !! wrong sign or factor, a rotation applied in the wrong frame or a
   !! preconditioner indexed the other way round does not get there: it stalls
   !! or lands somewhere else.
   !!
   !! **Converging is not enough to show the phase ran.** A request that
   !! quietly fell back to DIIS would converge too, to the same energy, so
   !! every case also reads `gdm_started_at` and `gdm_iterations`.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_scf_convergence, only: scf_convergence_t, CONV_METRIC_COMMUTATOR
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf, SCF_GUESS_CORE
   use mqc_diis, only: parse_accelerator_name, ACCEL_DIIS
   use mqc_scf_types, only: scf_numerics_t
   implicit none
   private

   public :: collect_mqc_czt_gdm_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   !! Hartree. Both runs stop on a commutator of 1e-9, and the energy's error
   !! is second order in that.
   real(dp), parameter :: TOL_ENERGY = 1.0e-10_dp
   !! A handover threshold no commutator fails, so the minimisation starts at
   !! the first opportunity.
   real(dp), parameter :: START_AT_ONCE = 1.0e9_dp

contains

   subroutine collect_mqc_czt_gdm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("gdm_lands_where_diis_lands", test_same_energy_as_diis), &
                  new_unittest("gdm_as_the_accelerator_lands_there_too", &
                               test_accelerator_route), &
                  new_unittest("gdm_from_a_core_guess_still_gets_there", test_early_start), &
                  new_unittest("gdm_leaves_canonical_orbitals", test_orbitals), &
                  new_unittest("the_accelerator_name_needs_a_caller_that_can_run_it", &
                               test_accelerator_name), &
                  new_unittest("gdm_and_second_order_together_are_refused", &
                               test_refuses_both), &
                  new_unittest("gdm_refuses_a_continuum_solvent", test_refuses_pcm), &
                  new_unittest("unrestricted_gdm_on_a_closed_shell_is_restricted_gdm", &
                               test_unrestricted_closed_shell), &
                  new_unittest("unrestricted_gdm_lands_where_diis_lands", &
                               test_unrestricted_same_energy), &
                  new_unittest("unrestricted_gdm_started_at_once_still_gets_there", &
                               test_unrestricted_early_start) &
                  ]
   end subroutine collect_mqc_czt_gdm_tests

   subroutine water_molecule(mol, err)
      !! Water in 6-31G, the geometry the second-order tests use
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
   end subroutine water_molecule

   subroutine water(mol, scf, err, gdm, accelerator, gdm_start, guess)
      !! Converge it, tightly and on the commutator alone, so two routes'
      !! energies are comparable at one threshold
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err
      logical, intent(in), optional :: gdm
      character(len=*), intent(in), optional :: accelerator
         !! Routed through an `scf_numerics_t`, which is how a deck reaches
         !! the SCF.
      real(dp), intent(in), optional :: gdm_start
      integer, intent(in), optional :: guess

      type(scf_convergence_t) :: conv
      type(scf_numerics_t) :: settings

      call water_molecule(mol, err)
      if (err%has_error()) return
      conv%metric = CONV_METRIC_COMMUTATOR
      conv%tolerance = 1.0e-9_dp
      if (present(accelerator)) then
         settings%accelerator = accelerator
         call run_czt_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                          in_core=.true., convergence=conv, scf=settings)
      else
         call run_czt_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                          in_core=.true., convergence=conv, gdm=gdm, &
                          gdm_start=gdm_start, guess=guess)
      end if
   end subroutine water

   subroutine water_unrestricted(mol, scf, err, nelec, multiplicity, gdm, gdm_start)
      !! The same molecule through the unrestricted driver, at the same
      !! threshold. Nine electrons as a doublet is the cation.
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err
      integer, intent(in) :: nelec, multiplicity
      logical, intent(in), optional :: gdm
      real(dp), intent(in), optional :: gdm_start

      type(scf_convergence_t) :: conv

      call water_molecule(mol, err)
      if (err%has_error()) return
      conv%metric = CONV_METRIC_COMMUTATOR
      conv%tolerance = 1.0e-9_dp
      call run_czt_uhf(mol, nelec, multiplicity, 300, 1.0e-12_dp, 1.0e-10_dp, .false., &
                       scf, err, in_core=.true., convergence=conv, gdm=gdm, &
                       gdm_start=gdm_start)
   end subroutine water_unrestricted

   subroutine test_unrestricted_closed_shell(error)
      !! Equal spins: the two-spin minimisation must find the restricted energy
      !!
      !! The per-spin gradient is `2 F_ai` where the restricted one is `4 F_ai`
      !! over half as many variables, and the preconditioner's factor follows
      !! it. A singlet is where a slip in either shows as a wrong energy or a
      !! stall rather than as "another open-shell solution".
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol_r, mol_u
      type(rhf_result_t) :: r, u
      type(error_t) :: err

      call water(mol_r, r, err)
      call check(error,.not. err%has_error() .and. r%converged, "RHF did not converge")
      if (allocated(error)) return
      call water_unrestricted(mol_u, u, err, 10, 1, gdm=.true.)
      call check(error,.not. err%has_error(), "the run failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, u%converged, "unrestricted direct minimisation did not converge")
      if (allocated(error)) return
      call check(error, u%gdm_started_at > 0 .and. u%gdm_iterations > 0, &
                 "direct minimisation never took a step")
      if (allocated(error)) return
      write (*, "(a,i0,a,i0,a,es10.2)") "    DIIS then ", u%gdm_iterations, &
         " steps, ", u%gdm_rejected, " rejected; E(ugdm) - E(rhf) = ", u%energy - r%energy
      call check(error, abs(u%energy - r%energy) < TOL_ENERGY, &
                 "unrestricted direct minimisation missed the restricted energy")
      call mol_r%destroy()
      call mol_u%destroy()
   end subroutine test_unrestricted_closed_shell

   subroutine test_unrestricted_same_energy(error)
      !! The water cation, a doublet: the same minimum as DIIS
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol_diis, mol_gdm
      type(rhf_result_t) :: diis, gdm
      type(error_t) :: err

      call water_unrestricted(mol_diis, diis, err, 9, 2)
      call check(error,.not. err%has_error() .and. diis%converged, "DIIS did not converge")
      if (allocated(error)) return
      call water_unrestricted(mol_gdm, gdm, err, 9, 2, gdm=.true.)
      call check(error,.not. err%has_error(), "the run failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, gdm%converged, "direct minimisation did not converge")
      if (allocated(error)) return
      call check(error, gdm%gdm_started_at > 0, "direct minimisation never took over")
      if (allocated(error)) return
      call check(error, gdm%gdm_iterations > 0, "the phase took no step")
      if (allocated(error)) return
      call check(error,.not. gdm%gdm_stalled, "the phase reported a stall")
      if (allocated(error)) return
      call check(error, gdm%gdm_off_aufbau == 0, "the phase ended off the aufbau occupation")
      if (allocated(error)) return
      write (*, "(a,i0,a,i0,a,es10.2)") "    DIIS then ", gdm%gdm_iterations, &
         " steps, ", gdm%gdm_rejected, " rejected; E(gdm) - E(diis) = ", &
         gdm%energy - diis%energy
      call check(error, abs(gdm%energy - diis%energy) < TOL_ENERGY, &
                 "direct minimisation and DIIS reached different energies")
      if (allocated(error)) return
      ! The orbital energies must be the canonical ones, in both spins: the
      ! phase semicanonicalises once at the end and reports those.
      call check(error, maxval(abs(gdm%orbital_energies - diis%orbital_energies)) < 1.0e-6_dp, &
                 "the alpha orbital energies are not the canonical ones")
      if (allocated(error)) return
      call check(error, maxval(abs(gdm%orbital_energies_beta &
                                   - diis%orbital_energies_beta)) < 1.0e-6_dp, &
                 "the beta orbital energies are not the canonical ones")
      call mol_diis%destroy()
      call mol_gdm%destroy()
   end subroutine test_unrestricted_same_energy

   subroutine test_unrestricted_early_start(error)
      !! Handed over after one DIIS iteration, far from convergence
      !!
      !! Here the history, the halving and the restart all get used, where a
      !! handover at the default threshold needs a handful of clean steps.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol_diis, mol_gdm
      type(rhf_result_t) :: diis, gdm
      type(error_t) :: err

      call water_unrestricted(mol_diis, diis, err, 9, 2)
      call check(error,.not. err%has_error() .and. diis%converged, "DIIS did not converge")
      if (allocated(error)) return
      call water_unrestricted(mol_gdm, gdm, err, 9, 2, gdm=.true., gdm_start=START_AT_ONCE)
      call check(error,.not. err%has_error(), "the run failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, gdm%gdm_started_at == 1, "the handover was not after iteration 1")
      if (allocated(error)) return
      call check(error, gdm%converged, "direct minimisation did not converge")
      if (allocated(error)) return
      write (*, "(a,i0,a,i0,a,es10.2)") "    from iteration 1: ", gdm%gdm_iterations, &
         " steps, ", gdm%gdm_rejected, " rejected; E(gdm) - E(diis) = ", &
         gdm%energy - diis%energy
      call check(error, abs(gdm%energy - diis%energy) < TOL_ENERGY, &
                 "an early start reached a different energy")
      call mol_diis%destroy()
      call mol_gdm%destroy()
   end subroutine test_unrestricted_early_start

   subroutine test_same_energy_as_diis(error)
      !! The same minimum, by a different route
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol_diis, mol_gdm
      type(rhf_result_t) :: diis, gdm
      type(error_t) :: err

      call water(mol_diis, diis, err)
      call check(error,.not. err%has_error() .and. diis%converged, "DIIS did not converge")
      if (allocated(error)) return
      call water(mol_gdm, gdm, err, gdm=.true.)
      call check(error,.not. err%has_error(), "the run failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, gdm%converged, "direct minimisation did not converge")
      if (allocated(error)) return
      call check(error, gdm%gdm_started_at > 0, "direct minimisation never took over")
      if (allocated(error)) return
      call check(error, gdm%gdm_iterations > 0, "the phase took no step")
      if (allocated(error)) return
      call check(error,.not. gdm%gdm_stalled, "the phase reported a stall")
      if (allocated(error)) return
      call check(error, gdm%gdm_off_aufbau == 0, "the phase ended off the aufbau occupation")
      if (allocated(error)) return
      write (*, "(a,i0,a,i0,a,es10.2)") "    DIIS then ", gdm%gdm_iterations, &
         " steps, ", gdm%gdm_rejected, " rejected; E(gdm) - E(diis) = ", &
         gdm%energy - diis%energy
      call check(error, abs(gdm%energy - diis%energy) < TOL_ENERGY, &
                 "direct minimisation and DIIS reached different energies")
      call mol_diis%destroy()
      call mol_gdm%destroy()
   end subroutine test_same_energy_as_diis

   subroutine test_accelerator_route(error)
      !! `accelerator: gdm` is the same run as `gdm: true`
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol_a, mol_b
      type(rhf_result_t) :: by_flag, by_name
      type(error_t) :: err

      call water(mol_a, by_flag, err, gdm=.true.)
      call water(mol_b, by_name, err, accelerator="gdm")
      call check(error,.not. err%has_error(), "a run failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, by_name%gdm_started_at > 0, "the name did not start the phase")
      if (allocated(error)) return
      call check(error, by_name%gdm_iterations == by_flag%gdm_iterations, &
                 "the name and the flag took different numbers of steps")
      if (allocated(error)) return
      call check(error, abs(by_name%energy - by_flag%energy) < TOL_ENERGY, &
                 "the name and the flag reached different energies")
      call mol_a%destroy()
      call mol_b%destroy()
   end subroutine test_accelerator_route

   subroutine test_early_start(error)
      !! Started from the core guess after one diagonalisation, where the
      !! commutator is of order one, it still reaches the DIIS energy
      !!
      !! This is the case the accept rule and the step-scale back-off exist
      !! for: the first quasi-Newton steps from there overshoot. It is also the
      !! only case here whose history fills and wraps.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol_diis, mol_gdm
      type(rhf_result_t) :: diis, gdm
      type(error_t) :: err

      call water(mol_diis, diis, err, guess=SCF_GUESS_CORE)
      call check(error,.not. err%has_error() .and. diis%converged, "DIIS did not converge")
      if (allocated(error)) return
      call water(mol_gdm, gdm, err, gdm=.true., gdm_start=START_AT_ONCE, guess=SCF_GUESS_CORE)
      call check(error,.not. err%has_error(), "the run failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, gdm%gdm_started_at == 1, "the phase should start after the "// &
                 "first diagonalisation")
      if (allocated(error)) return
      call check(error, gdm%converged, "direct minimisation did not converge")
      if (allocated(error)) return
      write (*, "(a,i0,a,i0,a,es10.2)") "    from the core guess: ", gdm%gdm_iterations, &
         " steps, ", gdm%gdm_rejected, " rejected; E(gdm) - E(diis) = ", &
         gdm%energy - diis%energy
      call check(error, abs(gdm%energy - diis%energy) < TOL_ENERGY, &
                 "started early, direct minimisation reached a different energy")
      call mol_diis%destroy()
      call mol_gdm%destroy()
   end subroutine test_early_start

   subroutine test_orbitals(error)
      !! The orbital energies that leave are DIIS's
      !!
      !! The phase never diagonalises a Fock matrix on the way; it
      !! semicanonicalises once at the end. At convergence that is the
      !! canonical set, so the eigenvalues have to be the ones a
      !! diagonalisation gives.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol_diis, mol_gdm
      type(rhf_result_t) :: diis, gdm
      type(error_t) :: err

      call water(mol_diis, diis, err)
      call water(mol_gdm, gdm, err, gdm=.true.)
      call check(error,.not. err%has_error(), "a run failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, size(gdm%orbital_energies) == size(diis%orbital_energies), &
                 "a different number of orbitals came back")
      if (allocated(error)) return
      call check(error, maxval(abs(gdm%orbital_energies - diis%orbital_energies)) &
                 < 1.0e-7_dp, "the orbital energies are not the canonical ones")
      if (allocated(error)) return
      call check(error, maxval(abs(gdm%density - diis%density)) < 1.0e-7_dp, &
                 "the density is not DIIS's")
      call mol_diis%destroy()
      call mol_gdm%destroy()
   end subroutine test_orbitals

   subroutine test_accelerator_name(error)
      !! `gdm` parses only for a caller that says it can run the phase
      type(error_type), allocatable, intent(out) :: error

      integer :: scheme
      logical :: ok, second_order, gdm

      call parse_accelerator_name("gdm", scheme, ok, second_order, gdm)
      call check(error, ok .and. gdm .and. .not. second_order .and. scheme == ACCEL_DIIS, &
                 "'gdm' should parse to a DIIS opening with the phase requested")
      if (allocated(error)) return
      call parse_accelerator_name("GDM", scheme, ok, second_order, gdm)
      call check(error, ok .and. gdm, "the name should not be case sensitive")
      if (allocated(error)) return
      call parse_accelerator_name("gdm", scheme, ok)
      call check(error,.not. ok, "a caller that cannot run the phase must not have "// &
                 "the name accepted on its behalf")
      if (allocated(error)) return
      call parse_accelerator_name("diis", scheme, ok, second_order, gdm)
      call check(error, ok .and. .not. gdm, "'diis' must not ask for the phase")
   end subroutine test_accelerator_name

   subroutine test_refuses_both(error)
      !! Two ways to finish the SCF cannot both be asked for
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err

      call water_molecule(mol, err)
      call check(error,.not. err%has_error(), "the molecule would not build")
      if (allocated(error)) return
      call run_czt_rhf(mol, 10, 20, 1.0e-8_dp, 1.0e-6_dp, .false., scf, err, &
                       in_core=.true., gdm=.true., second_order=.true.)
      call check(error, err%has_error(), "gdm with second_order should be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "second_order") > 0, &
                 "the refusal should name both keywords; it said: "//err%get_message())
      call mol%destroy()
   end subroutine test_refuses_both

   subroutine test_refuses_pcm(error)
      !! What it cannot do it declines, naming the keyword
      type(error_type), allocatable, intent(out) :: error

      use_pcm: block
         use mqc_czt_pcm, only: pcm_context_t
         type(czt_molecule_t) :: mol
         type(rhf_result_t) :: scf
         type(error_t) :: err
         type(pcm_context_t) :: pcm

         call water_molecule(mol, err)
         call check(error,.not. err%has_error(), "the molecule would not build")
         if (allocated(error)) return
         pcm%enabled = .true.
         call run_czt_rhf(mol, 10, 20, 1.0e-8_dp, 1.0e-6_dp, .false., scf, err, &
                          in_core=.true., pcm=pcm, gdm=.true.)
         call check(error, err%has_error(), "gdm in a continuum solvent should be refused")
         if (allocated(error)) return
         call check(error, index(err%get_message(), "gdm") > 0, &
                    "the refusal should name the keyword; it said: "//err%get_message())
         call mol%destroy()
      end block use_pcm
   end subroutine test_refuses_pcm

end module test_mqc_czt_gdm

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_gdm, only: collect_mqc_czt_gdm_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_gdm", collect_mqc_czt_gdm_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
