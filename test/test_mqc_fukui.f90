!! Where a molecule reacts, by difference in the electron count
module test_mqc_fukui
   !! The condensed Fukui indices are differences of atomic charges between the
   !! neutral molecule and its two ions, so what can be asserted without a
   !! reference number is quite a lot: they sum to one because exactly one
   !! electron moved, the ionisation potential is positive because removing an
   !! electron costs energy, and water's oxygen is the nucleophilic site because
   !! that is where the lone pairs are.
   !!
   !! Reference values are deliberately not used. Condensed Fukui indices depend
   !! on the population scheme and on the basis, so a stored number would pin
   !! this implementation to one combination and say nothing about whether the
   !! quantity is right.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_pcm, only: pcm_context_t
   use mqc_method_config, only: pcm_config_t
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_fukui, only: fukui_result_t, fukui_indices
   implicit none
   private

   public :: collect_mqc_fukui_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   !> Bohr, C2v, with the two hydrogens exactly equivalent
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])

contains

   subroutine collect_mqc_fukui_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("indices_sum_to_one_electron", indices_sum_to_one_electron), &
                  new_unittest("water_is_nucleophilic_at_oxygen", water_is_nucleophilic_at_oxygen), &
                  new_unittest("refuses_an_open_shell_neutral", refuses_an_open_shell_neutral), &
                  new_unittest("values_have_not_drifted", values_have_not_drifted), &
                  new_unittest("the_ions_feel_the_solvent", the_ions_feel_the_solvent), &
                  new_unittest("the_global_descriptors_agree", the_global_descriptors_agree), &
                  new_unittest("the_ions_start_from_the_neutral", the_ions_start_from_the_neutral) &
                  ]
   end subroutine collect_mqc_fukui_tests

   subroutine the_ions_start_from_the_neutral(error)
      !! The neutral's orbitals are a guess for its own ions, and get used
      !!
      !! A Fukui run solves the neutral first and then two ions that differ from
      !! it by one electron. The converged neutral is therefore the best guess
      !! available for both, and it is already in the routine's hands. Before
      !! this it was used only to condense charges: the ions fell through to the
      !! unrestricted default, Wolfsberg-Helmholz, which carries no two-electron
      !! information at all and starts the anion nearly two Hartree above its own
      !! answer.
      !!
      !! Pinned as a comparison rather than against a fixed iteration count,
      !! since the count depends on the basis, the tolerances and DIIS, and a
      !! magic number here would be re-tuned rather than believed. What must
      !! hold is the direction: seeded is cheaper, and never lands higher.
      type(error_type), allocatable, intent(out) :: error
      type(fukui_result_t) :: bare, seeded
      type(error_t) :: err
      logical :: ok

      call water_fukui("mulliken", bare, err, ok)
      call check(error, ok, "the unseeded Fukui run failed")
      if (allocated(error)) return
      call water_fukui("mulliken", seeded, err, ok, seed=.true.)
      call check(error, ok, "the seeded Fukui run failed")
      if (allocated(error)) return

      call check(error, seeded%iterations_anion < bare%iterations_anion, &
                 "seeding the anion from the neutral did not save an iteration")
      if (allocated(error)) return
      call check(error, seeded%iterations_cation <= bare%iterations_cation, &
                 "seeding made the cation worse")
      if (allocated(error)) return

      ! Never *worse*, rather than identical. On water/6-31g the two runs agree
      ! to 1e-12, but this test is in STO-3G, where the anion is a minimal-basis
      ! anion with nowhere to put the extra electron: it has more than one SCF
      ! solution and the guess decides which one is found. Seeding finds the
      ! lower of the two by 0.13 Hartree. That is the guess doing its job -- a
      ! variational method that lands lower has landed better -- but it means
      ! "the answer does not move" is the wrong invariant to pin. What must hold
      ! is that a better starting point never buys a worse stationary point.
      call check(error, seeded%energy_anion <= bare%energy_anion + 1.0e-10_dp, &
                 "the seeded anion converged above the unseeded one")
      if (allocated(error)) return
      call check(error, seeded%energy_cation <= bare%energy_cation + 1.0e-10_dp, &
                 "the seeded cation converged above the unseeded one")
   end subroutine the_ions_start_from_the_neutral

   subroutine the_ions_feel_the_solvent(error)
      !! A continuum reaches all three states, not only the neutral
      !!
      !! The ions used to be refused alongside a solvated neutral, and the
      !! refusal was right: solving them in the gas phase and differencing
      !! against a solvated neutral is an ionisation potential across two
      !! different physics, wrong by the solvation energy of a charged species.
      !! That is not a correction-sized error -- on this molecule it is 3.6 eV
      !! on the potential and 3.1 eV on the affinity, which is the size of the
      !! answers themselves.
      !!
      !! So the check is not that the solvated numbers are *right* -- there is
      !! no reference for a Fukui index -- but that they are *different*, and
      !! different in the direction a dielectric must push them. A continuum
      !! that quietly failed to reach the ions would reproduce the gas-phase
      !! numbers exactly, and nothing else in this file would notice.
      type(error_type), allocatable, intent(out) :: error

      type(fukui_result_t) :: gas, solvated
      type(error_t) :: err
      logical :: ok

      call water_fukui("chelpg", gas, err, ok)
      call check(error, ok, "the gas-phase analysis did not run: "//err%get_message())
      if (allocated(error)) return

      call water_fukui("chelpg", solvated, err, ok, solvate=.true.)
      call check(error, ok, "the solvated analysis did not run: "//err%get_message())
      if (allocated(error)) return

      ! Far more than any tolerance: an unsolvated ion differs from a solvated
      ! one by electronvolts, so a hundredth of a Hartree is a floor and not a
      ! threshold anyone has to tune.
      call check(error, abs(solvated%ionisation_potential - gas%ionisation_potential) > 1.0e-2_dp, &
                 "the continuum did not reach the cation -- the ionisation "// &
                 "potential is the gas-phase one")
      if (allocated(error)) return
      call check(error, abs(solvated%electron_affinity - gas%electron_affinity) > 1.0e-2_dp, &
                 "the continuum did not reach the anion -- the electron "// &
                 "affinity is the gas-phase one")
      if (allocated(error)) return

      ! Both ions are stabilised by the dielectric, so removing an electron
      ! costs less and adding one pays better. Either moving the wrong way
      ! would mean the solvent operator reached the ions with the wrong sign.
      call check(error, solvated%ionisation_potential < gas%ionisation_potential, &
                 "a solvated cation should make ionisation cheaper")
      if (allocated(error)) return
      call check(error, solvated%electron_affinity > gas%electron_affinity, &
                 "a solvated anion should make the electron affinity less negative")
   end subroutine the_ions_feel_the_solvent

   subroutine the_global_descriptors_agree(error)
      !! The three global descriptors are one definition, not three numbers
      !!
      !! `mu = -(IP + EA)/2` and `eta = (IP - EA)/2` are the first and second
      !! derivatives of the energy with respect to electron number, and
      !! `omega = mu^2 / 2 eta` is written for exactly those. The hardness once
      !! lost its half while the chemical potential kept its own, which left the
      !! hardness twice what it should be and the electrophilicity half -- and
      !! nothing caught it, because each number is plausible on its own and no
      !! reference for them is checked anywhere.
      !!
      !! So what is pinned here is the relation between them rather than any one
      !! value: recompute all three from the two energies the ions gave and
      !! require agreement. That survives a change of molecule, basis or scheme,
      !! and it is the only check that would have failed when the half went
      !! missing.
      type(error_type), allocatable, intent(out) :: error

      type(fukui_result_t) :: res
      type(error_t) :: err
      logical :: ok
      real(dp) :: mu, eta

      call water_fukui("chelpg", res, err, ok)
      call check(error, ok, "the analysis did not run: "//err%get_message())
      if (allocated(error)) return

      mu = -0.5_dp*(res%ionisation_potential + res%electron_affinity)
      eta = 0.5_dp*(res%ionisation_potential - res%electron_affinity)

      call check(error, res%chemical_potential, mu, thr=1.0e-12_dp, &
                 more="the chemical potential is not -(IP + EA)/2")
      if (allocated(error)) return
      call check(error, res%hardness, eta, thr=1.0e-12_dp, &
                 more="the hardness is not (IP - EA)/2 -- Parr and Pearson's half")
      if (allocated(error)) return
      call check(error, res%electrophilicity, mu**2/(2.0_dp*eta), thr=1.0e-12_dp, &
                 more="the electrophilicity is not mu^2 / 2 eta")
   end subroutine the_global_descriptors_agree

   subroutine water_fukui(scheme, res, err, ok, solvate, seed)
      !! Converge water and run the analysis on it
      character(len=*), intent(in) :: scheme
      type(fukui_result_t), intent(out) :: res
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok
      logical, intent(in), optional :: solvate
         !! Build a C-PCM water continuum and hand it to all three states.
      logical, intent(in), optional :: seed
         !! Hand the neutral's orbitals over, so the ions start from them.
         !! Absent is the old behaviour, where they fall back to GWH -- which
         !! is what makes the two runs comparable in the test below.

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(pcm_context_t) :: pcm_ctx
      type(pcm_config_t) :: pcm_cfg
      real(dp), allocatable :: orbs(:, :)
      logical :: wet

      ok = .false.
      wet = .false.
      if (present(solvate)) wet = solvate

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      if (err%has_error()) return

      if (wet) then
         pcm_cfg%enabled = .true.
         pcm_cfg%method = "cpcm"
         pcm_cfg%dielectric = 78.3553_dp
         pcm_cfg%angular_points = 110
         call pcm_ctx%build(mol, WATER_Z, pcm_cfg, err)
         if (err%has_error()) then
            call mol%destroy()
            return
         end if
      end if

      ! The neutral is solved in the same continuum the ions will be, which is
      ! the whole point: the three energies have to come from one Hamiltonian.
      call run_libcint_rhf(mol, 10, 100, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err, &
                           pcm=pcm_ctx)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if
      ! `orbs` stays unallocated when not seeding, and an unallocated
      ! allocatable is absent at an optional dummy -- so one call site serves
      ! both the seeded path and the bare one.
      if (present(seed)) then
         if (seed) orbs = scf%orbitals
      end if
      call fukui_indices(mol, 10, 1, scf%density, scf%energy, scheme, 100, &
                         1.0e-10_dp, 1.0e-8_dp, res, err, pcm=pcm_ctx, &
                         neutral_orbitals=orbs)
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine water_fukui

   subroutine indices_sum_to_one_electron(error)
      !! Exactly one electron was added, and exactly one removed
      !!
      !! The population scheme divides the whole density, so the difference
      !! between two states differing by one electron has to account for one
      !! electron. That holds for either scheme and for any molecule, and it is
      !! what catches an ion run at the wrong charge, a density that was only
      !! one spin, or a fit that lost part of the molecule.
      type(error_type), allocatable, intent(out) :: error

      type(fukui_result_t) :: res
      type(error_t) :: err
      logical :: ok

      call water_fukui("chelpg", res, err, ok)
      call check(error, ok, "the analysis did not run: "//err%get_message())
      if (allocated(error)) return

      call check(error, sum(res%f_plus), 1.0_dp, thr=1.0e-6_dp, &
                 more="f+ does not account for the electron that was added")
      if (allocated(error)) return
      call check(error, sum(res%f_minus), 1.0_dp, thr=1.0e-6_dp, &
                 more="f- does not account for the electron that was removed")
      if (allocated(error)) return
      ! The dual descriptor is a difference of two things that each sum to one.
      call check(error, sum(res%dual), 0.0_dp, thr=1.0e-6_dp, &
                 more="the dual descriptor should sum to zero")
   end subroutine indices_sum_to_one_electron

   subroutine water_is_nucleophilic_at_oxygen(error)
      !! The chemistry, which is the only thing worth checking without a reference
      !!
      !! Water gives up charge from the oxygen lone pairs and accepts it at the
      !! hydrogens, so the dual descriptor is negative on the oxygen and
      !! positive on both hydrogens. Any implementation that had `f+` and `f-`
      !! the wrong way round would reproduce every sum rule above and get this
      !! exactly backwards.
      !!
      !! The two hydrogens are symmetry-equivalent in this geometry, so they
      !! must agree -- loosely, since the electrostatic fit samples a grid that
      !! does not share the molecule's symmetry.
      type(error_type), allocatable, intent(out) :: error

      type(fukui_result_t) :: res
      type(error_t) :: err
      logical :: ok

      call water_fukui("chelpg", res, err, ok)
      call check(error, ok, "the analysis did not run: "//err%get_message())
      if (allocated(error)) return

      call check(error, res%dual(1) < 0.0_dp, &
                 "water's oxygen should be the nucleophilic site")
      if (allocated(error)) return
      call check(error, res%dual(2) > 0.0_dp .and. res%dual(3) > 0.0_dp, &
                 "water's hydrogens should be the electrophilic sites")
      if (allocated(error)) return
      call check(error, res%f_minus(1) > res%f_minus(2), &
                 "the oxygen should give up charge more readily than a hydrogen")
      if (allocated(error)) return
      call check(error, abs(res%dual(2) - res%dual(3)) < 0.05_dp, &
                 "the two equivalent hydrogens disagree")
      if (allocated(error)) return

      ! Removing an electron costs energy. Nothing else here would notice the
      ! two ion energies being swapped.
      call check(error, res%ionisation_potential > 0.0_dp, &
                 "the ionisation potential came out negative")
   end subroutine water_is_nucleophilic_at_oxygen

   subroutine refuses_an_open_shell_neutral(error)
      !! An open-shell neutral is refused rather than guessed at
      !!
      !! Both ions of a closed shell are doublets. For anything else the two
      !! multiplicities are a chemical choice, and picking one silently would
      !! produce numbers that look exactly like the ones above.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(fukui_result_t) :: res
      type(error_t) :: err

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call run_libcint_rhf(mol, 10, 100, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
      call check(error,.not. err%has_error(), "could not converge water")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call fukui_indices(mol, 10, 3, scf%density, scf%energy, "chelpg", 100, &
                         1.0e-10_dp, 1.0e-8_dp, res, err)
      call check(error, err%has_error(), "an open-shell neutral was accepted")
      call mol%destroy()
   end subroutine refuses_an_open_shell_neutral

   subroutine values_have_not_drifted(error)
      !! A regression pin, and not a claim that these are the right numbers
      !!
      !! **These are this code's own output, not a reference.** PySCF has no
      !! Fukui indices, so there is nothing external to check against, and the
      !! values depend on the population scheme and the basis besides -- water
      !! in 6-31G gives a quite different set. What they catch is drift: a
      !! change to the charge fitting, to the ion SCFs, or to the differencing
      !! that moves the answer without breaking the sum rule or reversing the
      !! chemistry, both of which the tests above already pin.
      !!
      !! Loose on purpose. The electrostatic fit samples a grid, so the last
      !! couple of figures move with anything that perturbs it; a tolerance
      !! tight enough to catch that would fail on a compiler change and teach
      !! everyone to update the number without looking.
      type(error_type), allocatable, intent(out) :: error

      type(fukui_result_t) :: res
      type(error_t) :: err
      logical :: ok
      real(dp), parameter :: PIN = 1.0e-4_dp

      call water_fukui("chelpg", res, err, ok)
      call check(error, ok, "the analysis did not run: "//err%get_message())
      if (allocated(error)) return

      call check(error, res%f_plus(1), 0.141096_dp, thr=PIN, more="O f+ has drifted")
      if (allocated(error)) return
      call check(error, res%f_minus(1), 0.562489_dp, thr=PIN, more="O f- has drifted")
      if (allocated(error)) return
      call check(error, res%f_plus(2), 0.429599_dp, thr=PIN, more="H f+ has drifted")
      if (allocated(error)) return
      call check(error, res%f_minus(2), 0.218800_dp, thr=PIN, more="H f- has drifted")
      if (allocated(error)) return
      call check(error, res%ionisation_potential, 0.307154_dp, thr=PIN, &
                 more="the ionisation potential has drifted")
      if (allocated(error)) return
      call check(error, res%electron_affinity, -0.713407_dp, thr=PIN, &
                 more="the electron affinity has drifted")
   end subroutine values_have_not_drifted

end module test_mqc_fukui

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fukui, only: collect_mqc_fukui_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fukui", collect_mqc_fukui_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
