module test_mqc_libcint_mp2
   !! Pins conventional MP2 against what a reference energy cannot check.
   !!
   !! The validation manifest already compares total energies to PySCF across
   !! bases, molecules and frozen-core counts, so this suite is for the
   !! structural properties:
   !!
   !!   * the four-index transform has to preserve the integrals' permutational
   !!     symmetry, (ia|jb) == (jb|ia). An index transposed in one of the four
   !!     quarter-steps breaks that while still producing a plausible energy;
   !!   * the correlation energy has to be negative for a stable closed-shell
   !!     reference, and the two spin components separately so;
   !!   * freezing orbitals has to raise it, monotonically, and by roughly what
   !!     a core contributes rather than by an arbitrary amount;
   !!   * the spin scaling has to be applied by `total` and not folded into the
   !!     components, so a scaled result still reports what it scaled;
   !!   * the frozen-core *gradient* has to difference the frozen-core energy,
   !!     match `pyscf.grad.mp2` element by element, and come out the same from
   !!     every memory path -- the occupied-frozen block of the relaxed density
   !!     is built from the Lagrangian, and no energy comparison can see it.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_result_types, only: mp2_energy_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_mp2, run_libcint_ri_mp2
   use mqc_libcint_mp2_gradient, only: libcint_mp2_gradient
   use mqc_libcint_ri_mp2_gradient, only: libcint_ri_mp2_gradient
   implicit none
   private
   public :: collect_mqc_libcint_mp2_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: E_TOL = 1.0e-11_dp
   real(dp), parameter :: D_TOL = 1.0e-9_dp

   ! The gradient tests' geometry, in Bohr already: the same water the PySCF
   ! reference below was computed at, so the pinned numbers mean what they say.
   real(dp), parameter :: WATER_B(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])

contains

   subroutine collect_mqc_libcint_mp2_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mp2_correlation_is_negative", test_negative), &
                  new_unittest("mp2_components_sum_to_the_correlation", test_sum), &
                  new_unittest("mp2_freezing_core_raises_the_energy", test_frozen), &
                  new_unittest("mp2_rejects_freezing_everything", test_bad_frozen), &
                  new_unittest("mp2_scaling_is_applied_by_total_only", test_scaling), &
                  new_unittest("mp2_scs_matches_its_published_factors", test_scs), &
                  new_unittest("ri_mp2_approaches_the_conventional_answer", test_ri), &
                  new_unittest("mp2_frozen_gradient_matches_pyscf", &
                               test_frozen_gradient_pyscf), &
                  new_unittest("mp2_frozen_gradient_differences_the_energy", &
                               test_frozen_gradient_fd), &
                  new_unittest("mp2_frozen_gradient_memory_paths_agree", &
                               test_frozen_gradient_paths), &
                  new_unittest("mp2_gradient_rejects_freezing_everything", &
                               test_frozen_gradient_refusal), &
                  new_unittest("ri_mp2_frozen_gradient_matches_the_numpy_reference", &
                               test_ri_frozen_gradient_reference), &
                  new_unittest("ri_mp2_frozen_gradient_differences_the_energy", &
                               test_ri_frozen_gradient_fd), &
                  new_unittest("ri_mp2_frozen_gradient_integral_paths_agree", &
                               test_ri_frozen_gradient_paths), &
                  new_unittest("ri_mp2_gradient_rejects_freezing_everything", &
                               test_ri_frozen_gradient_refusal) &
                  ]
   end subroutine collect_mqc_libcint_mp2_tests

   subroutine water_mp2(basis, frozen, mp2, err)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: frozen
      type(mp2_result_t), intent(out) :: mp2
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.757_dp*ANG, 0.587_dp*ANG, &
                   0.0_dp, -0.757_dp*ANG, 0.587_dp*ANG], [3, 3])
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 200, E_TOL, D_TOL, .false., scf, err)
      if (err%has_error()) return
      call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, 5, scf%energy, &
                           mp2, err, n_frozen=frozen)
      call mol%destroy()
   end subroutine water_mp2

   subroutine test_negative(error)
      !! Both spin components are negative for a stable reference
      !!
      !! Every denominator is e_occ - e_virt, which is negative below the gap,
      !! and both numerators are squares or near-squares. A positive component
      !! means either the reference is not a minimum or a sign is wrong, and
      !! the two are worth being able to tell apart.
      type(error_type), allocatable, intent(out) :: error
      type(mp2_result_t) :: mp2
      type(error_t) :: err

      call water_mp2("cc-pvdz", 0, mp2, err)
      call check(error,.not. err%has_error(), "water MP2 must run: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, mp2%opposite_spin < 0.0_dp, "E_OS must be negative")
      if (allocated(error)) return
      call check(error, mp2%same_spin < 0.0_dp, "E_SS must be negative")
      if (allocated(error)) return
      call check(error, mp2%total < mp2%scf_energy, &
                 "correlation must lower the reference energy")
   end subroutine test_negative

   subroutine test_sum(error)
      !! The components add up, and opposite-spin dominates
      !!
      !! The ratio is the useful half. E_OS is around three times E_SS for a
      !! closed-shell molecule, and a transform that mixed the two -- swapping
      !! (ib|ja) for (ia|jb) in the same-spin term, say -- would still sum to
      !! something negative and still look like a correlation energy.
      type(error_type), allocatable, intent(out) :: error
      type(mp2_result_t) :: mp2
      type(error_t) :: err

      call water_mp2("cc-pvdz", 0, mp2, err)
      call check(error, abs(mp2%correlation - (mp2%same_spin + mp2%opposite_spin)) &
                 < 1.0e-14_dp, "the components must sum to the correlation energy")
      if (allocated(error)) return
      call check(error, abs(mp2%total - (mp2%scf_energy + mp2%correlation)) < 1.0e-14_dp, &
                 "the total must be the reference plus the correlation")
      if (allocated(error)) return
      call check(error, mp2%opposite_spin < mp2%same_spin, &
                 "opposite-spin correlation must be the larger of the two")
      if (allocated(error)) return
      call check(error, mp2%opposite_spin/mp2%same_spin > 2.0_dp .and. &
                 mp2%opposite_spin/mp2%same_spin < 4.0_dp, &
                 "the OS/SS ratio should be around three for a closed shell")
   end subroutine test_sum

   subroutine test_frozen(error)
      !! Correlating fewer orbitals recovers less correlation
      type(error_type), allocatable, intent(out) :: error
      type(mp2_result_t) :: all_e, frozen
      type(error_t) :: err

      call water_mp2("cc-pvdz", 0, all_e, err)
      call water_mp2("cc-pvdz", 1, frozen, err)
      call check(error,.not. err%has_error(), "both must run: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, frozen%n_occupied, all_e%n_occupied - 1)
      if (allocated(error)) return
      call check(error, frozen%n_virtual, all_e%n_virtual, &
                 "freezing the core must not change the virtual space")
      if (allocated(error)) return
      call check(error, frozen%correlation > all_e%correlation, &
                 "a frozen core must recover less correlation")
      if (allocated(error)) return
      ! Oxygen's 1s is worth a couple of millihartree here, not tens.
      call check(error, abs(frozen%correlation - all_e%correlation) < 0.01_dp, &
                 "freezing one core orbital should cost only a few millihartree")
   end subroutine test_frozen

   subroutine test_bad_frozen(error)
      !! Freezing every occupied orbital leaves nothing to correlate
      type(error_type), allocatable, intent(out) :: error
      type(mp2_result_t) :: mp2
      type(error_t) :: err

      call water_mp2("cc-pvdz", 5, mp2, err)
      call check(error, err%has_error(), "freezing all five occupied must be refused")
      if (allocated(error)) return

      call err%clear()
      call water_mp2("cc-pvdz", -1, mp2, err)
      call check(error, err%has_error(), "a negative frozen count must be refused")
   end subroutine test_bad_frozen

   subroutine test_scaling(error)
      !! The scaling belongs to `total`, not to the components
      !!
      !! Folding the factors into ss and os would give the same total and lose
      !! the numbers it was scaled from, so a result could never be rescaled or
      !! compared against an unscaled one. This is what keeps them separate.
      type(error_type), allocatable, intent(out) :: error
      type(mp2_energy_t) :: e

      e%ss = -0.05_dp
      e%os = -0.15_dp

      ! Default is plain MP2: both factors one.
      call check(error, abs(e%total() - (-0.20_dp)) < 1.0e-14_dp, &
                 "unscaled, the total is the plain sum")
      if (allocated(error)) return

      e%ss_scale = 1.0_dp/3.0_dp
      e%os_scale = 1.2_dp
      call check(error, abs(e%total() - (-0.15_dp*1.2_dp - 0.05_dp/3.0_dp)) < 1.0e-14_dp, &
                 "scaled, the total applies the factors")
      if (allocated(error)) return
      call check(error, abs(e%ss - (-0.05_dp)) < 1.0e-14_dp, &
                 "the same-spin component must be untouched by scaling")
      if (allocated(error)) return
      call check(error, abs(e%os - (-0.15_dp)) < 1.0e-14_dp, &
                 "the opposite-spin component must be untouched by scaling")
   end subroutine test_scaling

   subroutine test_scs(error)
      !! `scs` reports Grimme's factors whatever the run used
      type(error_type), allocatable, intent(out) :: error
      type(mp2_energy_t) :: e

      e%ss = -0.05_dp
      e%os = -0.15_dp
      ! Deliberately set to something else: scs answers "what would SCS-MP2
      ! give", which is a different question from "what did this run give".
      e%ss_scale = 0.0_dp
      e%os_scale = 1.3_dp

      call check(error, abs(e%scs() - (-0.15_dp*1.2_dp - 0.05_dp/3.0_dp)) < 1.0e-14_dp, &
                 "scs must use the published factors, not the run's")
      if (allocated(error)) return
      call check(error, abs(e%total() - (-0.15_dp*1.3_dp)) < 1.0e-14_dp, &
                 "and total must still use the run's")
   end subroutine test_scs

   subroutine test_ri(error)
      !! The fitted answer has to be close to the exact one, and not too close
      !!
      !! Two bounds rather than one. The upper bound is the obvious check: a
      !! transform that lost an index would be wrong by millihartree, not by
      !! the tens of microhartree a fit costs. The lower bound is the one worth
      !! having -- if the fitted and conventional energies agreed to 1e-12, the
      !! fitted path would not be fitting anything, which is what a silent
      !! fallback to the four-index route would look like from outside.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: exact, ri
      type(error_t) :: err
      real(dp) :: c(3, 3), gap

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.757_dp*ANG, 0.587_dp*ANG, &
                   0.0_dp, -0.757_dp*ANG, 0.587_dp*ANG], [3, 3])
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, "cc-pvdz", mol, err)
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, "cc-pvdz-rifit", aux, err)
      call check(error,.not. err%has_error(), "water and its auxiliary must build: "//err%get_full_trace())
      if (allocated(error)) return

      call run_libcint_rhf(mol, 10, 200, E_TOL, D_TOL, .false., scf, err)
      call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, 5, scf%energy, &
                           exact, err, n_frozen=1)
      call run_libcint_ri_mp2(mol, aux, scf%orbitals, scf%orbital_energies, 5, &
                              scf%energy, ri, err, n_frozen=1)
      call check(error,.not. err%has_error(), "both routes must run: "//err%get_full_trace())
      if (allocated(error)) return

      gap = abs(ri%correlation - exact%correlation)
      call check(error, gap < 1.0e-3_dp, &
                 "the fitting error must be small; a lost index would be millihartree")
      if (allocated(error)) return
      call check(error, gap > 1.0e-9_dp, &
                 "the fitted route must actually fit; exact agreement means it did not")
      if (allocated(error)) return

      ! The components have to be fitted too, not just their sum.
      call check(error, abs(ri%opposite_spin - exact%opposite_spin) > 1.0e-9_dp, &
                 "the opposite-spin component must be fitted as well")
      if (allocated(error)) return
      call check(error, ri%n_occupied, exact%n_occupied)
      if (allocated(error)) return
      call check(error, ri%n_virtual, exact%n_virtual)
   end subroutine test_ri

   subroutine water_gradient(basis, frozen, gradient, err, block_bytes, force_blocked)
      !! One converged SCF at `WATER_B`, then the frozen-core MP2 gradient
      character(len=*), intent(in) :: basis
      integer, intent(in) :: frozen
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: err
      real(dp), intent(in), optional :: block_bytes
      logical, intent(in), optional :: force_blocked

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], WATER_B, basis, mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err)
      if (err%has_error()) return
      call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, 5, &
                                gradient, err, n_frozen=frozen, &
                                block_bytes=block_bytes, force_blocked=force_blocked)
      call mol%destroy()
   end subroutine water_gradient

   subroutine test_frozen_gradient_pyscf(error)
      !! Water/6-31G with the oxygen core frozen, against `pyscf.grad.mp2`
      !!
      !! The reference was computed from this repository's own basis JSON --
      !! PySCF's internal tables differ in the eighth decimal and would look
      !! exactly like a bug here -- and with PySCF's Z-vector solved densely
      !! rather than by its default Krylov solver, which carries ~6e-8 of
      !! unconverged residual that its `tol` does not control. Against that
      !! reference the two codes agree to 2e-11; the bound leaves room for the
      !! SCF landing differently, not for any term being wrong.
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: REFERENCE(3, 3) = reshape( &
                             [0.0_dp, 0.0_dp, 9.513625697403e-03_dp, &
                              0.0_dp, 2.247664933890e-02_dp, -4.756812848701e-03_dp, &
                              0.0_dp, -2.247664933891e-02_dp, -4.756812848701e-03_dp], &
                             [3, 3])
      real(dp), allocatable :: gradient(:, :)
      type(error_t) :: err

      call water_gradient("6-31g", 1, gradient, err)
      call check(error,.not. err%has_error(), &
                 "the frozen-core gradient must run: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(gradient - REFERENCE)) < 5.0e-10_dp, &
                 "the frozen-core MP2 gradient must match pyscf.grad.mp2")
   end subroutine test_frozen_gradient_pyscf

   subroutine test_frozen_gradient_fd(error)
      !! The frozen-core gradient differences the frozen-core energy
      !!
      !! Independent of PySCF: both numbers come from this code, so this pins
      !! "the gradient is the derivative of the energy it claims" -- with the
      !! *same* frozen count on both sides, since differencing an all-electron
      !! energy against a frozen gradient disagrees in the third decimal for a
      !! reason that has nothing to do with either. The tolerance is the
      !! central difference's own O(h^2) truncation, not the gradient's error.
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: STEP = 2.0e-3_dp
      real(dp), allocatable :: gradient(:, :)
      real(dp) :: moved(3, 3), plus, minus, worst
      type(error_t) :: err
      integer :: iatom, comp

      call water_gradient("sto-3g", 1, gradient, err)
      call check(error,.not. err%has_error(), &
                 "the frozen-core gradient must run: "//err%get_full_trace())
      if (allocated(error)) return

      worst = 0.0_dp
      do iatom = 1, 3
         do comp = 1, 3
            moved = WATER_B
            moved(comp, iatom) = WATER_B(comp, iatom) + STEP
            call frozen_energy_at(moved, plus, err)
            if (err%has_error()) exit
            moved(comp, iatom) = WATER_B(comp, iatom) - STEP
            call frozen_energy_at(moved, minus, err)
            if (err%has_error()) exit
            worst = max(worst, abs((plus - minus)/(2.0_dp*STEP) - gradient(comp, iatom)))
         end do
      end do
      call check(error,.not. err%has_error(), &
                 "the displaced energies must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, worst < 5.0e-6_dp, &
                 "the frozen-core gradient must difference the frozen-core energy")
   end subroutine test_frozen_gradient_fd

   subroutine frozen_energy_at(coords, energy, err)
      !! One frozen-core MP2 total energy at a displaced geometry
      real(dp), intent(in) :: coords(3, 3)
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2

      energy = 0.0_dp
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], coords, "sto-3g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 300, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err)
      if (err%has_error()) return
      call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, 5, scf%energy, &
                           mp2, err, n_frozen=1)
      if (err%has_error()) return
      energy = mp2%total
      call mol%destroy()
   end subroutine frozen_energy_at

   subroutine test_frozen_gradient_paths(error)
      !! Dense, dense split over `j`, and blocked agree with a core frozen
      !!
      !! The three build the same quantities from different amounts of memory,
      !! and with a frozen core every occupied loop and offset they carry runs
      !! over the *active* space -- a byte target of one forces several blocks
      !! of both the first index and the occupied one, which is what makes the
      !! offsets load-bearing. Not bitwise, for the same reason as the
      !! validation program: the blocked path screens where the dense one reads
      !! a tensor, and both land around 1e-12.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: dense(:, :), split(:, :), blocked(:, :)
      type(error_t) :: err

      call water_gradient("6-31g", 1, dense, err)
      call check(error,.not. err%has_error(), &
                 "the dense path must run: "//err%get_full_trace())
      if (allocated(error)) return
      call water_gradient("6-31g", 1, split, err, block_bytes=1.0_dp)
      call check(error,.not. err%has_error(), &
                 "the split path must run: "//err%get_full_trace())
      if (allocated(error)) return
      call water_gradient("6-31g", 1, blocked, err, block_bytes=1.0_dp, &
                          force_blocked=.true.)
      call check(error,.not. err%has_error(), &
                 "the blocked path must run: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, maxval(abs(dense - split)) < 1.0e-10_dp, &
                 "splitting the occupied loop must not move a frozen-core gradient")
      if (allocated(error)) return
      call check(error, maxval(abs(dense - blocked)) < 1.0e-10_dp, &
                 "the blocked path must agree with the dense one under a frozen core")
   end subroutine test_frozen_gradient_paths

   subroutine test_frozen_gradient_refusal(error)
      !! Freezing every occupied orbital leaves nothing to differentiate
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gradient(:, :)
      type(error_t) :: err

      call water_gradient("sto-3g", 5, gradient, err)
      call check(error, err%has_error(), &
                 "the gradient must refuse a core that freezes everything")
   end subroutine test_frozen_gradient_refusal

   subroutine water_ri_gradient(basis, aux_basis, frozen, gradient, err, force_direct)
      !! One converged SCF at `WATER_B`, then the frozen-core RI-MP2 gradient
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: frozen
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(in), optional :: force_direct

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf

      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], WATER_B, basis, mol, err)
      if (err%has_error()) return
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], WATER_B, aux_basis, &
                                  aux, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err)
      if (err%has_error()) return
      call libcint_ri_mp2_gradient(mol, aux, scf%orbitals, scf%orbital_energies, 5, &
                                   gradient, err, n_frozen=frozen, &
                                   force_direct=force_direct)
      call mol%destroy()
      call aux%destroy()
   end subroutine water_ri_gradient

   subroutine test_ri_frozen_gradient_reference(error)
      !! Water/6-31G, oxygen core frozen, against the numpy RI-MP2 reference
      !!
      !! PySCF has no fitted MP2 nuclear gradient -- `mp.dfmp2.DFMP2`'s
      !! `nuc_grad_method` raises NotImplementedError -- so the machine-precision
      !! reference is `tools/cpu_validation/ri_mp2_gradient.py` with the same
      !! frozen count, fed this repository's basis JSON and solving its
      !! Z-vector densely. Against that the two codes agree to 7e-11, the same
      !! floor the frozen=0 control shows, so the bound leaves room for the two
      !! SCFs landing differently, not for any block being wrong. PySCF is
      !! still in the chain where it can be: `pyscf.mp.dfmp2` reproduces the
      !! frozen fitted correlation energy to 7e-15 from the same orbitals, and
      !! central differences of that energy meet this gradient at O(h^2).
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: REFERENCE(3, 3) = reshape( &
                             [0.0_dp, 0.0_dp, 9.512385746353e-03_dp, &
                              0.0_dp, 2.247565931183e-02_dp, -4.756192873177e-03_dp, &
                              0.0_dp, -2.247565931183e-02_dp, -4.756192873176e-03_dp], &
                             [3, 3])
      real(dp), allocatable :: gradient(:, :)
      type(error_t) :: err

      call water_ri_gradient("6-31g", "cc-pvtz-rifit", 1, gradient, err)
      call check(error,.not. err%has_error(), &
                 "the frozen-core RI gradient must run: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(gradient - REFERENCE)) < 2.0e-9_dp, &
                 "the frozen-core RI-MP2 gradient must match the numpy reference")
   end subroutine test_ri_frozen_gradient_reference

   subroutine test_ri_frozen_gradient_fd(error)
      !! The frozen-core RI gradient differences the frozen-core RI energy
      !!
      !! Independent of every external reference: both numbers come from this
      !! code, with the *same* frozen count and the same auxiliary basis on
      !! both sides. The tolerance is the central difference's own O(h^2)
      !! truncation, not the gradient's error.
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: STEP = 2.0e-3_dp
      real(dp), allocatable :: gradient(:, :)
      real(dp) :: moved(3, 3), plus, minus, worst
      type(error_t) :: err
      integer :: iatom, comp

      call water_ri_gradient("sto-3g", "cc-pvdz-rifit", 1, gradient, err)
      call check(error,.not. err%has_error(), &
                 "the frozen-core RI gradient must run: "//err%get_full_trace())
      if (allocated(error)) return

      worst = 0.0_dp
      do iatom = 1, 3
         do comp = 1, 3
            moved = WATER_B
            moved(comp, iatom) = WATER_B(comp, iatom) + STEP
            call ri_frozen_energy_at(moved, plus, err)
            if (err%has_error()) exit
            moved(comp, iatom) = WATER_B(comp, iatom) - STEP
            call ri_frozen_energy_at(moved, minus, err)
            if (err%has_error()) exit
            worst = max(worst, abs((plus - minus)/(2.0_dp*STEP) - gradient(comp, iatom)))
         end do
      end do
      call check(error,.not. err%has_error(), &
                 "the displaced energies must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, worst < 5.0e-6_dp, &
                 "the frozen-core RI gradient must difference the frozen-core RI energy")
   end subroutine test_ri_frozen_gradient_fd

   subroutine ri_frozen_energy_at(coords, energy, err)
      !! One frozen-core RI-MP2 total energy at a displaced geometry
      real(dp), intent(in) :: coords(3, 3)
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2

      energy = 0.0_dp
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], coords, "sto-3g", mol, err)
      if (err%has_error()) return
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], coords, &
                                  "cc-pvdz-rifit", aux, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 300, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err)
      if (err%has_error()) return
      call run_libcint_ri_mp2(mol, aux, scf%orbitals, scf%orbital_energies, 5, &
                              scf%energy, mp2, err, n_frozen=1)
      if (err%has_error()) return
      energy = mp2%total
      call mol%destroy()
      call aux%destroy()
   end subroutine ri_frozen_energy_at

   subroutine test_ri_frozen_gradient_paths(error)
      !! Stored and recomputed reference integrals agree with a core frozen
      !!
      !! The two paths differ only in where the reference's integrals come
      !! from, and neither touches the fitted correlation -- so under a frozen
      !! core they still have to agree to rounding, which is what shows the
      !! new occupied-frozen block feeds both the same numbers.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: stored(:, :), direct(:, :)
      type(error_t) :: err

      call water_ri_gradient("6-31g", "cc-pvdz-rifit", 1, stored, err)
      call check(error,.not. err%has_error(), &
                 "the stored path must run: "//err%get_full_trace())
      if (allocated(error)) return
      call water_ri_gradient("6-31g", "cc-pvdz-rifit", 1, direct, err, &
                             force_direct=.true.)
      call check(error,.not. err%has_error(), &
                 "the direct path must run: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, maxval(abs(stored - direct)) < 1.0e-10_dp, &
                 "recomputing the reference integrals must not move a frozen-core "// &
                 "RI gradient")
   end subroutine test_ri_frozen_gradient_paths

   subroutine test_ri_frozen_gradient_refusal(error)
      !! Freezing every occupied orbital leaves nothing to differentiate
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gradient(:, :)
      type(error_t) :: err

      call water_ri_gradient("sto-3g", "cc-pvdz-rifit", 5, gradient, err)
      call check(error, err%has_error(), &
                 "the RI gradient must refuse a core that freezes everything")
   end subroutine test_ri_frozen_gradient_refusal

end module test_mqc_libcint_mp2

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_mp2, only: collect_mqc_libcint_mp2_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_mp2", collect_mqc_libcint_mp2_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
