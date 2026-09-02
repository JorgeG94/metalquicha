module test_mqc_libcint_cc
   !! Pins coupled cluster against what a reference energy cannot see.
   !!
   !! `validation/check_cc.f90` already compares CCSD and (T) to PySCF, which is
   !! the check that matters and the one that would catch a wrong equation. It
   !! needs the basis files and takes seconds, so it is a manual program rather
   !! than a unit test. What is here instead are the properties that hold by
   !! construction, and whose failure would mean something structurally wrong
   !! rather than numerically off:
   !!
   !!   * `<pq||rs>` is antisymmetric under swapping either pair, and symmetric
   !!     under swapping both. The conversion from chemists' to physicists'
   !!     notation is the likeliest place in the module for a sign error, and
   !!     these identities are exact -- they cannot be satisfied approximately by
   !!     a plausible-looking tensor.
   !!   * a spin-forbidden integral vanishes. Spin integration leaves Kronecker
   !!     deltas, and dropping them gives an answer that still converges.
   !!   * MP2 out of the spin-orbital machinery equals conventional MP2 exactly.
   !!     Two entirely separate routes to one number.
   !!   * CCSD correlation lies below MP2, and (T) below that, for a closed-shell
   !!     molecule at its minimum.
   !!   * freezing a core orbital raises the correlation energy by millihartree,
   !!     not by tens -- which is what a frozen-core indexing error looks like.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf, run_libcint_uhf
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_mp2
   use mqc_libcint_cc, only: cc_result_t, run_libcint_ccsd, spin_orbital_integrals, &
                             so_map_t, build_so_map, so_source_t
   implicit none
   private
   public :: collect_mqc_libcint_cc_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_libcint_cc_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("antisymmetrised_integrals_are_antisymmetric", test_antisymmetry), &
                  new_unittest("spin_forbidden_integrals_vanish", test_spin_blocks), &
                  new_unittest("spin_orbital_mp2_equals_conventional", test_mp2_identity), &
                  new_unittest("ccsd_lies_below_mp2_and_triples_below_that", test_ordering), &
                  new_unittest("freezing_a_core_orbital_raises_correlation", test_frozen), &
                  new_unittest("uccsd_on_a_closed_shell_equals_ccsd", test_uccsd_closed_shell), &
                  new_unittest("uccsd_open_shell_matches_pyscf", test_uccsd_open_shell), &
                  new_unittest("density_fitting_is_refused_when_unrestricted", &
                               test_uccsd_refuses_fitting), &
                  new_unittest("frozen_core_beyond_the_beta_space_is_refused", &
                               test_uccsd_frozen_beta), &
                  new_unittest("an_empty_active_beta_space_still_runs", &
                               test_uccsd_empty_beta) &
                  ]
   end subroutine collect_mqc_libcint_cc_tests

   subroutine water(mol, err, basis)
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
   end subroutine water

   subroutine converged_water(mol, scf, err, basis)
      type(libcint_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis

      call water(mol, err, basis)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                           in_core=.true.)
   end subroutine converged_water

   subroutine test_antisymmetry(error)
      !! <pq||rs> = -<qp||rs> = -<pq||sr> = <qp||sr>
      !!
      !! Checked on every element of a small case rather than on a sample. The
      !! tensor is 14^4 in water/STO-3G, which is cheap enough to check whole, and
      !! a sampled check would pass on a tensor that is antisymmetric almost
      !! everywhere -- which is exactly what one wrong spin delta produces.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: eri(:, :), mo(:, :, :, :), w(:, :, :, :), c_act(:, :)
      integer :: n_so, p, q, r, s
      type(so_map_t) :: map
      type(so_source_t) :: src
      real(dp) :: worst_bra, worst_ket, worst_both

      call converged_water(mol, scf, err, "sto-3g")
      call check(error,.not. err%has_error() .and. scf%converged, "water SCF must converge")
      if (allocated(error)) return

      call build_mo_tensor(mol, scf%orbitals, mo)
      n_so = 2*size(scf%orbitals, 2)
      ! Restricted, so this is the interleaving the module used to compute.
      call build_so_map(size(scf%orbitals, 2), size(scf%orbitals, 2), 0, 0, map)
      ! Restricted: one tensor, read for every spin case.
      src%restricted = .true.
      call move_alloc(mo, src%aa)
      call spin_orbital_integrals(src, map, n_so, w)
      call move_alloc(src%aa, mo)

      worst_bra = 0.0_dp
      worst_ket = 0.0_dp
      worst_both = 0.0_dp
      do s = 1, n_so
         do r = 1, n_so
            do q = 1, n_so
               do p = 1, n_so
                  worst_bra = max(worst_bra, abs(w(p, q, r, s) + w(q, p, r, s)))
                  worst_ket = max(worst_ket, abs(w(p, q, r, s) + w(p, q, s, r)))
                  worst_both = max(worst_both, abs(w(p, q, r, s) - w(q, p, s, r)))
               end do
            end do
         end do
      end do

      call check(error, worst_bra < 1.0e-13_dp, "swapping the bra pair must flip the sign")
      if (allocated(error)) return
      call check(error, worst_ket < 1.0e-13_dp, "swapping the ket pair must flip the sign")
      if (allocated(error)) return
      call check(error, worst_both < 1.0e-13_dp, "swapping both pairs must change nothing")

      call mol%destroy()
   end subroutine test_antisymmetry

   subroutine test_spin_blocks(error)
      !! An integral coupling opposite spins on both pairings is exactly zero
      !!
      !! With interleaved spins, `p` and `p+1` are the same spatial orbital with
      !! opposite spin. <pq||rs> then vanishes identically whenever neither
      !! pairing survives spin integration, and "identically" means bitwise --
      !! this is a term that was never added, not one that cancelled.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: mo(:, :, :, :), w(:, :, :, :)
      integer :: n_so, p, q, r, s
      type(so_map_t) :: map
      type(so_source_t) :: src
      logical :: found
      real(dp) :: worst

      call converged_water(mol, scf, err, "sto-3g")
      call check(error,.not. err%has_error() .and. scf%converged, "water SCF must converge")
      if (allocated(error)) return

      call build_mo_tensor(mol, scf%orbitals, mo)
      n_so = 2*size(scf%orbitals, 2)
      ! Restricted, so this is the interleaving the module used to compute.
      call build_so_map(size(scf%orbitals, 2), size(scf%orbitals, 2), 0, 0, map)
      ! Restricted: one tensor, read for every spin case.
      src%restricted = .true.
      call move_alloc(mo, src%aa)
      call spin_orbital_integrals(src, map, n_so, w)
      call move_alloc(src%aa, mo)

      worst = 0.0_dp
      found = .false.
      do s = 1, n_so
         do r = 1, n_so
            do q = 1, n_so
               do p = 1, n_so
                  ! Neither (p,r)(q,s) nor (p,s)(q,r) is spin-allowed.
                  if (mod(p, 2) /= mod(r, 2) .and. mod(p, 2) /= mod(s, 2)) then
                     found = .true.
                     worst = max(worst, abs(w(p, q, r, s)))
                  end if
               end do
            end do
         end do
      end do

      call check(error, found, "there must be spin-forbidden combinations to check")
      if (allocated(error)) return
      call check(error, worst == 0.0_dp, "a spin-forbidden integral must be exactly zero")

      call mol%destroy()
   end subroutine test_spin_blocks

   subroutine test_mp2_identity(error)
      !! Spin-orbital MP2 must reproduce the conventional path exactly
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2
      type(cc_result_t) :: cc
      type(error_t) :: err

      call converged_water(mol, scf, err, "sto-3g")
      call check(error,.not. err%has_error() .and. scf%converged, "water SCF must converge")
      if (allocated(error)) return

      call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, 5, scf%energy, &
                           mp2, err)
      call check(error,.not. err%has_error(), "conventional MP2 must run: "//err%get_full_trace())
      if (allocated(error)) return

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, 0, &
                            60, 1.0e-11_dp, .false., .false., cc, err)
      call check(error,.not. err%has_error() .and. cc%converged, "CCSD must converge")
      if (allocated(error)) return

      ! Two separate transforms, two separate energy expressions, one number.
      call check(error, abs(cc%e_mp2 - mp2%correlation) < 1.0e-12_dp, &
                 "spin-orbital MP2 must equal the conventional MP2 it is built beside")

      call mol%destroy()
   end subroutine test_mp2_identity

   subroutine test_ordering(error)
      !! E_CCSD < E_MP2 < 0, and (T) below CCSD, at a closed-shell minimum
      !!
      !! Not a law in general -- MP2 can overshoot, and (T) can come out positive
      !! for a stretched bond or a poor reference -- which is why this is asserted
      !! on water at its equilibrium geometry, where it is reliably true. It
      !! catches a sign error in the energy expression that a converged
      !! calculation would otherwise report happily.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(cc_result_t) :: cc
      type(error_t) :: err
      real(dp) :: e_ccsd

      call converged_water(mol, scf, err, "sto-3g")
      call check(error,.not. err%has_error() .and. scf%converged, "water SCF must converge")
      if (allocated(error)) return

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, 0, &
                            60, 1.0e-11_dp, .true., .false., cc, err)
      call check(error,.not. err%has_error() .and. cc%converged, "CCSD must converge")
      if (allocated(error)) return

      e_ccsd = cc%e_singles + cc%e_doubles
      call check(error, cc%e_mp2 < 0.0_dp, "MP2 correlation must be negative")
      if (allocated(error)) return
      call check(error, e_ccsd < cc%e_mp2, "CCSD must recover more correlation than MP2")
      if (allocated(error)) return
      call check(error, cc%e_triples < 0.0_dp, "(T) must be negative here")
      if (allocated(error)) return
      ! A guard on magnitude as well as sign: (T) is a correction, and one the
      ! size of the doubles would mean the 1/36 or a permutation sign is wrong.
      call check(error, abs(cc%e_triples) < 0.1_dp*abs(e_ccsd), &
                 "(T) must be small beside the CCSD correlation energy")

      call mol%destroy()
   end subroutine test_ordering

   subroutine test_frozen(error)
      !! Freezing the oxygen core raises the correlation energy, by millihartree
      !!
      !! The magnitude is the content. A frozen-core off-by-one that dropped a
      !! valence orbital instead of the core, or shifted the orbital energies
      !! against the coefficients, changes this by tens of millihartree -- and
      !! still converges.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(cc_result_t) :: all_e, frozen_e
      type(error_t) :: err
      real(dp) :: e_all, e_frozen, gap

      call converged_water(mol, scf, err, "sto-3g")
      call check(error,.not. err%has_error() .and. scf%converged, "water SCF must converge")
      if (allocated(error)) return

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, 0, &
                            60, 1.0e-11_dp, .false., .false., all_e, err)
      call check(error,.not. err%has_error() .and. all_e%converged, "CCSD must converge")
      if (allocated(error)) return
      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, 1, &
                            60, 1.0e-11_dp, .false., .false., frozen_e, err)
      call check(error,.not. err%has_error() .and. frozen_e%converged, &
                 "frozen-core CCSD must converge")
      if (allocated(error)) return

      e_all = all_e%e_singles + all_e%e_doubles
      e_frozen = frozen_e%e_singles + frozen_e%e_doubles
      gap = e_frozen - e_all

      call check(error, gap > 0.0_dp, "freezing a core orbital must raise the energy")
      if (allocated(error)) return
      call check(error, gap < 0.05_dp, &
                 "the frozen-core shift must be millihartree, not tens of them")

      call mol%destroy()
   end subroutine test_frozen

   subroutine build_mo_tensor(mol, coeff, mo)
      !! The active MO tensor, the way the CC driver builds it
      use mqc_libcint_mp2, only: transform_block
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)
      real(dp), allocatable, intent(out) :: mo(:, :, :, :)

      real(dp), allocatable :: eri(:, :)

      call mol%eris_packed(eri)
      call transform_block(eri, coeff, coeff, coeff, coeff, mo)
      deallocate (eri)
   end subroutine build_mo_tensor

   subroutine test_uccsd_closed_shell(error)
      !! An unrestricted reference on a closed shell must give the restricted answer
      !!
      !! This is the identity the whole unrestricted path rests on. UHF on a
      !! closed-shell singlet converges to the RHF solution, so the two spins end
      !! up with the same orbitals and the same energies -- but they arrive at the
      !! amplitude equations through different machinery: three spatial tensors
      !! instead of one, and orbital energies drawn per spin. Any wrong spin case
      !! in `chem` and any transposition in the mixed `(aa|bb)` tensor land here.
      !!
      !! The *ordering* is not what this case exercises, and saying otherwise
      !! would be a comfortable mistake. `build_so_map` branches on the counts
      !! rather than on the reference, and a closed-shell UHF has the counts
      !! matching, so this takes the interleaving -- the same one the restricted
      !! path takes. The blocked ordering is covered by the open-shell case
      !! below, which is the only one where the counts differ.
      !!
      !! Not asked bit for bit. UHF and RHF are separate SCFs and converge to the
      !! same solution only to their own tolerance, so the reference orbitals
      !! differ in the last few digits before coupled cluster ever runs. The
      !! correlation energy is what is compared, and 1e-9 is well inside the gap
      !! between "the same method" and "a spin case is wrong", which shows up in
      !! the third decimal or worse.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: rscf, uscf
      type(error_t) :: err
      type(cc_result_t) :: rcc, ucc

      call water(mol, err, "sto-3g")
      call check(error,.not. err%has_error(), "water must build")
      if (allocated(error)) return

      call run_libcint_rhf(mol, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., rscf, err, &
                           in_core=.true.)
      call check(error,.not. err%has_error() .and. rscf%converged, "RHF must converge")
      if (allocated(error)) return

      call run_libcint_uhf(mol, 10, 1, 200, 1.0e-11_dp, 1.0e-9_dp, .false., uscf, err, &
                           in_core=.true.)
      call check(error,.not. err%has_error() .and. uscf%converged, "UHF must converge")
      if (allocated(error)) return

      call run_libcint_ccsd(mol, rscf%orbitals, rscf%orbital_energies, 5, 0, &
                            60, 1.0e-10_dp, .false., .false., rcc, err)
      call check(error,.not. err%has_error() .and. rcc%converged, "CCSD must converge")
      if (allocated(error)) return

      call run_libcint_ccsd(mol, uscf%orbitals, uscf%orbital_energies, 5, 0, &
                            60, 1.0e-10_dp, .false., .false., ucc, err, &
                            coeff_b=uscf%orbitals_beta, &
                            orbital_energies_b=uscf%orbital_energies_beta, n_occ_b=5)
      call mol%destroy()
      call check(error,.not. err%has_error() .and. ucc%converged, "UCCSD must converge")
      if (allocated(error)) return

      call check(error, abs(ucc%e_mp2 - rcc%e_mp2) < 1.0e-9_dp, &
                 "the MP2 checkpoint disagrees between the restricted and "// &
                 "unrestricted paths, so the three-tensor transform is wrong "// &
                 "before coupled cluster ever sees it")
      if (allocated(error)) return

      call check(error, abs(ucc%e_correlation - rcc%e_correlation) < 1.0e-9_dp, &
                 "UCCSD on a closed shell does not reproduce CCSD")
   end subroutine test_uccsd_closed_shell

   subroutine test_uccsd_open_shell(error)
      !! A doublet against PySCF's GCCSD
      !!
      !! The closed-shell identity above cannot see the thing this is for. With
      !! `n_alpha == n_beta` the blocked ordering and the interleaved one place
      !! occupied in the same span, so a virtual index that lands past the wrong
      !! spin's occupied count still works out. A doublet is the smallest case
      !! where they differ, and where the mixed `(aa|bb)` tensor carries a
      !! genuinely different number from `(aa|aa)`.
      !!
      !! Checked against **GCCSD** rather than UCCSD. PySCF's GCCSD is general
      !! spin-orbital coupled cluster, which is what this module is; UCCSD is
      !! spin-blocked. They agree with each other to 1.8e-16, so the choice costs
      !! nothing and buys a comparison between two implementations of the same
      !! formulation rather than of two different ones.
      !!
      !! OH at 0.9697 A in STO-3G: five alpha and four beta electrons in six
      !! spatial orbitals.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(cc_result_t) :: ucc
      real(dp) :: c(3, 2)
      integer, parameter :: Z(2) = [8, 1]
      character(len=2), parameter :: SYM(2) = ["O ", "H "]
      !! PySCF 2.14 GCCSD from a UHF reference, this repository's own STO-3G.
      real(dp), parameter :: PYSCF_UHF = -74.362637545612_dp
      real(dp), parameter :: PYSCF_CORR = -0.024496069957_dp
      !! And its (T), which is spin-orbital and so should need no new equations.
      real(dp), parameter :: PYSCF_T = -0.000000231291_dp

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9697_dp*ANG], [3, 2])
      call build_libcint_molecule(Z, SYM, c, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "OH must build")
      if (allocated(error)) return

      call run_libcint_uhf(mol, 9, 2, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                           in_core=.true.)
      call check(error,.not. err%has_error() .and. scf%converged, "UHF must converge")
      if (allocated(error)) return
      call check(error, abs(scf%energy - PYSCF_UHF) < 1.0e-9_dp, &
                 "the UHF reference disagrees with PySCF, so nothing downstream "// &
                 "is being tested")
      if (allocated(error)) return

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, 0, &
                            60, 1.0e-10_dp, .true., .false., ucc, err, &
                            coeff_b=scf%orbitals_beta, &
                            orbital_energies_b=scf%orbital_energies_beta, n_occ_b=4)
      call mol%destroy()
      call check(error,.not. err%has_error() .and. ucc%converged, "UCCSD must converge")
      if (allocated(error)) return

      call check(error, abs(ucc%e_correlation - ucc%e_triples - PYSCF_CORR) < 1.0e-8_dp, &
                 "UCCSD on a doublet disagrees with PySCF GCCSD")
      if (allocated(error)) return

      ! (T) is tiny here -- 2.3e-7 -- so an absolute tolerance would pass on a
      ! zero. Compared to a tenth of its own size instead.
      call check(error, abs(ucc%e_triples - PYSCF_T) < 0.1_dp*abs(PYSCF_T), &
                 "UCCSD(T) on a doublet disagrees with PySCF G(T)")
   end subroutine test_uccsd_open_shell

   subroutine test_uccsd_refuses_fitting(error)
      !! RI over an unrestricted reference is refused, not silently run wrong
      !!
      !! The fitted blocks and the fitted ladder both still assume alpha and beta
      !! share a spatial orbital. Handed unrestricted orbitals they would not
      !! fail: they would converge, to the restricted answer built from the alpha
      !! orbitals alone. That is the failure mode worth a test -- an error is
      !! recoverable and a plausible wrong number is not.
      !!
      !! Delete this when URI-CCSD exists.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(cc_result_t) :: ucc
      real(dp) :: c(3, 2)
      integer, parameter :: Z(2) = [8, 1]
      character(len=2), parameter :: SYM(2) = ["O ", "H "]

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9697_dp*ANG], [3, 2])
      call build_libcint_molecule(Z, SYM, c, "sto-3g", mol, err)
      call build_libcint_molecule(Z, SYM, c, "cc-pvdz-rifit", aux, err)
      call check(error,.not. err%has_error(), "OH and its auxiliary basis must build")
      if (allocated(error)) return

      call run_libcint_uhf(mol, 9, 2, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                           in_core=.true.)
      call check(error,.not. err%has_error() .and. scf%converged, "UHF must converge")
      if (allocated(error)) return

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, 0, &
                            60, 1.0e-10_dp, .false., .false., ucc, err, &
                            coeff_b=scf%orbitals_beta, &
                            orbital_energies_b=scf%orbital_energies_beta, n_occ_b=4, &
                            aux=aux)
      call mol%destroy()
      call aux%destroy()

      call check(error, err%has_error(), &
                 "density fitting over an unrestricted reference was accepted; it "// &
                 "would have returned the restricted answer for the alpha orbitals")
   end subroutine test_uccsd_refuses_fitting

   subroutine test_uccsd_frozen_beta(error)
      !! A frozen core deeper than the beta occupation is refused
      !!
      !! The alpha guard alone does not catch it. OH here has five alpha and
      !! four beta electrons, so a frozen count of five clears `frozen < n_occ`
      !! on the alpha side and leaves the beta active space at `-1`.
      !!
      !! What that costs is not a wrong number. `build_so_map` sizes its arrays
      !! `n_o + n_ob + n_v + n_vb`, which a negative `n_ob` shrinks, and then
      !! fills `n_o + n_v + n_vb` of them -- the beta occupied loop does not run
      !! at a negative count. It writes past the end of its own allocation, so
      !! without this guard the symptom is whatever that memory happened to be.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(cc_result_t) :: ucc
      real(dp) :: c(3, 2)
      integer, parameter :: Z(2) = [8, 1]
      character(len=2), parameter :: SYM(2) = ["O ", "H "]

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9697_dp*ANG], [3, 2])
      call build_libcint_molecule(Z, SYM, c, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "OH must build")
      if (allocated(error)) return

      call run_libcint_uhf(mol, 9, 2, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                           in_core=.true.)
      call check(error,.not. err%has_error() .and. scf%converged, "UHF must converge")
      if (allocated(error)) return

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, 5, &
                            60, 1.0e-10_dp, .false., .false., ucc, err, &
                            coeff_b=scf%orbitals_beta, &
                            orbital_energies_b=scf%orbital_energies_beta, n_occ_b=4)
      call mol%destroy()

      call check(error, err%has_error(), &
                 "a frozen core deeper than the beta occupation was accepted; the "// &
                 "spin-orbital map would have been written past its end")
   end subroutine test_uccsd_frozen_beta

   subroutine test_uccsd_empty_beta(error)
      !! No *active* beta orbitals is a calculation, not an error
      !!
      !! The boundary the obvious guard gets wrong. Refusing at
      !! `frozen >= n_occ_b` would reject frozen-core lithium, whose one valence
      !! electron is alpha and whose beta electron is entirely core -- a real
      !! calculation with an empty active beta space. Only `frozen > n_occ_b` is
      !! out of bounds.
      !!
      !! Lithium at sto-3g: two alpha, one beta, one frozen core orbital. So the
      !! active space is one alpha occupied and none beta, and what is asserted
      !! is that it runs and produces a correlation energy at all -- with a
      !! single active electron there are no double excitations, so the number
      !! itself is small and the point is the shape.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(cc_result_t) :: ucc
      real(dp) :: c(3, 1)
      integer, parameter :: Z(1) = [3]
      character(len=2), parameter :: SYM(1) = ["Li"]

      c = 0.0_dp
      call build_libcint_molecule(Z, SYM, c, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "lithium must build")
      if (allocated(error)) return

      call run_libcint_uhf(mol, 3, 2, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                           in_core=.true.)
      call check(error,.not. err%has_error() .and. scf%converged, "UHF must converge")
      if (allocated(error)) return

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 2, 1, &
                            60, 1.0e-10_dp, .false., .false., ucc, err, &
                            coeff_b=scf%orbitals_beta, &
                            orbital_energies_b=scf%orbital_energies_beta, n_occ_b=1)
      call mol%destroy()

      call check(error,.not. err%has_error(), &
                 "an empty active beta space was refused, but frozen-core lithium "// &
                 "is a calculation: "//err%get_message())
      if (allocated(error)) return
      call check(error, ucc%converged, "it should converge")
      if (allocated(error)) return
      call check(error, ucc%e_correlation <= 0.0_dp, &
                 "a correlation energy should come back, and not be positive")
   end subroutine test_uccsd_empty_beta

end module test_mqc_libcint_cc

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_cc, only: collect_mqc_libcint_cc_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_cc", collect_mqc_libcint_cc_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
