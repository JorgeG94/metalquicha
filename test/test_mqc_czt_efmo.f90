!! The EFMO orchestrator, at the four limits where its total is known
module test_mqc_czt_efmo
   !! `run_efmo` assembles eq 6 of the EFMO paper (Sattasathuchana et al., JCTC
   !! 20, 2445 (2024)):
   !!
   !!     E = sum_I E_I^0
   !!       + sum_{R_IJ <= R_cut} (E_IJ^0 - E_I^0 - E_J^0 - E_IJ^pol)
   !!       + sum_{R_IJ >  R_cut} (E_IJ^Coul + E_IJ^disp + E_IJ^ExRep + E_IJ^CT)
   !!       + E_pol^total
   !!
   !! No reference for an EFMO energy exists yet -- Phase 3 gets those from
   !! GAMESS -- so what is checked here is not a number but the **identities the
   !! expression has by construction**, each of which fails if a term is
   !! dropped, double counted, or assembled with the wrong sign:
   !!
   !!   1. `R_cut` huge: every pair is quantum, so the expression collapses to
   !!      FMO2 in vacuo plus the *many-body* part of the induction.
   !!   2. `R_cut` zero: no pair is quantum, so it collapses to the in-vacuo
   !!      monomers plus the whole EFP-EFP interaction energy.
   !!   3. Two fragments, both quantum: their pair induction *is* the total
   !!      induction, the two cancel exactly, and every other term telescopes
   !!      away -- so the answer is the dimer's own RHF energy.
   !!   4. `R_cut` between the pair separations of a trimer: one quantum dimer
   !!      and two effective ones, against a total assembled here by hand from
   !!      the Phase 1 pieces.
   !!   5. Level two through the general many-body machinery is the pair form
   !!      **to the last bit**, not to a tolerance: the difference operator
   !!      applied to a pair has to be the subtraction that was written out
   !!      before it existed, or every reference in the tree moves.
   !!   6. Level equal to the fragment count with `R_cut` huge is the
   !!      *unfragmented* energy. The in-vacuo series telescopes to the
   !!      supersystem's own SCF and the induction series telescopes to
   !!      `E_pol^total`, which then cancels the last term of the expression
   !!      exactly -- so the whole polarization correction disappears and what
   !!      is left is one RHF energy. It is the sharpest check there is on the
   !!      subset-induction bookkeeping: a wrong sign, a missing subset or a
   !!      group induction solved on the wrong fragments all survive the
   !!      level-two limits and fail here.
   !!   7. The induction series alone: `sum_S dE_S^pol` over every subset of a
   !!      trimer is `E_pol^total`, asserted apart from the energy so that a
   !!      failure in six names which half moved.
   !!
   !! **Limit 1 is not "EFMO is nearly FMO2".** The induction is strongly
   !! non-additive -- three waters at four Angstrom carry 44 per cent of their
   !! total induction in terms no pair has -- so the difference in test one is a
   !! term of the method, of the same size as the pair energies. Only the
   !! identity is asserted, and the measured difference is reported by the run
   !! itself.
   !!
   !! Water in 6-31G, one trimer geometry for every case, and the reference
   !! potentials built once: MAKEFP is the whole cost of this file, and
   !! `run_efmo` builds its own set on each call.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_efmo, only: efmo_options_t, efmo_result_t, run_efmo, &
                           EFMO_CORR_NONE, EFMO_CORR_RI_MP2
   use mqc_czt_mp2, only: mp2_result_t, run_czt_ri_mp2
   use mqc_elements, only: core_orbital_count
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   use mqc_czt_efp_potential, only: efp_potential_t, make_efp_potential
   use mqc_czt_efp_read, only: efp_fragment_t
   use mqc_czt_efp_convert, only: potential_to_fragment
   use mqc_czt_efp_energy, only: efp_energy_t, efp_interaction_energy, &
                                 efp_pair_energy_t, efp_pair_terms, &
                                 pair_polarization_energy, subset_polarization_energy
   use mqc_czt_efp_interaction, only: efp_system_t, build_efp_system, polarization_energy
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_atomic_guess, only: build_restricted_guess
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_scf_types, only: scf_numerics_t
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_efmo_tests

   real(dp), parameter :: ANG = ANGSTROM_TO_BOHR

   character(len=*), parameter :: BASIS = "6-31g"
      !! s and p only, so the Cartesian form `make_efp_potential` forces and the
      !! spherical form `run_fmo2` builds are the same basis. Above d they are
      !! not, and test one would then compare two different models.

   character(len=*), parameter :: AUX = "cc-pvdz-rifit"
      !! The fitting set the correlated cases use. It does not match the orbital
      !! basis, and does not have to: what the two correlated tests assert is
      !! that `run_efmo` correlates the same orbitals with the same fitting
      !! space that the reference here does, which is true of any auxiliary
      !! basis. A matched set would make the number closer to conventional MP2
      !! and the test no sharper.

   real(dp), parameter :: TOL = 1.0e-9_dp
      !! Every identity below is exact, so what this has to clear is the SCF
      !! convergence the pieces are held to -- 1e-10 on the energy -- summed
      !! over the handful of SCFs a trimer needs.

   real(dp), parameter :: SPACING_12 = 4.0_dp
   real(dp), parameter :: SPACING_13 = 12.0_dp
      !! Where waters two and three sit along `x`, in Angstrom, from water one.
      !! Oxygen's Bondi radius is 1.52, so `R_IJ` is the separation over 3.04:
      !! 1.32 for the first pair, 3.95 and 2.63 for the other two. `R_cut = 2.0`
      !! therefore splits them one to two, which is what test four needs, and
      !! the two limits sweep past both ends of that range.

   ! Built once. Every test wants the same three potentials, and they are
   ! identical every time.
   type(efp_fragment_t), save :: cached_frag(3)
   real(dp), save :: cached_mono(3) = 0.0_dp
   real(dp), save :: cached_dimer_12 = 0.0_dp
   real(dp), save :: cached_mono_mp2(3) = 0.0_dp
      !! Each monomer's RI-MP2 correlation energy, on the same orbitals its
      !! potential was built from -- which is where `run_efmo` takes its own.
   logical, save :: cached_ready = .false.

contains

   subroutine collect_mqc_czt_efmo_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("efmo_all_pairs_quantum_is_fmo2_plus_many_body_induction", &
                               test_all_quantum), &
                  new_unittest("efmo_no_pair_quantum_is_monomers_plus_efp", test_no_quantum), &
                  new_unittest("efmo_two_fragments_is_the_dimer_energy", test_two_fragments), &
                  new_unittest("efmo_trimer_split_one_quantum_two_effective", test_mixed), &
                  new_unittest("efmo_rimp2_two_fragments_is_the_dimer_rimp2_energy", &
                               test_rimp2_dimer), &
                  new_unittest("efmo_rimp2_all_quantum_is_the_correlated_pair_sum", &
                               test_rimp2_trimer), &
                  new_unittest("efmo_level_two_reduces_to_the_pair_form_exactly", &
                               test_level_two_is_the_pair_form), &
                  new_unittest("efmo_full_level_is_the_unfragmented_energy", &
                               test_full_level), &
                  new_unittest("efmo_induction_series_sums_to_the_total", &
                               test_induction_series) &
                  ]
   end subroutine collect_mqc_czt_efmo_tests

   subroutine water_geometry(z, symbols, coords)
      !! One water, in Bohr, in the `yz` plane
      integer, intent(out) :: z(3)
      character(len=2), intent(out) :: symbols(3)
      real(dp), intent(out) :: coords(3, 3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      coords = reshape([0.00000000000000_dp, 0.00000000009155_dp, 0.10077199490609_dp, &
                        0.00000000000000_dp, 0.77250895271063_dp, -0.46780199741728_dp, &
                        0.00000000000000_dp, -0.77250895280218_dp, -0.46780199748881_dp], &
                       [3, 3])*ANG
   end subroutine water_geometry

   subroutine water_chain(n, z, symbols, xyz, owner)
      !! `n` copies of that water along `x`, at the spacings above
      !!
      !! `n = 2` is exactly waters one and two of the trimer, at the same
      !! coordinates, so the two-fragment test and the trimer tests share a
      !! geometry and a dimer energy.
      integer, intent(in) :: n
      integer, intent(out) :: z(:), owner(:)
      character(len=2), intent(out) :: symbols(:)
      real(dp), intent(out) :: xyz(:, :)

      integer :: zw(3), k, f, at
      character(len=2) :: sw(3)
      real(dp) :: coords(3, 3), offset(3)

      call water_geometry(zw, sw, coords)
      at = 0
      do f = 1, n
         offset = 0.0_dp
         if (f == 2) offset(1) = SPACING_12*ANG
         if (f == 3) offset(1) = SPACING_13*ANG
         do k = 1, 3
            at = at + 1
            z(at) = zw(k)
            symbols(at) = sw(k)
            xyz(:, at) = coords(:, k) + offset
            owner(at) = f
         end do
      end do
   end subroutine water_chain

   subroutine efmo_settings(opts)
      !! The options every case runs with
      !!
      !! MAKEFP's own SCF defaults, which `run_efmo` uses for its dimers too:
      !! the reference energies computed here have to come from the same SCF, or
      !! an identity exact in exact arithmetic fails on convergence.
      type(efmo_options_t), intent(out) :: opts

      opts%basis = BASIS
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-10_dp
      opts%scf_density_tol = 1.0e-8_dp
      opts%scf_grad_tol = 1.0e-8_dp
      opts%scf%grad_tol = 1.0e-8_dp
   end subroutine efmo_settings

   subroutine build_reference(err)
      !! The three potentials, their `E_I^0`, and the 1-2 dimer's RHF energy
      type(error_t), intent(inout) :: err

      integer :: z(9), owner(9)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9)
      type(efp_potential_t) :: pot
      type(efmo_options_t) :: opts
      type(rhf_result_t) :: scf
      integer :: k
      integer, allocatable :: idx(:)

      if (cached_ready) return
      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)

      do k = 1, 3
         idx = atoms_of(owner, k)
         call make_efp_potential(z(idx), symbols(idx), xyz(:, idx), BASIS, "FRAG", pot, &
                                 err, charge=0, &
                                 energy_tol=opts%scf_energy_tol, &
                                 density_tol=opts%scf_density_tol, &
                                 grad_tol_in=opts%scf_grad_tol, scf_in=opts%scf, &
                                 max_iter_in=opts%scf_max_iter, scf_out=scf)
         if (err%has_error()) return
         cached_mono(k) = pot%scf_energy
         cached_mono_mp2(k) = ri_mp2_on(z(idx), symbols(idx), xyz(:, idx), &
                                        sum(z(idx)), scf, opts, err)
         if (err%has_error()) return
         call potential_to_fragment(pot, cached_frag(k), err)
         call pot%destroy()
         if (err%has_error()) return
      end do

      call dimer_rhf(z(1:6), symbols(1:6), xyz(:, 1:6), 0, opts, cached_dimer_12, err)
      if (err%has_error()) return
      cached_ready = .true.
   end subroutine build_reference

   pure function atoms_of(owner, k) result(idx)
      !! The indices of fragment `k`'s atoms
      integer, intent(in) :: owner(:)
      integer, intent(in) :: k
      integer, allocatable :: idx(:)

      integer :: i, n

      n = count(owner == k)
      allocate (idx(n))
      n = 0
      do i = 1, size(owner)
         if (owner(i) == k) then
            n = n + 1
            idx(n) = i
         end if
      end do
   end function atoms_of

   subroutine dimer_rhf(z, symbols, xyz, charge, opts, energy, err)
      !! `E_IJ^0` computed here, the way `run_efmo` computes it
      !!
      !! Deliberately a second implementation rather than a call into the
      !! module: what test three and test four assert is that the orchestrator's
      !! dimer is *this* number, and reusing its own routine would assert
      !! nothing. Cartesian, because the monomer SCFs behind the potentials are.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: charge
      type(efmo_options_t), intent(in) :: opts
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: guess_density(:, :)
      integer :: guess_kind

      energy = 0.0_dp
      call build_czt_molecule(z, symbols, xyz, BASIS, mol, err, force_cartesian=.true.)
      if (err%has_error()) return
      call build_restricted_guess(mol, "auto", guess_kind, guess_density, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, sum(z) - charge, opts%scf_max_iter, opts%scf_energy_tol, &
                       opts%scf_density_tol, .false., scf, err, &
                       guess=guess_kind, guess_density=guess_density, &
                       grad_tol=opts%scf_grad_tol, scf=opts%scf)
      call mol%destroy()
      if (err%has_error()) return
      energy = scf%energy
   end subroutine dimer_rhf

   function total_induction(frags, err) result(energy)
      !! `E_pol^total` over the fragments given
      type(efp_fragment_t), intent(in) :: frags(:)
      type(error_t), intent(inout) :: err
      real(dp) :: energy

      type(efp_system_t) :: system
      real(dp), allocatable :: shifts(:, :)

      allocate (shifts(3, size(frags)), source=0.0_dp)
      energy = 0.0_dp
      call build_efp_system(frags, shifts, system, err)
      if (err%has_error()) return
      energy = polarization_energy(system, frags, err)
      call system%destroy()
   end function total_induction

   function pair_induction_sum(frags, err) result(energy)
      !! `sum_{I<J} E_IJ^pol` over every pair of the fragments given
      type(efp_fragment_t), intent(in) :: frags(:)
      type(error_t), intent(inout) :: err
      real(dp) :: energy

      real(dp) :: zero(3)
      integer :: a, b

      zero = 0.0_dp
      energy = 0.0_dp
      do a = 1, size(frags) - 1
         do b = a + 1, size(frags)
            energy = energy + pair_polarization_energy(frags(a), frags(b), zero, zero, err)
            if (err%has_error()) return
         end do
      end do
   end function pair_induction_sum

   subroutine test_all_quantum(error)
      !! Limit one: `R_cut` huge, so EFMO is FMO2 in vacuo plus many-body induction
      !!
      !! With every pair quantum, eq 6 has no far sum and reads
      !!
      !!     sum_I E_I^0 + sum_IJ (E_IJ^0 - E_I^0 - E_J^0)  -  sum_IJ E_IJ^pol
      !!                                                    +  E_pol^total
      !!
      !! whose first line is exactly the many-body expansion at level two over
      !! in-vacuo fragments -- `run_fmo2` with the embedding off and the MBE
      !! assembly, which is a completely independent implementation of that sum.
      !! So the difference between the two totals must be the induction
      !! remainder and nothing else.
      !!
      !! **That remainder is not small.** It is the part of the induction no
      !! pair carries, and for this trimer it is a sizeable fraction of the
      !! total; the plan's earlier expectation that the two energies nearly
      !! agree was wrong. Only the identity is asserted.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      type(fmo_options_t) :: fmo_opts
      type(fmo_result_t) :: fmo_res
      integer :: z(9), owner(9)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9), e_pol_total, e_pol_pairs, expected

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 1.0e6_dp
      call run_efmo(z, symbols, xyz, owner, [0, 0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%n_qm_pairs, 3, message="a huge cutoff left a pair effective")
      if (allocated(error)) return
      call check(error, res%n_efp_pairs, 0, message="a huge cutoff kept an EFP pair")
      if (allocated(error)) return

      ! FMO2 with no embedding and the MBE assembly: sum E_I + sum (E_IJ - E_I - E_J),
      ! every fragment solved in vacuo. The same first line as above, from other code.
      fmo_opts%basis = BASIS
      fmo_opts%esp = "none"
      fmo_opts%expansion = "mbe"
      fmo_opts%level = 2
      fmo_opts%scf_max_iter = 200
      fmo_opts%scf_energy_tol = 1.0e-10_dp
      fmo_opts%scf_density_tol = 1.0e-8_dp
      fmo_opts%scf%grad_tol = 1.0e-8_dp
      call run_fmo2(z, symbols, xyz, owner, fmo_opts, fmo_res, err)
      call check(error,.not. err%has_error(), "run_fmo2 failed: "//err%get_full_trace())
      if (allocated(error)) return

      e_pol_total = total_induction(cached_frag, err)
      e_pol_pairs = pair_induction_sum(cached_frag, err)
      call check(error,.not. err%has_error(), "the induction terms failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      expected = fmo_res%energy + (e_pol_total - e_pol_pairs)
      call check(error, res%energy, expected, thr=TOL, &
                 message="EFMO with every pair quantum is not FMO2 in vacuo plus the "// &
                 "many-body induction")
      if (allocated(error)) return

      ! And that the orchestrator's own two induction sums are the ones just
      ! computed, so a failure above says which half moved.
      call check(error, res%polarization_total, e_pol_total, thr=1.0e-12_dp, &
                 message="E_pol^total")
      if (allocated(error)) return
      call check(error, res%induction_correction, e_pol_pairs, thr=1.0e-12_dp, &
                 message="sum E_IJ^pol")
      if (allocated(error)) return
      ! The remainder is a term of the method, not a residue: assert it is not
      ! negligible, so a version that silently dropped it would fail here.
      call check(error, abs(e_pol_total - e_pol_pairs) > 1.0e-6_dp, &
                 "the many-body induction vanished, so this geometry cannot tell a "// &
                 "dropped remainder from a correct one")
   end subroutine test_all_quantum

   subroutine test_no_quantum(error)
      !! Limit two: `R_cut` zero, so EFMO is in-vacuo monomers plus EFP-EFP
      !!
      !! With no quantum dimer there is no near sum at all, and eq 6 reads
      !! `sum_I E_I^0` plus the four pair terms over every pair plus
      !! `E_pol^total` -- which is exactly `sum_I E_I^0` plus what
      !! `efp_interaction_energy` returns for the same fragments, since that
      !! routine's five terms are those four and that induction.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      type(efp_energy_t) :: efp
      real(dp) :: shifts(3, 3), expected
      integer :: z(9), owner(9)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9)

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 0.0_dp
      call run_efmo(z, symbols, xyz, owner, [0, 0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%n_qm_pairs, 0, message="a zero cutoff left a pair quantum")
      if (allocated(error)) return
      call check(error, res%n_efp_pairs, 3, message="a zero cutoff lost an EFP pair")
      if (allocated(error)) return
      call check(error, res%nmer_correction, 0.0_dp, thr=0.0_dp, &
                 message="there is no quantum dimer, so there is no dimer correction")
      if (allocated(error)) return
      call check(error, res%induction_correction, 0.0_dp, thr=0.0_dp, &
                 message="there is no quantum dimer, so nothing subtracts pair induction")
      if (allocated(error)) return

      shifts = 0.0_dp
      efp = efp_interaction_energy(cached_frag, shifts, err)
      call check(error,.not. err%has_error(), "the EFP energy failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      expected = sum(cached_mono) + efp%total
      call check(error, res%energy, expected, thr=TOL, &
                 message="EFMO with no quantum dimer is not the in-vacuo monomers "// &
                 "plus the EFP-EFP interaction")
      if (allocated(error)) return
      call check(error, res%monomer_sum, sum(cached_mono), thr=TOL, &
                 message="the monomer sum is not the sum of the potentials' own SCFs")
   end subroutine test_no_quantum

   function ri_mp2_on(z, symbols, xyz, nelec, scf, opts, err) result(energy)
      !! The RI-MP2 correlation energy on orbitals already converged
      !!
      !! A second implementation of what `run_efmo` does after each of its
      !! SCFs, written out here for the same reason `dimer_rhf` is: calling the
      !! module's own routine would assert nothing about it. The frozen core is
      !! counted from the elements, which is `run_efmo`'s default, and the
      !! molecule is Cartesian because the orbitals came from a Cartesian one.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: nelec
      type(rhf_result_t), intent(in) :: scf
      type(efmo_options_t), intent(in) :: opts
      type(error_t), intent(inout) :: err
      real(dp) :: energy

      type(czt_molecule_t) :: mol, fit
      type(mp2_result_t) :: mp2

      energy = 0.0_dp
      call build_czt_molecule(z, symbols, xyz, BASIS, mol, err, force_cartesian=.true.)
      if (err%has_error()) return
      call build_czt_molecule(z, symbols, xyz, AUX, fit, err, force_cartesian=.true.)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if
      call run_czt_ri_mp2(mol, fit, scf%orbitals, scf%orbital_energies, nelec/2, &
                          scf%energy, mp2, err, n_frozen=core_orbital_count(z))
      call mol%destroy()
      call fit%destroy()
      if (err%has_error()) return
      energy = mp2%same_spin + mp2%opposite_spin
      if (opts%scf_max_iter < 0) energy = 0.0_dp   ! never taken; keeps `opts` used
   end function ri_mp2_on

   subroutine dimer_rimp2(z, symbols, xyz, charge, opts, energy, correlation, err)
      !! `E_IJ^0` at RI-MP2: one RHF here, then the correlation on its orbitals
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: charge
      type(efmo_options_t), intent(in) :: opts
      real(dp), intent(out) :: energy, correlation
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: guess_density(:, :)
      integer :: guess_kind

      energy = 0.0_dp
      correlation = 0.0_dp
      call build_czt_molecule(z, symbols, xyz, BASIS, mol, err, force_cartesian=.true.)
      if (err%has_error()) return
      call build_restricted_guess(mol, "auto", guess_kind, guess_density, err)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if
      call run_czt_rhf(mol, sum(z) - charge, opts%scf_max_iter, opts%scf_energy_tol, &
                       opts%scf_density_tol, .false., scf, err, &
                       guess=guess_kind, guess_density=guess_density, &
                       grad_tol=opts%scf_grad_tol, scf=opts%scf)
      call mol%destroy()
      if (err%has_error()) return
      correlation = ri_mp2_on(z, symbols, xyz, sum(z) - charge, scf, opts, err)
      if (err%has_error()) return
      energy = scf%energy + correlation
   end subroutine dimer_rimp2

   subroutine test_rimp2_dimer(error)
      !! Two fragments at RI-MP2: EFMO is the dimer's own correlated energy
      !!
      !! Limit three again, with `model.method: ri-mp2`. It stays exact for the
      !! same reason it was exact at Hartree-Fock -- on two fragments
      !! `E_IJ^pol` *is* `E_pol^total`, so the induction cancels and the
      !! monomer energies telescope out of `E_IJ^0 - E_I^0 - E_J^0` -- and the
      !! correlation rides along, because a correlated `E_I^0` is still just
      !! `E_I^0`.
      !!
      !! **It is the sharpest test of the correlation plumbing there is.** The
      !! monomer correlation is computed on the orbitals `make_efp_potential`
      !! handed back and the dimer's on a fresh SCF's, so a run that correlated
      !! the wrong determinant, froze a different core on one side, or fitted
      !! against a differently built auxiliary molecule fails here by a
      !! millihartree rather than passing with a plausible number.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      integer :: z(6), owner(6)
      character(len=2) :: symbols(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: reference, reference_corr

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call water_chain(2, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 1.0e6_dp
      opts%correlation = EFMO_CORR_RI_MP2
      opts%corr_aux_basis = AUX
      call run_efmo(z, symbols, xyz, owner, [0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      call dimer_rimp2(z, symbols, xyz, 0, opts, reference, reference_corr, err)
      call check(error,.not. err%has_error(), "the reference RI-MP2 failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%energy, reference, thr=TOL, &
                 message="EFMO/RI-MP2 on two fragments is not the dimer's own "// &
                 "in-vacuo RI-MP2 energy")
      if (allocated(error)) return
      ! The correlation is reported broken in two, and the two have to add up to
      ! the dimer's own: the monomers' plus (dimer - monomers).
      call check(error, res%monomer_correlation + res%nmer_correlation, &
                 reference_corr, thr=TOL, &
                 message="the reported correlation does not add up to the dimer's")
      if (allocated(error)) return
      ! And switching it off has to give the Hartree-Fock answer back exactly,
      ! not nearly: the correlated path must not have moved an SCF.
      opts%correlation = EFMO_CORR_NONE
      call run_efmo(z, symbols, xyz, owner, [0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, res%energy, cached_dimer_12, thr=TOL, &
                 message="switching the correlation off did not give the RHF total back")
      if (allocated(error)) return
      call check(error, res%monomer_correlation, 0.0_dp, thr=0.0_dp, &
                 message="a Hartree-Fock run reported a correlation energy")
   end subroutine test_rimp2_dimer

   subroutine test_rimp2_trimer(error)
      !! Three fragments, every pair quantum: the correlated two-body sum
      !!
      !! With `R_cut` huge eq 6 has no effective-fragment half, so
      !!
      !!     E = sum_I E_I + sum_{I<J} (E_IJ - E_I - E_J)
      !!         + (E_pol^total - sum_IJ E_IJ^pol)
      !!
      !! -- a many-body-expansion pair sum at RI-MP2, plus the part of the
      !! induction no pair holds. Every monomer and every dimer energy on the
      !! right is computed here, independently, from an SCF plus a fitted MP2,
      !! so what is asserted is that the orchestrator's correlated `E_I^0` and
      !! `E_IJ^0` are those numbers and that they enter the expansion with the
      !! signs Hartree-Fock's do.
      !!
      !! The induction terms are taken from the run rather than recomputed: they
      !! are Hartree-Fock quantities built from the potentials and are the same
      !! whatever the correlation is, which is itself the assertion made by
      !! comparing them against `test_all_quantum`'s.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      integer :: z(9), owner(9)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9)
      real(dp) :: dimer(3), dimer_corr(3), expected, pair_sum
      integer :: pairs(2, 3), k
      integer, allocatable :: idx(:)

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 1.0e6_dp
      opts%correlation = EFMO_CORR_RI_MP2
      opts%corr_aux_basis = AUX
      call run_efmo(z, symbols, xyz, owner, [0, 0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, res%n_qm_pairs, 3, message="not every pair is quantum")
      if (allocated(error)) return

      ! The monomer sum is the three cached Hartree-Fock energies plus their
      ! three correlations, each computed on the potential's own orbitals.
      call check(error, res%monomer_sum, sum(cached_mono) + sum(cached_mono_mp2), &
                 thr=TOL, message="sum E_I^0 is not the correlated monomer sum")
      if (allocated(error)) return
      call check(error, res%monomer_correlation, sum(cached_mono_mp2), thr=TOL, &
                 message="the reported monomer correlation is not the sum of the three")
      if (allocated(error)) return

      pairs = reshape([1, 2, 1, 3, 2, 3], [2, 3])
      pair_sum = 0.0_dp
      do k = 1, 3
         idx = [atoms_of(owner, pairs(1, k)), atoms_of(owner, pairs(2, k))]
         call dimer_rimp2(z(idx), symbols(idx), xyz(:, idx), 0, opts, dimer(k), &
                          dimer_corr(k), err)
         if (err%has_error()) exit
         pair_sum = pair_sum + dimer(k) &
                    - (cached_mono(pairs(1, k)) + cached_mono_mp2(pairs(1, k))) &
                    - (cached_mono(pairs(2, k)) + cached_mono_mp2(pairs(2, k)))
      end do
      call check(error,.not. err%has_error(), "a reference dimer failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%nmer_correction, pair_sum, thr=TOL, &
                 message="the quantum dimer correction is not the correlated pair sum")
      if (allocated(error)) return

      expected = sum(cached_mono) + sum(cached_mono_mp2) + pair_sum &
                 - res%induction_correction + res%polarization_total
      call check(error, res%energy, expected, thr=TOL, &
                 message="EFMO/RI-MP2 with every pair quantum is not the correlated "// &
                 "pair sum plus the many-body induction")
   end subroutine test_rimp2_trimer

   subroutine test_level_two_is_the_pair_form(error)
      !! Level two through the general machinery is the pair form, bit for bit
      !!
      !! The energy is now assembled by a difference operator over subsets, and
      !! at level two that operator has to reproduce the subtraction the pair
      !! expression was written as -- `E_IJ^0 - E_I^0 - E_J^0` and `E_IJ^pol` --
      !! in the same floating-point order, or every EFMO reference in this
      !! repository moves in its last digits for no physical reason.
      !!
      !! **So the tolerance is zero.** Both sums are recomputed here from the
      !! run's own reported per-pair numbers, accumulated in the order the pair
      !! table is in, and required to be exactly equal. A tolerance here would
      !! test nothing that the four limit tests above do not already test.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      integer :: z(9), owner(9)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9), pair_form, pol_form, total
      integer :: k

      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 1.0e6_dp
      opts%level = 2
      call run_efmo(z, symbols, xyz, owner, [0, 0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%n_qm_groups, 3, &
                 message="level two on three fragments is three groups")
      if (allocated(error)) return
      call check(error, size(res%level_vacuum), 2, &
                 message="a level-two run reported more than two levels")
      if (allocated(error)) return

      pair_form = 0.0_dp
      pol_form = 0.0_dp
      do k = 1, size(res%pairs)
         if (.not. res%pairs(k)%qm) cycle
         ! Parenthesised on purpose: the many-body difference is formed per
         ! pair and *then* accumulated, so the reference here has to associate
         ! the same way to be comparable at zero tolerance.
         pair_form = pair_form + (res%pairs(k)%e_dimer &
                                  - res%monomer_energy(res%pairs(k)%i) &
                                  - res%monomer_energy(res%pairs(k)%j))
         pol_form = pol_form + res%pairs(k)%e_pair_pol
      end do

      call check(error, res%nmer_correction, pair_form, thr=0.0_dp, &
                 message="the many-body difference at level two is not exactly the "// &
                 "pair subtraction")
      if (allocated(error)) return
      call check(error, res%level_vacuum(2), pair_form, thr=0.0_dp, &
                 message="the level-two vacuum sum is not exactly the pair subtraction")
      if (allocated(error)) return
      call check(error, res%induction_correction, pol_form, thr=0.0_dp, &
                 message="the induction difference at level two is not exactly the "// &
                 "sum of the pair inductions")
      if (allocated(error)) return
      call check(error, res%level_vacuum(1), res%monomer_sum, thr=0.0_dp, &
                 message="level one is not the fragment sum")
      if (allocated(error)) return
      call check(error, res%level_induction(1), 0.0_dp, thr=0.0_dp, &
                 message="a single fragment was given an induction energy")
      if (allocated(error)) return

      ! And that the total is those sums and nothing else, again exactly.
      total = res%monomer_sum + res%nmer_correction - res%induction_correction &
              + res%far_electrostatics + res%far_dispersion &
              + res%far_exchange_repulsion + res%far_charge_transfer &
              + res%polarization_total
      call check(error, res%energy, total, thr=0.0_dp, &
                 message="the total is not the reported sums")
   end subroutine test_level_two_is_the_pair_form

   subroutine test_full_level(error)
      !! Level = N with `R_cut` huge is the unfragmented energy of the cluster
      !!
      !! **The strongest statement the method makes, and the sharpest test of
      !! the subset-induction bookkeeping.** With every group near and the level
      !! at the fragment count, two series each telescope:
      !!
      !!   * `sum_S dE_S^0` over every subset is the in-vacuo energy of the
      !!     whole cluster -- one ordinary RHF, computed here independently.
      !!   * `sum_S dE_S^pol` over every subset is `E_pol^total`, the last term
      !!     of the energy, which it therefore cancels **entirely**.
      !!
      !! So EFMO at full level is not "close to" the supersystem: it *is* the
      !! supersystem's RHF energy, to SCF convergence, with no induction
      !! correction left over at all. A missing subset, a sign, or a group
      !! induction solved over the wrong fragments all leave the level-two
      !! limits intact and break this.
      !!
      !! Three waters, so the run is three MAKEFPs, three dimer SCFs and one
      !! trimer SCF.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      integer :: z(9), owner(9)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9), supersystem

      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 1.0e6_dp
      opts%level = 3
      call run_efmo(z, symbols, xyz, owner, [0, 0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      ! Three pairs and one trimer: the whole subset lattice above the monomers.
      call check(error, res%n_qm_groups, 4, &
                 message="full level on three fragments is not four groups")
      if (allocated(error)) return
      call check(error, res%level_count(3), 1, message="the trimer was not enumerated")
      if (allocated(error)) return

      ! The polarization correction cancels the total induction exactly, which
      ! is the induction series telescoping. Held to the induced-dipole solve's
      ! own tolerance rather than to the SCF's: no SCF enters either side.
      call check(error, res%induction_correction, res%polarization_total, &
                 thr=1.0e-12_dp, &
                 message="the induction series does not sum to E_pol^total, so the "// &
                 "polarization correction does not cancel at full level")
      if (allocated(error)) return
      ! And that the three-body induction is not zero, so the cancellation above
      ! is an identity being satisfied and not two zeros agreeing.
      call check(error, abs(res%level_induction(3)) > 1.0e-6_dp, &
                 "the three-body induction vanished, so this geometry cannot tell a "// &
                 "telescoping series from a dropped one")
      if (allocated(error)) return

      call dimer_rhf(z, symbols, xyz, 0, opts, supersystem, err)
      call check(error,.not. err%has_error(), "the supersystem RHF failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call check(error, res%energy, supersystem, thr=TOL, &
                 message="EFMO at level = N with a huge cutoff is not the "// &
                 "unfragmented RHF energy of the cluster")
   end subroutine test_full_level

   subroutine test_induction_series(error)
      !! `sum_S dE_S^pol` over every subset of a trimer is `E_pol^total`
      !!
      !! The induction half of the identity above, asserted on the potentials
      !! alone -- no SCF, no energy, no cutoff -- so that a failure in the
      !! full-level test names which of the two series moved. Written out here
      !! against the subset solver directly:
      !!
      !!     dE_I^pol = 0
      !!     dE_IJ^pol = E_IJ^pol
      !!     dE_IJK^pol = E_IJK^pol - E_IJ^pol - E_IK^pol - E_JK^pol
      !!
      !! whose sum is `E_IJK^pol`, the induction over all three at once, because
      !! every pair term appears once with each sign.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: shifts(3, 3), series, three_body, total
      integer :: a, b

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      shifts = 0.0_dp
      series = 0.0_dp
      three_body = subset_polarization_energy(cached_frag, shifts, [1, 2, 3], err)
      do a = 1, 2
         do b = a + 1, 3
            series = series + subset_polarization_energy(cached_frag, shifts, [a, b], err)
            three_body = three_body &
                         - subset_polarization_energy(cached_frag, shifts, [a, b], err)
         end do
      end do
      series = series + three_body
      total = total_induction(cached_frag, err)
      call check(error,.not. err%has_error(), "an induction solve failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call check(error, series, total, thr=1.0e-12_dp, &
                 message="the induction expansion does not telescope to the "// &
                 "induction over all three fragments")
      if (allocated(error)) return

      ! And that the pair entry of the series is the two-fragment routine to the
      ! last bit, so the general solver has not changed what a pair means.
      call check(error, subset_polarization_energy(cached_frag, shifts, [1, 2], err), &
                 pair_polarization_energy(cached_frag(1), cached_frag(2), shifts(:, 1), &
                                          shifts(:, 2), err), thr=0.0_dp, &
                 message="the subset induction on two fragments is not exactly the "// &
                 "pair induction")
   end subroutine test_induction_series

   subroutine test_two_fragments(error)
      !! Limit three: two fragments, so EFMO is the dimer's own RHF energy
      !!
      !! `E_IJ^pol` on a two-fragment system *is* `E_pol^total` -- the same
      !! solver on the same system -- so the last term of eq 6 and the
      !! subtraction inside the near sum cancel exactly, and what is left is
      !! `E_I^0 + E_J^0 + (E_IJ^0 - E_I^0 - E_J^0) = E_IJ^0`.
      !!
      !! The sharpest of the four: it says the induction bookkeeping is right
      !! *and* that the monomer energies entering the difference are the ones
      !! the potentials were built from. Off by a wrong monomer energy, this
      !! fails by the size of an SCF energy rather than of an interaction.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      integer :: z(6), owner(6)
      character(len=2) :: symbols(6)
      real(dp) :: xyz(3, 6)

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call water_chain(2, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 1.0e6_dp
      call run_efmo(z, symbols, xyz, owner, [0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%n_qm_pairs, 1, message="the one pair is not quantum")
      if (allocated(error)) return
      call check(error, res%induction_correction, res%polarization_total, thr=1.0e-12_dp, &
                 message="on two fragments the pair induction is not the total")
      if (allocated(error)) return
      call check(error, res%energy, cached_dimer_12, thr=TOL, &
                 message="EFMO on two fragments is not the dimer's in-vacuo RHF energy")
   end subroutine test_two_fragments

   subroutine test_mixed(error)
      !! Limit four: one quantum dimer and two effective ones, assembled by hand
      !!
      !! The trimer's separations are 1.32, 3.95 and 2.63 in contact units, so
      !! `R_cut = 2.0` puts the close pair in the quantum list and the other two
      !! in the effective one. This is the only case where both halves of eq 6
      !! are non-empty, and it is where a sign or a double count that cancels in
      !! the limits shows up.
      !!
      !! Every piece of the expected total comes from a Phase 1 routine called
      !! directly here: the monomer energies from the cached potentials, the
      !! dimer from an independent RHF, the pair induction and the four far
      !! terms from their own entries, and the total induction over all three.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      type(efp_pair_energy_t), allocatable :: far(:)
      real(dp) :: shifts(3, 3), zero(3), e_pair_pol, e_pol_total, expected
      integer :: z(9), owner(9), far_pairs(2, 2)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9)

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 2.0_dp
      call run_efmo(z, symbols, xyz, owner, [0, 0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%n_qm_pairs, 1, message="the split did not keep one QM dimer")
      if (allocated(error)) return
      call check(error, res%n_efp_pairs, 2, message="the split did not leave two EFP dimers")
      if (allocated(error)) return
      call check(error, res%pairs(1)%i == 1 .and. res%pairs(1)%j == 2, &
                 "the quantum dimer is not the close pair")
      if (allocated(error)) return
      ! The separations this geometry was built for, so a change to the radii or
      ! to the chain would fail here rather than silently reclassify a pair.
      call check(error, res%pairs(1)%r, SPACING_12/(2.0_dp*1.40_dp), thr=1.0e-10_dp, &
                 message="R_12 is not the oxygen separation over twice the GAMESS FMO oxygen radius")
      if (allocated(error)) return

      shifts = 0.0_dp
      zero = 0.0_dp
      far_pairs = reshape([1, 3, 2, 3], [2, 2])
      far = efp_pair_terms(cached_frag, shifts, far_pairs, err, charge_transfer_on=.true.)
      e_pair_pol = pair_polarization_energy(cached_frag(1), cached_frag(2), zero, zero, err)
      e_pol_total = total_induction(cached_frag, err)
      call check(error,.not. err%has_error(), "the reference terms failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      expected = sum(cached_mono) &
                 + (cached_dimer_12 - cached_mono(1) - cached_mono(2) - e_pair_pol) &
                 + sum(far%total) + e_pol_total
      call check(error, res%energy, expected, thr=TOL, &
                 message="the trimer total is not the sum assembled from the pieces")
      if (allocated(error)) return

      ! And each of the six sums separately, so a failure names its term.
      call check(error, res%monomer_sum, sum(cached_mono), thr=TOL, message="sum E_I^0")
      if (allocated(error)) return
      call check(error, res%nmer_correction, &
                 cached_dimer_12 - cached_mono(1) - cached_mono(2), thr=TOL, &
                 message="the QM dimer correction")
      if (allocated(error)) return
      call check(error, res%induction_correction, e_pair_pol, thr=1.0e-12_dp, &
                 message="sum E_IJ^pol")
      if (allocated(error)) return
      call check(error, res%far_electrostatics, sum(far%electrostatics), thr=1.0e-12_dp, &
                 message="the far Coulomb sum")
      if (allocated(error)) return
      call check(error, res%far_dispersion, sum(far%dispersion), thr=1.0e-12_dp, &
                 message="the far dispersion sum")
      if (allocated(error)) return
      call check(error, res%far_exchange_repulsion, sum(far%exchange_repulsion), &
                 thr=1.0e-12_dp, message="the far exchange repulsion sum")
      if (allocated(error)) return
      call check(error, res%far_charge_transfer, sum(far%charge_transfer), thr=1.0e-12_dp, &
                 message="the far charge transfer sum")
      if (allocated(error)) return
      call check(error, res%polarization_total, e_pol_total, thr=1.0e-12_dp, &
                 message="E_pol^total")
   end subroutine test_mixed

end module test_mqc_czt_efmo

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_efmo, only: collect_mqc_czt_efmo_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_efmo", collect_mqc_czt_efmo_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
