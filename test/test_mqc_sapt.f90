!! SAPT0's molecules: the dimer basis, and ghosts in it
module test_mqc_sapt
   !! M1 of the SAPT ladder. What is asserted here is everything that can be
   !! checked before a single term is computed, and each assertion corresponds to
   !! a mistake that is otherwise silent:
   !!
   !!   * a ghost that kept its nuclear charge -- caught by the charge sum
   !!   * a ghost that lost its basis functions -- caught by the AO count
   !!   * a monomer built in its own atom order rather than the dimer's -- caught
   !!     by `Tr(D S)`, which is the occupied count only when the AO orderings
   !!     agree. This one is not hypothetical; the PySCF prototype was written
   !!     that way first and every AO quantity was quietly permuted.
   !!   * ghosting that changed the nuclear repulsion of the real atoms -- caught
   !!     by comparing against the monomer alone
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_sapt, only: sapt_molecules_t, build_sapt_molecules, &
                       sapt_cache_t, build_sapt_cache, sapt_elst10, sapt_exch10_s2, sapt_exch10, &
                       sapt_induction, sapt_terms_t, sapt_disp20, sapt_exch_disp20, run_sapt0, &
                       sapt2_cache_t, build_sapt2_cache, sapt2_amps_t, build_sapt2_amps, &
                       sapt2_zero_amps, sapt2_amps_mp2_energy, sapt2_k2f, sapt_elst12, &
                       sapt_exch11, sapt_exch12, sapt_ind22, run_sapt2
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_mp2, only: run_libcint_mp2, mp2_result_t
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_sapt_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_sapt_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("sapt_dimer_basis_and_ghosts", test_ghosts), &
                  new_unittest("sapt_monomer_ao_ordering", test_ordering), &
                  new_unittest("sapt_bsse_is_negative", test_bsse), &
                  new_unittest("sapt_elst10", test_elst10), &
                  new_unittest("sapt_exch10_s2", test_exch10_s2), &
                  new_unittest("sapt_exch10_sinf", test_exch10), &
                  new_unittest("sapt_induction", test_induction), &
                  new_unittest("sapt_disp20", test_disp20), &
                  new_unittest("sapt_exch_disp20", test_exch_disp20), &
                  new_unittest("sapt_total_and_dhf", test_total), &
                  new_unittest("sapt2_amplitudes_mp2", test_sapt2_amplitudes), &
                  new_unittest("sapt2_terms", test_sapt2_terms), &
                  new_unittest("sapt2_total", test_sapt2_total) &
                  ]
   end subroutine collect_mqc_sapt_tests

   subroutine water_dimer(mols, err)
      !! Two waters 3.0 Angstrom apart along x -- the geometry every EFP number
      !! in this repository is pinned to, so the two decompositions can one day
      !! be read against each other on a case where the EFP side is exact.
      type(sapt_molecules_t), intent(out) :: mols
      type(error_t), intent(inout) :: err

      integer :: z(3)
      character(len=2) :: sym(3)
      real(dp) :: a(3, 3), b(3, 3)

      z = [8, 1, 1]
      sym = ["O ", "H ", "H "]
      a = reshape([0.0_dp, 0.0_dp, 0.10077199_dp, &
                   0.0_dp, 0.77250895_dp, -0.46780200_dp, &
                   0.0_dp, -0.77250895_dp, -0.46780200_dp], [3, 3])*ANG
      b = a
      b(1, :) = b(1, :) + 3.0_dp*ANG
      call build_sapt_molecules(z, sym, a, z, sym, b, "6-31g", mols, err)
   end subroutine water_dimer

   subroutine test_ghosts(error)
      !! Both monomers span the dimer basis, and only their own atoms are charged
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(error_t) :: err

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return

      call check(error, mols%mono_a%nao == mols%dimer%nao, &
                 "monomer A must span the whole dimer basis")
      if (allocated(error)) return
      call check(error, mols%mono_b%nao == mols%dimer%nao, &
                 "monomer B must span the whole dimer basis")
      if (allocated(error)) return

      ! Charges: A's three atoms carry 10, B's ghosts carry nothing.
      call check(error, abs(sum(mols%mono_a%charges(1:3)) - 10.0_dp) < 1.0e-12_dp, &
                 "monomer A's own atoms lost their nuclear charge")
      if (allocated(error)) return
      call check(error, abs(sum(mols%mono_a%charges(4:6))) < 1.0e-12_dp, &
                 "monomer A's ghosts kept a nuclear charge")
      if (allocated(error)) return
      call check(error, abs(sum(mols%mono_b%charges(1:3))) < 1.0e-12_dp, &
                 "monomer B's ghosts kept a nuclear charge")
      if (allocated(error)) return

      ! And ghosting must not disturb the real atoms' mutual repulsion.
      call check(error, mols%mono_a%nuclear_repulsion() > 0.0_dp, &
                 "monomer A should still have a nuclear repulsion")
      if (allocated(error)) return
      call check(error, mols%dimer%nuclear_repulsion() > &
                 mols%mono_a%nuclear_repulsion() + mols%mono_b%nuclear_repulsion(), &
                 "the dimer's nuclear repulsion must exceed the monomers' sum")

      call mols%destroy()
   end subroutine test_ghosts

   subroutine test_ordering(error)
      !! `Tr(D S) = n_occ`, which holds only if the AO orderings agree
      !!
      !! `D` comes from a monomer's SCF and `S` from the dimer. If the monomer
      !! were built in its own atom order the two would be permuted relative to
      !! each other, the trace would come out wrong, and every SAPT term after it
      !! would be quietly meaningless.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(error_t) :: err
      type(rhf_result_t) :: scf
      real(dp), allocatable :: s(:, :), d(:, :), ds(:, :)
      real(dp) :: trace
      integer :: i, nocc

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return

      call run_libcint_rhf(mols%mono_a, mols%n_elec_a, 200, 1.0e-11_dp, 1.0e-9_dp, &
                           .false., scf, err, in_core=.true.)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "monomer A's SCF failed in the dimer basis")
      if (allocated(error)) return

      nocc = mols%n_elec_a/2
      call mols%dimer%overlap(s)
      allocate (d(mols%dimer%nao, mols%dimer%nao))
      allocate (ds(mols%dimer%nao, mols%dimer%nao))
      call pic_gemm(scf%orbitals(:, 1:nocc), scf%orbitals(:, 1:nocc), d, transb="T")
      call pic_gemm(d, s, ds)
      trace = 0.0_dp
      do i = 1, mols%dimer%nao
         trace = trace + ds(i, i)
      end do

      call check(error, abs(trace - real(nocc, dp)) < 1.0e-8_dp, &
                 "Tr(D S) is not the occupied count -- the monomer's AO ordering "// &
                 "does not match the dimer's")

      deallocate (s, d, ds)
      call mols%destroy()
   end subroutine test_ordering

   subroutine test_bsse(error)
      !! The dimer basis must lower a monomer's energy, and by a sane amount
      !!
      !! No wrong ghost has this property: one that kept its nuclear charge gives
      !! a wildly lower energy, and one that lost its basis functions gives
      !! exactly the monomer's own energy.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(libcint_molecule_t) :: alone
      type(rhf_result_t) :: dcbs, own
      type(error_t) :: err
      integer :: z(3)
      character(len=2) :: sym(3)
      real(dp) :: a(3, 3), bsse

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return

      z = [8, 1, 1]
      sym = ["O ", "H ", "H "]
      a = reshape([0.0_dp, 0.0_dp, 0.10077199_dp, &
                   0.0_dp, 0.77250895_dp, -0.46780200_dp, &
                   0.0_dp, -0.77250895_dp, -0.46780200_dp], [3, 3])*ANG
      call build_libcint_molecule(z, sym, a, "6-31g", alone, err)
      call check(error,.not. err%has_error(), "building the lone monomer failed")
      if (allocated(error)) return

      call run_libcint_rhf(mols%mono_a, mols%n_elec_a, 200, 1.0e-11_dp, 1.0e-9_dp, &
                           .false., dcbs, err, in_core=.true.)
      if (.not. err%has_error()) then
         call run_libcint_rhf(alone, mols%n_elec_a, 200, 1.0e-11_dp, 1.0e-9_dp, &
                              .false., own, err, in_core=.true.)
      end if
      call check(error,.not. err%has_error(), "an SCF failed")
      if (allocated(error)) return

      bsse = dcbs%energy - own%energy
      call check(error, bsse < 0.0_dp, &
                 "the dimer basis must lower the monomer's energy")
      if (allocated(error)) return
      ! Loose upper bound: a ghost that kept its charge is worth whole Hartrees,
      ! not milli-Hartrees, so this separates the two cases by orders of magnitude.
      call check(error, bsse > -0.1_dp, &
                 "the basis-set superposition error is far too large, which is "// &
                 "what a ghost that kept its nuclear charge looks like")

      call alone%destroy()
      call mols%destroy()
   end subroutine test_bsse

   subroutine test_elst10(error)
      !! `Elst10,r` against the conventional four-index reference
      !!
      !! The reference is `validation/check_sapt0.py`, which computes the same
      !! quantity in PySCF with exact integrals. psi4 cannot serve here: its
      !! closed-shell SAPT is density-fitted by construction, so its number is a
      !! different quantity in the same basis and agrees only to ~1e-6.
      !!
      !! Regenerate with:  python validation/check_sapt0.py --json <path>
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(error_t) :: err
      real(dp) :: e

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      e = sapt_elst10(cache)
      ! Conventional four-index, water dimer 3.0 Angstrom along x, 6-31G.
      call check(error, e, 0.006384228637_dp, thr=1.0e-9_dp, &
                 message="Elst10,r disagrees with the PySCF reference")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_elst10

   subroutine test_exch10_s2(error)
      !! `Exch10(S^2)` against the conventional four-index reference
      !!
      !! The term that exposed psi4's `vector_dot` being elementwise rather than
      !! a trace: two of its operands are non-symmetric, and taking `Tr(X Y)`
      !! there puts the answer 55% high while leaving every shape valid.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(error_t) :: err
      real(dp) :: e

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      e = sapt_exch10_s2(cache)
      call check(error, e, 0.003235393978_dp, thr=1.0e-9_dp, &
                 message="Exch10(S^2) disagrees with the PySCF reference")
      if (allocated(error)) return
      ! First-order exchange is repulsive. A sign slip here still produces a
      ! plausible-looking interaction energy once it is summed with the rest.
      call check(error, e > 0.0_dp, "first-order exchange must be repulsive")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_exch10_s2

   subroutine test_exch10(error)
      !! `Exch10` at S^inf -- the form the SAPT0 total actually uses
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(error_t) :: err
      real(dp) :: e, e_s2

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      e = sapt_exch10(cache, err)
      call check(error,.not. err%has_error(), "Exch10 failed")
      if (allocated(error)) return
      call check(error, e, 0.003241500408_dp, thr=1.0e-9_dp, &
                 message="Exch10 disagrees with the PySCF reference")
      if (allocated(error)) return

      ! The S^2 form is the leading term of this one, so at a separation where
      ! the overlap is small the two must be close and S^inf the larger. They
      ! part company as the monomers approach -- 19% at 2 Angstrom on LiH.
      e_s2 = sapt_exch10_s2(cache)
      call check(error, e > e_s2, "S^inf exchange should exceed the S^2 form here")
      if (allocated(error)) return
      call check(error, abs(e - e_s2) < 0.05_dp*abs(e), &
                 "S^2 and S^inf differ by far more than the overlap warrants")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_exch10

   subroutine test_induction(error)
      !! `Ind20` and `Exch-Ind20`, uncoupled and coupled
      !!
      !! The uncoupled pair is asserted as well as the coupled one, because it
      !! needs no response solver: if it passes and the coupled one does not, the
      !! fault is in `cphf_solve`'s conventions rather than in the contraction,
      !! and the two are otherwise indistinguishable from one failing number.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(sapt_terms_t) :: t
      type(error_t) :: err

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      call sapt_induction(mols, cache, t, err)
      call check(error,.not. err%has_error(), "induction failed")
      if (allocated(error)) return

      call check(error, t%ind20_u, -0.000954737579_dp, thr=1.0e-9_dp, &
                 message="Ind20,u disagrees with the PySCF reference")
      if (allocated(error)) return
      call check(error, t%exch_ind20_u, 0.000805018637_dp, thr=1.0e-9_dp, &
                 message="Exch-Ind20,u disagrees with the PySCF reference")
      if (allocated(error)) return
      call check(error, t%ind20_r, -0.001134897260_dp, thr=1.0e-8_dp, &
                 message="Ind20,r disagrees with the PySCF reference")
      if (allocated(error)) return
      call check(error, t%exch_ind20_r, 0.000949574150_dp, thr=1.0e-8_dp, &
                 message="Exch-Ind20,r disagrees with the PySCF reference")
      if (allocated(error)) return

      ! Induction is stabilising and its exchange counterpart repulsive, and
      ! letting the monomer relax can only lower the energy further.
      call check(error, t%ind20_r < t%ind20_u, &
                 "the response must be more stabilising than the uncoupled form")
      if (allocated(error)) return
      call check(error, t%exch_ind20_r > 0.0_dp, &
                 "exchange-induction must be repulsive")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_induction

   subroutine test_disp20(error)
      !! `Disp20` against the conventional four-index reference
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(error_t) :: err
      real(dp) :: e

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      e = sapt_disp20(cache)
      call check(error, e, -0.000317465444_dp, thr=1.0e-9_dp, &
                 message="Disp20 disagrees with the PySCF reference")
      if (allocated(error)) return
      ! Dispersion is attractive. The denominator is negative throughout, so a
      ! positive answer means the occupied and virtual energies were swapped.
      call check(error, e < 0.0_dp, "dispersion must be attractive")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_disp20

   subroutine test_exch_disp20(error)
      !! `Exch-Disp20` against the conventional four-index reference
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(error_t) :: err
      real(dp) :: e, d

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      e = sapt_exch_disp20(cache)
      call check(error, e, 0.000117840542_dp, thr=1.0e-9_dp, &
                 message="Exch-Disp20 disagrees with the PySCF reference")
      if (allocated(error)) return

      ! Repulsive, and a sizeable fraction of dispersion -- 10-25% across psi4's
      ! own tests, always opposite in sign. A term that quietly evaluated to zero
      ! would leave a total that is wrong in a systematic, plausible direction.
      d = sapt_disp20(cache)
      call check(error, e > 0.0_dp, "exchange-dispersion must be repulsive")
      if (allocated(error)) return
      call check(error, abs(e) > 0.05_dp*abs(d) .and. abs(e) < 0.6_dp*abs(d), &
                 "exchange-dispersion is not a plausible fraction of dispersion")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_exch_disp20

   subroutine test_total(error)
      !! `dHF(2)`, the total, and the identity that validates four terms at once
      !!
      !! `dHF` is *defined* as the residual, so
      !!
      !!     Elst10 + Exch10 + Ind20,r + Exch-Ind20,r + dHF == E_int^HF(CP)
      !!
      !! holds by construction. That makes it a check on those four *collectively*
      !! against a supermolecular Hartree-Fock number, to machine precision, with
      !! no SAPT reference involved -- the sharpest thing available here. It costs
      !! one dimer SCF, the monomers having already been run in the dimer basis.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_terms_t) :: t
      type(error_t) :: err
      real(dp) :: four, assembled

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return

      call run_sapt0(mols, t, err)
      call check(error,.not. err%has_error(), "SAPT0 failed")
      if (allocated(error)) return

      ! The identity, first: it needs no reference and localises a fault to the
      ! four terms as a set before any of the stored numbers are consulted.
      four = t%elst10 + t%exch10 + t%ind20_r + t%exch_ind20_r
      call check(error, abs(four + t%delta_hf - t%e_int_hf) < 1.0e-12_dp, &
                 "the four terms plus dHF must be the supermolecular HF energy")
      if (allocated(error)) return

      call check(error, t%e_int_hf, 0.009315543651_dp, thr=1.0e-8_dp, &
                 message="E_int^HF disagrees with the PySCF reference")
      if (allocated(error)) return
      call check(error, t%delta_hf, -0.000124862284_dp, thr=1.0e-8_dp, &
                 message="dHF(2) disagrees with the PySCF reference")
      if (allocated(error)) return
      call check(error, t%total, 0.009115918750_dp, thr=1.0e-8_dp, &
                 message="the SAPT0 total disagrees with the PySCF reference")
      if (allocated(error)) return

      ! And that the total really is the sum of what it claims, with Exch10 and
      ! not Exch10(S^2) -- a total built from the S^2 form is comparable to
      ! nothing published, and the two differ by only 6e-6 here so the wrong one
      ! looks entirely reasonable.
      assembled = four + t%delta_hf + t%disp20 + t%exch_disp20
      call check(error, abs(assembled - t%total) < 1.0e-14_dp, &
                 "the total is not the sum of its stated parts")
      if (allocated(error)) return
      call check(error, abs(t%exch10 - t%exch10_s2) > 1.0e-7_dp, &
                 "S^2 and S^inf exchange should differ measurably, so that using "// &
                 "the wrong one in the total would be detectable")

      call mols%destroy()
   end subroutine test_total

   subroutine test_sapt2_amplitudes(error)
      !! The SAPT2 amplitude machinery, checked against code already trusted
      !!
      !! Two identities, neither needing a SAPT2 reference:
      !!
      !! * the monomer MP2 correlation energy recovered from the amplitudes
      !!   must equal `run_libcint_mp2` on the same ghosted monomer to machine
      !!   precision (the plan's 0.2) -- this pins the amplitude conventions
      !!   (denominator in, nothing antisymmetrized) before any term consumes
      !!   them;
      !! * the `vAR`/`vBS` potentials assembled through psi4's thirteen dressed
      !!   contractions must reproduce `Exch-Ind20,r` when contracted with the
      !!   CPHF coefficients -- which validates the entire dressed three-index
      !!   machinery against the USAPT0-factorized potential SAPT0 already
      !!   pinned to PySCF.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(sapt2_cache_t) :: c2
      type(sapt2_amps_t) :: amps_a, amps_b
      type(sapt_terms_t) :: t
      type(mp2_result_t) :: mp2
      type(error_t) :: err
      real(dp), allocatable :: chf_a(:, :), chf_b(:, :)
      real(dp), allocatable :: u_ar(:, :), v_ar(:, :), u_bs(:, :), v_bs(:, :)
      real(dp) :: e2, exind

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return
      call build_sapt2_cache(cache, c2, err)
      call check(error,.not. err%has_error(), "the SAPT2 cache failed")
      if (allocated(error)) return

      call build_sapt2_amps(c2, "A", amps_a)
      call build_sapt2_amps(c2, "B", amps_b)

      ! The MP2 identity, against an implementation that shares none of the
      ! amplitude machinery. Conventional water dimer basis, monomer A.
      e2 = sapt2_amps_mp2_energy(c2, amps_a, "A")
      call run_libcint_mp2(mols%mono_a, cache%c_a, cache%eps_a, cache%nocc_a, &
                           cache%e_scf_a, mp2, err)
      call check(error,.not. err%has_error(), "the reference MP2 failed")
      if (allocated(error)) return
      call check(error, abs(e2 - mp2%correlation) < 1.0e-11_dp, &
                 "the amplitudes' MP2 energy must match run_libcint_mp2")
      if (allocated(error)) return
      e2 = sapt2_amps_mp2_energy(c2, amps_b, "B")
      call check(error, abs(e2 - mp2%correlation) < 1.0e-11_dp, &
                 "monomer B is monomer A translated, so its E2 is the same")
      if (allocated(error)) return

      ! The dressed-machinery identity: -2 x . vAR summed over both
      ! directions is Exch-Ind20,r.
      call sapt_induction(mols, cache, t, err, chf_a=chf_a, chf_b=chf_b)
      call check(error,.not. err%has_error(), "induction failed")
      if (allocated(error)) return
      call sapt2_k2f(c2, u_ar, v_ar, u_bs, v_bs)
      exind = -2.0_dp*sum(chf_a*v_ar) - 2.0_dp*sum(chf_b*v_bs)
      call check(error, abs(exind - t%exch_ind20_r) < 1.0e-11_dp, &
                 "the K2f route to Exch-Ind20,r disagrees with the USAPT0 "// &
                 "factorization -- a dressing column is wrong somewhere")

      call amps_a%destroy()
      call amps_b%destroy()
      call c2%destroy()
      call cache%destroy()
      call mols%destroy()
   end subroutine test_sapt2_amplitudes

   subroutine test_sapt2_terms(error)
      !! The four SAPT2 terms against the conventional four-index reference
      !!
      !! The reference is `validation/check_sapt2.py`, the same conventional
      !! PySCF route the SAPT0 numbers come from; psi4 cannot serve, being
      !! density-fitted (and, in its own test suite, natural-orbital
      !! truncated) by construction.
      !!
      !! Regenerate with:  python validation/check_sapt2.py --json <path>
      !!
      !! The `t2 -> 0` identity is asserted alongside each value: with the
      !! amplitudes zeroed, `Elst12`, `Exch11` and `Exch12` vanish
      !! identically, and `Ind22`'s pieces 2-7 do -- its first piece builds
      !! the perturbed doubles from bare integrals times the uncoupled
      !! response and survives, which is a correction this project made to
      !! its own plan.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(sapt2_cache_t) :: c2
      type(sapt2_amps_t) :: amps_a, amps_b
      type(sapt_terms_t) :: t
      type(error_t) :: err
      real(dp), allocatable :: chf_a(:, :), chf_b(:, :)
      real(dp), allocatable :: u_ar(:, :), v_ar(:, :), u_bs(:, :), v_bs(:, :)
      real(dp) :: e, pieces(7, 2)

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return
      call build_sapt2_cache(cache, c2, err)
      call check(error,.not. err%has_error(), "the SAPT2 cache failed")
      if (allocated(error)) return
      call build_sapt2_amps(c2, "A", amps_a)
      call build_sapt2_amps(c2, "B", amps_b)
      call sapt_induction(mols, cache, t, err, chf_a=chf_a, chf_b=chf_b)
      call check(error,.not. err%has_error(), "induction failed")
      if (allocated(error)) return
      call sapt2_k2f(c2, u_ar, v_ar, u_bs, v_bs)

      ! Elst12,r. The tolerance is set by the iterative CPHF behind the
      ! orbital-relaxation contraction, as for Ind20,r above.
      e = sapt_elst12(c2, amps_a, amps_b, chf_a, chf_b)
      call check(error, e, -0.000735464482_dp, thr=1.0e-8_dp, &
                 message="Elst12,r disagrees with the PySCF reference")
      if (allocated(error)) return

      e = sapt_exch11(c2, amps_a, amps_b)
      call check(error, e, -0.000218219611_dp, thr=1.0e-9_dp, &
                 message="Exch11 disagrees with the PySCF reference")
      if (allocated(error)) return

      e = sapt_exch12(c2, amps_a, amps_b, u_ar, u_bs)
      call check(error, e, 0.001170593809_dp, thr=1.0e-9_dp, &
                 message="Exch12 disagrees with the PySCF reference")
      if (allocated(error)) return

      call sapt_ind22(c2, amps_a, amps_b, e, pieces)
      call check(error, e, -0.000467708358_dp, thr=1.0e-9_dp, &
                 message="Ind22 disagrees with the PySCF reference")
      if (allocated(error)) return
      ! Piece 5 is where the second-order doubles enter -- the one place t2
      ! is consumed, and worth 13% of Ind22 when the first-order t is used
      ! there by mistake.
      call check(error, pieces(5, 1), 0.000000546805_dp, thr=1.0e-11_dp, &
                 message="Ind22's t2 piece disagrees with the reference")
      if (allocated(error)) return

      ! The t2 -> 0 fallback, piecewise.
      call sapt2_zero_amps(amps_a)
      call sapt2_zero_amps(amps_b)
      e = sapt_elst12(c2, amps_a, amps_b, chf_a, chf_b)
      call check(error, e == 0.0_dp, "Elst12 must vanish with the amplitudes")
      if (allocated(error)) return
      e = sapt_exch11(c2, amps_a, amps_b)
      call check(error, e == 0.0_dp, "Exch11 must vanish with the amplitudes")
      if (allocated(error)) return
      e = sapt_exch12(c2, amps_a, amps_b, u_ar, u_bs)
      call check(error, e == 0.0_dp, "Exch12 must vanish with the amplitudes")
      if (allocated(error)) return
      call sapt_ind22(c2, amps_a, amps_b, e, pieces)
      call check(error, all(pieces(2:7, :) == 0.0_dp), &
                 "Ind22's pieces 2-7 must vanish with the amplitudes")
      if (allocated(error)) return
      call check(error, abs(pieces(1, 1)) > 1.0e-8_dp, &
                 "Ind22's first piece is integral-driven and must NOT vanish "// &
                 "-- if it does, its perturbed doubles lost their source term")

      call amps_a%destroy()
      call amps_b%destroy()
      call c2%destroy()
      call cache%destroy()
      call mols%destroy()
   end subroutine test_sapt2_terms

   subroutine test_sapt2_total(error)
      !! `run_sapt2`: the SAPT0 part is bit-compatible with `run_sapt0`, and
      !! the SAPT2 total is the sum it claims to be
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_terms_t) :: t0, t2
      type(error_t) :: err
      real(dp) :: assembled

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return

      call run_sapt0(mols, t0, err)
      call check(error,.not. err%has_error(), "SAPT0 failed")
      if (allocated(error)) return
      call run_sapt2(mols, t2, err)
      call check(error,.not. err%has_error(), "SAPT2 failed")
      if (allocated(error)) return

      ! The SAPT0 part of SAPT2 goes through exactly run_sapt0's code -- the
      ! same subroutine, on a cache built the same way -- so any difference
      ! here is the run-to-run reproducibility of the monomer SCFs, measured
      ! at ~1e-14 on the terms and ~1e-12 on the supermolecular residual
      ! (76-hartree totals differenced). A wiring error is orders of
      ! magnitude above this threshold.
      call check(error, abs(t2%total - t0%total) < 1.0e-11_dp, &
                 "SAPT2's SAPT0 part must reproduce run_sapt0")
      if (allocated(error)) return
      call check(error, abs(t2%elst10 - t0%elst10) < 1.0e-11_dp .and. &
                 abs(t2%exch10 - t0%exch10) < 1.0e-11_dp .and. &
                 abs(t2%ind20_r - t0%ind20_r) < 1.0e-11_dp .and. &
                 abs(t2%exch_ind20_r - t0%exch_ind20_r) < 1.0e-11_dp .and. &
                 abs(t2%disp20 - t0%disp20) < 1.0e-11_dp .and. &
                 abs(t2%exch_disp20 - t0%exch_disp20) < 1.0e-11_dp, &
                 "a SAPT0 term differs between run_sapt0 and run_sapt2")
      if (allocated(error)) return

      call check(error, t2%total_sapt2, 0.009256454004_dp, thr=1.0e-8_dp, &
                 message="the SAPT2 total disagrees with the PySCF reference")
      if (allocated(error)) return
      call check(error, t2%exch_ind22, 0.000391333896_dp, thr=1.0e-9_dp, &
                 message="Exch-Ind22 disagrees with the PySCF reference")
      if (allocated(error)) return

      ! Exch-Ind22 is defined as a scaling, not computed (ind22.cc:52).
      call check(error, abs(t2%exch_ind22 &
                            - t2%ind22*(t2%exch_ind20_r/t2%ind20_r)) &
                 < 1.0e-14_dp, &
                 "Exch-Ind22 must be Ind22 scaled by the second-order ratio")
      if (allocated(error)) return

      assembled = t2%total + t2%elst12 + t2%exch11 + t2%exch12 &
                  + t2%ind22 + t2%exch_ind22
      call check(error, abs(assembled - t2%total_sapt2) < 1.0e-14_dp, &
                 "the SAPT2 total is not the sum of its stated parts")

      call mols%destroy()
   end subroutine test_sapt2_total

end module test_mqc_sapt

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sapt, only: collect_mqc_sapt_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sapt", collect_mqc_sapt_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
