!! Manual check that the CPU coupled cluster reproduces PySCF
!!
!!     cmake -B build -DMQC_ENABLE_CZT=ON && ./build/check_cc
!!
!! Runs RHF, then CCSD and (T), on a ladder of systems chosen so a failure
!! localises. The milestone order of CC_PLAN.md, in one program:
!!
!!   * **MP2 out of the spin-orbital integrals must equal `run_czt_mp2`.**
!!     Checked first and to 1e-10, because it is free -- MP2 is the first CCSD
!!     iteration's energy -- and because it isolates the three things most
!!     likely to be wrong (the antisymmetrisation, the spin-orbital index map,
!!     the denominators) with no amplitude equations in the way. If this
!!     disagrees, nothing after it can be right.
!!   * **CCSD against PySCF `cc.CCSD`**, and **(T) against `.ccsd_t()`**.
!!
!! References are PySCF 2.14 fed this repository's own basis JSON through
!! `bse_to_pyscf`, not PySCF's internal tables: on Pople sets those differ in
!! the eighth decimal, which is enough to look exactly like a bug here. See
!! tools/cpu_validation/gen_cpu_validation.py.
program check_cc
   use pic_types, only: dp
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mp2, only: mp2_result_t, run_czt_mp2
   use mqc_czt_cc, only: cc_result_t, run_czt_ccsd
   use mqc_error, only: error_t
   implicit none

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer :: failures

   failures = 0

   ! H2O / STO-3G, no frozen core. Small enough to debug by hand: seven spatial
   ! orbitals, fourteen spin orbitals, four virtual.
   call one_case("H2O sto-3g    ", "sto-3g", 0, &
                 -74.962005712275_dp, &      ! PySCF RHF
                 -0.035266770320_dp, &       ! PySCF MP2 correlation
                 -0.049173927119_dp, &       ! PySCF CCSD correlation
                 -0.000072742387_dp)         ! PySCF (T)

   ! The same molecule in a basis with a d shell, which is where an angular
   ! momentum error in the transform would first show.
   call one_case("H2O cc-pvdz   ", "cc-pvdz", 0, &
                 -76.026435906898_dp, &
                 -0.203838686392_dp, &
                 -0.213190638318_dp, &
                 -0.003037892387_dp)

   ! The frozen-core path, on the same molecule so only the freezing changes.
   call one_case("H2O cc-pvdz f1", "cc-pvdz", 1, &
                 -76.026435906898_dp, &
                 -0.201503456922_dp, &
                 -0.211097952217_dp, &
                 -0.003015723855_dp)

   ! RI. References are pyscf.cc.dfccsd.RCCSD on a *conventional* mean field, so
   ! the comparison isolates the CC-side fitting rather than mixing in a fitted
   ! Hartree-Fock. The fitting error against conventional CCSD is 1.4e-4 -- three
   ! orders above the 4e-10 the exact path agrees to -- so conventional CCSD is
   ! not the yardstick here and cannot be.
   call ri_case("RI H2O sto-3g    ", "sto-3g", "cc-pvdz-rifit", 0, &
                -0.049173482895_dp, -0.000072751429_dp)
   call ri_case("RI H2O sto-3g  f1", "sto-3g", "cc-pvdz-rifit", 1, &
                -0.049095167549_dp, -0.000072820877_dp)
   call ri_case("RI H2O cc-pvdz   ", "cc-pvdz", "cc-pvdz-rifit", 0, &
                -0.213326680500_dp, -0.003041992009_dp)
   call ri_case("RI H2O cc-pvdz f1", "cc-pvdz", "cc-pvdz-rifit", 1, &
                -0.211233849569_dp, -0.003019775484_dp)

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[cc] all ok -- CCSD(T) reproduces PySCF, fitted and exact"
   else
      write (*, "(A,I0,A)") "[cc] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine one_case(label, basis, frozen, e_scf_ref, e_mp2_ref, e_ccsd_ref, e_t_ref)
      character(len=*), intent(in) :: label, basis
      integer, intent(in) :: frozen
      real(dp), intent(in) :: e_scf_ref, e_mp2_ref, e_ccsd_ref, e_t_ref

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2
      type(cc_result_t) :: cc
      type(error_t) :: err
      real(dp) :: c(3, 3)
      real(dp) :: e_ccsd
      integer :: n_occ

      ! The geometry every CPU water case in the validation suite uses.
      c = reshape([0.0_dp, 0.00000000009155_dp*ANG, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[cc] ", label, " basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      n_occ = 5
      call run_czt_rhf(mol, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                       in_core=.true.)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A,A,A)") "[cc] ", label, " SCF failed"
         failures = failures + 1
         return
      end if
      call report(label//" SCF ", scf%energy, e_scf_ref, 1.0e-8_dp)

      ! Conventional MP2 through the existing path, as the thing the
      ! spin-orbital MP2 has to reproduce exactly.
      call run_czt_mp2(mol, scf%orbitals, scf%orbital_energies, n_occ, &
                       scf%energy, mp2, err, n_frozen=frozen)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[cc] ", label, " MP2 failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call run_czt_ccsd(mol, scf%orbitals, scf%orbital_energies, n_occ, frozen, &
                        100, 1.0e-11_dp, .true., .false., cc, err)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[cc] ", label, " CCSD failed: ", err%get_message()
         failures = failures + 1
         call mol%destroy()
         return
      end if

      ! The free checkpoint: same MP2, two entirely different routes to it.
      call report(label//" MP2=", cc%e_mp2, mp2%correlation, 1.0e-10_dp)
      call report(label//" MP2 ", cc%e_mp2, e_mp2_ref, 1.0e-7_dp)
      e_ccsd = cc%e_singles + cc%e_doubles
      call report(label//" CCSD", e_ccsd, e_ccsd_ref, 1.0e-8_dp)
      call report(label//" (T) ", cc%e_triples, e_t_ref, 1.0e-8_dp)

      call mol%destroy()
   end subroutine one_case

   subroutine ri_case(label, basis, auxbasis, frozen, e_ccsd_ref, e_t_ref)
      !! RI-CCSD and RI-CCSD(T) against pyscf.cc.dfccsd.RCCSD
      character(len=*), intent(in) :: label, basis, auxbasis
      integer, intent(in) :: frozen
      real(dp), intent(in) :: e_ccsd_ref, e_t_ref

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(cc_result_t) :: cc
      type(error_t) :: err
      real(dp) :: c(3, 3)
      real(dp) :: e_ccsd

      c = reshape([0.0_dp, 0.00000000009155_dp*ANG, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[cc] ", label, " basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, auxbasis, aux, err)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[cc] ", label, " aux basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      ! Exact SCF on purpose, matching how the reference was produced.
      call run_czt_rhf(mol, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                       in_core=.true.)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A,A,A)") "[cc] ", label, " SCF failed"
         failures = failures + 1
         return
      end if

      call run_czt_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, frozen, &
                        100, 1.0e-11_dp, .true., .false., cc, err, aux=aux)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[cc] ", label, " RI-CCSD failed: ", err%get_message()
         failures = failures + 1
         call mol%destroy()
         call aux%destroy()
         return
      end if

      e_ccsd = cc%e_singles + cc%e_doubles
      call report(label//" CCSD", e_ccsd, e_ccsd_ref, 1.0e-8_dp)
      call report(label//" (T) ", cc%e_triples, e_t_ref, 1.0e-8_dp)

      call mol%destroy()
      call aux%destroy()
   end subroutine ri_case

   subroutine report(what, got, want, tol)
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: got, want, tol

      if (abs(got - want) <= tol) then
         write (*, "(A,A,A,F18.10,A,ES10.3)") "  ok   ", what, " ", got, "  diff ", &
            abs(got - want)
      else
         write (*, "(A,A,A,F18.10,A,F18.10,A,ES10.3)") &
            "  FAIL ", what, " got ", got, " want ", want, " diff ", abs(got - want)
         failures = failures + 1
      end if
   end subroutine report

end program check_cc
