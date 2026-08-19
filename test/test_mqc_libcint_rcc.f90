module test_mqc_libcint_rcc
   !! Ties the spatial-orbital coupled cluster to the spin-orbital one.
   !!
   !! For a closed shell the two formulations are the same theory written twice,
   !! so they must agree to machine precision -- not to a tolerance, and not
   !! merely to the accuracy of some external reference. That identity is the
   !! whole reason `mqc_libcint_cc` was kept when this module was added, and it
   !! is a far sharper instrument than a PySCF number: a reference energy checks
   !! the total, while this checks the total, the singles/doubles split and the
   !! MP2 starting point, against an implementation that shares none of this
   !! one's index conventions.
   !!
   !! That last point is what makes it worth the runtime. The translation from
   !! PySCF's row-major einsum strings to column-major Fortran is where a
   !! spin-adapted implementation actually goes wrong, and a transposed index
   !! does not produce garbage -- it produces a converged, plausible, wrong
   !! number. The spin-orbital path indexes everything differently and cannot
   !! share such a mistake.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_mp2
   use mqc_libcint_cc, only: cc_result_t, run_libcint_ccsd
   use mqc_libcint_rcc, only: rcc_result_t, run_libcint_rccsd, ri_ladder_prefers_direct
   implicit none
   private
   public :: collect_mqc_libcint_rcc_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_libcint_rcc_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("spatial_mp2_equals_spin_orbital_mp2", test_mp2_identity), &
                  new_unittest("spatial_ccsd_equals_spin_orbital_ccsd", test_ccsd_identity), &
                  new_unittest("spatial_ccsd_matches_pyscf_cc_pvdz", test_pyscf_reference), &
                  new_unittest("fitted_path_agrees_with_spin_orbital", test_fitted_identity), &
                  new_unittest("direct_ri_ladder_agrees_with_assembled", test_ri_ladder_direct), &
                  new_unittest("frozen_core_agrees_between_formulations", test_frozen_identity), &
                  new_unittest("spatial_triples_equal_spin_orbital_triples", test_triples_identity) &
                  ]
   end subroutine collect_mqc_libcint_rcc_tests

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

   subroutine both_paths(basis, frozen, spin, spatial, ok, triples)
      !! Run the two formulations over one set of converged orbitals
      !!
      !! One SCF, not two: the amplitudes are being compared, so any difference
      !! in the reference would show up as a difference in them and would be
      !! indistinguishable from an error in the equations.
      character(len=*), intent(in) :: basis
      integer, intent(in) :: frozen
      type(cc_result_t), intent(out) :: spin
      type(rcc_result_t), intent(out) :: spatial
      logical, intent(out) :: ok
      logical, intent(in), optional :: triples

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      logical :: want_t

      want_t = .false.
      if (present(triples)) want_t = triples

      ok = .false.
      call converged_water(mol, scf, err, basis)
      if (err%has_error() .or. .not. scf%converged) return

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, frozen, &
                            60, 1.0e-10_dp, want_t, .false., spin, err)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if

      call run_libcint_rccsd(mol, scf%orbitals, scf%orbital_energies, 5, frozen, &
                             60, 1.0e-10_dp, want_t, .false., spatial, err)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if

      call mol%destroy()
      ok = spin%converged .and. spatial%converged
   end subroutine both_paths

   subroutine test_mp2_identity(error)
      !! The MP2 starting points must agree exactly
      !!
      !! Checked before CCSD because it isolates the half of the work that can
      !! fail on its own: the block layout, the chemists'-to-spatial index
      !! mapping and the frozen-core offset. If this passes and CCSD does not,
      !! the fault is in an amplitude equation and not in the integrals.
      type(error_type), allocatable, intent(out) :: error
      type(cc_result_t) :: spin
      type(rcc_result_t) :: spatial
      logical :: ok

      call both_paths("sto-3g", 0, spin, spatial, ok)
      call check(error, ok, "both coupled cluster paths must converge")
      if (allocated(error)) return

      call check(error, abs(spatial%e_mp2 - spin%e_mp2) < 1.0e-10_dp, &
                 "spatial and spin-orbital MP2 must agree")
   end subroutine test_mp2_identity

   subroutine test_ccsd_identity(error)
      !! The converged correlation energies must agree to machine precision
      !!
      !! And the singles and doubles parts separately. A total can agree while
      !! the two halves are individually wrong in opposite directions -- which
      !! is exactly what a misplaced t1 t1 product in tau does, since it moves
      !! weight between them without changing the sum.
      type(error_type), allocatable, intent(out) :: error
      type(cc_result_t) :: spin
      type(rcc_result_t) :: spatial
      logical :: ok

      call both_paths("sto-3g", 0, spin, spatial, ok)
      call check(error, ok, "both coupled cluster paths must converge")
      if (allocated(error)) return

      call check(error, abs(spatial%e_correlation - spin%e_singles - spin%e_doubles) &
                 < 1.0e-9_dp, "spatial and spin-orbital CCSD correlation must agree")
      if (allocated(error)) return

      call check(error, abs(spatial%e_doubles - spin%e_doubles) < 1.0e-9_dp, &
                 "the doubles parts must agree, not just the total")
      if (allocated(error)) return

      call check(error, abs(spatial%e_singles - spin%e_singles) < 1.0e-9_dp, &
                 "the singles parts must agree, not just the total")
      if (allocated(error)) return

      ! Sanity, in case both paths were wrong in the same direction: CCSD must
      ! lie below MP2 for a closed-shell molecule at its minimum.
      call check(error, spatial%e_correlation < spatial%e_mp2, &
                 "CCSD correlation must lie below MP2")
   end subroutine test_ccsd_identity

   subroutine test_pyscf_reference(error)
      !! cc-pVDZ against the PySCF number the spin-orbital path is validated on
      !!
      !! STO-3G has four virtual orbitals, which is few enough that several
      !! terms are nearly degenerate and a wrong contraction can hide. cc-pVDZ
      !! has nineteen and a d shell, so the particle-particle ladder does real
      !! work and an angular-momentum error in the transform would show.
      !!
      !! The reference is `validation/check_cc.f90`'s, quoted rather than
      !! recomputed: pyscf.cc.CCSD on this geometry, and the value the
      !! spin-orbital path already reproduces to 1e-8. Two implementations
      !! hitting one external number is a stronger statement than either hitting
      !! it alone.
      type(error_type), allocatable, intent(out) :: error
      type(cc_result_t) :: spin
      type(rcc_result_t) :: spatial
      logical :: ok

      real(dp), parameter :: E_CCSD_PYSCF = -0.213190638318_dp
      real(dp), parameter :: E_MP2_PYSCF = -0.203838686392_dp

      call both_paths("cc-pvdz", 0, spin, spatial, ok)
      call check(error, ok, "both paths must converge in cc-pVDZ")
      if (allocated(error)) return

      call check(error, abs(spatial%e_mp2 - E_MP2_PYSCF) < 1.0e-7_dp, &
                 "spatial MP2 must reproduce PySCF in cc-pVDZ")
      if (allocated(error)) return

      call check(error, abs(spatial%e_correlation - E_CCSD_PYSCF) < 1.0e-8_dp, &
                 "spatial CCSD must reproduce PySCF in cc-pVDZ")
      if (allocated(error)) return

      ! And still equal to the spin-orbital path, which is the tighter of the
      ! two statements -- an external reference is quoted to eight decimals and
      ! these two share a machine.
      call check(error, abs(spatial%e_correlation - spin%e_singles - spin%e_doubles) &
                 < 1.0e-9_dp, "the two formulations must agree in cc-pVDZ too")
   end subroutine test_pyscf_reference

   subroutine fitted_both_paths(basis, frozen, triples, spin, spatial, ok)
      !! Run both formulations density-fitted, over one set of exact orbitals
      !!
      !! The reference this is checked against is pyscf.cc.dfccsd on a
      !! *conventional* mean field, so the SCF here is exact too -- fitting it
      !! as well would mix a second approximation into a comparison meant to
      !! isolate the first.
      character(len=*), intent(in) :: basis
      integer, intent(in) :: frozen
      logical, intent(in) :: triples
      type(cc_result_t), intent(out) :: spin
      type(rcc_result_t), intent(out) :: spatial
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err

      ok = .false.
      call converged_water(mol, scf, err, basis)
      if (err%has_error() .or. .not. scf%converged) return

      call water(aux, err, "cc-pvdz-rifit")
      if (err%has_error()) then
         call mol%destroy()
         return
      end if

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, frozen, &
                            60, 1.0e-10_dp, triples, .false., spin, err, aux=aux)
      if (.not. err%has_error()) then
         call run_libcint_rccsd(mol, scf%orbitals, scf%orbital_energies, 5, frozen, &
                                60, 1.0e-10_dp, triples, .false., spatial, err, aux=aux)
      end if

      call mol%destroy()
      call aux%destroy()
      if (err%has_error()) return
      ok = spin%converged .and. spatial%converged
   end subroutine fitted_both_paths

   subroutine test_fitted_identity(error)
      !! The density-fitted spatial path, CCSD and (T), against every reference
      !!
      !! Its own case because the fitted route reaches `(ac|bd)` a completely
      !! different way -- a gemm over the auxiliary index inside the ladder,
      !! against a permuted read of a stored tensor -- and that is the one part
      !! of this module the conventional tests never exercise.
      !!
      !! The triples needed no fitted route of their own, which is what the
      !! second half of this asserts rather than assumes: the particle half
      !! wants (ia|bc) and the ring half (ij|ka), and `build_rcc_eris_fitted`
      !! already produces both. Nothing in `triples_correction` reads a block
      !! that only the conventional path builds.
      !!
      !! Compared against the fitted spin-orbital path, and against
      !! pyscf.cc.dfccsd -- never against conventional CCSD. The fitting error
      !! is 1.4e-4, three orders above the agreement asserted here, so the
      !! conventional answer is not the yardstick and cannot be.
      type(error_type), allocatable, intent(out) :: error
      type(cc_result_t) :: spin
      type(rcc_result_t) :: spatial
      logical :: ok

      ! validation/check_cc.f90's ri_case references, pyscf.cc.dfccsd.RCCSD
      real(dp), parameter :: E_CC_STO3G = -0.049173482895_dp
      real(dp), parameter :: E_T_STO3G = -0.000072751429_dp
      real(dp), parameter :: E_CC_PVDZ = -0.213326680500_dp
      real(dp), parameter :: E_T_PVDZ = -0.003041992009_dp
      real(dp), parameter :: E_CC_PVDZ_F1 = -0.211233849569_dp
      real(dp), parameter :: E_T_PVDZ_F1 = -0.003019775484_dp

      call fitted_both_paths("sto-3g", 0, .true., spin, spatial, ok)
      call check(error, ok, "both fitted paths must converge in STO-3G")
      if (allocated(error)) return
      call check(error, abs(spatial%e_correlation - spatial%e_triples - E_CC_STO3G) &
                 < 1.0e-9_dp, "fitted spatial CCSD must reproduce dfccsd, STO-3G")
      if (allocated(error)) return
      call check(error, abs(spatial%e_triples - E_T_STO3G) < 1.0e-9_dp, &
                 "fitted spatial (T) must reproduce dfccsd, STO-3G")
      if (allocated(error)) return
      call check(error, abs(spatial%e_triples - spin%e_triples) < 1.0e-10_dp, &
                 "the two fitted (T) corrections must agree, STO-3G")
      if (allocated(error)) return

      ! cc-pVDZ, where the ladder's auxiliary-index gemm does real work and a
      ! d shell is present.
      call fitted_both_paths("cc-pvdz", 0, .true., spin, spatial, ok)
      call check(error, ok, "both fitted paths must converge in cc-pVDZ")
      if (allocated(error)) return
      call check(error, abs(spatial%e_correlation - spatial%e_triples - E_CC_PVDZ) &
                 < 1.0e-8_dp, "fitted spatial CCSD must reproduce dfccsd, cc-pVDZ")
      if (allocated(error)) return
      call check(error, abs(spatial%e_triples - E_T_PVDZ) < 1.0e-8_dp, &
                 "fitted spatial (T) must reproduce dfccsd, cc-pVDZ")
      if (allocated(error)) return
      call check(error, abs(spatial%e_correlation - spin%e_singles - spin%e_doubles &
                            - spin%e_triples) < 1.0e-9_dp, &
                 "the two fitted formulations must agree, cc-pVDZ")
      if (allocated(error)) return

      ! And with a core frozen, which the triples index three times over.
      call fitted_both_paths("cc-pvdz", 1, .true., spin, spatial, ok)
      call check(error, ok, "both fitted paths must converge, frozen core")
      if (allocated(error)) return
      call check(error, abs(spatial%e_correlation - spatial%e_triples - E_CC_PVDZ_F1) &
                 < 1.0e-8_dp, "fitted spatial CCSD must reproduce dfccsd, frozen core")
      if (allocated(error)) return
      call check(error, abs(spatial%e_triples - E_T_PVDZ_F1) < 1.0e-8_dp, &
                 "fitted spatial (T) must reproduce dfccsd, frozen core")
      if (allocated(error)) return
      call check(error, abs(spatial%e_triples - spin%e_triples) < 1.0e-9_dp, &
                 "the two fitted (T) corrections must agree, frozen core")
   end subroutine test_fitted_identity

   subroutine test_ri_ladder_direct(error)
      !! The fitted ladder's other algorithm, on a system that selects it
      !!
      !! `ri_ladder_prefers_direct` chooses between assembling (ac|bd) and
      !! contracting the fitted tensor against tau without ever forming it, on
      !! measured operation counts. Every other fitted case in this file lands
      !! on the assembling side -- water has five occupied orbitals and the
      !! crossover needs roughly n_vir > 2 n_occ^2 -- so without this the whole
      !! direct path would be dead code as far as the suite is concerned.
      !!
      !! H2 in cc-pVDZ has one occupied orbital and nine virtual, which selects
      !! direct by a wide margin and costs almost nothing to run.
      !!
      !! The choice is asserted, not assumed. If the crossover is ever retuned
      !! so that this case stops selecting direct, this fails rather than
      !! quietly going back to testing the branch that was already covered.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(cc_result_t) :: spin
      type(rcc_result_t) :: spatial
      real(dp) :: c(3, 2)

      ! Representative counts for this case: one occupied, nine virtual, and an
      ! auxiliary basis of a few tens. The exact naux does not move the answer
      ! anywhere near the crossover.
      call check(error, ri_ladder_prefers_direct(1, 9, 50), &
                 "H2/cc-pVDZ must select the direct fitted ladder")
      if (allocated(error)) return
      call check(error,.not. ri_ladder_prefers_direct(5, 19, 84), &
                 "water/cc-pVDZ must select the assembled fitted ladder")
      if (allocated(error)) return

      c = reshape([0.0_dp, 0.0_dp, -0.3707_dp*ANG, &
                   0.0_dp, 0.0_dp, 0.3707_dp*ANG], [3, 2])
      call build_libcint_molecule([1, 1], ["H ", "H "], c, "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "H2 basis must load")
      if (allocated(error)) return

      call build_libcint_molecule([1, 1], ["H ", "H "], c, "cc-pvdz-rifit", aux, err)
      call check(error,.not. err%has_error(), "H2 auxiliary basis must load")
      if (allocated(error)) return

      call run_libcint_rhf(mol, 2, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                           in_core=.true.)
      call check(error,.not. err%has_error() .and. scf%converged, "H2 SCF must converge")
      if (allocated(error)) return

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 1, 0, &
                            60, 1.0e-10_dp, .true., .false., spin, err, aux=aux)
      call check(error,.not. err%has_error(), "fitted spin-orbital CCSD(T) must converge")
      if (allocated(error)) return

      call run_libcint_rccsd(mol, scf%orbitals, scf%orbital_energies, 1, 0, &
                             60, 1.0e-10_dp, .true., .false., spatial, err, aux=aux)
      call check(error,.not. err%has_error(), "fitted spatial CCSD(T) must converge")
      if (allocated(error)) return

      ! The spin-orbital path has only the one ladder, so this compares the
      ! direct algorithm against an assembled one -- which is the whole point.
      call check(error, abs(spatial%e_correlation - spin%e_singles - spin%e_doubles &
                            - spin%e_triples) < 1.0e-10_dp, &
                 "the direct fitted ladder must reproduce the assembled one")
      if (allocated(error)) return

      ! Two electrons in one orbital: CCSD is exact for a two-electron system,
      ! so the triples must be zero to numerical noise. A ladder that dropped
      ! or double-counted a term would not respect that.
      call check(error, abs(spatial%e_triples) < 1.0e-12_dp, &
                 "(T) must vanish for two electrons, where CCSD is already exact")

      call mol%destroy()
      call aux%destroy()
   end subroutine test_ri_ladder_direct

   subroutine test_frozen_identity(error)
      !! And they must still agree with a core orbital frozen
      !!
      !! Its own case because the frozen-core offset enters the two paths
      !! differently -- one indexes orbital energies through a spin-to-spatial
      !! map, the other slices the coefficient matrix -- so an off-by-one in
      !! either is invisible until the counts differ from zero.
      type(error_type), allocatable, intent(out) :: error
      type(cc_result_t) :: spin
      type(rcc_result_t) :: spatial
      logical :: ok

      real(dp), parameter :: E_T_PVDZ_F1 = -0.003015723855_dp   ! PySCF, check_cc.f90

      call both_paths("sto-3g", 1, spin, spatial, ok)
      call check(error, ok, "both paths must converge with a frozen core")
      if (allocated(error)) return

      call check(error, abs(spatial%e_mp2 - spin%e_mp2) < 1.0e-10_dp, &
                 "frozen-core MP2 must agree between formulations")
      if (allocated(error)) return

      call check(error, abs(spatial%e_correlation - spin%e_singles - spin%e_doubles) &
                 < 1.0e-9_dp, "frozen-core CCSD must agree between formulations")
      if (allocated(error)) return

      ! And the triples with a core frozen, in cc-pVDZ where the reference is
      ! quoted. The triples index the occupied space three times over, so a
      ! frozen-core offset that survives CCSD can still fail here.
      call both_paths("cc-pvdz", 1, spin, spatial, ok, triples=.true.)
      call check(error, ok, "both paths must converge, frozen core with triples")
      if (allocated(error)) return

      call check(error, abs(spatial%e_triples - spin%e_triples) < 1.0e-9_dp, &
                 "frozen-core (T) must agree between formulations")
      if (allocated(error)) return

      call check(error, abs(spatial%e_triples - E_T_PVDZ_F1) < 1.0e-8_dp, &
                 "frozen-core (T) must reproduce PySCF")
   end subroutine test_frozen_identity

   subroutine test_triples_identity(error)
      !! The (T) corrections must agree, and match PySCF
      !!
      !! Its own case rather than folded into the CCSD comparison because (T)
      !! is a separate algorithm reading the same amplitudes: it can be wrong
      !! while CCSD is right, and a converged CCSD is the precondition for it
      !! meaning anything at all.
      !!
      !! Checked in both bases. STO-3G is where a wrong permutation table still
      !! gives a plausible magnitude -- four virtual orbitals leave few distinct
      !! triples -- and cc-pVDZ is where it cannot hide.
      type(error_type), allocatable, intent(out) :: error
      type(cc_result_t) :: spin
      type(rcc_result_t) :: spatial
      logical :: ok

      real(dp), parameter :: E_T_STO3G = -0.000072742387_dp   ! PySCF, check_cc.f90
      real(dp), parameter :: E_T_PVDZ = -0.003037892387_dp

      call both_paths("sto-3g", 0, spin, spatial, ok, triples=.true.)
      call check(error, ok, "both paths must converge with triples")
      if (allocated(error)) return

      call check(error, abs(spatial%e_triples - spin%e_triples) < 1.0e-10_dp, &
                 "spatial and spin-orbital (T) must agree in STO-3G")
      if (allocated(error)) return

      call check(error, abs(spatial%e_triples - E_T_STO3G) < 1.0e-9_dp, &
                 "spatial (T) must reproduce PySCF in STO-3G")
      if (allocated(error)) return

      call both_paths("cc-pvdz", 0, spin, spatial, ok, triples=.true.)
      call check(error, ok, "both paths must converge with triples in cc-pVDZ")
      if (allocated(error)) return

      call check(error, abs(spatial%e_triples - spin%e_triples) < 1.0e-9_dp, &
                 "spatial and spin-orbital (T) must agree in cc-pVDZ")
      if (allocated(error)) return

      call check(error, abs(spatial%e_triples - E_T_PVDZ) < 1.0e-8_dp, &
                 "spatial (T) must reproduce PySCF in cc-pVDZ")
      if (allocated(error)) return

      ! (T) must lower the energy further for a closed-shell molecule at its
      ! minimum. A sign error in the permutation table is the failure this
      ! catches that the magnitude comparisons above might not.
      call check(error, spatial%e_triples < 0.0_dp, &
                 "(T) must be negative here")
   end subroutine test_triples_identity

end module test_mqc_libcint_rcc

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_rcc, only: collect_mqc_libcint_rcc_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_rcc", collect_mqc_libcint_rcc_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
