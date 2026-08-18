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
   use mqc_libcint_rcc, only: rcc_result_t, run_libcint_rccsd
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

   subroutine test_fitted_identity(error)
      !! The density-fitted spatial path against the fitted spin-orbital one
      !!
      !! Its own case because the fitted route reaches `(ac|bd)` a completely
      !! different way -- a gemm over the auxiliary index inside the ladder,
      !! against a permuted read of a stored tensor -- and that is the one part
      !! of this module the conventional tests never exercise.
      !!
      !! Compared against the fitted spin-orbital path rather than against the
      !! conventional spatial one: the fitting error is 1.4e-4, three orders
      !! above the agreement being asserted, so conventional CCSD is not the
      !! yardstick here and cannot be.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(cc_result_t) :: spin
      type(rcc_result_t) :: spatial

      call converged_water(mol, scf, err, "sto-3g")
      call check(error,.not. err%has_error() .and. scf%converged, "water SCF must converge")
      if (allocated(error)) return

      call water(aux, err, "cc-pvdz-rifit")
      call check(error,.not. err%has_error(), "the auxiliary basis must load")
      if (allocated(error)) return

      call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, 5, 0, &
                            60, 1.0e-10_dp, .false., .false., spin, err, aux=aux)
      call check(error,.not. err%has_error(), "fitted spin-orbital CCSD must converge")
      if (allocated(error)) return

      call run_libcint_rccsd(mol, scf%orbitals, scf%orbital_energies, 5, 0, &
                             60, 1.0e-10_dp, .false., .false., spatial, err, aux=aux)
      call check(error,.not. err%has_error(), "fitted spatial CCSD must converge")
      if (allocated(error)) return

      call check(error, abs(spatial%e_correlation - spin%e_singles - spin%e_doubles) &
                 < 1.0e-9_dp, "the two fitted formulations must agree")

      call mol%destroy()
      call aux%destroy()
   end subroutine test_fitted_identity

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
