module test_mqc_czt_neo
   !! NEO-HF and NEO-DFT: quantum protons against PySCF-NEO
   !!
   !! The reference is the `pyscf/neo` module of Yang Yang's PySCF fork, at
   !! commit f9c0266. The two codes share nothing: PySCF builds the cross
   !! Coulomb terms from its own two-basis integrals and runs one DIIS over
   !! all components, this code makes the cross terms out of a combined
   !! molecule and macro-iterates. Agreement to the macro-iteration's
   !! tolerance is therefore agreement on the physics -- the mass, the sign of
   !! every coupling, the ghosting of the quantum nucleus and the block
   !! bookkeeping in the combined basis.
   !!
   !! Every reference has two values because the electronic basis can be
   !! built in spherical or Cartesian form and the proton basis follows it;
   !! which one applies is read off the result rather than assumed.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_neo, only: neo_result_t, run_czt_neo_hf
   implicit none
   private
   public :: collect_mqc_czt_neo_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_czt_neo_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("hcn_with_a_quantum_proton_matches_pyscf_neo", test_hcn), &
                  new_unittest("two_quantum_protons_in_h2_match_pyscf_neo", test_h2), &
                  new_unittest("neo_dft_without_epc_matches_pyscf_neo", test_hcn_ks), &
                  new_unittest("epc17_2_matches_pyscf_neo", test_hcn_epc), &
                  new_unittest("epc_needs_a_functional", test_epc_needs_dft), &
                  new_unittest("an_unknown_nuclear_basis_is_refused", test_bad_nuclear_basis), &
                  new_unittest("an_unknown_epc_form_is_refused", test_bad_epc), &
                  new_unittest("a_heavy_atom_is_refused_for_now", test_refusal) &
                  ]
   end subroutine collect_mqc_czt_neo_tests

   subroutine hcn(result, err, quantum, basis, nuclear_basis, functional, epc)
      !! HCN as PySCF-NEO's test has it: H at the origin, C and N on z, in Angstrom
      !!
      !! Grid level 6 when a functional is named: PySCF-NEO's own level 3 and
      !! level 6 differ by 1e-6 on this molecule, so a comparison at the
      !! microhartree needs both codes on a fine grid.
      type(neo_result_t), intent(out) :: result
      type(error_t), intent(inout) :: err
      logical, intent(in) :: quantum(3)
      character(len=*), intent(in), optional :: basis, nuclear_basis, functional, epc
      character(len=:), allocatable :: use_basis, use_nuclear
      real(dp) :: c(3, 3)
      use_basis = "cc-pvdz"
      if (present(basis)) use_basis = basis
      use_nuclear = "pb4-d"
      if (present(nuclear_basis)) use_nuclear = nuclear_basis
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 1.064_dp*ANG, &
                   0.0_dp, 0.0_dp, 2.220_dp*ANG], [3, 3])
      call run_czt_neo_hf([1, 6, 7], ["H ", "C ", "N "], c, use_basis, use_nuclear, quantum, &
                          14, 200, 1.0e-10_dp, 1.0e-8_dp, .false., result, err, in_core=.true., &
                          functional=functional, grid_level=6, epc=epc)
   end subroutine hcn

   subroutine test_hcn(error)
      !! The energy PySCF-NEO gives for the same molecule, basis and proton basis
      !!
      !! `-92.8437063565785` is the fork's own test value in spherical harmonics,
      !! `-92.8442210525` the same run with `cart=True`, both with PB4-D. The
      !! proton density is checked to be one particle on the same run.
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err
      real(dp) :: reference, trace

      call hcn(result, err, [.true., .false., .false.])
      call check(error,.not. err%has_error(), "NEO-HF failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, result%converged, "the macro-iteration did not converge")
      if (allocated(error)) return
      if (result%cartesian) then
         reference = -92.8442210525_dp
      else
         reference = -92.8437063565785_dp
      end if
      call check(error, abs(result%energy - reference) < 2.0e-6_dp, &
                 "the NEO-HF energy disagrees with PySCF-NEO")
      if (allocated(error)) return
      ! The proton sits in a well: its orbital energy is negative, and the
      ! electron-proton attraction is what holds it there.
      call check(error, result%nuclear_orbital_energies(1, 1) < 0.0_dp, &
                 "the proton orbital energy is not negative")
      if (allocated(error)) return
      call check(error, result%electron_nucleus < 0.0_dp, &
                 "the electron-proton energy is not attractive")
      if (allocated(error)) return
      trace = sum(result%nuclear_densities(:, :, 1)*result%nuclear_overlaps(:, :, 1))
      call check(error, abs(trace - 1.0_dp) < 1.0e-10_dp, "the proton density is not one particle")
   end subroutine test_hcn

   subroutine test_h2(error)
      !! Two quantum protons at once: H2 with both nuclei quantised
      !!
      !! The only case that exercises the other-proton terms -- the `q /= p`
      !! field in a proton's Fock matrix, the proton-proton repulsion in the
      !! energy, and the second proton block of the combined basis. From
      !! PySCF-NEO, H at 0 and 0.74 Angstrom on z, cc-pVDZ, PB4-D:
      !! `-1.0507815470` spherical, `-1.0508155378` Cartesian. The two protons
      !! are equivalent by symmetry, so their orbital energies must agree.
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err
      real(dp) :: c(3, 2), reference, trace
      integer :: p

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.74_dp*ANG], [3, 2])
      call run_czt_neo_hf([1, 1], ["H ", "H "], c, "cc-pvdz", "pb4-d", [.true., .true.], &
                          2, 200, 1.0e-10_dp, 1.0e-8_dp, .false., result, err, in_core=.true.)
      call check(error,.not. err%has_error(), "NEO-HF on H2 failed: "//err%get_message())
      if (allocated(error)) return
      if (result%cartesian) then
         reference = -1.0508155378_dp
      else
         reference = -1.0507815470_dp
      end if
      call check(error, abs(result%energy - reference) < 2.0e-6_dp, &
                 "the two-proton NEO-HF energy disagrees with PySCF-NEO")
      if (allocated(error)) return
      call check(error, result%n_quantum == 2, "two protons were asked for")
      if (allocated(error)) return
      call check(error, abs(result%nuclear_orbital_energies(1, 1) &
                            - result%nuclear_orbital_energies(1, 2)) < 1.0e-6_dp, &
                 "the two equivalent protons have different orbital energies")
      if (allocated(error)) return
      call check(error, result%nucleus_nucleus > 0.0_dp, &
                 "the proton-proton repulsion is not positive")
      if (allocated(error)) return
      do p = 1, 2
         trace = sum(result%nuclear_densities(:, :, p)*result%nuclear_overlaps(:, :, p))
         call check(error, abs(trace - 1.0_dp) < 1.0e-10_dp, "a proton density is not one particle")
         if (allocated(error)) return
      end do
   end subroutine test_h2

   subroutine test_hcn_ks(error)
      !! NEO-DFT with B3LYP5 electrons and no electron-proton functional
      !!
      !! Pins the Kohn-Sham coupling on its own, which is what caught the XC
      !! grid being built for "element 0" on the ghosted proton. In 6-31G to
      !! keep it quick: `-93.3001036793` (spherical) and `-93.3001188171`
      !! (Cartesian) from PySCF-NEO, `xc='HYB_GGA_XC_B3LYP5'`, grid level 6,
      !! which is the functional PySCF's `b3lyp5` names.
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err
      real(dp) :: reference

      call hcn(result, err, [.true., .false., .false.], basis="6-31g", &
               functional="hyb_gga_xc_b3lyp5")
      call check(error,.not. err%has_error(), "NEO-DFT failed: "//err%get_message())
      if (allocated(error)) return
      if (result%cartesian) then
         reference = -93.3001188171_dp
      else
         reference = -93.3001036793_dp
      end if
      call check(error, abs(result%energy - reference) < 2.0e-6_dp, &
                 "the NEO-DFT energy disagrees with PySCF-NEO")
      if (allocated(error)) return
      call check(error, result%kohn_sham, "the electrons did not run as Kohn-Sham")
   end subroutine test_hcn_ks

   subroutine test_hcn_epc(error)
      !! epc17-2 on top of B3LYP5: `-93.3670509407` (spherical) and
      !! `-93.3694593587` (Cartesian) from PySCF-NEO at grid level 6, cc-pVDZ
      !!
      !! The functional is integrated on the electronic grid where the proton
      !! density lives, its electronic potential rides frozen through each
      !! macro-iteration, and its energy replaces the frozen `Tr(D V)` the SCF
      !! counted. Each of those is a place to be wrong by a few millihartree,
      !! which is also the size of the whole correction here.
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err
      real(dp) :: reference

      call hcn(result, err, [.true., .false., .false.], functional="hyb_gga_xc_b3lyp5", &
               epc="17-2")
      call check(error,.not. err%has_error(), "NEO-DFT with epc failed: "//err%get_message())
      if (allocated(error)) return
      if (result%cartesian) then
         reference = -93.3694593587_dp
      else
         reference = -93.3670509407_dp
      end if
      call check(error, abs(result%energy - reference) < 2.0e-6_dp, &
                 "the epc17-2 energy disagrees with PySCF-NEO")
      if (allocated(error)) return
      call check(error, result%epc_energy < -0.01_dp .and. result%epc_energy > -0.1_dp, &
                 "the electron-proton correlation energy is not a few tens of millihartree")
   end subroutine test_hcn_epc

   subroutine test_epc_needs_dft(error)
      !! A correlation functional on a Hartree-Fock electron is refused
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err

      call hcn(result, err, [.true., .false., .false.], epc="17-2")
      call check(error, err%has_error(), "epc on Hartree-Fock electrons was accepted")
   end subroutine test_epc_needs_dft

   subroutine test_bad_nuclear_basis(error)
      !! A proton basis that does not exist is refused by name, not by a crash
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err

      call hcn(result, err, [.true., .false., .false.], nuclear_basis="pb9-z")
      call check(error, err%has_error(), "an unknown nuclear basis was accepted")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "pb9-z") > 0, &
                 "the refusal does not name the basis")
   end subroutine test_bad_nuclear_basis

   subroutine test_bad_epc(error)
      !! An epc form this code does not know is refused, naming it
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err

      call hcn(result, err, [.true., .false., .false.], functional="hyb_gga_xc_b3lyp5", &
               epc="19")
      call check(error, err%has_error(), "an unknown epc form was accepted")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "19") > 0, "the refusal does not name the form")
   end subroutine test_bad_epc

   subroutine test_refusal(error)
      !! Quantising carbon is not implemented yet and must say so, not run
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err

      call hcn(result, err, [.true., .true., .false.])
      call check(error, err%has_error(), "a quantum carbon was accepted")
   end subroutine test_refusal

end module test_mqc_czt_neo

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_neo, only: collect_mqc_czt_neo_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_neo", collect_mqc_czt_neo_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
