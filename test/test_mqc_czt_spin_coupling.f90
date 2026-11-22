!! Indirect nuclear spin-spin couplings against PySCF
module test_mqc_czt_spin_coupling
   !! The Fermi-contact reduced couplings of water in 6-31G, against
   !! `pyscf.prop.ssc`'s `make_fc` fed this repository's basis set: at
   !! Hartree-Fock, where the triplet response is exchange alone, and at
   !! B3LYP, where it also carries the functional's spin kernel. The numbers
   !! are the reduced couplings in atomic units, which no isotope enters.
   !!
   !! Beside them, what needs no other program: the matrix is symmetric with
   !! a zero diagonal, the conversion to hertz reproduces a value worked out
   !! by hand, and the isotope table gives carbon as 13C.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_czt_spin_coupling, only: czt_fermi_contact
   use mqc_nmr_isotopes, only: nmr_isotope, reduced_to_hertz, reduced_to_si
   implicit none
   private

   public :: collect_mqc_czt_spin_coupling_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   ! PySCF, water/6-31G: O-H(1), O-H(2), H-H. Atomic units.
   real(dp), parameter :: K_HF(3) = [3.587784200978e-8_dp, 3.594115077811e-8_dp, &
                                     -1.579497529373e-9_dp]
   real(dp), parameter :: K_B3LYP(3) = [2.113607620152e-8_dp, 2.124488044674e-8_dp, &
                                        -4.594923671731e-10_dp]
   real(dp), parameter :: REL_TOL = 1.0e-5_dp
      !! 4e-8 on the large couplings and 1.5e-6 on the small one, measured;
      !! the small one is where the two response solves' residuals show.

contains

   subroutine collect_mqc_czt_spin_coupling_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_hartree_fock_matches_pyscf", test_hf), &
                  new_unittest("water_b3lyp_matches_pyscf", test_b3lyp), &
                  new_unittest("isotopes_and_units", test_units) &
                  ]
   end subroutine collect_mqc_czt_spin_coupling_tests

   subroutine water(mol, err)
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
   end subroutine water

   subroutine compare(error, k, reference, label)
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: k(:, :), reference(3)
      character(len=*), intent(in) :: label
      real(dp) :: got(3)

      got = [k(1, 2), k(1, 3), k(2, 3)]
      write (*, "(4x,a,3es18.10)") label//"  ", got
      write (*, "(4x,a,3es18.10)") "PySCF           ", reference
      write (*, "(4x,a,es10.2)") "worst relative difference ", &
         maxval(abs(got - reference)/abs(reference))
      call check(error, maxval(abs(k - transpose(k))) < 1.0e-20_dp, &
                 "the coupling matrix is not symmetric")
      if (allocated(error)) return
      call check(error, abs(k(1, 1)) + abs(k(2, 2)) + abs(k(3, 3)) < tiny(1.0_dp), &
                 "a nucleus is coupled to itself")
      if (allocated(error)) return
      call check(error, maxval(abs(got - reference)/abs(reference)) < REL_TOL, &
                 label//": the couplings do not match PySCF's")
   end subroutine compare

   subroutine test_hf(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: k(:, :)

      call water(mol, err)
      call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "water must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call czt_fermi_contact(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                             k, err)
      call check(error,.not. err%has_error(), "couplings: "//err%get_full_trace())
      if (allocated(error)) return
      call compare(error, k, K_HF, "Hartree-Fock    ")
   end subroutine test_hf

   subroutine test_b3lyp(error)
      !! The same with the functional's spin kernel in the response. A
      !! triplet solve against the singlet kernel is off by a factor near
      !! two on the H-H coupling, so this is what says the triplet one ran.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc
      type(error_t) :: err
      real(dp), allocatable :: k(:, :)

      if (.not. xc_available()) return
      call water(mol, err)
      if (.not. err%has_error()) call xc_context_create(mol, "b3lyp", xc, err)
      if (.not. err%has_error()) &
         call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-10_dp, .false., scf, err, xc=xc)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "water must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call czt_fermi_contact(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                             k, err, xc=xc, density=scf%density)
      call check(error,.not. err%has_error(), "couplings: "//err%get_full_trace())
      if (allocated(error)) return
      call compare(error, k, K_B3LYP, "B3LYP           ")
   end subroutine test_b3lyp

   subroutine test_units(error)
      !! The table and the two conversions, against numbers worked by hand.
      type(error_type), allocatable, intent(out) :: error

      character(len=8) :: label
      real(dp) :: g_h, g_c, g_o, g_none, j, k_si

      call nmr_isotope(1, label, g_h)
      call check(error, trim(label) == "1H" .and. abs(g_h - 5.5856946893_dp) < 1.0e-9_dp, &
                 "hydrogen is not 1H")
      if (allocated(error)) return
      ! The abundant carbon has no spin; the one observed is 13C.
      call nmr_isotope(6, label, g_c)
      call check(error, trim(label) == "13C" .and. abs(g_c - 1.4048236_dp) < 1.0e-7_dp, &
                 "carbon is not 13C")
      if (allocated(error)) return
      call nmr_isotope(8, label, g_o)
      call nmr_isotope(26, label, g_none)
      call check(error, len_trim(label) == 0 .and. abs(g_none) < tiny(1.0_dp), &
                 "an element beyond the table has an isotope")
      if (allocated(error)) return

      ! J = (E_h / h) (m_e / 2 m_p)^2 g_A g_B K
      !   = 6.5796839e15 * 7.4151505e-8 * (-0.757516 * 5.5856947) * 3.5877842e-8
      !   = -74.067 Hz, water's one-bond 17O-1H at Hartree-Fock.
      j = reduced_to_hertz(K_HF(1), g_o, g_h)
      call check(error, abs(j - (-74.067_dp)) < 2.0e-3_dp, "J(17O-1H) in hertz is wrong")
      if (allocated(error)) return
      ! K = E_h (m_e / 2 m_p)^2 / mu_N^2 * K_au = 1.26727e28 * 3.5877842e-8 T^2/J
      k_si = reduced_to_si(K_HF(1))
      call check(error, abs(k_si - 45.467_dp) < 2.0e-3_dp, &
                 "the reduced coupling in 1e19 T^2/J is wrong")
   end subroutine test_units

end module test_mqc_czt_spin_coupling

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_spin_coupling, only: collect_mqc_czt_spin_coupling_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_spin_coupling", collect_mqc_czt_spin_coupling_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
