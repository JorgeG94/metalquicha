module test_mqc_mbis
   !! Pins MBIS charges and spin populations against denspart.
   !!
   !! The reference numbers come from `tools/mbis/mbis_reference.py`, which runs
   !! denspart -- the MBIS code of the method's author -- on a PySCF density
   !! integrated on PySCF's own 150x974 grid. denspart minimises the extended
   !! Kullback-Leibler divergence from its own table of starting parameters;
   !! this code iterates the MBIS fixed point from shells filled in aufbau
   !! order with Slater-screened exponents. The two meet at the same pro-model,
   !! so agreement checks the iteration and its start as well as the
   !! quadrature. What is checked:
   !!
   !!   * water's charges and fitted shells, RHF/cc-pVDZ;
   !!   * its cation's charges and spin populations, UHF;
   !!   * hydrogen sulfide, whose sulfur carries three shells;
   !!   * lithium fluoride, ionic, where the start decides the answer: an
   !!     earlier start collapsed lithium's valence shell onto its core;
   !!   * a lone Li+, which needs no fit at all;
   !!   * the charges sum to the molecular charge to round-off, renormalised
   !!     to Tr(D S), including an anion in aug-cc-pVDZ whose diffuse tail the
   !!     grid misses 5.8e-6 electrons of;
   !!   * the shell count per row, and the refusal of an effective core.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf
   use mqc_czt_mbis, only: mbis_model_t, mbis_populations, mbis_shell_count
   implicit none
   private
   public :: collect_mqc_mbis_tests

   character(len=2), parameter :: SYM(3) = ["X ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])
      !! Bohr, the water of test_mqc_fukui and of the reference script
   real(dp), parameter :: H2S(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.8187_dp, 1.7419_dp, &
                           0.0_dp, -1.8187_dp, 1.7419_dp], [3, 3])

   real(dp), parameter :: Q_WATER(3) = [-0.861711465820_dp, 0.430855733573_dp, &
                                        0.430855733198_dp]
   real(dp), parameter :: SHELLS_WATER(2, 4) = reshape( &
                          [1.6062306175_dp, 17.9535165069_dp, &
                           7.2554808483_dp, 2.4898333644_dp, &
                           0.5691442664_dp, 2.8751690013_dp, &
                           0.5691442668_dp, 2.8751689988_dp], [2, 4])
      !! (population, exponent) per shell: oxygen's two, then each hydrogen's
   real(dp), parameter :: Q_CATION(3) = [-0.067655869303_dp, 0.533827928566_dp, &
                                         0.533827931540_dp]
   real(dp), parameter :: S_CATION(3) = [1.008953830950_dp, -0.004476915392_dp, &
                                         -0.004476915432_dp]
   real(dp), parameter :: Q_H2S(3) = [-0.333305106574_dp, 0.166652546902_dp, &
                                      0.166652547024_dp]

   real(dp), parameter :: Q_LIF(2) = [0.940409209270_dp, -0.940409390994_dp]
   real(dp), parameter :: LIF_TOL = 1.0e-6_dp
      !! denspart's own LiF is converged no further than this: its two charges
      !! differ in magnitude by 1.8e-7.
   real(dp), parameter :: CHARGE_TOL = 1.0e-7_dp
      !! Measured at 2e-8 or better on charges, spin populations and shell
      !! populations, from grid level 4 up, and exponents to 6e-9 relative.
      !! Unchanged by the DIIS that cut water's fit from 72 iterations to 19.
   real(dp), parameter :: SUM_TOL = 1.0e-10_dp
      !! Round-off: the populations are renormalised to Tr(D S). Before, the
      !! grid's own error showed here, 1.4e-8 at the default level.

contains

   subroutine collect_mqc_mbis_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_matches_denspart", test_water), &
                  new_unittest("cation_spins_match_denspart", test_cation), &
                  new_unittest("h2s_matches_denspart", test_h2s), &
                  new_unittest("shells_per_row", test_shell_count), &
                  new_unittest("effective_core_is_refused", test_ecp_refused), &
                  new_unittest("bare_cation_needs_no_fit", test_bare_cation), &
                  new_unittest("ionic_lif_matches_denspart", test_lif), &
                  new_unittest("diffuse_anion_sums_exactly", test_diffuse_anion) &
                  ]
   end subroutine collect_mqc_mbis_tests

   subroutine converge(z, coords, charge, multiplicity, mol, scf, err)
      integer, intent(in) :: z(3)
      real(dp), intent(in) :: coords(3, 3)
      integer, intent(in) :: charge, multiplicity
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      character(len=2) :: symbols(3)

      symbols = SYM
      symbols(1) = merge("O ", "S ", z(1) == 8)
      call build_czt_molecule(z, symbols, coords, "cc-pvdz", mol, err)
      if (err%has_error()) return
      if (multiplicity == 1) then
         call run_czt_rhf(mol, sum(z) - charge, 100, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err)
      else
         call run_czt_uhf(mol, sum(z) - charge, multiplicity, 100, 1.0e-11_dp, 1.0e-9_dp, &
                          .false., scf, err)
      end if
   end subroutine converge

   subroutine test_water(error)
      !! Charges, and every fitted shell
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(mbis_model_t) :: model
      real(dp), allocatable :: q(:)
      integer :: i

      call converge([8, 1, 1], WATER, 0, 1, mol, scf, err)
      call mbis_populations(mol, scf%density, q, err, model=model)
      call check(error,.not. err%has_error(), "water: "//err%get_full_trace())
      if (allocated(error)) return

      do i = 1, 3
         call check(error, q(i), Q_WATER(i), "water charge", thr=CHARGE_TOL)
         if (allocated(error)) return
      end do
      call check(error, model%n_shells, 4, "oxygen two shells, each hydrogen one")
      if (allocated(error)) return
      call check(error, model%retired + model%merged, 0, &
                 "water reduces nothing, so this is denspart's model")
      if (allocated(error)) return
      do i = 1, 4
         call check(error, model%population(i), SHELLS_WATER(1, i), "shell population", &
                    thr=CHARGE_TOL)
         if (allocated(error)) return
         call check(error, model%exponent(i), SHELLS_WATER(2, i), "shell exponent", &
                    thr=CHARGE_TOL*SHELLS_WATER(2, i))
         if (allocated(error)) return
      end do
      call check(error, abs(sum(q)) < SUM_TOL, "water's charges must sum to zero")
      call mol%destroy()
   end subroutine test_water

   subroutine test_cation(error)
      !! Open-shell charges, and spin populations from the converged weights
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: q(:), s(:)
      integer :: i

      call converge([8, 1, 1], WATER, 1, 2, mol, scf, err)
      call mbis_populations(mol, scf%density + scf%density_beta, q, err, &
                            spin_density=scf%density - scf%density_beta, spin_populations=s)
      call check(error,.not. err%has_error(), "cation: "//err%get_full_trace())
      if (allocated(error)) return

      do i = 1, 3
         call check(error, q(i), Q_CATION(i), "cation charge", thr=CHARGE_TOL)
         if (allocated(error)) return
         call check(error, s(i), S_CATION(i), "cation spin population", thr=CHARGE_TOL)
         if (allocated(error)) return
      end do
      call check(error, abs(sum(q) - 1.0_dp) < SUM_TOL, "the cation's charges must sum to one")
      if (allocated(error)) return
      call check(error, abs(sum(s) - 1.0_dp) < SUM_TOL, "one unpaired electron")
      call mol%destroy()
   end subroutine test_cation

   subroutine test_h2s(error)
      !! A third-row atom, fitted with three shells
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(mbis_model_t) :: model
      real(dp), allocatable :: q(:)
      integer :: i

      call converge([16, 1, 1], H2S, 0, 1, mol, scf, err)
      call mbis_populations(mol, scf%density, q, err, model=model)
      call check(error,.not. err%has_error(), "h2s: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, model%n_shells, 5, "sulfur three shells, each hydrogen one")
      if (allocated(error)) return
      do i = 1, 3
         call check(error, q(i), Q_H2S(i), "h2s charge", thr=CHARGE_TOL)
         if (allocated(error)) return
      end do
      call mol%destroy()
   end subroutine test_h2s

   subroutine test_shell_count(error)
      !! One shell per row reached
      type(error_type), allocatable, intent(out) :: error

      call check(error, all([mbis_shell_count(1), mbis_shell_count(2), mbis_shell_count(3), &
                             mbis_shell_count(10), mbis_shell_count(11), &
                             mbis_shell_count(36), mbis_shell_count(37), &
                             mbis_shell_count(118)] == [1, 1, 2, 2, 3, 4, 5, 7]), &
                 "shells per row")
   end subroutine test_shell_count

   subroutine test_ecp_refused(error)
      !! An effective core removes electrons the shell count does not know of
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: q(:), d(:, :)
      real(dp) :: c(3, 2)

      c = 0.0_dp
      c(3, 2) = 3.0_dp
      call build_czt_molecule([1, 53], ["H ", "I "], c, "def2-svp", mol, err, &
                              ecp_name="def2-ecp")
      call check(error,.not. err%has_error(), "HI must build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, any(mol%core_electrons > 0), "iodine must carry a core")
      if (allocated(error)) return
      allocate (d(mol%nao, mol%nao))
      d = 0.0_dp
      call mbis_populations(mol, d, q, err)
      call check(error, err%has_error(), "an effective core must be refused")
      call mol%destroy()
   end subroutine test_ecp_refused

   subroutine test_lif(error)
      !! An ionic molecule, whose lithium keeps a valence shell of 0.08 electrons
      !!
      !! The case the lone-atom shortcut does not cover: two atoms, so a real
      !! fit, with one shell close to empty but not empty enough to retire.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(mbis_model_t) :: model
      real(dp), allocatable :: q(:)
      real(dp) :: c(3, 2)
      integer :: i

      c = 0.0_dp
      c(3, 2) = 2.955_dp
      call build_czt_molecule([3, 9], ["Li", "F "], c, "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "LiF must build: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_rhf(mol, 12, 100, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err)
      call mbis_populations(mol, scf%density, q, err, model=model)
      call check(error,.not. err%has_error(), "LiF: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, model%n_shells, 4, "two shells each, none retired or merged")
      if (allocated(error)) return
      do i = 1, 2
         call check(error, q(i), Q_LIF(i), "LiF charge", thr=LIF_TOL)
         if (allocated(error)) return
      end do
      call mol%destroy()
   end subroutine test_lif

   subroutine test_diffuse_anion(error)
      !! Water's anion in aug-cc-pVDZ, whose diffuse tail the grid misses
      !!
      !! 5.8e-6 electrons short at the default level. Renormalised to Tr(D S),
      !! the charges sum to -1 and the spins to one unpaired electron anyway.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: q(:), s(:)

      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], WATER, "aug-cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "water must build: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_uhf(mol, 11, 2, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err)
      call mbis_populations(mol, scf%density + scf%density_beta, q, err, &
                            spin_density=scf%density - scf%density_beta, spin_populations=s)
      call check(error,.not. err%has_error(), "anion: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, sum(q), -1.0_dp, "the anion's charges sum to -1", thr=SUM_TOL)
      if (allocated(error)) return
      call check(error, sum(s), 1.0_dp, "and its spins to one unpaired electron", thr=SUM_TOL)
      call mol%destroy()
   end subroutine test_diffuse_anion

   subroutine test_bare_cation(error)
      !! A lone Li+ is charged +1 without a fit
      !!
      !! With one atom every weight is one, so no pro-model can change the
      !! answer -- and the fit would be the flattest there is, two shells for
      !! one 1s pair, which crawled along the split between them to the
      !! 1000-iteration cap before lone atoms skipped it. This is what a bare
      !! cation monomer of an ion-water cluster looks like.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(mbis_model_t) :: model
      real(dp), allocatable :: q(:)
      real(dp) :: c(3, 1)

      c = 0.0_dp
      call build_czt_molecule([3], ["Li"], c, "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "Li+ must build: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_rhf(mol, 2, 100, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err)
      call mbis_populations(mol, scf%density, q, err, model=model)
      call check(error,.not. err%has_error(), "Li+: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, model%converged .and. model%iterations == 0, "a lone atom is not fitted")
      if (allocated(error)) return
      call check(error, q(1), 1.0_dp, "a lone Li+ carries its whole charge", thr=SUM_TOL)
      if (allocated(error)) return

      call model%destroy()
      call check(error,.not. allocated(model%population) .and. model%n_shells == 0, &
                 "destroy must release the shells")
      call mol%destroy()
   end subroutine test_bare_cation

end module test_mqc_mbis

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mbis, only: collect_mqc_mbis_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_mbis", collect_mqc_mbis_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
