module test_mqc_term_list
   !! `generate_mbe_term_list`, which a frozen optimization rests on
   !!
   !! The driver and `mqc_geometry_optimizer` both call this, and a geometry
   !! optimization freezes whatever it returns and then uses it for every step
   !! of the run. Two properties therefore have to hold, and neither is
   !! self-evident from reading the code:
   !!
   !!   * the same geometry gives the same list, or the frozen list is not the
   !!     list the run would otherwise have used
   !!   * a moved geometry can give a *different* list, which is the whole
   !!     reason freezing exists
   !!
   !! No DL-FIND needed: this is the fragment layer, so it runs in CI whatever
   !! the optimizer backend is set to.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_frag_utils, only: generate_mbe_term_list, get_nfrags
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_adapter, only: driver_config_t
   use pic_types, only: dp, int64
   implicit none
   private
   public :: collect_mqc_term_list_tests

contains

   subroutine collect_mqc_term_list_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("all_terms_without_cutoffs", test_all_terms_without_cutoffs), &
                  new_unittest("deterministic", test_deterministic), &
                  new_unittest("screening_removes_terms", test_screening_removes_terms), &
                  new_unittest("screened_list_keeps_monomers", test_screened_keeps_monomers), &
                  new_unittest("list_moves_with_geometry", test_list_moves_with_geometry) &
                  ]
   end subroutine collect_mqc_term_list_tests

   subroutine make_chain(sys_geom, spacing)
      !! Four one-atom monomers in a line, `spacing` Bohr apart
      !!
      !! One atom per monomer keeps the inter-monomer distance equal to the
      !! atom separation, so a cutoff in the test means exactly what it says.
      type(system_geometry_t), intent(out) :: sys_geom
      real(dp), intent(in) :: spacing

      integer :: i

      sys_geom%n_monomers = 4
      sys_geom%atoms_per_monomer = 1
      sys_geom%total_atoms = 4
      sys_geom%charge = 0
      sys_geom%multiplicity = 1
      allocate (sys_geom%element_numbers(4))
      allocate (sys_geom%coordinates(3, 4))
      sys_geom%element_numbers = 10  ! neon: closed shell, never bonded
      sys_geom%coordinates = 0.0_dp
      do i = 1, 4
         sys_geom%coordinates(1, i) = real(i - 1, dp)*spacing
      end do
   end subroutine make_chain

   subroutine test_all_terms_without_cutoffs(error)
      !! With no cutoffs the list is the full combinatorial one
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms

      call make_chain(sys_geom, 4.0_dp)
      config%nlevel = 2

      call generate_mbe_term_list(sys_geom, config, 2, polymers, n_terms)

      ! 4 monomers + 6 pairs
      call check(error, n_terms == 10_int64, &
                 "MBE(2) over 4 monomers is 4 monomers and 6 dimers")
      if (allocated(error)) return

      call check(error, n_terms == get_nfrags(4, 2), &
                 "the count should match get_nfrags, which sizes the array")
   end subroutine test_all_terms_without_cutoffs

   subroutine test_deterministic(error)
      !! The same geometry gives the same list, term for term
      !!
      !! This is what makes freezing sound: the list captured once at the start
      !! is the list every later step would have generated for that geometry.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: first(:, :), second(:, :)
      integer(int64) :: n_first, n_second

      call make_chain(sys_geom, 4.0_dp)
      config%nlevel = 2
      allocate (config%fragment_cutoffs(2))
      config%fragment_cutoffs = [0.0_dp, 5.0_dp]

      call generate_mbe_term_list(sys_geom, config, 2, first, n_first)
      call generate_mbe_term_list(sys_geom, config, 2, second, n_second)

      call check(error, n_first == n_second, "the term count should not change between calls")
      if (allocated(error)) return

      call check(error, all(first(1:n_first, :) == second(1:n_second, :)), &
                 "the terms should be identical, in the same order")
   end subroutine test_deterministic

   subroutine test_screening_removes_terms(error)
      !! A cutoff shorter than the chain drops the distant pairs
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms

      ! Neighbours 4 Bohr apart, so the pairs are at 4, 8 and 12 Bohr. The
      ! cutoff is in Angstrom, and 4 Bohr is about 2.117 Angstrom, so a cutoff
      ! of 3 Angstrom keeps only the three adjacent pairs.
      call make_chain(sys_geom, 4.0_dp)
      config%nlevel = 2
      allocate (config%fragment_cutoffs(2))
      config%fragment_cutoffs = [0.0_dp, 3.0_dp]

      call generate_mbe_term_list(sys_geom, config, 2, polymers, n_terms)

      call check(error, n_terms == 7_int64, &
                 "4 monomers and the 3 adjacent dimers should survive a 3 Angstrom cutoff")
   end subroutine test_screening_removes_terms

   subroutine test_screened_keeps_monomers(error)
      !! Screening never removes a monomer
      !!
      !! Monomers carry the leading term of the expansion, so a screen that
      !! dropped one would not be an approximation but a different system.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms, iterm
      integer :: n_monomer_terms

      call make_chain(sys_geom, 20.0_dp)  ! far enough that every pair is screened
      config%nlevel = 2
      allocate (config%fragment_cutoffs(2))
      config%fragment_cutoffs = [0.0_dp, 2.0_dp]

      call generate_mbe_term_list(sys_geom, config, 2, polymers, n_terms)

      n_monomer_terms = 0
      do iterm = 1, n_terms
         if (count(polymers(iterm, :) > 0) == 1) n_monomer_terms = n_monomer_terms + 1
      end do

      call check(error, n_monomer_terms == 4, &
                 "every monomer should survive however tight the cutoff")
      if (allocated(error)) return

      call check(error, n_terms == 4_int64, &
                 "and with every pair beyond the cutoff, nothing else should")
   end subroutine test_screened_keeps_monomers

   subroutine test_list_moves_with_geometry(error)
      !! The same system at a different geometry can give a different list
      !!
      !! The reason `keywords.optimization.freeze_terms` exists. Move the
      !! monomers closer and a pair that was screened out comes back, which
      !! mid-optimization is a step in the energy the optimizer reads as real.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: spread_out, close_up
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_spread, n_close

      config%nlevel = 2
      allocate (config%fragment_cutoffs(2))
      config%fragment_cutoffs = [0.0_dp, 3.0_dp]

      call make_chain(spread_out, 4.0_dp)
      call generate_mbe_term_list(spread_out, config, 2, polymers, n_spread)

      call make_chain(close_up, 2.0_dp)
      call generate_mbe_term_list(close_up, config, 2, polymers, n_close)

      call check(error, n_close > n_spread, &
                 "bringing the monomers together should bring screened pairs back")
   end subroutine test_list_moves_with_geometry

end module test_mqc_term_list

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_term_list, only: collect_mqc_term_list_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_term_list", collect_mqc_term_list_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
