!! The accurate atomic minimal basis and the dimensions it fixes
module test_mqc_aambs
   !! Everything here is checkable without an integral, a basis set of the
   !! working kind, or an SCF -- which is the point. The quasi-atomic
   !! construction defines each of its orbital spaces by an exact dimension
   !! rather than by a threshold, so the counting is load-bearing: get
   !! `N_VVO = N_mbs - N_occ` wrong and the SVD selects a different number of
   !! vectors and the whole analysis is of a different space, silently.
   !!
   !! Paper I (West, Schmidt, Gordon, Ruedenberg, J. Chem. Phys. 139, 234107
   !! (2013)) tabulates these counts for eight molecules in its Table I, which
   !! makes them a published reference for a stage that has no other one.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_aambs, only: aambs_element_counts, aambs_dimensions, aambs_dimensions_t
   implicit none
   private

   public :: collect_mqc_aambs_tests

contains

   subroutine collect_mqc_aambs_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("element_counts_match_gamess", test_element_counts), &
                  new_unittest("paper_one_table_one", test_paper_one_table_one), &
                  new_unittest("core_is_semicore_aware", test_semicore), &
                  new_unittest("refusals", test_refusals) &
                  ]
   end subroutine collect_mqc_aambs_tests

   subroutine test_element_counts(error)
      !! Spot values against GAMESS's own tables
      !!
      !! `NVVOS_NUMCOR` (vvos.src) for the core and `LOCAL_NUMVAL` (locsvd.src)
      !! for the valence. Hydrogen has no core at all; carbon is the textbook
      !! 1s core with four valence; potassium is the case where the 3d shell
      !! exists in the basis but the valence is a single 4s.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer :: core, valence

      call aambs_element_counts(1, core, valence, err)
      call check(error,.not. err%has_error(), "hydrogen should be in the table")
      if (allocated(error)) return
      call check(error, core, 0, "hydrogen has no chemical core")
      if (allocated(error)) return
      call check(error, valence, 1, "hydrogen has one valence orbital")
      if (allocated(error)) return

      call aambs_element_counts(6, core, valence, err)
      call check(error, core, 1, "carbon has a 1s core")
      if (allocated(error)) return
      call check(error, valence, 4, "carbon has 2s2p valence")
      if (allocated(error)) return

      call aambs_element_counts(19, core, valence, err)
      call check(error, core, 9, "potassium has a nine-orbital core")
      if (allocated(error)) return
      call check(error, valence, 1, "potassium has a single 4s valence orbital")
      if (allocated(error)) return

      call aambs_element_counts(54, core, valence, err)
      call check(error,.not. err%has_error(), "xenon is the last element in the table")
      if (allocated(error)) return
      call check(error, core, 23, "xenon has a 23-orbital core")
      if (allocated(error)) return
      call check(error, valence, 4, "xenon has 5s5p valence")
   end subroutine test_element_counts

   subroutine test_semicore(error)
      !! The core is the *chemical* core, which is not the frozen core
      !!
      !! Scandium through zinc keep 3d in the valence, so the core stops at the
      !! argon shell: nine orbitals. Gallium through krypton have a filled and
      !! chemically inert 3d, which is counted as core: fourteen. A frozen-core
      !! count would say ten for both. This distinction is why the two counts
      !! are shipped as data rather than derived from the electron
      !! configuration, and getting it wrong moves four orbitals between the
      !! core and valence spaces for every fourth-row element.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer :: core, valence

      call aambs_element_counts(30, core, valence, err)   ! zinc
      call check(error, core, 9, "zinc keeps 3d in the valence")
      if (allocated(error)) return
      call check(error, valence, 6, "zinc has 3d and 4s in the valence")
      if (allocated(error)) return

      call aambs_element_counts(31, core, valence, err)   ! gallium
      call check(error, core, 14, "gallium counts its filled 3d as core")
      if (allocated(error)) return
      call check(error, valence, 4, "gallium has 4s4p valence")
   end subroutine test_semicore

   subroutine test_paper_one_table_one(error)
      !! Table I of Paper I, all eight molecules
      !!
      !! Columns are the inner-shell count, the filled valence count, the number
      !! of valence-virtual orbitals to be recovered, and the total minimal-basis
      !! dimension. The electron count is supplied rather than derived so the
      !! test states its own premises; the point of the check is the partition,
      !! not the arithmetic that produced the electron count.
      type(error_type), allocatable, intent(out) :: error

      ! HNO
      call one_row(error, [1, 7, 8], 16, 2, 6, 3, 11, "HNO, nitroxyl")
      if (allocated(error)) return
      ! NaCl
      call one_row(error, [11, 17], 28, 10, 4, 1, 15, "NaCl, table salt")
      if (allocated(error)) return
      ! HNCO
      call one_row(error, [1, 7, 6, 8], 22, 3, 8, 5, 16, "HNCO, isocyanic acid")
      if (allocated(error)) return
      ! HOCH=O, formic acid
      call one_row(error, [1, 8, 6, 1, 8], 24, 3, 9, 5, 17, "HCOOH, formic acid")
      if (allocated(error)) return
      ! H2Si=CH2, silene
      call one_row(error, [1, 1, 14, 6, 1, 1], 24, 6, 6, 6, 18, "H2SiCH2, silene")
      if (allocated(error)) return
      ! SO2
      call one_row(error, [16, 8, 8], 32, 7, 9, 3, 19, "SO2, sulfur dioxide")
      if (allocated(error)) return
      ! MnO4-, permanganate. One extra electron, hence 57 rather than 56.
      call one_row(error, [25, 8, 8, 8, 8], 58, 13, 16, 6, 35, "MnO4-, permanganate")
      if (allocated(error)) return
      ! As4S4, realgar
      call one_row(error, [33, 33, 33, 33, 16, 16, 16, 16], 196, 76, 22, 10, 108, &
                   "As4S4, realgar")
   end subroutine test_paper_one_table_one

   subroutine one_row(error, atomic_numbers, n_electrons, n_core, n_valocc, &
                      n_vvo, n_mbs, label)
      !! One row of Paper I Table I
      type(error_type), allocatable, intent(inout) :: error
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: n_electrons
      integer, intent(in) :: n_core, n_valocc, n_vvo, n_mbs
      character(len=*), intent(in) :: label

      type(aambs_dimensions_t) :: dims
      type(error_t) :: err

      call aambs_dimensions(atomic_numbers, n_electrons, dims, err)
      call check(error,.not. err%has_error(), label//": dimensions should be available")
      if (allocated(error)) return

      call check(error, dims%n_core, n_core, label//": inner-shell count")
      if (allocated(error)) return
      call check(error, dims%n_valocc, n_valocc, label//": filled valence count")
      if (allocated(error)) return
      call check(error, dims%n_vvo, n_vvo, label//": valence-virtual count")
      if (allocated(error)) return
      call check(error, dims%n_mbs, n_mbs, label//": minimal-basis dimension")
      if (allocated(error)) return

      ! The identity the whole partition rests on. Stated separately because it
      ! is what a wrong core or valence table breaks first, and because Table I
      ! is internally consistent only if it holds.
      call check(error, dims%n_core + dims%n_valocc + dims%n_vvo, dims%n_mbs, &
                 label//": core + filled valence + valence-virtual must be the "// &
                 "whole minimal basis")
   end subroutine one_row

   subroutine test_refusals(error)
      !! The three ways the counting can be asked something it cannot answer
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(aambs_dimensions_t) :: dims
      integer :: core, valence

      ! Beyond xenon there is no non-relativistic table.
      call aambs_element_counts(55, core, valence, err)
      call check(error, err%has_error(), "caesium is past the end of the table and "// &
                 "should be refused rather than silently given zero orbitals")
      if (allocated(error)) return
      call err%clear()

      ! An odd electron count is not a closed shell.
      call aambs_dimensions([8, 1], 9, dims, err)
      call check(error, err%has_error(), "an odd electron count should be refused")
      if (allocated(error)) return
      call err%clear()

      ! More occupied orbitals than the minimal basis can hold. Two hydrogens
      ! have two minimal-basis orbitals between them, so six electrons cannot
      ! fit -- n_vvo would come out negative.
      call aambs_dimensions([1, 1], 6, dims, err)
      call check(error, err%has_error(), "an occupied space larger than the minimal "// &
                 "basis should be refused rather than produce a negative dimension")
   end subroutine test_refusals

end module test_mqc_aambs

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_aambs, only: collect_mqc_aambs_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_aambs", collect_mqc_aambs_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
