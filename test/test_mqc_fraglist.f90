module test_mqc_fraglist
   !! Tests the fragment term list, and the round trip the Python interface
   !! will be made of: generate a list, filter it, hand the survivors back.
   !!
   !! The list exists as an object so that the screening criterion can live
   !! outside it. Distance screening used to be inseparable from generation,
   !! which is why energy-based fragmentation was not possible; these tests are
   !! mostly about that separation holding -- that a filtered list is a
   !! well-formed list, and that what a caller puts back is what comes out.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_fraglist, only: fraglist_t
   use mqc_error, only: error_t
   use pic_types, only: dp, default_int, int64
   implicit none
   private
   public :: collect_mqc_fraglist_tests

contains

   subroutine collect_mqc_fraglist_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("generate_dimers", test_generate_dimers), &
                  new_unittest("generate_up_to_trimers", test_generate_trimers), &
                  new_unittest("terms_are_zero_padded", test_padding), &
                  new_unittest("keep_selects_and_reorders_nothing", test_keep), &
                  new_unittest("keep_everything_is_a_no_op", test_keep_all), &
                  new_unittest("keep_nothing_leaves_an_empty_list", test_keep_none), &
                  new_unittest("replace_round_trips", test_replace), &
                  new_unittest("replace_drops_stale_distances", test_replace_clears), &
                  new_unittest("rejects_impossible_levels", test_bad_level), &
                  new_unittest("rejects_a_mismatched_mask", test_bad_mask) &
                  ]
   end subroutine collect_mqc_fraglist_tests

   subroutine test_generate_dimers(error)
      !! Five monomers at level 2 is C(5,2) = 10 pairs
      type(error_type), allocatable, intent(out) :: error
      type(fraglist_t) :: list
      type(error_t) :: err

      call list%create(5_default_int, 2_default_int, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, int(list%n_terms), 10, "C(5,2) = 10 dimers")
      if (allocated(error)) return
      call check(error, int(list%max_level), 2)
      if (allocated(error)) return
      call check(error, int(list%level_of(1_int64)), 2, "a dimer has two monomers")

      call list%destroy()
   end subroutine test_generate_dimers

   subroutine test_generate_trimers(error)
      !! Level 3 gives pairs and triples together: C(5,2) + C(5,3) = 10 + 10
      type(error_type), allocatable, intent(out) :: error
      type(fraglist_t) :: list
      type(error_t) :: err

      call list%create(5_default_int, 3_default_int, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, int(list%n_terms), 20, "C(5,2) + C(5,3) = 20 terms")

      call list%destroy()
   end subroutine test_generate_trimers

   subroutine test_padding(error)
      !! A dimer inside a level-3 list occupies two slots and zeroes the third
      !!
      !! Everything downstream reads the order off the non-zero count, so a
      !! stray value in the padding would silently promote a dimer to a trimer.
      type(error_type), allocatable, intent(out) :: error
      type(fraglist_t) :: list
      type(error_t) :: err
      integer(int64) :: iterm
      integer :: n_dimers, n_trimers

      call list%create(5_default_int, 3_default_int, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return

      n_dimers = 0
      n_trimers = 0
      do iterm = 1, list%n_terms
         select case (int(list%level_of(iterm)))
         case (2)
            n_dimers = n_dimers + 1
            call check(error, list%terms(iterm, 3), 0, &
                       "a dimer's third slot must be zero")
            if (allocated(error)) return
         case (3)
            n_trimers = n_trimers + 1
         case default
            call check(error, .false., "a term of unexpected order appeared")
            return
         end select
         ! Monomer indices are 1-based, so an occupied slot is never zero.
         call check(error, all(list%terms(iterm, 1:list%level_of(iterm)) > 0), &
                    "occupied slots must hold 1-based monomer indices")
         if (allocated(error)) return
      end do

      call check(error, n_dimers, 10, "ten dimers")
      if (allocated(error)) return
      call check(error, n_trimers, 10, "ten trimers")

      call list%destroy()
   end subroutine test_padding

   subroutine test_keep(error)
      !! The operation every screen reduces to
      type(error_type), allocatable, intent(out) :: error
      type(fraglist_t) :: list
      type(error_t) :: err
      logical, allocatable :: mask(:)
      integer(default_int) :: first_kept(3), third_kept(3)

      call list%create(5_default_int, 3_default_int, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return

      ! Keep terms 2, 5 and 9, and remember what they were.
      allocate (mask(list%n_terms))
      mask = .false.
      mask(2) = .true.
      mask(5) = .true.
      mask(9) = .true.
      first_kept = list%terms(2, :)
      third_kept = list%terms(9, :)

      call list%keep(mask, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, int(list%n_terms), 3, "three terms survive")
      if (allocated(error)) return

      ! Survivors keep their contents and their relative order.
      call check(error, all(list%terms(1, :) == first_kept), &
                 "the first survivor should be the old term 2")
      if (allocated(error)) return
      call check(error, all(list%terms(3, :) == third_kept), &
                 "the last survivor should be the old term 9")

      deallocate (mask)
      call list%destroy()
   end subroutine test_keep

   subroutine test_keep_all(error)
      type(error_type), allocatable, intent(out) :: error
      type(fraglist_t) :: list
      type(error_t) :: err
      logical, allocatable :: mask(:)
      integer(int64) :: before

      call list%create(6_default_int, 2_default_int, err)
      before = list%n_terms
      allocate (mask(before))
      mask = .true.

      call list%keep(mask, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, int(list%n_terms), int(before), "keeping everything changes nothing")

      deallocate (mask)
      call list%destroy()
   end subroutine test_keep_all

   subroutine test_keep_none(error)
      !! An empty list is a legitimate outcome of a strict screen
      type(error_type), allocatable, intent(out) :: error
      type(fraglist_t) :: list
      type(error_t) :: err
      logical, allocatable :: mask(:)

      call list%create(6_default_int, 2_default_int, err)
      allocate (mask(list%n_terms))
      mask = .false.

      call list%keep(mask, err)
      call check(error, .not. err%has_error(), &
                 "screening everything away should not be an error: "//err%get_message())
      if (allocated(error)) return
      call check(error, int(list%n_terms), 0, "no terms survive")

      deallocate (mask)
      call list%destroy()
   end subroutine test_keep_none

   subroutine test_replace(error)
      !! What a caller puts back is what comes out
      type(error_type), allocatable, intent(out) :: error
      type(fraglist_t) :: list
      type(error_t) :: err
      integer(default_int) :: mine(3, 2)

      ! Three dimers of a caller's own choosing, in a caller's own order.
      mine(1, :) = [4, 7]
      mine(2, :) = [1, 2]
      mine(3, :) = [3, 9]

      call list%replace(mine, 3_int64, 2_default_int, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, int(list%n_terms), 3)
      if (allocated(error)) return
      call check(error, int(list%max_level), 2)
      if (allocated(error)) return
      call check(error, all(list%terms(1, :) == [4, 7]), "term 1 should be [4, 7]")
      if (allocated(error)) return
      call check(error, all(list%terms(3, :) == [3, 9]), "term 3 should be [3, 9]")
      if (allocated(error)) return
      ! Not reordered or sorted on the way in -- a caller's order is meaningful,
      ! and a restart depends on it.
      call check(error, all(list%terms(2, :) == [1, 2]), &
                 "the caller's ordering should be preserved")

      call list%destroy()
   end subroutine test_replace

   subroutine test_replace_clears(error)
      !! Distances describe the old list, so they must not survive the new one
      type(error_type), allocatable, intent(out) :: error
      type(fraglist_t) :: list
      type(error_t) :: err
      integer(default_int) :: mine(1, 2)

      call list%create(4_default_int, 2_default_int, err)
      ! Pretend a measure has happened; the point is that replace clears it.
      allocate (list%distances(list%n_terms))
      list%distances = 1.0_dp
      list%has_distances = .true.

      mine(1, :) = [1, 2]
      call list%replace(mine, 1_int64, 2_default_int, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, .not. list%has_distances, &
                 "distances from the old list must not be read as the new one's")
      if (allocated(error)) return
      call check(error, .not. allocated(list%distances), &
                 "the stale distance array should be gone")

      call list%destroy()
   end subroutine test_replace_clears

   subroutine test_bad_level(error)
      !! Levels that cannot produce terms are rejected, not silently emptied
      type(error_type), allocatable, intent(out) :: error
      type(fraglist_t) :: list
      type(error_t) :: err

      call list%create(5_default_int, 1_default_int, err)
      call check(error, err%has_error(), "level 1 has no n-mer terms and should be rejected")
      if (allocated(error)) return

      err = error_t()
      call list%create(3_default_int, 5_default_int, err)
      call check(error, err%has_error(), &
                 "a level above the monomer count should be rejected")
      if (allocated(error)) return

      err = error_t()
      call list%create(0_default_int, 2_default_int, err)
      call check(error, err%has_error(), "a system with no monomers should be rejected")

      call list%destroy()
   end subroutine test_bad_level

   subroutine test_bad_mask(error)
      type(error_type), allocatable, intent(out) :: error
      type(fraglist_t) :: list
      type(error_t) :: err
      logical :: mask(3)

      call list%create(5_default_int, 2_default_int, err)
      mask = .true.
      call list%keep(mask, err)
      call check(error, err%has_error(), &
                 "a mask that is not one entry per term should be rejected")

      call list%destroy()
   end subroutine test_bad_mask

end module test_mqc_fraglist

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fraglist, only: collect_mqc_fraglist_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_fraglist", collect_mqc_fraglist_tests) &
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
