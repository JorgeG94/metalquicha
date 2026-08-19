module test_mqc_gpu_team
   !! What the multi-GPU mode gets right without a GPU in the room.
   !!
   !! Two things here are worth more than they look. The first is the radial
   !! partition: the whole claim of the mode is that slicing the quadrature
   !! across ranks and adding the pieces back reproduces the single-GPU answer
   !! *exactly*, not nearly. That claim is arithmetic on the mesh and needs no
   !! device, so it is checked here rather than discovered later as a two
   !! millihartree disagreement between one GPU and four.
   !!
   !! The second is the memory-policy name table. `parse_df_derivative_memory`
   !! repeats cuEST's enumerator values by hand, because the module it lives in
   !! is compiled on the fpm path where the cuEST bindings do not exist. A
   !! hand-copied enum is exactly the kind of thing that rots silently when the
   !! vendor renumbers, so the values are pinned here; the companion assertion
   !! against the real bindings lives in test_mqc_cuest_memory_policy, which
   !! only builds when they are present.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_gpu_team, only: gpu_team_active, gpu_team_rank, gpu_team_size, &
                           gpu_team_reduce, gpu_team_reduce_scalar, gpu_team_agree, &
                           gpu_team_disable
   use mqc_cuest_iface, only: parse_df_derivative_memory
   implicit none
   private
   public :: collect_mqc_gpu_team_tests

contains

   subroutine collect_mqc_gpu_team_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("radial shards partition the mesh", test_shards_partition), &
                  new_unittest("radial shards are balanced", test_shards_balanced), &
                  new_unittest("one shard is the whole mesh", test_single_shard), &
                  new_unittest("dormant team is the identity", test_dormant_team), &
                  new_unittest("memory policy names", test_memory_policy_names), &
                  new_unittest("memory policy rejects nonsense", test_memory_policy_bad) &
                  ]
   end subroutine collect_mqc_gpu_team_tests

   subroutine test_shards_partition(error)
      !! Every radial point belongs to exactly one shard, and the shards
      !! together are the whole mesh
      !!
      !! This is the multi-GPU mode's correctness argument in one assertion.
      !! The Becke weight at a point does not depend on which other points
      !! exist, so if the point sets partition, the quadrature sums do too --
      !! and the reduced result is bit-for-bit the single-GPU one rather than
      !! an approximation to it.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N_RADIAL = 75
      integer, parameter :: N_SHARDS = 4
      logical :: seen(N_RADIAL)
      integer :: shard, i, total

      seen = .false.
      total = 0
      do shard = 0, N_SHARDS - 1
         do i = shard + 1, N_RADIAL, N_SHARDS
            call check(error,.not. seen(i), "radial point claimed by two shards")
            if (allocated(error)) return
            seen(i) = .true.
            total = total + 1
         end do
      end do

      call check(error, total, N_RADIAL)
      if (allocated(error)) return
      call check(error, all(seen), "a radial point belongs to no shard")
   end subroutine test_shards_partition

   subroutine test_shards_balanced(error)
      !! No shard is more than one point larger than any other
      !!
      !! Striding rather than blocking is what buys this. It matters for more
      !! than tidiness: the ranks run in lockstep and reduce every iteration,
      !! so the slowest shard sets the pace, and a blocked split would hand one
      !! rank the dense core of every atom and another the vacuum.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N_RADIAL = 75
      integer, parameter :: N_SHARDS = 8
      integer :: shard, counts(0:N_SHARDS - 1), formula

      do shard = 0, N_SHARDS - 1
         counts(shard) = count_shard(N_RADIAL, shard, N_SHARDS)
         ! The closed form the backend uses to size its arrays, checked against
         ! the loop that actually fills them -- a mismatch here is an
         ! out-of-bounds write on the GPU path.
         formula = (N_RADIAL - shard - 1)/N_SHARDS + 1
         call check(error, counts(shard), formula)
         if (allocated(error)) return
      end do

      call check(error, sum(counts), N_RADIAL)
      if (allocated(error)) return
      call check(error, maxval(counts) - minval(counts) <= 1, &
                 "radial shards differ by more than one point")
   end subroutine test_shards_balanced

   subroutine test_single_shard(error)
      !! One shard keeps every point, so a single-rank run is untouched
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N_RADIAL = 75

      call check(error, count_shard(N_RADIAL, 0, 1), N_RADIAL)
      if (allocated(error)) return
      call check(error, (N_RADIAL - 0 - 1)/1 + 1, N_RADIAL)
   end subroutine test_single_shard

   subroutine test_dormant_team(error)
      !! With no team live, every collective is the identity
      !!
      !! The reductions are called unconditionally from the SCF, so a serial
      !! run passes through them on every iteration. If they were not exactly
      !! the identity, single-GPU answers would change the day this mode
      !! landed -- which is the regression this catches.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: array(3), scalar
      logical :: flag

      call gpu_team_disable()

      call check(error,.not. gpu_team_active(), "a disabled team reports itself active")
      if (allocated(error)) return
      call check(error, gpu_team_size(), 1)
      if (allocated(error)) return
      call check(error, gpu_team_rank(), 0)
      if (allocated(error)) return

      array = [1.0_dp, -2.5_dp, 3.25_dp]
      call gpu_team_reduce(array)
      call check(error, array(1), 1.0_dp)
      if (allocated(error)) return
      call check(error, array(2), -2.5_dp)
      if (allocated(error)) return
      call check(error, array(3), 3.25_dp)
      if (allocated(error)) return

      scalar = -7.5_dp
      call gpu_team_reduce_scalar(scalar)
      call check(error, scalar, -7.5_dp)
      if (allocated(error)) return

      ! Both branches: a dormant `agree` must not turn a false into a true, and
      ! must not turn a true into a false either -- it would stall or shorten
      ! every serial SCF respectively.
      flag = .true.
      call gpu_team_agree(flag)
      call check(error, flag, "dormant agree changed a true")
      if (allocated(error)) return

      flag = .false.
      call gpu_team_agree(flag)
      call check(error,.not. flag, "dormant agree changed a false")
   end subroutine test_dormant_team

   subroutine test_memory_policy_names(error)
      !! Each accepted name maps to cuEST's enumerator, and 'auto' to -1
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: parse_error
      integer :: policy

      call parse_df_derivative_memory("auto", policy, parse_error)
      call check(error,.not. parse_error%has_error(), "'auto' was rejected")
      if (allocated(error)) return
      call check(error, policy, -1)
      if (allocated(error)) return

      ! -1 is not a cuEST value; it is this code's "set nothing at all", which
      ! has to stay distinguishable from DEVICECACHE = 0. Conflating them would
      ! mean every gradient silently pinned to one policy.
      call parse_df_derivative_memory("", policy, parse_error)
      call check(error, policy, -1)
      if (allocated(error)) return

      call parse_df_derivative_memory("devicecache", policy, parse_error)
      call check(error, policy, 0)
      if (allocated(error)) return

      call parse_df_derivative_memory("hostcache", policy, parse_error)
      call check(error, policy, 1)
      if (allocated(error)) return

      call parse_df_derivative_memory("overwrite", policy, parse_error)
      call check(error, policy, 2)
      if (allocated(error)) return

      call parse_df_derivative_memory("recompute", policy, parse_error)
      call check(error, policy, 3)
      if (allocated(error)) return

      ! Case and surrounding space are a deck author's business, not cuEST's.
      call parse_df_derivative_memory("  OverWrite  ", policy, parse_error)
      call check(error,.not. parse_error%has_error(), "a spelled-out name was rejected")
      if (allocated(error)) return
      call check(error, policy, 2)
   end subroutine test_memory_policy_names

   subroutine test_memory_policy_bad(error)
      !! An unknown policy is refused rather than defaulted
      !!
      !! Falling back to 'auto' would mean a typo bought cuEST's default and
      !! said nothing -- on the one keyword whose entire purpose is to be
      !! reached for when the default has already run out of memory.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: parse_error
      integer :: policy

      call parse_df_derivative_memory("hostcahce", policy, parse_error)
      call check(error, parse_error%has_error(), "a misspelled policy was accepted")
      if (allocated(error)) return
      call check(error, index(parse_error%get_message(), "hostcahce") > 0, &
                 "the refusal does not quote the name it refused")
   end subroutine test_memory_policy_bad

   pure integer function count_shard(n_radial, shard, n_shards)
      !! How many radial points a shard actually receives, by the same loop the
      !! backend uses to fill them
      integer, intent(in) :: n_radial, shard, n_shards
      integer :: i

      count_shard = 0
      do i = shard + 1, n_radial, n_shards
         count_shard = count_shard + 1
      end do
   end function count_shard

end module test_mqc_gpu_team

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_gpu_team, only: collect_mqc_gpu_team_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_gpu_team", collect_mqc_gpu_team_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
