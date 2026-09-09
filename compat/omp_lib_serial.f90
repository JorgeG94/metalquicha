module omp_lib
   !! The OpenMP runtime API, answering as one thread, for `MQC_ENABLE_SERIAL`.
   !!
   !! A serial build compiles without the OpenMP flag, so every `!$omp`
   !! directive is a comment and every `!$` sentinel line disappears. What does
   !! not disappear is the twenty-odd `use omp_lib` statements outside those
   !! sentinels: the module is still needed at compile time, and its symbols at
   !! link time. gfortran ships `omp_lib.mod` unconditionally but leaves
   !! `omp_get_max_threads_` in libgomp, so a serial link fails without it.
   !!
   !! Hence this. Linking libgomp instead would work and would be one line, but
   !! it puts back the dependency the serial build exists to remove, and it
   !! answers `omp_get_max_threads` with the machine's core count in a build
   !! that has no parallel region to run on them. One is the true answer here.
   !!
   !! Named `omp_lib` on purpose -- the point is that no call site changes. It
   !! lives outside `src/` because fpm globs `src/` and would compile it into
   !! every build, shadowing the real runtime; CMake adds it only when
   !! `MQC_ENABLE_SERIAL` is on, and the real `omp_lib` is used otherwise.
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none
   private

   public :: omp_get_max_threads, omp_get_num_threads, omp_get_thread_num
   public :: omp_set_num_threads, omp_set_max_active_levels
   public :: omp_lock_kind, omp_init_lock, omp_set_lock, omp_unset_lock, &
             omp_destroy_lock

   integer, parameter :: omp_lock_kind = int64
      !! Kind of the opaque lock handle, matching the real runtime's width.

contains

   function omp_get_max_threads() result(n)
      !! Threads a parallel region would use: one, there being no such region.
      integer :: n
      n = 1
   end function omp_get_max_threads

   function omp_get_num_threads() result(n)
      !! Threads in the current team. One, always, outside a parallel region.
      integer :: n
      n = 1
   end function omp_get_num_threads

   function omp_get_thread_num() result(n)
      !! This thread's index in its team. Zero, and it is the only member.
      integer :: n
      n = 0
   end function omp_get_thread_num

   subroutine omp_set_num_threads(n)
      !! Accepted and ignored. A serial build has one thread to give.
      integer, intent(in) :: n
      integer :: unused
      unused = n
   end subroutine omp_set_num_threads

   subroutine omp_set_max_active_levels(n)
      !! Accepted and ignored, as above: nesting needs a region to nest in.
      integer, intent(in) :: n
      integer :: unused
      unused = n
   end subroutine omp_set_max_active_levels

   ! The four lock routines below are reached only from `!$`-sentinel lines,
   ! which a serial build drops, so nothing here calls them today. They are
   ! carried anyway because the sources do name them, and a lock with one thread
   ! contending for it is trivially satisfied -- which is what these do.

   subroutine omp_init_lock(lock)
      !! Initialize a lock. Held by no one, and no one else can ask.
      integer(omp_lock_kind), intent(out) :: lock
      lock = 0_int64
   end subroutine omp_init_lock

   subroutine omp_set_lock(lock)
      !! Acquire a lock. Uncontended by construction, so this always succeeds.
      integer(omp_lock_kind), intent(inout) :: lock
      lock = 1_int64
   end subroutine omp_set_lock

   subroutine omp_unset_lock(lock)
      !! Release a lock.
      integer(omp_lock_kind), intent(inout) :: lock
      lock = 0_int64
   end subroutine omp_unset_lock

   subroutine omp_destroy_lock(lock)
      !! Discard a lock.
      integer(omp_lock_kind), intent(inout) :: lock
      lock = 0_int64
   end subroutine omp_destroy_lock

end module omp_lib
