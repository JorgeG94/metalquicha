!! Manual check that mqc can set CREST's sampling up and run it
!!
!!     cmake -B build -DMQC_ENABLE_CREST=ON
!!     cmake --build build --target check_crest_search
!!     ./build/check_crest_search
!!
!! `check_crest_callback` proved one gradient crosses the boundary. This is the
!! next question: whether a whole sampling run can be *configured* from mqc,
!! rather than from CREST's command line, and driven with mqc supplying every
!! gradient.
!!
!! **The settings come from CREST's own parser, not from here.** `systemdata`
!! has on the order of a hundred fields and the algorithms read far more of
!! them than the entry point names; filling it by hand would be a slow way of
!! reimplementing `confparse` and getting the defaults subtly wrong. So this
!! builds an argument vector *in memory* -- no command line, no file -- and
!! hands it to `parseflags`, which is exactly what `crest_main` does with the
!! real argv. The calculator is then overridden to point at this program.
!!
!! `parseflags` is a bare external subroutine with a deferred-length array
!! dummy, so it is called through an implicit interface here. That is what
!! `crest_main.f90` does too; there is no module wrapping it to use instead.
program check_crest_search
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use crest_data, only: systemdata, timer
   use crest_calculator, only: jobtype
   use lj, only: lj_engrad
   implicit none

   !> Explicit interfaces for two bare external subroutines. `crest_main.f90`
   !> calls both through an implicit interface; declaring them properly here
   !> costs nothing and makes the compiler check the call, which matters most
   !> for `arg`, an assumed-length explicit-shape array receiving a
   !> deferred-length allocatable.
   interface
      subroutine parseflags(env, arg, nra)
         import :: systemdata
         implicit none
         type(systemdata), intent(inout) :: env
         integer, intent(in) :: nra
         !> No intent, deliberately: the definition in confparse.f90 declares
         !> it without one, and intent is part of a procedure's
         !> characteristics -- adding it here would make this interface
         !> disagree with the procedure it describes.
         ! allow(missing-intent)
         character(len=*) :: arg(nra)
      end subroutine parseflags

      subroutine crest_search_imtdgc(env, tim)
         import :: systemdata, timer
         implicit none
         type(systemdata), intent(inout) :: env
         type(timer), intent(inout) :: tim
      end subroutine crest_search_imtdgc
   end interface

   type(systemdata) :: env
   type(timer) :: tim
   character(len=:), allocatable :: argv(:)
   integer :: nargs, unit, i, j
   integer :: ncalls = 0
   real(wp), allocatable :: r0(:, :)
   real(wp) :: d(3)

   !> The starting structure is a file, and that is not the file-based
   !> gradient exchange this design exists to avoid: it is read once, at
   !> setup, the way any program reads its input geometry. Every gradient
   !> after this crosses in memory.
   open (newunit=unit, file="crest_probe.xyz", status="replace", action="write")
   write (unit, "(a)") "3"
   write (unit, "(a)") "water"
   write (unit, "(a)") "O   0.0000000   0.0000000   0.1173000"
   write (unit, "(a)") "H   0.0000000   0.7572000  -0.4692000"
   write (unit, "(a)") "H   0.0000000  -0.7572000  -0.4692000"
   close (unit)

   call tim%init(20)

   nargs = 1
   allocate (character(len=64) :: argv(nargs))
   argv(1) = "crest_probe.xyz"

   call parseflags(env, argv, nargs)

   write (*, "(a)") ""
   write (*, "(a)") "  parseflags returned; what it built:"
   write (*, "(a,i0)") "    calculation levels   ", env%calc%ncalculations
   write (*, "(a,i0)") "    runver               ", env%runver
   write (*, "(a,i0)") "    threads              ", env%threads
   write (*, "(a,l1)") "    legacy               ", env%legacy

   if (env%calc%ncalculations < 1) then
      write (*, "(a)") "  FAIL: no calculation level was set up to override"
      error stop 1
   end if

   write (*, "(a,i0)") "    level 1 jobtype id   ", env%calc%calcs(1)%id

   !> Swap whatever level parseflags chose for this program's own hook. This
   !> is the whole point: CREST keeps its algorithm, its defaults and its
   !> bookkeeping, and every energy it asks for comes from here.
   env%calc%calcs(1)%id = jobtype%callback
   env%calc%calcs(1)%external_engrad => probe_engrad
   call env%calc%calcs(1)%printid(1, 1)

   write (*, "(a)") ""
   !> spring lengths from the structure parseflags read, so the model's
   !> minimum is where CREST starts
   allocate (r0(env%ref%nat, env%ref%nat), source=0.0_wp)
   do i = 1, env%ref%nat
      do j = 1, env%ref%nat
         d = env%ref%xyz(:, i) - env%ref%xyz(:, j)
         r0(i, j) = sqrt(sum(d*d))
      end do
   end do

   write (*, "(a)") "  running crest_search_imtdgc with this program as the calculator"
   ncalls = 0
   call crest_search_imtdgc(env, tim)

   write (*, "(a)") ""
   write (*, "(a,i0)") "  search returned, iostatus_meta = ", env%iostatus_meta
   write (*, "(a,i0)") "  callback invocations          = ", ncalls
   write (*, "(a)") "  OK: CREST sampled, and every gradient came from here"

contains

   subroutine probe_engrad(nat, at, xyz, energy, grad, iostatus)
      !! A pairwise harmonic model, standing in for mqc's SCF
      !!
      !! Cheap on purpose: a full iMTD-GC run asks for tens of thousands of
      !! gradients, and what is under test here is whether the run can be
      !! configured from memory and driven through the hook -- not the physics.
      !! `check_crest_callback` already establishes that mqc's own numbers
      !! survive the crossing.
      !!
      !! Every pair is a spring at its distance in the input geometry, so the
      !! minimum is the structure CREST started from. That matters: a
      !! Lennard-Jones gas was the first attempt and it optimised 0 of 686
      !! sampled structures, because sigma = 3.4 A against an O-H bond of
      !! 0.96 A is deep in the repulsive wall. The plumbing was fine -- 71692
      !! callback invocations -- and the ensemble was still empty. A potential
      !! a molecule can actually sit in is the difference between testing the
      !! integration and testing the toy.
      use, intrinsic :: iso_fortran_env, only: wp => real64
      integer, intent(in) :: nat
      integer, intent(in) :: at(nat)
      real(wp), intent(in) :: xyz(3, nat)
      real(wp), intent(out) :: energy
      real(wp), intent(out) :: grad(3, nat)
      integer, intent(out) :: iostatus

      real(wp), parameter :: K_BOND = 0.5_wp
      real(wp) :: dx(3)
      real(wp) :: r, dr, pref
      integer :: i, j

      iostatus = 0
      energy = 0.0_wp
      grad = 0.0_wp
      !$omp atomic
      ncalls = ncalls + 1
      if (nat < 1 .or. at(1) < 1) then
         iostatus = 1
         return
      end if
      if (.not. allocated(r0)) then
         iostatus = 2
         return
      end if

      do i = 1, nat
         do j = i + 1, nat
            dx = xyz(:, i) - xyz(:, j)
            r = sqrt(sum(dx*dx))
            if (r < 1.0e-8_wp) cycle
            dr = r - r0(i, j)
            energy = energy + 0.5_wp*K_BOND*dr*dr
            pref = K_BOND*dr/r
            grad(:, i) = grad(:, i) + pref*dx
            grad(:, j) = grad(:, j) - pref*dx
         end do
      end do
   end subroutine probe_engrad

end program check_crest_search
