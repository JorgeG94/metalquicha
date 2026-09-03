!! That `incremental_fock` reaches the SCF, counted rather than timed
module test_mqc_incremental_fock
   !! **The one SCF setting with no observable effect on any answer.**
   !!
   !! An incremental Fock build is exact to the convergence threshold: it
   !! rebuilds `G` from the density *change* rather than the density, and since
   !! `G` is linear in the density the two agree. So the energy, the iteration
   !! count, the orbitals and every property computed from them are identical
   !! whether the setting arrived or was dropped on the floor. Six other
   !! settings on this path were dropped for months and this is the one that
   !! could not have been caught by comparing numbers.
   !!
   !! That left timing, and timing turned out to be a poor thing to assert.
   !! Measured on the development machine: an SCF small enough for a unit test
   !! reports 0.000 s per Fock build either way, at STO-3G, 6-31G and cc-pVDZ
   !! alike. One large enough to separate them -- adenine in 6-31G, where the
   !! per-build time falls 0.48 s to 0.18 s with incremental building on and
   !! stays flat at 0.48 s with it off -- costs about nine seconds a run, and
   !! wall-clock on a shared machine is exactly the kind of assertion that
   !! fails for reasons that have nothing to do with the code.
   !!
   !! Capping the iteration count to make such a test affordable does not work
   !! either, and the reason is worth stating: the saving accrues *late*. The
   !! incremental path screens on the largest density element a quartet will
   !! multiply, and only once the SCF settles does the change become small
   !! enough to discard most of them. Truncating the run at ten iterations
   !! removes the effect being measured -- adenine capped that way takes 0.30 s
   !! of Fock time both ways.
   !!
   !! So the SCF counts its builds instead. `full_fock_builds` and
   !! `incremental_updates` are exact, free, and the same on every machine.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   implicit none
   private

   public :: collect_incremental_fock

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

   !! Iterations between full rebuilds when incremental building is on. Must
   !! match `INCREMENTAL_RESET` in `mqc_czt_rhf`; the tests below stay well
   !! inside one window so the arithmetic does not depend on it.
   integer, parameter :: RESET_WINDOW = 16

contains

   subroutine collect_incremental_fock(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("incremental_off_builds_every_fock_in_full", test_off), &
                  new_unittest("incremental_on_builds_one_and_updates", test_on), &
                  new_unittest("incremental_does_not_move_the_energy", test_same_energy) &
                  ]
   end subroutine collect_incremental_fock

   subroutine converge(scf, error, incremental)
      !! One water SCF, with incremental building on or off
      type(rhf_result_t), intent(out) :: scf
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: incremental

      type(czt_molecule_t) :: mol
      type(error_t) :: err

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return

      call run_czt_rhf(mol, 10, 50, 1.0e-9_dp, 1.0e-7_dp, .false., scf, err, &
                       incremental_fock=incremental)
      call check(error,.not. err%has_error(), "the SCF should run: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, scf%converged, "the SCF should converge")
   end subroutine converge

   subroutine test_off(error)
      !! Off: every iteration builds a Fock matrix in full, and none updates
      type(error_type), allocatable, intent(out) :: error

      type(rhf_result_t) :: scf

      call converge(scf, error, incremental=.false.)
      if (allocated(error)) return

      call check(error, scf%full_fock_builds == scf%iterations, &
                 "with incremental building off, every iteration is a full build")
      if (allocated(error)) return
      call check(error, scf%incremental_updates == 0, &
                 "and none of them is an update")
   end subroutine test_off

   subroutine test_on(error)
      !! On: one full build to seed the reference, then updates
      !!
      !! The count is one plus a rebuild every `RESET_WINDOW` iterations, and
      !! this SCF converges well inside one window, so it is exactly one. That
      !! is the assertion the timing could not make: a setting that changes no
      !! number still changes something countable.
      type(error_type), allocatable, intent(out) :: error

      type(rhf_result_t) :: scf

      call converge(scf, error, incremental=.true.)
      if (allocated(error)) return

      call check(error, scf%iterations < RESET_WINDOW, &
                 "this SCF must stay inside one reset window for the count below "// &
                 "to be exact; if it no longer does, the expectation needs the "// &
                 "periodic rebuilds added back")
      if (allocated(error)) return
      call check(error, scf%full_fock_builds == 1, &
                 "with incremental building on, only the first build is full")
      if (allocated(error)) return
      call check(error, scf%incremental_updates == scf%iterations - 1, &
                 "and every iteration after it is an update")
      if (allocated(error)) return
      ! The point of the pair: the setting arrived and did something, and the
      ! something is invisible in every number the SCF reports as a result.
      call check(error, scf%full_fock_builds < scf%iterations, &
                 "incremental building must skip full builds, or it did not arrive")
   end subroutine test_on

   subroutine test_same_energy(error)
      !! And it is exact -- the whole reason the count was needed
      !!
      !! An incremental build is `G(D) = G(D_ref) + G(D - D_ref)`, which is an
      !! identity rather than an approximation, so the two paths must agree far
      !! beyond the convergence threshold. They differ only in how much
      !! screening discards, which is why the counts differ and the energies
      !! must not.
      type(error_type), allocatable, intent(out) :: error

      type(rhf_result_t) :: on, off

      call converge(on, error, incremental=.true.)
      if (allocated(error)) return
      call converge(off, error, incremental=.false.)
      if (allocated(error)) return

      call check(error, abs(on%energy - off%energy) < 1.0e-10_dp, &
                 "incremental building must not move the energy")
      if (allocated(error)) return
      call check(error, on%iterations == off%iterations, &
                 "nor the iteration count")
   end subroutine test_same_energy

end module test_mqc_incremental_fock

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_incremental_fock, only: collect_incremental_fock
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_incremental_fock", collect_incremental_fock)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
