!! The many-body expansion contexts, and what they refuse to run without
module test_mqc_mbe_contexts
   !! Each expansion -- MBE, GMBE, FMO -- is a context object that a caller
   !! fills in before asking it to run: a geometry, a partition, and for the
   !! distributed path a set of MPI resources. Every one of those is optional
   !! at the type level and mandatory in fact, so the run entry points check
   !! them and refuse.
   !!
   !! **A refusal is not a number, so no validation deck can express it.** A
   !! deck that forgot its geometry is not a deck; the case being guarded
   !! against is a *caller* -- the driver, the C interface, a future
   !! fragmentation scheme -- assembling a context and missing a field, which
   !! is why the checks exist and why nothing was exercising them. What each
   !! test asserts is that the refusal happens and that the object survives it
   !! well enough to be destroyed, since a run that stops halfway through
   !! allocating is the other way this goes wrong.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, int32
   use mqc_many_body_expansion, only: fmo_context_t
   use mqc_method_config, only: method_config_t
   use mqc_calc_types, only: CALC_TYPE_ENERGY
   use mqc_json_output_types, only: json_output_data_t
   implicit none
   private

   public :: collect_mqc_mbe_contexts_tests

contains

   subroutine collect_mqc_mbe_contexts_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("a_fresh_context_has_neither_geometry_nor_mpi", test_fresh), &
                  new_unittest("serial_run_refuses_without_a_geometry", test_serial_refuses), &
                  new_unittest("distributed_run_refuses_without_resources", test_distributed_refuses), &
                  new_unittest("init_carries_the_method_settings", test_init_settings) &
                  ]
   end subroutine collect_mqc_mbe_contexts_tests

   subroutine fresh_context(context)
      type(fmo_context_t), intent(out) :: context
      type(method_config_t) :: config

      config%basis_set = "6-31g"
      call context%init(config, int(CALC_TYPE_ENERGY, int32))
   end subroutine fresh_context

   subroutine test_fresh(error)
      type(error_type), allocatable, intent(out) :: error

      type(fmo_context_t) :: context

      call fresh_context(context)
      ! Both are pointers or allocatables that a caller fills in afterwards, and
      ! the two queries are how every run entry point decides whether it may
      ! proceed. A fresh context answering yes to either would let a run start
      ! on a geometry that is not there.
      call check(error,.not. context%has_geometry(), &
                 "a fresh context claims to have a geometry")
      if (allocated(error)) return
      call check(error,.not. context%has_mpi(), &
                 "a fresh context claims to have MPI resources")
      if (allocated(error)) return
      call context%destroy()
   end subroutine test_fresh

   subroutine test_serial_refuses(error)
      type(error_type), allocatable, intent(out) :: error

      type(fmo_context_t) :: context
      type(json_output_data_t) :: json_data

      call fresh_context(context)
      ! No geometry and no partition: the run must come back having done
      ! nothing rather than dereferencing what it was not given. The energy
      ! staying at its initial zero is what says it did not proceed -- an FMO
      ! that ran would have put something there.
      call context%run_serial(json_data)
      call check(error, abs(context%energy) < 1.0e-14_dp, &
                 "a context with no geometry produced an energy")
      if (allocated(error)) return
      call check(error,.not. json_data%has_energy, &
                 "a refused run still wrote an energy into the output document")
      if (allocated(error)) return
      call context%destroy()
   end subroutine test_serial_refuses

   subroutine test_distributed_refuses(error)
      type(error_type), allocatable, intent(out) :: error

      type(fmo_context_t) :: context
      type(json_output_data_t) :: json_data

      call fresh_context(context)
      ! The distributed path needs MPI resources as well, and it is the one a
      ! caller reaches by mistake: a serial build calling the distributed entry
      ! point gets a null pointer, not a slower calculation.
      call context%run_distributed(json_data)
      call check(error, abs(context%energy) < 1.0e-14_dp, &
                 "a context with no resources produced an energy")
      if (allocated(error)) return
      call check(error,.not. json_data%has_energy, &
                 "a refused distributed run wrote an energy anyway")
      if (allocated(error)) return
      call context%destroy()
   end subroutine test_distributed_refuses

   subroutine test_init_settings(error)
      !! **The context carries the basis twice, and `init` fills only one.**
      !!
      !! `init` copies the whole `method_config`, whose `basis_set` is where a
      !! deck's basis lands -- and leaves `context%basis`, which is what the FMO
      !! run actually reads, at the type's default of 6-31g. The driver bridges
      !! them by hand, one line after the `init` call, and any other caller
      !! that forgets gets an FMO in 6-31g that converges and reports a number.
      !!
      !! This test does not assert that they agree, because they do not. It
      !! pins the shape as it is, so that a change making `init` carry the
      !! basis through breaks here and gets noticed, rather than silently
      !! double-setting a field the driver is already setting.
      type(error_type), allocatable, intent(out) :: error

      type(fmo_context_t) :: context
      type(method_config_t) :: config

      config%basis_set = "cc-pvdz"
      call context%init(config, int(CALC_TYPE_ENERGY, int32))

      call check(error, trim(context%method_config%basis_set) == "cc-pvdz", &
                 "init did not copy the method configuration")
      if (allocated(error)) return
      call check(error, trim(context%basis) == "6-31g", &
                 "init now sets the context basis; the driver's copy of that "// &
                 "line is redundant and this test should say so")
      if (allocated(error)) return

      ! `init` is `intent(out)`, so everything it does not touch is the type's
      ! own default: FMO2 rather than a monomer-only expansion, and the exact
      ! embedding rather than point charges.
      call check(error, context%level == 2, "the default expansion level is not FMO2")
      if (allocated(error)) return
      call check(error, trim(context%esp) == "exact", &
                 "the default embedding is not the exact one")
      if (allocated(error)) return
      call context%destroy()
   end subroutine test_init_settings

end module test_mqc_mbe_contexts

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mbe_contexts, only: collect_mqc_mbe_contexts_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_mbe_contexts", collect_mqc_mbe_contexts_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
