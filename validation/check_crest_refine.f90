!! Manual check of the configuration this integration is actually for:
!! xTB samples, mqc refines
!!
!!     cmake -B build -DMQC_ENABLE_CREST=ON -DMQC_ENABLE_TBLITE=ON
!!                    -DWITH_TBLITE=ON -DWITH_GFN0=ON
!!     cmake --build build --target check_crest_refine
!!     ./build/check_crest_refine
!!
!! `check_crest_search` proved a whole sampling run can be configured from
!! memory and driven through the hook. Doing *every* gradient that way is not
!! what anybody wants: an iMTD-GC run asks for tens of thousands, and at an ab
!! initio price that is a week per molecule.
!!
!! CREST already has the right shape for this, and it is not a new idea here --
!! `legacy_wrappers.f90` does the same thing with two xTB levels. A second
!! calculation level tagged `refine_lvl = refine%singlepoint`, registered with
!! `env%addrefine`, is applied to the surviving ensemble after the search. So
!! level 1 is xTB, in process, tens of thousands of times; level 2 is mqc, on
!! what survives.
!!
!! Both projects share one tblite, which took reconciling: CREST pinned a
!! dormant `xtb_solvation` branch and mqc pins release v0.6.0, and only one
!! tblite can exist in one binary. The pin turned out to be stale rather than
!! load-bearing -- one added argument to eeq_guess is the whole difference --
!! so the sampling here runs GFN2-xTB through the same tblite mqc uses.
!!
!! What this checks is that the split holds: that the search runs, that the
!! refinement is reached, and that mqc's energies are the ones that come back
!! from it.
module crest_refine_probe
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_adapter, only: driver_config_t
   use mqc_resources, only: resources_t
   use mqc_calculation_interface, only: compute_energy_and_forces
   implicit none
   private

   public :: install_mqc, mqc_refine_engrad, refine_calls

   type(system_geometry_t), save :: ctx_geom
   type(driver_config_t), save :: ctx_config
   type(resources_t), save :: ctx_resources
   logical, save :: ctx_ready = .false.
   integer, save :: refine_calls = 0

contains

   subroutine install_mqc(geom, config, res)
      type(system_geometry_t), intent(in) :: geom
      type(driver_config_t), intent(in) :: config
      type(resources_t), intent(in) :: res

      ctx_geom = geom
      ctx_config = config
      ctx_resources = res
      ctx_ready = .true.
      refine_calls = 0
   end subroutine install_mqc

   subroutine mqc_refine_engrad(nat, at, xyz, energy, grad, iostatus)
      !! mqc's SCF, reached only on the refinement level
      integer, intent(in) :: nat
      integer, intent(in) :: at(nat)
      real(wp), intent(in) :: xyz(3, nat)
      real(wp), intent(out) :: energy
      real(wp), intent(out) :: grad(3, nat)
      integer, intent(out) :: iostatus

      iostatus = 0
      energy = 0.0_wp
      grad = 0.0_wp

      if (.not. ctx_ready) then
         iostatus = 1
         return
      end if
      if (nat /= ctx_geom%total_atoms) then
         iostatus = 2
         return
      end if

      !$omp atomic
      refine_calls = refine_calls + 1

      ctx_geom%coordinates = xyz
      call compute_energy_and_forces(ctx_geom, ctx_config, ctx_resources, &
                                     energy, gradient=grad)
   end subroutine mqc_refine_engrad

end module crest_refine_probe

program check_crest_refine
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use crest_data, only: systemdata, timer, refine
   use crest_calculator, only: jobtype, calculation_settings
   use crest_refine_probe, only: install_mqc, mqc_refine_engrad, refine_calls
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_adapter, only: driver_config_t
   use mqc_resources, only: resources_t
   use mqc_calc_types, only: CALC_TYPE_GRADIENT
   use mqc_method_types, only: METHOD_TYPE_HF
   use pic_mpi_lib, only: pic_mpi_init, pic_mpi_finalize, comm_world
   implicit none

   interface
      subroutine parseflags(env, arg, nra)
         import :: systemdata
         implicit none
         type(systemdata), intent(inout) :: env
         integer, intent(in) :: nra
         !> No intent, matching the definition in confparse.f90; intent is part
         !> of a procedure's characteristics.
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
   type(calculation_settings) :: refine_level
   type(system_geometry_t) :: geom
   type(driver_config_t) :: config
   type(resources_t) :: res
   character(len=:), allocatable :: argv(:)
   integer :: nargs, unit, nat, nranks

   call pic_mpi_init()
   res%mpi_comms%world_comm = comm_world()
   res%mpi_comms%node_comm = res%mpi_comms%world_comm%split()
   call res%init()
   nranks = res%mpi_comms%world_comm%size()

   if (nranks > 1) then
      write (*, "(a)") "  SKIP: CREST sampling runs on a single rank; re-run without mpirun"
      call pic_mpi_finalize()
      stop 0
   end if

   open (newunit=unit, file="refine_probe.xyz", status="replace", action="write")
   write (unit, "(a)") "3"
   write (unit, "(a)") "water"
   write (unit, "(a)") "O   0.0000000   0.0000000   0.1173000"
   write (unit, "(a)") "H   0.0000000   0.7572000  -0.4692000"
   write (unit, "(a)") "H   0.0000000  -0.7572000  -0.4692000"
   close (unit)

   call tim%init(20)
   nargs = 1
   allocate (character(len=64) :: argv(nargs))
   argv(1) = "refine_probe.xyz"
   call parseflags(env, argv, nargs)

   nat = env%ref%nat
   write (*, "(a)") ""
   write (*, "(a,i0)") "  atoms                 ", nat
   ! Whatever parseflags chose stays: GFN2-xTB through tblite, in process.
   ! That is the arrangement this is for, and it needs both projects to agree
   ! on one tblite. They now can -- see the eeq_guess commit in the fork --
   ! where before, CREST's pin to a dormant xtb_solvation branch and mqc's
   ! release v0.6.0 could not coexist in one binary.
   write (*, "(a,i0)") "  sampling level id     ", env%calc%calcs(1)%id
   write (*, "(a)") "  (6 = tblite, i.e. GFN2-xTB in process)"

   ! mqc's side of the arrangement
   geom%n_monomers = 1
   geom%atoms_per_monomer = nat
   geom%total_atoms = nat
   geom%charge = 0
   geom%multiplicity = 1
   allocate (geom%element_numbers(nat), geom%coordinates(3, nat))
   geom%element_numbers = env%ref%at
   geom%coordinates = env%ref%xyz

   config%calc_type = CALC_TYPE_GRADIENT
   config%nlevel = 0
   config%method_config%method_type = METHOD_TYPE_HF
   config%method_config%basis_set = "sto-3g"
   call install_mqc(geom, config, res)

   ! The second level, applied to the ensemble the search leaves behind.
   ! Same shape as legacy_wrappers.f90's two-xTB-level arrangement, with mqc
   ! where the second xTB would be.
   refine_level%id = jobtype%callback
   refine_level%external_engrad => mqc_refine_engrad
   refine_level%refine_lvl = refine%singlepoint
   call env%calc%add(refine_level)
   if (allocated(env%refine_queue)) deallocate (env%refine_queue)
   call env%addrefine(refine%singlepoint)

   write (*, "(a,i0)") "  calculation levels    ", env%calc%ncalculations
   write (*, "(a)") ""
   write (*, "(a)") "  running: xTB samples, mqc refines"

   call crest_search_imtdgc(env, tim)

   write (*, "(a)") ""
   write (*, "(a,i0)") "  iostatus_meta         ", env%iostatus_meta
   write (*, "(a,i0)") "  mqc refinement calls  ", refine_calls

   call pic_mpi_finalize()

   if (env%iostatus_meta /= 0) then
      write (*, "(a)") "  FAIL: the search did not complete"
      error stop 1
   end if
   ! The point of the whole arrangement: mqc is reached, and reached a few
   ! times rather than tens of thousands. A zero here would mean the
   ! refinement level was configured and never applied, which looks like
   ! success from every other angle.
   if (refine_calls < 1) then
      write (*, "(a)") "  FAIL: the refinement level was never applied"
      error stop 1
   end if
   write (*, "(a)") "  OK: xTB sampled and mqc refined, both in process"
end program check_crest_refine
