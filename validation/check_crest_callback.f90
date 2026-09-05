!! Manual check that mqc can serve CREST its gradients in memory
!!
!!     cmake -B build -DMQC_ENABLE_CREST=ON
!!     cmake --build build --target check_crest_callback
!!     ./build/check_crest_callback
!!
!! This is the seam the whole CREST integration rests on, and it is worth a
!! check of its own because everything about it is quiet when it goes wrong.
!! CREST holds a procedure pointer; this program installs one and asks CREST
!! for an energy and gradient. Nothing is written to disk, no process is
!! spawned, and if the pointer were never called the answer would come back
!! zero rather than as an error.
!!
!! The potential is harmonic -- E = K/2 * sum(x^2), dE/dx = K*x -- for the same
!! reason the HDF5 layout check uses index-encoding values: an analytic result
!! that is wrong in an obvious way if the geometry arrives transposed, scaled
!! or not at all. A real ab initio gradient would prove less, not more, because
!! nothing here could tell a wrong one from a right one.
module crest_callback_probe
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_adapter, only: driver_config_t
   use mqc_resources, only: resources_t
   use mqc_calculation_interface, only: compute_energy_and_forces
   implicit none
   private

   public :: harmonic_engrad, K_FORCE
   public :: install_mqc_engrad, mqc_engrad

   real(wp), parameter :: K_FORCE = 0.25_wp

   !! Context for the mqc-backed callback below.
   !!
   !! Held as module state for the same reason mqc_geometry_optimizer holds
   !! ctx_config and ctx_resources that way: CREST's hook takes a geometry and
   !! nothing else, so everything a calculation needs besides the geometry has
   !! to be waiting for it. A driver installs this once and CREST then calls
   !! the hook as many times as its sampling wants.
   type(system_geometry_t), save :: ctx_geom
   type(driver_config_t), save :: ctx_config
   type(resources_t), save :: ctx_resources
   logical, save :: ctx_ready = .false.

contains

   subroutine harmonic_engrad(nat, at, xyz, energy, grad, iostatus)
      !! Stands in for mqc's own energy and gradient
      !!
      !! The signature is CREST's `engrad_callback` abstract interface: plain
      !! arrays, coordinates in bohr, gradient in Eh/bohr. Deliberately not
      !! mqc's `compute_energy_and_forces` yet -- that one needs a
      !! `system_geometry_t`, a `driver_config_t` and a communicator, and
      !! whether the pointer is reached at all should be established before any
      !! of that is in the way.
      integer, intent(in) :: nat
      integer, intent(in) :: at(nat)
      real(wp), intent(in) :: xyz(3, nat)
      real(wp), intent(out) :: energy
      real(wp), intent(out) :: grad(3, nat)
      integer, intent(out) :: iostatus

      iostatus = 0
      energy = 0.0_wp
      grad = 0.0_wp
      if (nat < 1 .or. at(1) < 1) then
         iostatus = 1
         return
      end if

      energy = 0.5_wp*K_FORCE*sum(xyz*xyz)
      grad = K_FORCE*xyz
   end subroutine harmonic_engrad

   subroutine install_mqc_engrad(geom, config, res)
      type(system_geometry_t), intent(in) :: geom
      type(driver_config_t), intent(in) :: config
      type(resources_t), intent(in) :: res

      ctx_geom = geom
      ctx_config = config
      ctx_resources = res
      ctx_ready = .true.
   end subroutine install_mqc_engrad

   subroutine mqc_engrad(nat, at, xyz, energy, grad, iostatus)
      !! CREST's engrad_callback, answered by mqc's own SCF
      !!
      !! Only the coordinates move. Everything else -- the method, the basis,
      !! the charge, the communicator -- was fixed when the driver installed
      !! this, which is exactly the contract a sampling loop wants: it varies a
      !! geometry and nothing else.
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
         ! CREST may not change the atom count under a sampling run, and if it
         ! ever does the installed context is the wrong one rather than a
         ! resizable one. Refuse instead of reallocating into a mismatch.
         iostatus = 2
         return
      end if
      if (any(at /= ctx_geom%element_numbers)) then
         iostatus = 3
         return
      end if

      ctx_geom%coordinates = xyz
      call compute_energy_and_forces(ctx_geom, ctx_config, ctx_resources, &
                                     energy, gradient=grad)
   end subroutine mqc_engrad

end module crest_callback_probe

program check_crest_callback
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use crest_calculator, only: calcdata, calculation_settings, jobtype, engrad
   use strucrd, only: coord
   use crest_callback_probe, only: harmonic_engrad, K_FORCE, &
                                   install_mqc_engrad, mqc_engrad
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_adapter, only: driver_config_t
   use mqc_resources, only: resources_t
   use mqc_calculation_interface, only: compute_energy_and_forces
   use mqc_calc_types, only: CALC_TYPE_GRADIENT
   use mqc_method_types, only: METHOD_TYPE_HF
   use pic_mpi_lib, only: pic_mpi_init, pic_mpi_finalize, comm_world
   implicit none

   integer, parameter :: NAT = 3
   real(wp), parameter :: CALLBACK_TOL = 1.0e-12_wp

   type(calcdata) :: calc, calc2
   type(calculation_settings) :: sett, sett2
   type(coord) :: mol
   real(wp) :: energy, e_ref, gmax
   real(wp), allocatable :: grad(:, :), g_ref(:, :)
   integer :: io, failures

   type(system_geometry_t) :: geom
   type(driver_config_t) :: config
   type(resources_t) :: res
   real(wp) :: e_direct, e_crest
   real(wp), allocatable :: g_direct(:, :), g_crest(:, :)
   integer :: nranks

   failures = 0

   ! a geometry with no symmetry at all, so a transposed or partially
   ! copied array cannot accidentally give the right answer
   mol%nat = NAT
   allocate (mol%at(NAT), mol%xyz(3, NAT))
   mol%at = [8, 1, 1]
   mol%xyz(:, 1) = [0.10_wp, 0.20_wp, 0.30_wp]
   mol%xyz(:, 2) = [1.40_wp, -0.50_wp, 0.70_wp]
   mol%xyz(:, 3) = [-1.10_wp, 0.90_wp, -1.30_wp]

   allocate (grad(3, NAT), g_ref(3, NAT))
   e_ref = 0.5_wp*K_FORCE*sum(mol%xyz*mol%xyz)
   g_ref = K_FORCE*mol%xyz

   sett%id = jobtype%callback
   sett%external_engrad => harmonic_engrad
   call calc%add(sett)

   call engrad(mol, calc, energy, grad, io)

   write (*, "(a)") "  mqc -> CREST -> mqc, in memory"
   write (*, "(a,i0)") "    iostatus            ", io
   write (*, "(a,f20.12)") "    energy from CREST   ", energy
   write (*, "(a,f20.12)") "    analytic reference  ", e_ref
   gmax = maxval(abs(grad - g_ref))
   write (*, "(a,es20.8)") "    max gradient error  ", gmax

   if (io /= 0) then
      write (*, "(a)") "  FAIL: CREST reported a non-zero status"
      failures = failures + 1
   end if

   ! Exact, not approximate. Both sides evaluate the same expression on the
   ! same numbers; the only thing between them is the pointer call, so any
   ! difference at all is a transport fault rather than arithmetic.
   if (energy /= e_ref) then
      write (*, "(a)") "  FAIL: energy did not survive the callback"
      failures = failures + 1
   end if
   if (any(grad /= g_ref)) then
      write (*, "(a)") "  FAIL: gradient did not survive the callback"
      failures = failures + 1
   end if

   ! and the reference must not be zero, or the two comparisons above would
   ! pass without the callback ever having been reached
   if (e_ref == 0.0_wp) then
      write (*, "(a)") "  FAIL: the reference is zero; this check proves nothing"
      failures = failures + 1
   end if

   ! ------------------------------------------------------------------
   ! and now the same path with mqc's own SCF behind the pointer
   ! ------------------------------------------------------------------

   ! The communicator has to be built the way app/main.f90 builds it.
   ! resources_init only sets thread and GPU defaults -- it does not touch
   ! MPI -- and compute_energy_and_forces returns its results only under
   ! `world_comm%rank() == 0`, so an uninitialised comm makes the whole
   ! calculation run and hand back nothing. It ran, printed a converged SCF,
   ! and returned zero. Whatever drives CREST has to own this.
   call pic_mpi_init()
   res%mpi_comms%world_comm = comm_world()
   res%mpi_comms%node_comm = res%mpi_comms%world_comm%split()
   call res%init()
   nranks = res%mpi_comms%world_comm%size()
   write (*, "(a)") ""
   write (*, "(a,i0)") "  ranks in this run     ", nranks

   ! Conformer sampling runs on one rank: CREST loops on the master with
   ! OpenMP underneath, while compute_energy_and_forces expects every rank to
   ! call it. Refusing here is the guard the scoping settled on -- at run time
   ! rather than at build time, so one binary serves both and `mpirun -n 1`
   ! is still allowed.
   if (nranks > 1) then
      write (*, "(a)") "  SKIP: CREST sampling runs on a single rank; re-run without mpirun"
   else
      geom%n_monomers = 1
      geom%atoms_per_monomer = NAT
      geom%total_atoms = NAT
      geom%charge = 0
      geom%multiplicity = 1
      allocate (geom%element_numbers(NAT), geom%coordinates(3, NAT))
      geom%element_numbers = mol%at
      geom%coordinates = mol%xyz

      config%calc_type = CALC_TYPE_GRADIENT
      config%nlevel = 0
      config%method_config%method_type = METHOD_TYPE_HF
      config%method_config%basis_set = "sto-3g"

      allocate (g_direct(3, NAT), g_crest(3, NAT))

      ! mqc on its own
      call compute_energy_and_forces(geom, config, res, e_direct, gradient=g_direct)

      ! the same calculation, reached through CREST's pointer
      call install_mqc_engrad(geom, config, res)
      sett2%id = jobtype%callback
      sett2%external_engrad => mqc_engrad
      call calc2%add(sett2)
      call engrad(mol, calc2, e_crest, g_crest, io)

      write (*, "(a,f20.12)") "    E, mqc direct       ", e_direct
      write (*, "(a,f20.12)") "    E, through CREST    ", e_crest
      write (*, "(a,es20.8)") "    max gradient error  ", maxval(abs(g_crest - g_direct))

      if (io /= 0) then
         write (*, "(a,i0)") "  FAIL: CREST reported status ", io
         failures = failures + 1
      end if
      ! mqc against mqc, so the callback is the only thing in between. It says
      ! nothing about whether the gradient is right -- only that it arrives
      ! unaltered, which is what this checks. To rounding and not to the bit:
      ! the Fock and gradient builds merge per-thread partial sums in the
      ! order the threads finish, so two runs on more than one thread differ
      ! in the last place (2e-15 on two threads, 7e-15 on sixteen). A
      ! transposed, scaled or missing gradient is off by 1e-3 or more.
      if (abs(e_crest - e_direct) > CALLBACK_TOL) then
         write (*, "(a)") "  FAIL: the SCF energy changed crossing the callback"
         failures = failures + 1
      end if
      if (any(abs(g_crest - g_direct) > CALLBACK_TOL)) then
         write (*, "(a)") "  FAIL: the gradient changed crossing the callback"
         failures = failures + 1
      end if
      if (e_direct == 0.0_wp) then
         write (*, "(a)") "  FAIL: mqc returned zero; the comparison is vacuous"
         failures = failures + 1
      end if
   end if

   call pic_mpi_finalize()

   if (failures == 0) then
      write (*, "(a)") ""
      write (*, "(a)") "  OK: CREST called into this program and took its numbers"
   else
      write (*, "(a,i0,a)") "  ", failures, " failure(s)"
      error stop 1
   end if
end program check_crest_callback
