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
   implicit none
   private

   public :: harmonic_engrad, K_FORCE

   real(wp), parameter :: K_FORCE = 0.25_wp

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

end module crest_callback_probe

program check_crest_callback
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use crest_calculator, only: calcdata, calculation_settings, jobtype, engrad
   use strucrd, only: coord
   use crest_callback_probe, only: harmonic_engrad, K_FORCE
   implicit none

   integer, parameter :: NAT = 3

   type(calcdata) :: calc
   type(calculation_settings) :: sett
   type(coord) :: mol
   real(wp) :: energy, e_ref, gmax
   real(wp), allocatable :: grad(:, :), g_ref(:, :)
   integer :: io, failures

   failures = 0

   !> a geometry with no symmetry at all, so a transposed or partially
   !> copied array cannot accidentally give the right answer
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

   !> Exact, not approximate. Both sides evaluate the same expression on the
   !> same numbers; the only thing between them is the pointer call, so any
   !> difference at all is a transport fault rather than arithmetic.
   if (energy /= e_ref) then
      write (*, "(a)") "  FAIL: energy did not survive the callback"
      failures = failures + 1
   end if
   if (any(grad /= g_ref)) then
      write (*, "(a)") "  FAIL: gradient did not survive the callback"
      failures = failures + 1
   end if

   !> and the reference must not be zero, or the two comparisons above would
   !> pass without the callback ever having been reached
   if (e_ref == 0.0_wp) then
      write (*, "(a)") "  FAIL: the reference is zero; this check proves nothing"
      failures = failures + 1
   end if

   if (failures == 0) then
      write (*, "(a)") "  OK: CREST called into this program and took its numbers"
   else
      write (*, "(a,i0,a)") "  ", failures, " failure(s)"
      error stop 1
   end if
end program check_crest_callback
