!! Stand-in for DL-FIND when the build has no libdlfind
module mqc_dlfind_bridge
   !! Same name and same entry points as the real bridge, declining.
   !!
   !! A stub rather than a preprocessor guard at the call site, matching
   !! mqc_libcint_bridge_stub and mqc_cuest_bridge_stub: the optimizer reads
   !! the same either way, and `dlfind_available` is how it asks in advance so
   !! the refusal can name the build option.
   use pic_types, only: dp
   use mqc_optimizer_types, only: optimizer_settings_t, energy_gradient_i, step_callback_i, &
                                  hessian_i
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: dlfind_available
   public :: dlfind_optimize
   public :: dlfind_connected_minima

contains

   pure function dlfind_available() result(available)
      !! Whether this build can optimize a geometry
      logical :: available
      available = .false.
   end function dlfind_available

   subroutine dlfind_connected_minima(a, energy_first, b, energy_second, found, e_saddle)
      !! Nothing ran, so nothing was connected
      real(dp), allocatable, intent(out) :: a(:, :), b(:, :)
      real(dp), intent(out) :: energy_first, energy_second
      logical, intent(out) :: found
      real(dp), intent(out) :: e_saddle

      found = .false.
      e_saddle = 0.0_dp
      energy_first = 0.0_dp
      energy_second = 0.0_dp
      if (.false.) then
         allocate (a(0, 0), b(0, 0))
      end if
   end subroutine dlfind_connected_minima

   subroutine dlfind_optimize(opt_settings, natoms, znuc, residues, coords, &
                              energy_gradient, step_taken, final_energy, error, &
                              hessian, endpoint)
      !! No-op stand-in: optimizing a geometry needs the DL-FIND backend
      !!
      !! The argument list has to track the real bridge's exactly, including the
      !! arguments nothing here reads: the caller is compiled once against
      !! whichever of the two is present, so a signature that has drifted is a
      !! compile error in a build nobody makes locally.
      type(optimizer_settings_t), intent(in) :: opt_settings
      integer, intent(in) :: natoms
      integer, intent(in) :: znuc(natoms)
      integer, intent(in) :: residues(natoms)
      real(dp), intent(inout) :: coords(3, natoms)
      procedure(energy_gradient_i) :: energy_gradient
      procedure(step_callback_i) :: step_taken
      real(dp), intent(out) :: final_energy
      type(error_t), intent(inout) :: error
      procedure(hessian_i), optional :: hessian
      real(dp), intent(in), optional :: endpoint(:, :)

      final_energy = 0.0_dp
      call error%set(ERROR_VALIDATION, &
                     "This build has no geometry optimizer. Configure with "// &
                     "-DMQC_ENABLE_DLFIND=ON to fetch and link DL-FIND.")

   end subroutine dlfind_optimize

end module mqc_dlfind_bridge
