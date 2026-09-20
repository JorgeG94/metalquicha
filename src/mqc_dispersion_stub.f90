!! Stand-in for s-dftd3 when the build has no dispersion library
module mqc_dispersion
   !! The fpm build, and any CMake build with MQC_ENABLE_DFTD3=OFF, compiles
   !! this. It names the option and computes nothing, rather than failing to
   !! link or -- far worse -- returning zero and letting a deck that asked for
   !! dispersion quietly get a bare Kohn-Sham energy.
   !!
   !! The real implementation is `backends/dftd3/mqc_dispersion.f90`: same
   !! module name, same two public procedures, exactly one of the pair ever
   !! compiled. It lives outside `src/` because fpm globs `src/` and cannot
   !! fetch or link s-dftd3, so fpm gets this one; and this one lives here,
   !! beside the DFT it belongs to, rather than in `src/methods/stubs`, because
   !! nothing in its interface mentions a backend type.
   use pic_types, only: dp, default_int
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: dispersion_available
   public :: dispersion_correction

contains

   pure function dispersion_available() result(available)
      !! .false. -- this build has no s-dftd3
      !!
      !! Asked by `read_dft_dispersion`, so a deck naming
      !! `keywords.dft.dispersion` on a build without the library is refused
      !! while it is still a deck, and by the version banner, so
      !! `run_validation.py` skips the gated cases rather than failing them.
      logical :: available

      available = .false.
   end function dispersion_available

   subroutine dispersion_correction(kind, functional, atomic_numbers, coordinates, &
                                    energy, gradient, error)
      !! No-op stand-in: report the missing library, compute nothing
      character(len=*), intent(in) :: kind
      character(len=*), intent(in) :: functional
      integer(default_int), intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), intent(out) :: energy
      real(dp), intent(out), optional :: gradient(:, :)
      type(error_t), intent(out) :: error

      energy = 0.0_dp
      if (present(gradient)) gradient = 0.0_dp
      call error%set(ERROR_VALIDATION, "an empirical dispersion correction ('"// &
                     trim(adjustl(kind))//"', for functional '"//trim(adjustl(functional))// &
                     "') was asked for, and this build has no dispersion library. "// &
                     "Build with CMake and -DMQC_ENABLE_DFTD3=ON.")
   end subroutine dispersion_correction

end module mqc_dispersion
