!! Stand-in for dftd4 when the build has no D4 library
module mqc_dispersion_d4
   !! The fpm build, and any CMake build with MQC_ENABLE_DFTD4=OFF, compiles
   !! this. It names the option and computes nothing, rather than failing to
   !! link or -- far worse -- returning zero and letting a deck that asked for
   !! dispersion quietly get a bare Kohn-Sham energy.
   !!
   !! The real implementation is `backends/dftd4/mqc_dispersion_d4.f90`: same
   !! module name, same two public procedures, exactly one of the pair ever
   !! compiled. It lives outside `src/` because fpm globs `src/` and cannot
   !! fetch or link dftd4, so fpm gets this one; and this one lives here, beside
   !! the DFT it belongs to, rather than in `src/methods/stubs`, because nothing
   !! in its interface mentions a backend type.
   !!
   !! It carries `charge` like its twin does, unused. A stub whose signature
   !! differed from the real one would compile on whichever build was tried
   !! first and fail on the other.
   use pic_types, only: dp, default_int
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: dispersion_d4_available
   public :: dispersion_d4_correction

contains

   pure function dispersion_d4_available() result(available)
      !! .false. -- this build has no dftd4
      !!
      !! Asked through `mqc_dispersion_apply`, so a deck naming
      !! `keywords.dft.dispersion: "d4"` on a build without the library is
      !! refused while it is still a deck, and by the version banner, so
      !! `run_validation.py` skips the gated cases rather than failing them.
      logical :: available

      available = .false.
   end function dispersion_d4_available

   subroutine dispersion_d4_correction(kind, functional, charge, atomic_numbers, coordinates, &
                                       energy, gradient, error)
      !! No-op stand-in: report the missing library, compute nothing
      character(len=*), intent(in) :: kind
      character(len=*), intent(in) :: functional
      real(dp), intent(in) :: charge
      integer(default_int), intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), intent(out) :: energy
      real(dp), intent(out), optional :: gradient(:, :)
      type(error_t), intent(out) :: error

      energy = 0.0_dp
      if (present(gradient)) gradient = 0.0_dp
      call error%set(ERROR_VALIDATION, "an empirical dispersion correction ('"// &
                     trim(adjustl(kind))//"', for functional '"//trim(adjustl(functional))// &
                     "') was asked for, and this build has no D4 dispersion library. "// &
                     "Build with CMake and -DMQC_ENABLE_DFTD4=ON.")
   end subroutine dispersion_d4_correction

end module mqc_dispersion_d4
