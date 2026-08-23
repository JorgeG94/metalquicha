!! How this binary was built -- the fallback for builds that are not CMake's
module mqc_build_info
   !! **There are two of this module and only ever one in a build.**
   !!
   !! CMake generates the real one from `cmake/templates/mqc_build_info.f90.in`
   !! into its build tree, filled with the values the configuration actually
   !! used, and names that file explicitly in `target_sources`. It never
   !! compiles this one, because `src/core/CMakeLists.txt` lists its sources by
   !! hand and this is not among them.
   !!
   !! fpm has no configure step to fill anything in, and discovers everything
   !! under `src/` by convention, so this is what it gets. The values say they
   !! are unknown rather than guessing: an fpm build genuinely does not know
   !! which BLAS it was handed, and a benchmark recorded against "Release" that
   !! was not one is worse than a benchmark recorded against "unrecorded".
   !!
   !! Adding this file to `src/core/CMakeLists.txt` would break the CMake build
   !! with a duplicate module, which is the failure mode to expect if the two
   !! ever drift.
   implicit none
   private

   public :: MQC_BUILD_TYPE, MQC_BLAS_VENDOR, MQC_ARCH_FLAGS
   public :: MQC_FORTRAN_COMPILER, MQC_FORTRAN_FLAGS, MQC_OPENMP
   public :: build_info_line

   character(len=*), parameter :: MQC_BUILD_TYPE = "unrecorded (fpm)"
   character(len=*), parameter :: MQC_BLAS_VENDOR = "unrecorded (fpm)"
   character(len=*), parameter :: MQC_ARCH_FLAGS = "unrecorded (fpm)"
   character(len=*), parameter :: MQC_FORTRAN_COMPILER = "unrecorded (fpm)"
   character(len=*), parameter :: MQC_FORTRAN_FLAGS = "unrecorded (fpm)"
   character(len=*), parameter :: MQC_OPENMP = "unrecorded (fpm)"

contains

   pure function build_info_line() result(line)
      !! One line, in the same `key=value` shape `--version` already uses for
      !! features, so whatever parses that can parse this.
      character(len=:), allocatable :: line

      line = "build: type="//trim(MQC_BUILD_TYPE)// &
             " compiler="//trim(MQC_FORTRAN_COMPILER)// &
             " blas="//trim(MQC_BLAS_VENDOR)// &
             " openmp="//trim(MQC_OPENMP)// &
             " arch="//trim(MQC_ARCH_FLAGS)
   end function build_info_line

end module mqc_build_info
