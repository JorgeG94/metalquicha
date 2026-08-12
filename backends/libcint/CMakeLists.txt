# CPU Gaussian-integral backend, kept out of src/ so the fpm build never sees it
# -- the same reason backends/cuest/ and backends/hdf5/ live here.
#
# These modules `use libcint_fortran`, which only exists once CMake has fetched
# libcint. fpm globs src/ and has no such dependency declared, so leaving them
# there fails the fpm build outright:
#
# <ERROR> *cmd_build* Target error: Unable to find source for module dependency:
# "libcint_fortran" used by "src/methods/mqc_libcint_integrals.f90"
#
# and there is no CMake option that helps, because fpm never reads CMake. So fpm
# builds the xTB path and CMake adds this directory when MQC_ENABLE_LIBCINT is
# ON.
target_sources(
  ${main_lib}
  PRIVATE mqc_libcint_atomic_guess.f90
          mqc_libcint_bridge.f90
          mqc_libcint_direct.f90
          mqc_libcint_integrals.f90
          mqc_libcint_rhf.f90
          mqc_libcint_mp2.f90
          mqc_libcint_cc.f90
          mqc_libcint_ao_data.f90
          mqc_libcint_ao.f90)
