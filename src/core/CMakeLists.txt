target_sources(
  ${main_lib}
  PRIVATE # Core modules
          mqc_calculation_defaults.f90
          mqc_calculation_keywords.f90
          mqc_physical_constants.F90
          mqc_program_limits.f90
          mqc_version.f90
          mqc_elements.f90
          mqc_nuclear_repulsion.f90
          mqc_atomic_radii.f90
          mqc_atomic_guess_common.f90
          mqc_pcm_radii.f90
          mqc_geometry.f90
          mqc_logo.f90
          mqc_acknowledgements.F90
          mqc_mpi_comms.f90
          mqc_resources.f90
          mqc_quao_rows.f90
          mqc_result_types.f90
          mqc_fingerprint.f90
          mqc_json_output_types.f90)

# What this binary was built with, recorded where CMake still knows it.
#
# The BLAS choice alone moves a density functional run by a factor of five and
# cannot be recovered from inside the program, so a benchmark that does not
# report it is not reproducible. Emitted into the build tree rather than the
# source tree: it is a product of the configuration, and two build directories
# off one checkout legitimately disagree. `mqc_build_info.f90` beside this file
# is the fpm build's stand-in and is not compiled here.
set(MQC_BUILD_INFO_BLAS "${MQC_BLA_VENDOR}")
if(NOT MQC_BUILD_INFO_BLAS)
  set(MQC_BUILD_INFO_BLAS "cmake-default")
endif()
if(OpenMP_Fortran_FOUND)
  set(MQC_BUILD_INFO_OPENMP "yes")
else()
  set(MQC_BUILD_INFO_OPENMP "no")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" _mqc_cfg)
set(MQC_BUILD_INFO_FFLAGS
    "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_mqc_cfg}}")
string(STRIP "${MQC_BUILD_INFO_FFLAGS}" MQC_BUILD_INFO_FFLAGS)
# The architecture flags in effect, read from the flags themselves rather than
# from MQC_ARCH_FLAGS, whose default is the sentinel "DEFAULT" and says nothing
# about whether -march=native was applied.
string(REGEX MATCHALL "(-march=|-mtune=|-mno-|-x|-mcpu=)[^ ]+" _mqc_arch
             "${MQC_BUILD_INFO_FFLAGS}")
list(REMOVE_DUPLICATES _mqc_arch)
string(JOIN " " MQC_BUILD_INFO_ARCH ${_mqc_arch})
if(NOT MQC_BUILD_INFO_ARCH)
  set(MQC_BUILD_INFO_ARCH "none")
endif()
configure_file("${PROJECT_SOURCE_DIR}/cmake/templates/mqc_build_info.f90.in"
               "${PROJECT_BINARY_DIR}/generated/mqc_build_info.f90" @ONLY)
target_sources(${main_lib}
               PRIVATE "${PROJECT_BINARY_DIR}/generated/mqc_build_info.f90")
