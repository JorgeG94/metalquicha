# cuEST GPU backend, kept out of src/ so the fpm build never sees it -- fpm
# globs src/ and cannot link the cuEST library, so it builds the CPU/xTB path
# only. This directory is pulled in by the top-level CMakeLists when
# MQC_ENABLE_CUEST is ON.

# --- bindings: vendored copy or fetched from upstream mod_cuest --------------
#
# The sources are identical either way; the vendored copy is a straight copy.
# Fetching needs network at configure time, so configure on a login node.
if(MQC_CUEST_BINDINGS STREQUAL "fetch")
  include(FetchContent)
  FetchContent_Declare(
    mod_cuest
    GIT_REPOSITORY https://github.com/JorgeG94/mod_cuest.git
    GIT_TAG ${MQC_CUEST_BINDINGS_TAG}
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(mod_cuest)
  message(
    STATUS "cuEST bindings: fetched from mod_cuest@${MQC_CUEST_BINDINGS_TAG}")
  # mod_cuest has no CMakeLists of its own -- it is a source drop, not a library
  # -- so the files are added directly.
  target_sources(
    ${main_lib}
    PRIVATE ${mod_cuest_SOURCE_DIR}/cuest.f90
            ${mod_cuest_SOURCE_DIR}/cuest_helpers.f90
            ${mod_cuest_SOURCE_DIR}/cudafort/cuda_runtime.F90
            ${mod_cuest_SOURCE_DIR}/cudafort/cuda_helpers.f90
            # Ours, not mod_cuest's: hand-written, covering the handful of BLAS
            # calls the SCF needs and the one dense eigensolver. Upstream has no
            # equivalent, so this path has to take them from the vendored
            # directory whichever way the generated bindings arrive -- without
            # them the build stops at `Cannot open module file 'cublas.mod'`,
            # which is what it did.
            ${CMAKE_CURRENT_SOURCE_DIR}/bindings/cublas.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/bindings/cusolver.f90)
  # `.F90` upstream as well as here, so the capital F asks for the preprocessing
  # the CUDA 13 guards need on both paths. The property is set for the same
  # belt-and-braces reason as on the vendored copy.
  set_source_files_properties(${mod_cuest_SOURCE_DIR}/cudafort/cuda_runtime.F90
                              PROPERTIES Fortran_PREPROCESS ON)
elseif(MQC_CUEST_BINDINGS STREQUAL "vendored")
  message(STATUS "cuEST bindings: vendored under backends/cuest/bindings")
  add_subdirectory(bindings)
else()
  message(
    FATAL_ERROR
      "MQC_CUEST_BINDINGS must be 'vendored' or 'fetch', got '${MQC_CUEST_BINDINGS}'"
  )
endif()

add_subdirectory(backend)
