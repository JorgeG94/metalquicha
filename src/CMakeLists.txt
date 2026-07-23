add_subdirectory(basis)
add_subdirectory(cli)
add_subdirectory(core)
add_subdirectory(fragmentation)
add_subdirectory(interface)
add_subdirectory(io)
add_subdirectory(methods)
add_subdirectory(parallel)
add_subdirectory(utils)
add_subdirectory(vibrational)
# The cuEST/CUDA bindings must compile before anything that USEs them. They can
# come from the vendored copy under cuest_bindings/, or be fetched from the
# upstream mod_cuest repository -- the sources are identical either way, since
# the vendored copy is a straight copy of them.
#
# Fetching needs network access at configure time. Gadi compute nodes have
# none, so configure on the login node (which is where you build anyway).
if(MQC_ENABLE_CUEST)
  if(MQC_CUEST_BINDINGS STREQUAL "fetch")
    include(FetchContent)
    FetchContent_Declare(
      mod_cuest
      GIT_REPOSITORY https://github.com/JorgeG94/mod_cuest.git
      GIT_TAG ${MQC_CUEST_BINDINGS_TAG}
      GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(mod_cuest)
    message(STATUS "cuEST bindings: fetched from mod_cuest@${MQC_CUEST_BINDINGS_TAG}")
    # mod_cuest has no CMakeLists of its own at the root -- it is a source
    # drop, not a library -- so the files are added directly.
    target_sources(
      ${main_lib}
      PRIVATE ${mod_cuest_SOURCE_DIR}/cuest.f90
              ${mod_cuest_SOURCE_DIR}/cuest_helpers.f90
              ${mod_cuest_SOURCE_DIR}/cudafort/cuda_runtime.f90
              ${mod_cuest_SOURCE_DIR}/cudafort/cuda_helpers.f90)
  elseif(MQC_CUEST_BINDINGS STREQUAL "vendored")
    message(STATUS "cuEST bindings: vendored under src/cuest_bindings")
    add_subdirectory(cuest_bindings)
  else()
    message(FATAL_ERROR
            "MQC_CUEST_BINDINGS must be 'vendored' or 'fetch', got '${MQC_CUEST_BINDINGS}'")
  endif()
endif()

target_sources(${main_lib} PRIVATE mqc_calc_types.f90 mqc_method_types.f90
                                   mqc_driver.f90)
