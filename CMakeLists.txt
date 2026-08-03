cmake_minimum_required(VERSION 3.22)
cmake_policy(SET CMP0144 NEW)

project(
  "metalquicha"
  LANGUAGES Fortran C
  VERSION 0.1.0
  DESCRIPTION "A hastily put together fragmentation quantum chemistry program")

set(project_name metalquicha)
set(main_lib ${project_name})
set(exe_name mqc)
set(all_targets ${main_lib} ${exe_name})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# see cmake/ subdirectory for this
add_subdirectory(cmake)
include(demo_status)

option(MQC_ENABLE_TESTING "Enable testing for ${project_name}" ON)
option(MQC_ENABLE_TBLITE "Link to tblite for xTB calculations" ON)
option(MQC_ENABLE_CUEST "Link to NVIDIA cuEST for GPU integrals (HF/DFT)" OFF)
set(CUEST_ROOT
    ""
    CACHE PATH "Root of the unpacked libcuest-<ver>-archive (needs lib/)")
set(MQC_CUEST_BINDINGS
    "vendored"
    CACHE STRING "Where the cuEST Fortran bindings come from: vendored | fetch")
set_property(CACHE MQC_CUEST_BINDINGS PROPERTY STRINGS vendored fetch)
set(MQC_CUEST_BINDINGS_TAG
    "main"
    CACHE STRING "Git tag/branch of mod_cuest when MQC_CUEST_BINDINGS=fetch")

add_library(${main_lib} STATIC)
set_target_properties(${main_lib} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set(libraries_to_link)

if(MQC_ENABLE_TESTING)
  enable_testing()
  if(NOT TARGET test-drive::test-drive)
    find_package("test-drive" REQUIRED)
  endif()
endif()

if(NOT TARGET pic::pic)
  find_package("pic" REQUIRED)
endif()

if(NOT TARGET pic-mpi::pic-mpi)
  find_package("pic-mpi" REQUIRED)
endif()

if(NOT TARGET pic-blas::pic-blas)
  find_package("pic-blas" REQUIRED)
endif()

if(MQC_ENABLE_TBLITE)
  if(NOT TARGET "tblite::tblite")
    find_package("tblite" REQUIRED)
  endif()
else()
  target_compile_definitions(${main_lib} PRIVATE MQC_WITHOUT_TBLITE)
endif()

# NVIDIA cuEST -- GPU integrals behind the HF/DFT methods.
#
# Only the shared library is needed to build: the Fortran bindings under
# src/cuest_bindings are pre-generated and vendored, so the cuEST headers are
# not required here. Linking libcuest.so by full path makes CMake set the rpath
# for us.
#
# Runtime note: cuEST ships GPU code for sm_80 and newer only. It compiles and
# links anywhere, but on an older card (Volta/V100, e.g. the Gadi login node)
# cuestCreate returns CUEST_STATUS_UNSUPPORTED_ARCHITECTURE.
if(MQC_ENABLE_CUEST)
  if(NOT CUEST_ROOT)
    message(
      FATAL_ERROR "MQC_ENABLE_CUEST=ON requires -DCUEST_ROOT=<cuest archive>")
  endif()
  find_library(
    CUEST_LIBRARY
    NAMES cuest
    PATHS "${CUEST_ROOT}/lib"
    NO_DEFAULT_PATH REQUIRED)
  find_package(CUDAToolkit REQUIRED)
  list(APPEND libraries_to_link ${CUEST_LIBRARY} CUDA::cudart)
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_CUEST)
  message(STATUS "cuEST enabled: ${CUEST_LIBRARY}")
endif()

if(NOT TARGET "jsonfortran::jsonfortran")
  find_package("jsonfortran" REQUIRED)
  list(APPEND libraries_to_link jsonfortran::jsonfortran)
endif()

find_package(
  OpenMP
  COMPONENTS Fortran
  REQUIRED)

add_subdirectory(src)

# The cuEST backend lives outside src/ on purpose: the fpm build globs all of
# src/, and fpm cannot link the cuEST binary library, so fpm builds the CPU/xTB
# path only. CMake adds the backend explicitly here.
if(MQC_ENABLE_CUEST)
  add_subdirectory(cuest_backend)
endif()

if(MQC_ENABLE_TESTING)
  add_subdirectory(test)
endif()

set_target_properties(${main_lib} PROPERTIES Fortran_MODULE_DIRECTORY
                                             ${CMAKE_BINARY_DIR}/modules)

list(APPEND libraries_to_link pic::pic pic-mpi::pic-mpi pic-blas::pic-blas
     OpenMP::OpenMP_Fortran)

if(MQC_ENABLE_TBLITE)
  list(APPEND libraries_to_link tblite::tblite)
endif()

target_link_libraries(${main_lib} PUBLIC ${libraries_to_link})

# Add json-fortran module directory for Fortran USE statements and ensure
# jsonfortran is built before metalquicha (for .mod files)
target_include_directories(${main_lib} PRIVATE "${jsonfortran_BINARY_DIR}")
add_dependencies(${main_lib} jsonfortran)

add_executable(${exe_name} ${PROJECT_SOURCE_DIR}/app/main.f90)

# Finite-difference check on the analytic gradient. Not a ctest case: it needs a
# GPU of sm_80 or newer, so it cannot run in CI or on a login node.
if(MQC_ENABLE_CUEST)
  add_executable(check_gradient
                 ${PROJECT_SOURCE_DIR}/validation/check_gradient.f90)
  target_include_directories(check_gradient
                             PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_gradient PRIVATE ${main_lib})
endif()

target_include_directories(${exe_name} PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(${exe_name} PRIVATE ${main_lib})

foreach(tgt ${all_targets})
  message(STATUS "${tgt}")
  install(
    TARGETS ${tgt}
    EXPORT ${project_name}Targets
    ARCHIVE DESTINATION lib
    INCLUDES
    DESTINATION include/${project_name})
endforeach()

install(
  DIRECTORY ${CMAKE_BINARY_DIR}/modules/
  DESTINATION include/${project_name}
  FILES_MATCHING
  PATTERN "*.mod")

install(
  EXPORT ${project_name}Targets
  FILE ${project_name}Targets.cmake
  NAMESPACE ${project_name}::
  DESTINATION lib/cmake/${project_name})

include(CMakePackageConfigHelpers)

# RENAME YOUR sampleConfig.cmake.in to match ${project_name}Config.cmake
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${project_name}Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${project_name}Config.cmake"
  INSTALL_DESTINATION lib/cmake/${project_name})

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${project_name}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${project_name}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${project_name}ConfigVersion.cmake"
        DESTINATION lib/cmake/${project_name})

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/${project_name}Config.cmake"
  DESTINATION lib/cmake/${project_name}
  RENAME ${project_name}-config.cmake)

print_summary()

if(CMAKE_BUILD_TYPE STREQUAL "Coverage-mqc")
  add_custom_target(
    coverage
    COMMAND ${CMAKE_CTEST_COMMAND} -R "mqc" --output-on-failure
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/validation python3
            run_validation.py
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/validation python3
            run_validation.py --mpi --np 2
    COMMAND ${CMAKE_COMMAND} -E remove -f coverage.info
    COMMAND lcov --directory . --capture --output-file coverage.info
    COMMAND lcov --ignore-errors unused --remove coverage.info '/usr/*'
            '*/build/*' --output-file coverage_filtered.info
    COMMAND genhtml coverage_filtered.info --output-directory coverage_report
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating code coverage report")
endif()
