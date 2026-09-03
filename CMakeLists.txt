cmake_minimum_required(VERSION 3.22)
cmake_policy(SET CMP0144 NEW)

project(
  "metalquicha"
  LANGUAGES Fortran C
  VERSION 0.2.0
  DESCRIPTION "A hastily put together fragmentation quantum chemistry program")

set(project_name metalquicha)
set(main_lib ${project_name})
set(exe_name mqc)
set(all_targets ${main_lib} ${exe_name})

# Set before the dependencies are declared, so they are built PIC too: the
# ctypes shared library links every one of them into a .so, and a dependency
# built without it fails the link outright rather than degrading.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Sequential MKL, not threaded, when MKL is what will be found.
#
# Every path in this program that computes a fragment pins itself to one OpenMP
# thread -- the serial and distributed MBE processors, the GMBE processor, the
# unfragmented driver. Threaded MKL underneath that is thread creation for work
# that will never use it, and on the coordinator it is threads competing with
# the ones handing out fragments. Parallelism here is MPI over fragments, not
# BLAS over one small matrix.
#
# This is NOT a fix for the tblite corruption above >1 OpenMP thread. That was
# measured against MKL_NUM_THREADS=1 and MKL_THREADING_LAYER=SEQUENTIAL and
# fails identically, so it is tblite's own OpenMP and not the BLAS. The pins
# remain load-bearing.
#
# Defaulted rather than forced: MKLROOT is the reliable signal that oneAPI is in
# the environment, and where it is not -- CI uses OpenBLAS -- CMake's own
# detection is left alone. Override with -DMQC_BLA_VENDOR=... , or set it empty
# to get the previous behaviour back.
if(DEFINED ENV{MKLROOT} AND NOT DEFINED MQC_BLA_VENDOR)
  set(MQC_BLA_VENDOR
      "Intel10_64lp_seq"
      CACHE STRING "BLA_VENDOR passed to find_package(BLAS/LAPACK)")
else()
  set(MQC_BLA_VENDOR
      ""
      CACHE STRING "BLA_VENDOR passed to find_package(BLAS/LAPACK)")
endif()

if(MQC_BLA_VENDOR)
  # Read by FindBLAS/FindLAPACK, which run inside the dependencies, so this has
  # to be set before they are declared.
  set(BLA_VENDOR ${MQC_BLA_VENDOR})
  message(STATUS "BLAS vendor: ${BLA_VENDOR}")
endif()

# Threaded MKL under gfortran must find *gfortran's* libgomp.
#
# `FindBLAS` completes the GNU threading layer with a bare `find_library(gomp)`,
# and on a box with a stray `/usr/local/lib/libgomp.so` -- which on at least one
# of ours is a symlink to LLVM's `libomp.so` -- that is what it finds. The link
# then records a DT_NEEDED on LLVM's runtime while `-fopenmp` records one on
# GNU's, and both load. Two runtimes, two thread pools, two notions of thread
# identity, and threads from one writing into state owned by the other: wrong
# answers above one thread, different every run, with nothing in the build log
# to say so.
#
# Asking the compiler where its own libgomp is and searching there first costs
# nothing and cannot pick the wrong one. It matters only for the threaded Intel
# layers -- the sequential one links no OpenMP at all, which is why the default
# above has never tripped over this.
if(BLA_VENDOR MATCHES "^Intel"
   AND NOT BLA_VENDOR MATCHES "_seq$"
   AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  execute_process(
    COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgomp.so
    OUTPUT_VARIABLE mqc_gomp
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(IS_ABSOLUTE "${mqc_gomp}" AND EXISTS "${mqc_gomp}")
    get_filename_component(mqc_gomp_dir "${mqc_gomp}" DIRECTORY)
    list(PREPEND CMAKE_LIBRARY_PATH "${mqc_gomp_dir}")
    message(STATUS "Threaded MKL: OpenMP runtime pinned to ${mqc_gomp}")
  else()
    message(
      WARNING
        "Threaded MKL with gfortran, but gfortran's own libgomp could not be "
        "located. If a foreign libgomp is on the link path this build will "
        "load two OpenMP runtimes and give wrong answers above one thread. "
        "Check `ldd` on the binary for libomp.so.")
  endif()
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(GNUInstallDirs)

# cmake/ holds this project's own modules, cmake/modules/ the Find modules for
# its dependencies. Both are on the path so `include(Mqc...)` and
# `find_package(...)` each resolve without a path spelled out at the call.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
     "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# see cmake/ subdirectory for this
add_subdirectory(cmake)
include(demo_status)

# ---------------------------------------------------------------------------
# What this file is, and what it is not.
#
# Four includes and the targets they surround. Everything that used to be
# spelled out here has moved to the file whose subject it is, so that a change
# lands in one place and a reader looking for one thing does not scroll past the
# other three:
#
# cmake/MqcOptions.cmake        every build option, with its default
# cmake/MqcValidation.cmake     how a validation program is registered
# cmake/MqcDependencies.cmake   resolving each dependency and wiring it in
# cmake/MqcInstall.cmake        install, export, and the coverage target
#
# Dependencies themselves are declared one per file under cmake/modules/, all of
# them through `mqc_fetch`; nothing is fetched inline any more.
# ---------------------------------------------------------------------------
include(MqcOptions)
include(ExtractBasisSets)
mqc_extract_basis_sets()

include(MqcValidation)

add_library(${main_lib} STATIC)
set_target_properties(${main_lib} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set(libraries_to_link)

include(MqcDependencies)

add_executable(${exe_name} ${PROJECT_SOURCE_DIR}/app/main.f90)
target_include_directories(${exe_name} PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(${exe_name} PRIVATE ${main_lib})

# The backend-independent validation programs. The ones that need CPU integrals
# are registered by backends/cenzontle/validation.cmake instead, which the
# dependency resolution above has already reached if that backend is on.
add_subdirectory(validation)

# Shared library for ctypes. The C entry points live in the static archive and
# nothing in this translation unit calls them, so the linker would drop every
# one; the whole archive has to be forced in.
add_library(mqc_shared SHARED
            ${PROJECT_SOURCE_DIR}/src/interface/mqc_capi_export.f90)
set_target_properties(
  mqc_shared PROPERTIES OUTPUT_NAME mqc Fortran_MODULE_DIRECTORY
                                        "${CMAKE_BINARY_DIR}/modules_shared")
# Forcing the whole archive in is spelled differently per linker: GNU ld (Linux)
# brackets the archive with --whole-archive/--no-whole-archive, while Apple's ld
# takes the archive path directly via -force_load and has no closing flag.
# Linking ${main_lib} again in both cases pulls in its transitive dependencies
# (BLAS, MPI, tblite, ...); the force-load only covers the archive itself.
if(APPLE)
  target_link_libraries(
    mqc_shared PRIVATE "-Wl,-force_load,$<TARGET_FILE:${main_lib}>" ${main_lib})
else()
  target_link_libraries(
    mqc_shared PRIVATE "-Wl,--whole-archive" "$<TARGET_FILE:${main_lib}>"
                       "-Wl,--no-whole-archive" ${main_lib})
endif()
add_dependencies(mqc_shared ${main_lib})

include(MqcInstall)
