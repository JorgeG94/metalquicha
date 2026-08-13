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

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# see cmake/ subdirectory for this
add_subdirectory(cmake)
include(demo_status)

option(MQC_ENABLE_TESTING "Enable testing for ${project_name}" ON)
option(MQC_ENABLE_TBLITE "Link to tblite for xTB calculations" ON)
option(
  MQC_ENABLE_HDF5
  "Link to HDF5 for binary checkpoints (needed for gradient/Hessian restart)"
  OFF)
# On by default, and the pair of them is what makes this program able to do ab
# initio quantum chemistry on a machine without a GPU. They were opt-in, which
# had two costs worth naming: `cmake -B build` produced a program that could run
# xTB and nothing else, and libcint-on/libxc-off -- the natural thing to type if
# you only wanted integrals -- produced a build where every deck naming a
# functional is refused. Neither is a good default.
#
# The price is that a default configure reaches the network: libxc is a tarball
# and libcint a clone. For a machine without one, point
# FETCHCONTENT_SOURCE_DIR_LIBCINT and FETCHCONTENT_SOURCE_DIR_LIBXC at local
# copies, or turn these off and build the xTB path alone.
option(MQC_ENABLE_LIBCINT
       "Fetch libcint for CPU Gaussian integrals (real ab initio without a GPU)"
       ON)
option(MQC_ENABLE_LIBXC
       "Fetch libxc for exchange-correlation functionals (CPU DFT)" ON)

# A coverage build measures this project's own Fortran, not its dependencies.
# libxc will not even configure under a build type it does not recognise: its
# autocmake safeguard allows only Debug/Release/RelWithDebInfo and aborts on
# "Coverage-mqc", so a coverage configure that fetched it dies before a line of
# mqc is compiled. Drop both CPU backends for coverage -- which is also how
# coverage ran before they became default-on. A normal (non-cache) set so it
# shadows the option here and in the subdirectories without rewriting the cache,
# leaving a later reconfigure of the same tree under another build type intact.
if(CMAKE_BUILD_TYPE STREQUAL "Coverage-mqc")
  set(MQC_ENABLE_LIBXC OFF)
  set(MQC_ENABLE_LIBCINT OFF)
  message(STATUS "Coverage-mqc build: libxc and libcint disabled")
endif()

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

# Basis sets come out of the Basis Set Exchange bundle at configure time. Add a
# name here (lowercase, as BSE spells it) to have it extracted; the bundle has
# every basis BSE publishes, so nothing needs downloading. basis_sets/*.json is
# derived entirely from this list -- anything else found there is deleted -- so
# removing a name here removes it from the build.
set(MQC_BASIS_BUNDLE
    "basis_sets-json-0.12.tar.bz2"
    CACHE STRING "Basis Set Exchange bundle under basis_sets/")
set(MQC_BASIS_SETS
    # Pople, Dunning, and the auxiliaries the validation inputs name
    "sto-3g;3-21g;6-31g;6-31g_st_;6-31g_st__st_;cc-pvdz;cc-pvtz;aug-cc-pvdz;aug-cc-pvtz;cc-pvtz-jkfit;cc-pvdz-rifit;cc-pvtz-rifit"
    # every def2 the bundle carries: orbital, diffuse, RI/JK auxiliaries, ECP
    "def2-sv(p);def2-sv(p)-jkfit;def2-sv(p)-rifit;def2-svp;def2-svp-rifit;def2-svpd;def2-svpd-rifit"
    "def2-tzvp;def2-tzvp-rifit;def2-tzvpd;def2-tzvpd-rifit;def2-tzvpp;def2-tzvpp-rifit;def2-tzvppd;def2-tzvppd-rifit"
    "def2-qzvp;def2-qzvp-rifit;def2-qzvpd;def2-qzvpp;def2-qzvpp-rifit;def2-qzvppd;def2-qzvppd-rifit"
    "def2-mtzvp;def2-mtzvpp;def2-mtzvpp-rij;def2-universal-jfit;def2-universal-jkfit;def2-ecp"
    CACHE STRING "Basis sets to extract from ${MQC_BASIS_BUNDLE}")

include(ExtractBasisSets)
mqc_extract_basis_sets()

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

# tblite does not build with every Fortran compiler, and it is on by default, so
# an unsupported one fails deep inside somebody else's source rather than here.
# Refused at configure time instead, naming the flag that gets you the rest of
# the program: the ab initio path builds fine on all of these.
if(MQC_ENABLE_TBLITE)
  if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Intel|IntelLLVM)$")
    message(
      FATAL_ERROR
        "MQC_ENABLE_TBLITE=ON needs gfortran, ifort or ifx; this is "
        "${CMAKE_Fortran_COMPILER_ID}. tblite does not build with it. Configure "
        "with -DMQC_ENABLE_TBLITE=OFF to build everything else -- the CPU ab "
        "initio path (libcint/libxc) and the GPU path both work here.")
  endif()
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
  # cuBLAS backs the SCF's own linear algebra (the Fock assembly, the DIIS
  # commutator and the orthogonal-basis transform) and cuSOLVER the Fock
  # diagonalization, both reached through the hand-written bindings in
  # backends/cuest/bindings/. cuEST links them internally too, but that is its
  # private dependency -- ours has to be requested explicitly.
  list(APPEND libraries_to_link ${CUEST_LIBRARY} CUDA::cudart CUDA::cublas
       CUDA::cusolver)
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
  add_subdirectory(backends/cuest)
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

# libxc, for exchange-correlation functionals.
#
# FetchContent rather than ExternalProject, which is what the sample beside this
# file used. The difference is when the subproject is configured:
# ExternalProject builds at *build* time, so its targets do not exist while this
# project is being configured and a caller has to name the .a files by hand,
# guess the module directory and add_dependencies to force an order.
# FetchContent configures it now, so `xc` and `xcf03` are ordinary targets --
# their include directories, their link order and their build order all come
# along on their own.
#
# The tarball rather than the git repository: libxc's history carries the
# generated functional sources and is far larger than any release. URL_HASH pins
# it, which a tag would not -- a tag can be moved and a hash cannot.
add_executable(check_grid ${PROJECT_SOURCE_DIR}/validation/check_grid.f90)
target_include_directories(check_grid PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(check_grid PRIVATE ${main_lib})

add_executable(check_partition
               ${PROJECT_SOURCE_DIR}/validation/check_partition.f90)
target_include_directories(check_partition
                           PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(check_partition PRIVATE ${main_lib})

add_executable(check_radial ${PROJECT_SOURCE_DIR}/validation/check_radial.f90)
target_include_directories(check_radial PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(check_radial PRIVATE ${main_lib})

add_executable(check_lebedev ${PROJECT_SOURCE_DIR}/validation/check_lebedev.f90)
target_include_directories(check_lebedev PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(check_lebedev PRIVATE ${main_lib})

if(MQC_ENABLE_LIBXC)
  include(FetchContent)
  # Static, so nothing has to be found on the library path at run time; the
  # global POSITION_INDEPENDENT_CODE above keeps it usable from libmqc.so.
  set(BUILD_SHARED_LIBS
      OFF
      CACHE BOOL "" FORCE)
  set(ENABLE_FORTRAN
      ON
      CACHE BOOL "" FORCE)
  set(BUILD_TESTING
      OFF
      CACHE BOOL "" FORCE)
  set(ENABLE_PYTHON
      OFF
      CACHE BOOL "" FORCE)
  FetchContent_Declare(
    libxc
    URL https://gitlab.com/libxc/libxc/-/archive/7.1.2/libxc-7.1.2.tar.gz
    URL_HASH
      SHA256=c517ce61820ea8114664a4280b6a6bc74a4f22f1fd1ea4ddecd6df0caeeae4f4
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(libxc)
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_LIBXC)
  # xcf03 first, then xc: the Fortran interface depends on the C library, and a
  # static link resolves left to right.
  target_link_libraries(${main_lib} PRIVATE $<BUILD_INTERFACE:xcf03>
                                            $<BUILD_INTERFACE:xc>)
  target_include_directories(${main_lib} PRIVATE ${libxc_BINARY_DIR})
  message(STATUS "libxc enabled: exchange-correlation functionals available")
endif()

# libcint, for Gaussian integrals on the CPU.
#
# The point is not speed. Everything ab initio here runs through cuEST, which
# needs an A100 -- so the HF and DFT paths are compile-checked and never
# executed anywhere the developer is sitting, which is how they got this far
# unverified. A CPU backend makes them runnable on a laptop, and gives cuEST a
# second implementation to disagree with, which is the check it has never had.
#
# JorgeG94's fork rather than upstream sunqm: it adds `libcint_fortran`, a
# Fortran API over the C one, so this program needs no hand-written bindings at
# all -- including libcint_gto_norm, which is the easy thing to get wrong
# (libcint wants contraction coefficients pre-multiplied by the primitive
# normalisation, and an overlap diagonal that is not 1 is the symptom).
#
# libcint rather than libint, measured rather than assumed: clean clone to
# libcint.so in 9 seconds against 46 MB of pre-generated source and a long
# build. For a backend whose job is to be correct and always available, build
# time is the property that matters; if throughput ever does, libint can sit
# behind the same interface later.
if(MQC_ENABLE_LIBCINT)
  include(FetchContent)
  FetchContent_Declare(
    libcint
    GIT_REPOSITORY https://github.com/JorgeG94/libcint.git
    # A commit, not a branch. libcint develops on `develop`, so tracking that
    # would change what this program builds against on every push there, with
    # nothing here recording that it changed. This hash is reachable from
    # develop and carries libcint_3c2e_sph and libcint_2c2e_sph, without which
    # the density-fitting path will not compile -- and their Cartesian twins
    # libcint_3c2e_cart, libcint_2c2e_cart and libcint_tot_cgto_cart, without
    # which a Cartesian basis such as 6-31G* can reach the exact four-index path
    # but not the density-fitted one. It also carries the three- and two-centre
    # optimizers: libcint exports them, the Fortran interface had no way to make
    # one, so every fitting integral rebuilt its shell-pair data per call. Worth
    # replacing with a tag once the fork has one -- a tag says what it is, a
    # hash does not, and a tag also lets GIT_SHALLOW come back, which a bare
    # hash cannot resolve.
    GIT_TAG 525cab503f6134502cf69de42234b5897276dbf8)
  # Offline builds: point FETCHCONTENT_SOURCE_DIR_LIBCINT at a local clone and
  # no network is touched. Compute nodes on most clusters have none, and a
  # configure step that reaches for GitHub fails there in a way that reads like
  # a compiler problem.
  FetchContent_MakeAvailable(libcint)
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_LIBCINT)
  # BUILD_INTERFACE, not just PRIVATE: a static library records even its private
  # dependencies in its interface as $<LINK_ONLY:...>, so CMake wants them in
  # this project's export set -- which a fetched third-party library is never
  # going to join. Scoping them to the build keeps the integrals backend
  # internal, which it is.
  target_link_libraries(${main_lib} PRIVATE $<BUILD_INTERFACE:cint_fortran>
                                            $<BUILD_INTERFACE:cint>)
  add_subdirectory(backends/libcint)
  message(
    STATUS "libcint enabled: CPU Gaussian integrals via the Fortran interface")
endif()

# HDF5, for checkpoints that carry derivatives. Only the C library is needed:
# the bindings are hand-written against the C ABI precisely so that no
# compiler-specific hdf5.mod has to exist. Optional, because a laptop build
# doing energy-only screening has no use for it.
if(MQC_ENABLE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C)
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_HDF5)
  target_include_directories(${main_lib} PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${main_lib} PUBLIC ${HDF5_C_LIBRARIES})
  add_subdirectory(backends/hdf5)
  message(STATUS "HDF5 checkpoints enabled: ${HDF5_C_LIBRARIES}")
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

# Manual check that a session splits rank 0 from the workers. Not a ctest case:
# it is only meaningful under mpirun with more than one rank, and ctest here has
# no MPI harness.
add_executable(check_run ${PROJECT_SOURCE_DIR}/validation/check_run.f90)
target_include_directories(check_run PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(check_run PRIVATE ${main_lib})

if(MQC_ENABLE_HDF5)
  add_executable(check_hdf5 ${PROJECT_SOURCE_DIR}/validation/check_hdf5.f90)
  target_include_directories(check_hdf5 PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_hdf5 PRIVATE ${main_lib})

  add_executable(check_hdf5_checkpoint
                 ${PROJECT_SOURCE_DIR}/validation/check_hdf5_checkpoint.f90)
  target_include_directories(check_hdf5_checkpoint
                             PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_hdf5_checkpoint PRIVATE ${main_lib})
endif()

# Needs both: libcint for the density, libxc for the functional.
if(MQC_ENABLE_LIBXC AND MQC_ENABLE_LIBCINT)
  add_executable(check_dft ${PROJECT_SOURCE_DIR}/validation/check_dft.f90)
  target_include_directories(check_dft PRIVATE "${CMAKE_BINARY_DIR}/modules"
                                               ${libxc_BINARY_DIR})
  target_link_libraries(check_dft PRIVATE ${main_lib} cint_fortran cint xcf03
                                          xc)
endif()

if(MQC_ENABLE_LIBXC)
  add_executable(check_libxc ${PROJECT_SOURCE_DIR}/validation/check_libxc.f90)
  target_include_directories(check_libxc PRIVATE "${CMAKE_BINARY_DIR}/modules"
                                                 ${libxc_BINARY_DIR})
  target_link_libraries(check_libxc PRIVATE ${main_lib} xcf03 xc)
endif()

if(MQC_ENABLE_LIBCINT)
  add_executable(bench_df_k ${PROJECT_SOURCE_DIR}/validation/bench_df_k.f90)
  target_include_directories(bench_df_k PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(bench_df_k PRIVATE ${main_lib} cint_fortran cint)

  add_executable(check_direct ${PROJECT_SOURCE_DIR}/validation/check_direct.f90)
  target_include_directories(check_direct PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_direct PRIVATE ${main_lib} cint_fortran cint)

  add_executable(check_libcint
                 ${PROJECT_SOURCE_DIR}/validation/check_libcint.f90)
  target_include_directories(check_libcint
                             PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_libcint PRIVATE ${main_lib} cint_fortran cint)

  add_executable(check_df ${PROJECT_SOURCE_DIR}/validation/check_df.f90)
  target_include_directories(check_df PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_df PRIVATE ${main_lib} cint_fortran cint)

  add_executable(check_rhf ${PROJECT_SOURCE_DIR}/validation/check_rhf.f90)
  target_include_directories(check_rhf PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_rhf PRIVATE ${main_lib} cint_fortran cint)

  add_executable(check_cc ${PROJECT_SOURCE_DIR}/validation/check_cc.f90)
  target_include_directories(check_cc PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_cc PRIVATE ${main_lib} cint_fortran cint)

  add_executable(check_ao ${PROJECT_SOURCE_DIR}/validation/check_ao.f90)
  target_include_directories(check_ao PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_ao PRIVATE ${main_lib} cint_fortran cint)

  # Multipole integrals against PySCF's, and the first dipole the CPU path can
  # produce. Milestone 1 of the EFP work.
  add_executable(check_multipole
                 ${PROJECT_SOURCE_DIR}/validation/check_multipole.f90)
  target_include_directories(check_multipole
                             PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_multipole PRIVATE ${main_lib} cint_fortran cint)

  # Boys localization: the invariants first, then against pyscf.lo.Boys.
  add_executable(check_localize
                 ${PROJECT_SOURCE_DIR}/validation/check_localize.f90)
  target_include_directories(check_localize
                             PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_localize PRIVATE ${main_lib} cint_fortran cint)

  # CPHF and the static polarizability against PySCF's own response solver. The
  # finite-field check lives in the unit tests, needing no reference.
  add_executable(check_cphf ${PROJECT_SOURCE_DIR}/validation/check_cphf.f90)
  target_include_directories(check_cphf PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_cphf PRIVATE ${main_lib} cint_fortran cint)

  # The screening fit itself, in Fortran: grid, objective and search.
  add_executable(check_screening_fit
                 ${PROJECT_SOURCE_DIR}/validation/check_screening_fit.f90)
  target_include_directories(check_screening_fit
                             PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_screening_fit PRIVATE ${main_lib} cint_fortran
                                                    cint)

  # Electrostatic potential integrals at arbitrary points, against PySCF. The
  # enabling piece for a charge-penetration screening fit.
  add_executable(check_esp ${PROJECT_SOURCE_DIR}/validation/check_esp.f90)
  target_include_directories(check_esp PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_esp PRIVATE ${main_lib} cint_fortran cint)

  # The twelve imaginary frequencies, their static limit, and GAMESS's own
  # dynamic polarizable points.
  add_executable(check_dynamic
                 ${PROJECT_SOURCE_DIR}/validation/check_dynamic.f90)
  target_include_directories(check_dynamic
                             PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_dynamic PRIVATE ${main_lib} cint_fortran cint)

  # The projection data a potential carries: the LMO Fock matrix, which needs no
  # basis ordering, against GAMESS's own.
  add_executable(check_projection
                 ${PROJECT_SOURCE_DIR}/validation/check_projection.f90)
  target_include_directories(check_projection
                             PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_projection PRIVATE ${main_lib} cint_fortran cint)

  # Distributed multipoles: the sum rule in Fortran, needing no reference, then
  # the distribution itself against GAMESS in Python.
  add_executable(check_dma ${PROJECT_SOURCE_DIR}/validation/check_dma.f90)
  target_include_directories(check_dma PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_dma PRIVATE ${main_lib} cint_fortran cint)

  # Distributed LMO polarizabilities against a GAMESS MAKEFP potential. The
  # first milestone whose reference is GAMESS rather than PySCF.
  add_executable(
    check_distributed_polarizability
    ${PROJECT_SOURCE_DIR}/validation/check_distributed_polarizability.f90)
  target_include_directories(check_distributed_polarizability
                             PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(check_distributed_polarizability
                        PRIVATE ${main_lib} cint_fortran cint)
endif()

add_executable(check_session_run
               ${PROJECT_SOURCE_DIR}/validation/check_session_run.f90)
target_include_directories(check_session_run
                           PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(check_session_run PRIVATE ${main_lib})

add_executable(check_bcast_system
               ${PROJECT_SOURCE_DIR}/validation/check_bcast_system.f90)
target_include_directories(check_bcast_system
                           PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(check_bcast_system PRIVATE ${main_lib})

add_executable(check_bcast ${PROJECT_SOURCE_DIR}/validation/check_bcast.f90)
target_include_directories(check_bcast PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(check_bcast PRIVATE ${main_lib})

add_executable(check_session ${PROJECT_SOURCE_DIR}/validation/check_session.f90)
target_include_directories(check_session PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(check_session PRIVATE ${main_lib})

target_include_directories(${exe_name} PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(${exe_name} PRIVATE ${main_lib})

# Shared library for ctypes. The C entry points live in the static archive and
# nothing in it references them -- they are called from outside -- so a normal
# archive link would take none of them and produce an empty .so.
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
            '*/build/*' '*/validation/*' --output-file coverage_filtered.info
    COMMAND genhtml coverage_filtered.info --output-directory coverage_report
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating code coverage report")
endif()
