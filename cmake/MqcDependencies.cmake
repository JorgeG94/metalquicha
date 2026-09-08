# Resolving every dependency, and wiring each one to ${main_lib}.
#
# The fetching itself is not here: each package is declared in its own
# cmake/modules/Find<pkg>.cmake, all of which call `mqc_fetch` and read the
# same. What is here is the part that is specific to *this* project -- which
# compile definitions a dependency turns on, which targets it adds to the link,
# and the order the answers have to be worked out in.
#
# Order is load-bearing in two places, both noted where they happen: BLA_VENDOR
# has to be settled before CREST is resolved, and libxc before libcint, because
# a validation program in the backend directory links both.
include_guard(GLOBAL)

# Dependencies do not bring their test suites with them.
#
# Set before the first fetch, because it is what each subproject's own
# `include(CTest)` / `option(BUILD_TESTING)` sees, and CTest defaults it ON.
# Left alone, `ctest` in this build tree runs 167 cases belonging to pic,
# tblite, mctc-lib, json-fortran, s-dftd3, test-drive and toml-f before it
# reaches one of ours -- which is why every instruction here says `ctest -R
# "mqc"`. toml-f's two have never passed, and a red square that is not ours is
# worse than no square.
#
# This project's own tests are unaffected: they go through `enable_testing()`
# and `add_test` directly and never consult BUILD_TESTING.
foreach(
  _mqc_dep_tests
  BUILD_TESTING # toml-f, and test-drive's fallback
  TESTDRIVE_BUILD_TESTING # test-drive
  PIC_ENABLE_TESTING # pic
  ENABLE_TESTING # pic-blas, pic-mpi
  WITH_TESTS # tblite, s-dftd3, CREST
  JSONFORTRAN_ENABLE_TESTS) # json-fortran
  set(${_mqc_dep_tests}
      OFF
      CACHE BOOL "Build a dependency's own test suite" FORCE)
endforeach()
unset(_mqc_dep_tests)

# mctc-lib and dftd4 are not on that list because neither can be told: both
# `add_subdirectory("test")` unconditionally, with no option in front of it.
# Their forty-odd cases are the remainder, and removing them means patching them
# upstream rather than configuring them away.

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

target_compile_definitions(${main_lib} PRIVATE CODATA_YEAR=${MQC_CODATA_YEAR})
message(STATUS "CODATA revision for the Bohr radius: ${MQC_CODATA_YEAR}")

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
  # CUDA 13 removed every _v2 entry point from the runtime, so the vendored
  # bindings have to know which major version they are being built against --
  # see the guard on cudaGetDeviceProperties in bindings/cuda_runtime.f90.
  #
  # 13.2 is the floor on that major version because it is the only CUDA 13 this
  # has been built against; 13.0 and 13.1 are refused rather than assumed to
  # behave like it. Nothing below 13 is gated -- 12.x is what the bindings were
  # generated from and is the well-trodden path.
  if(CUDAToolkit_VERSION_MAJOR EQUAL 13 AND CUDAToolkit_VERSION VERSION_LESS
                                            13.2)
    message(
      FATAL_ERROR
        "cuEST needs CUDA 12.x or CUDA >= 13.2, found ${CUDAToolkit_VERSION}")
  endif()
  target_compile_definitions(
    ${main_lib} PRIVATE MQC_CUDA_VERSION_MAJOR=${CUDAToolkit_VERSION_MAJOR})
  message(STATUS "CUDA toolkit: ${CUDAToolkit_VERSION}")
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

# LFortran: answer FindOpenMP's questions, because it does not know how to ask.
#
# CMake's FindOpenMP carries a table of per-compiler OpenMP flags and has no
# entry for LFortran, so the flag search comes back empty and the REQUIRED find
# below fails before anything is compiled. LFortran does support OpenMP -- the
# flag is `--openmp` -- so the fix is to fill in the row rather than to make
# OpenMP optional: `omp_lib` is used in twenty-odd modules here and the `!$omp`
# directives are not decoration.
#
# The runtime is the other half. `--openmp` makes LFortran's driver link
# `-lomp`, LLVM's runtime, which the conda toolchain does not ship; libgomp is
# what is there, and it resolves the same symbols. So the flag is attached to
# compilation only and the library is named explicitly, which is exactly the
# shape FindOpenMP already uses for every other compiler:
# `OpenMP_Fortran_LIB_NAMES` plus one `OpenMP_<name>_LIBRARY` per entry. Set
# them as cache variables and FindOpenMP skips its own detection and believes
# these.
#
# Override `OpenMP_Fortran_FLAGS` or `OpenMP_omp_LIBRARY` on the command line
# for an LFortran built against LLVM's runtime instead.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "LFortran" AND NOT MQC_ENABLE_SERIAL)
  if(NOT DEFINED OpenMP_Fortran_FLAGS)
    set(OpenMP_Fortran_FLAGS
        "--openmp"
        CACHE STRING "Fortran compiler flags for OpenMP parallelization")
  endif()
  find_library(
    OpenMP_gomp_LIBRARY
    NAMES gomp omp
    HINTS ENV CONDA_PREFIX
    PATH_SUFFIXES lib lib64
    DOC "OpenMP runtime linked into LFortran builds")
  if(NOT DEFINED OpenMP_Fortran_LIB_NAMES)
    set(OpenMP_Fortran_LIB_NAMES
        "gomp"
        CACHE STRING "Fortran compiler libraries for OpenMP parallelization")
  endif()
  # Told, not detected: the version-probe program FindOpenMP would compile to
  # decide these does not survive LFortran, and `use omp_lib` does work.
  set(OpenMP_Fortran_HAVE_OMPLIB_MODULE
      TRUE
      CACHE BOOL INTERNAL "")
  set(OpenMP_Fortran_HAVE_OMPLIB_HEADER
      FALSE
      CACHE BOOL INTERNAL "")
  mark_as_advanced(OpenMP_Fortran_FLAGS OpenMP_Fortran_LIB_NAMES
                   OpenMP_gomp_LIBRARY)
endif()

if(MQC_ENABLE_SERIAL)
  # No OpenMP at all -- not even the search, which is REQUIRED and would stop a
  # configure on a toolchain that has none. `compat/` supplies the `omp_lib` the
  # sources still `use`, so nothing under src/ or backends/ changes; see
  # compat/omp_lib_serial.f90.
  add_subdirectory(compat)
else()
  find_package(
    OpenMP
    COMPONENTS Fortran
    REQUIRED)
endif()

add_subdirectory(src)

# The cuEST backend lives outside src/ on purpose: the fpm build globs all of
# src/, and fpm cannot link the cuEST binary library, so fpm builds the CPU/xTB
# path only. CMake adds the backend explicitly here.
if(MQC_ENABLE_CUEST)
  add_subdirectory(backends/cuest)
endif()

# terco is built by nvfortran and reached through its C ABI, so mqc needs the
# shared library and terco's own interface declaration -- never a copy of it.
# Configure terco with TERCO_BUILD_SHARED=ON; a static libterco cannot be linked
# into this build because its device code was compiled by a different compiler
# than the one building mqc.
if(MQC_ENABLE_TERCO)
  if(NOT TERCO_ROOT)
    message(
      FATAL_ERROR
        "MQC_ENABLE_TERCO=ON requires -DTERCO_ROOT=<terco checkout or install>")
  endif()
  find_library(
    TERCO_LIBRARY
    NAMES terco
    PATHS "${TERCO_ROOT}/lib" "${TERCO_ROOT}/build-shared-gpu"
          "${TERCO_ROOT}/build-xc-gpu" "${TERCO_ROOT}"
    NO_DEFAULT_PATH REQUIRED)
  # A terco built in place rather than installed is the common case, so its
  # build directories are searched too. Naming the file outright with
  # -DTERCO_LIBRARY=<path> wins over all of them: find_library leaves a cache
  # entry that is already set alone. The authoritative declaration, compiled
  # from terco's tree. See backends/terco/CMakeLists.txt for why this is not
  # vendored. `include/terco` is where `make install` puts it and `include` is
  # where a source checkout keeps it, so both layouts are reachable from one
  # TERCO_ROOT.
  find_file(
    TERCO_INTERFACES
    NAMES trc_c_interfaces.f90
    PATHS "${TERCO_ROOT}/include/terco" "${TERCO_ROOT}/include"
    NO_DEFAULT_PATH REQUIRED)
  list(APPEND libraries_to_link ${TERCO_LIBRARY})
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_TERCO)
  message(STATUS "terco enabled: ${TERCO_LIBRARY}")
  add_subdirectory(backends/terco)
endif()

if(MQC_ENABLE_TESTING)
  add_subdirectory(test)
endif()

set_target_properties(${main_lib} PROPERTIES Fortran_MODULE_DIRECTORY
                                             ${CMAKE_BINARY_DIR}/modules)

list(APPEND libraries_to_link pic::pic pic-mpi::pic-mpi pic-blas::pic-blas)
if(NOT MQC_ENABLE_SERIAL)
  list(APPEND libraries_to_link OpenMP::OpenMP_Fortran)
endif()

if(MQC_ENABLE_TBLITE)
  list(APPEND libraries_to_link tblite::tblite)
endif()

# libxc, for exchange-correlation functionals.
#
# FetchContent rather than ExternalProject, and worth knowing why the difference
# shows up here: ExternalProject builds at *build* time, so its targets would
# not exist while this project is being configured, and the lines below would
# have to name the .a files by hand, guess the module directory and
# add_dependencies to force an order. Because it is fetched now, `xc` and
# `xcf03` are ordinary targets and their include directories, link order and
# build order all come along on their own.
#
# Resolved before libcint, and that order is load-bearing: check_dft is
# registered from backends/cenzontle/validation.cmake and links both.
if(MQC_ENABLE_LIBXC)
  find_package(libxc REQUIRED)
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_LIBXC)
  # xcf03 first, then xc: the Fortran interface depends on the C library, and a
  # static link resolves left to right.
  target_link_libraries(${main_lib} PRIVATE $<BUILD_INTERFACE:xcf03>
                                            $<BUILD_INTERFACE:xc>)
  # PUBLIC, and for the same reason the libcint include below is: nvfortran
  # needs every transitively used module file on the include path, not only the
  # ones a source names itself. `app/main.f90` uses no libxc module -- it uses
  # `mqc_driver`, which reaches `xc_f03_lib_m` several modules down -- and
  # nvfortran still has to open `xc_f03_lib_m.mod` to read them. PRIVATE keeps
  # it off the executable's include path and the build dies at 65 per cent with
  #
  # NVFORTRAN-F-0004-Unable to open MODULE file xc_f03_lib_m.mod
  #
  # gfortran does not need it, which is why this stood for so long.
  #
  # BUILD_INTERFACE because a fetched third-party library never joins this
  # project's export set.
  target_include_directories(${main_lib}
                             PUBLIC $<BUILD_INTERFACE:${libxc_BINARY_DIR}>)
  message(STATUS "libxc enabled: exchange-correlation functionals available")
endif()

# The CPU integrals backend: libfint by default, libcint behind
# -DMQC_USE_LIBFINT=OFF. The pins, and why each is pinned the way it is, are in
# cmake/modules/Findlibfint.cmake and cmake/modules/Findlibcint.cmake.
#
# The point of having it at all is not speed. Everything ab initio started out
# going through cuEST, which needs an A100 -- so the HF and DFT paths were
# compile-checked and never executed anywhere the developer was sitting, which
# is how they got that far unverified. A CPU backend makes them runnable on a
# laptop, and gives cuEST a second implementation to disagree with, which is the
# check it never had.
#
# libcint rather than libint, measured rather than assumed: clean clone to
# libcint.so in 9 seconds against 46 MB of pre-generated source and a long
# build. For a backend whose job is to be correct and always available, build
# time is the property that matters; if throughput ever does, libint can sit
# behind the same interface later.
#
# Nothing under backends/cenzontle/ changes with the choice: libfint defines the
# same `libcint_fortran` module with the same 25 procedures, so this program
# cannot tell which one it linked. That is the whole point of a port that is
# bit-identical -- the swap is a build option, not a migration.
#
# One thing the swap is *not* is purely a link: libfint carries L (sp) shells
# and libcint cannot represent one, so selecting it also defines
# MQC_WITH_SP_SHELLS and the fused-sp view in mqc_czt_integrals.F90 is compiled
# in. A libcint build takes the split-shell path instead. The two configurations
# therefore run different source, which is why both are built in CI rather than
# only the default.
if(MQC_ENABLE_CZT AND MQC_USE_LIBFINT)
  find_package(libfint REQUIRED)
  message(STATUS "libfint: ${MQC_LIBFINT_REPOSITORY}@${MQC_LIBFINT_TAG}")
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_LIBCINT)
  # libfint carries L (sp) shells; libcint cannot represent one at all.
  #
  # This is the first place the two backends are not interchangeable, so it is
  # the first thing here that has to be compiled differently rather than merely
  # linked differently. `build_sp_view` marks a fused sp shell with KAPPA_OF,
  # which libcint would read as a relativistic kappa and act on -- so the fused
  # view can only be built when libfint is what will read it.
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_SP_SHELLS)
  # Which of the two actually linked. MQC_WITH_LIBCINT above is set by both
  # branches -- it means "this build has CPU integrals" -- so it cannot answer
  # that question, and the startup acknowledgement has to name what ran rather
  # than what the interface is called.
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_LIBFINT)
  # One target, not two: libfint has no C library beside it.
  set(MQC_INTEGRALS_LIBS fint)
  # Answer to libcint's target names as well as our own.
  #
  # Most targets link ${MQC_INTEGRALS_LIBS} and so follow whichever backend is
  # selected. A few name `cint_fortran cint` directly, and every rebase onto
  # main brings more of them, because a contributor adding a test has no reason
  # to know this branch exists. Routing them one at a time is what three of the
  # four commits on this branch have been, and it will be due again after the
  # next merge.
  #
  # Aliases end that. With libfint selected libcint is never fetched, so the
  # names are free, and pointing both at `fint` makes a target that asks for
  # libcint by name get libfint without knowing. Nothing on main has to be aware
  # of this branch, and a rebase stops carrying link fixes.
  if(NOT TARGET cint)
    add_library(cint ALIAS fint)
  endif()
  if(NOT TARGET cint_fortran)
    add_library(cint_fortran ALIAS fint)
  endif()
  target_link_libraries(${main_lib} PRIVATE $<BUILD_INTERFACE:fint>)
  # PUBLIC, for the reason the libxc include above is and the libcint one below
  # is: nvfortran opens every transitively used module file, so anything linking
  # this library has to be able to find `cint_workspace.mod` even though no
  # source outside backends/cenzontle/ names it. The `fint` link is PRIVATE, so
  # the target's own interface include does not carry through -- which is why
  # this has to be stated rather than inherited.
  target_include_directories(
    ${main_lib} PUBLIC $<BUILD_INTERFACE:${libfint_BINARY_DIR}/include>)
  add_subdirectory(backends/cenzontle)
  message(
    STATUS
      "libfint enabled: CPU Gaussian integrals, all-Fortran, no C and no BLAS")

elseif(MQC_ENABLE_CZT)
  find_package(libcint REQUIRED)
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_LIBCINT)

  # BUILD_INTERFACE, not just PRIVATE: a static library records even its private
  # dependencies in its interface as $<LINK_ONLY:...>, so CMake wants them in
  # this project's export set -- which a fetched third-party library is never
  # going to join. Scoping them to the build keeps the integrals backend
  # internal, which it is.
  set(MQC_INTEGRALS_LIBS cint_fortran cint)
  target_link_libraries(${main_lib} PRIVATE $<BUILD_INTERFACE:cint_fortran>
                                            $<BUILD_INTERFACE:cint>)

  # The module directory has to propagate even though the library does not.
  #
  # nvfortran needs the `.mod` of every module that a used module itself uses,
  # transitively. So `use mqc_resources` in app/main.f90 makes it reach for
  # `libcint_fortran.mod`, which is two levels down and behind a PRIVATE link.
  # gfortran answers that from its own `.mod` and never asks, which is why this
  # was not needed until the first NVHPC build.
  #
  # BUILD_INTERFACE for the same reason the link above is: a fetched third-party
  # library is never going to join this project's export set.
  target_include_directories(
    ${main_lib} PUBLIC $<BUILD_INTERFACE:${libcint_BINARY_DIR}/include>)

  add_subdirectory(backends/cenzontle)
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

# CREST, for conformer and ensemble sampling.
#
# Resolved here rather than earlier because it must come after BLA_VENDOR is set
# and after tblite: CREST calls find_package(LAPACK/BLAS) unguarded and picking
# a different BLAS than the rest of this build links fine and computes wrong,
# and it guards tblite with `if(NOT TARGET "tblite::tblite")`, which this
# project has already defined above so that CREST links the same one. What it
# does with its other dependencies is in cmake/modules/Findcrest.cmake.
if(MQC_ENABLE_CREST)
  find_package(crest REQUIRED)
  target_link_libraries(${main_lib} PUBLIC lib-crest-static)
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_CREST)
  add_subdirectory(backends/crest)
  message(STATUS "CREST enabled: ${MQC_CREST_REPOSITORY}@${MQC_CREST_TAG}")

  # The seam the integration rests on, checked on its own because every way it
  # can fail is quiet: an uncalled pointer returns zero rather than an error.
endif()

# DL-FIND, behind libdlfind, for geometry optimization.
#
# Off by default, and the reason is licensing rather than size: DL-FIND is
# LGPL-3 and this program is MIT. libdlfind builds it as a shared library, which
# is what keeps the two separable -- the .so stays LGPL, mqc stays MIT, and a
# user retains the right to relink. Compiling those sources into libmqc.a
# instead would pull the relink obligation onto every binary shipped from here,
# so this deliberately does not add_subdirectory the DL-FIND sources.
#
# The interface is C, not Fortran: `api_dl_find` is bind(c) and takes the seven
# callbacks as function pointers, so nothing here needs libdlfind's .mod files
# and the two projects can be built by different compilers. That is also why
# there is no include-directory line below -- backends/dlfind declares the
# interface itself.
#
# A commit, not a branch, for the same reason libcint is pinned: what this
# program builds against should not change because somebody pushed to main.
if(MQC_ENABLE_DLFIND)
  find_package(libdlfind REQUIRED)
  target_compile_definitions(${main_lib} PRIVATE MQC_WITH_DLFIND)
  target_link_libraries(${main_lib} PRIVATE $<BUILD_INTERFACE:dlfind>)
  add_subdirectory(backends/dlfind)
  message(STATUS "DL-FIND enabled: geometry optimization via libdlfind")
endif()

# Cray: settle which BLAS the binary binds to, because the wrapper does not.
#
# `find_package(MPI)` learns what to link by asking the `ftn` wrapper for its
# link line, and that line names libsci -- `libsci_gnu`, `libsci_gnu_mpi` --
# whenever the cray-libsci module is loaded. `find_package(OpenMP)` asks with
# `-fopenmp` and gets the threaded pair, `libsci_gnu_mp` and
# `libsci_gnu_mpi_mp`. The wrapper also appends libsci to every link it performs
# itself. So libsci reaches the link line however the BLAS was chosen,
# `-DBLAS_LIBRARIES=` does not keep it out, and the dynamic linker then binds
# each `dgemm` and `dgetrf` to whichever BLAS-bearing library comes first in the
# executable's DT_NEEDED order. Reordering `libraries_to_link` does not settle
# that: CMake sorts the link line by dependency and MPI's copy lands first
# whatever order is written here. Rewriting MPI's imported target does not reach
# it either -- `MPI::MPI_Fortran` is created in pic-mpi's directory and is not
# visible from this one.
#
# Link options are emitted before any library, so naming a BLAS there is what
# puts it first. Two cases:
#
# A vendor was requested (`MQC_BLA_VENDOR`, defaulted to sequential MKL when
# MKLROOT is set): that BLAS goes first. libsci may still be loaded for MPI's
# sake, but nothing binds to it. This is what the `perlmutter-cpu` preset relies
# on, and it holds whether or not the module is loaded when a rebuild relinks --
# which is what broke the previous version of this block, a configure-time check
# of the OpenMP libraries that a later relink in another shell walked straight
# past. `mqc_check_blas_binding` reads the built executable to make sure.
#
# No vendor: libsci it is, and the threaded copy goes first. That is a
# correctness fix, not a performance one. Cray's serial libsci is not re-entrant
# and it fails quietly: two threads inside `dsygvd` at once corrupt each other's
# workspace, and on Perlmutter 190 of 200 concurrent calls on one matrix came
# back with different eigenvalues, the rest with `info` in the forties.
# Downstream that reads as an SCF that will not converge, or an eigensolver
# failing on a perfectly ordinary structure -- never as a threading problem, and
# never in a traceback naming the BLAS. The threaded copy is re-entrant only
# when called from one thread at a time, and it is no faster: see the
# response-solver block below.
set(mqc_libsci_mp)
foreach(mqc_omp_lib IN LISTS OpenMP_Fortran_LIBRARIES)
  if(mqc_omp_lib MATCHES "/libsci_[A-Za-z0-9_]*_mp\\.so")
    list(APPEND mqc_libsci_mp "${mqc_omp_lib}")
  endif()
endforeach()
set(mqc_blas_first)
if(MQC_BLA_VENDOR AND CMAKE_Fortran_COMPILER_WRAPPER STREQUAL "CrayPrgEnv")
  # The same answer pic-blas already got; FindBLAS caches it.
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  set(mqc_blas_first ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  list(REMOVE_DUPLICATES mqc_blas_first)
  target_link_options(${main_lib} PUBLIC ${mqc_blas_first})
  message(
    STATUS "Cray: ${BLA_VENDOR} linked ahead of libsci; the built executable "
           "is checked for it")
elseif(mqc_libsci_mp)
  target_link_options(${main_lib} PUBLIC ${mqc_libsci_mp})
  message(STATUS "Cray libsci: threaded copy linked first (${mqc_libsci_mp}); "
                 "the serial libsci the MPI wrapper names is not thread-safe")
  if(MQC_ENABLE_CZT)
    # The CPU backend calls the BLAS from inside its own parallel regions -- the
    # fitted Coulomb build, the ESP grid, every split GEMM -- and libsci,
    # threaded copy included, is not re-entrant under that: a 956-function
    # fitted SCF segfaulted in `B = (mn|Q) J^-1/2` where the MKL build ran.
    message(
      WARNING "Cray libsci is the BLAS of a CPU ab initio build. It is not "
              "re-entrant and its LAPACK slows down with every thread; expect "
              "segfaults in threaded fitted builds. Use the perlmutter-cpu "
              "preset (MKLROOT set, MQC_BLA_VENDOR=Intel10_64lp_seq) instead.")
  endif()
endif()

# Attach the post-link check to an executable. A no-op off a Cray or without a
# requested vendor, so every other platform's build is untouched.
function(mqc_check_blas_binding target)
  if(NOT (MQC_BLA_VENDOR AND CMAKE_Fortran_COMPILER_WRAPPER STREQUAL
                             "CrayPrgEnv"))
    return()
  endif()
  add_custom_command(
    TARGET ${target}
    POST_BUILD
    COMMAND
      ${CMAKE_COMMAND} -DBINARY=$<TARGET_FILE:${target}>
      -DVENDOR=${MQC_BLA_VENDOR} -P
      ${PROJECT_SOURCE_DIR}/cmake/MqcCheckBlasBinding.cmake
    VERBATIM)
endfunction()

# Whether the BLAS underneath threads itself, which decides where the response
# solver puts its parallelism.
#
# The twelve frequency solves are independent, and each is one `getrf` over the
# whole occupied-virtual space. Under a threaded BLAS that factorization already
# fills the machine, so the loop around it should stay serial. Under a
# sequential one it does not, and the solves have to be threaded here or the run
# spends its time on one core: at 8350 pairs, 175 seconds against 15.
#
# They do not compose. Nesting them is *slower* than a threaded BLAS alone --
# measured at 107 s against 102 s on a glycine tripeptide -- because MKL
# serialises a call made from inside a parallel region, leaving twelve-way
# concurrency where there was forty-way. Cray's libsci does not serialise: its
# `dgetrf` is not re-entrant, and thirteen concurrent calls spawned thirteen
# thousand threads and segfaulted inside `dtrsm`, or came back with "parameter
# number 1 had an illegal value". Serial is the only layout it survives, and it
# is slow there too -- one 6475 LU takes 3.8 s on one thread and gets slower
# with every thread added, 70 s on 64 -- so a Perlmutter build that cares about
# the response solve should unload cray-libsci and link MKL, where the same
# thirteen solves take 4 s at once. Left in for the builds that link libsci
# anyway, so that they finish rather than crash.
#
# Outside the backend branches above: it belongs to the response solver, which
# is the same code under libfint and libcint, and it once sat in the libcint
# branch alone -- so the default build, libfint, factorized every frequency on
# one core with the other thirty-nine idle. After the Cray block above, because
# that block is what settles which libsci the calls bind to, whatever
# `BLAS_LIBRARIES` was set to on the command line.
#
# Read off the vendor rather than probed at runtime, because there is no
# portable way to ask a BLAS how many threads it intends to use.
#
# An unknown vendor -- CMake's own detection, which is what CI gets -- is
# treated as sequential, because the two mistakes do not cost the same. Guess
# sequential wrongly and the nesting costs 5%; guess threaded wrongly and every
# solve runs on one core, which is the 175-second case. The exception is a Cray
# that binds to libsci, whose LAPACK cannot be called from two threads at once
# whatever the vendor string says.
if(MQC_ENABLE_CZT)
  if(mqc_libsci_mp AND NOT mqc_blas_first)
    message(
      STATUS
        "Response solver: frequencies serial (Cray libsci is not re-entrant)")
  elseif(NOT BLA_VENDOR OR BLA_VENDOR MATCHES "_seq$|_SEQ$")
    target_compile_definitions(${main_lib} PRIVATE MQC_SEQUENTIAL_BLAS)
    message(STATUS "Response solver: frequencies threaded (BLAS is sequential)")
  else()
    message(STATUS "Response solver: frequencies serial (BLAS threads itself)")
  endif()
endif()
unset(mqc_omp_lib)
unset(mqc_libsci_mp)
unset(mqc_blas_first)

target_link_libraries(${main_lib} PUBLIC ${libraries_to_link})

# Add json-fortran module directory for Fortran USE statements and ensure
# jsonfortran is built before metalquicha (for .mod files)
target_include_directories(
  ${main_lib} PUBLIC "$<BUILD_INTERFACE:${jsonfortran_BINARY_DIR}>")
add_dependencies(${main_lib} jsonfortran)
