# Every build option this project defines, and the data two of them drive.
#
# One file so that `cmake -LH` has a written counterpart: what a flag is for,
# and what it costs, next to its default rather than three hundred lines from
# it. Included from the top-level before anything is resolved, because the
# dependency modules read these.
include_guard(GLOBAL)

option(MQC_ENABLE_TESTING "Enable testing for ${project_name}" ON)
option(MQC_ENABLE_TBLITE "Link to tblite for xTB calculations" ON)

#
# Build without MPI.
#
# mqc parallelises over fragments with MPI and that is how it is meant to run.
# The point of turning it off is not to run serially where MPI would do -- it is
# that a build should be possible where MPI is not: a CI runner without one, a
# laptop, or a compiler whose vendor MPI is awkward to obtain. Today each of
# those means no build at all, since `find_package(pic-mpi)` fails before
# anything else is tried.
#
# OFF selects pic-mpi's single-rank backend, which presents the same interface
# with one rank underneath. No mqc source changes: everything still goes through
# `pic_mpi_lib`.
#
# The flag propagated to pic-mpi is PIC_ENABLE_MPI, not this one -- pic-mpi has
# other consumers and none of them should have to know what MQC_ is. See
# cmake/modules/Findpic-mpi.cmake.
#
option(MQC_ENABLE_MPI "Build with MPI. OFF uses pic-mpi's single-rank backend"
       ON)
#
# Build with no parallelism and no semi-empirical engine: BLAS/LAPACK, and the
# libraries this project cannot be built without, and nothing else.
#
# This is one switch rather than three because the three are not independent.
# tblite is built with OpenMP and cannot be configured without it, so "no
# OpenMP" already implies "no tblite"; and a build with neither has no reason to
# carry MPI either. Asking for them one at a time is how you discover that in
# the wrong order.
#
# What it is for: a compiler that cannot yet build the parallel forms. LFortran
# is the case that prompted it -- it compiles this project's sources and links
# `libmetalquicha.a`, but its OpenMP support rejects a continued `!$omp`
# directive, which is how nearly every directive here is written. Serial builds
# what such a compiler can build instead of failing at the first `find_package`.
#
# It is not a performance mode. Nothing is threaded, nothing is distributed, and
# the fragment loop runs one fragment at a time. Do not benchmark against it.
#
# Normal (non-cache) sets so this shadows the two options without rewriting the
# cache, leaving a later reconfigure of the same tree without it intact -- the
# same reason the Coverage-mqc block below does it that way.
option(MQC_ENABLE_SERIAL
       "Build with no OpenMP, MPI or tblite. Only BLAS/LAPACK remain" OFF)
if(MQC_ENABLE_SERIAL)
  set(MQC_ENABLE_MPI OFF)
  set(MQC_ENABLE_TBLITE OFF)
  # A dependency asks for OpenMP itself and gets it whatever this project
  # decided -- libfint does, to make its shell-quartet workspace `!$omp
  # threadprivate`, and its `find_package` is its own. That would put the flag
  # back on the one library whose OpenMP a serial build most needs off: libfint
  # is where LFortran's named-`critical` gap bites. Disabling the find is how
  # CMake lets a parent answer for a subproject it does not own.
  #
  # It costs nothing here. libfint warns that shared workspaces are "only safe
  # to call from one thread at a time", which is the serial build exactly.
  set(CMAKE_DISABLE_FIND_PACKAGE_OpenMP ON)
  message(
    STATUS "Serial build: OpenMP, MPI and tblite disabled. Threading and rank "
           "parallelism are both absent -- this is a portability mode, not a "
           "performance one")
endif()

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
option(
  MQC_ENABLE_CZT
  "Build cenzontle, the CPU ab initio backend (real ab initio without a GPU)"
  ON)

# This option was `MQC_ENABLE_LIBCINT` until the backend was named. Honour the
# old spelling rather than ignoring it: CMake accepts an unknown -D silently, so
# a script still passing it would otherwise configure a build with no CPU ab
# initio path at all and say nothing about why every deck was refused.
if(DEFINED MQC_ENABLE_LIBCINT)
  message(
    WARNING
      "MQC_ENABLE_LIBCINT is now MQC_ENABLE_CZT (the backend is cenzontle, and "
      "it is not a libcint wrapper -- libfint backs it by default). Using the "
      "old value, ${MQC_ENABLE_LIBCINT}, for MQC_ENABLE_CZT.")
  set(MQC_ENABLE_CZT
      "${MQC_ENABLE_LIBCINT}"
      CACHE BOOL "" FORCE)
endif()
option(MQC_ENABLE_LIBXC
       "Fetch libxc for exchange-correlation functionals (CPU DFT)" ON)

# Which of the two CPU integral libraries provides them. They define the same
# `libcint_fortran` module with the same procedures, so no source under
# backends/cenzontle/ knows which one it linked -- but they are not the same
# build: libfint carries L (sp) shells and turns on a different code path for
# them. See cmake/MqcDependencies.cmake, which is where that happens.
option(MQC_USE_LIBFINT
       "Use libfint (the Fortran port) instead of libcint for CPU integrals" ON)

# A coverage build measures this project's own Fortran, not its dependencies --
# and `backends/cenzontle/` is this project's own Fortran. It used to be dropped
# here alongside libxc, which meant every line of the CPU backend -- the SCF,
# MP2, coupled cluster, MCSCF, SAPT, EFP, the analyses -- was invisible to
# coverage rather than uncovered by it. Half the program was not being measured,
# and a number that does not include it is not this program's coverage.
#
# libxc is a different case and stays off: its autocmake safeguard allows only
# Debug/Release/RelWithDebInfo and aborts on "Coverage-mqc", so a coverage
# configure that fetched it dies before a line of mqc is compiled. That costs
# the Kohn-Sham paths, which the build already refuses without libxc anyway.
#
# A normal (non-cache) set so it shadows the option here and in the
# subdirectories without rewriting the cache, leaving a later reconfigure of the
# same tree under another build type intact.
if(CMAKE_BUILD_TYPE STREQUAL "Coverage-mqc")
  set(MQC_ENABLE_LIBXC OFF)
  message(STATUS "Coverage-mqc build: libxc disabled (it refuses this build "
                 "type); libcint measured as part of the project")
endif()

# Off by default: DL-FIND is LGPL-3 and this program is MIT, so linking it is a
# choice the person building makes rather than one made for them. See the
# FetchContent block below for how the two licenses are kept separable.
option(MQC_ENABLE_DLFIND "Fetch libdlfind (DL-FIND) for geometry optimization"
       OFF)

# Which CODATA revision the Bohr radius in mqc_physical_constants follows. Only
# the length conversion is selected this way; see that module for why.
set(MQC_CODATA_YEAR
    "2018"
    CACHE STRING "CODATA revision for the Bohr radius: 2018 | 2010")
set_property(CACHE MQC_CODATA_YEAR PROPERTY STRINGS 2018 2010)

option(MQC_ENABLE_CUEST "Link to NVIDIA cuEST for GPU integrals (HF/DFT)" OFF)
option(MQC_ENABLE_CREST "Link to CREST for conformer and ensemble sampling" OFF)
option(MQC_ENABLE_TERCO "Link to terco for a device-resident SCF (HF/DFT)" OFF)
set(MQC_LIBFINT_REPOSITORY
    "https://github.com/JorgeG94/libfint.git"
    CACHE STRING "Where to fetch libfint from")
set(MQC_LIBFINT_TAG
    "v0.3.1"
    CACHE STRING "libfint revision to build against")

set(MQC_CREST_REPOSITORY
    "https://github.com/JorgeG94/crest.git"
    CACHE STRING "Where to fetch CREST from")
set(MQC_CREST_TAG
    "feat/engrad-callback"
    CACHE STRING "CREST revision to build against")
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
    "sto-3g;3-21g;6-31g;6-31g_st_;6-31g_st__st_;cc-pvdz;cc-pvtz;cc-pvqz;6-31g(2df,p);aug-cc-pvdz;aug-cc-pvtz;cc-pvtz-jkfit;cc-pvdz-rifit;cc-pvtz-rifit"
    # The 6-311 family. Two things want it. MAKEFP, because an effective
    # fragment potential wants a large Pople set with diffuse functions and the
    # f shells in (3df,3pd) are what the polarizability and dispersion tensors
    # are built from; note the bundle has no 6-311++G(3df,2p), the set the EFMO
    # literature names, so (3df,3pd) is the nearest thing it carries. And the
    # quasi-atomic bonding papers, whose numbers are quoted in that same set and
    # so can be reproduced rather than approached in a smaller one.
    "6-311g;6-311g_st_;6-311g_st__st_;6-311++g;6-311++g_st__st_;6-311++g(2d,2p);6-311++g(3df,3pd)"
    # every def2 the bundle carries: orbital, diffuse, RI/JK auxiliaries, ECP
    "def2-sv(p);def2-sv(p)-jkfit;def2-sv(p)-rifit;def2-svp;def2-svp-rifit;def2-svpd;def2-svpd-rifit"
    "def2-tzvp;def2-tzvp-rifit;def2-tzvpd;def2-tzvpd-rifit;def2-tzvpp;def2-tzvpp-rifit;def2-tzvppd;def2-tzvppd-rifit"
    "def2-qzvp;def2-qzvp-rifit;def2-qzvpd;def2-qzvpp;def2-qzvpp-rifit;def2-qzvppd;def2-qzvppd-rifit"
    "def2-mtzvp;def2-mtzvpp;def2-mtzvpp-rij;def2-universal-jfit;def2-universal-jkfit;def2-ecp"
    CACHE STRING "Basis sets to extract from ${MQC_BASIS_BUNDLE}")

include(ExtractBasisSets)
mqc_extract_basis_sets()
