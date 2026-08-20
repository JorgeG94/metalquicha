# cmake/CMakeLists.txt

# Release, not RelWithDebInfo.
#
# The two differ by more than a `-g`. RelWithDebInfo is `-O2` and picks up none
# of the architecture tuning below, which is attached to the Release
# configuration alone -- so the default build was leaving both `-O3` and
# `-march=native` on the table. It also propagates: an integral library fetched
# into this build inherits *this* configuration, and libfint's own tuning is
# expressed for Release, so a RelWithDebInfo parent quietly compiled it at `-O2`
# and undid it.
#
# `-march=native` comes with the Release configuration, and on a cluster that is
# the wrong default -- see the MQC_ARCH_FLAGS comment below, and pass
# `-DMQC_ARCH_FLAGS=` for anything built somewhere other than where it runs. The
# knob was already there; this only changes which configuration a build lands in
# when the caller says nothing.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Build type" FORCE)
  message(
    STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified.")
endif()

# Module installation
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# =============================================================================
# Manual Flag Override
# =============================================================================
# Use -DMQC_FC_EXTRA_FLAGS="..." to override all automatic flag detection
# =============================================================================
set(MQC_FC_EXTRA_FLAGS
    ""
    CACHE STRING "Override all Fortran flags (leave empty for auto-detection)")

if(NOT MQC_FC_EXTRA_FLAGS STREQUAL "")
  message(STATUS "Using manual Fortran flags: ${MQC_FC_EXTRA_FLAGS}")
  set(CMAKE_Fortran_FLAGS
      "${MQC_FC_EXTRA_FLAGS}"
      PARENT_SCOPE)
  return()
endif()

# =============================================================================
# Compiler Flags Configuration (auto-detected)
# =============================================================================
# Edit the flags below for your compiler. Each compiler has: - STANDARD  :
# Always applied - DEBUG     : Added for Debug builds - RELEASE   : Added for
# Release builds - WARNINGS  : Warning flags (optional) - COVERAGE  : For
# coverage builds (GCC only typically)
# =============================================================================

# Architecture tuning, overridable.
#
# `-march=native` tunes for the machine doing the *build*, which on a cluster is
# rarely the machine doing the running. Two ways that bites:
#
# * A build node newer than the compute nodes produces a binary that dies with
#   SIGILL on an instruction the compute node has never heard of.
# * A compiler newer than the assembler produces instructions the assembler has
#   never heard of. On Gadi, `-march=native` on a Sapphire Rapids node enables
#   AVX512-FP16, GCC then uses `vmovw` for ordinary 16-bit moves -- in tblite's
#   zip writer, whose headers are full of 2-byte fields -- and binutils older
#   than 2.38 rejects it. The fix there is a newer binutils if the module
#   exists, and this variable if it does not.
#
# Kept as the default because it is what this project has always done and it is
# right for a workstation. Override for anything shared:
#
# -DMQC_ARCH_FLAGS=                       portable, no tuning
# -DMQC_ARCH_FLAGS=-march=x86-64-v3       AVX2, runs on anything since 2015
# -DMQC_ARCH_FLAGS="-march=native -mno-avx512fp16"   native minus the bit the
# assembler lacks
set(MQC_ARCH_FLAGS
    "DEFAULT"
    CACHE STRING "Architecture flags for Release builds; empty for portable")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(FLAGS_STANDARD "-ffree-line-length-none -fbacktrace")
  set(FLAGS_DEBUG "-ffpe-trap=invalid,zero,overflow -fcheck=all")
  set(FLAGS_RELEASE "-march=native")
  set(FLAGS_WARNINGS "-Wall -Wextra -Wpedantic")
  set(FLAGS_COVERAGE "-O0 -g --coverage")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(FLAGS_STANDARD "-traceback")
  set(FLAGS_DEBUG "-check all -fpe0")
  set(FLAGS_RELEASE "-xHost")
  set(FLAGS_WARNINGS "-warn all")
  set(FLAGS_COVERAGE "")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  set(FLAGS_STANDARD "-traceback")
  set(FLAGS_DEBUG "-check all -fpe0")
  set(FLAGS_RELEASE "-xHost")
  set(FLAGS_WARNINGS "-warn all")
  set(FLAGS_COVERAGE "")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  set(FLAGS_STANDARD "-traceback -Minfo=mp")
  set(FLAGS_DEBUG "-Mbounds -Mchkptr")
  set(FLAGS_RELEASE "-fast")
  set(FLAGS_WARNINGS "")
  set(FLAGS_COVERAGE "")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  set(FLAGS_STANDARD "")
  set(FLAGS_DEBUG "")
  set(FLAGS_RELEASE "")
  set(FLAGS_WARNINGS "")
  set(FLAGS_COVERAGE "")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  # CrayPE's ftn wrapper driving the Cray Fortran compiler (CCE). CMake already
  # supplies the module flags (-em -J) and, for a PIC build, -h PIC.
  # Architecture tuning is the job of the loaded craype-x86-* target module
  # rather than a -march flag, so RELEASE carries none and Cray's own -O2
  # default stands.
  #
  # -homp is set here rather than left to find_package(OpenMP): CMake's
  # FindOpenMP detects the Cray compiler as OpenMP-capable with an *empty* flag
  # (older PrgEnv turned OpenMP on by default), so OpenMP_Fortran_FLAGS comes
  # back blank, the !$omp directives are ignored and the OpenMP runtime is never
  # linked. Keeping it in the standard set applies it to every compile, and
  # because CMAKE_Fortran_FLAGS also reaches the Cray link line the runtime is
  # pulled in.
  set(FLAGS_STANDARD "-homp")
  set(FLAGS_DEBUG "-g -Rbcps")
  set(FLAGS_RELEASE "")
  set(FLAGS_WARNINGS "")
  set(FLAGS_COVERAGE "")

else()
  message(WARNING "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
  set(FLAGS_STANDARD "")
  set(FLAGS_DEBUG "")
  set(FLAGS_RELEASE "")
  set(FLAGS_WARNINGS "")
  set(FLAGS_COVERAGE "")
endif()

# =============================================================================
# Apply flags to CMake variables
# =============================================================================
#
# **Everything that adds to CMAKE_Fortran_FLAGS accumulates into this local
# first, and the parent scope is written exactly once, at the end.** `set(...
# PARENT_SCOPE)` does not touch the local variable, so two of them in a row
# means the second reads the *unmodified* local and silently discards whatever
# the first added. That is not hypothetical: the coverage build lost
# `-ffree-line-length-none` this way and stopped being able to compile the
# banner, whose box-drawing characters push those lines past 132 bytes.
set(MQC_FORTRAN_FLAGS "${CMAKE_Fortran_FLAGS} ${FLAGS_STANDARD}")
set(CMAKE_Fortran_FLAGS_DEBUG
    "${CMAKE_Fortran_FLAGS_DEBUG} ${FLAGS_DEBUG}"
    PARENT_SCOPE)
# An explicit MQC_ARCH_FLAGS replaces the per-compiler default, including when
# it is set to the empty string -- which is how a portable build is asked for.
if(NOT MQC_ARCH_FLAGS STREQUAL "DEFAULT")
  set(FLAGS_RELEASE "${MQC_ARCH_FLAGS}")
  message(STATUS "Architecture flags overridden: '${MQC_ARCH_FLAGS}'")
endif()

set(CMAKE_Fortran_FLAGS_RELEASE
    "${CMAKE_Fortran_FLAGS_RELEASE} ${FLAGS_RELEASE}"
    PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO
    "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} ${FLAGS_RELEASE}"
    PARENT_SCOPE)

# Optional: Enable warnings via option
option(MQC_ENABLE_WARNINGS "Enable compiler warnings" OFF)
option(MQC_ENABLE_WERROR "Enable -Werror" OFF)
if(MQC_ENABLE_WARNINGS)
  set(MQC_FORTRAN_FLAGS "${MQC_FORTRAN_FLAGS} ${FLAGS_WARNINGS}")
  if(MQC_ENABLE_WERROR)
    set(MQC_FORTRAN_FLAGS "${MQC_FORTRAN_FLAGS} -Werror")
  endif()
endif()

# Coverage build type. `Coverage-mqc` is not one of CMake's own build types, so
# nothing else supplies flags for it -- these are all of them, and they have to
# sit alongside the standard set rather than replace it.
if(CMAKE_BUILD_TYPE STREQUAL "Coverage-mqc")
  set(MQC_FORTRAN_FLAGS "${MQC_FORTRAN_FLAGS} ${FLAGS_COVERAGE}")
endif()

set(CMAKE_Fortran_FLAGS
    "${MQC_FORTRAN_FLAGS}"
    PARENT_SCOPE)
