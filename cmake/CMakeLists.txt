# cmake/CMakeLists.txt

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "RelWithDebInfo"
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
set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} ${FLAGS_STANDARD}"
    PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS_DEBUG
    "${CMAKE_Fortran_FLAGS_DEBUG} ${FLAGS_DEBUG}"
    PARENT_SCOPE)
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
  set(CMAKE_Fortran_FLAGS
      "${CMAKE_Fortran_FLAGS} ${FLAGS_WARNINGS}"
      PARENT_SCOPE)
  if(MQC_ENABLE_WERROR)
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} -Werror"
        PARENT_SCOPE)
  endif()
endif()

# Coverage build type
if(CMAKE_BUILD_TYPE STREQUAL "Coverage-mqc")
  set(CMAKE_Fortran_FLAGS
      "${CMAKE_Fortran_FLAGS} ${FLAGS_COVERAGE}"
      PARENT_SCOPE)
endif()
