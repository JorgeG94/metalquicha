set(_lib "pic-mpi")
set(_url "https://github.com/JorgeG94/pic-mpi/")

# Turn MQC_ENABLE_MPI=OFF into the flag pic-mpi actually understands.
#
# pic-mpi is used by more than mqc, so it cannot be asked to recognise an
# MQC_-prefixed option. It has its own, PIC_ENABLE_MPI, and the translation
# belongs here -- pre-set as a cache entry so that when pic-mpi is fetched as a
# subproject its `option()` sees this value instead of applying its default.
if(DEFINED MQC_ENABLE_MPI AND NOT MQC_ENABLE_MPI)
  set(PIC_ENABLE_MPI
      OFF
      CACHE BOOL "Link against MPI (set by MQC_ENABLE_MPI=OFF)" FORCE)
endif()

# Pass PIC_USE_LEGACY_MPI option to the fetched package if set
if(DEFINED PIC_USE_LEGACY_MPI)
  set(PIC_USE_LEGACY_MPI
      ${PIC_USE_LEGACY_MPI}
      CACHE BOOL "Use legacy MPI module" FORCE)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

# Use the main branch
set(_rev "v0.5.2")

# The first pic-mpi release carrying the single-rank backend, i.e. the first
# one that understands PIC_ENABLE_MPI.
set(_serial_min "0.6.0")

# Checked BEFORE the fetch, on purpose. FetchContent configures pic-mpi as soon
# as it is made available, so a pic-mpi that predates the option fails inside
# its own find_package(MPI) -- with an error naming FindMPI and nothing about
# why MPI was wanted. Anything checked afterwards runs too late to help.
# `FETCHCONTENT_SOURCE_DIR_PIC-MPI` means a developer has pointed the build at
# a local checkout, where the pinned tag says nothing about what is in it. That
# is how the serial backend gets tested before it is released, so the version
# check does not apply.
if(DEFINED MQC_ENABLE_MPI
   AND NOT MQC_ENABLE_MPI
   AND NOT FETCHCONTENT_SOURCE_DIR_PIC-MPI)
  string(REGEX REPLACE "^v" "" _rev_number "${_rev}")
  if(_rev_number VERSION_LESS _serial_min)
    message(
      FATAL_ERROR
        "MQC_ENABLE_MPI=OFF needs a pic-mpi with the single-rank backend "
        "(>= v${_serial_min}), and this pins ${_rev}. Bump _rev in "
        "cmake/modules/Findpic-mpi.cmake, or configure with "
        "MQC_ENABLE_MPI=ON.")
  endif()
  unset(_rev_number)
endif()
my_fetch_package("${_lib}" "${_url}" "${_rev}")

unset(_lib)
unset(_url)
unset(_rev)
unset(_serial_min)
