# pic-mpi -- the MPI wrapper layer, and the one dependency with a build option
# of its own to translate.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

# TEMPORARY: the branch, not a release. v0.6.5 predates the 1-D `irecv` fix
# (pic-mpi#26) and the removal of the distributed-arrays module (pic-mpi#27),
# and this build exists to put both in front of this project's CI before they
# are tagged. Put a `vX.Y.Z` back the moment there is one -- `MqcFetch.cmake`
# documents GIT_TAG as a tag or a SHA precisely so a build is reproducible, and
# a branch is not.
set(_rev "main")

# The first pic-mpi release carrying the single-rank backend, i.e. the first one
# that understands PIC_ENABLE_MPI.
set(_serial_min "0.6.0")

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

# pic-mpi's own switches, forwarded when the configure named one. Nothing here
# chooses them; they are how a machine that needs a workaround asks for it, and
# `CMakePresets.json` is where the combinations that machines need are written
# down -- `perlmutter` sets the first two.
#
# PIC_USE_LEGACY_MPI     the old `mpi` module instead of `mpi_f08`
# PIC_MPICH_PERLMUTTER   drop the MPI_Win_allocate wrappers, whose generic does
# not resolve under nvfortran with Cray MPICH. Safe here because nothing in this
# program uses RMA PIC_USE_VAPAA          reach MPI through its C ABI. The
# bigger hammer: if a second generic-resolution mismatch like the above turns
# up, this fixes the class rather than the instance
foreach(_pic_mpi_opt PIC_USE_LEGACY_MPI PIC_MPICH_PERLMUTTER PIC_USE_VAPAA)
  if(DEFINED ${_pic_mpi_opt})
    set(${_pic_mpi_opt}
        ${${_pic_mpi_opt}}
        CACHE BOOL "Forwarded to pic-mpi" FORCE)
  endif()
endforeach()
unset(_pic_mpi_opt)

# Checked BEFORE the fetch, on purpose. FetchContent configures pic-mpi as soon
# as it is made available, so a pic-mpi that predates the option fails inside
# its own find_package(MPI) -- with an error naming FindMPI and nothing about
# why MPI was wanted. Anything checked afterwards runs too late to help.
# `FETCHCONTENT_SOURCE_DIR_PIC-MPI` means a developer has pointed the build at a
# local checkout, where the pinned tag says nothing about what is in it. That is
# how the serial backend gets tested before it is released, so the version check
# does not apply.
if(DEFINED MQC_ENABLE_MPI
   AND NOT MQC_ENABLE_MPI
   AND NOT FETCHCONTENT_SOURCE_DIR_PIC-MPI)
  # Only a `vX.Y.Z` pin can be compared. A branch name is not a version, and
  # CMake reads one as 0 -- so `main`, which is ahead of every release, would
  # fail this check claiming to be older than v0.6.0. A branch is pinned
  # deliberately and by someone who knows what is on it, which is the same
  # reason FETCHCONTENT_SOURCE_DIR_PIC-MPI is exempt above.
  set(_rev_number "")
  if(_rev MATCHES "^v?[0-9]+\\.[0-9]+\\.[0-9]+$")
    string(REGEX REPLACE "^v" "" _rev_number "${_rev}")
  endif()
  if(_rev_number AND _rev_number VERSION_LESS _serial_min)
    message(
      FATAL_ERROR
        "MQC_ENABLE_MPI=OFF needs a pic-mpi with the single-rank backend "
        "(>= v${_serial_min}), and this pins ${_rev}. Bump _rev in "
        "cmake/modules/Findpic-mpi.cmake, or configure with "
        "MQC_ENABLE_MPI=ON.")
  endif()
  unset(_rev_number)
endif()

mqc_fetch(
  NAME
  pic-mpi
  GIT_REPOSITORY
  "https://github.com/JorgeG94/pic-mpi/"
  GIT_TAG
  "${_rev}"
  NAMESPACED_TARGET)

unset(_rev)
unset(_serial_min)
