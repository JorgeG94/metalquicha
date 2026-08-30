# CREST -- conformer and ensemble sampling.
#
# Resolved after BLA_VENDOR is set, because CREST calls
# find_package(LAPACK/BLAS) unguarded -- the only two of its dependencies not
# behind an `if(NOT TARGET ...)` -- and picking a different BLAS than the rest
# of this build links fine and computes wrong. MqcDependencies.cmake is what
# keeps that order.
#
# Nothing here declares toml-f, gfnff, gfn0, libpvol or lwoniom, and that is
# deliberate: CREST's own config/modules/crest-utils.cmake resolves each of them
# through "subproject, cmake, fetch, pkgconf" in order and, when it has to
# fetch, creates the `pkg::pkg` alias itself. Those projects export their
# namespaced target only from install(), so without that alias a FetchContent
# build of them would not satisfy CREST's own `if(NOT TARGET ...)` guards.
# Duplicating any of it here would fight the mechanism rather than help it.
#
# tblite is the exception worth naming: CREST guards it with `if(NOT TARGET
# "tblite::tblite")`, and this project defines that target already, so CREST
# links the same tblite rather than building a second one.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

# Features this integration does not use, off by default so they are neither
# built nor fetched at an unpinned HEAD. A deck that wants one can turn it back
# on, and will then be choosing to resolve its version too.
foreach(_crest_off WITH_GFNFF WITH_LIBPVOL WITH_LWONIOM)
  if(NOT DEFINED ${_crest_off})
    set(${_crest_off}
        OFF
        CACHE BOOL "CREST optional backend, off unless asked for" FORCE)
  endif()
endforeach()
unset(_crest_off)

# GIT_SUBMODULES names gfn0 and nothing else, and both halves matter.
#
# CREST carries gfn0, gfnff, toml-f, lwoniom, pvol, tblite and test-drive as git
# submodules under subprojects/, and FetchContent updates every one of them
# unless told otherwise. CREST's resolver then prefers a populated subprojects/
# directory over everything else, so it add_subdirectory()s its own toml-f --
# which tblite has already fetched into _deps by the time we get here. Two
# directories defining `toml-f-lib` is a hard CMake error that names CMP0002 and
# nothing about CREST, so it does not read as a dependency conflict at all.
# Excluding toml-f and tblite lets CREST's own `if(NOT TARGET ...)` guards see
# the targets this project already made and skip both, which is the sharing the
# arrangement wanted.
#
# But the submodules are *version pins*, so excluding one is not free: with none
# checked out CREST fetched HEAD of gfnff and stopped on a `print=` keyword its
# own source still passes. gfn0 is the one dependency it needs that nothing here
# provides -- the MTD-length heuristic calls it for Wiberg bond orders -- so it
# keeps its pin. The rest are switched off above rather than fetched at whatever
# HEAD happens to be.
mqc_fetch(
  NAME
  crest
  GIT_REPOSITORY
  "${MQC_CREST_REPOSITORY}"
  GIT_TAG
  "${MQC_CREST_TAG}"
  GIT_SUBMODULES
  "subprojects/gfn0"
  PROVIDES
  lib-crest-static)
