# OpenTrustRegion -- the second-order trust-region orbital optimizer.
#
# Pinned by SHA rather than by branch, which is the stronger form: a tag can be
# moved and a SHA cannot. 8fa7769 is the commit whose `stability_check` this
# project was written against, and the pin is overridable through
# MQC_OTR_REPOSITORY and MQC_OTR_TAG so a bisect over an upstream change has a
# way back. Note that MQC_OTR_TAG is a cache variable: bumping the literal here
# does not reach a build tree that has already been configured, which has to be
# done with `-DMQC_OTR_TAG=` or a fresh tree.
#
# ## The three options forced below, and why each one has to be
#
# `INTEGER_SIZE` selects the width of OpenTrustRegion's own `ip` kind, and
# everything crossing the interface -- the error codes, the settings counts --
# is declared with it. It has to equal pic's `default_int`, which is `int32`
# unless pic itself is built with USE_INT8. Nothing in this project builds pic
# that way, so 4 is the matching width; there is no way to ask the library for
# both, and a mismatch is a silent one -- an `integer(int64)` error code read
# through an `integer(int32)` dummy is the low half of it, which is zero for
# every error code the library defines.
#
# `OpenTrustRegion_HOST_PROVIDES_BLAS` stops the library finding and linking a
# BLAS of its own, leaving it to resolve `ddot`, `dnrm2`, `dgemm`, `dgemv` and
# `dsyev` against whatever this project already linked. That is what keeps one
# BLAS in the binary: the library's own `find_package(BLAS)` would pick up
# threaded MKL here, and a second threaded BLAS in the same process is the
# failure mode this project spends cmake/MqcCheckBlasBinding.cmake detecting. It
# also needs the same integer width, which is why the option is documented
# upstream as requiring INTEGER_SIZE to be set -- the two are one decision.
#
# `OpenTrustRegion_ENABLE_XHOST` defaults ON upstream and appends -march=native
# (or -xHost) to its Release flags. This project settles architecture flags
# once, in MQC_ARCH_FLAGS, precisely so that a build on a login node runs on the
# compute nodes; a dependency quietly compiling for the build host breaks every
# cross-compiled cluster build with an illegal instruction at run time rather
# than at link time. Forced OFF.
#
# `OpenTrustRegion_BUILD_TESTING` would build a shared test library that links a
# BLAS at build time, which the option above has just declined to find. Upstream
# already defaults it to PROJECT_IS_TOP_LEVEL, false here; forced for the same
# reason every other dependency's suite is forced off in
# cmake/MqcDependencies.cmake.
#
# These are set before the fetch because they are what the subproject's own
# `option()` and `set(... CACHE)` calls see, the same ordering constraint
# Findlibxc.cmake documents.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

set(INTEGER_SIZE
    "4"
    CACHE STRING "" FORCE)
set(OpenTrustRegion_HOST_PROVIDES_BLAS
    ON
    CACHE BOOL "" FORCE)
set(OpenTrustRegion_ENABLE_XHOST
    OFF
    CACHE BOOL "" FORCE)
set(OpenTrustRegion_BUILD_TESTING
    OFF
    CACHE BOOL "" FORCE)

mqc_fetch(
  NAME
  opentrustregion
  GIT_REPOSITORY
  "${MQC_OTR_REPOSITORY}"
  GIT_TAG
  "${MQC_OTR_TAG}"
  PROVIDES
  opentrustregion)
