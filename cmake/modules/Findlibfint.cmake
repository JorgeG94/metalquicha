# libfint -- the all-Fortran port of libcint, and the default CPU integrals.
#
# Pinned by tag now that it is the default. Tracking `main` was right while this
# was the alternative backend and one job built it: an upstream commit showed up
# as a red square here, which is the drift detection the two repositories wanted
# from each other. As the default it is every job, so the same commit would
# redden the whole matrix on a metalquicha change that touched nothing, and the
# signal would be indistinguishable from a real break.
#
# v0.1.3 is 9a9926aacf8a4fe2fb364454dd45fbb39cf82dbf today. It is the first
# release carrying the scalar ECP integrals, which nothing else provides:
# libcint has no ECP code at all, so `mqc_libcint_ecp` refuses `model.ecp`
# outright in a `-DMQC_USE_LIBFINT=OFF` build rather than failing at the linker.
# It also succeeds v0.1.1, the first release whose C ABI carries the second
# derivatives -- which is what lets `mqc_libcint_hess_abi` declare one set of
# entry points for both backends -- and v0.1.2, the first whose F12 integrals
# build. Both are ancestors of it, so neither is given up here.
#
# libcint pins by SHA rather than by name, which is the stronger form: a tag can
# be moved and a SHA cannot. This one is a tag because it was a *branch* for the
# length of the ECP work, and a moving pin reddens the whole matrix with a
# signal indistinguishable from a real break. A release is the fixed point that
# ends that, and the SHA is written above so a moved tag is detectable.
#
# Overridable through MQC_LIBFINT_REPOSITORY and MQC_LIBFINT_TAG, which is the
# point of doing it that way rather than editing the literal:
# `-DMQC_LIBFINT_TAG=v0.1.2` builds against the previous release and turns the
# ECP support off, so a bisect over an unrelated failure has a way back.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME
  libfint
  GIT_REPOSITORY
  "${MQC_LIBFINT_REPOSITORY}"
  GIT_TAG
  "${MQC_LIBFINT_TAG}"
  PROVIDES
  fint)
