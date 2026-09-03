# libfint -- the all-Fortran port of libcint, and the default CPU integrals.
#
# Pinned by tag now that it is the default. Tracking `main` was right while this
# was the alternative backend and one job built it: an upstream commit showed up
# as a red square here, which is the drift detection the two repositories wanted
# from each other. As the default it is every job, so the same commit would
# redden the whole matrix on a metalquicha change that touched nothing, and the
# signal would be indistinguishable from a real break.
#
# v0.2.1 is d02d7dd00bd14d2230e0ece7a42faa9738aa8128 today. It is the first
# release whose `libcint_fortran` carries the scalar ECP -- `libcint_ecp_sph`
# and `libcint_ecp_cart` -- so `mqc_libcint_ecp` reaches it the way every other
# integral here is reached, instead of declaring a `bind(C)` interface to
# libcint's own names and marshalling `atm`, `bas` and `env` into C-typed copies
# to call it. libcint has no ECP code at all, so this is the one integral the
# two backends cannot both answer and the `#ifdef` around it stays: what it
# marks now is a gap in what the libcint cross-check can cover, not a
# portability constraint.
#
# It succeeds v0.1.4, the first release whose C ABI exposes `int1e_irp`, `<i| r
# nabla |j>` -- the nuclear derivative of the dipole integrals, and so what an
# analytic infrared intensity needs. The integral was in the library before
# this; what it lacked was a C entry point, and reaching the Fortran module
# directly for it would have made the caller unlinkable against libcint, which
# is the one property this indirection exists to keep.
#
# It succeeds v0.1.3, the first release carrying the scalar ECP integrals, which
# nothing else provides: libcint has no ECP code at all, so `mqc_czt_ecp`
# refuses `model.ecp` outright in a `-DMQC_USE_LIBFINT=OFF` build rather than
# failing at the linker. It also succeeds v0.1.1, the first release whose C ABI
# carries the second derivatives -- which is what lets `mqc_czt_hess_abi`
# declare one set of entry points for both backends -- and v0.1.2, the first
# whose F12 integrals build. Both are ancestors of it, so neither is given up
# here.
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
