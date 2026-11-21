# libfint -- the all-Fortran port of libcint, and the CPU integrals.
#
# Pinned by tag, because the integrals are every job: tracking `main` would
# redden the whole matrix on an upstream commit that touched nothing here, with
# a signal indistinguishable from a real break. A tag can be moved and a SHA
# cannot, so the SHA is written below to make a moved tag detectable.
#
# v0.9.0 is 57fccec6014117a968ec7d0fae72eac135f135f7 today. It is the first
# release whose `libcint_fortran` carries everything `backends/cenzontle/` calls
# directly: the gauge-including and nuclear-moment integrals a shielding needs
# (`libcint_1e_igovlp` through `libcint_2e_ig1`), the multipoles, the
# point-charge field, the grid integral, the one- and two-electron second
# derivatives, `int1e_irp`, the scalar ECP (`libcint_ecp_sph` and
# `libcint_ecp_cart`) and its first derivatives (`libcint_ecp_ipnuc_*` and
# `libcint_ecp_iprinv_*`), and the fitting integrals' second derivatives
# (`libcint_3c2e_ipip1` and the rest), which the RI-MP2 Hessian used to reach
# through libfint's internal modules. Its module is called `libcint_fortran` and
# its procedures and constants keep libcint's `libcint_*` and `LIBCINT_*`
# spelling.
#
# Overridable through MQC_LIBFINT_REPOSITORY and MQC_LIBFINT_TAG, which is the
# point of doing it that way rather than editing the literal:
# `-DMQC_LIBFINT_TAG=<tag>` builds against another release, so a bisect over an
# unrelated failure has a way back. Anything older than the pin lacks procedures
# this backend calls and will not link: the ECP derivatives are in no release up
# to v0.7.0.
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
