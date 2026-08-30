# libcint -- the C Gaussian integral library, selected by MQC_USE_LIBFINT=OFF.
#
# JorgeG94's fork rather than upstream sunqm: it adds `libcint_fortran`, a
# Fortran API over the C one, so this program needs no hand-written bindings at
# all -- including libcint_gto_norm, which is the easy thing to get wrong
# (libcint wants contraction coefficients pre-multiplied by the primitive
# normalisation, and an overlap diagonal that is not 1 is the symptom).
#
# A commit, not a branch. libcint develops on `develop`, so tracking that would
# change what this program builds against on every push there, with nothing here
# recording that it changed. This hash is reachable from develop and carries
# libcint_3c2e_sph and libcint_2c2e_sph, without which the density-fitting path
# will not compile -- and their Cartesian twins libcint_3c2e_cart,
# libcint_2c2e_cart and libcint_tot_cgto_cart, without which a Cartesian basis
# such as 6-31G* can reach the exact four-index path but not the density-fitted
# one. It also carries the three- and two-centre optimizers: libcint exports
# them, the Fortran interface had no way to make one, so every fitting integral
# rebuilt its shell-pair data per call. This hash also carries the
# nuclear-derivative integrals -- ipkin, ipnuc and iprinv, plus the three-centre
# and two-centre derivative pairs -- without which there is no analytic gradient
# on this backend at all. Worth replacing with a tag once the fork has one -- a
# tag says what it is, a hash does not, and a tag also lets GIT_SHALLOW come
# back, which a bare hash cannot resolve.
#
# Offline builds: point FETCHCONTENT_SOURCE_DIR_LIBCINT at a local clone and no
# network is touched. Compute nodes on most clusters have none, and a configure
# step that reaches for GitHub fails there in a way that reads like a compiler
# problem.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME
  libcint
  GIT_REPOSITORY
  "https://github.com/JorgeG94/libcint.git"
  GIT_TAG
  "3c78069d6c92e9d7f4db792ea2a33779b21539cc"
  PROVIDES
  cint
  cint_fortran)
