# libxc -- the exchange-correlation functionals.
#
# A release tarball rather than a clone: libxc's history is large and nothing
# here needs it, and the hash pins the contents in a way a git tag does not.
#
# Static, so nothing has to be found on the library path at run time; the global
# POSITION_INDEPENDENT_CODE keeps it usable from libmqc.so. These have to be set
# before the fetch, because they are what the subproject's own `option()` calls
# see. BUILD_TESTING is handled once for every dependency in
# cmake/MqcDependencies.cmake.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

set(BUILD_SHARED_LIBS
    OFF
    CACHE BOOL "" FORCE)
set(ENABLE_FORTRAN
    ON
    CACHE BOOL "" FORCE)
set(ENABLE_PYTHON
    OFF
    CACHE BOOL "" FORCE)

# Third functional derivatives, which libxc does not build by default.
#
# `MAXORDER` is 2 upstream, and asking a functional for `kxc` in such a build
# does not return an error: libxc prints "does not provide an implementation of
# kxc" and *aborts the process*. So this is not an optimisation to be skipped
# when unneeded -- without it a double-hybrid Hessian kills the run.
#
# That Hessian is what needs them. Its perturbative term is not stationary in
# the orbitals, so the perturbed Z-vector's operator is the differentiated
# kernel, and the fixed-density skeleton term carries `g_xc rho^X rho^Y`;
# omitting that one is a 400 per cent error on it rather than a correction to
# it. Nothing below a double hybrid needs this -- a Kohn-Sham gradient stops at
# `v_xc` and a Kohn-Sham Hessian at `f_xc` -- and it costs about 40 per cent of
# libxc's build time, measured: 970 s against 696 s.
#
# Here rather than beside the fetch for the reason the header above gives: this
# is what the subproject's own `option()` sees, and it only sees it if it is set
# first.
set(MAXORDER
    3
    CACHE STRING "" FORCE)

mqc_fetch(
  NAME
  libxc
  URL
  "https://gitlab.com/libxc/libxc/-/archive/7.1.2/libxc-7.1.2.tar.gz"
  URL_HASH
  "SHA256=c517ce61820ea8114664a4280b6a6bc74a4f22f1fd1ea4ddecd6df0caeeae4f4"
  PROVIDES
  xc
  xcf03)
