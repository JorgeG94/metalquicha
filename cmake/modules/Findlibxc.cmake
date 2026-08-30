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
