# libdlfind -- DL-FIND, behind `driver: optimize`.
#
# OFF by default: DL-FIND carries its own licence, and whether to pull it in is
# a choice the person building makes rather than one made for them.
#
# Offline builds: point FETCHCONTENT_SOURCE_DIR_LIBDLFIND at a local clone, as
# with libcint. Compute nodes on most clusters have no network.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME
  libdlfind
  GIT_REPOSITORY
  "https://github.com/JorgeG94/libdlfind.git"
  GIT_TAG
  "4167998d16d8dac4a484ba9305f27d6325a7a28d"
  PROVIDES
  dlfind)
