# dftd4 -- the charge-dependent D4 correction, behind `keywords.dft.dispersion`.
#
# OFF by default: dftd4 is LGPL-3.0-or-later and this program is MIT, so whether
# to pull it in is a choice the person building makes. It is used through its C
# API (include/dftd4.h) and built shared, which is what keeps the two licences
# separable -- the same treatment s-dftd3 and libdlfind get.
#
# The pin has to agree with the one tblite carries in
# config/cmake/Finddftd4.cmake, because with both enabled there is exactly one
# copy and it is this one. See the block in cmake/MqcDependencies.cmake.
#
# dftd4 resolves mctc-lib, multicharge and a BLAS of its own at its top level;
# all three are already in this build tree for other reasons, so nothing extra
# is declared here.
#
# Offline builds: point FETCHCONTENT_SOURCE_DIR_DFTD4 at a local clone.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME
  dftd4
  GIT_REPOSITORY
  "${MQC_DFTD4_REPOSITORY}"
  GIT_TAG
  "${MQC_DFTD4_TAG}"
  NAMESPACED_TARGET
  PROVIDES
  dftd4-lib)
