# s-dftd3 -- simple-dftd3, behind `keywords.dft.dispersion`.
#
# OFF by default: s-dftd3 is LGPL-3-or-later and this program is MIT, so whether
# to pull it in is a choice the person building makes. It is used through its C
# API (include/s-dftd3.h) and built shared, which is what keeps the two licences
# separable -- the same treatment libdlfind gets.
#
# The pin has to agree with the one tblite carries in
# config/cmake/Finds-dftd3.cmake, because with both enabled there is exactly one
# copy and it is this one. See the block in cmake/MqcDependencies.cmake.
#
# Offline builds: point FETCHCONTENT_SOURCE_DIR_S-DFTD3 at a local clone.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME
  s-dftd3
  GIT_REPOSITORY
  "${MQC_DFTD3_REPOSITORY}"
  GIT_TAG
  "${MQC_DFTD3_TAG}"
  NAMESPACED_TARGET
  PROVIDES
  s-dftd3-lib)
