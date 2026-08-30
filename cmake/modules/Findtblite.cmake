# tblite -- the xTB engine behind GFN1 and GFN2.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME
  tblite
  GIT_REPOSITORY
  "https://github.com/tblite/tblite"
  GIT_TAG
  "v0.6.0"
  NAMESPACED_TARGET)
