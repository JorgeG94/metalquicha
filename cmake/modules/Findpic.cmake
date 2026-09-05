# pic -- sorting, logging, timers.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME
  pic
  GIT_REPOSITORY
  "https://github.com/JorgeG94/pic/"
  GIT_TAG
  "v0.6.1"
  NAMESPACED_TARGET)
