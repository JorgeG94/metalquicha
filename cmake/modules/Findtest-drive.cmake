# test-drive -- the unit test framework.
#
# `main` rather than a tag, and the one dependency here where that is right: it
# is a test harness, so an upstream change reddens the suite and nothing that
# ships.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME
  test-drive
  GIT_REPOSITORY
  "https://github.com/fortran-lang/test-drive"
  GIT_TAG
  "v0.6.1"
  NAMESPACED_TARGET)
