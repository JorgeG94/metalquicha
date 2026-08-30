# json-fortran -- reads every input deck and writes every result file.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME
  jsonfortran
  GIT_REPOSITORY
  "https://github.com/jacobwilliams/json-fortran.git"
  GIT_TAG
  "9.3.1"
  NAMESPACED_TARGET)
