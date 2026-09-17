# pic-blas -- the BLAS and LAPACK wrappers.
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME
  pic-blas
  GIT_REPOSITORY
  "https://github.com/JorgeG94/pic-blas/"
  GIT_TAG
  "main"
  NAMESPACED_TARGET)
