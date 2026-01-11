set(_lib "jsonfortran")
set(_pkg "JSONFORTRAN")
set(_url "https://github.com/jacobwilliams/json-fortran.git")
set(_rev "9.2.0")

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

my_fetch_package("${_lib}" "${_url}" "${_rev}")

unset(_lib)
unset(_pkg)
unset(_url)
unset(_rev)
