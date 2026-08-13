# based on finding tblite inside xtb

set(_lib "tblite")
set(_pkg "TBLITE")
set(_url "https://github.com/tblite/tblite")
set(_rev "v0.6.0")

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

my_fetch_package("${_lib}" "${_url}" "${_rev}")

unset(_lib)
unset(_pkg)
unset(_url)
unset(_rev)
