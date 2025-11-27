set(_lib "test-drive")
set(_pkg "test-drive")
set(_url "https://github.com/JorgeG94/test-drive")
set(_rev "main")

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

my_fetch_package("${_lib}" "${_url}" "${_rev}")

unset(_lib)
unset(_pkg)
unset(_url)
unset(_rev)
