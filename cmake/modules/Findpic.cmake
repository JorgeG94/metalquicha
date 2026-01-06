set(_lib "pic")
set(_url "https://github.com/JorgeG94/pic/")

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

# Use the main branch
set(_rev "v0.2.0")
my_fetch_package("${_lib}" "${_url}" "${_rev}")

unset(_lib)
unset(_url)
unset(_rev)
