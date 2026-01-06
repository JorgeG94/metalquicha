set(_lib "pic-blas")
set(_url "https://github.com/JorgeG94/pic-blas/")

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

# Use the main branch
set(_rev "v0.1.0")
my_fetch_package("${_lib}" "${_url}" "${_rev}")

unset(_lib)
unset(_url)
unset(_rev)
