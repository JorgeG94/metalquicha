set(_lib "pic-mpi")
set(_url "https://github.com/JorgeG94/pic-mpi/")

# Pass PIC_USE_LEGACY_MPI option to the fetched package if set
if(DEFINED PIC_USE_LEGACY_MPI)
  set(PIC_USE_LEGACY_MPI
      ${PIC_USE_LEGACY_MPI}
      CACHE BOOL "Use legacy MPI module" FORCE)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

# Use the main branch
set(_rev "v0.1.0")
my_fetch_package("${_lib}" "${_url}" "${_rev}")

unset(_lib)
unset(_url)
unset(_rev)
