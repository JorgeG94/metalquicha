macro("my_fetch_package" package url)
  # Optional third argument for git reference (branch, tag, or commit)
  set(_git_ref "${ARGV2}")
  if(NOT _git_ref)
    set(_git_ref "HEAD")
  endif()

  string(TOLOWER "${package}" _pkg_lc)
  string(TOUPPER "${package}" _pkg_uc)

  message(STATUS "Retrieving ${package} from ${url} (ref: ${_git_ref})")
  include(FetchContent)

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(FETCHCONTENT_QUIET FALSE)
  endif()

  FetchContent_Declare(
    "${_pkg_lc}"
    GIT_REPOSITORY "${url}"
    GIT_TAG "${_git_ref}")
  FetchContent_MakeAvailable("${_pkg_lc}")

  add_library("${package}::${package}" INTERFACE IMPORTED)
  target_link_libraries("${package}::${package}" INTERFACE "${package}")

  if(NOT EXISTS "${${_pkg_lc}_BINARY_DIR}/include")
    file(MAKE_DIRECTORY "${${_pkg_lc}_BINARY_DIR}/include")
  endif()

  unset(_pkg_lc)
  unset(_pkg_uc)
  unset(_git_ref)

  # sanity check
  if(NOT TARGET "${package}::${package}")
    message(FATAL_ERROR "Could not find dependency ${package}")
  endif()
endmacro()
