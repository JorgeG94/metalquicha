# One way to bring a dependency in.
#
# Every `Find<pkg>.cmake` beside this file calls `mqc_fetch` and nothing else,
# so a dependency is declared in exactly one place and every declaration reads
# the same. What differs between packages -- a git tag against a release
# tarball, a namespaced target against a bare one -- is a keyword here rather
# than a different shape of file.
#
# cmake-format: off
#   mqc_fetch(
#     NAME               <pkg>            # FetchContent name, lower case
#     [GIT_REPOSITORY    <url>]
#     [GIT_TAG           <tag|sha>]       # a tag or a SHA, never a branch
#     [GIT_SUBMODULES    <path>...]       # "" checks out none
#     [URL               <archive url>]
#     [URL_HASH          <ALGO=hex>]
#     [NAMESPACED_TARGET]                 # also define <pkg>::<pkg>
#     [PROVIDES          <target>...])    # assert these exist afterwards
# cmake-format: on
#
# `PROVIDES` is the part worth having. A fetched project that renames a target
# otherwise fails at link time, in another directory, with a message naming
# neither the dependency nor the pin. Asserting here fails during configure and
# says which package moved.
#
# A macro rather than a function, and deliberately: `FetchContent_MakeAvailable`
# configures the dependency in the calling scope, so a function would swallow
# every variable the subproject sets at its top level. The `_mqc_fetch_*` locals
# are unset at the end for the same reason -- in a macro they are the caller's.
include_guard(GLOBAL)

include(FetchContent)

macro(mqc_fetch)
  cmake_parse_arguments(
    _mqc_fetch "NAMESPACED_TARGET" "NAME;GIT_REPOSITORY;GIT_TAG;URL;URL_HASH"
    "PROVIDES;GIT_SUBMODULES" ${ARGN})

  if(NOT _mqc_fetch_NAME)
    message(FATAL_ERROR "mqc_fetch: NAME is required")
  endif()
  if(_mqc_fetch_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "mqc_fetch(${_mqc_fetch_NAME}): unrecognised arguments "
                        "'${_mqc_fetch_UNPARSED_ARGUMENTS}'")
  endif()

  string(TOLOWER "${_mqc_fetch_NAME}" _mqc_fetch_lc)

  # A single provenance line per dependency, so a configure log says what was
  # fetched and from which pin without reading any of these files.
  if(_mqc_fetch_GIT_REPOSITORY)
    if(NOT _mqc_fetch_GIT_TAG)
      set(_mqc_fetch_GIT_TAG "HEAD")
    endif()
    message(
      STATUS
        "Dependency ${_mqc_fetch_NAME}: ${_mqc_fetch_GIT_REPOSITORY}@${_mqc_fetch_GIT_TAG}"
    )
    if(DEFINED _mqc_fetch_GIT_SUBMODULES)
      FetchContent_Declare(
        "${_mqc_fetch_lc}"
        GIT_REPOSITORY "${_mqc_fetch_GIT_REPOSITORY}"
        GIT_TAG "${_mqc_fetch_GIT_TAG}"
        GIT_SUBMODULES "${_mqc_fetch_GIT_SUBMODULES}")
    else()
      FetchContent_Declare(
        "${_mqc_fetch_lc}"
        GIT_REPOSITORY "${_mqc_fetch_GIT_REPOSITORY}"
        GIT_TAG "${_mqc_fetch_GIT_TAG}")
    endif()
  elseif(_mqc_fetch_URL)
    message(STATUS "Dependency ${_mqc_fetch_NAME}: ${_mqc_fetch_URL}")
    if(_mqc_fetch_URL_HASH)
      FetchContent_Declare(
        "${_mqc_fetch_lc}"
        URL "${_mqc_fetch_URL}"
        URL_HASH "${_mqc_fetch_URL_HASH}"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    else()
      FetchContent_Declare(
        "${_mqc_fetch_lc}" URL "${_mqc_fetch_URL}" DOWNLOAD_EXTRACT_TIMESTAMP
                               TRUE)
    endif()
  else()
    message(
      FATAL_ERROR "mqc_fetch(${_mqc_fetch_NAME}): needs GIT_REPOSITORY or URL")
  endif()

  # Progress output only where someone is watching a Debug configure.
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(FETCHCONTENT_QUIET FALSE)
  endif()

  FetchContent_MakeAvailable("${_mqc_fetch_lc}")

  if(_mqc_fetch_NAMESPACED_TARGET)
    if(NOT TARGET "${_mqc_fetch_NAME}::${_mqc_fetch_NAME}")
      add_library("${_mqc_fetch_NAME}::${_mqc_fetch_NAME}" INTERFACE IMPORTED)
      target_link_libraries("${_mqc_fetch_NAME}::${_mqc_fetch_NAME}"
                            INTERFACE "${_mqc_fetch_NAME}")
    endif()
    # Some of these projects name an include directory they never create.
    if(NOT EXISTS "${${_mqc_fetch_lc}_BINARY_DIR}/include")
      file(MAKE_DIRECTORY "${${_mqc_fetch_lc}_BINARY_DIR}/include")
    endif()
    list(APPEND _mqc_fetch_PROVIDES "${_mqc_fetch_NAME}::${_mqc_fetch_NAME}")
  endif()

  foreach(_mqc_fetch_target IN LISTS _mqc_fetch_PROVIDES)
    if(NOT TARGET "${_mqc_fetch_target}")
      message(
        FATAL_ERROR
          "${_mqc_fetch_NAME} was fetched but defines no target "
          "'${_mqc_fetch_target}'. The pin in "
          "cmake/modules/Find${_mqc_fetch_NAME}.cmake probably needs updating, "
          "or the project renamed it.")
    endif()
  endforeach()

  unset(_mqc_fetch_NAME)
  unset(_mqc_fetch_GIT_REPOSITORY)
  unset(_mqc_fetch_GIT_TAG)
  unset(_mqc_fetch_GIT_SUBMODULES)
  unset(_mqc_fetch_URL)
  unset(_mqc_fetch_URL_HASH)
  unset(_mqc_fetch_PROVIDES)
  unset(_mqc_fetch_NAMESPACED_TARGET)
  unset(_mqc_fetch_UNPARSED_ARGUMENTS)
  unset(_mqc_fetch_KEYWORDS_MISSING_VALUES)
  unset(_mqc_fetch_lc)
  unset(_mqc_fetch_target)
endmacro()
