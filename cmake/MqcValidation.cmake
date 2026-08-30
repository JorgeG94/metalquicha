# One way to register a validation program.
#
# These are the `check_*` and `bench_*` programs under validation/: standalone
# executables that compare against a second code, or time something. Forty-seven
# of them were spelled out three lines at a time in the top-level file, followed
# by a foreach listing forty-seven names again to mark them EXCLUDE_FROM_ALL --
# a list that had to be kept in step with the definitions by hand, and quietly
# was not.
#
# cmake-format: off
#   mqc_add_validation_program(<name>
#     [SOURCE               <file>]     # default validation/<name>.f90
#     [LIBRARIES            <lib>...]   # beyond ${main_lib}
#     [INCLUDE_DIRECTORIES  <dir>...]   # beyond the module directory
#     [IN_DEFAULT_BUILD])               # a ctest case drives it
# cmake-format: on
#
# **EXCLUDE_FROM_ALL by default.** They are one-off checks, and several predate
# things they still name: the reference is another program's output rather than
# anything in this tree, so they cannot run in CI anyway. Building them by
# default bought nothing and cost a broken build every time one of them drifted
# out of step with a guard it depended on -- which is how a libcint-off job
# failed on a file nobody had asked for.
#
# EXCLUDE_FROM_ALL rather than deleting them: each is still a target, so
#
# cmake --build build --target check_rhf
#
# builds exactly the one you want, and the ones still worth keeping can be
# converted into test-drive tests -- where CI does run them -- one at a time
# rather than in a single sweep. The cost to know about: a plain `cmake --build
# build` does not rebuild these, so an out-of-date one sits on disk printing
# whatever it printed last time.
#
# `IN_DEFAULT_BUILD` is for the three that `add_test` drives, where excluding
# them would leave ctest with entries and no executables to run.
#
# The output directory is forced back to the top of the build tree. These used
# to be defined in the top-level file and land in `build/`; registering them
# from a subdirectory would otherwise scatter them, and `run_validation.py`,
# every doc page and everyone's fingers expect `build/check_thing`.
include_guard(GLOBAL)

function(mqc_add_validation_program name)
  cmake_parse_arguments(MV "IN_DEFAULT_BUILD" "SOURCE"
                        "LIBRARIES;INCLUDE_DIRECTORIES" ${ARGN})
  if(MV_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "mqc_add_validation_program(${name}): unrecognised "
                        "arguments '${MV_UNPARSED_ARGUMENTS}'")
  endif()

  if(NOT MV_SOURCE)
    set(MV_SOURCE "${PROJECT_SOURCE_DIR}/validation/${name}.f90")
  endif()
  if(NOT EXISTS "${MV_SOURCE}")
    message(FATAL_ERROR "mqc_add_validation_program(${name}): no such source "
                        "'${MV_SOURCE}'")
  endif()

  add_executable(${name} "${MV_SOURCE}")
  target_include_directories(${name} PRIVATE "${CMAKE_BINARY_DIR}/modules"
                                             ${MV_INCLUDE_DIRECTORIES})
  target_link_libraries(${name} PRIVATE ${main_lib} ${MV_LIBRARIES})
  set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                           "${CMAKE_BINARY_DIR}")
  if(NOT MV_IN_DEFAULT_BUILD)
    set_target_properties(${name} PROPERTIES EXCLUDE_FROM_ALL TRUE)
  endif()
endfunction()
