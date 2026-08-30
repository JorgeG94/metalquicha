# Installing, exporting, and the coverage target.
#
# Last thing the top-level does: everything here needs every target to exist.
include_guard(GLOBAL)

foreach(tgt ${all_targets})
  message(STATUS "${tgt}")
  install(
    TARGETS ${tgt}
    EXPORT ${project_name}Targets
    ARCHIVE DESTINATION lib
    INCLUDES
    DESTINATION include/${project_name})
endforeach()

install(
  DIRECTORY ${CMAKE_BINARY_DIR}/modules/
  DESTINATION include/${project_name}
  FILES_MATCHING
  PATTERN "*.mod")

install(
  EXPORT ${project_name}Targets
  FILE ${project_name}Targets.cmake
  NAMESPACE ${project_name}::
  DESTINATION lib/cmake/${project_name})

include(CMakePackageConfigHelpers)

# RENAME YOUR sampleConfig.cmake.in to match ${project_name}Config.cmake
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${project_name}Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${project_name}Config.cmake"
  INSTALL_DESTINATION lib/cmake/${project_name})

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${project_name}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${project_name}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${project_name}ConfigVersion.cmake"
        DESTINATION lib/cmake/${project_name})

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/${project_name}Config.cmake"
  DESTINATION lib/cmake/${project_name}
  RENAME ${project_name}-config.cmake)

print_summary()

if(CMAKE_BUILD_TYPE STREQUAL "Coverage-mqc")
  # Both manifests, every time. A case naming a backend this build does not have
  # is skipped by the runner rather than failed -- that is what the `requires`
  # markers are for -- so the xTB manifest costs nothing without tblite and the
  # CPU one costs nothing without libcint. Running only the manifest someone
  # remembered to enable is how half the program stopped being measured in the
  # first place.
  #
  # The CPU manifest is where the driver, the fragmentation schemes and the JSON
  # writer get exercised: those are reached by *running* the program rather than
  # by calling into it, so no unit test touches them, and the coverage report
  # showed them near zero until this ran.
  add_custom_target(
    coverage
    COMMAND ${CMAKE_CTEST_COMMAND} -R "mqc" --output-on-failure
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/validation python3
            run_validation.py
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/validation python3
            run_validation.py --mpi --np 2
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/validation python3
            run_validation.py --manifest validation_tests_cpu.json
    # And the Python interface, which is the only thing that drives the
    # program through its C entry points. `src/parallel` is 0% without this
    # and it is not untested code: `mqc_bcast` and `mqc_bcast_system` are used
    # by `mqc_session` and by nothing else, so no deck reaches them however
    # many ranks it runs on. The examples assert against PySCF references and
    # skip whatever this build cannot do, so they are a test rather than a
    # warm-up lap.
    COMMAND
      ${CMAKE_COMMAND} -E env MQC_LIBRARY=$<TARGET_FILE:mqc_shared>
      PYTHONPATH=${CMAKE_SOURCE_DIR}/python python3
      ${CMAKE_SOURCE_DIR}/python/examples/backends.py
    COMMAND
      ${CMAKE_COMMAND} -E env MQC_LIBRARY=$<TARGET_FILE:mqc_shared>
      PYTHONPATH=${CMAKE_SOURCE_DIR}/python python3
      ${CMAKE_SOURCE_DIR}/python/examples/methods.py
    COMMAND ${CMAKE_COMMAND} -E remove -f coverage.info
    # `negative` because gcov's counters are not atomic and this program is
    # threaded: two threads incrementing the same counter can leave it below
    # zero, which lcov 1.x warned about and lcov 2.x -- what the runner
    # installs -- refuses outright, failing the whole target over a count
    # rather than over coverage. The alternative is -fprofile-update=atomic,
    # which costs time in a build that is already fifty times slower than
    # release. A raced line was still executed, so the report is not the
    # thing in doubt here.
    COMMAND lcov --directory . --capture --ignore-errors negative --output-file
            coverage.info
    COMMAND lcov --ignore-errors unused --remove coverage.info '/usr/*'
            '*/build/*' '*/validation/*' --output-file coverage_filtered.info
    COMMAND genhtml coverage_filtered.info --output-directory coverage_report
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating code coverage report")
  # The examples load the shared library through ctypes, so it has to exist
  # before the target runs rather than being built by whoever remembered.
  add_dependencies(coverage mqc_shared)
endif()
