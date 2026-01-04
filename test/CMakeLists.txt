# Unit testing

find_package(Threads REQUIRED)

# Macro to add a test executable for .f90 files
macro(ADDTEST name)
  add_executable(test_${name} test_${name}.f90)
  set_target_properties(
    test_${name} PROPERTIES Fortran_MODULE_DIRECTORY
                            "${CMAKE_BINARY_DIR}/tests/modules/${name}")
  target_include_directories(test_${name} PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(test_${name} PRIVATE "${main_lib}"
                                             "test-drive::test-drive")
  add_test(
    NAME ${project_name}/${name}
    COMMAND $<TARGET_FILE:test_${name}>
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADDTEST)

# Macro to add a test executable for preprocessed .F90 files
macro(ADDTESTPP name)
  add_executable(test_${name} test_${name}.F90)
  set_target_properties(
    test_${name} PROPERTIES Fortran_MODULE_DIRECTORY
                            "${CMAKE_BINARY_DIR}/tests/modules/${name}")
  target_include_directories(test_${name} PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(test_${name} PRIVATE "${main_lib}"
                                             "test-drive::test-drive")
  add_test(
    NAME ${project_name}/${name}
    COMMAND $<TARGET_FILE:test_${name}>
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADDTESTPP)

# Add individual test executables
message(STATUS "Adding individual test executables for ${project_name}...")

# Non-MPI tests
addtest(mqc_xyz_reader)
addtest(mqc_xyz_reader_detailed)
addtest(mqc_basis_reader)
addtest(mqc_basis_file_reader)
addtest(mqc_basis_utils)
addtest(mqc_elements)
addtest(mqc_cgto)
addtest(mqc_physical_fragment)
addtest(mqc_frag_utils)
addtest(mqc_gmbe_intersection)
addtest(mqc_fragment_from_atoms)
addtest(mqc_gmbe_energy)
addtest(mqc_mbe)
addtest(mqc_geometry)
addtest(mqc_config_parser)
addtest(mqc_config_adapter)
addtest(mqc_cli_parser)
addtest(mqc_split_lines)
addtest(mqc_result_types)
addtest(mqc_finite_differences)

if(MQC_ENABLE_TBLITE)
  addtest(mqc_method_api)
  addtest(mqc_hessian_finite_diff)
endif()
