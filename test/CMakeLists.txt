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
addtest(mqc_config_roundtrip)
addtest(mqc_json_config)
addtest(mqc_json_reader)
addtest(mqc_json_schema)
addtest(mqc_ecp_reader)
addtest(mqc_fraglist)
addtest(mqc_bond_perception)
addtest(mqc_fingerprint)
addtest(mqc_checkpoint)

# Parses basis_sets/*.json, which mqc_extract_basis_sets() unpacks from the BSE
# bundle. Only the CMake build does that, so the test compiles to a skip stub
# elsewhere (notably under fpm, where the files are absent and gitignored).
addtestpp(mqc_json_basis_reader)
target_compile_definitions(test_mqc_json_basis_reader
                           PRIVATE MQC_HAVE_BASIS_SETS)

# Only meaningful with the cuEST backend compiled in.
if(MQC_ENABLE_CUEST)
  addtestpp(mqc_cuest_spin)
  target_compile_definitions(test_mqc_cuest_spin PRIVATE MQC_WITH_CUEST)
endif()
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
addtest(mqc_config_adapter)
addtest(mqc_cli_parser)
addtest(mqc_split_lines)
addtest(mqc_result_types)
addtest(mqc_finite_differences)
addtest(mqc_diis)
addtest(mqc_lebedev)
addtest(mqc_dft_radial)
addtest(mqc_dft_partition)
addtest(mqc_dft_grid)
addtest(mqc_vibrational_analysis)
addtest(mqc_thermochemistry)
addtest(mqc_method_types)
addtest(mqc_method_placeholders)
addtest(mqc_work_queue)

# C-API / external-interface layer. These call the bind(C) entry points the way
# a C or Python caller would; none of them needs MPI or a QC backend.
addtest(mqc_validate)
addtest(mqc_capi_system)
addtest(mqc_capi_fraglist)
addtest(mqc_capi_session)
addtest(mqc_capi_run)

if(MQC_ENABLE_TBLITE)
  addtest(mqc_method_api)
  addtest(mqc_hessian_finite_diff)
endif()
