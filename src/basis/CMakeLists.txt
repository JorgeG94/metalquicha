target_sources(
  ${main_lib}
  PRIVATE mqc_basis_utils.F90 mqc_cgto.f90 mqc_ecp.f90
          mqc_json_basis_reader.f90 mqc_json_ecp_reader.f90
          mqc_basis_normalization.f90)

# Where find_basis_file falls back to when neither $MQC_BASIS_PATH nor
# ./basis_sets has the file -- the tree this build was configured from, which is
# what lets `mqc` run from any directory.
#
# On the target, not on the source file: SOURCE properties are scoped to the
# directory that sets them, and ${main_lib} is defined in the top-level
# CMakeLists, so a source-scoped definition here would never reach the compiler.
# Only mqc_basis_utils.F90 reads the macro.
target_compile_definitions(
  ${main_lib}
  PRIVATE "MQC_DEFAULT_BASIS_DIR=\"${CMAKE_SOURCE_DIR}/basis_sets\"")
