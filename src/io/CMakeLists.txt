target_sources(
  ${main_lib}
  PRIVATE mqc_basis_file_reader.f90
          mqc_config_parser.F90
          mqc_io_helpers.f90
          mqc_json.f90
          mqc_config_adapter.f90
          mqc_xyz_reader.f90
          mqc_json_writer.f90)
add_subdirectory(parsers)
