target_sources(
  ${main_lib}
  PRIVATE mqc_basis_file_reader.f90
          mqc_config_parser.F90
          mqc_config_parser_basic_sections.f90
          mqc_config_parser_structure.f90
          mqc_config_parser_fragments.f90
          mqc_config_parser_calc_settings.f90
          mqc_config_parser_molecules.f90
          mqc_io_helpers.f90
          mqc_json.f90
          mqc_config_adapter.f90
          mqc_xyz_reader.f90
          mqc_json_writer.f90)
