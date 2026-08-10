target_sources(
  ${main_lib}
  PRIVATE mqc_config_types.f90
          mqc_io_helpers.f90
          mqc_json.f90
          mqc_config_adapter.f90
          mqc_xyz_reader.f90
          mqc_json_schema.f90
          mqc_json_config_reader.f90
          mqc_json_writer.f90
          mqc_fragment_table_writer.f90
          mqc_checkpoint.f90)
