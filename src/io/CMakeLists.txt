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

# The real backend lives in backends/hdf5/, outside src/, because fpm globs src/
# and cannot choose between two files defining mqc_hdf5_checkpoint -- having
# both here is a hard error there whatever CMake is told. Only one of the two is
# ever compiled, the same arrangement as mqc_cuest_bridge_stub.
if(NOT MQC_ENABLE_HDF5)
  target_sources(${main_lib} PRIVATE mqc_hdf5_checkpoint_stub.f90
                                     mqc_hdf5_amplitudes_stub.f90)
endif()
