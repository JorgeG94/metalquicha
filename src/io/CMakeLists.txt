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

# Compiled only when HDF5 is on. The module declares the library's globals as
# bind(C) variables, which the linker treats as definitions to be merged with
# the real ones -- so pulling this object into a build without libhdf5 leaves
# H5open undefined. That does not fail the link of a shared library, it fails
# the dlopen, which is to say it fails in Python rather than in CMake.
if(MQC_ENABLE_HDF5)
  target_sources(${main_lib} PRIVATE mqc_hdf5_bindings.f90)
endif()
