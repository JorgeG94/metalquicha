target_sources(
  ${main_lib}
  PRIVATE mqc_basis_file_reader.f90
          mqc_basis_reader.f90
          mqc_cgto.f90
          mqc_elements.f90
          mqc_cli_parser.f90
          mqc_geometry.f90
          mqc_xyz_reader.f90
          mqc_mpi_tags.f90
          mqc_chemistry_algorithms.f90
          mqc_mbe.f90
          mqc_physical_fragment.f90
          mqc_input_parser.f90
          mqc_driver.f90)
