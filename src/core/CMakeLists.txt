target_sources(
  ${main_lib}
  PRIVATE # Core modules
          mqc_calculation_keywords.f90
          mqc_physical_constants.f90
          mqc_program_limits.f90
          mqc_elements.f90
          mqc_geometry.f90
          mqc_logo.f90
          mqc_mpi_comms.f90
          mqc_resources.f90
          mqc_result_types.f90)
