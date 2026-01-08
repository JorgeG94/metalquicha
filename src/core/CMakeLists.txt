target_sources(
  ${main_lib}
  PRIVATE # Core modules
          mqc_calculation_keywords.f90 mqc_elements.f90 mqc_geometry.f90
          mqc_logo.f90 mqc_mpi_context.f90 mqc_result_types.f90)
