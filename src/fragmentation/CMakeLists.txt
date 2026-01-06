target_sources(
  ${main_lib}
  PRIVATE mqc_mbe.f90
          mqc_combinatorics.f90
          mqc_fragment_lookup.f90
          mqc_gmbe_utils.f90
          mqc_frag_utils.f90
          mqc_mbe_fragment_distribution_scheme.F90
          mqc_gmbe_fragment_distribution_scheme.f90
          mqc_mbe_fragment_distribution_scheme_hessian.F90
          mqc_mbe_mpi_fragment_distribution_scheme.F90
          mqc_unfragmented_workflow.f90
          mqc_serial_fragment_processor.f90
          mqc_mbe_io.f90
          mqc_physical_fragment.f90)
