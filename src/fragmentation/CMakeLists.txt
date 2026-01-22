# Common types and utilities
target_sources(
  ${main_lib}
  PRIVATE common/mqc_physical_fragment.f90 common/mqc_fragment_lookup.f90
          common/mqc_combinatorics.f90 common/mqc_frag_utils.f90
          common/mqc_many_body_expansion.f90)

# Standard MBE
target_sources(
  ${main_lib}
  PRIVATE mbe/mqc_mbe.f90
          mbe/mqc_mbe_fragment_distribution_scheme.F90
          mbe/mqc_mbe_mpi_fragment_distribution_scheme.F90
          mbe/mqc_mbe_fragment_distribution_scheme_hessian.F90
          mbe/mqc_serial_fragment_processor.f90
          mbe/mqc_unfragmented_workflow.f90
          mbe/mqc_mbe_io.f90)

# Generalized MBE (overlapping fragments)
target_sources(
  ${main_lib} PRIVATE gmbe/mqc_gmbe_utils.f90
                      gmbe/mqc_gmbe_fragment_distribution_scheme.f90)
