# The exchange-correlation quadrature: the grid, the atomic partition, and the
# functional specification a deck's `functional` parses into. Where a run's time
# actually goes -- 89 per cent of an LDA one, 98 of a meta-GGA.
target_sources(
  ${main_lib}
  PRIVATE mqc_lebedev_data.f90
          mqc_lebedev.f90
          mqc_dft_radial_data.f90
          mqc_dft_radial.f90
          mqc_dft_partition.f90
          mqc_dft_prune.f90
          mqc_dft_grid.f90
          mqc_xc_spec.f90)
