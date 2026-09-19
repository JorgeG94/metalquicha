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
          mqc_xc_spec.f90
          mqc_dispersion_names.f90)

# The dispersion correction, in exactly one of its two forms. The real wrapper
# is backends/dftd3/mqc_dispersion.f90 and defines the same module; this stub
# refuses, naming the option. The name mapping above is compiled either way --
# it is a fact about two vocabularies, not about a link line, and it is what the
# unit tests can reach without the library.
if(NOT MQC_ENABLE_DFTD3)
  target_sources(${main_lib} PRIVATE mqc_dispersion_stub.f90)
endif()
