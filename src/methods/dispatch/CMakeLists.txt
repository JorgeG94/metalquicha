# The dispatch layer: what a deck's `method` resolves to, and the type every
# implementation extends. It decides which backend runs; it does not do physics.
target_sources(
  ${main_lib}
  PRIVATE mqc_method_base.f90 mqc_method_config.f90 mqc_method_hf.F90
          mqc_method_dft.F90 mqc_method_mcscf.f90 mqc_cuest_iface.f90)

if(MQC_ENABLE_TBLITE)
  target_sources(${main_lib} PRIVATE mqc_method_xtb.f90)
endif()

# Factory must come after method types are defined
target_sources(${main_lib} PRIVATE mqc_method_factory.F90)
