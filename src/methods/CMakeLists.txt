target_sources(
  ${main_lib}
  PRIVATE # Method implementations
          mqc_method_base.f90 mqc_method_config.f90 mqc_method_hf.F90
          mqc_method_dft.F90 mqc_method_mcscf.f90)

if(MQC_ENABLE_TBLITE)
  target_sources(${main_lib} PRIVATE mqc_method_xtb.f90)
endif()

if(MQC_ENABLE_CUEST)
  add_subdirectory(cuest)
endif()

# Factory must come after method types are defined
target_sources(${main_lib} PRIVATE mqc_method_factory.F90)
