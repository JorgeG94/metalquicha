target_sources(${main_lib} PRIVATE # Method implementations
                                   mqc_method_base.f90 mqc_method_hf.f90)

if(MQC_ENABLE_TBLITE)
  target_sources(${main_lib} PRIVATE mqc_method_xtb.f90)
endif()
