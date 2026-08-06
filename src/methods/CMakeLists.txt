target_sources(
  ${main_lib}
  PRIVATE # Method implementations
          mqc_method_base.f90
          mqc_semi_numerical_hessian.f90
          mqc_method_config.f90
          mqc_method_hf.F90
          mqc_method_dft.F90
          mqc_method_mcscf.f90
          mqc_cuest_iface.f90)

# The cuEST bridge has two forms. The real one, in cuest_backend/, is compiled
# only with the backend; the stub here stands in otherwise, so the method files
# always resolve `mqc_cuest_bridge` without a preprocessor conditional.
if(NOT MQC_ENABLE_CUEST)
  target_sources(${main_lib} PRIVATE mqc_cuest_bridge_stub.f90)
endif()

if(MQC_ENABLE_TBLITE)
  target_sources(${main_lib} PRIVATE mqc_method_xtb.f90)
endif()

# Factory must come after method types are defined
target_sources(${main_lib} PRIVATE mqc_method_factory.F90)
