target_sources(
  ${main_lib} PRIVATE mqc_optimizer_types.f90 mqc_optimization_output.f90
                      mqc_geometry_optimizer.f90)

# The DL-FIND bridge has two forms. The real one, in backends/dlfind/, is
# compiled only with the backend; the stub here stands in otherwise, so
# mqc_geometry_optimizer resolves `mqc_dlfind_bridge` without a preprocessor
# conditional. The same arrangement as mqc_cuest_bridge and mqc_libcint_bridge
# under src/methods. Both defining that module at once would be a duplicate
# symbol, which is why this is an either/or and not an addition.
if(NOT MQC_ENABLE_DLFIND)
  target_sources(${main_lib} PRIVATE mqc_dlfind_bridge_stub.f90)
endif()
