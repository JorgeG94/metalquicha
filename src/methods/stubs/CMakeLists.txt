# Each of these defines a module name that a real backend also defines, and
# exactly one of the pair is ever compiled. That is what lets the dispatch layer
# `use mqc_cuest_bridge` with no preprocessor guard at the call site, and it is
# also why the real ones cannot live under src/: fpm globs src/ and would reject
# the duplicate module.
if(NOT MQC_ENABLE_CUEST)
  target_sources(${main_lib} PRIVATE mqc_cuest_bridge_stub.f90)
endif()

if(NOT MQC_ENABLE_CZT)
  target_sources(${main_lib} PRIVATE mqc_czt_bridge_stub.f90)
endif()

if(NOT MQC_ENABLE_TERCO)
  target_sources(${main_lib} PRIVATE mqc_terco_bridge_stub.f90)
endif()

if(NOT MQC_ENABLE_CREST)
  target_sources(${main_lib} PRIVATE mqc_crest_driver_stub.f90)
endif()
