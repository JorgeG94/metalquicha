# The mqc side of terco: marshalling only. Built when MQC_ENABLE_TERCO is ON.
target_sources(${main_lib} PRIVATE mqc_terco_driver.f90 mqc_terco_bridge.f90)
