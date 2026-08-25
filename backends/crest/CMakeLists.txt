# Conformer sampling through CREST, kept out of src/ so the fpm build never sees
# it.
#
# fpm globs src/ and has no way to choose between two files defining the same
# module, so mqc_crest_driver.f90 and mqc_crest_driver_stub.f90 sitting side by
# side there is a hard error regardless of any CMake option. Only the stub stays
# in src/, and fpm therefore builds a program that declines conformer sampling,
# the same arrangement as mqc_hdf5_checkpoint_stub and mqc_cuest_bridge_stub.
#
# Pulled in by the top-level CMakeLists when MQC_ENABLE_CREST is ON.
target_sources(${main_lib} PRIVATE mqc_crest_driver.f90)
