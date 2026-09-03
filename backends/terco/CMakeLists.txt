# terco device-resident SCF backend, kept out of src/ for the same reason the
# cuEST one is: fpm globs src/ and cannot link a binary library, so fpm builds
# the CPU/xTB path only. This directory is pulled in by MqcDependencies when
# MQC_ENABLE_TERCO is ON.
#
# `trc_c_interfaces.f90` is compiled from the terco checkout rather than copied
# here. That is the whole point of the arrangement: a copied interface block
# still compiles and still links after terco changes an argument, because a C
# symbol carries no signature, and fails at run time by reading whatever was in
# the register. Compiling terco's own declaration means a drifted ABI is a
# compile error in this build.
target_sources(${main_lib} PRIVATE ${TERCO_INTERFACES})

add_subdirectory(backend)
