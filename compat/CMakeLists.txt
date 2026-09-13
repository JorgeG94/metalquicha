# Compiled only for MQC_ENABLE_SERIAL, and only ever reached from there: the
# module inside defines `omp_lib`, which would shadow the compiler's own in any
# build that actually wants OpenMP. MqcDependencies.cmake adds this directory
# under that condition and no other. See omp_lib_serial.f90 for the why.
target_sources(${main_lib} PRIVATE omp_lib_serial.f90)
