# Vendored Fortran bindings to NVIDIA cuEST and the CUDA runtime.
#
# GENERATED SOURCES -- do not hand-edit. cuest.f90 and cuda_runtime.f90 are
# machine-produced from the cuEST headers and the CUDA headers respectively by
# the generators in the upstream mod_cuest project. To pick up a new cuEST
# release, regenerate there and re-copy, rather than patching here.
#
# These are copied in rather than fetched so the build has one fewer moving part
# while the integration is young; switching to FetchContent later means
# replacing this file, not touching any metalquicha source.
#
# They do not follow the mqc_ naming convention or the fortitude `only` rule for
# the same reason -- they are third-party output, not project code.
target_sources(${main_lib} PRIVATE cuest.f90 cuest_helpers.f90 cuda_runtime.f90
                                   cuda_helpers.f90)
