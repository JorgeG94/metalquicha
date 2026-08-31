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
#
# cublas.f90 and cusolver.f90 are the exceptions: hand-written, covering only
# the handful of BLAS calls the SCF needs and the single dense eigensolver it
# needs. They sit here because they are the same kind of thing -- thin
# iso_c_binding layers over a vendor C API -- and share the convention of
# binding device pointers as TYPE(c_ptr).
target_sources(${main_lib} PRIVATE cuest.f90 cuest_helpers.f90 cuda_runtime.f90
                                   cuda_helpers.f90 cublas.f90 cusolver.f90)

# cuda_runtime.f90 carries one #if on MQC_CUDA_VERSION_MAJOR, for the entry
# point CUDA 13 renamed out from under it. Preprocessing is asked for on the
# file rather than by renaming it to .F90 -- the convention everywhere else in
# this project -- so the name still matches the generated original and a
# re-copy from mod_cuest stays a plain overwrite.
set_source_files_properties(cuda_runtime.f90 PROPERTIES Fortran_PREPROCESS ON)
