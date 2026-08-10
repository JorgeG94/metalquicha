# HDF5 checkpoint backend, kept out of src/ so the fpm build never sees it.
#
# fpm globs src/ and has no way to choose between two files defining the same
# module, so mqc_hdf5_checkpoint.f90 and mqc_hdf5_checkpoint_stub.f90 sitting
# side by side there is a hard error -- "one or more duplicate module names
# found" -- regardless of any CMake option. Only the stub stays in src/, and
# fpm therefore builds the text-checkpoint path, exactly as it builds the
# CPU/xTB path and not cuEST.
#
# mqc_hdf5_bindings.f90 moves for a second reason. It declares the library's
# globals as bind(C) variables, which the linker treats as definitions to merge
# with the real ones, so compiling it without libhdf5 leaves H5open undefined.
# That does not fail the link of a shared library; it fails the dlopen, which
# is to say it fails in Python rather than in CMake.
#
# Pulled in by the top-level CMakeLists when MQC_ENABLE_HDF5 is ON.
target_sources(${main_lib} PRIVATE mqc_hdf5_bindings.f90 mqc_hdf5_checkpoint.f90)
