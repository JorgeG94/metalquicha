# The validation programs that do not need a particular integrals backend.
#
# The cenzontle ones live with the backend they exercise, in
# backends/cenzontle/validation.cmake, and are registered from there so that a
# build without that backend never mentions them.
#
# All of these are EXCLUDE_FROM_ALL -- see cmake/MqcValidation.cmake for why
# that matters more than it looks.

# The grid, and the three pieces it is built from. No backend: these are this
# project's own quadrature.
mqc_add_validation_program(check_grid)
mqc_add_validation_program(check_partition)
mqc_add_validation_program(check_radial)
mqc_add_validation_program(check_lebedev)

# Driving the program the way a user does, through its own entry points.
mqc_add_validation_program(check_run)
mqc_add_validation_program(check_session)
# Manual check that a session splits rank 0 from the workers. Not a ctest case:
# it is only meaningful under mpirun with more than one rank, and ctest here has
# no MPI harness.
mqc_add_validation_program(check_session_run)
mqc_add_validation_program(check_bcast)
mqc_add_validation_program(check_bcast_system)

# Bond orders come from an xTB single point, so this needs tblite, not libcint.
# It was registered inside the cenzontle block once, so a cenzontle-on
# tblite-off build compiled it against an mqc_method_xtb that CMake had never
# built.
if(MQC_ENABLE_TBLITE)
  mqc_add_validation_program(check_bond_orders)
endif()

# libxc on its own -- the functional values, with no integrals involved.
if(MQC_ENABLE_LIBXC)
  mqc_add_validation_program(check_libxc LIBRARIES xcf03 xc INCLUDE_DIRECTORIES
                             ${libxc_BINARY_DIR})
endif()

if(MQC_ENABLE_HDF5)
  mqc_add_validation_program(check_hdf5)
  # Writes the file that validation/check_amplitude_layout.py reads. Split in
  # two because the check cannot be made from Fortran: a block stored under the
  # wrong declared shape round trips through our own bindings exactly.
  mqc_add_validation_program(check_amplitude_layout)
endif()

# The seam the CREST integration rests on, checked on its own because every way
# it can fail is quiet: an uncalled pointer returns zero rather than an error.
if(MQC_ENABLE_CREST)
  mqc_add_validation_program(check_crest_callback)
  mqc_add_validation_program(check_crest_search)
  mqc_add_validation_program(check_crest_refine)
endif()

# Finite-difference check on the analytic gradient. Not a ctest case: it needs a
# GPU of sm_80 or newer, so it cannot run in CI or on a login node.
if(MQC_ENABLE_CUEST)
  mqc_add_validation_program(check_gradient)
endif()

# The ECP layout and the ECP matrix, against PySCF. Not a ctest case yet: the
# reference comes from PySCF rather than from a number this program printed
# earlier, so it belongs with the other validation programs that need a second
# code to compare against. libfint rather than libcint -- the ECP integrals are
# libfint's and libcint has none, so this cannot link in a -DMQC_USE_LIBFINT=OFF
# build.
if(MQC_ENABLE_CZT AND MQC_USE_LIBFINT)
  mqc_add_validation_program(check_ecp)
endif()
