# cuEST-backed integrals, grids and SCF. Built only when MQC_ENABLE_CUEST is ON.
#
# Listed in dependency order for readability; CMake derives the real Fortran
# compile order from the MODULE/USE graph.
target_sources(
  ${main_lib}
  PRIVATE mqc_cuest_runtime.f90
          mqc_cuest_context.f90
          mqc_cuest_basis.f90
          mqc_cuest_grid.f90
          mqc_cuest_functionals.f90
          mqc_cuest_integrals.f90
          mqc_diis_device.f90
          mqc_cuest_scf.f90
          mqc_cuest_atomic_guess.f90
          mqc_cuest_gradient.f90
          mqc_cuest_driver.f90
          mqc_cuest_bridge.f90)
