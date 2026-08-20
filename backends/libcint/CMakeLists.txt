# CPU Gaussian-integral backend, kept out of src/ so the fpm build never sees it
# -- the same reason backends/cuest/ and backends/hdf5/ live here.
#
# These modules `use libcint_fortran`, which only exists once CMake has fetched
# libcint. fpm globs src/ and has no such dependency declared, so leaving them
# there fails the fpm build outright:
#
# <ERROR> *cmd_build* Target error: Unable to find source for module dependency:
# "libcint_fortran" used by "src/methods/mqc_libcint_integrals.f90"
#
# and there is no CMake option that helps, because fpm never reads CMake. So fpm
# builds the xTB path and CMake adds this directory when MQC_ENABLE_LIBCINT is
# ON.
target_sources(
  ${main_lib}
  PRIVATE mqc_libcint_atomic_guess.f90
          mqc_libcint_gradient.f90
          mqc_libcint_response.f90
          mqc_libcint_bridge.f90
          mqc_libcint_direct.f90
          mqc_libcint_integrals.F90
          mqc_libcint_multipole.f90
          mqc_libcint_avas.f90
          mqc_libcint_bonding.f90
          mqc_libcint_ieda.f90
          mqc_libcint_quao.f90
          mqc_libcint_quao_report.f90
          mqc_libcint_casci.f90
          mqc_libcint_mcscf.f90
          mqc_libcint_localize.f90
          mqc_libcint_cphf.F90
          mqc_libcint_dma.f90
          mqc_libcint_esp.f90
          mqc_libcint_charges.f90
          mqc_libcint_fukui.f90
          mqc_libcint_fmo.f90
          mqc_libcint_screening.f90
          mqc_libcint_projection.f90
          mqc_libcint_rhf.f90
          mqc_libcint_mp2.f90
          mqc_libcint_mp2_gradient.f90
          mqc_libcint_ri_mp2_gradient.f90
          mqc_libcint_cc.f90
          mqc_libcint_rcc.f90
          mqc_libcint_ao_data.f90
          mqc_libcint_ao.f90
          mqc_libcint_xc.F90
          mqc_efp_potential.f90
          mqc_efp_read.f90
          mqc_efp_interaction.f90
          mqc_efp_pair.f90
          mqc_efp_rotate.f90
          mqc_efp_energy.f90
          mqc_sapt.f90)

# The analytic Hessian is libfint's or it is nothing. Its second derivatives
# come from `cint_gen_hess` and the workspace type beside it, which a libcint
# build has no module files for -- so that build compiles a stub that reports
# the absence, and `mqc_method_hf` falls back to finite differences exactly as
# it does for an unrestricted or density-fitted reference.
#
# Split at the file rather than with `#ifdef` because the whole of both modules
# is libfint-only: preprocessing them would leave two files that are one guard
# and nothing else in half the configurations.
if(MQC_USE_LIBFINT)
  target_sources(${main_lib} PRIVATE mqc_libcint_hess_ints.f90
                                     mqc_libcint_hessian.f90)
else()
  target_sources(${main_lib} PRIVATE mqc_libcint_hessian_stub.f90)
endif()
