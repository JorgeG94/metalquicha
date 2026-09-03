# CPU Gaussian-integral backend, kept out of src/ so the fpm build never sees it
# -- the same reason backends/cuest/ and backends/hdf5/ live here.
#
# These modules `use libcint_fortran`, which only exists once CMake has fetched
# libcint. fpm globs src/ and has no such dependency declared, so leaving them
# there fails the fpm build outright:
#
# <ERROR> *cmd_build* Target error: Unable to find source for module dependency:
# "libcint_fortran" used by "src/methods/mqc_czt_integrals.f90"
#
# and there is no CMake option that helps, because fpm never reads CMake. So fpm
# builds the xTB path and CMake adds this directory when MQC_ENABLE_CZT is ON.
target_sources(
  ${main_lib}
  PRIVATE mqc_czt_afo.f90
          mqc_czt_atomic_guess.f90
          mqc_czt_gradient.f90
          mqc_czt_response.f90
          mqc_czt_hess_abi.f90
          mqc_czt_hess_ints.f90
          mqc_czt_hessian.f90
          mqc_czt_xc_hessian.f90
          mqc_czt_bridge.f90
          mqc_czt_direct.f90
          mqc_czt_ecp.F90
          mqc_czt_integrals.F90
          mqc_czt_multipole.f90
          mqc_czt_avas.f90
          mqc_czt_bonding.f90
          mqc_czt_ieda.f90
          mqc_czt_quao.f90
          mqc_czt_quao_report.f90
          mqc_czt_casci.f90
          mqc_czt_mcscf.f90
          mqc_czt_localize.f90
          mqc_czt_cphf.F90
          mqc_czt_dma.f90
          mqc_czt_esp.f90
          mqc_czt_pcm.f90
          mqc_czt_charges.f90
          mqc_czt_fukui.f90
          mqc_czt_fmo.f90
          mqc_czt_screening.f90
          mqc_czt_projection.f90
          mqc_czt_rhf.f90
          mqc_czt_mp2.f90
          mqc_czt_mp2_gradient.f90
          mqc_czt_mp2_hessian.f90
          mqc_czt_ri_mp2_gradient.f90
          mqc_czt_mcscf_gradient.f90
          mqc_czt_cc.f90
          mqc_czt_rcc.f90
          mqc_czt_ao_data.f90
          mqc_czt_ao.f90
          mqc_czt_xc.F90
          mqc_czt_vv10.f90
          mqc_czt_efp_potential.f90
          mqc_czt_efp_read.f90
          mqc_czt_efp_interaction.f90
          mqc_czt_efp_pair.f90
          mqc_czt_efp_rotate.f90
          mqc_czt_efp_energy.f90
          mqc_czt_sapt.f90)

# The validation and benchmark programs that exercise this backend. Registered
# from here rather than from the top-level file so that a build without CPU
# integrals never defines them: this directory only exists when the backend
# does.
include("${CMAKE_CURRENT_LIST_DIR}/validation.cmake")
