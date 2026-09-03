target_sources(
  ${main_lib}
  PRIVATE mqc_czt_ao.f90
          mqc_czt_ao_data.f90
          mqc_czt_integrals.F90
          mqc_czt_ecp.F90
          mqc_czt_screening.f90
          mqc_czt_direct.f90
          mqc_czt_hess_ints.f90
          mqc_czt_multipole.f90
          mqc_czt_esp.f90)
