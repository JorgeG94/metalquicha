target_sources(
  ${main_lib}
  PRIVATE mqc_czt_mcscf.f90
          mqc_czt_mcscf_gradient.f90
          mqc_czt_casci.f90
          mqc_czt_avas.f90
          mqc_czt_sa_hessian.f90
          mqc_czt_sa_effective_density.f90
          mqc_czt_sa_gradient.f90
          mqc_czt_sa_nac.f90)
