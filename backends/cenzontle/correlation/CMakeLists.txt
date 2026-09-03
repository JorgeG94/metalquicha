target_sources(
  ${main_lib}
  PRIVATE mqc_czt_mp2.f90 mqc_czt_cc.f90 mqc_czt_rcc.f90
          mqc_czt_mp2_gradient.f90 mqc_czt_mp2_hessian.f90
          mqc_czt_ri_mp2_gradient.f90)
