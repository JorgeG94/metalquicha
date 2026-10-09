target_sources(
  ${main_lib}
  PRIVATE mqc_czt_density_contraction.f90
          mqc_czt_effective_density.f90
          mqc_czt_gradient.f90
          mqc_czt_hessian.f90
          mqc_czt_hess_abi.f90
          mqc_czt_cphf.F90
          mqc_czt_response.f90
          mqc_czt_response_product.f90
          mqc_czt_rpa_solver.f90
          mqc_czt_tddft.f90
          mqc_czt_tddft_properties.f90)
