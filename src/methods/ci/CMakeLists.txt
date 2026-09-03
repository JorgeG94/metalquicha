# Configuration interaction: determinants, the Davidson solver, the reduced
# density matrices and the ORMAS active-space machinery. Used by MCSCF and by
# CASCI, so it belongs to neither.
target_sources(
  ${main_lib}
  PRIVATE mqc_determinants.f90
          mqc_ci.f90
          mqc_ci_transform.f90
          mqc_davidson.f90
          mqc_rdm.f90
          mqc_ormas_space.f90
          mqc_ormas_ci.f90)
