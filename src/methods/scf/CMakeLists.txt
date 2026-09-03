# The SCF as a shared thing: what convergence means, what the settings are, and
# the accelerators. The SCF loops themselves live in the backends, which is why
# these have to be here rather than in one of them -- cenzontle and cuEST both
# use them, and the point is that they cannot drift apart.
target_sources(
  ${main_lib}
  PRIVATE mqc_scf_types.f90 mqc_scf_common.f90 mqc_scf_convergence.f90
          mqc_fock_projector.f90 mqc_diis.f90 mqc_ediis.f90)
