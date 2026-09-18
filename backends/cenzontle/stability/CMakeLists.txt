# The wavefunction stability analysis.
#
# Two modules, split on their dependencies rather than on their subject.
# `mqc_czt_ov_hessian` is the electronic Hessian over the non-redundant orbital
# rotations and needs nothing but this backend, so it is always compiled and
# always available to a test. `mqc_czt_stability` is the part that calls
# OpenTrustRegion to diagonalise it, so it is compiled only when that library
# was fetched -- and when it was not, the stub beside it supplies a module of
# the same name that declines.
#
# The stub lives here rather than in src/methods/stubs with the others because
# it names this backend's types in its signature and so can only be compiled
# where they are, and because fpm globs `src/` and would try to resolve those
# modules from a directory it never compiles.
target_sources(${main_lib} PRIVATE mqc_czt_ov_hessian.f90)

if(MQC_ENABLE_OTR)
  target_sources(${main_lib} PRIVATE mqc_czt_stability.f90)
else()
  target_sources(${main_lib} PRIVATE mqc_czt_stability_stub.f90)
endif()
