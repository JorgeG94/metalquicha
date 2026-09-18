# The wavefunction stability analysis.
#
# Three modules, split on their dependencies rather than on their subject.
# `mqc_czt_ov_hessian` is the electronic Hessian over the non-redundant orbital
# rotations and needs nothing but this backend. `mqc_czt_native_stability`
# diagonalises it with this program's own Davidson, so it needs nothing either.
# Both are core functionality and are always compiled -- a stability analysis
# has to work in a default build.
#
# `mqc_czt_stability` is the part that calls OpenTrustRegion to diagonalise the
# same operator, so it is compiled only when that library was fetched; when it
# was not, the stub beside it supplies a module of the same name that declines.
# It is the independent implementation the native one is checked against, not
# what a deck gets by default.
#
# That stub lives here rather than in src/methods/stubs with the others for two
# reasons: it names this backend's types in its signature, so it can only be
# compiled where they are, and fpm globs `src/` and would try to resolve those
# modules from a directory it never compiles.
target_sources(${main_lib} PRIVATE mqc_czt_ov_hessian.f90
                                   mqc_czt_native_stability.f90)

if(MQC_ENABLE_OTR)
  target_sources(${main_lib} PRIVATE mqc_czt_stability.f90)
else()
  target_sources(${main_lib} PRIVATE mqc_czt_stability_stub.f90)
endif()
