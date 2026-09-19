# The wavefunction stability analysis.
#
# Two modules, split on their dependencies rather than on their subject.
# `mqc_czt_ov_hessian` is the electronic Hessian over the non-redundant orbital
# rotations and needs nothing but this backend. `mqc_czt_native_stability`
# diagonalises it with this program's own Davidson, so it needs nothing either.
# Both are core functionality and are always compiled -- a stability analysis
# has to work in a default build, and the second-order SCF applies the same
# operator.
target_sources(${main_lib} PRIVATE mqc_czt_ov_hessian.f90
                                   mqc_czt_native_stability.f90)
