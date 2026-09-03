# Everything a method is, one level down from the dispatch that chooses one.
#
# dispatch/  what `model.method` resolves to, and which backend runs it scf/
# convergence, settings and accelerators, shared by both backends dft/       the
# exchange-correlation grid and the functional specification ci/ determinants,
# Davidson, RDMs, ORMAS -- MCSCF and CASCI both stubs/     the compiled-out twin
# of each optional backend
#
# scf/ and ci/ were `src/scf/` and `src/mcscf/`, beside `src/methods/` rather
# than inside it, which said they were peers of the whole method layer instead
# of parts of it.
add_subdirectory(dispatch)
add_subdirectory(scf)
add_subdirectory(dft)
add_subdirectory(ci)
add_subdirectory(stubs)

# Finite-difference second derivatives for any method that has a gradient and no
# analytic Hessian. It belongs to no single method, which is the whole idea, so
# it sits here rather than being filed under one of them.
target_sources(${main_lib} PRIVATE mqc_semi_numerical_hessian.f90)
