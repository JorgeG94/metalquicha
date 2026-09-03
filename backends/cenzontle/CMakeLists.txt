# Cenzontle, the CPU ab initio backend, kept out of src/ so the fpm build never
# sees it -- fpm globs src/ and cannot fetch the integral library, so it builds
# the xTB path and CMake adds this directory when MQC_ENABLE_CZT is ON.
#
# These modules `use libcint_fortran`, which only exists once CMake has resolved
# libfint (the default) or libcint. fpm has no such dependency declared, so
# leaving them under src/ fails the fpm build outright:
#
# <ERROR> *cmd_build* Target error: Unable to find source for module dependency:
# "libcint_fortran" used by "src/mqc_czt_integrals.f90"
#
# and no CMake option helps, because fpm never reads CMake.
#
# ## The layout
#
# One directory per thing the backend does, because 49 files in a flat list is
# not a structure and cannot be given a dependency list. Roughly in the order a
# calculation moves through them:
#
# integrals/    AO integrals, ECPs, screening, the direct and multipole paths
# scf/          the SCF itself, its guesses, and continuum solvation dft/
# exchange-correlation quadrature and VV10 correlation/  MP2 and coupled
# cluster, and their derivatives mcscf/        CASCI, MCSCF, AVAS derivatives/
# gradients, Hessians and the coupled-perturbed solver analysis/     charges,
# multipoles, localisation, QUAO, IEDA, Fukui fragments/    FMO and the adjusted
# frozen orbitals that cut a bond efp/          the effective fragment potential
# sapt/         symmetry-adapted perturbation theory
#
# `mqc_czt_bridge` stays here, alone, because it is the only way in: every
# caller outside this directory reaches the backend through it and nothing
# reaches past it. Keeping it at the top says that out loud.
add_subdirectory(integrals)
add_subdirectory(scf)
add_subdirectory(dft)
add_subdirectory(correlation)
add_subdirectory(mcscf)
add_subdirectory(derivatives)
add_subdirectory(analysis)
add_subdirectory(fragments)
add_subdirectory(efp)
add_subdirectory(sapt)

target_sources(${main_lib} PRIVATE mqc_czt_bridge.f90)

# The validation and benchmark programs that exercise this backend. Registered
# from here rather than from the top-level file so that a build without CPU
# integrals never defines them: this directory only exists when the backend
# does.
include("${CMAKE_CURRENT_LIST_DIR}/validation.cmake")
