# Unit testing

find_package(Threads REQUIRED)

# Macro to add a test executable for .f90 files
macro(ADDTEST name)
  add_executable(test_${name} test_${name}.f90)
  set_target_properties(
    test_${name} PROPERTIES Fortran_MODULE_DIRECTORY
                            "${CMAKE_BINARY_DIR}/tests/modules/${name}")
  target_include_directories(test_${name} PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(test_${name} PRIVATE "${main_lib}"
                                             "test-drive::test-drive")
  add_test(
    NAME ${project_name}/${name}
    COMMAND $<TARGET_FILE:test_${name}>
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADDTEST)

# Macro to add a test executable for preprocessed .F90 files
macro(ADDTESTPP name)
  add_executable(test_${name} test_${name}.F90)
  set_target_properties(
    test_${name} PROPERTIES Fortran_MODULE_DIRECTORY
                            "${CMAKE_BINARY_DIR}/tests/modules/${name}")
  target_include_directories(test_${name} PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(test_${name} PRIVATE "${main_lib}"
                                             "test-drive::test-drive")
  add_test(
    NAME ${project_name}/${name}
    COMMAND $<TARGET_FILE:test_${name}>
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADDTESTPP)

# Add individual test executables
message(STATUS "Adding individual test executables for ${project_name}...")

# Non-MPI tests
addtest(mqc_xyz_reader)
addtest(mqc_xyz_reader_detailed)
addtest(mqc_config_roundtrip)
addtest(mqc_json_config)
addtest(mqc_json_reader)
addtest(mqc_json_schema)
addtest(mqc_ecp_reader)
addtest(mqc_fraglist)
addtest(mqc_bond_perception)
addtest(mqc_fingerprint)

# The output document itself: every consumer of this program reads that file,
# and nothing checked its shape. Reached only by whole calculations before,
# which assert on energies rather than on which keys exist and what units they
# are in.
addtest(mqc_json_writer)

# What every backend is handed. `reset` restates the defaults the type already
# declares, so the two are a pair that can drift; the solvation summary is the
# only place a run says which continuum it used.
addtest(mqc_method_config)

# What each expansion context refuses to run without. A refusal is not a number,
# so no validation deck can express one: the case guarded against is a caller
# assembling a context and missing a field, not a deck missing a geometry.
addtest(mqc_mbe_contexts)
addtest(mqc_checkpoint)

# Parses basis_sets/*.json, which mqc_extract_basis_sets() unpacks from the BSE
# bundle. Only the CMake build does that, so the test compiles to a skip stub
# elsewhere (notably under fpm, where the files are absent and gitignored).
addtestpp(mqc_json_basis_reader)
target_compile_definitions(test_mqc_json_basis_reader
                           PRIVATE MQC_HAVE_BASIS_SETS)

# Which angular form the reader takes each basis to be in. Reads the same
# extracted files, so it is gated the same way.
addtestpp(mqc_basis_cartesian)
target_compile_definitions(test_mqc_basis_cartesian PRIVATE MQC_HAVE_BASIS_SETS)

# That the angular form reaches the integrals, counted in basis functions. Needs
# libcint, since libcint_molecule_t is the thing it inspects.
if(MQC_ENABLE_LIBCINT)
  # The overlap between two different bases on the same atoms, which is the one
  # new integral the quasi-atomic construction needs. Gated because
  # libcint_molecule_t is what it builds.
  addtest(mqc_mixed_overlap)
  target_link_libraries(test_mqc_mixed_overlap PRIVATE cint_fortran cint)

  # The valence-virtual space: that the projection separates cleanly, and that
  # the count does not move with the orbital basis.
  addtest(mqc_quao_vvo)
  target_link_libraries(test_mqc_quao_vvo PRIVATE cint_fortran cint)

  # The quasi-atomic orbitals themselves, and the density in their basis.
  addtest(mqc_quao)
  target_link_libraries(test_mqc_quao PRIVATE cint_fortran cint)
  # Second derivatives, the response solve and the assembled Hessian. Runs on
  # either backend now that the entry points are reached through the C ABI.
  addtest(mqc_hess_ints)
  target_link_libraries(test_mqc_hess_ints PRIVATE cint_fortran cint)
  # Fukui indices by difference in the electron count.
  addtest(mqc_fukui)
  target_link_libraries(test_mqc_fukui PRIVATE cint_fortran cint)
  # Unrestricted MP2, whose closed-shell limit is the restricted one.
  addtest(mqc_ump2)
  target_link_libraries(test_mqc_ump2 PRIVATE cint_fortran cint)
  # Regrouping the kinetic energy onto atoms and atom pairs.
  addtest(mqc_ieda)
  target_link_libraries(test_mqc_ieda PRIVATE cint_fortran cint)
  # CASCI on real molecules: the seam between the CI machinery and a basis set.
  addtest(mqc_casci)
  target_link_libraries(test_mqc_casci PRIVATE cint_fortran cint)

  addtest(mqc_mcscf)
  target_link_libraries(test_mqc_mcscf PRIVATE cint_fortran cint)

  # Choosing an active space from atomic orbital labels.
  addtest(mqc_avas)
  target_link_libraries(test_mqc_avas PRIVATE cint_fortran cint)

  addtest(mqc_libcint_cartesian)
  target_link_libraries(test_mqc_libcint_cartesian
                        PRIVATE ${MQC_INTEGRALS_LIBS})

  # That a general contraction reaches libcint as one shell rather than one per
  # coefficient column, and that merging them moves no basis function.
  addtest(mqc_libcint_contraction)
  target_link_libraries(test_mqc_libcint_contraction
                        PRIVATE ${MQC_INTEGRALS_LIBS})

  # The unrestricted SCF, against the properties a reference energy cannot see:
  # closed-shell UHF must equal RHF, direct must equal in-core.
  addtest(mqc_libcint_uhf)
  target_link_libraries(test_mqc_libcint_uhf PRIVATE ${MQC_INTEGRALS_LIBS})

  # A fragment potential written and read back. The format is GAMESS's, so
  # writer and reader are two halves of one thing and a round trip is what pins
  # them.
  addtest(mqc_efp_read)
  target_link_libraries(test_mqc_efp_read PRIVATE ${MQC_INTEGRALS_LIBS})

  # The electrostatic interaction between two fragments, against the energy
  # GAMESS reports for the same dimer. Chiefly a guard on the charge-dipole
  # sign, which both ways round gives a plausible-looking answer.
  addtest(mqc_efp_interaction)
  target_link_libraries(test_mqc_efp_interaction PRIVATE ${MQC_INTEGRALS_LIBS})

  # A molecule spanning two fragments, checked by the block structure of its own
  # overlap: the diagonal blocks must be each fragment's overlap exactly.
  addtest(mqc_efp_pair)
  target_link_libraries(test_mqc_efp_pair PRIVATE ${MQC_INTEGRALS_LIBS})

  # All five EFP2 terms assembled into one energy. Each was pinned against
  # GAMESS on its own, so what this checks is the bookkeeping: no term dropped,
  # none counted twice, dispersion totalled damped rather than undamped.
  addtest(mqc_efp_energy)
  target_link_libraries(test_mqc_efp_energy PRIVATE ${MQC_INTEGRALS_LIBS})

  # Placing a fragment in an arbitrary orientation, piece by piece: the frame
  # construction, the Cartesian shell rotation as a group representation, and
  # the invariants a rotated fragment must keep. The energy test asserts the
  # same thing end to end but cannot say which piece broke.
  addtest(mqc_efp_rotate)
  target_link_libraries(test_mqc_efp_rotate PRIVATE ${MQC_INTEGRALS_LIBS})

  # The free-atom initial guesses. Chiefly that every guess reaches the same
  # energy, and that SAD is spherically symmetric where SAC is not -- without
  # which the two would be one guess under two names.
  addtest(mqc_libcint_guess)
  target_link_libraries(test_mqc_libcint_guess PRIVATE ${MQC_INTEGRALS_LIBS})

  # Unrestricted Kohn-Sham, on the closed-shell limit: with equal spin densities
  # it must reproduce the restricted energy exactly, for every rung. That is
  # what catches libxc's interleaved arrays being read with the wrong stride and
  # the cross-spin gradient term being dropped -- both of which converge.
  addtest(mqc_libcint_uks)
  target_link_libraries(test_mqc_libcint_uks PRIVATE ${MQC_INTEGRALS_LIBS})

  # Range separation, on the two limits of the erf kernel. The validation
  # manifest checks the energies; this checks that omega reaches libcint at all,
  # which is the part that failed silently -- an ignored parameter returns
  # full-range exchange from the long-range pass and still converges.
  addtest(mqc_libcint_rsh)
  target_link_libraries(test_mqc_libcint_rsh PRIVATE ${MQC_INTEGRALS_LIBS})

  # Basis-function values on the grid, against the identity that the overlap
  # matrix is an integral of them. libcint cannot check these -- it does not
  # evaluate basis functions -- so the quadrature is the check.
  addtest(mqc_libcint_ao)
  target_link_libraries(test_mqc_libcint_ao PRIVATE ${MQC_INTEGRALS_LIBS})

  # The potential integral against its own centre-gradient. This is what fixes
  # the sign of libcint's `nabla-rinv` for every caller, charge transfer's
  # dipole rank being the first.
  addtest(mqc_libcint_esp)
  target_link_libraries(test_mqc_libcint_esp PRIVATE ${MQC_INTEGRALS_LIBS})
  # The continuum: cavity geometry, the Born ion in closed form, and one
  # solvated SCF pinned against PySCF.
  addtest(mqc_pcm)
  target_link_libraries(test_mqc_pcm PRIVATE ${MQC_INTEGRALS_LIBS})
  # SAPT0's molecules: ghosts in the dimer basis. Everything checkable before a
  # term is computed -- that ghosts keep their basis and lose their charge, and
  # that a monomer's AO ordering is the dimer's, which Tr(D S) catches.
  addtest(mqc_sapt)
  target_link_libraries(test_mqc_sapt PRIVATE ${MQC_INTEGRALS_LIBS})

  # Coupled cluster, against the identities a reference energy cannot see:
  # <pq||rs> antisymmetry, spin-forbidden integrals vanishing exactly, and the
  # spin-orbital MP2 equalling the conventional one. PySCF comparison lives in
  # validation/check_cc.f90, which needs the basis files.
  addtest(mqc_libcint_cc)
  target_link_libraries(test_mqc_libcint_cc PRIVATE ${MQC_INTEGRALS_LIBS})

  # Spin-adapted coupled cluster, tied to the spin-orbital path above. For a
  # closed shell the two are the same theory written twice and must agree to
  # machine precision, which catches the transposed index that a reference
  # energy alone would not -- see the module header for why that is the failure
  # mode worth designing against.
  addtest(mqc_libcint_rcc)
  target_link_libraries(test_mqc_libcint_rcc PRIVATE cint_fortran cint)

  # Conventional MP2: the structural properties a reference energy cannot see,
  # chiefly that the four-index transform keeps (ia|jb) == (jb|ia).
  addtest(mqc_libcint_mp2)
  target_link_libraries(test_mqc_libcint_mp2 PRIVATE ${MQC_INTEGRALS_LIBS})

  # CPHF, which is the one response property needing no external reference: the
  # analytical polarizability is differenced out of the same SCF in a field, so
  # this belongs in the suite rather than in validation/.
  addtest(mqc_libcint_cphf)
  target_link_libraries(test_mqc_libcint_cphf PRIVATE ${MQC_INTEGRALS_LIBS})

  # Ghost centres and the superposition error they expose. Everything a
  # counterpoise-corrected expansion rests on, checked before any expansion is
  # assembled: that a ghost drops its nucleus and keeps its functions, that the
  # AO space is untouched, and that a monomer really is lowered by its partner's
  # basis -- which is the quantity the correction removes. Which partitions FMO
  # refuses. A refusal is not a number, so no validation deck can express it --
  # and the refusals happen before any SCF, so this is milliseconds rather than
  # the seconds check_fmo spends on the physics.
  addtest(mqc_fmo_partitions)
  target_link_libraries(test_mqc_fmo_partitions PRIVATE ${MQC_INTEGRALS_LIBS})

  addtest(mqc_counterpoise)
  target_link_libraries(test_mqc_counterpoise PRIVATE ${MQC_INTEGRALS_LIBS})

  # Which density symmetry each direct two-electron build may be given. The fast
  # one folds permutations into a multiplicity factor and is valid only for a
  # symmetric density; asserting that it is genuinely wrong otherwise is what
  # stops the general build being deleted as redundant.
  addtest(mqc_libcint_direct)
  target_link_libraries(test_mqc_libcint_direct PRIVATE ${MQC_INTEGRALS_LIBS})

  # Atomic charges, against the arithmetic every population analysis owes: both
  # schemes sum to the molecular charge, every basis function belongs to exactly
  # one atom, and the CHELPG grid stays outside the molecule it fits. The
  # numbers themselves are compared against other programs in validation/; none
  # of that catches a constraint row assembled wrongly, which still fits.
  addtest(mqc_libcint_charges)
  target_link_libraries(test_mqc_libcint_charges PRIVATE ${MQC_INTEGRALS_LIBS})

  # The multipole matrices every distributed multipole and every localization
  # centroid is built from. Origin shifts, symmetry, and the redundancy in the
  # full Cartesian tensor -- the three things a transposed component breaks
  # while leaving the reported dipole intact.
  addtest(mqc_libcint_multipole)
  target_link_libraries(test_mqc_libcint_multipole
                        PRIVATE ${MQC_INTEGRALS_LIBS})

  # Foster-Boys localization, against what a rotation among occupied orbitals
  # cannot change: the density, orthonormality, and the direction its own
  # functional moves.
  addtest(mqc_libcint_localize)
  target_link_libraries(test_mqc_libcint_localize PRIVATE ${MQC_INTEGRALS_LIBS})

  # Distributed multipoles, against the molecular moments they are a
  # redistribution of: the site charges must account for every electron, and the
  # site charges and dipoles together must reproduce the molecular dipole.
  # Everything an effective fragment potential's electrostatics rests on.
  addtest(mqc_libcint_dma)
  target_link_libraries(test_mqc_libcint_dma PRIVATE ${MQC_INTEGRALS_LIBS})

  # Carrying a density between basis sets, which is what the projection guess
  # does one rung at a time. Checked against the properties that make the result
  # a density rather than a matrix of the right shape -- idempotency, the
  # electron count, and the identity when both bases are the same.
  addtest(mqc_libcint_projection)
  target_link_libraries(test_mqc_libcint_projection
                        PRIVATE ${MQC_INTEGRALS_LIBS})
endif()

# Only meaningful with the cuEST backend compiled in.
if(MQC_ENABLE_CUEST)
  addtestpp(mqc_cuest_spin)
  target_compile_definitions(test_mqc_cuest_spin PRIVATE MQC_WITH_CUEST)
endif()
# The continuum cavity radii. Unconditional: the table is ours, a wrong entry
# gives a wrong solvation energy with nothing to object, and neither cuEST nor
# libxc is needed to check it.
addtest(mqc_pcm_radii)
addtest(mqc_atomic_radii)
addtest(mqc_scf_common)
addtest(mqc_shared_helpers)

# What a functional name means as a composition. Unconditional: the
# double-hybrid coefficients are ours, and libxc is not needed to check them.
addtest(mqc_xc_spec)

addtest(mqc_basis_utils)
addtest(mqc_elements)
addtest(mqc_cgto)
addtest(mqc_physical_fragment)
addtest(mqc_frag_utils)
addtest(mqc_gmbe_intersection)
addtest(mqc_fragment_from_atoms)
addtest(mqc_gmbe_energy)
addtest(mqc_mbe)
addtest(mqc_geometry)
addtest(mqc_config_adapter)
addtest(mqc_cli_parser)
addtest(mqc_split_lines)
addtest(mqc_result_types)
addtest(mqc_finite_differences)
addtest(mqc_diis)
addtest(mqc_fock_projector)
addtest(mqc_afo_model)
addtest(mqc_aambs)
addtest(mqc_determinants)
addtest(mqc_ci)
addtest(mqc_ci_transform)
addtest(mqc_davidson)
addtest(mqc_rdm)
addtest(mqc_ormas_space)
addtest(mqc_ormas_ci)
addtest(mqc_lebedev)
addtest(mqc_dft_radial)
addtest(mqc_dft_partition)
addtest(mqc_dft_grid)
addtest(mqc_dft_prune)
addtest(mqc_vibrational_analysis)
addtest(mqc_thermochemistry)
addtest(mqc_method_types)
addtest(mqc_method_placeholders)
# Which backend a run credits. Resolved from the method, the request and what
# linked, so this runs on any build -- with and without cuEST it just expects a
# different answer for `auto`.
addtest(mqc_acknowledgements)
# The optimizer's vocabulary, not the DL-FIND bridge: mqc_optimizer_types is
# compiled whether or not MQC_ENABLE_DLFIND is on, so this runs in CI, where the
# optimizer backend is off.
addtest(mqc_optimizer_types)
# The term list a frozen optimization is built on. Fragment layer, so this runs
# regardless of MQC_ENABLE_DLFIND.
addtest(mqc_term_list)
addtest(mqc_work_queue)

# C-API / external-interface layer. These call the bind(C) entry points the way
# a C or Python caller would; none of them needs MPI or a QC backend.
addtest(mqc_validate)
addtest(mqc_capi_system)
addtest(mqc_capi_fraglist)
addtest(mqc_capi_session)
addtest(mqc_capi_run)

# The checkpoint format only exists with HDF5 on, so the test that exercises its
# ragged packing and sort ordering follows the same guard the format does.
if(MQC_ENABLE_HDF5)
  addtest(mqc_hdf5_checkpoint)
endif()

if(MQC_ENABLE_TBLITE)
  addtest(mqc_method_api)
  addtest(mqc_hessian_finite_diff)
endif()
