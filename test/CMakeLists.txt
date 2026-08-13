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
  addtest(mqc_libcint_cartesian)
  target_link_libraries(test_mqc_libcint_cartesian PRIVATE cint_fortran cint)

  # That a general contraction reaches libcint as one shell rather than one per
  # coefficient column, and that merging them moves no basis function.
  addtest(mqc_libcint_contraction)
  target_link_libraries(test_mqc_libcint_contraction PRIVATE cint_fortran cint)

  # The unrestricted SCF, against the properties a reference energy cannot see:
  # closed-shell UHF must equal RHF, direct must equal in-core.
  addtest(mqc_libcint_uhf)
  target_link_libraries(test_mqc_libcint_uhf PRIVATE cint_fortran cint)

  # A fragment potential written and read back. The format is GAMESS's, so
  # writer and reader are two halves of one thing and a round trip is what pins
  # them.
  addtest(mqc_efp_read)
  target_link_libraries(test_mqc_efp_read PRIVATE cint_fortran cint)

  # The free-atom initial guesses. Chiefly that every guess reaches the same
  # energy, and that SAD is spherically symmetric where SAC is not -- without
  # which the two would be one guess under two names.
  addtest(mqc_libcint_guess)
  target_link_libraries(test_mqc_libcint_guess PRIVATE cint_fortran cint)

  # Unrestricted Kohn-Sham, on the closed-shell limit: with equal spin densities
  # it must reproduce the restricted energy exactly, for every rung. That is
  # what catches libxc's interleaved arrays being read with the wrong stride and
  # the cross-spin gradient term being dropped -- both of which converge.
  addtest(mqc_libcint_uks)
  target_link_libraries(test_mqc_libcint_uks PRIVATE cint_fortran cint)

  # Range separation, on the two limits of the erf kernel. The validation
  # manifest checks the energies; this checks that omega reaches libcint at all,
  # which is the part that failed silently -- an ignored parameter returns
  # full-range exchange from the long-range pass and still converges.
  addtest(mqc_libcint_rsh)
  target_link_libraries(test_mqc_libcint_rsh PRIVATE cint_fortran cint)

  # Basis-function values on the grid, against the identity that the overlap
  # matrix is an integral of them. libcint cannot check these -- it does not
  # evaluate basis functions -- so the quadrature is the check.
  addtest(mqc_libcint_ao)
  target_link_libraries(test_mqc_libcint_ao PRIVATE cint_fortran cint)

  # Coupled cluster, against the identities a reference energy cannot see:
  # <pq||rs> antisymmetry, spin-forbidden integrals vanishing exactly, and the
  # spin-orbital MP2 equalling the conventional one. PySCF comparison lives in
  # validation/check_cc.f90, which needs the basis files.
  addtest(mqc_libcint_cc)
  target_link_libraries(test_mqc_libcint_cc PRIVATE cint_fortran cint)

  # Conventional MP2: the structural properties a reference energy cannot see,
  # chiefly that the four-index transform keeps (ia|jb) == (jb|ia).
  addtest(mqc_libcint_mp2)
  target_link_libraries(test_mqc_libcint_mp2 PRIVATE cint_fortran cint)

  # CPHF, which is the one response property needing no external reference: the
  # analytical polarizability is differenced out of the same SCF in a field, so
  # this belongs in the suite rather than in validation/.
  addtest(mqc_libcint_cphf)
  target_link_libraries(test_mqc_libcint_cphf PRIVATE cint_fortran cint)

  # Which density symmetry each direct two-electron build may be given. The fast
  # one folds permutations into a multiplicity factor and is valid only for a
  # symmetric density; asserting that it is genuinely wrong otherwise is what
  # stops the general build being deleted as redundant.
  addtest(mqc_libcint_direct)
  target_link_libraries(test_mqc_libcint_direct PRIVATE cint_fortran cint)
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
addtest(mqc_lebedev)
addtest(mqc_dft_radial)
addtest(mqc_dft_partition)
addtest(mqc_dft_grid)
addtest(mqc_dft_prune)
addtest(mqc_vibrational_analysis)
addtest(mqc_thermochemistry)
addtest(mqc_method_types)
addtest(mqc_method_placeholders)
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

if(MQC_ENABLE_TBLITE)
  addtest(mqc_method_api)
  addtest(mqc_hessian_finite_diff)
endif()
