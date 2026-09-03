# The validation and benchmark programs for the CPU integrals backend.
#
# Here rather than in the top-level file because every one of them links
# ${MQC_INTEGRALS_LIBS} and means nothing without it: a build configured with
# -DMQC_ENABLE_CZT=OFF should not mention these targets at all, and the way to
# guarantee that is to define them in the directory that only exists when the
# backend does.
#
# Included from backends/cenzontle/CMakeLists.txt, which the dependency
# resolution add_subdirectory()s once it has settled which of libfint and
# libcint is providing the integrals. `${MQC_INTEGRALS_LIBS}` therefore already
# names the right one, and nothing here has to know which.

# --- integrals, and the SCF over them ---------------------------------------
mqc_add_validation_program(check_ao LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_libcint LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_rhf LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_direct LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_df LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_cc LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_sapt LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_charges LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_projection LIBRARIES ${MQC_INTEGRALS_LIBS})

# --- derivatives -------------------------------------------------------------
mqc_add_validation_program(check_ao_deriv2 LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_df_deriv LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_df_cphf LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_partition_deriv LIBRARIES
                           ${MQC_INTEGRALS_LIBS})
# The analytic gradient against finite differences of its own energy, and
# against translational invariance. Not a ctest case: it converges an SCF per
# displaced coordinate, which is seconds rather than the milliseconds the unit
# suite is meant to take.
mqc_add_validation_program(check_scf_gradient LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_mp2_gradient LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_ri_mp2_gradient LIBRARIES
                           ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_dh_gradient LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_xc_kernel LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(check_xc_potential_gradient LIBRARIES
                           ${MQC_INTEGRALS_LIBS})

# These two name libcint's targets rather than ${MQC_INTEGRALS_LIBS}. With
# libfint selected libcint is never fetched, so `cint` and `cint_fortran` are
# aliases onto `fint` -- see the alias block in cmake/MqcDependencies.cmake --
# and both spellings reach the same library. Kept as they were written rather
# than normalised, because the aliases exist precisely so that a target naming
# libcint directly does not have to be found and changed.
mqc_add_validation_program(check_df_ref_gradient LIBRARIES cint_fortran cint)
mqc_add_validation_program(check_fitted_reference_gradient LIBRARIES
                           cint_fortran cint)

# --- benchmarks --------------------------------------------------------------
mqc_add_validation_program(bench_df_k LIBRARIES ${MQC_INTEGRALS_LIBS})
# What an EFP interaction energy costs at cluster size. A dimer is too small to
# time, so this builds a lattice of fragments.
mqc_add_validation_program(bench_efp LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(bench_scf_gradient LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(bench_mp2_gradient LIBRARIES ${MQC_INTEGRALS_LIBS})
mqc_add_validation_program(bench_ri_mp2_gradient LIBRARIES
                           ${MQC_INTEGRALS_LIBS})

# --- needs libxc as well as the integrals ------------------------------------
#
# Needs both: the integrals backend for the density, libxc for the functional.
if(MQC_ENABLE_LIBXC)
  mqc_add_validation_program(
    check_dft
    LIBRARIES
    ${MQC_INTEGRALS_LIBS}
    xcf03
    xc
    INCLUDE_DIRECTORIES
    ${libxc_BINARY_DIR})
endif()

# --- FMO, the three that ctest drives ----------------------------------------
#
# IN_DEFAULT_BUILD, unlike everything above: `add_test` runs these, so excluding
# them from the default build would leave ctest with three entries and no
# executables to run.
mqc_add_validation_program(check_fmo LIBRARIES ${MQC_INTEGRALS_LIBS}
                           IN_DEFAULT_BUILD)
mqc_add_validation_program(check_fmo_mpi LIBRARIES ${MQC_INTEGRALS_LIBS}
                           IN_DEFAULT_BUILD)
mqc_add_validation_program(probe_esp LIBRARIES ${MQC_INTEGRALS_LIBS}
                           IN_DEFAULT_BUILD)
mqc_add_validation_program(sweep_fmo LIBRARIES ${MQC_INTEGRALS_LIBS})

if(MQC_ENABLE_TESTING)
  # The physics, on one rank.
  add_test(NAME ${project_name}/mqc_czt_fmo COMMAND $<TARGET_FILE:check_fmo>)

  # The embedding is local, which is what lets the cutoff exist.
  add_test(NAME ${project_name}/mqc_czt_fmo_locality
           COMMAND $<TARGET_FILE:probe_esp>)

  # Same answer however many ranks. Run on two by default, which is enough to
  # catch a missing reduction or a rank computing someone else's fragment;
  # oversubscribe so it does not depend on the machine having free cores.
  find_program(MQC_MPIEXEC NAMES mpirun mpiexec)
  if(MQC_MPIEXEC)
    add_test(NAME ${project_name}/mqc_czt_fmo_mpi
             COMMAND ${MQC_MPIEXEC} --oversubscribe -np 2
                     $<TARGET_FILE:check_fmo_mpi>)
  else()
    add_test(NAME ${project_name}/mqc_czt_fmo_mpi
             COMMAND $<TARGET_FILE:check_fmo_mpi>)
  endif()
endif()
