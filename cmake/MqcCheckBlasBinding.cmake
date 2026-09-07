# Run after linking an executable on a Cray: refuse one whose BLAS would bind to
# libsci when another vendor was asked for.
#
# `cmake -DBINARY=<exe> -DVENDOR=<name> -P MqcCheckBlasBinding.cmake`
#
# The dynamic linker resolves each symbol against the executable's DT_NEEDED
# entries in order, so the first BLAS-bearing library in that list is the one
# every `dgemm` and `dgetrf` will run in. `readelf -d` prints the list; the
# check reads it rather than trusting the link line, because the `ftn` wrapper
# appends libsci to every link it performs whenever the cray-libsci module is
# loaded, and that can differ between the configure and any later rebuild.
find_program(MQC_READELF readelf)
if(NOT MQC_READELF)
  message(STATUS "readelf not found; the BLAS binding of ${BINARY} was not checked")
  return()
endif()
execute_process(
  COMMAND ${MQC_READELF} -d "${BINARY}"
  OUTPUT_VARIABLE mqc_dynamic
  RESULT_VARIABLE mqc_status
  ERROR_QUIET)
if(NOT mqc_status EQUAL 0)
  message(STATUS "readelf failed on ${BINARY}; the BLAS binding was not checked")
  return()
endif()
string(REGEX MATCHALL "NEEDED[^\n]*\\[[^]]*\\]" mqc_needed "${mqc_dynamic}")
set(mqc_first "")
foreach(mqc_entry IN LISTS mqc_needed)
  if(mqc_entry MATCHES "\\[(libsci[^]]*|libmkl[^]]*|libopenblas[^]]*|libblas[^]]*|libflexiblas[^]]*)\\]")
    set(mqc_first "${CMAKE_MATCH_1}")
    break()
  endif()
endforeach()
if(mqc_first MATCHES "^libsci")
  message(
    FATAL_ERROR
      "${BINARY} would bind its BLAS and LAPACK to Cray libsci (${mqc_first}) "
      "although ${VENDOR} was requested: libsci precedes it in the executable's "
      "dependency order. Its dgetrf is not re-entrant and slows down with every "
      "thread it is given. Reconfigure, or run `module unload cray-libsci` and "
      "link again.")
endif()
message(STATUS "BLAS binding of ${BINARY}: ${mqc_first}")
