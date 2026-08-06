# Unpack the Basis Set Exchange bundle into basis_sets/ at configure time.
#
# The bundle ships every basis BSE knows about -- 922 of them, 572 MB unpacked
# -- which is far more than any build needs and far too much to keep in git.
# Only the names in MQC_BASIS_SETS are extracted.
#
# Bundle members are "<bse-name>.<revision>.json", with several revisions of the
# same basis side by side. Only the revision suffix is stripped: the file keeps
# its Basis Set Exchange name, which is exactly what
# mqc_basis_utils::normalize_basis_name produces, so there is no mapping between
# the two languages to keep in step.
#
# The highest revision of each basis wins. Extraction is skipped when every
# wanted file is already present and newer than the tarball, so a reconfigure
# costs nothing.

# The on-disk stem for a bundle name. Identity apart from case, but kept as a
# function so the one place that defines it stays obvious.
function(mqc_canonical_basis_name bse_name out_var)
  string(TOLOWER "${bse_name}" name)
  set(${out_var}
      "${name}"
      PARENT_SCOPE)
endfunction()

function(mqc_extract_basis_sets)
  set(bundle "${CMAKE_CURRENT_SOURCE_DIR}/basis_sets/${MQC_BASIS_BUNDLE}")
  set(dest "${CMAKE_CURRENT_SOURCE_DIR}/basis_sets")

  if(NOT EXISTS "${bundle}")
    message(
      FATAL_ERROR
        "Basis set bundle not found: ${bundle}\n"
        "Download basis_sets-json-<ver>.tar.bz2 from https://www.basissetexchange.org/downloads "
        "and put it in basis_sets/, or point -DMQC_BASIS_BUNDLE= at the one you have."
    )
  endif()

  # ---- do two entries want the same file? -----------------------------------
  #
  # Only reachable by listing a name twice in different cases now that the stem
  # is the bundle name, but a duplicate would still have one extraction silently
  # overwrite the other, so say so rather than pick.
  set(seen_stems "")
  foreach(bse_name ${MQC_BASIS_SETS})
    mqc_canonical_basis_name("${bse_name}" stem)
    if("${stem}" IN_LIST seen_stems)
      message(
        FATAL_ERROR
          "MQC_BASIS_SETS lists '${bse_name}' and '${owner_of_${stem}}', which "
          "both resolve to ${stem}.json. Keep one.")
    endif()
    list(APPEND seen_stems "${stem}")
    set(owner_of_${stem} "${bse_name}")
  endforeach()

  # ---- drop anything the bundle is no longer being asked for ----------------
  #
  # basis_sets/*.json is entirely derived -- only the bundle and PROVENANCE.md
  # are tracked -- so the directory should hold exactly what MQC_BASIS_SETS asks
  # for and nothing else. Without this, a name dropped from the list, or left
  # over from an older naming scheme, stays on disk and keeps resolving: a basis
  # that quietly outlives the decision to stop shipping it.
  file(GLOB present "${dest}/*.json")
  foreach(path ${present})
    # NAME, not NAME_WE: the latter strips from the first dot, which would maim
    # any basis name that ever contains one.
    get_filename_component(stem "${path}" NAME)
    string(REGEX REPLACE "\\.json$" "" stem "${stem}")
    string(TOLOWER "${stem}" stem)
    if(NOT "${stem}" IN_LIST seen_stems)
      message(STATUS "Basis sets: removing stale ${stem}.json")
      file(REMOVE "${path}")
    endif()
  endforeach()

  # ---- is there anything to do? ---------------------------------------------
  set(missing "")
  foreach(bse_name ${MQC_BASIS_SETS})
    mqc_canonical_basis_name("${bse_name}" stem)
    if(NOT EXISTS "${dest}/${stem}.json" OR "${bundle}" IS_NEWER_THAN
                                            "${dest}/${stem}.json")
      list(APPEND missing "${bse_name}")
    endif()
  endforeach()

  if(NOT missing)
    list(LENGTH MQC_BASIS_SETS n_total)
    message(STATUS "Basis sets: ${n_total} present in ${dest}")
    return()
  endif()

  list(LENGTH missing n_missing)
  message(STATUS "Basis sets: extracting ${n_missing} from ${MQC_BASIS_BUNDLE}")

  # ---- list the bundle once -------------------------------------------------
  #
  # One tar pass to learn which revisions exist, rather than one probe per
  # basis: the archive is bzip2 and a full scan takes seconds, which is not a
  # thing to do 900 times.
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar tjf "${bundle}"
    OUTPUT_VARIABLE listing
    ERROR_VARIABLE tar_error
    RESULT_VARIABLE tar_status
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(NOT tar_status EQUAL 0)
    message(FATAL_ERROR "Could not read ${bundle}: ${tar_error}")
  endif()

  string(REPLACE "\n" ";" members "${listing}")

  set(wanted_members "")
  set(not_in_bundle "")
  foreach(bse_name ${missing})
    string(TOLOWER "${bse_name}" lower)
    # Regex-escape the name: "6-31+g_st_" is full of metacharacters, and an
    # unescaped "+" would match a different basis or nothing at all.
    string(REGEX REPLACE "([][+*()^$.?|\\\\])" "\\\\\\1" escaped "${lower}")

    set(best "")
    set(best_revision -1)
    foreach(member ${members})
      if(member MATCHES "/${escaped}\\.([0-9]+)\\.json$")
        if(CMAKE_MATCH_1 GREATER best_revision)
          set(best_revision "${CMAKE_MATCH_1}")
          set(best "${member}")
        endif()
      endif()
    endforeach()

    if(best)
      list(APPEND wanted_members "${best}")
      set(member_for_${lower} "${best}")
    else()
      list(APPEND not_in_bundle "${bse_name}")
    endif()
  endforeach()

  if(not_in_bundle)
    string(REPLACE ";" ", " pretty "${not_in_bundle}")
    message(
      FATAL_ERROR
        "Not in ${MQC_BASIS_BUNDLE}: ${pretty}\n"
        "MQC_BASIS_SETS takes Basis Set Exchange names -- lowercase, with '*' "
        "written '_st_' (so 6-31G** is 6-31g_st__st_). Check the spelling at "
        "https://www.basissetexchange.org.")
  endif()

  # ---- extract to a scratch dir, then rename into place ---------------------
  set(staging "${CMAKE_BINARY_DIR}/basis_sets_staging")
  file(REMOVE_RECURSE "${staging}")
  file(MAKE_DIRECTORY "${staging}")

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${bundle}" -- ${wanted_members}
    WORKING_DIRECTORY "${staging}"
    ERROR_VARIABLE tar_error
    RESULT_VARIABLE tar_status)

  if(NOT tar_status EQUAL 0)
    message(FATAL_ERROR "Could not extract from ${bundle}: ${tar_error}")
  endif()

  foreach(bse_name ${missing})
    string(TOLOWER "${bse_name}" lower)
    mqc_canonical_basis_name("${bse_name}" stem)
    set(member "${member_for_${lower}}")
    if(NOT EXISTS "${staging}/${member}")
      message(FATAL_ERROR "Expected ${member} in the extracted bundle")
    endif()
    file(RENAME "${staging}/${member}" "${dest}/${stem}.json")
    # tar restores each member's archived mtime, which predates the bundle file
    # itself -- so without this the freshness check above fires on every
    # configure and re-extracts all 40. Stamp them as of now instead.
    file(TOUCH "${dest}/${stem}.json")
  endforeach()

  file(REMOVE_RECURSE "${staging}")
  message(STATUS "Basis sets: ${dest}")
endfunction()
