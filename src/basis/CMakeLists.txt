target_sources(
  ${main_lib}
  PRIVATE mqc_basis_utils.f90 mqc_cgto.f90 mqc_basis_reader.f90
          mqc_gbs_reader.f90 mqc_json_basis_reader.f90
          mqc_basis_normalization.f90)
