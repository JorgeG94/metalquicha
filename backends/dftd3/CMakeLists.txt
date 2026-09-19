# The real dispersion wrapper, compiled only when s-dftd3 is being fetched.
#
# Outside src/ for the same reason the cuEST, cenzontle and DL-FIND bridges are:
# the fpm build globs all of src/ and cannot fetch or link s-dftd3, so fpm
# compiles the stub in src/methods/dft/ and gets a build that refuses
# `keywords.dft.dispersion` rather than one that fails to link.
target_sources(${main_lib} PRIVATE mqc_dispersion.f90)
