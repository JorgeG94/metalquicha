# The real D4 wrapper, compiled only when dftd4 is being fetched.
#
# Outside src/ for the same reason backends/dftd3 is: the fpm build globs all of
# src/ and cannot fetch or link dftd4, so fpm compiles the stub in
# src/methods/dft/ and gets a build that refuses `"d4"` by name rather than one
# that fails to link.
target_sources(${main_lib} PRIVATE mqc_dispersion_d4.f90)
