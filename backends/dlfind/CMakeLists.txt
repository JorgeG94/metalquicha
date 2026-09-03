# The real DL-FIND bridge, compiled only when libdlfind is being fetched.
#
# It lives outside src/ for the same reason the cuEST and cenzontle backends do:
# the fpm build globs all of src/, and fpm cannot fetch and link libdlfind, so
# fpm gets the stub and a build that refuses `driver: "Optimize"` rather than
# one that fails to link.
target_sources(${main_lib} PRIVATE mqc_dlfind_bridge.f90)
