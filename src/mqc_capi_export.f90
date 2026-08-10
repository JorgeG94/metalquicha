!! Anchor for the shared library
module mqc_capi_export
   !! Exists so `libmqc.so` has a source file of its own.
   !!
   !! The C entry points are in a static archive and nothing in the program
   !! references them -- they are called from outside it -- so a linker asked
   !! for the archive in the ordinary way would take none of them. The shared
   !! library is linked with `--whole-archive` for that reason, and needs one
   !! source of its own to be a target at all. This is that source.
   implicit none
   private
end module mqc_capi_export
