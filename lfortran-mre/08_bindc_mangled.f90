! LFortran 0.64.0: `--mangle-underscore-external` appends an underscore to a
! `bind(C)` symbol that has no explicit `name=`, which breaks C interop.
!
!   $ lfortran --mangle-underscore-external -c 08_bindc_mangled.f90 -o m.o
!   $ nm m.o | grep my_c_fn
!                    U my_c_fn_          <- wrong, the C symbol is my_c_fn
!
! Without the flag, or with `bind(C, name="my_c_fn")` spelled out, the symbol
! is `my_c_fn` as the standard requires: an interoperable procedure with no
! `name=` takes the lowercased Fortran name as its binding label, and no
! trailing underscore.
!
! The flag itself is needed here -- LFortran calls a plain external procedure
! by its bare name while every conventional BLAS exports `dgemm_`, so without
! it nothing links against libopenblas. The two needs collide: this is what
! makes libxc's example programs fail to link in a metalquicha build, with
! undefined references to `xc_func_end_`, `xc_func_init_` and the rest. libxc's
! Fortran interface is written `bind(c)` with no `name=` throughout.
!
! Note this MRE is checked by inspecting the symbol, not by compiling -- see
! the nm line above. run.sh reports it as "ok" because it does compile.
module bindc_mangled
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   interface
      function my_c_fn(x) bind(C)
         import :: c_int
         integer(c_int), value :: x
         integer(c_int) :: my_c_fn
      end function my_c_fn
   end interface
end module bindc_mangled

program call_it
   use bindc_mangled, only: my_c_fn
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   print *, my_c_fn(1_c_int)
end program call_it
