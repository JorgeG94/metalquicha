!! Writes an amplitude file for a reader that is not Fortran to check
!!
!!     cmake -B build -DMQC_ENABLE_HDF5=ON
!!     cmake --build build --target check_amplitude_layout
!!     ./build/check_amplitude_layout && python3 validation/check_amplitude_layout.py
!!
!! This program asserts nothing about the file it writes, and that is the
!! point. The failure it exists to catch is invisible from inside Fortran:
!! HDF5 declares extents in C order while Fortran's leftmost index varies
!! fastest, so a block written with `shape()` passed straight through holds
!! the right bytes under the wrong declared shape. Reading it back through
!! the same bindings is exact, because the same mistake is made twice. Only a
!! reader with the opposite convention can see it, so the assertions live in
!! the Python script and all this does is produce the evidence.
!!
!! **The extents are deliberately all different.** `t2(2,3,4,5)` rather than
!! anything square: a reversal bug in a block whose dimensions match is
!! undetectable, and `nocc`/`nvir` being equal is exactly the accident that
!! would let one through on a small test case.
program check_amplitude_layout
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_hdf5_amplitudes, only: hdf5_amplitudes_t
   implicit none

   character(len=*), parameter :: PATH = "amplitude_layout.h5"
   character(len=*), parameter :: FP = "1234abcd5678ef90"
   integer, parameter :: N1 = 2, N2 = 3, N3 = 4, N4 = 5
   integer, parameter :: ITERATION = 3
   real(dp), parameter :: ENERGY = -0.5_dp

   type(hdf5_amplitudes_t) :: store
   type(error_t) :: err
   real(dp), allocatable :: t1(:, :), t2(:, :, :, :)
   integer :: i, a, j, b, unit
   logical :: exists

   ! Every element encodes its own indices, so a block that arrives
   ! transposed is caught by its values and not only by its shape.
   allocate (t1(N1, N2))
   do a = 1, N2
      do i = 1, N1
         t1(i, a) = real(100*i + a, dp)
      end do
   end do

   allocate (t2(N1, N2, N3, N4))
   do b = 1, N4
      do j = 1, N3
         do a = 1, N2
            do i = 1, N1
               t2(i, a, j, b) = real(1000*i + 100*a + 10*j + b, dp)
            end do
         end do
      end do
   end do

   inquire (file=PATH, exist=exists)
   if (exists) then
      open (newunit=unit, file=PATH, status="old", action="readwrite")
      close (unit, status="delete")
   end if

   call store%open(PATH, FP, err)
   if (err%has_error()) then
      write (*, "(a)") "could not open "//PATH//": "//err%get_message()
      error stop 1
   end if

   call store%put("t1", t1)
   call store%put("t2", t2)
   call store%commit(ITERATION, ENERGY)
   call store%close()

   write (*, "(a)") "wrote "//PATH
   write (*, "(a,4(1x,i0))") "  t1 Fortran extents:", N1, N2
   write (*, "(a,4(1x,i0))") "  t2 Fortran extents:", N1, N2, N3, N4
   write (*, "(a)") "  a C-order reader must see these reversed"
end program check_amplitude_layout
