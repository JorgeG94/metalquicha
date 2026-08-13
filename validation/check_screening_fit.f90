!! Fit the charge-penetration screening exponents, in Fortran
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON
!!     ./build/check_screening_fit
!!
!! The two screening blocks of a potential, produced rather than validated: the grid,
!! the objective and the search all live in `mqc_libcint_screening`, so this is the
!! path a real potential would take. `tools/efp_validation/check_screening.py` is what
!! established that the objective is GAMESS's -- by showing every one of its published
!! exponents is stationary under it -- and this reproduces that objective's minimum
!! from the Fortran side.
!!
!! Expect the exponents *not* to equal a reference's. At the upper bound the damping
!! is off and the objective is flat, so an exponent that wanders there sticks; water's
!! reference has one bond midpoint at 10.0 and the other at 2.48 with nothing physical
!! separating them. What is reproducible is the objective, and this reaches 0.641
!! kcal/mol for the exponential form where GAMESS's own exponents give 0.649.
program check_screening_fit
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_dma, only: dma_result_t, distributed_multipoles
   use mqc_libcint_screening, only: fit_screening, SCREEN_EXPONENTIAL, SCREEN_GAUSSIAN
   use mqc_error, only: error_t
   implicit none
   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp
   type(libcint_molecule_t) :: mol
   type(rhf_result_t) :: scf
   type(dma_result_t) :: dma
   type(error_t) :: err
   real(dp), allocatable :: alpha(:)
   real(dp) :: c(3, 3)
   real(dp) :: resid
   integer :: z(3)
   integer :: k, npts, kind
   character(len=2) :: sym(3)
   z = [8, 1, 1]
   sym = ["O ", "H ", "H "]
   c = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
   call build_libcint_molecule(z, sym, c, "6-31g*", mol, err)
   call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
   call distributed_multipoles(mol, scf%density, z, dma, err)
   if (err%has_error()) then
      write (*, *) "dma failed: ", err%get_message()
      stop 1
   end if
   do kind = 1, 2
      call fit_screening(mol, scf%density, dma, z, kind, alpha, err, &
                         residual=resid, grid_size=npts)
      if (err%has_error()) then
         write (*, *) "fit failed: ", err%get_message()
         stop 1
      end if
      write (*, "(a,i0,a,i0,a,f10.5,a)") " kind ", kind, "  grid ", npts, &
         "  residual ", resid, " kcal/mol"
      do k = 1, size(alpha)
         write (*, "(a,a8,f14.6)") "    ", dma%labels(k), alpha(k)
      end do
      deallocate (alpha)
   end do
   call mol%destroy()
end program check_screening_fit
