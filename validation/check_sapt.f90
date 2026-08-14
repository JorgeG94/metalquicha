!! SAPT0 on a pair of waters from the prism cluster
!!
!!     cmake --build build --target check_sapt
!!     ./build/check_sapt
!!
!! SAPT0 decomposes the interaction of *two* monomers, so a six-water cluster is
!! fifteen SAPT calculations rather than one. This runs the pairs of
!! `validation/inputs/sample_inputs/prism.xyz` and prints the breakdown, which is
!! the thing SAPT is for: not the total, which a supermolecular calculation also
!! gives, but which physics it is made of.
!!
!! The prism is a harder test than the reference dimer in `test_mqc_sapt`. Those
!! two waters share an orientation and sit 3 Angstrom apart on an axis; these are
!! hydrogen-bonded at assorted angles and distances, so every term is exercised
!! at a geometry nobody chose to be convenient.
program check_sapt
   use pic_types, only: dp
   use mqc_sapt, only: sapt_molecules_t, build_sapt_molecules, sapt_terms_t, run_sapt0
   use mqc_error, only: error_t
   implicit none

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer, parameter :: N_MONOMER = 6
   integer, parameter :: N_ATOM = 3                  !! Atoms in a water
   integer, parameter :: N_PRISM_ATOMS = N_MONOMER*N_ATOM

   real(dp) :: prism(3, 3, N_MONOMER)
   integer :: z(3)
   character(len=2) :: sym(3)
   integer :: i, j, pair_count
   type(sapt_molecules_t) :: mols
   type(sapt_terms_t) :: t
   type(error_t) :: err
   real(dp) :: r, total_sum

   z = [8, 1, 1]
   sym = ["O ", "H ", "H "]
   call load_prism(prism)

   write (*, "(A)") "SAPT0 / 6-31G on the water prism, pair by pair"
   write (*, "(A)") "geometry: validation/inputs/sample_inputs/prism.xyz"
   write (*, "(A)") ""
   write (*, "(A)") "  pair    R(O-O)      Elst10      Exch10     Ind20,r  "// &
      "ExchInd20      Disp20  ExchDisp20       dHF       TOTAL"
   write (*, "(A)") repeat("-", 118)

   pair_count = 0
   total_sum = 0.0_dp
   do i = 1, N_MONOMER
      do j = i + 1, N_MONOMER
         call build_sapt_molecules(z, sym, prism(:, :, i), z, sym, prism(:, :, j), &
                                   "6-31g", mols, err)
         if (err%has_error()) then
            write (*, "(A,I0,A,I0,A,A)") "  ", i, "-", j, "  FAILED: ", err%get_message()
            call err%clear()
            cycle
         end if
         call run_sapt0(mols, t, err)
         if (err%has_error()) then
            write (*, "(A,I0,A,I0,A,A)") "  ", i, "-", j, "  FAILED: ", err%get_message()
            call err%clear()
            call mols%destroy()
            cycle
         end if

         r = sqrt(sum((prism(:, 1, i) - prism(:, 1, j))**2))/ANG
         write (*, "(2X,I1,A,I1,4X,F7.3,7F12.8)") i, "-", j, r, &
            t%elst10, t%exch10, t%ind20_r, t%exch_ind20_r, &
            t%disp20, t%exch_disp20, t%delta_hf
         write (*, "(A,F12.8)") repeat(" ", 106), t%total
         pair_count = pair_count + 1
         total_sum = total_sum + t%total
         call mols%destroy()
      end do
   end do

   write (*, "(A)") repeat("-", 118)
   write (*, "(A,I0,A)") "  ", pair_count, " pairs"
   write (*, "(A,F14.8,A)") "  sum of pair totals  ", total_sum, " Hartree"
   write (*, "(A,F14.8,A)") "                      ", total_sum*627.509474_dp, " kcal/mol"
   write (*, "(A)") ""
   write (*, "(A)") "  The sum is the two-body part of the cluster's binding energy."
   write (*, "(A)") "  It is not the binding energy: three-body induction in a"
   write (*, "(A)") "  hydrogen-bonded cluster is worth several kcal/mol, and no"
   write (*, "(A)") "  sum over pairs contains it."

contains

   subroutine load_prism(geom)
      !! The six monomers, in Bohr. Hard-coded rather than read, so this runs
      !! from any directory; the source is the xyz named in the header.
      real(dp), intent(out) :: geom(3, 3, N_MONOMER)

      real(dp) :: raw(3, N_PRISM_ATOMS)
      integer :: m, a, k

      raw = reshape([ &
                    0.8974_dp, -1.285111_dp, 1.375674_dp, &
                    0.93366_dp, -1.620249_dp, 0.461291_dp, &
                    1.538424_dp, -0.56327_dp, 1.325981_dp, &
                    -1.559218_dp, -0.241778_dp, 1.423474_dp, &
                    -2.064102_dp, -0.501966_dp, 2.197464_dp, &
                    -0.677656_dp, -0.678646_dp, 1.525237_dp, &
                    2.090789_dp, 0.984452_dp, -0.045693_dp, &
                    1.258568_dp, 1.498421_dp, -0.018646_dp, &
                    2.788216_dp, 1.642098_dp, -0.093886_dp, &
                    -0.447032_dp, 2.02126_dp, 0.012114_dp, &
                    -0.903998_dp, 1.539985_dp, 0.716909_dp, &
                    -0.911081_dp, 1.700472_dp, -0.771056_dp, &
                    -1.776681_dp, -0.211411_dp, -1.351176_dp, &
                    -2.560698_dp, -0.516179_dp, -1.814061_dp, &
                    -1.949801_dp, -0.387242_dp, -0.40953_dp, &
                    0.9244_dp, -1.339405_dp, -1.432395_dp, &
                    0.044626_dp, -0.992399_dp, -1.635267_dp, &
                    1.466581_dp, -0.545087_dp, -1.34041_dp], [3, N_PRISM_ATOMS])

      k = 0
      do m = 1, N_MONOMER
         do a = 1, 3
            k = k + 1
            geom(:, a, m) = raw(:, k)*ANG
         end do
      end do
   end subroutine load_prism

end program check_sapt
