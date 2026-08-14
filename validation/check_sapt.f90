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
   character(len=32) :: mode
   type(sapt_molecules_t) :: mols
   type(sapt_terms_t) :: t
   type(error_t) :: err
   real(dp) :: r, total_sum

   z = [8, 1, 1]
   sym = ["O ", "H ", "H "]

   call get_command_argument(1, mode)
   if (trim(mode) == "--scan") then
      call long_range_scan()
      stop
   end if

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

   subroutine long_range_scan()
      !! Every term's decay with separation, against what the physics demands
      !!
      !! This is the check that catches a term which is right at one geometry by
      !! coincidence. Each SAPT term has a known long-range form, and none of
      !! them can be got right by accident across a decade of separation:
      !!
      !!   Elst10   ~ R^-3   two dipoles. Positive here: the monomers' dipoles
      !!                     are parallel and perpendicular to the separation,
      !!                     which is the repulsive arrangement
      !!   Ind20    ~ R^-6   dipole inducing a dipole
      !!   Disp20   ~ R^-6   the London term, and the one whose coefficient is
      !!                     the Casimir-Polder C6 an EFP potential also carries
      !!   Exch10   exponential, not a power law -- it is an overlap effect, so
      !!                     `ln|E|` against `R` is the straight line, not
      !!                     `ln|E|` against `ln R`
      !!
      !! The exponent is measured locally between consecutive points, so it is a
      !! sequence converging on the exact value rather than a single fit that
      !! could be dragged into place by one bad point.
      integer, parameter :: N_SEP = 7          !! Separations in the scan
      integer, parameter :: N_TRACKED = 4      !! Terms whose decay is followed
      real(dp), parameter :: NOISE_FLOOR = 1.0e-16_dp
      real(dp), parameter :: SEP(N_SEP) = [4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, &
                                           8.0_dp, 10.0_dp, 12.0_dp]
      real(dp) :: a(3, 3), b(3, 3)
      real(dp) :: e(N_SEP, N_TRACKED)
      integer :: k, t_i

      a = reshape([0.0_dp, 0.0_dp, 0.10077199_dp, &
                   0.0_dp, 0.77250895_dp, -0.46780200_dp, &
                   0.0_dp, -0.77250895_dp, -0.46780200_dp], [3, 3])*ANG

      write (*, "(A)") "SAPT0 / 6-31G, water dimer, separation scan along x"
      write (*, "(A)") ""
      write (*, "(A)") "     R/Ang        Elst10        Exch10       Ind20,r"// &
         "        Disp20"
      write (*, "(A)") repeat("-", 74)
      do k = 1, size(SEP)
         b = a
         b(1, :) = b(1, :) + SEP(k)*ANG
         call build_sapt_molecules(z, sym, a, z, sym, b, "6-31g", mols, err)
         if (err%has_error()) then
            write (*, "(F10.2,A)") SEP(k), "   FAILED: "//err%get_message()
            call err%clear()
            cycle
         end if
         call run_sapt0(mols, t, err)
         if (err%has_error()) then
            write (*, "(F10.2,A)") SEP(k), "   FAILED: "//err%get_message()
            call err%clear()
            call mols%destroy()
            cycle
         end if
         e(k, 1) = t%elst10
         e(k, 2) = t%exch10
         e(k, 3) = t%ind20_r
         e(k, 4) = t%disp20
         write (*, "(F10.2,4ES14.5)") SEP(k), e(k, 1), e(k, 2), e(k, 3), e(k, 4)
         call mols%destroy()
      end do

      write (*, "(A)") ""
      write (*, "(A)") "Local decay exponent  n  in  |E| ~ R^-n , between consecutive points"
      write (*, "(A)") "     R range      Elst10        Exch10       Ind20,r"// &
         "        Disp20      expected"
      write (*, "(A)") repeat("-", 88)
      do k = 1, size(SEP) - 1
         write (*, "(F6.1,A,F5.1,4F14.3,A)") SEP(k), " -", SEP(k + 1), &
            (exponent_between(SEP(k), SEP(k + 1), e(k, t_i), e(k + 1, t_i)), &
             t_i=1, N_TRACKED), &
            "     3 / exp / 6 / 6"
      end do

      write (*, "(A)") ""
      write (*, "(A)") "Exchange is not a power law. `ln|Exch10|` against R, whose"
      write (*, "(A)") "slope is minus the decay constant and should be near constant"
      write (*, "(A)") "-- until it is not, see the note below."
      write (*, "(A)") "     R range    d ln|Exch10| / dR   (1/Ang)"
      write (*, "(A)") repeat("-", 48)
      do k = 1, size(SEP) - 1
         write (*, "(F6.1,A,F5.1,F18.3,A)") SEP(k), " -", SEP(k + 1), &
            (log(abs(e(k + 1, 2))) - log(abs(e(k, 2))))/(SEP(k + 1) - SEP(k)), &
            merge("   <- noise", "           ", abs(e(k + 1, 2)) < NOISE_FLOOR)
      end do

      write (*, "(A)") ""
      write (*, "(A)") "  Exchange falls below 1e-16 past about 7 Angstrom, at which"
      write (*, "(A)") "  point it is the numerical floor rather than physics: it is"
      write (*, "(A)") "  assembled from traces of matrices of order 1e+2, so double"
      write (*, "(A)") "  precision cannot carry it further. The decay constant is"
      write (*, "(A)") "  meaningless there and the rows say so. The physics -- that"
      write (*, "(A)") "  it decays exponentially rather than as any power of R -- is"
      write (*, "(A)") "  established well before then: seven orders of magnitude"
      write (*, "(A)") "  between 4 and 7 Angstrom, where the power-law fit returns"
      write (*, "(A)") "  a nonsensical and steadily growing exponent."
   end subroutine long_range_scan

   pure function exponent_between(r1, r2, e1, e2) result(n)
      !! `n` such that `|E| ~ R^-n` between two points
      real(dp), intent(in) :: r1, r2, e1, e2
      real(dp) :: n

      if (abs(e1) <= 0.0_dp .or. abs(e2) <= 0.0_dp) then
         n = 0.0_dp
      else
         n = -(log(abs(e2)) - log(abs(e1)))/(log(r2) - log(r1))
      end if
   end function exponent_between

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
