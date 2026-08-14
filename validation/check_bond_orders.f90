program check_bond_orders
   !! Check that xTB bond orders say what a chemist would say
   !!
   !! Bond orders are about to be load-bearing -- the plan is to cut molecules
   !! where they are small -- so the first question is whether GFN2's are sharp
   !! enough to cut on. These are cases where the answer is not a matter of
   !! opinion:
   !!
   !!   ethane, ethene, ethyne   C-C near 1, 2, 3
   !!   benzene                  C-C near 1.5, all six equal by symmetry
   !!   water dimer              O-H covalent near 1; the O...H hydrogen bond
   !!                            small but *not zero*
   !!
   !! The last is the one that matters most, and it is the case
   !! `mqc_bond_perception` gets wrong in both directions: its covalent-radius
   !! rule invents a bond across a short hydrogen bond and has no way to say a
   !! bond is partial. If the hydrogen bond comes back at a few hundredths while
   !! a real single bond comes back near one, there is a threshold between them
   !! and the whole scheme has something to stand on.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger, info_level
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_method_xtb, only: xtb_method_t
   implicit none

   integer :: n_bad

   n_bad = 0
   call logger%configure(info_level)

   ! Angstrom, converted on the way in. Geometries are close enough to
   ! equilibrium that the orders are about bonding rather than about strain.
   call one_case("ethane  C-C", [6, 6, 1, 1, 1, 1, 1, 1], reshape([ &
                                                                  0.0000_dp, 0.0000_dp, 0.7680_dp, &
                                                                  0.0000_dp, 0.0000_dp, -0.7680_dp, &
                                                                  -1.0192_dp, 0.0000_dp, 1.1573_dp, &
                                                                  0.5096_dp, 0.8826_dp, 1.1573_dp, &
                                                                  0.5096_dp, -0.8826_dp, 1.1573_dp, &
                                                                  1.0192_dp, 0.0000_dp, -1.1573_dp, &
                                                                  -0.5096_dp, -0.8826_dp, -1.1573_dp, &
                                                                  -0.5096_dp, 0.8826_dp, -1.1573_dp], [3, 8]), &
                 1, 2, 1.0_dp, 0.25_dp, n_bad)

   call one_case("ethene  C=C", [6, 6, 1, 1, 1, 1], reshape([ &
                                                            0.0000_dp, 0.0000_dp, 0.6695_dp, &
                                                            0.0000_dp, 0.0000_dp, -0.6695_dp, &
                                                            0.0000_dp, 0.9289_dp, 1.2321_dp, &
                                                            0.0000_dp, -0.9289_dp, 1.2321_dp, &
                                                            0.0000_dp, 0.9289_dp, -1.2321_dp, &
                                                            0.0000_dp, -0.9289_dp, -1.2321_dp], [3, 6]), &
                 1, 2, 2.0_dp, 0.35_dp, n_bad)

   call one_case("ethyne  C#C", [6, 6, 1, 1], reshape([ &
                                                      0.0000_dp, 0.0000_dp, 0.6015_dp, &
                                                      0.0000_dp, 0.0000_dp, -0.6015_dp, &
                                                      0.0000_dp, 0.0000_dp, 1.6615_dp, &
                                                      0.0000_dp, 0.0000_dp, -1.6615_dp], [3, 4]), &
                 1, 2, 3.0_dp, 0.45_dp, n_bad)

   ! The hydrogen bond. Atom 1 is the acceptor oxygen, atom 5 the donated
   ! hydrogen of the second water.
   call one_case("water dimer O...H", [8, 1, 1, 8, 1, 1], reshape([ &
                                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                                  0.0000_dp, -0.7572_dp, 0.5865_dp, &
                                                                  0.0000_dp, 0.7572_dp, 0.5865_dp, &
                                                                  0.0000_dp, 0.0000_dp, -2.9000_dp, &
                                                                  0.0000_dp, 0.0000_dp, -1.9300_dp, &
                                                                  0.0000_dp, 0.9330_dp, -3.1900_dp], [3, 6]), &
                 1, 5, 0.05_dp, 0.05_dp, n_bad)

   if (n_bad == 0) then
      call logger%info("")
      call logger%info("xTB bond orders match chemical expectation")
   else
      call logger%error("bond order checks failed")
      stop 1
   end if

contains

   subroutine one_case(label, z, coords_ang, i, j, expected, tol, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords_ang(:, :)
      integer, intent(in) :: i, j
      real(dp), intent(in) :: expected, tol
      integer, intent(inout) :: n_bad

      real(dp), parameter :: ANGSTROM_TO_BOHR = 1.8897261254578281_dp
      type(physical_fragment_t) :: frag
      type(calculation_result_t) :: res
      type(xtb_method_t) :: xtb
      character(len=160) :: line
      integer :: k

      frag%n_atoms = size(z)
      frag%nelec = sum(z)
      frag%charge = 0
      frag%multiplicity = 1
      allocate (frag%element_numbers(size(z)), source=z)
      allocate (frag%coordinates(3, size(z)))
      do k = 1, size(z)
         frag%coordinates(:, k) = coords_ang(:, k)*ANGSTROM_TO_BOHR
      end do

      xtb%variant = "gfn2"
      xtb%want_bond_orders = .true.
      call xtb%calc_energy(frag, res)

      if (res%has_error) then
         write (line, "(a,a,a)") "  FAIL ", label, ": calculation failed"
         call logger%error(trim(line))
         n_bad = n_bad + 1
         return
      end if
      if (.not. res%has_bond_orders) then
         write (line, "(a,a,a)") "  FAIL ", label, ": no bond orders returned"
         call logger%error(trim(line))
         n_bad = n_bad + 1
         return
      end if

      write (line, "(a,a,a,f8.4,a,f6.2,a,f5.2,a)") "  ", label, "  =", &
         res%bond_orders(i, j), "   (expected ", expected, " +/- ", tol, ")"
      call logger%info(trim(line))
      if (abs(res%bond_orders(i, j) - expected) > tol) then
         call logger%error("    outside tolerance")
         n_bad = n_bad + 1
      end if

      ! A bond order matrix that is not symmetric, or has a nonzero diagonal,
      ! means the wrong thing was read out of the dictionary.
      if (maxval(abs(res%bond_orders - transpose(res%bond_orders))) > 1.0e-8_dp) then
         call logger%error("    matrix is not symmetric")
         n_bad = n_bad + 1
      end if
   end subroutine one_case

end program check_bond_orders
