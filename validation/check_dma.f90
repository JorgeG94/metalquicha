!! Manual check of the distributed multipoles
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON
!!     ./build/check_dma
!!     python3 validation/check_dma.py
!!
!! Milestone 4. Two things are established here and they fail differently.
!!
!! **The sum rule, in Fortran, needing nothing external.** Translate every
!! distributed moment to a common origin, add, and the total molecular moments
!! must come back -- the same moments `check_multipole` gets straight from
!! `Tr(D r)`. That is an exact identity and it catches a partition that loses or
!! double-counts density, which is the failure the primitive-pair assignment is
!! most likely to have. It is checked here rather than in Python because it needs
!! no reference at all.
!!
!! **The distribution itself, against GAMESS.** The sum rule is blind to where
!! the density went: move a dipole from one atom to another and the total is
!! unchanged. Only the `.efp` comparison in the Python constrains that, and only
!! it constrains the component ordering -- a transposed quadrupole sums to a
!! transposed total, which is exactly how the polarizability transpose survived
!! until a per-point comparison found it.
program check_dma
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_dma, only: dma_result_t, distributed_multipoles, N_QUAD
   use mqc_error, only: error_t
   implicit none

   !> Bohr per Angstrom **as GAMESS converts it**, not as CODATA now has it.
   !>
   !> GAMESS uses 0.52917724924 Angstrom per Bohr, an older CODATA value; the
   !> current one is 0.529177210903, which is what the rest of this repository
   !> uses. The difference is 7e-8 relative -- invisible in an energy, and not
   !> invisible here. Converting water's geometry with our constant put every
   !> coordinate 1.3e-7 Bohr from GAMESS's, which is a *different molecule* by the
   !> standards of a parameter comparison and would have been charged against the
   !> partition instead of against arithmetic.
   !>
   !> So a program that compares against a GAMESS potential converts the way
   !> GAMESS does. `0.9584 * this` reproduces their printed 1.8111133866 exactly.
   !> Programs that compare against PySCF must *not* use this -- PySCF is on the
   !> current value, and this is the same trap as feeding it a different basis
   !> table.
   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp
   integer :: failures

   failures = 0
   call one_case(1, "6-31g*", "water")
   call one_case(2, "6-31g*", "hcn")

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[dma] wrote every case; run validation/check_dma.py"
   else
      write (*, "(A,I0,A)") "[dma] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine one_case(case_index, basis, molecule)
      integer, intent(in) :: case_index
      character(len=*), intent(in) :: basis
      character(len=*), intent(in) :: molecule

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(dma_result_t) :: dma
      type(error_t) :: err
      real(dp), allocatable :: dip(:, :, :), quad(:, :, :)
      real(dp) :: c(3, 3), total_dipole(3), total_second(3, 3)
      real(dp) :: summed_dipole(3), summed_second(3, 3)
      real(dp) :: electrons, net_charge, d_gap, s_gap
      integer :: z(3)
      integer :: nelec, unit, k, a, b, i
      character(len=2) :: symbols(3)
      character(len=64) :: name
      ! Which (row, column) of the second-moment tensor each packed
      ! quadrupole component is, in the file's XX YY ZZ XY XZ YZ order.
      integer, parameter :: PACK_A(N_QUAD) = [1, 2, 3, 1, 1, 2]
      integer, parameter :: PACK_B(N_QUAD) = [1, 2, 3, 2, 3, 3]

      select case (molecule)
      case ("water")
         z = [8, 1, 1]
         symbols = ["O ", "H ", "H "]
         nelec = 10
         c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                      0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                      0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      case default
         z = [7, 6, 1]
         symbols = ["N ", "C ", "H "]
         nelec = 14
         c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                      0.0_dp, 0.0_dp, 1.1560_dp*ANG, &
                      0.0_dp, 0.0_dp, 2.2230_dp*ANG], [3, 3])
      end select

      call build_libcint_molecule(z, symbols, c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[dma] basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A)") "[dma] SCF failed"
         failures = failures + 1
         call mol%destroy()
         return
      end if

      call distributed_multipoles(mol, scf%density, z, dma, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[dma] failed: ", err%get_message()
         failures = failures + 1
         call mol%destroy()
         return
      end if

      ! The totals, straight from the contracted basis about the origin.
      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, err)
      if (.not. err%has_error()) then
         call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 2, quad, err)
      end if
      if (err%has_error()) then
         write (*, "(A)") "[dma] multipole integrals failed"
         failures = failures + 1
         call mol%destroy()
         return
      end if
      do k = 1, 3
         total_dipole(k) = -sum(scf%density*dip(:, :, k)) &
                           + sum(mol%charges*mol%coords(k, :))
      end do
      do b = 1, 3
         do a = 1, 3
            total_second(a, b) = -sum(scf%density*quad(:, :, 3*(a - 1) + b)) &
                                 + sum(mol%charges*mol%coords(a, :)*mol%coords(b, :))
         end do
      end do

      ! And the same totals assembled from the distributed set, which is the
      ! identity being tested.
      electrons = -sum(dma%electronic)
      net_charge = sum(dma%electronic + dma%nuclear)
      summed_dipole = 0.0_dp
      summed_second = 0.0_dp
      do k = 1, size(dma%points, 2)
         summed_dipole = summed_dipole &
                         + (dma%electronic(k) + dma%nuclear(k))*dma%points(:, k) &
                         + dma%dipole(:, k)
         do i = 1, 6
            a = PACK_A(i)
            b = PACK_B(i)
            summed_second(a, b) = summed_second(a, b) + dma%quadrupole(i, k)
            if (a /= b) summed_second(b, a) = summed_second(b, a) + dma%quadrupole(i, k)
         end do
         do b = 1, 3
            do a = 1, 3
               summed_second(a, b) = summed_second(a, b) &
                                     + (dma%electronic(k) + dma%nuclear(k)) &
                                     *dma%points(a, k)*dma%points(b, k) &
                                     + dma%dipole(a, k)*dma%points(b, k) &
                                     + dma%dipole(b, k)*dma%points(a, k)
            end do
         end do
      end do

      d_gap = maxval(abs(summed_dipole - total_dipole))
      s_gap = maxval(abs(summed_second - total_second))

      write (name, "(A,I0,A)") "/tmp/mqc_dma_", case_index, ".txt"
      open (newunit=unit, file=trim(name), status="replace", action="write")
      write (unit, "(A)") trim(basis)
      write (unit, "(A)") trim(molecule)
      write (unit, "(I0)") size(dma%points, 2)
      do k = 1, size(dma%points, 2)
         write (unit, "(A)") trim(dma%labels(k))
      end do
      write (unit, "(3es25.16e3)") dma%points
      write (unit, "(es25.16e3)") dma%electronic
      write (unit, "(es25.16e3)") dma%nuclear
      write (unit, "(3es25.16e3)") dma%dipole
      write (unit, "(es25.16e3)") dma%quadrupole
      write (unit, "(es25.16e3)") dma%octopole
      close (unit)

      write (*, "(A,A6,1X,A8,A,I0,A,F14.9,A,ES9.2)") &
         "  ", trim(molecule), trim(basis), "  points ", size(dma%points, 2), &
         "   electrons ", electrons, "   net charge ", net_charge
      write (*, "(A,ES10.2,A,ES10.2)") &
         "        sum rule: dipole ", d_gap, "   second moment ", s_gap

      if (abs(electrons - real(nelec, dp)) > 1.0e-9_dp) then
         write (*, "(A)") "    FAIL the partition does not account for every electron"
         failures = failures + 1
      end if
      if (d_gap > 1.0e-9_dp) then
         write (*, "(A)") "    FAIL the distributed dipoles do not sum to the total"
         failures = failures + 1
      end if
      if (s_gap > 1.0e-8_dp) then
         write (*, "(A)") "    FAIL the distributed quadrupoles do not sum to the total"
         failures = failures + 1
      end if

      call mol%destroy()
      deallocate (dip, quad)
   end subroutine one_case

end program check_dma
