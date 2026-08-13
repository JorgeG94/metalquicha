!! Does the MAKEFP machinery converge as the basis grows
!!
!!     ./build/check_makefp_scaling
!!
!! Everything a fragment potential needs *computed*, at a series of basis sets, with
!! the writing left out. The write needs a map onto GAMESS's orbital ordering that
!! exists for Cartesian s, p and d only, so it would refuse a Dunning set before any
!! of the physics ran -- and whether the physics holds up is the question here.
!!
!! **Why this is worth its own program.** The dynamic response is solved at twelve
!! imaginary frequencies for up to nine perturbations, each an iterative solve with
!! another iterative solve inside it. That is where the cost is and where convergence
!! fails first: a preconditioner mismatched at the top of the frequency range passed
!! water in 6-31G* and failed the same molecule in cc-pVTZ. Watching one basis prove
!! nothing about the next is the lesson; this walks the series so a regression shows
!! up as the size where it stops rather than as a surprise later.
!!
!! **Where this stands against GAMESS.** With the Hessian stored, water in cc-pVQZ
!! takes 5139 CPU seconds and about 190 seconds of wall time on forty threads. GAMESS
!! does the same potential in 1108 seconds of wall time on one core plus a spinning
!! helper, 1665 CPU seconds. So we finish 5.8 times sooner and use 3.1 times more
!! processor to do it: the parallel scaling is carrying us, and its integrals are
!! genuinely better per cycle.
!!
!! Two attempts at closing the CPU gap failed on measurement and are worth recording
!! so they are not retried. Solving the response in product form, which is what GAMESS
!! does, is 16 times *slower* here -- squaring the operator squares its condition
!! number, and a Krylov method pays for that immediately, whereas forming the same
!! product once as a matrix costs nothing. And building both `(A+B)` and `(A-B)` from
!! one pass over the integrals, by routing the symmetric half through the build that
!! writes all eight permutations, cost 1.6 times more rather than less: it halves
!! integral generation and more than gives that back in scatter, because the folded
!! build does six updates with a degeneracy factor where the explicit one does eight.
!! Getting that win needs a build that folds one half and expands the other in the
!! same quartet loop.
!!
!! Reported per basis: orbital count, SCF energy, the isotropic static and
!! highest-frequency polarizabilities, and the wall time. No reference values -- the
!! polarizabilities should fall smoothly towards a limit and the run should finish,
!! which is what a convergence check is.
program check_makefp_scaling
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_localize, only: boys_localize
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_dma, only: dma_result_t, distributed_multipoles
   use mqc_libcint_cphf, only: distributed_polarizability, &
                               distributed_dynamic_polarizability, &
                               distributed_dynamic_cross, &
                               casimir_polder_frequencies, N_CASIMIR_POLDER
   use mqc_libcint_screening, only: fit_screening, SCREEN_EXPONENTIAL
   use mqc_error, only: error_t
   implicit none

   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp
   !> libcint's full-Cartesian quadrupole slots, xx,xy,xz,yx,...,zz.
   integer, parameter :: QXX = 1, QXY = 2, QXZ = 3, QYY = 5, QYZ = 6, QZZ = 9

   character(len=16), parameter :: BASES(4) = [character(len=16) :: &
                                               "6-31g*", "cc-pvdz", "cc-pvtz", "cc-pvqz"]
   integer :: failures, ib
   real(dp) :: c(3, 3)
   integer :: z(3)
   character(len=2) :: symbols(3)

   z = [8, 1, 1]
   symbols = ["O ", "H ", "H "]
   c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])

   failures = 0
   write (*, "(A)") "  water, the geometry the GAMESS references use"
   write (*, "(A)") ""
   write (*, "(A8,A6,A18,A13,A13,A9)") "basis", "nao", "SCF energy", &
      "alpha(0)", "alpha(nu_max)", "seconds"
   do ib = 1, size(BASES)
      call one_basis(trim(BASES(ib)))
   end do

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[scaling] every basis completed the whole response set"
   else
      write (*, "(A,I0,A)") "[scaling] ", failures, " basis set(s) failed"
      error stop 1
   end if

contains

   subroutine one_basis(basis)
      character(len=*), intent(in) :: basis

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(dma_result_t) :: dma
      type(error_t) :: err
      real(dp), allocatable :: loc(:, :), cen(:, :), stat_pol(:, :, :)
      real(dp), allocatable :: dyn(:, :, :, :), cross(:, :, :, :)
      real(dp), allocatable :: dip(:, :, :), quad(:, :, :), buck(:, :, :)
      real(dp), allocatable :: alpha_scr(:)
      real(dp) :: nu(N_CASIMIR_POLDER), com(3), mass_total
      real(dp) :: iso0, isomax, t0, t1
      integer :: n_lmo, k, i

      call cpu_time(t0)
      call build_libcint_molecule(z, symbols, c, basis, mol, err)
      if (bail(err, basis, "basis")) return

      call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (bail(err, basis, "SCF")) return
      if (.not. scf%converged) then
         write (*, "(A8,A)") basis, "   SCF did not converge"
         failures = failures + 1
         call mol%destroy()
         return
      end if

      n_lmo = scf%n_occupied - 1
      call boys_localize(mol, scf%orbitals(:, 2:scf%n_occupied), n_lmo, loc, cen, err)
      if (bail(err, basis, "localization")) return

      call distributed_multipoles(mol, scf%density, z, dma, err)
      if (bail(err, basis, "distributed multipoles")) return

      call distributed_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                      scf%n_occupied, stat_pol, cen, err, n_core=1)
      if (bail(err, basis, "static polarizability")) return

      nu = casimir_polder_frequencies()
      call distributed_dynamic_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                              scf%n_occupied, nu, dyn, cen, err, &
                                              n_core=1)
      if (bail(err, basis, "dynamic dipole polarizability")) return

      ! The mixed and quadrupole-quadrupole blocks, which are the expensive ones:
      ! nine driving operators rather than three, at every frequency.
      com = 0.0_dp
      mass_total = 0.0_dp
      do i = 1, 3
         com = com + real(z(i), dp)*c(:, i)
         mass_total = mass_total + real(z(i), dp)
      end do
      com = com/mass_total
      call multipole_matrices(mol, com, 1, dip, err)
      if (bail(err, basis, "dipole integrals")) return
      call multipole_matrices(mol, com, 2, quad, err)
      if (bail(err, basis, "quadrupole integrals")) return
      allocate (buck(mol%nao, mol%nao, 9))
      buck(:, :, QXX) = 0.5_dp*(2.0_dp*quad(:, :, QXX) - quad(:, :, QYY) - quad(:, :, QZZ))
      buck(:, :, QYY) = 0.5_dp*(2.0_dp*quad(:, :, QYY) - quad(:, :, QXX) - quad(:, :, QZZ))
      buck(:, :, QZZ) = 0.5_dp*(2.0_dp*quad(:, :, QZZ) - quad(:, :, QXX) - quad(:, :, QYY))
      buck(:, :, QXY) = 1.5_dp*quad(:, :, QXY)
      buck(:, :, QXZ) = 1.5_dp*quad(:, :, QXZ)
      buck(:, :, QYZ) = 1.5_dp*quad(:, :, QYZ)
      buck(:, :, 4) = buck(:, :, QXY)
      buck(:, :, 7) = buck(:, :, QXZ)
      buck(:, :, 8) = buck(:, :, QYZ)

      call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                     scf%n_occupied, nu, buck, dip, cross, cen, err, &
                                     n_core=1)
      if (bail(err, basis, "dipole-quadrupole response")) return
      deallocate (cross)
      call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                     scf%n_occupied, nu, buck, buck, cross, cen, err, &
                                     n_core=1)
      if (bail(err, basis, "quadrupole-quadrupole response")) return

      call fit_screening(mol, scf%density, dma, z, SCREEN_EXPONENTIAL, alpha_scr, err)
      if (bail(err, basis, "screening fit")) return

      iso0 = 0.0_dp
      isomax = 0.0_dp
      do k = 1, n_lmo
         iso0 = iso0 + (stat_pol(1, 1, k) + stat_pol(2, 2, k) + stat_pol(3, 3, k))/3.0_dp
         isomax = isomax + (dyn(1, 1, k, N_CASIMIR_POLDER) &
                            + dyn(2, 2, k, N_CASIMIR_POLDER) &
                            + dyn(3, 3, k, N_CASIMIR_POLDER))/3.0_dp
      end do
      call cpu_time(t1)

      write (*, "(A8,I6,F18.8,F13.6,F13.6,F9.1)") basis, mol%nao, scf%energy, &
         iso0, isomax, t1 - t0
      flush (6)

      call mol%destroy()
      deallocate (loc, cen, stat_pol, dyn, cross, dip, quad, buck, alpha_scr)
   end subroutine one_basis

   logical function bail(err, basis, stage)
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis, stage

      bail = err%has_error()
      if (bail) then
         write (*, "(A8,A,A,A,A)") basis, "   FAILED in ", stage, ": ", &
            trim(err%get_message())
         flush (6)
         failures = failures + 1
      end if
   end function bail

end program check_makefp_scaling
