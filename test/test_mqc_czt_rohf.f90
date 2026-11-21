!! Restricted open-shell Hartree-Fock, against PySCF and against itself
module test_mqc_czt_rohf
   !! The pinned energies are PySCF's `scf.ROHF` on water at this geometry in
   !! 6-31G, fed this repository's basis: the cation as a doublet and the
   !! neutral molecule as a triplet.
   !!
   !! Beside them, what needs no other program:
   !!
   !!   * the two spins share their orbital spaces -- beta's occupied space
   !!     lies inside alpha's -- so `<S^2>` is exact, which an unrestricted
   !!     solution of the same state is not;
   !!   * the energy is above the unrestricted one for the same state, since
   !!     it is the same functional minimised under a constraint;
   !!   * at convergence the three off-diagonal blocks of Roothaan's
   !!     effective Fock matrix vanish in the orbital basis, which is the
   !!     stationary condition itself rather than a consequence of it;
   !!   * on a closed shell the effective matrix is the Fock matrix and the
   !!     run is restricted Hartree-Fock;
   !!   * each spin's reported orbitals diagonalise that spin's own Fock
   !!     matrix inside its occupied block, which is what lets a gradient
   !!     build its energy-weighted density from them.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf, roothaan_fock, &
                          assemble_fock_uhf
   implicit none
   private

   public :: collect_mqc_czt_rohf_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: E_TOL = 1.0e-12_dp
   real(dp), parameter :: D_TOL = 1.0e-10_dp
   real(dp), parameter :: E_CATION_ROHF = -75.5784700861_dp
   real(dp), parameter :: E_CATION_UHF = -75.5805935893_dp
   real(dp), parameter :: E_TRIPLET_ROHF = -75.7236726674_dp

contains

   subroutine collect_mqc_czt_rohf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("doublet_energy_matches_pyscf", test_doublet), &
                  new_unittest("triplet_energy_matches_pyscf", test_triplet), &
                  new_unittest("spins_share_their_orbital_spaces", test_shared_spaces), &
                  new_unittest("effective_fock_is_block_diagonal_at_convergence", &
                               test_stationary), &
                  new_unittest("occupied_blocks_are_semicanonical", test_semicanonical), &
                  new_unittest("closed_shell_is_restricted_hartree_fock", test_closed_shell), &
                  new_unittest("direct_minimisation_is_refused", test_refusals) &
                  ]
   end subroutine collect_mqc_czt_rohf_tests

   subroutine water(mol, err)
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
   end subroutine water

   subroutine cation(mol, r, err)
      !! The water cation, restricted open shell, converged tightly
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: r
      type(error_t), intent(inout) :: err

      call water(mol, err)
      if (err%has_error()) return
      call run_czt_uhf(mol, 9, 2, 200, E_TOL, D_TOL, .false., r, err, &
                       grad_tol=1.0e-9_dp, restricted_open=.true.)
   end subroutine cation

   subroutine test_doublet(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: r, u
      type(error_t) :: err

      call cation(mol, r, err)
      call check(error,.not. err%has_error() .and. r%converged, &
                 "the cation must converge: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(a,f18.10)") "    ROHF energy ", r%energy
      call check(error, abs(r%energy - E_CATION_ROHF) < 1.0e-8_dp, &
                 "the doublet's energy does not match PySCF's")
      if (allocated(error)) return
      call check(error, abs(r%spin_squared - 0.75_dp) < 1.0e-10_dp, &
                 "a restricted open-shell doublet must have <S^2> = 3/4 exactly")
      if (allocated(error)) return

      ! The same functional without the constraint cannot be higher. Pinned
      ! to PySCF's too, so "lower" is not satisfied by a wrong state.
      call run_czt_uhf(mol, 9, 2, 200, E_TOL, D_TOL, .false., u, err)
      call check(error, u%converged .and. abs(u%energy - E_CATION_UHF) < 1.0e-8_dp, &
                 "the unrestricted cation does not match PySCF's")
      if (allocated(error)) return
      call check(error, r%energy > u%energy + 1.0e-4_dp, &
                 "the restricted open-shell energy must lie above the unrestricted one")
      if (allocated(error)) return
      call check(error, u%spin_squared > 0.75_dp + 1.0e-4_dp, &
                 "the unrestricted cation should be spin contaminated")
      call mol%destroy()
   end subroutine test_doublet

   subroutine test_triplet(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: r
      type(error_t) :: err

      call water(mol, err)
      call run_czt_uhf(mol, 10, 3, 200, E_TOL, D_TOL, .false., r, err, &
                       restricted_open=.true.)
      call check(error,.not. err%has_error() .and. r%converged, &
                 "triplet water must converge: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(a,f18.10)") "    triplet ROHF energy ", r%energy
      call check(error, abs(r%energy - E_TRIPLET_ROHF) < 1.0e-8_dp, &
                 "the triplet's energy does not match PySCF's")
      if (allocated(error)) return
      call check(error, r%n_occupied == 6 .and. r%n_occupied_beta == 4, &
                 "a triplet of ten electrons is six alpha and four beta")
      if (allocated(error)) return
      call check(error, abs(r%spin_squared - 2.0_dp) < 1.0e-10_dp, &
                 "a restricted open-shell triplet must have <S^2> = 2 exactly")
      call mol%destroy()
   end subroutine test_triplet

   subroutine test_shared_spaces(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: r
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), leak(:, :)

      call cation(mol, r, err)
      call check(error,.not. err%has_error() .and. r%converged, "the cation must converge")
      if (allocated(error)) return
      call mol%overlap(s)

      ! D_a S D_b = D_b says beta's occupied space is inside alpha's. The
      ! orbitals themselves differ between the spins -- each set diagonalises
      ! its own Fock matrix -- so it is the spaces that are compared.
      leak = matmul(r%density, matmul(s, r%density_beta)) - r%density_beta
      call check(error, maxval(abs(leak)) < 1.0e-9_dp, &
                 "beta's occupied space is not inside alpha's")
      if (allocated(error)) return
      ! And the virtual spaces the other way round: every alpha virtual is a
      ! beta virtual.
      leak = matmul(transpose(r%orbitals(:, 6:)), matmul(s, r%density_beta))
      call check(error, maxval(abs(leak)) < 1.0e-9_dp, &
                 "an alpha virtual orbital overlaps beta's occupied space")
      call mol%destroy()
   end subroutine test_shared_spaces

   subroutine test_stationary(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: r
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), h(:, :), fa(:, :), fb(:, :), eff(:, :), copy(:, :)
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :)
      real(dp) :: e_elec
      integer :: n

      call cation(mol, r, err)
      call check(error,.not. err%has_error() .and. r%converged, "the cation must converge")
      if (allocated(error)) return
      n = mol%nao
      call mol%overlap(s)
      call mol%core_hamiltonian(h)
      call mol%eris(eri)
      allocate (fa(n, n), fb(n, n))
      call assemble_fock_uhf(mol, h, r%density, r%density_beta, eri, bounds, fock_a=fa, &
                             fock_b=fb, e_elec=e_elec, error=err)
      call check(error,.not. err%has_error(), "Fock build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, abs(e_elec + mol%nuclear_repulsion() - r%energy) < 1.0e-9_dp, &
                 "the reported energy is not that of the reported densities")
      if (allocated(error)) return

      eff = fa
      copy = fb
      call roothaan_fock(eff, copy, r%density, r%density_beta, s)
      call check(error, maxval(abs(eff - transpose(eff))) < 1.0e-12_dp, &
                 "the effective Fock matrix is not symmetric")
      if (allocated(error)) return
      call check(error, maxval(abs(eff - copy)) < 1.0e-14_dp, &
                 "both spins must be handed the same effective matrix")
      if (allocated(error)) return

      ! Between spaces rather than between orbitals: alpha's reported
      ! occupied orbitals mix the doubly and the singly occupied ones, since
      ! they diagonalise alpha's Fock matrix across both. `Db` projects onto
      ! the closed space, `Da - Db` onto the open one, and alpha's virtuals
      ! (orbital 6 on) are the virtuals of both spins.
      block
         real(dp), allocatable :: open_d(:, :), virt(:, :)
         real(dp) :: closed_open, closed_virtual, open_virtual

         open_d = r%density - r%density_beta
         virt = r%orbitals(:, 6:)
         closed_open = maxval(abs(matmul(r%density_beta, matmul(eff, open_d))))
         closed_virtual = maxval(abs(matmul(r%density_beta, matmul(eff, virt))))
         open_virtual = maxval(abs(matmul(open_d, matmul(eff, virt))))
         write (*, "(a,3es10.2)") "    closed-open, closed-virtual, open-virtual ", &
            closed_open, closed_virtual, open_virtual
         call check(error, closed_open < 1.0e-7_dp, &
                    "the closed-open block of the effective Fock matrix does not vanish")
         if (allocated(error)) return
         call check(error, closed_virtual < 1.0e-7_dp, &
                    "the closed-virtual block of the effective Fock matrix does not vanish")
         if (allocated(error)) return
         call check(error, open_virtual < 1.0e-7_dp, &
                    "the open-virtual block of the effective Fock matrix does not vanish")
         if (allocated(error)) return
         ! And the blocks are not empty: the open space holds one electron.
         call check(error, abs(sum(open_d*s) - 1.0_dp) < 1.0e-9_dp, &
                    "the open shell of a doublet must hold one electron")
      end block
      call mol%destroy()
   end subroutine test_stationary

   subroutine test_semicanonical(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: r
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), h(:, :), fa(:, :), fb(:, :), fmo(:, :)
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), w(:, :), dfd(:, :)
      real(dp) :: e_elec
      integer :: n, i

      call cation(mol, r, err)
      call check(error,.not. err%has_error() .and. r%converged, "the cation must converge")
      if (allocated(error)) return
      n = mol%nao
      call mol%overlap(s)
      call mol%core_hamiltonian(h)
      call mol%eris(eri)
      allocate (fa(n, n), fb(n, n))
      call assemble_fock_uhf(mol, h, r%density, r%density_beta, eri, bounds, fock_a=fa, &
                             fock_b=fb, e_elec=e_elec, error=err)
      call check(error,.not. err%has_error(), "Fock build: "//err%get_full_trace())
      if (allocated(error)) return

      ! Beta's four occupied orbitals diagonalise beta's Fock matrix among
      ! themselves, with the reported energies on the diagonal.
      fmo = matmul(transpose(r%orbitals_beta(:, 1:4)), matmul(fb, r%orbitals_beta(:, 1:4)))
      do i = 1, 4
         fmo(i, i) = fmo(i, i) - r%orbital_energies_beta(i)
      end do
      call check(error, maxval(abs(fmo)) < 1.0e-8_dp, &
                 "beta's occupied block is not diagonal with the reported energies")
      if (allocated(error)) return

      ! Which is the point: C eps C^T over the occupied orbitals is then
      ! D F D, the energy-weighted density a gradient contracts with the
      ! overlap derivative. Alpha here, five orbitals.
      allocate (w(n, n))
      w = 0.0_dp
      do i = 1, 5
         w = w + r%orbital_energies(i)*matmul(r%orbitals(:, i:i), &
                                              transpose(r%orbitals(:, i:i)))
      end do
      dfd = matmul(r%density, matmul(fa, r%density))
      call check(error, maxval(abs(w - dfd)) < 1.0e-8_dp, &
                 "alpha's energy-weighted density is not D F D")
      call mol%destroy()
   end subroutine test_semicanonical

   subroutine test_closed_shell(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: r, ro
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), fa(:, :), fb(:, :), half(:, :)

      call water(mol, err)
      call run_czt_rhf(mol, 10, 200, E_TOL, D_TOL, .false., r, err)
      call run_czt_uhf(mol, 10, 1, 200, E_TOL, D_TOL, .false., ro, err, &
                       restricted_open=.true.)
      call check(error,.not. err%has_error() .and. r%converged .and. ro%converged, &
                 "both closed-shell runs must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, abs(ro%energy - r%energy) < 1.0e-10_dp, &
                 "a closed shell run restricted open is not restricted Hartree-Fock")
      if (allocated(error)) return

      ! With no open shell there is no open projector and the two spins'
      ! matrices are one: the effective matrix is that matrix, on the
      ! occupied-virtual blocks that decide the orbitals and everywhere else.
      call mol%overlap(s)
      call mol%core_hamiltonian(fa)
      fb = fa
      half = 0.5_dp*r%density
      call roothaan_fock(fa, fb, half, half, s)
      call mol%core_hamiltonian(fb)
      call check(error, maxval(abs(fa - fb)) < 1.0e-10_dp, &
                 "with equal spins the effective matrix must be the Fock matrix itself")
      call mol%destroy()
   end subroutine test_closed_shell

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: r
      type(error_t) :: err

      call water(mol, err)
      call run_czt_uhf(mol, 9, 2, 50, E_TOL, D_TOL, .false., r, err, gdm=.true., &
                       restricted_open=.true.)
      call check(error, err%has_error(), &
                 "direct minimisation over a restricted open shell must be refused")
      if (allocated(error)) return
      call err%clear()
      call run_czt_uhf(mol, 9, 2, 50, E_TOL, D_TOL, .false., r, err, second_order=.true., &
                       restricted_open=.true.)
      call check(error, err%has_error(), &
                 "a second-order SCF over a restricted open shell must be refused")
      call mol%destroy()
   end subroutine test_refusals

end module test_mqc_czt_rohf

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_rohf, only: collect_mqc_czt_rohf_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_rohf", collect_mqc_czt_rohf_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
