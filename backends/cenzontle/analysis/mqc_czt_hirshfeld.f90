!! Hirshfeld atomic charges, by stockholder partition of the molecular density
module mqc_czt_hirshfeld
   !! Each point of space is shared between the atoms in proportion to how much
   !! density each would put there on its own:
   !!
   !!     w_A(r) = rho_A^0(|r - R_A|) / sum_B rho_B^0(|r - R_B|)
   !!     q_A    = Z_A - integral w_A(r) rho(r) dr
   !!
   !! (Hirshfeld, Theor. Chim. Acta 44, 129 (1977)). The weights depend on the
   !! geometry and the elements only, never on the density being divided, so
   !! two densities of one molecule -- a neutral and its ion -- are split by the
   !! same partition. That is what makes the Hirshfeld-condensed Fukui index a
   !! difference of densities rather than of two partitions.
   !!
   !! **The free atoms are one fixed recipe, whatever the molecule ran.** Each
   !! element is the spherically averaged, spin-polarised atom of
   !! `mqc_czt_spherical_atom`, solved with wB97X in def2-TZVP on grid level 7.
   !! A charge is therefore a property of the molecular density alone and does
   !! not move with the molecule's basis or method. A fractional-occupation
   !! Hartree-Fock atom would not do: its fractional electrons keep part of
   !! their self-repulsion, and oxygen sits 0.21 hartree above the UHF atom
   !! where PBE's is 0.013 above. The free atoms follow the approach of the
   !! PySCF extension github.com/frobnitzem/hirshfeld.
   !!
   !! Each free atom is solved once per process and kept as a table of its
   !! radial density, interpolated by cubic Hermite in (ln r, ln rho) with the
   !! slope taken from the analytic density gradient.
   !!
   !! Needs libxc, for the functional. Covers H-Kr, the elements the free-atom
   !! configurations cover. Refuses an effective core potential, whose
   !! missing core electrons an all-electron free atom would still count.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_number_to_symbol
   use mqc_string_utils, only: int_to_text
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_ao, only: eval_ao_block, eval_rho, shell_extents, block_significant_aos, &
                         AO_POINT_BLOCK
   use mqc_czt_xc, only: xc_available
   use mqc_czt_spherical_atom, only: spherical_atom_settings_t, spherical_atom_t, &
                                     solve_spherical_atom
   implicit none
   private

   public :: hirshfeld_populations     !! Charges, and spin populations for an open shell
   public :: free_atom_density         !! A cached free atom's density at given radii
   public :: clear_free_atom_cache     !! Drop every cached free atom
   public :: HIRSHFELD_GRID_LEVEL

   integer, parameter :: HIRSHFELD_GRID_LEVEL = 5
      !! Quadrature level for the molecular integral, when the caller names none.

   character(len=*), parameter :: FREE_ATOM_BASIS = "def2-tzvp"
   character(len=*), parameter :: FREE_ATOM_FUNCTIONAL = "wb97x"
   integer, parameter :: FREE_ATOM_GRID_LEVEL = 7
      !! wB97X's stiff semilocal part is the slowest functional here to settle on
      !! a grid; at level 7 its free atoms agree with PySCF's on 200x1202 to
      !! 2e-11 hartree.

   integer, parameter :: Z_MAX = 36
      !! Krypton, the last element with a free-atom configuration.

   integer, parameter :: N_RADIAL = 2000
   real(dp), parameter :: R_MIN = 1.0e-6_dp
   real(dp), parameter :: R_MAX = 60.0_dp
      !! The table's span, Bohr, uniform in ln r. Below `R_MIN` a Gaussian
      !! density is flat to far better than the table's precision.
   real(dp), parameter :: RHO_FLOOR = 1.0e-30_dp
      !! A free atom is zero beyond the radius where it falls below this.

   real(dp), parameter :: AO_SCREEN = 1.0e-12_dp
      !! AO value below which a shell is dropped from a grid block, as for the
      !! exchange-correlation quadrature.

   type :: free_atom_table_t
      ! One element's radial density, uniform in ln r.
      logical :: valid = .false.
      real(dp) :: log_r0 = 0.0_dp
      real(dp) :: step = 0.0_dp
      real(dp) :: r_cut = 0.0_dp
      real(dp), allocatable :: log_rho(:)   ! (N_RADIAL)
      real(dp), allocatable :: slope(:)     ! d ln rho / d ln r
   end type free_atom_table_t

   type(free_atom_table_t), save :: tables(Z_MAX)

contains

   subroutine hirshfeld_populations(mol, density, charges, error, spin_density, &
                                    spin_populations, grid_level)
      !! q_A = Z_A - N_A, with N_A the stockholder share of `density`
      !!
      !! `density` is the total, both spins. With `spin_density` (alpha minus
      !! beta) the same weights give `spin_populations`. A ghost centre carries
      !! no free atom, so it takes no share and its charge is zero.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)              !! (nao, nao)
      real(dp), allocatable, intent(out) :: charges(:)   !! (natm), e
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: spin_density(:, :)
      real(dp), allocatable, intent(out), optional :: spin_populations(:)
      integer, intent(in), optional :: grid_level        !! Default `HIRSHFELD_GRID_LEVEL`

      type(dft_grid_t) :: grid
      real(dp), allocatable :: electrons(:), spins(:), extents(:)
      logical, allocatable :: real_atom(:)
      integer :: iatom, level
      logical :: want_spin

      if (error%has_error()) return
      want_spin = present(spin_density) .and. present(spin_populations)
      if (.not. xc_available()) then
         call error%set(ERROR_VALIDATION, "Hirshfeld charges need libxc, which their free "// &
                        "atoms are solved with: configure with -DMQC_ENABLE_LIBXC=ON, or "// &
                        "choose another scheme ('chelpg' or 'mulliken')")
         return
      end if
      if (allocated(mol%core_electrons)) then
         if (any(mol%core_electrons > 0)) then
            call error%set(ERROR_VALIDATION, "Hirshfeld charges: not available with an "// &
                           "effective core potential")
            return
         end if
      end if

      allocate (real_atom(mol%natm))
      do iatom = 1, mol%natm
         real_atom(iatom) = mol%charges(iatom) > 0.0_dp
         if (.not. real_atom(iatom)) cycle
         call ensure_free_atom(mol%atomic_numbers(iatom), error)
         if (error%has_error()) then
            call error%add_context("hirshfeld charges")
            return
         end if
      end do

      level = HIRSHFELD_GRID_LEVEL
      if (present(grid_level)) level = grid_level
      call build_dft_grid(mol%coords, mol%atomic_numbers, grid, error, level=level)
      if (error%has_error()) return
      call shell_extents(mol, AO_SCREEN, extents)

      allocate (electrons(mol%natm), spins(mol%natm))
      if (want_spin) then
         call integrate_shares(mol, grid, extents, real_atom, density, electrons, error, &
                               spin_density, spins)
      else
         call integrate_shares(mol, grid, extents, real_atom, density, electrons, error)
      end if
      call grid%destroy()
      if (error%has_error()) return

      allocate (charges(mol%natm))
      charges = 0.0_dp
      where (real_atom) charges = mol%charges - electrons
      if (want_spin) spin_populations = spins
   end subroutine hirshfeld_populations

   subroutine integrate_shares(mol, grid, extents, real_atom, density, electrons, error, &
                               spin_density, spins)
      ! Each atom's stockholder share of the density, and of the spin density
      ! when one is given, block by block over the grid.
      type(czt_molecule_t), intent(in) :: mol
      type(dft_grid_t), intent(in) :: grid
      real(dp), intent(in) :: extents(:)
      logical, intent(in) :: real_atom(:)
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(out) :: electrons(:)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: spin_density(:, :)
      real(dp), intent(out), optional :: spins(:)

      real(dp), allocatable :: ao(:, :), rho(:), spin_rho(:), d_sig(:, :), s_sig(:, :)
      real(dp), allocatable :: share(:, :), local_e(:), local_s(:), r_cut(:)
      integer, allocatable :: ao_list(:), ao_offset(:), near(:)
      logical, allocatable :: shell_mask(:)
      type(error_t) :: local_error
      integer :: g0, g1, n_sig, ia, ja, n_near
      logical :: failed, want_spin

      want_spin = present(spin_density) .and. present(spins)
      electrons = 0.0_dp
      if (want_spin) spins = 0.0_dp
      allocate (r_cut(mol%natm))
      r_cut = 0.0_dp
      do ia = 1, mol%natm
         if (real_atom(ia)) r_cut(ia) = tables(mol%atomic_numbers(ia))%r_cut
      end do
      failed = .false.

      ! `local_error` firstprivate rather than private, for the reason
      ! `xc_add_potential` gives: a private copy need not be default-initialised.
      !$omp parallel default(shared) &
      !$omp    private(g0, g1, n_sig, ia, ja, n_near, ao, rho, spin_rho, d_sig, &
      !$omp            s_sig, share, local_e, local_s, ao_list, ao_offset, near, shell_mask) &
      !$omp    firstprivate(local_error)
      allocate (d_sig(mol%nao, mol%nao), shell_mask(mol%nbas), ao_offset(mol%nbas), &
                ao_list(mol%nao), near(mol%natm), local_e(mol%natm), local_s(mol%natm))
      if (want_spin) allocate (s_sig(mol%nao, mol%nao))
      local_e = 0.0_dp
      local_s = 0.0_dp

      !$omp do schedule(dynamic)
      do g0 = 1, grid%n_points, AO_POINT_BLOCK
         if (failed) cycle
         g1 = min(g0 + AO_POINT_BLOCK - 1, grid%n_points)

         call block_significant_aos(mol, grid%coords(:, g0:g1), extents, shell_mask, &
                                    ao_list, ao_offset, n_sig)
         if (n_sig == 0) cycle
         call atoms_reaching(mol, grid%coords(:, g0:g1), r_cut, near, n_near)
         if (n_near == 0) cycle

         do ja = 1, n_sig
            do ia = 1, n_sig
               d_sig(ia, ja) = density(ao_list(ia), ao_list(ja))
            end do
         end do
         call eval_ao_block(mol, grid%coords(:, g0:g1), ao, local_error, &
                            shell_mask=shell_mask, ao_offset=ao_offset, n_ao_out=n_sig)
         if (local_error%has_error()) then
            !$omp critical (mqc_hirshfeld_error)
            if (.not. failed) error = local_error
            failed = .true.
            !$omp end critical (mqc_hirshfeld_error)
            cycle
         end if
         call eval_rho(ao, d_sig(1:n_sig, 1:n_sig), rho)

         call stockholder_shares(mol, grid%coords(:, g0:g1), near(1:n_near), share)
         do ia = 1, n_near
            local_e(near(ia)) = local_e(near(ia)) + &
                                sum(grid%weights(g0:g1)*rho*share(:, ia))
         end do

         if (want_spin) then
            do ja = 1, n_sig
               do ia = 1, n_sig
                  s_sig(ia, ja) = spin_density(ao_list(ia), ao_list(ja))
               end do
            end do
            call eval_rho(ao, s_sig(1:n_sig, 1:n_sig), spin_rho)
            do ia = 1, n_near
               local_s(near(ia)) = local_s(near(ia)) + &
                                   sum(grid%weights(g0:g1)*spin_rho*share(:, ia))
            end do
         end if
      end do
      !$omp end do

      !$omp critical (mqc_hirshfeld_sum)
      electrons = electrons + local_e
      if (want_spin) spins = spins + local_s
      !$omp end critical (mqc_hirshfeld_sum)
      !$omp end parallel
   end subroutine integrate_shares

   subroutine atoms_reaching(mol, coords, r_cut, near, n_near)
      ! The atoms whose free density reaches any point of the block, by the
      ! block's bounding sphere: one distance per atom rather than per point.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)
      real(dp), intent(in) :: r_cut(:)
      integer, intent(out) :: near(:)
      integer, intent(out) :: n_near

      real(dp) :: centre(3)
      real(dp) :: radius
      integer :: ia, ig

      centre = sum(coords, dim=2)/real(size(coords, 2), dp)
      radius = 0.0_dp
      do ig = 1, size(coords, 2)
         radius = max(radius, norm2(coords(:, ig) - centre))
      end do

      n_near = 0
      do ia = 1, mol%natm
         if (r_cut(ia) <= 0.0_dp) cycle
         if (norm2(mol%coords(:, ia) - centre) > r_cut(ia) + radius) cycle
         n_near = n_near + 1
         near(n_near) = ia
      end do
   end subroutine atoms_reaching

   subroutine stockholder_shares(mol, coords, near, share)
      ! w_A at every point of the block, for the atoms in `near`. Zero at a
      ! point no free atom reaches, which the molecular density cannot reach
      ! either to any precision that matters.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: near(:)
      real(dp), allocatable, intent(out) :: share(:, :)   ! (n_points, size(near))

      real(dp) :: total
      integer :: ig, ia

      allocate (share(size(coords, 2), size(near)))
      do ig = 1, size(coords, 2)
         do ia = 1, size(near)
            share(ig, ia) = table_density(tables(mol%atomic_numbers(near(ia))), &
                                          norm2(coords(:, ig) - mol%coords(:, near(ia))))
         end do
         total = sum(share(ig, :))
         if (total > 0.0_dp) then
            share(ig, :) = share(ig, :)/total
         else
            share(ig, :) = 0.0_dp
         end if
      end do
   end subroutine stockholder_shares

   subroutine free_atom_density(atomic_number, radii, rho, error)
      !! The free atom's total density at each radius, Bohr in, e/Bohr^3 out
      !!
      !! Solves and caches the atom on first use, as a charge calculation would.
      integer, intent(in) :: atomic_number
      real(dp), intent(in) :: radii(:)
      real(dp), allocatable, intent(out) :: rho(:)
      type(error_t), intent(inout) :: error

      integer :: k

      call ensure_free_atom(atomic_number, error)
      if (error%has_error()) return
      allocate (rho(size(radii)))
      do k = 1, size(radii)
         rho(k) = table_density(tables(atomic_number), radii(k))
      end do
   end subroutine free_atom_density

   subroutine ensure_free_atom(atomic_number, error)
      ! Solve and tabulate one element unless it is already cached. The check
      ! and the fill are one critical section, so two threads missing on the
      ! same element cannot both write it.
      integer, intent(in) :: atomic_number
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      if (atomic_number < 1 .or. atomic_number > Z_MAX) then
         call error%set(ERROR_VALIDATION, "no free atom for Z="// &
                        int_to_text(atomic_number)//"; Hirshfeld covers H-Kr")
         return
      end if

      !$omp critical (mqc_hirshfeld_cache)
      if (.not. tables(atomic_number)%valid) then
         call tabulate_free_atom(atomic_number, tables(atomic_number), error)
      end if
      !$omp end critical (mqc_hirshfeld_cache)
   end subroutine ensure_free_atom

   subroutine tabulate_free_atom(atomic_number, table, error)
      ! Converge the free atom and record its radial density and log slope.
      integer, intent(in) :: atomic_number
      type(free_atom_table_t), intent(inout) :: table
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: atom_mol
      type(spherical_atom_settings_t) :: settings
      type(spherical_atom_t) :: atom
      real(dp), allocatable :: points(:, :), ao(:, :), ao_grad(:, :, :), rho(:), rho_grad(:, :)
      real(dp) :: origin(3, 1)
      integer :: k, last

      origin = 0.0_dp
      call build_czt_molecule([atomic_number], [element_number_to_symbol(atomic_number)], &
                              origin, FREE_ATOM_BASIS, atom_mol, error)
      if (error%has_error()) then
         ! The likeliest cause is a build whose MQC_BASIS_SETS left the free
         ! atoms' basis out, which nothing else would point at.
         call error%add_context("Hirshfeld's free atoms are solved in "//FREE_ATOM_BASIS// &
                                ", which has to be among the build's MQC_BASIS_SETS")
         return
      end if

      settings%functional = FREE_ATOM_FUNCTIONAL
      settings%grid_level = FREE_ATOM_GRID_LEVEL
      settings%polarized = .true.
      call solve_spherical_atom(atom_mol, settings, atom, error)
      if (error%has_error()) then
         call atom_mol%destroy()
         return
      end if
      if (.not. atom%converged) then
         call error%set(ERROR_VALIDATION, "the free atom Z="//int_to_text(atomic_number)// &
                        " did not converge")
         call atom_mol%destroy()
         return
      end if

      table%log_r0 = log(R_MIN)
      table%step = (log(R_MAX) - log(R_MIN))/real(N_RADIAL - 1, dp)
      allocate (points(3, N_RADIAL))
      points = 0.0_dp
      do k = 1, N_RADIAL
         points(3, k) = exp(table%log_r0 + real(k - 1, dp)*table%step)
      end do
      call eval_ao_block(atom_mol, points, ao, error, grad=ao_grad)
      call atom_mol%destroy()
      if (error%has_error()) return
      call eval_rho(ao, atom%d_alpha + atom%d_beta, rho, ao_grad=ao_grad, rho_grad=rho_grad)

      ! Out to the last radius above the floor; the density falls monotonically
      ! there, so nothing beyond it comes back.
      last = 1
      do k = 1, N_RADIAL
         if (rho(k) > RHO_FLOOR) last = k
      end do
      table%r_cut = points(3, last)

      allocate (table%log_rho(N_RADIAL), table%slope(N_RADIAL))
      do k = 1, N_RADIAL
         if (k <= last) then
            table%log_rho(k) = log(rho(k))
            ! d ln rho / d ln r = r rho'/rho, rho' along the axis the points lie on
            table%slope(k) = points(3, k)*rho_grad(k, 3)/rho(k)
         else
            table%log_rho(k) = log(RHO_FLOOR)
            table%slope(k) = 0.0_dp
         end if
      end do
      table%valid = .true.
   end subroutine tabulate_free_atom

   pure function table_density(table, r) result(rho)
      ! Cubic Hermite in (ln r, ln rho) between the two bracketing radii.
      type(free_atom_table_t), intent(in) :: table
      real(dp), intent(in) :: r
      real(dp) :: rho

      real(dp) :: x, t, h00, h10, h01, h11
      integer :: k

      rho = 0.0_dp
      if (r >= table%r_cut) return
      x = (log(max(r, R_MIN)) - table%log_r0)/table%step
      k = min(int(x) + 1, N_RADIAL - 1)
      t = x - real(k - 1, dp)
      h00 = (1.0_dp + 2.0_dp*t)*(1.0_dp - t)**2
      h10 = t*(1.0_dp - t)**2
      h01 = t*t*(3.0_dp - 2.0_dp*t)
      h11 = t*t*(t - 1.0_dp)
      rho = exp(h00*table%log_rho(k) + h10*table%step*table%slope(k) + &
                h01*table%log_rho(k + 1) + h11*table%step*table%slope(k + 1))
   end function table_density

   subroutine clear_free_atom_cache()
      !! Drop every cached free atom, for the tests
      integer :: z

      do z = 1, Z_MAX
         if (allocated(tables(z)%log_rho)) deallocate (tables(z)%log_rho)
         if (allocated(tables(z)%slope)) deallocate (tables(z)%slope)
         tables(z)%valid = .false.
      end do
   end subroutine clear_free_atom_cache

end module mqc_czt_hirshfeld
