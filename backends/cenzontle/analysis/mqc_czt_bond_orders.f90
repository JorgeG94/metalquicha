!! Mayer bond orders over a converged ab initio density
module mqc_czt_bond_orders
   !! What this is, and what it is not.
   !!
   !! Mayer's bond order is a property of the density and the overlap: the
   !! block of `D S` spanning two atoms, squared and summed. It needs no
   !! parameters, no reference molecule and no distance criterion, so it says
   !! what the wave function that was just converged thinks the bonding is --
   !! in the basis it was converged in, which is the one thing it does depend
   !! on.
   !!
   !! Three things in this tree are called bond orders and none of them agree
   !! beyond a trend:
   !!
   !!   * **these**, from the ab initio density in the AO basis;
   !!   * the **xTB Wiberg-Mayer** orders (`mqc_method_xtb`), from a
   !!     semi-empirical Hamiltonian in its own minimal basis. Cheap enough to
   !!     drive a fragmentation search, which is what they are for;
   !!   * the **QUAO kinetic** bond orders (`mqc_czt_quao`), Ruedenberg's
   !!     definition in the orthonormal quasi-atomic basis, which is a
   !!     different quantity rather than a different approximation to this one.
   !!
   !! Comparing the first two on a system small enough to afford both is the
   !! reason this exists: it says whether the cheap ranking can be trusted.
   !!
   !! **No Wiberg here.** Wiberg's order is this sum of squares in an
   !! *orthonormal* basis, where `S` is the identity and the product collapses.
   !! Computed over non-orthogonal AOs it is not Wiberg's quantity and not
   !! anybody else's, so it is not offered. Over the quasi-atomic basis, which
   !! is orthonormal, `mqc_czt_quao`'s `population_bond_order` already is it.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_charges, only: ao_to_atom
   use mqc_population_analysis, only: mayer_atomic_bond_orders, &
                                      mayer_atomic_bond_orders_open_shell, &
                                      mayer_atomic_valences
   implicit none
   private

   public :: mayer_bond_orders
   public :: mayer_bond_orders_open_shell
   public :: mayer_valences
   public :: print_mayer_report

   ! Pairs weaker than this are not printed. Not a cutoff on the numbers --
   ! the JSON carries the whole matrix -- only on the table, which is
   ! unreadable at N^2 rows and where everything below a tenth is either a
   ! weak interaction the table cannot rank or basis-set noise.
   real(dp), parameter :: REPORT_THRESHOLD = 0.1_dp

contains

   subroutine mayer_bond_orders(mol, density, overlap, orders, error)
      !! B_AB from a closed-shell total density
      !!
      !! `density` is the *total* density, which for a restricted build
      !! already carries its factor of two. An unrestricted density must not
      !! come through here even summed -- see `mayer_bond_orders_open_shell`.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :), overlap(:, :)
      real(dp), allocatable, intent(out) :: orders(:, :)
      type(error_t), intent(inout) :: error

      integer, allocatable :: owner(:)

      call ao_to_atom(mol, owner)
      call mayer_atomic_bond_orders(owner, mol%natm, density, overlap, orders, error)
   end subroutine mayer_bond_orders

   subroutine mayer_bond_orders_open_shell(mol, density_alpha, density_beta, overlap, &
                                           orders, error)
      !! B_AB from two spin densities, which is not the closed-shell formula
      !!
      !! It reduces to it when the two are equal, so this is the general entry
      !! and the other is the cheaper special case rather than a different
      !! definition.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density_alpha(:, :), density_beta(:, :)
      real(dp), intent(in) :: overlap(:, :)
      real(dp), allocatable, intent(out) :: orders(:, :)
      type(error_t), intent(inout) :: error

      integer, allocatable :: owner(:)

      call ao_to_atom(mol, owner)
      call mayer_atomic_bond_orders_open_shell(owner, mol%natm, density_alpha, &
                                               density_beta, overlap, orders, error)
   end subroutine mayer_bond_orders_open_shell

   subroutine mayer_valences(orders, valences)
      !! V_A = sum_{B /= A} B_AB
      !!
      !! Re-exported from `mqc_population_analysis` so a caller that has the
      !! matrix from here does not have to reach past this module for the one
      !! line that sums it.
      real(dp), intent(in) :: orders(:, :)
      real(dp), allocatable, intent(out) :: valences(:)

      call mayer_atomic_valences(orders, valences)
   end subroutine mayer_valences

   subroutine print_mayer_report(orders, valences, element_symbols, threshold)
      !! The bonded pairs, then the valence each atom's orders add up to
      !!
      !! Printed whatever the verbosity, like every other analysis a deck asks
      !! for by name: a run that was asked for bond orders and says nothing has
      !! reported that it found no bonds.
      real(dp), intent(in) :: orders(:, :)
      real(dp), intent(in) :: valences(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in), optional :: threshold   !! Weakest pair to print

      character(len=160) :: line
      character(len=16) :: label, partner
      real(dp) :: cutoff
      integer :: natm, iatom, jatom, shown

      cutoff = REPORT_THRESHOLD
      if (present(threshold)) cutoff = threshold
      natm = size(valences)

      call logger%info("")
      call logger%info("  ==== Mayer bond orders =================================")
      call logger%info("")
      ! The basis is not named here because this routine does not know it; the
      ! header of the calculation does. It matters: a Mayer order moves with
      ! the basis, if far less violently than a Mulliken charge.
      write (line, "(4x,a,f6.2,a)") "pairs above ", cutoff, ":"
      call logger%info(trim(line))
      shown = 0
      do iatom = 1, natm
         do jatom = iatom + 1, natm
            if (orders(iatom, jatom) < cutoff) cycle
            write (label, "(a,i0)") trim(adjustl(element_symbols(iatom)))//" ", iatom
            write (partner, "(a,i0)") trim(adjustl(element_symbols(jatom)))//" ", jatom
            write (line, "(6x,a8,a3,a8,f12.4)") label, "--", partner, &
               orders(iatom, jatom)
            call logger%info(trim(line))
            shown = shown + 1
         end do
      end do
      if (shown == 0) call logger%info("      none")

      call logger%info("")
      call logger%info("     atom         valence")
      do iatom = 1, natm
         write (label, "(a,i0)") trim(adjustl(element_symbols(iatom)))//" ", iatom
         write (line, "(4x,a8,f12.4)") label, valences(iatom)
         call logger%info(trim(line))
      end do
      call logger%info("")
   end subroutine print_mayer_report

end module mqc_czt_bond_orders
