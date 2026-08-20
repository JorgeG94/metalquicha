!! Point-charge nuclear repulsion and its gradient
module mqc_nuclear_repulsion
   !! `E_NN = sum_{A<B} Z_A Z_B / R_AB` and its derivative, which depend on
   !! nothing but charges and coordinates.
   !!
   !! Here because both were written twice, once per integral backend: the
   !! energy as `nuclear_repulsion_energy` in the cuEST SCF and as
   !! `molecule_nuclear_repulsion` in the libcint molecule, the gradient as a
   !! `nuclear_repulsion_gradient` in each of the two gradient modules. Four
   !! copies of a formula with no backend content in it at all.
   !!
   !! The copies had drifted in the small ways copies do -- one guarded against
   !! coincident atoms and the others did not, one took integer atomic numbers
   !! and the others real charges, one overwrote the gradient where the other
   !! accumulated. None of those differences was a decision anybody made twice.
   use pic_types, only: dp
   implicit none
   private

   public :: nuclear_repulsion
   public :: add_nuclear_repulsion_gradient

contains

   pure function nuclear_repulsion(charges, coordinates) result(energy)
      !! `sum_{A<B} Z_A Z_B / R_AB`, Hartree
      !!
      !! Charges are real rather than integer so a ghost centre, which carries
      !! basis functions and no nucleus, costs nothing to express.
      real(dp), intent(in) :: charges(:)         !! Nuclear charges
      real(dp), intent(in) :: coordinates(:, :)  !! (3, n_atoms), Bohr
      real(dp) :: energy

      integer :: a, b
      real(dp) :: r

      energy = 0.0_dp
      do a = 1, size(charges)
         do b = a + 1, size(charges)
            r = norm2(coordinates(:, a) - coordinates(:, b))
            ! Two atoms at the same point are a broken geometry rather than an
            ! infinite energy, and the caller is better placed to say so than a
            ! division here is.
            if (r > 0.0_dp) energy = energy + charges(a)*charges(b)/r
         end do
      end do
   end function nuclear_repulsion

   pure subroutine add_nuclear_repulsion_gradient(charges, coordinates, gradient)
      !! Add `dE_NN/dR` to `gradient`, Hartree/Bohr
      !!
      !! The derivative on A points away from every other nucleus: repulsion
      !! pushes atoms apart, and a gradient is the direction of steepest
      !! *ascent*, which is what sets the sign below.
      !!
      !! Accumulates rather than assigns, because every caller but one was
      !! adding this to a gradient that already had electronic terms in it. The
      !! one that was not now zeroes its own array first, where that intent is
      !! visible.
      real(dp), intent(in) :: charges(:)         !! Nuclear charges
      real(dp), intent(in) :: coordinates(:, :)  !! (3, n_atoms), Bohr
      real(dp), intent(inout) :: gradient(:, :)  !! (3, n_atoms)

      integer :: a, b
      real(dp) :: rvec(3)
      real(dp) :: r

      do a = 1, size(charges)
         do b = 1, size(charges)
            if (a == b) cycle
            rvec = coordinates(:, a) - coordinates(:, b)
            r = norm2(rvec)
            if (r <= 0.0_dp) cycle
            gradient(:, a) = gradient(:, a) - charges(a)*charges(b)*rvec/r**3
         end do
      end do
   end subroutine add_nuclear_repulsion_gradient

end module mqc_nuclear_repulsion
