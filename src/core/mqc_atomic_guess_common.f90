!! Free-atom facts an SCF initial guess needs, independent of backend
module mqc_atomic_guess_common
   !! What a superposition-of-atomic-densities guess needs to know about a free
   !! atom before any integrals exist.
   !!
   !! The libcint and cuEST guesses build their atomic densities with completely
   !! different machinery, but they agree on which atom they are building: the
   !! occupancy follows from the periodic table, not from the backend. Both
   !! carried identical copies of this; one copy means the two paths cannot
   !! start from different atoms and blame each other for the difference.
   use pic_types, only: dp
   implicit none
   private

   public :: hund_multiplicity

contains

   pure function hund_multiplicity(atomic_number) result(multiplicity)
      !! Ground-state spin multiplicity of a neutral free atom
      !!
      !! Aufbau filling in Madelung order, with Hund's first rule inside each
      !! subshell: electrons occupy distinct orbitals with parallel spin before
      !! pairing. Good enough for a guess. It gets the transition metals
      !! nominally wrong wherever the real ground state defies Madelung -- Cr and
      !! Cu are the standard examples -- which costs iterations rather than
      !! correctness, because the molecular SCF is what decides the answer.
      integer, intent(in) :: atomic_number
      integer :: multiplicity

      ! Subshell capacities in Madelung (n+l, then n) order: 1s 2s 2p 3s 3p 4s
      ! 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s. Enough for the whole periodic table, which
      ! is more than the basis sets here cover.
      integer, parameter :: N_SUBSHELLS = 16
      integer, parameter :: CAPACITY(N_SUBSHELLS) = &
                            [2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2]
      integer :: remaining, i, in_shell, degeneracy, unpaired

      remaining = atomic_number
      unpaired = 0
      do i = 1, size(CAPACITY)
         if (remaining <= 0) exit
         degeneracy = CAPACITY(i)/2
         in_shell = min(remaining, CAPACITY(i))
         remaining = remaining - in_shell
         ! Hund: singly occupy every orbital of the subshell, then pair up.
         if (in_shell <= degeneracy) then
            unpaired = in_shell
         else
            unpaired = 2*degeneracy - in_shell
         end if
      end do
      multiplicity = unpaired + 1
   end function hund_multiplicity

end module mqc_atomic_guess_common
