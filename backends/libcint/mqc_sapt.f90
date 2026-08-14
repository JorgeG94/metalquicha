!! Symmetry-adapted perturbation theory, SAPT0
module mqc_sapt
   !! The interaction energy as a sum of named physical terms -- electrostatics,
   !! exchange, induction, dispersion -- rather than as a difference of two totals.
   !!
   !! **Everything is in the dimer-centred basis.** Both monomers carry every basis
   !! function of the dimer and only their own nuclei, which is what makes the terms
   !! counterpoise-corrected by construction. SAPT in a monomer-centred basis is a
   !! different method with different numbers, and the `dHF` correction is not even
   !! defined there.
   !!
   !! **The three molecules share one atom ordering, and that is load-bearing.**
   !! A monomer built in its own atom order rather than the dimer's spans the same
   !! space with a *permuted* AO basis, so every matrix built from one molecule and
   !! contracted against another is silently wrong. It is not a hypothetical: the
   !! PySCF prototype this is checked against was written that way first, and the
   !! symptom was `Tr(D_B S) = 5.46` where it has to be the occupied count. That
   !! trace is asserted here for the same reason.
   !!
   !! Reference values come from `validation/check_sapt0.py`, a conventional
   !! four-index SAPT0 in PySCF. psi4 cannot serve as the reference: its
   !! closed-shell SAPT module is density-fitted by construction and, in its own
   !! documentation, "cannot be run with exact integrals".
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   implicit none
   private

   public :: sapt_molecules_t
   public :: build_sapt_molecules

   type :: sapt_molecules_t
      !! The dimer and its two counterpoise-corrected monomers
      type(libcint_molecule_t) :: dimer
      type(libcint_molecule_t) :: mono_a   !! A real, B ghosted
      type(libcint_molecule_t) :: mono_b   !! B real, A ghosted
      integer :: n_atoms_a = 0
      integer :: n_atoms_b = 0
      integer :: n_elec_a = 0   !! Electrons on A, from A's nuclei alone
      integer :: n_elec_b = 0
   contains
      procedure :: destroy => sapt_molecules_destroy
   end type sapt_molecules_t

contains

   subroutine build_sapt_molecules(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, &
                                   basis_name, mols, error, charge_a, charge_b)
      !! The dimer and both ghosted monomers, all in one AO ordering
      !!
      !! The atoms are laid out A then B, once, and all three molecules use that
      !! layout -- the monomers differ from the dimer only in which atoms carry a
      !! nuclear charge. So an overlap, a density or a Fock matrix from any of the
      !! three is directly contractable with one from any other.
      integer, intent(in) :: z_a(:), z_b(:)
      character(len=*), intent(in) :: sym_a(:), sym_b(:)
      real(dp), intent(in) :: xyz_a(:, :), xyz_b(:, :)   !! (3, n), Bohr
      character(len=*), intent(in) :: basis_name
      type(sapt_molecules_t), intent(out) :: mols
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: charge_a, charge_b

      integer, allocatable :: z(:)
      real(dp), allocatable :: xyz(:, :)
      character(len=len(sym_a)), allocatable :: sym(:)
      logical, allocatable :: ghost_a(:), ghost_b(:)
      integer :: na, nb, n, k

      na = size(z_a)
      nb = size(z_b)
      n = na + nb
      if (size(sym_a) /= na .or. size(xyz_a, 2) /= na) then
         call error%set(ERROR_VALIDATION, "sapt: monomer A's arrays disagree on "// &
                        "how many atoms it has")
         return
      end if
      if (size(sym_b) /= nb .or. size(xyz_b, 2) /= nb) then
         call error%set(ERROR_VALIDATION, "sapt: monomer B's arrays disagree on "// &
                        "how many atoms it has")
         return
      end if

      allocate (z(n), xyz(3, n), sym(n), ghost_a(n), ghost_b(n))
      z(1:na) = z_a
      z(na + 1:n) = z_b
      xyz(:, 1:na) = xyz_a
      xyz(:, na + 1:n) = xyz_b
      do k = 1, na
         sym(k) = sym_a(k)
      end do
      do k = 1, nb
         sym(na + k) = sym_b(k)
      end do

      ! A monomer is the dimer with the other half's charges switched off. Note
      ! both masks are over the *same* n atoms in the *same* order.
      ghost_a = .false.
      ghost_a(na + 1:n) = .true.      ! monomer A: B is ghosted
      ghost_b = .false.
      ghost_b(1:na) = .true.          ! monomer B: A is ghosted

      call build_libcint_molecule(z, sym, xyz, basis_name, mols%dimer, error)
      if (error%has_error()) return
      call build_libcint_molecule(z, sym, xyz, basis_name, mols%mono_a, error, &
                                  ghost=ghost_a)
      if (error%has_error()) return
      call build_libcint_molecule(z, sym, xyz, basis_name, mols%mono_b, error, &
                                  ghost=ghost_b)
      if (error%has_error()) return

      if (mols%mono_a%nao /= mols%dimer%nao .or. mols%mono_b%nao /= mols%dimer%nao) then
         call error%set(ERROR_VALIDATION, "sapt: a ghosted monomer does not span "// &
                        "the dimer basis, so ghosting dropped basis functions")
         return
      end if

      mols%n_atoms_a = na
      mols%n_atoms_b = nb
      mols%n_elec_a = sum(z_a)
      mols%n_elec_b = sum(z_b)
      if (present(charge_a)) mols%n_elec_a = mols%n_elec_a - charge_a
      if (present(charge_b)) mols%n_elec_b = mols%n_elec_b - charge_b

      deallocate (z, xyz, sym, ghost_a, ghost_b)
   end subroutine build_sapt_molecules

   subroutine sapt_molecules_destroy(self)
      class(sapt_molecules_t), intent(inout) :: self

      call self%dimer%destroy()
      call self%mono_a%destroy()
      call self%mono_b%destroy()
      self%n_atoms_a = 0
      self%n_atoms_b = 0
   end subroutine sapt_molecules_destroy

end module mqc_sapt
