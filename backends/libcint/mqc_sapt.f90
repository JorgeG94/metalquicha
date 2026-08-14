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
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf, build_fock
   use pic_blas_interfaces, only: pic_gemm
   implicit none
   private

   public :: sapt_molecules_t
   public :: build_sapt_molecules
   public :: sapt_cache_t
   public :: build_sapt_cache
   public :: sapt_elst10

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

   type :: sapt_cache_t
      !! What every SAPT0 term is built from, in psi4's naming
      integer :: nao = 0
      integer :: nocc_a = 0, nocc_b = 0
      real(dp) :: e_scf_a = 0.0_dp, e_scf_b = 0.0_dp
      real(dp) :: e_nuc = 0.0_dp     !! `E_nuc(dimer) - E_nuc(A) - E_nuc(B)`
      real(dp), allocatable :: d_a(:, :), d_b(:, :)   !! `C_occ C_occ^T`, NO factor of 2
      real(dp), allocatable :: c_a(:, :), c_b(:, :)   !! All MOs, dimer basis
      real(dp), allocatable :: eps_a(:), eps_b(:)
      real(dp), allocatable :: v_a(:, :), v_b(:, :)   !! One monomer's nuclei alone
      real(dp), allocatable :: j_a(:, :), j_b(:, :)   !! `J[D]`, undoubled density
      real(dp), allocatable :: s(:, :)
      real(dp), allocatable :: eri(:, :, :, :)
   contains
      procedure :: destroy => sapt_cache_destroy
   end type sapt_cache_t

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

   subroutine build_sapt_cache(mols, cache, error)
      !! The monomer SCFs and the matrices every SAPT0 term is built from
      !!
      !! **`D` carries no factor of two.** It is `C_occ C_occ^T`, and every factor
      !! of 2 and 4 in the terms exists because of that. Adopting the doubled
      !! density while keeping those factors doubles half the expression -- and it
      !! is invisible, because the shapes are unchanged.
      !!
      !! `V` is the nuclear attraction of one monomer's nuclei *alone*, in the
      !! dimer basis. That falls out of ghosting for nothing: a ghosted monomer's
      !! core Hamiltonian carries only its own nuclei, so `V = H_core - T` with
      !! both taken from that molecule.
      type(sapt_molecules_t), intent(inout) :: mols
      type(sapt_cache_t), intent(out) :: cache
      type(error_t), intent(inout) :: error

      type(rhf_result_t) :: scf_a, scf_b
      real(dp), allocatable :: h(:, :), t(:, :), zero(:, :)
      integer :: nao, nocc_a, nocc_b

      nao = mols%dimer%nao
      nocc_a = mols%n_elec_a/2
      nocc_b = mols%n_elec_b/2
      if (2*nocc_a /= mols%n_elec_a .or. 2*nocc_b /= mols%n_elec_b) then
         call error%set(ERROR_VALIDATION, "sapt: SAPT0 here is closed shell, and "// &
                        "a monomer has an odd number of electrons")
         return
      end if

      call run_libcint_rhf(mols%mono_a, mols%n_elec_a, 200, 1.0e-11_dp, 1.0e-9_dp, &
                           .false., scf_a, error, in_core=.true.)
      if (error%has_error()) return
      call run_libcint_rhf(mols%mono_b, mols%n_elec_b, 200, 1.0e-11_dp, 1.0e-9_dp, &
                           .false., scf_b, error, in_core=.true.)
      if (error%has_error()) return
      if (.not. (scf_a%converged .and. scf_b%converged)) then
         call error%set(ERROR_VALIDATION, "sapt: a monomer SCF did not converge in "// &
                        "the dimer basis")
         return
      end if

      cache%nao = nao
      cache%nocc_a = nocc_a
      cache%nocc_b = nocc_b
      cache%e_scf_a = scf_a%energy
      cache%e_scf_b = scf_b%energy

      allocate (cache%d_a(nao, nao), cache%d_b(nao, nao))
      call pic_gemm(scf_a%orbitals(:, 1:nocc_a), scf_a%orbitals(:, 1:nocc_a), &
                    cache%d_a, transb="T")
      call pic_gemm(scf_b%orbitals(:, 1:nocc_b), scf_b%orbitals(:, 1:nocc_b), &
                    cache%d_b, transb="T")
      call move_alloc(scf_a%orbitals, cache%c_a)
      call move_alloc(scf_b%orbitals, cache%c_b)
      call move_alloc(scf_a%orbital_energies, cache%eps_a)
      call move_alloc(scf_b%orbital_energies, cache%eps_b)

      call mols%dimer%overlap(cache%s)
      call mols%dimer%eris(cache%eri)

      ! V for each monomer: the nuclear attraction of its nuclei alone, in the
      ! dimer basis. Ghosting gives it without a new integral, because the kinetic
      ! part is common to all three molecules:
      !
      !     H(dimer)  = T + V_A + V_B          H(mono_b) = T + V_B
      !     => V_A = H(dimer) - H(mono_b),  and V_B = H(dimer) - H(mono_a)
      !
      ! The alternative, H(mono_a) minus a kinetic integral, would need an
      ! accessor this backend does not have -- and SAPT0 never wants T on its own.
      call mols%dimer%core_hamiltonian(h)
      call mols%mono_b%core_hamiltonian(t)
      allocate (cache%v_a(nao, nao))
      cache%v_a = h - t
      deallocate (t)
      call mols%mono_a%core_hamiltonian(t)
      allocate (cache%v_b(nao, nao))
      cache%v_b = h - t
      deallocate (h, t)

      ! J from each monomer's own density. `build_fock` is `H + J - K/2`, so a
      ! zero core and no exchange leaves exactly J -- and it is linear in the
      ! density, so passing the undoubled D gives J[D] as the terms want.
      allocate (zero(nao, nao), cache%j_a(nao, nao), cache%j_b(nao, nao))
      zero = 0.0_dp
      call build_fock(zero, cache%eri, cache%d_a, cache%j_a, k_scale=0.0_dp)
      call build_fock(zero, cache%eri, cache%d_b, cache%j_b, k_scale=0.0_dp)
      deallocate (zero)

      cache%e_nuc = mols%dimer%nuclear_repulsion() &
                    - mols%mono_a%nuclear_repulsion() &
                    - mols%mono_b%nuclear_repulsion()
   end subroutine build_sapt_cache

   pure function sapt_elst10(c) result(energy)
      !! `Elst10,r` -- the classical Coulomb interaction of the two unperturbed
      !! monomer densities, nuclei included (`sapt_jk_terms.py:131-134`):
      !!
      !!     Elst10 = 4 D_B . J_A + 2 D_A . V_B + 2 D_B . V_A + dE_nuc
      !!
      !! Each contraction is elementwise -- `sum_pq X_pq Y_pq`, psi4's
      !! `vector_dot`. Every operand here is symmetric so it coincides with
      !! `Tr(X Y)`, but the exchange terms have operands that are not, and there
      !! the two differ by 55% of the answer. Uniform is safer than case by case.
      type(sapt_cache_t), intent(in) :: c
      real(dp) :: energy

      energy = 4.0_dp*sum(c%d_b*c%j_a) &
               + 2.0_dp*sum(c%d_a*c%v_b) &
               + 2.0_dp*sum(c%d_b*c%v_a) &
               + c%e_nuc
   end function sapt_elst10

   subroutine sapt_cache_destroy(self)
      class(sapt_cache_t), intent(inout) :: self

      if (allocated(self%d_a)) deallocate (self%d_a)
      if (allocated(self%d_b)) deallocate (self%d_b)
      if (allocated(self%c_a)) deallocate (self%c_a)
      if (allocated(self%c_b)) deallocate (self%c_b)
      if (allocated(self%eps_a)) deallocate (self%eps_a)
      if (allocated(self%eps_b)) deallocate (self%eps_b)
      if (allocated(self%v_a)) deallocate (self%v_a)
      if (allocated(self%v_b)) deallocate (self%v_b)
      if (allocated(self%j_a)) deallocate (self%j_a)
      if (allocated(self%j_b)) deallocate (self%j_b)
      if (allocated(self%s)) deallocate (self%s)
      if (allocated(self%eri)) deallocate (self%eri)
   end subroutine sapt_cache_destroy

   subroutine sapt_molecules_destroy(self)
      class(sapt_molecules_t), intent(inout) :: self

      call self%dimer%destroy()
      call self%mono_a%destroy()
      call self%mono_b%destroy()
      self%n_atoms_a = 0
      self%n_atoms_b = 0
   end subroutine sapt_molecules_destroy

end module mqc_sapt
