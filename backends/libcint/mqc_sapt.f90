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
   use mqc_libcint_mp2, only: transform_block
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   implicit none
   private

   public :: sapt_molecules_t
   public :: build_sapt_molecules
   public :: sapt_cache_t
   public :: build_sapt_cache
   public :: sapt_elst10
   public :: sapt_exch10_s2
   public :: sapt_exch10
   public :: sapt_induction
   public :: sapt_disp20
   public :: sapt_exch_disp20
   public :: run_sapt0
   public :: sapt_terms_t

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

   type :: sapt_terms_t
      !! The named pieces of a SAPT0 interaction energy
      real(dp) :: elst10 = 0.0_dp
      real(dp) :: exch10_s2 = 0.0_dp
      real(dp) :: exch10 = 0.0_dp
      real(dp) :: ind20_u = 0.0_dp, ind20_r = 0.0_dp
      real(dp) :: exch_ind20_u = 0.0_dp, exch_ind20_r = 0.0_dp
      real(dp) :: disp20 = 0.0_dp
      real(dp) :: exch_disp20 = 0.0_dp
      real(dp) :: e_int_hf = 0.0_dp   !! Counterpoise-corrected supermolecular HF
      real(dp) :: delta_hf = 0.0_dp
      real(dp) :: total = 0.0_dp
   end type sapt_terms_t

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
      real(dp), allocatable :: k_a(:, :), k_b(:, :)   !! `K[D]`
      real(dp), allocatable :: w_a(:, :), w_b(:, :)   !! `V + 2J`, electrostatic potential
      real(dp), allocatable :: h_a(:, :), h_b(:, :)   !! `V + 2J - K`
      real(dp), allocatable :: k_o(:, :)
         !! `K[D_A S D_B]`, from the inter-monomer transition density.
         !!
         !! **Not symmetric**, and that matters everywhere it is used. psi4 builds
         !! it with the JK pair reversed and transposes in place afterwards
         !! (`sapt_jk_terms.py:96-111`); it is built here in FISAPT's orientation
         !! (`fisapt.cc:3665`) so no transpose is needed and the most bug-prone
         !! line in that file has no counterpart here.
      real(dp), allocatable :: s(:, :)
      real(dp), allocatable :: eri(:, :, :, :)
      real(dp), allocatable :: eri_packed(:, :)   !! For the four-index transforms
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
      real(dp), allocatable :: h(:, :), t(:, :), zero(:, :), tmp(:, :)
      integer :: nao, nocc_a, nocc_b

      nao = mols%dimer%nao
      nocc_a = mols%n_elec_a/2
      nocc_b = mols%n_elec_b/2
      if (2*nocc_a /= mols%n_elec_a .or. 2*nocc_b /= mols%n_elec_b) then
         call error%set(ERROR_VALIDATION, "sapt: SAPT0 here is closed shell, and "// &
                        "a monomer has an odd number of electrons")
         return
      end if

      ! Tight on purpose. A SAPT term is a small difference of large quantities,
      ! and Exch10 in particular tracks the monomer density closely enough that a
      ! 1e-9 density threshold shows up in the ninth decimal of the answer.
      call run_libcint_rhf(mols%mono_a, mols%n_elec_a, 200, 1.0e-12_dp, 1.0e-11_dp, &
                           .false., scf_a, error, in_core=.true.)
      if (error%has_error()) return
      call run_libcint_rhf(mols%mono_b, mols%n_elec_b, 200, 1.0e-12_dp, 1.0e-11_dp, &
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
      call mols%dimer%eris_packed(cache%eri_packed)

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
      allocate (cache%j_a(nao, nao), cache%j_b(nao, nao))
      allocate (cache%k_a(nao, nao), cache%k_b(nao, nao))
      allocate (cache%k_o(nao, nao))
      call coulomb_exchange(cache%eri, cache%d_a, cache%j_a, cache%k_a)
      call coulomb_exchange(cache%eri, cache%d_b, cache%j_b, cache%k_b)

      ! The inter-monomer transition density, and the exchange over it. It is
      ! NOT symmetric -- one index runs through an occupied orbital of A and the
      ! other through one of B -- so K_O and its transpose are different matrices
      ! and both appear in the terms below.
      allocate (zero(nao, nao), tmp(nao, nao))
      call pic_gemm(cache%d_a, cache%s, tmp)
      call pic_gemm(tmp, cache%d_b, zero)          ! zero now holds D_A S D_B
      call coulomb_exchange(cache%eri, zero, tmp, cache%k_o)
      deallocate (zero, tmp)

      allocate (cache%w_a(nao, nao), cache%w_b(nao, nao))
      allocate (cache%h_a(nao, nao), cache%h_b(nao, nao))
      cache%w_a = cache%v_a + 2.0_dp*cache%j_a
      cache%w_b = cache%v_b + 2.0_dp*cache%j_b
      cache%h_a = cache%w_a - cache%k_a
      cache%h_b = cache%w_b - cache%k_b

      cache%e_nuc = mols%dimer%nuclear_repulsion() &
                    - mols%mono_a%nuclear_repulsion() &
                    - mols%mono_b%nuclear_repulsion()
   end subroutine build_sapt_cache

   subroutine coulomb_exchange(eri, d, j, k)
      !! `J[D]` and `K[D]` for an arbitrary, possibly non-symmetric, density
      !!
      !! `build_fock` is `H + J - K/2`, so with a zero core it gives `J` at
      !! `k_scale = 0` and `J - K/2` at 1, and the difference is `K/2`. Two calls
      !! rather than a new integral routine, and correct for a non-symmetric
      !! density because `build_fock` assumes neither symmetry nor antisymmetry.
      real(dp), intent(in) :: eri(:, :, :, :)
      real(dp), intent(in) :: d(:, :)
      real(dp), intent(out) :: j(:, :), k(:, :)

      real(dp), allocatable :: zero(:, :), jk(:, :)
      integer :: n

      n = size(d, 1)
      allocate (zero(n, n), jk(n, n))
      zero = 0.0_dp
      call build_fock(zero, eri, d, j, k_scale=0.0_dp)
      call build_fock(zero, eri, d, jk, k_scale=1.0_dp)
      k = 2.0_dp*(j - jk)
      deallocate (zero, jk)
   end subroutine coulomb_exchange

   function sapt_exch10_s2(c) result(energy)
      !! `Exch10(S^2)`, FISAPT's MCBS form -- `fisapt.cc:2450-2465`
      !!
      !!     -2 D_A.K_B  -2 (D_A S D_B).h_A  -2 (D_B S D_A).h_B
      !!     +2 (D_B S D_A S D_B).w_A  +2 (D_A S D_B S D_A).w_B
      !!     -2 (D_A S D_B).K_O
      !!
      !! Preferred over the DCBS form of `sapt_jk_terms.py:215-222`, which needs
      !! the virtual density `P` and is valid only while `D + P = S^-1`. That
      !! identity fails silently if the monomer SCF drops linearly dependent MOs,
      !! and this form never mentions `P` at all.
      !!
      !! **Every contraction is elementwise**, `sum_pq X_pq Y_pq` -- psi4's
      !! `vector_dot`. The last term has two non-symmetric operands, where that is
      !! *not* `Tr(X Y)`; taking the trace convention there puts the answer 55%
      !! high, which is how it was found in the PySCF reference.
      type(sapt_cache_t), intent(in) :: c
      real(dp) :: energy

      real(dp), allocatable :: ab(:, :), ba(:, :), t(:, :), u(:, :)
      integer :: n

      n = c%nao
      allocate (ab(n, n), ba(n, n), t(n, n), u(n, n))

      call pic_gemm(c%d_a, c%s, t)
      call pic_gemm(t, c%d_b, ab)                ! D_A S D_B
      call pic_gemm(c%d_b, c%s, t)
      call pic_gemm(t, c%d_a, ba)                ! D_B S D_A

      energy = -2.0_dp*sum(c%d_a*c%k_b)
      energy = energy - 2.0_dp*sum(ab*c%h_a)
      energy = energy - 2.0_dp*sum(ba*c%h_b)

      call pic_gemm(ba, c%s, t)
      call pic_gemm(t, c%d_b, u)                 ! D_B S D_A S D_B
      energy = energy + 2.0_dp*sum(u*c%w_a)

      call pic_gemm(ab, c%s, t)
      call pic_gemm(t, c%d_a, u)                 ! D_A S D_B S D_A
      energy = energy + 2.0_dp*sum(u*c%w_b)

      energy = energy - 2.0_dp*sum(ab*c%k_o)

      deallocate (ab, ba, t, u)
   end function sapt_exch10_s2

   function sapt_exch10(c, error) result(energy)
      !! `Exch10` at `S^inf` -- `sapt_jk_terms.py:169-237`
      !!
      !! **This is the one the SAPT0 total uses** (`sapt0.cc:231`).
      !! `Exch10(S^2)` is computed and printed but only forms the ratio that
      !! sSAPT0's exchange scaling needs, so a total built from the `S^2` form is
      !! not comparable to any published SAPT0 number.
      !!
      !! The single-exchange approximation is lifted by inverting the overlap
      !! metric over both monomers' occupied spaces:
      !!
      !!     Sab = [[1, S_AB], [S_AB^T, 1]]     Tmo = Sab^-1 - 1
      !!
      !! psi4 spells the inverse `Matrix::power(-1.0, 1e-14)`, an eigendecomposition
      !! with a *relative* cutoff that zeroes small eigenvalues rather than
      !! inverting them; that is what is done here.
      type(sapt_cache_t), intent(in) :: c
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), parameter :: EIGEN_FLOOR = 1.0e-14_dp
      real(dp), allocatable :: sab(:, :), vals(:), work(:, :)
      real(dp), allocatable :: t_a(:, :), t_b(:, :), t_ab(:, :), tmp(:, :)
      real(dp), allocatable :: jt_a(:, :), kt_a(:, :), jt_ab(:, :), kt_ab(:, :)
      integer :: na, nb, n, no, i, j, info
      real(dp) :: biggest

      energy = 0.0_dp
      na = c%nocc_a
      nb = c%nocc_b
      no = na + nb
      n = c%nao

      ! S_AB over the occupied blocks, then the bordered metric.
      allocate (tmp(n, na), sab(no, no))
      call pic_gemm(c%s, c%c_a(:, 1:na), tmp)
      sab = 0.0_dp
      call pic_gemm(c%c_b(:, 1:nb), tmp, sab(na + 1:no, 1:na), transa="T")
      sab(1:na, na + 1:no) = transpose(sab(na + 1:no, 1:na))
      deallocate (tmp)
      do i = 1, no
         sab(i, i) = sab(i, i) + 1.0_dp
      end do

      ! Invert by eigendecomposition, mirroring psi4's `power(-1.0, cutoff)`.
      allocate (vals(no), work(no, no))
      work = sab
      call pic_syev(work, vals, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "sapt: the occupied overlap metric "// &
                        "could not be diagonalised")
         return
      end if
      biggest = maxval(abs(vals))
      do i = 1, no
         if (abs(vals(i)) < EIGEN_FLOOR*biggest) then
            vals(i) = 0.0_dp
         else
            vals(i) = 1.0_dp/vals(i)
         end if
      end do
      do j = 1, no
         do i = 1, no
            sab(i, j) = sum(work(i, :)*vals(:)*work(j, :))
         end do
      end do
      ! `- 1` on the diagonal. Note the induction path's S^inf branch builds the
      ! same object WITHOUT this subtraction (`sapt_jk_terms.py:387-388` against
      ! `:175-177`) under the same variable names -- do not share a helper.
      do i = 1, no
         sab(i, i) = sab(i, i) - 1.0_dp
      end do
      deallocate (vals)

      ! Back to the AO basis.
      allocate (t_a(n, n), t_b(n, n), t_ab(n, n), tmp(n, no))
      call pic_gemm(c%c_a(:, 1:na), sab(1:na, 1:na), tmp(:, 1:na))
      call pic_gemm(tmp(:, 1:na), c%c_a(:, 1:na), t_a, transb="T")
      call pic_gemm(c%c_b(:, 1:nb), sab(na + 1:no, na + 1:no), tmp(:, 1:nb))
      call pic_gemm(tmp(:, 1:nb), c%c_b(:, 1:nb), t_b, transb="T")
      call pic_gemm(c%c_a(:, 1:na), sab(1:na, na + 1:no), tmp(:, 1:nb))
      call pic_gemm(tmp(:, 1:nb), c%c_b(:, 1:nb), t_ab, transb="T")
      deallocate (tmp, sab, work)

      ! **`JT_AB`/`KT_AB` come from T_AB TRANSPOSED.** psi4 feeds the pair as
      ! (C_left = Cocc_B, C_right = Cocc_A Tmo_AB), so the implied density is
      ! `T_AB^T` (`sapt_jk_terms.py:201-202`), and the `.T` on `KT_AB` in the last
      ! term below puts it back. Building `T_AB` here instead is a silent 55%
      ! error and changes nothing else.
      allocate (jt_a(n, n), kt_a(n, n), jt_ab(n, n), kt_ab(n, n))
      call coulomb_exchange(c%eri, t_a, jt_a, kt_a)
      call coulomb_exchange(c%eri, transpose(t_ab), jt_ab, kt_ab)

      energy = -2.0_dp*sum(c%d_a*c%k_b)
      energy = energy + 2.0_dp*sum(t_a*c%h_b)
      energy = energy + 2.0_dp*sum(t_b*c%h_a)
      energy = energy + 2.0_dp*sum(t_ab*(c%h_a + c%h_b))
      energy = energy + 4.0_dp*sum(t_b*(jt_ab - 0.5_dp*kt_ab))
      energy = energy + 4.0_dp*sum(t_a*(jt_ab - 0.5_dp*kt_ab))
      energy = energy + 4.0_dp*sum(t_b*(jt_a - 0.5_dp*kt_a))
      energy = energy + 4.0_dp*sum(t_ab*(jt_ab - 0.5_dp*transpose(kt_ab)))

      deallocate (t_a, t_b, t_ab, jt_a, kt_a, jt_ab, kt_ab)
   end function sapt_exch10

   subroutine sapt_induction(mols, c, terms, error)
      !! `Ind20` and `Exch-Ind20`, uncoupled and coupled
      !!
      !! **Uncoupled and coupled differ only in how `x` is obtained**; the
      !! contraction after is identical (`sapt_jk_terms.py:350` against `:523`).
      !! So both come out of one code path, and the uncoupled half is computable
      !! with no solver at all -- which is how the sign and scale conventions here
      !! were pinned before `cphf_solve` was involved.
      !!
      !!     w_B_MOA = C_occ_A^T w_B C_vir_A
      !!     uncoupled: x = w / (eps_occ - eps_vir)
      !!     coupled:   x from the response equations
      !!     Ind20      = 2 x . w
      !!     Exch-Ind20 = 2 x . EX
      !!
      !! Only the response form is reported as SAPT0 (`COUPLED_INDUCTION` defaults
      !! true, `read_options.cc:1068`).
      type(sapt_molecules_t), intent(in) :: mols
      type(sapt_cache_t), intent(in) :: c
      type(sapt_terms_t), intent(inout) :: terms
      type(error_t), intent(inout) :: error

      real(dp) :: u_ab, u_ba, r_ab, r_ba, xu_ab, xu_ba, xr_ab, xr_ba

      call one_direction(mols%mono_a, c, "A", "B", u_ab, r_ab, xu_ab, xr_ab, error)
      if (error%has_error()) return
      call one_direction(mols%mono_b, c, "B", "A", u_ba, r_ba, xu_ba, xr_ba, error)
      if (error%has_error()) return

      terms%ind20_u = u_ab + u_ba
      terms%ind20_r = r_ab + r_ba
      terms%exch_ind20_u = xu_ab + xu_ba
      terms%exch_ind20_r = xr_ab + xr_ba
   end subroutine sapt_induction

   subroutine one_direction(mol, c, this, other, ind_u, ind_r, exch_u, exch_r, error)
      !! One monomer polarised by the other's field
      use mqc_libcint_cphf, only: cphf_solve
      type(libcint_molecule_t), intent(in) :: mol
      type(sapt_cache_t), intent(in) :: c
      character(len=1), intent(in) :: this, other
      real(dp), intent(out) :: ind_u, ind_r, exch_u, exch_r
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: coeff(:, :), eps(:), w_ao(:, :), w_mo(:, :)
      real(dp), allocatable :: ex(:, :), x_u(:, :), x_r(:, :), pert(:, :, :)
      real(dp), allocatable :: resp(:, :, :), tmp(:, :)
      integer :: nocc, nvir, nmo, a, r

      if (this == "A") then
         nocc = c%nocc_a
         coeff = c%c_a
         eps = c%eps_a
         w_ao = c%w_b
      else
         nocc = c%nocc_b
         coeff = c%c_b
         eps = c%eps_b
         w_ao = c%w_a
      end if
      nmo = size(coeff, 2)
      nvir = nmo - nocc

      ! The partner's electrostatic potential in this monomer's occupied-virtual
      ! block. `w` and not `h`: no exchange in the perturbing field.
      allocate (w_mo(nocc, nvir), tmp(c%nao, nvir))
      call pic_gemm(w_ao, coeff(:, nocc + 1:nmo), tmp)
      call pic_gemm(coeff(:, 1:nocc), tmp, w_mo, transa="T")
      deallocate (tmp)

      call exch_ind_potential(c, this, other, ex)

      ! Uncoupled: the orbital-energy denominator, negative throughout.
      allocate (x_u(nocc, nvir))
      do r = 1, nvir
         do a = 1, nocc
            x_u(a, r) = w_mo(a, r)/(eps(a) - eps(nocc + r))
         end do
      end do

      ! Coupled: the same contraction, with `x` from the response equations.
      allocate (pert(c%nao, c%nao, 1))
      pert(:, :, 1) = w_ao
      call cphf_solve(mol, coeff, eps, nocc, pert, resp, error, in_core=.true.)
      if (error%has_error()) return
      allocate (x_r(nocc, nvir))
      do r = 1, nvir
         do a = 1, nocc
            x_r(a, r) = resp(r, a, 1)
         end do
      end do

      ind_u = 2.0_dp*sum(x_u*w_mo)
      ind_r = 2.0_dp*sum(x_r*w_mo)
      exch_u = 2.0_dp*sum(x_u*ex)
      exch_r = 2.0_dp*sum(x_r*ex)

      deallocate (coeff, eps, w_ao, w_mo, ex, x_u, x_r, pert, resp)
   end subroutine one_direction

   subroutine exch_ind_potential(c, this, other, ex)
      !! The exchange-induction potential, USAPT0's factorisation
      !!
      !! `usapt0.cc:1261-1315`. Two triplet products where `sapt_jk_terms.py`
      !! spends eighteen matrix chains, and verified equal to them to 4.5e-17.
      !!
      !!     W  = -K_o - 2 J_O + K_O + 2 J[D_o S D_t S D_o]
      !!     T1 = -h_t + S D_t w_o + w_t D_o S - K_O^T
      !!     W += S D_o T1
      !!     T2 = -h_o + w_o D_t S - K_O
      !!     W += T2 D_o S
      !!     EX = C_occ_t^T W C_vir_t
      !!
      !! **`T1` contains `S D_t w_o` -- `D_t`, this monomer's density, not the
      !! partner's.** The outer `S D_o` supplies the other projector, and getting
      !! it wrong is dimensionally invisible.
      !!
      !! One routine serves both directions: called the other way round the only
      !! further change is `K_O -> K_O^T`, which is why it is written this way
      !! rather than twice.
      type(sapt_cache_t), intent(in) :: c
      character(len=1), intent(in) :: this, other
      real(dp), allocatable, intent(out) :: ex(:, :)

      real(dp), allocatable :: d_t(:, :), d_o(:, :), h_t(:, :), h_o(:, :)
      real(dp), allocatable :: w_t(:, :), w_o(:, :), k_o_mat(:, :), k_tr(:, :)
      real(dp), allocatable :: w(:, :), t1(:, :), t2(:, :), p(:, :), q(:, :)
      real(dp), allocatable :: j_o(:, :), k_junk(:, :), tmp(:, :)
      integer :: n, nocc, nmo

      n = c%nao
      if (this == "A") then
         d_t = c%d_a
         d_o = c%d_b
         h_t = c%h_a
         h_o = c%h_b
         w_t = c%w_a
         w_o = c%w_b
         k_o_mat = c%k_b
         k_tr = c%k_o
         nocc = c%nocc_a
         nmo = size(c%c_a, 2)
      else
         d_t = c%d_b
         d_o = c%d_a
         h_t = c%h_b
         h_o = c%h_a
         w_t = c%w_b
         w_o = c%w_a
         k_o_mat = c%k_a
         k_tr = transpose(c%k_o)        ! the one asymmetry between the two calls
         nocc = c%nocc_b
         nmo = size(c%c_b, 2)
      end if

      allocate (w(n, n), t1(n, n), t2(n, n), p(n, n), q(n, n), tmp(n, n))
      allocate (j_o(n, n), k_junk(n, n))

      ! J_O, from the transition density in *this* call's orientation.
      call pic_gemm(d_t, c%s, tmp)
      call pic_gemm(tmp, d_o, p)                      ! D_t S D_o
      call coulomb_exchange(c%eri, p, j_o, k_junk)

      ! J of D_o S D_t S D_o
      call pic_gemm(d_o, c%s, tmp)
      call pic_gemm(tmp, d_t, q)
      call pic_gemm(q, c%s, tmp)
      call pic_gemm(tmp, d_o, p)
      call coulomb_exchange(c%eri, p, q, k_junk)      ! q now holds J[D_o S D_t S D_o]

      w = -k_o_mat - 2.0_dp*j_o + k_tr + 2.0_dp*q

      call pic_gemm(c%s, d_t, tmp)
      call pic_gemm(tmp, w_o, t1)                     ! S D_t w_o
      call pic_gemm(w_t, d_o, tmp)
      call pic_gemm(tmp, c%s, q)                      ! w_t D_o S
      t1 = -h_t + t1 + q - transpose(k_tr)

      call pic_gemm(c%s, d_o, tmp)
      call pic_gemm(tmp, t1, q)
      w = w + q

      call pic_gemm(w_o, d_t, tmp)
      call pic_gemm(tmp, c%s, t2)                     ! w_o D_t S
      t2 = -h_o + t2 - k_tr

      call pic_gemm(d_o, c%s, tmp)
      call pic_gemm(t2, tmp, q)
      w = w + q

      allocate (ex(nocc, nmo - nocc))
      deallocate (tmp)
      allocate (tmp(n, nmo - nocc))
      if (this == "A") then
         call pic_gemm(w, c%c_a(:, nocc + 1:nmo), tmp)
         call pic_gemm(c%c_a(:, 1:nocc), tmp, ex, transa="T")
      else
         call pic_gemm(w, c%c_b(:, nocc + 1:nmo), tmp)
         call pic_gemm(c%c_b(:, 1:nocc), tmp, ex, transa="T")
      end if

      deallocate (d_t, d_o, h_t, h_o, w_t, w_o, k_o_mat, k_tr)
      deallocate (w, t1, t2, p, q, tmp, j_o, k_junk)
   end subroutine exch_ind_potential

   function sapt_disp20(c) result(energy)
      !! `Disp20` -- `exch-disp20.cc:184-195`
      !!
      !!     Disp20 = 4 sum_{a,r in A} sum_{b,s in B}
      !!              (ar|bs)^2 / (eps_a + eps_b - eps_r - eps_s)
      !!
      !! Chemists' notation, no antisymmetrisation: the bra pair is (occ A, vir A)
      !! and the ket pair (occ B, vir B). Pairing `(ab)` with `(rs)` instead gives
      !! a plausible and entirely wrong number.
      !!
      !! **Not the routine in `disp20.cc`**, which is debug-only Laplace code
      !! called under `if (debug_)`, and whose neighbouring block computes a
      !! natural-orbital variant assigning a different variable.
      !!
      !! The denominator is negative, so `Disp20` is negative. Factor 4, not 2 --
      !! the closed-shell spin sum over both monomers. Its same-spin and
      !! opposite-spin halves are exactly equal by construction, so a decomposition
      !! that is not 50/50 is a bug rather than physics.
      type(sapt_cache_t), intent(in) :: c
      real(dp) :: energy

      real(dp), allocatable :: v(:, :, :, :)
      integer :: na, ra, nb, rb, a, r, b, s
      integer :: nmo_a, nmo_b
      real(dp) :: denom

      na = c%nocc_a
      nb = c%nocc_b
      nmo_a = size(c%c_a, 2)
      nmo_b = size(c%c_b, 2)
      ra = nmo_a - na
      rb = nmo_b - nb

      call transform_block(c%eri_packed, c%c_a(:, 1:na), c%c_a(:, na + 1:nmo_a), &
                           c%c_b(:, 1:nb), c%c_b(:, nb + 1:nmo_b), v)

      energy = 0.0_dp
      do s = 1, rb
         do b = 1, nb
            do r = 1, ra
               do a = 1, na
                  denom = c%eps_a(a) + c%eps_b(b) &
                          - c%eps_a(na + r) - c%eps_b(nb + s)
                  energy = energy + 4.0_dp*v(a, r, b, s)*v(a, r, b, s)/denom
               end do
            end do
         end do
      end do
      deallocate (v)
   end function sapt_disp20

   function sapt_exch_disp20(c) result(energy)
      !! `Exch-Disp20`, the S^2 form -- FISAPT's, `fisapt.cc:4351-4733`
      !!
      !! `libsapt_solver`'s version is factorised into some twenty `h`/`q` pieces
      !! and is not term-comparable; FISAPT's is the readable one and the two
      !! agree to 2.6e-8, that being libsapt's Laplace denominators.
      !!
      !! Four extra `(occ,vir)` transforms per monomer with modified coefficient
      !! sets, plus four rank-one terms built from `J`, `K` and `V`. The shapes
      !! are Disp20's throughout, so there is no new bottleneck -- about five
      !! times its integral cost and none of its scaling.
      !!
      !! **Not optional if a total is to be reported.** It runs 10-25% of `Disp20`
      !! across psi4's own tests and always with the opposite sign; omitting it
      !! makes a total systematically and plausibly wrong.
      type(sapt_cache_t), intent(in) :: c
      real(dp) :: energy

      real(dp), allocatable :: cr1(:, :), cs1(:, :), ca2(:, :), cb2(:, :)
      real(dp), allocatable :: cr3(:, :), cs3(:, :), ca4(:, :), cb4(:, :)
      real(dp), allocatable :: oa(:, :), va(:, :), ob(:, :), vb(:, :)
      real(dp), allocatable :: g1(:, :, :, :), g2(:, :, :, :), g3(:, :, :, :)
      real(dp), allocatable :: g4(:, :, :, :), g5(:, :, :, :), g6(:, :, :, :)
      real(dp), allocatable :: amp(:, :, :, :), v(:, :, :, :)
      real(dp), allocatable :: sas(:, :), sbr(:, :), sbar(:, :), sabs(:, :)
      real(dp), allocatable :: qbr(:, :), qas(:, :), qar(:, :), qbs(:, :)
      real(dp), allocatable :: t(:, :), u(:, :), eye(:, :)
      integer :: n, na, ra, nb, rb, nmo_a, nmo_b, a, r, b, s, i
      real(dp) :: denom

      n = c%nao
      na = c%nocc_a
      nb = c%nocc_b
      nmo_a = size(c%c_a, 2)
      nmo_b = size(c%c_b, 2)
      ra = nmo_a - na
      rb = nmo_b - nb

      oa = c%c_a(:, 1:na)
      va = c%c_a(:, na + 1:nmo_a)
      ob = c%c_b(:, 1:nb)
      vb = c%c_b(:, nb + 1:nmo_b)

      allocate (eye(n, n), t(n, n), u(n, n))
      eye = 0.0_dp
      do i = 1, n
         eye(i, i) = 1.0_dp
      end do

      ! The modified coefficient sets, fisapt.cc:4351-4370.
      call pic_gemm(c%d_b, c%s, t)
      cr1 = matmul(eye - t, va)
      ca2 = matmul(t, oa)
      call pic_gemm(c%d_a, c%s, u)
      cs1 = matmul(eye - u, vb)
      cb2 = matmul(u, ob)
      cr3 = 2.0_dp*matmul(t - matmul(u, t), va)
      cs3 = 2.0_dp*matmul(u - matmul(t, u), vb)
      ca4 = -2.0_dp*matmul(matmul(u, t), oa)
      cb4 = -2.0_dp*matmul(matmul(t, u), ob)

      ! v[a,r,b,s], assembled in FISAPT's order (fisapt.cc:4706-4716).
      call transform_block(c%eri_packed, ob, cr1, oa, cs1, g1)
      call transform_block(c%eri_packed, cb2, va, ca2, vb, g2)
      call transform_block(c%eri_packed, oa, va, ob, cs3, g3)
      call transform_block(c%eri_packed, oa, va, cb4, vb, g4)
      call transform_block(c%eri_packed, oa, cr3, ob, vb, g5)
      call transform_block(c%eri_packed, ca4, va, ob, vb, g6)

      allocate (v(na, ra, nb, rb))
      do s = 1, rb
         do b = 1, nb
            do r = 1, ra
               do a = 1, na
                  v(a, r, b, s) = g1(b, r, a, s) + g2(b, r, a, s) &
                                  + g3(a, r, b, s) + g4(a, r, b, s) &
                                  + g5(a, r, b, s) + g6(a, r, b, s)
               end do
            end do
         end do
      end do
      deallocate (g1, g2, g3, g4, g5, g6)

      ! Orbital-space overlap blocks and the Q matrices, fisapt.cc:4374-4452.
      sas = matmul(transpose(oa), matmul(c%s, vb))
      sbr = matmul(transpose(ob), matmul(c%s, va))
      sbar = matmul(transpose(oa), matmul(c%s, matmul(c%d_b, matmul(c%s, va))))
      sabs = matmul(transpose(ob), matmul(c%s, matmul(c%d_a, matmul(c%s, vb))))

      qbr = 2.0_dp*matmul(transpose(ob), matmul(c%j_a, va)) &
            - matmul(transpose(ob), matmul(c%k_a, va)) &
            + matmul(transpose(ob), matmul(transpose(c%k_o), va)) &
            - 2.0_dp*matmul(transpose(ob), matmul(c%s, matmul(c%d_a, matmul(c%j_b, va)))) &
            - 2.0_dp*matmul(transpose(ob), matmul(c%j_a, matmul(c%d_b, matmul(c%s, va)))) &
            - matmul(transpose(ob), matmul(c%s, matmul(c%d_a, matmul(c%v_b, va)))) &
            + matmul(transpose(ob), matmul(c%v_a, matmul(p_of(c, "B"), matmul(c%s, va))))
      qas = 2.0_dp*matmul(transpose(oa), matmul(c%j_b, vb)) &
            - matmul(transpose(oa), matmul(c%k_b, vb)) &
            + matmul(transpose(oa), matmul(c%k_o, vb)) &
            - 2.0_dp*matmul(transpose(oa), matmul(c%j_b, matmul(c%d_a, matmul(c%s, vb)))) &
            - 2.0_dp*matmul(transpose(oa), matmul(c%s, matmul(c%d_b, matmul(c%j_a, vb)))) &
            - matmul(transpose(oa), matmul(c%s, matmul(c%d_b, matmul(c%v_a, vb)))) &
            + matmul(transpose(oa), matmul(c%v_b, matmul(p_of(c, "A"), matmul(c%s, vb))))
      qar = 4.0_dp*matmul(transpose(oa), matmul(c%j_b, va)) &
            + 2.0_dp*matmul(transpose(oa), matmul(c%v_b, va))
      qbs = 4.0_dp*matmul(transpose(ob), matmul(c%j_a, vb)) &
            + 2.0_dp*matmul(transpose(ob), matmul(c%v_a, vb))

      do s = 1, rb
         do b = 1, nb
            do r = 1, ra
               do a = 1, na
                  v(a, r, b, s) = v(a, r, b, s) &
                                  + qbr(b, r)*sas(a, s) + sbr(b, r)*qas(a, s) &
                                  + qar(a, r)*sabs(b, s) + sbar(a, r)*qbs(b, s)
               end do
            end do
         end do
      end do

      ! `t` is the plain dispersion amplitude, reused (fisapt.cc:4697).
      call transform_block(c%eri_packed, oa, va, ob, vb, amp)
      energy = 0.0_dp
      do s = 1, rb
         do b = 1, nb
            do r = 1, ra
               do a = 1, na
                  denom = c%eps_a(a) + c%eps_b(b) &
                          - c%eps_a(na + r) - c%eps_b(nb + s)
                  energy = energy - 2.0_dp*(amp(a, r, b, s)/denom)*v(a, r, b, s)
               end do
            end do
         end do
      end do

      deallocate (v, amp, eye, t, u)
   end function sapt_exch_disp20

   function p_of(c, which) result(p)
      !! The virtual density of one monomer, `C_vir C_vir^T`
      type(sapt_cache_t), intent(in) :: c
      character(len=1), intent(in) :: which
      real(dp), allocatable :: p(:, :)

      integer :: nocc, nmo

      allocate (p(c%nao, c%nao))
      if (which == "A") then
         nocc = c%nocc_a
         nmo = size(c%c_a, 2)
         call pic_gemm(c%c_a(:, nocc + 1:nmo), c%c_a(:, nocc + 1:nmo), p, transb="T")
      else
         nocc = c%nocc_b
         nmo = size(c%c_b, 2)
         call pic_gemm(c%c_b(:, nocc + 1:nmo), c%c_b(:, nocc + 1:nmo), p, transb="T")
      end if
   end function p_of

   subroutine run_sapt0(mols, terms, error)
      !! Every SAPT0 term, and the total
      !!
      !!     E(SAPT0) = Elst10,r + Exch10 + Ind20,r + Exch-Ind20,r + dHF(2)
      !!                + Disp20 + Exch-Disp20
      !!
      !! from `sapt0.cc:229-231`. Three things in that are not guessable and each
      !! changes the answer:
      !!
      !! * **`dHF(2)` is in the total, unconditionally.** There is no option to
      !!   disable it, so every SAPT0 number in the literature includes it. On a
      !!   water dimer it is worth around a sixth of the total.
      !! * **The total uses `Exch10`, not `Exch10(S^2)`.** The `S^2` form is
      !!   computed and printed but only forms the ratio sSAPT0's exchange scaling
      !!   needs, so a total built from it is comparable to nothing.
      !! * **`Ind20,r`, the response form** -- `COUPLED_INDUCTION` defaults true.
      !!   But `Exch-Ind20` and `Exch-Disp20` are always the `S^2` approximation
      !!   in SAPT0; their `S^inf` variants are not implemented in psi4 at all.
      !!
      !! `dHF(2)` is *defined* as the residual, so the first four terms plus it
      !! are the counterpoise-corrected supermolecular Hartree-Fock interaction
      !! energy by construction. That identity is the sharpest check available on
      !! those four collectively, and it needs no SAPT reference -- see
      !! `test_mqc_sapt`.
      type(sapt_molecules_t), intent(inout) :: mols
      type(sapt_terms_t), intent(out) :: terms
      type(error_t), intent(inout) :: error

      type(sapt_cache_t) :: c
      type(rhf_result_t) :: scf_d
      real(dp) :: four

      call build_sapt_cache(mols, c, error)
      if (error%has_error()) return

      terms%elst10 = sapt_elst10(c)
      terms%exch10_s2 = sapt_exch10_s2(c)
      terms%exch10 = sapt_exch10(c, error)
      if (error%has_error()) return
      call sapt_induction(mols, c, terms, error)
      if (error%has_error()) return
      terms%disp20 = sapt_disp20(c)
      terms%exch_disp20 = sapt_exch_disp20(c)

      ! The supermolecular reference. Counterpoise-corrected for free: the
      ! monomer energies already came from SCFs in the dimer basis.
      call run_libcint_rhf(mols%dimer, mols%n_elec_a + mols%n_elec_b, 200, &
                           1.0e-12_dp, 1.0e-11_dp, .false., scf_d, error, &
                           in_core=.true.)
      if (error%has_error()) return
      if (.not. scf_d%converged) then
         call error%set(ERROR_VALIDATION, "sapt: the dimer SCF did not converge")
         return
      end if
      terms%e_int_hf = scf_d%energy - c%e_scf_a - c%e_scf_b

      four = terms%elst10 + terms%exch10 + terms%ind20_r + terms%exch_ind20_r
      terms%delta_hf = terms%e_int_hf - four
      terms%total = four + terms%delta_hf + terms%disp20 + terms%exch_disp20

      call c%destroy()
   end subroutine run_sapt0

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
      if (allocated(self%k_a)) deallocate (self%k_a)
      if (allocated(self%k_b)) deallocate (self%k_b)
      if (allocated(self%w_a)) deallocate (self%w_a)
      if (allocated(self%w_b)) deallocate (self%w_b)
      if (allocated(self%h_a)) deallocate (self%h_a)
      if (allocated(self%h_b)) deallocate (self%h_b)
      if (allocated(self%k_o)) deallocate (self%k_o)
      if (allocated(self%s)) deallocate (self%s)
      if (allocated(self%eri)) deallocate (self%eri)
      if (allocated(self%eri_packed)) deallocate (self%eri_packed)
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
