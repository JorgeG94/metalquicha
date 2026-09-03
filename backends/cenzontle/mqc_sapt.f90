!! Symmetry-adapted perturbation theory, SAPT0 and SAPT2
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
   !! contracted against another is silently wrong. `Tr(D_B S)` has to be the
   !! occupied count, and is asserted here for that reason.
   !!
   !! Reference values come from `validation/check_sapt0.py`, a conventional
   !! four-index SAPT0 in PySCF. psi4 cannot serve as the reference: its
   !! closed-shell SAPT module is density-fitted by construction and, in its own
   !! documentation, "cannot be run with exact integrals".
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, &
                                pair_index
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, build_fock
   use mqc_czt_mp2, only: transform_block
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
   public :: sapt2_cache_t
   public :: build_sapt2_cache
   public :: sapt2_amps_t
   public :: build_sapt2_amps
   public :: sapt2_zero_amps
   public :: sapt2_amps_mp2_energy
   public :: sapt2_k2f
   public :: sapt_elst12
   public :: sapt_exch11
   public :: sapt_exch12
   public :: sapt_ind22
   public :: run_sapt2

   type :: sapt_molecules_t
      !! The dimer and its two counterpoise-corrected monomers
      type(czt_molecule_t) :: dimer
      type(czt_molecule_t) :: mono_a   !! A real, B ghosted
      type(czt_molecule_t) :: mono_b   !! B real, A ghosted
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
      ! -- the SAPT2 terms, filled by `run_sapt2` only --
      real(dp) :: elst12 = 0.0_dp
      real(dp) :: exch11 = 0.0_dp
      real(dp) :: exch12 = 0.0_dp
      real(dp) :: ind22 = 0.0_dp
      real(dp) :: exch_ind22 = 0.0_dp
         !! Not computed directly: `Ind22` scaled by `Exch-Ind20,r / Ind20,r`
         !! (`ind22.cc:52`). Zero when `Ind20,r` is too small for the ratio
         !! to mean anything -- psi4 does not guard that division; we do.
      real(dp) :: total_sapt2 = 0.0_dp
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
         !! **Not symmetric**, and that matters everywhere it is used. Built in
         !! FISAPT's orientation (`fisapt.cc:3665`), not psi4's reversed one
         !! (`sapt_jk_terms.py:96-111`), so no transpose is needed.
      real(dp), allocatable :: s(:, :)
      real(dp), allocatable :: eri(:, :, :, :)
         !! Dense `(nao, nao, nao, nao)`, for the `J`/`K` builds
      real(dp), allocatable :: eri_packed(:, :)   !! For the four-index transforms
   contains
      procedure :: destroy => sapt_cache_destroy
   end type sapt_cache_t

   type :: sapt2_cache_t
      !! What the four SAPT2 terms are built from, mirroring psi4's SAPT2 state
      !!
      !! The three-index blocks `b_*` are MO blocks of an **exact**
      !! factorization of the conventional ERIs -- an eigendecomposition of the
      !! packed AO ERI matrix stands in for psi4's fitting basis, so
      !! `sum_P b_xy^P b_zw^P = (xy|zw)` to machine precision and every dressed
      !! contraction of `libsapt_solver` can be transcribed as written while
      !! the result stays conventional.
      !!
      !! Each block carries `ndf + 3` columns. The last three are psi4's
      !! "dressing" (`utils.cc`, `get_*_ints`), which folds the one-electron
      !! and nuclear parts of the intermolecular operator into the same
      !! contractions: an A-monomer pair dresses as
      !! `(delta, vB/NB, e*delta)`, a B-monomer pair as
      !! `(vA/NA, delta, e*delta)`, and a mixed pair carries `sAB` wherever a
      !! same-monomer pair carries `delta`, with `e = sqrt(enuc/(NA*NB))`.
      !! The raw blocks stored here have those columns zero; the `dressed_*`
      !! helpers fill them per psi4's dress flags.
      integer :: ndf = 0             !! Rank of the exact factorization
      integer :: nd3 = 0            !! `ndf + 3`, the stored column count
      integer :: nocc_a = 0, nvir_a = 0, nmo_a = 0
      integer :: nocc_b = 0, nvir_b = 0, nmo_b = 0
      integer :: nea = 0, neb = 0   !! Electron counts, psi4's `NA_`/`NB_`
      real(dp) :: enuc = 0.0_dp     !! Intermolecular nuclear repulsion
      real(dp) :: esq = 0.0_dp      !! `sqrt(enuc/(NA*NB))`
      real(dp), allocatable :: eps_a(:), eps_b(:)
      real(dp), allocatable :: s_ab(:, :)    !! `CA^T S CB`, full MO ranges
      real(dp), allocatable :: v_baa(:, :)   !! `CA^T V_B CA`
      real(dp), allocatable :: v_abb(:, :)   !! `CB^T V_A CB`
      real(dp), allocatable :: v_aab(:, :)   !! `CA^T V_A CB`
      real(dp), allocatable :: v_bab(:, :)   !! `CA^T V_B CB`
      real(dp), allocatable :: w_baa(:, :), w_bar(:, :), w_brr(:, :)
         !! `w_B = V_B + 2 J_B` in monomer A's occ-occ, occ-vir, vir-vir blocks
      real(dp), allocatable :: w_abb(:, :), w_abs(:, :), w_ass(:, :)
      real(dp), allocatable :: b_aa(:, :, :), b_ar(:, :, :), b_rr(:, :, :)
      real(dp), allocatable :: b_bb(:, :, :), b_bs(:, :, :), b_ss(:, :, :)
      real(dp), allocatable :: b_ab(:, :, :), b_as(:, :, :), b_rb(:, :, :)
      real(dp), allocatable :: diag_aa(:), diag_bb(:)
         !! `sum_a b_aa(a,a,P)` with the dressed tail, `sapt2.cc:815-833`
   contains
      procedure :: destroy => sapt2_cache_destroy
   end type sapt2_cache_t

   type :: sapt2_amps_t
      !! One monomer's MP2-level amplitudes, psi4's `amplitudes.cc` labels
      real(dp), allocatable :: t(:, :, :, :)
         !! `tARAR`: `(ar|a'r') / (e_a + e_a' - e_r - e_r')`. The denominator
         !! is included and nothing is antisymmetrized -- getting either wrong
         !! is a factor of 2 or a sign in every term at once.
      real(dp), allocatable :: theta(:, :, :, :)   !! `2 t - t(occ swapped)`
      real(dp), allocatable :: th_ov(:, :, :)
         !! "Theta AR Intermediates": `theta` against the dressed occ-vir block
      real(dp), allocatable :: p_oo(:, :), p_vv(:, :)
         !! The unrelaxed MP2 density blocks, WITHOUT their -2/+2 -- those
         !! factors live in the consuming terms, as in psi4
      real(dp), allocatable :: y2(:, :)       !! "Y2 AR Amplitudes"
      real(dp), allocatable :: t_singles(:, :)
         !! "T2 AR Amplitudes": Y2's theta part over `(e_occ - e_vir)`,
         !! written BEFORE the density parts are added (`amplitudes.cc:216`)
      real(dp), allocatable :: t2(:, :, :, :)  !! The second-order doubles
      real(dp), allocatable :: th2_ov(:, :, :)  !! "Theta 2": theta() on `t2`
   contains
      procedure :: destroy => sapt2_amps_destroy
   end type sapt2_amps_t

   integer, parameter :: SAPT2_PIECES = 7
      !! Both `Ind22` and the `Exch12` K2f energies come in seven pieces,
      !! following psi4's `ind220_1..7` and `e1..e7` decompositions

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

      call build_czt_molecule(z, sym, xyz, basis_name, mols%dimer, error)
      if (error%has_error()) return
      call build_czt_molecule(z, sym, xyz, basis_name, mols%mono_a, error, &
                              ghost=ghost_a)
      if (error%has_error()) return
      call build_czt_molecule(z, sym, xyz, basis_name, mols%mono_b, error, &
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
      !! density while keeping those factors doubles half the expression, and the
      !! shapes are unchanged so nothing objects.
      !!
      !! `V` is the nuclear attraction of one monomer's nuclei *alone*, in the
      !! dimer basis, which ghosting supplies with no new integral.
      ! TODO(mqc): fills both `eri` and `eri_packed` and holds them for the whole
      ! run, so the AO integrals are stored twice over in two different layouts
      ! -- the dense copy is four times the packed one.
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
         call error%set(ERROR_VALIDATION, "sapt: the monomer references here are "// &
                        "RHF, and a monomer has an odd number of electrons. If "// &
                        "the monomer is an ion, give it its charge -- an "// &
                        "unstated one is counted as neutral.")
         return
      end if

      ! Tight on purpose: a SAPT term is a small difference of large quantities,
      ! and a 1e-9 density threshold shows up in Exch10's ninth decimal.
      call run_czt_rhf(mols%mono_a, mols%n_elec_a, 200, 1.0e-12_dp, 1.0e-11_dp, &
                       .false., scf_a, error, in_core=.true.)
      if (error%has_error()) return
      call run_czt_rhf(mols%mono_b, mols%n_elec_b, 200, 1.0e-12_dp, 1.0e-11_dp, &
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

      ! V for each monomer, using that the kinetic part is common to all three
      ! molecules:
      !
      !     H(dimer)  = T + V_A + V_B          H(mono_b) = T + V_B
      !     => V_A = H(dimer) - H(mono_b),  and V_B = H(dimer) - H(mono_a)
      call mols%dimer%core_hamiltonian(h)
      call mols%mono_b%core_hamiltonian(t)
      allocate (cache%v_a(nao, nao))
      cache%v_a = h - t
      deallocate (t)
      call mols%mono_a%core_hamiltonian(t)
      allocate (cache%v_b(nao, nao))
      cache%v_b = h - t
      deallocate (h, t)

      ! J from each monomer's own density. `build_fock` is linear in the density,
      ! so passing the undoubled D gives J[D] as the terms want.
      allocate (cache%j_a(nao, nao), cache%j_b(nao, nao))
      allocate (cache%k_a(nao, nao), cache%k_b(nao, nao))
      allocate (cache%k_o(nao, nao))
      call coulomb_exchange(cache%eri, cache%d_a, cache%j_a, cache%k_a)
      call coulomb_exchange(cache%eri, cache%d_b, cache%j_b, cache%k_b)

      ! The inter-monomer transition density, and the exchange over it. NOT
      ! symmetric, so `K_O` and its transpose are different matrices and both
      ! appear in the terms below.
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
      !! `k_scale = 0` and `J - K/2` at 1, and the difference is `K/2`. Correct
      !! for a non-symmetric density, because `build_fock` assumes neither
      !! symmetry nor antisymmetry.
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
      !! the virtual density `P` and is valid only while `D + P = S^-1` -- an
      !! identity that fails silently if the monomer SCF drops linearly dependent
      !! MOs.
      !!
      !! **Every contraction is elementwise**, `sum_pq X_pq Y_pq` -- psi4's
      !! `vector_dot`. The last term has two non-symmetric operands, where that
      !! is *not* `Tr(X Y)`.
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
      !! **This is the one the SAPT0 total uses** (`sapt0.cc:231`), not
      !! `Exch10(S^2)`, which only forms the ratio sSAPT0's exchange scaling
      !! needs.
      !!
      !! The single-exchange approximation is lifted by inverting the overlap
      !! metric over both monomers' occupied spaces:
      !!
      !!     Sab = [[1, S_AB], [S_AB^T, 1]]     Tmo = Sab^-1 - 1
      !!
      !! The inverse is psi4's `Matrix::power(-1.0, 1e-14)`: an
      !! eigendecomposition with a *relative* cutoff that zeroes small
      !! eigenvalues rather than inverting them.
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
      ! `- 1` on the diagonal. The induction path's S^inf branch builds the same
      ! object WITHOUT this subtraction (`sapt_jk_terms.py:387-388` against
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
      ! `T_AB^T` (`sapt_jk_terms.py:201-202`), and the transpose on `KT_AB` in
      ! the last term below puts it back. Building `T_AB` here instead is a
      ! silent error and changes nothing else.
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

   subroutine sapt_induction(mols, c, terms, error, chf_a, chf_b)
      !! `Ind20` and `Exch-Ind20`, uncoupled and coupled
      !!
      !! **Uncoupled and coupled differ only in how `x` is obtained**; the
      !! contraction after is identical (`sapt_jk_terms.py:350` against `:523`),
      !! so both come out of one code path.
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
      real(dp), allocatable, intent(out), optional :: chf_a(:, :), chf_b(:, :)
         !! The CPHF coefficients behind `Ind20,r`, `(occ, vir)`. SAPT2's
         !! `Elst12` reuses them so that no Z-vector is ever solved
         !! (`elst12.cc:83`, the `CHFA` contraction).

      real(dp) :: u_ab, u_ba, r_ab, r_ba, xu_ab, xu_ba, xr_ab, xr_ba

      call one_direction(mols%mono_a, c, "A", "B", u_ab, r_ab, xu_ab, xr_ab, &
                         error, chf=chf_a)
      if (error%has_error()) return
      call one_direction(mols%mono_b, c, "B", "A", u_ba, r_ba, xu_ba, xr_ba, &
                         error, chf=chf_b)
      if (error%has_error()) return

      terms%ind20_u = u_ab + u_ba
      terms%ind20_r = r_ab + r_ba
      terms%exch_ind20_u = xu_ab + xu_ba
      terms%exch_ind20_r = xr_ab + xr_ba
   end subroutine sapt_induction

   subroutine one_direction(mol, c, this, other, ind_u, ind_r, exch_u, exch_r, &
                            error, chf)
      !! One monomer polarised by the other's field
      use mqc_czt_cphf, only: cphf_solve
      type(czt_molecule_t), intent(in) :: mol
      type(sapt_cache_t), intent(in) :: c
      character(len=1), intent(in) :: this, other
      real(dp), intent(out) :: ind_u, ind_r, exch_u, exch_r
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: chf(:, :)

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

      if (present(chf)) chf = x_r

      deallocate (coeff, eps, w_ao, w_mo, ex, x_u, x_r, pert, resp)
   end subroutine one_direction

   subroutine exch_ind_potential(c, this, other, ex)
      !! The exchange-induction potential, USAPT0's factorisation
      !!
      !! `usapt0.cc:1261-1315`, which is two triplet products where
      !! `sapt_jk_terms.py` spends eighteen matrix chains.
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
      !! One routine serves both directions; the only difference between them is
      !! `K_O -> K_O^T`.
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
      !! **Not the routine in `disp20.cc`**, which is debug-only Laplace code.
      !!
      !! The denominator is negative, so `Disp20` is negative. Factor 4, not 2 --
      !! the closed-shell spin sum over both monomers. Its same-spin and
      !! opposite-spin halves are exactly equal by construction, so a
      !! decomposition that is not 50/50 is a bug rather than physics.
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
      !! and is not term-comparable; FISAPT's is the readable one, and the two
      !! differ only by libsapt's Laplace denominators.
      !!
      !! Four extra `(occ,vir)` transforms per monomer with modified coefficient
      !! sets, plus four rank-one terms built from `J`, `K` and `V`. The shapes
      !! are Disp20's throughout, so there is no new bottleneck.
      !!
      !! **Not optional if a total is to be reported.** It runs at 10-25% of
      !! `Disp20` and always with the opposite sign.
      ! TODO(mqc): built out of chained `matmul` where every other routine in
      ! this module goes through `pic_gemm`, so these `nao` by `nao` chains miss
      ! BLAS entirely and allocate a temporary per link.
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
      !! * **`dHF(2)` is in the total, unconditionally**, so every SAPT0 number
      !!   in the literature includes it.
      !! * **The total uses `Exch10`, not `Exch10(S^2)`.** A total built from the
      !!   `S^2` form is comparable to nothing.
      !! * **`Ind20,r`, the response form** -- `COUPLED_INDUCTION` defaults true.
      !!   But `Exch-Ind20` and `Exch-Disp20` are always the `S^2` approximation
      !!   in SAPT0; their `S^inf` variants are not in psi4 at all.
      !!
      !! `dHF(2)` is *defined* as the residual, so the first four terms plus it
      !! are the counterpoise-corrected supermolecular Hartree-Fock interaction
      !! energy by construction -- an identity `test_mqc_sapt` checks with no
      !! SAPT reference.
      type(sapt_molecules_t), intent(inout) :: mols
      type(sapt_terms_t), intent(out) :: terms
      type(error_t), intent(inout) :: error

      type(sapt_cache_t) :: c

      call build_sapt_cache(mols, c, error)
      if (error%has_error()) return
      call sapt0_from_cache(mols, c, terms, error)
      call c%destroy()
   end subroutine run_sapt0

   subroutine sapt0_from_cache(mols, c, terms, error, chf_a, chf_b)
      !! The SAPT0 terms from an already-built cache
      !!
      !! Split out of `run_sapt0` so that `run_sapt2` computes its SAPT0 part
      !! through exactly this code: with the amplitudes zeroed, SAPT2 falls back
      !! to these numbers identically.
      type(sapt_molecules_t), intent(inout) :: mols
      type(sapt_cache_t), intent(in) :: c
      type(sapt_terms_t), intent(out) :: terms
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: chf_a(:, :), chf_b(:, :)

      type(rhf_result_t) :: scf_d
      real(dp) :: four

      terms%elst10 = sapt_elst10(c)
      terms%exch10_s2 = sapt_exch10_s2(c)
      terms%exch10 = sapt_exch10(c, error)
      if (error%has_error()) return
      call sapt_induction(mols, c, terms, error, chf_a=chf_a, chf_b=chf_b)
      if (error%has_error()) return
      terms%disp20 = sapt_disp20(c)
      terms%exch_disp20 = sapt_exch_disp20(c)

      ! The supermolecular reference, counterpoise-corrected for free because
      ! the monomer energies came from SCFs in the dimer basis.
      call run_czt_rhf(mols%dimer, mols%n_elec_a + mols%n_elec_b, 200, &
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
   end subroutine sapt0_from_cache

   pure function sapt_elst10(c) result(energy)
      !! `Elst10,r` -- the classical Coulomb interaction of the two unperturbed
      !! monomer densities, nuclei included (`sapt_jk_terms.py:131-134`):
      !!
      !!     Elst10 = 4 D_B . J_A + 2 D_A . V_B + 2 D_B . V_A + dE_nuc
      !!
      !! Each contraction is elementwise -- `sum_pq X_pq Y_pq`, psi4's
      !! `vector_dot`. Every operand here is symmetric so it coincides with
      !! `Tr(X Y)`, but the exchange terms have operands that are not, so the
      !! elementwise convention is used uniformly.
      type(sapt_cache_t), intent(in) :: c
      real(dp) :: energy

      energy = 4.0_dp*sum(c%d_b*c%j_a) &
               + 2.0_dp*sum(c%d_a*c%v_b) &
               + 2.0_dp*sum(c%d_b*c%v_a) &
               + c%e_nuc
   end function sapt_elst10

   ! ------------------------------------------------------------------
   ! SAPT2: four intramonomer-correlation corrections on top of SAPT0.
   !
   ! Everything below transcribes psi4's libsapt_solver -- amplitudes.cc,
   ! elst12.cc, exch11.cc, exch12.cc, ind22.cc and the uAR/uBS potentials
   ! of exch-ind20.cc -- against the conventional reference implementation
   ! in validation/check_sapt2.py, which carries the line-by-line citations.
   ! psi4's density-fitting factors are replaced by an exact factorization
   ! of the conventional ERIs, so the transcription is verbatim while the
   ! numbers are what a four-index code must produce. No frozen core.
   ! ------------------------------------------------------------------

   subroutine build_sapt2_cache(c, c2, error)
      !! The three-index factors, one-electron MO matrices and omega blocks
      !!
      !! The factorization diagonalizes the packed AO ERI matrix and keeps the
      !! positive spectrum: `(pq|rs) = sum_k B^k_pq B^k_rs` to machine precision,
      !! because the ERI matrix is a Gram matrix and its negative eigenvalues are
      !! roundoff. Rank is about `nao(nao+1)/2`, so this is an in-core object on
      !! the same footing as the cache's full `eri`, and the eigensolve needs a
      !! second copy of the packed matrix on top of it.
      type(sapt_cache_t), intent(in) :: c
      type(sapt2_cache_t), intent(out) :: c2
      type(error_t), intent(inout) :: error

      real(dp), parameter :: EIGEN_KEEP = 1.0e-15_dp
      real(dp), allocatable :: work(:, :), vals(:), m(:, :)
      real(dp), allocatable :: half_a(:, :), half_b(:, :)
      real(dp), allocatable :: full_aa(:, :), full_ab(:, :), full_bb(:, :)
      real(dp), allocatable :: tmp(:, :)
      integer :: nao, n_pair, no_a, no_b, nv_a, nv_b, nmo_a, nmo_b
      integer :: i, j, k, kk, info, nd3
      real(dp) :: biggest, scal

      nao = c%nao
      no_a = c%nocc_a
      no_b = c%nocc_b
      nmo_a = size(c%c_a, 2)
      nmo_b = size(c%c_b, 2)
      nv_a = nmo_a - no_a
      nv_b = nmo_b - no_b
      n_pair = nao*(nao + 1)/2

      allocate (work(n_pair, n_pair), vals(n_pair))
      work = c%eri_packed
      call pic_syev(work, vals, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "sapt2: the AO ERI matrix could "// &
                        "not be diagonalised for the exact factorization")
         return
      end if
      biggest = maxval(abs(vals))
      c2%ndf = count(vals > biggest*EIGEN_KEEP)
      nd3 = c2%ndf + 3
      c2%nd3 = nd3
      c2%nocc_a = no_a
      c2%nvir_a = nv_a
      c2%nmo_a = nmo_a
      c2%nocc_b = no_b
      c2%nvir_b = nv_b
      c2%nmo_b = nmo_b
      c2%nea = 2*no_a
      c2%neb = 2*no_b
      c2%enuc = c%e_nuc
      ! The intermolecular nuclear repulsion is positive whenever both monomers
      ! carry nuclei, so the square root is safe.
      c2%esq = sqrt(c%e_nuc/real(c2%nea*c2%neb, dp))
      c2%eps_a = c%eps_a
      c2%eps_b = c%eps_b

      allocate (c2%b_aa(no_a, no_a, nd3), c2%b_ar(no_a, nv_a, nd3), &
                c2%b_rr(nv_a, nv_a, nd3), c2%b_bb(no_b, no_b, nd3), &
                c2%b_bs(no_b, nv_b, nd3), c2%b_ss(nv_b, nv_b, nd3), &
                c2%b_ab(no_a, no_b, nd3), c2%b_as(no_a, nv_b, nd3), &
                c2%b_rb(nv_a, no_b, nd3))
      c2%b_aa = 0.0_dp
      c2%b_ar = 0.0_dp
      c2%b_rr = 0.0_dp
      c2%b_bb = 0.0_dp
      c2%b_bs = 0.0_dp
      c2%b_ss = 0.0_dp
      c2%b_ab = 0.0_dp
      c2%b_as = 0.0_dp
      c2%b_rb = 0.0_dp

      allocate (m(nao, nao), half_a(nao, nmo_a), half_b(nao, nmo_b))
      allocate (full_aa(nmo_a, nmo_a), full_ab(nmo_a, nmo_b), &
                full_bb(nmo_b, nmo_b))
      kk = 0
      do k = 1, n_pair
         if (vals(k) <= biggest*EIGEN_KEEP) cycle
         kk = kk + 1
         scal = sqrt(vals(k))
         do j = 1, nao
            do i = 1, nao
               m(i, j) = work(pair_index(i, j), k)*scal
            end do
         end do
         call pic_gemm(m, c%c_a, half_a)
         call pic_gemm(m, c%c_b, half_b)
         call pic_gemm(c%c_a, half_a, full_aa, transa="T")
         call pic_gemm(c%c_a, half_b, full_ab, transa="T")
         call pic_gemm(c%c_b, half_b, full_bb, transa="T")
         c2%b_aa(:, :, kk) = full_aa(1:no_a, 1:no_a)
         c2%b_ar(:, :, kk) = full_aa(1:no_a, no_a + 1:nmo_a)
         c2%b_rr(:, :, kk) = full_aa(no_a + 1:nmo_a, no_a + 1:nmo_a)
         c2%b_bb(:, :, kk) = full_bb(1:no_b, 1:no_b)
         c2%b_bs(:, :, kk) = full_bb(1:no_b, no_b + 1:nmo_b)
         c2%b_ss(:, :, kk) = full_bb(no_b + 1:nmo_b, no_b + 1:nmo_b)
         c2%b_ab(:, :, kk) = full_ab(1:no_a, 1:no_b)
         c2%b_as(:, :, kk) = full_ab(1:no_a, no_b + 1:nmo_b)
         c2%b_rb(:, :, kk) = full_ab(no_a + 1:nmo_a, 1:no_b)
      end do
      deallocate (work, vals, m, half_a, half_b, full_aa, full_ab, full_bb)

      ! One-electron matrices over the full MO ranges (sapt.cc:228-271).
      allocate (tmp(nao, nmo_b))
      call pic_gemm(c%s, c%c_b, tmp)
      allocate (c2%s_ab(nmo_a, nmo_b))
      call pic_gemm(c%c_a, tmp, c2%s_ab, transa="T")
      call pic_gemm(c%v_a, c%c_b, tmp)
      allocate (c2%v_aab(nmo_a, nmo_b))
      call pic_gemm(c%c_a, tmp, c2%v_aab, transa="T")
      call pic_gemm(c%v_b, c%c_b, tmp)
      allocate (c2%v_bab(nmo_a, nmo_b))
      call pic_gemm(c%c_a, tmp, c2%v_bab, transa="T")
      deallocate (tmp)
      allocate (tmp(nao, nmo_a))
      call pic_gemm(c%v_b, c%c_a, tmp)
      allocate (c2%v_baa(nmo_a, nmo_a))
      call pic_gemm(c%c_a, tmp, c2%v_baa, transa="T")
      deallocate (tmp)
      allocate (tmp(nao, nmo_b))
      call pic_gemm(c%v_a, c%c_b, tmp)
      allocate (c2%v_abb(nmo_b, nmo_b))
      call pic_gemm(c%c_b, tmp, c2%v_abb, transa="T")
      deallocate (tmp)

      ! Omega integrals: w = v + 2J of the partner, in MO blocks
      ! (sapt2.cc:815-905). The cache already holds w in the AO basis.
      allocate (tmp(nao, nmo_a))
      call pic_gemm(c%w_b, c%c_a, tmp)
      allocate (full_aa(nmo_a, nmo_a))
      call pic_gemm(c%c_a, tmp, full_aa, transa="T")
      c2%w_baa = full_aa(1:no_a, 1:no_a)
      c2%w_bar = full_aa(1:no_a, no_a + 1:nmo_a)
      c2%w_brr = full_aa(no_a + 1:nmo_a, no_a + 1:nmo_a)
      deallocate (tmp, full_aa)
      allocate (tmp(nao, nmo_b))
      call pic_gemm(c%w_a, c%c_b, tmp)
      allocate (full_bb(nmo_b, nmo_b))
      call pic_gemm(c%c_b, tmp, full_bb, transa="T")
      c2%w_abb = full_bb(1:no_b, 1:no_b)
      c2%w_abs = full_bb(1:no_b, no_b + 1:nmo_b)
      c2%w_ass = full_bb(no_b + 1:nmo_b, no_b + 1:nmo_b)
      deallocate (tmp, full_bb)

      ! The dressed diagonal vectors (sapt2.cc:815-833).
      allocate (c2%diag_aa(nd3), c2%diag_bb(nd3))
      c2%diag_aa = 0.0_dp
      c2%diag_bb = 0.0_dp
      do k = 1, c2%ndf
         do i = 1, no_a
            c2%diag_aa(k) = c2%diag_aa(k) + c2%b_aa(i, i, k)
         end do
         do i = 1, no_b
            c2%diag_bb(k) = c2%diag_bb(k) + c2%b_bb(i, i, k)
         end do
      end do
      c2%diag_aa(c2%ndf + 1) = real(no_a, dp)
      do i = 1, no_a
         c2%diag_aa(c2%ndf + 2) = c2%diag_aa(c2%ndf + 2) &
                                  + c2%v_baa(i, i)/real(c2%neb, dp)
      end do
      c2%diag_aa(c2%ndf + 3) = real(no_a, dp)*c2%esq
      do i = 1, no_b
         c2%diag_bb(c2%ndf + 1) = c2%diag_bb(c2%ndf + 1) &
                                  + c2%v_abb(i, i)/real(c2%nea, dp)
      end do
      c2%diag_bb(c2%ndf + 2) = real(no_b, dp)
      c2%diag_bb(c2%ndf + 3) = real(no_b, dp)*c2%esq
   end subroutine build_sapt2_cache

   ! -- the dressed blocks, transcribing utils.cc get_*_ints ---------------
   ! A-monomer pairs dress as (delta, vB/NB, e*delta), B-monomer pairs as
   ! (vA/NA, delta, e*delta); mixed pairs carry sAB where a same-monomer pair
   ! carries delta, and the dress flag selects which side of the convention
   ! the block plays.

   function dressed_aa(c2) result(a)
      type(sapt2_cache_t), intent(in) :: c2
      real(dp), allocatable :: a(:, :, :)

      integer :: i, nf

      a = c2%b_aa
      nf = c2%ndf
      a(:, :, nf + 2) = c2%v_baa(1:c2%nocc_a, 1:c2%nocc_a)/real(c2%neb, dp)
      do i = 1, c2%nocc_a
         a(i, i, nf + 1) = 1.0_dp
         a(i, i, nf + 3) = c2%esq
      end do
   end function dressed_aa

   function dressed_bb(c2) result(a)
      type(sapt2_cache_t), intent(in) :: c2
      real(dp), allocatable :: a(:, :, :)

      integer :: i, nf

      a = c2%b_bb
      nf = c2%ndf
      a(:, :, nf + 1) = c2%v_abb(1:c2%nocc_b, 1:c2%nocc_b)/real(c2%nea, dp)
      do i = 1, c2%nocc_b
         a(i, i, nf + 2) = 1.0_dp
         a(i, i, nf + 3) = c2%esq
      end do
   end function dressed_bb

   function dressed_ar(c2) result(a)
      type(sapt2_cache_t), intent(in) :: c2
      real(dp), allocatable :: a(:, :, :)

      a = c2%b_ar
      a(:, :, c2%ndf + 2) = c2%v_baa(1:c2%nocc_a, c2%nocc_a + 1:c2%nmo_a) &
                            /real(c2%neb, dp)
   end function dressed_ar

   function dressed_bs(c2) result(a)
      type(sapt2_cache_t), intent(in) :: c2
      real(dp), allocatable :: a(:, :, :)

      a = c2%b_bs
      a(:, :, c2%ndf + 1) = c2%v_abb(1:c2%nocc_b, c2%nocc_b + 1:c2%nmo_b) &
                            /real(c2%nea, dp)
   end function dressed_bs

   function dressed_rr(c2) result(a)
      type(sapt2_cache_t), intent(in) :: c2
      real(dp), allocatable :: a(:, :, :)

      integer :: i, nf

      a = c2%b_rr
      nf = c2%ndf
      a(:, :, nf + 2) = c2%v_baa(c2%nocc_a + 1:c2%nmo_a, &
                                 c2%nocc_a + 1:c2%nmo_a)/real(c2%neb, dp)
      do i = 1, c2%nvir_a
         a(i, i, nf + 1) = 1.0_dp
         a(i, i, nf + 3) = c2%esq
      end do
   end function dressed_rr

   function dressed_ss(c2) result(a)
      type(sapt2_cache_t), intent(in) :: c2
      real(dp), allocatable :: a(:, :, :)

      integer :: i, nf

      a = c2%b_ss
      nf = c2%ndf
      a(:, :, nf + 1) = c2%v_abb(c2%nocc_b + 1:c2%nmo_b, &
                                 c2%nocc_b + 1:c2%nmo_b)/real(c2%nea, dp)
      do i = 1, c2%nvir_b
         a(i, i, nf + 2) = 1.0_dp
         a(i, i, nf + 3) = c2%esq
      end do
   end function dressed_ss

   function dressed_ab(c2, dress) result(a)
      type(sapt2_cache_t), intent(in) :: c2
      integer, intent(in) :: dress
      real(dp), allocatable :: a(:, :, :)

      integer :: nf

      a = c2%b_ab
      nf = c2%ndf
      associate (s => c2%s_ab(1:c2%nocc_a, 1:c2%nocc_b))
         if (dress == 1) then
            a(:, :, nf + 1) = s
            a(:, :, nf + 2) = c2%v_bab(1:c2%nocc_a, 1:c2%nocc_b) &
                              /real(c2%neb, dp)
         else
            a(:, :, nf + 1) = c2%v_aab(1:c2%nocc_a, 1:c2%nocc_b) &
                              /real(c2%nea, dp)
            a(:, :, nf + 2) = s
         end if
         a(:, :, nf + 3) = c2%esq*s
      end associate
   end function dressed_ab

   function dressed_as(c2, dress) result(a)
      type(sapt2_cache_t), intent(in) :: c2
      integer, intent(in) :: dress
      real(dp), allocatable :: a(:, :, :)

      integer :: nf

      a = c2%b_as
      nf = c2%ndf
      associate (s => c2%s_ab(1:c2%nocc_a, c2%nocc_b + 1:c2%nmo_b))
         if (dress == 1) then
            a(:, :, nf + 1) = s
            a(:, :, nf + 2) = c2%v_bab(1:c2%nocc_a, c2%nocc_b + 1:c2%nmo_b) &
                              /real(c2%neb, dp)
         else
            a(:, :, nf + 1) = c2%v_aab(1:c2%nocc_a, c2%nocc_b + 1:c2%nmo_b) &
                              /real(c2%nea, dp)
            a(:, :, nf + 2) = s
         end if
         a(:, :, nf + 3) = c2%esq*s
      end associate
   end function dressed_as

   function dressed_rb(c2, dress) result(a)
      type(sapt2_cache_t), intent(in) :: c2
      integer, intent(in) :: dress
      real(dp), allocatable :: a(:, :, :)

      integer :: nf

      a = c2%b_rb
      nf = c2%ndf
      associate (s => c2%s_ab(c2%nocc_a + 1:c2%nmo_a, 1:c2%nocc_b))
         if (dress == 1) then
            a(:, :, nf + 1) = c2%v_aab(c2%nocc_a + 1:c2%nmo_a, 1:c2%nocc_b) &
                              /real(c2%nea, dp)
            a(:, :, nf + 2) = s
         else
            a(:, :, nf + 1) = s
            a(:, :, nf + 2) = c2%v_bab(c2%nocc_a + 1:c2%nmo_a, 1:c2%nocc_b) &
                              /real(c2%neb, dp)
         end if
         a(:, :, nf + 3) = c2%esq*s
      end associate
   end function dressed_rb

   subroutine build_sapt2_amps(c2, side, amps)
      !! One monomer's amplitude set -- `SAPT2::amplitudes()`, the labels the
      !! four SAPT2 terms consume: `tARAR`, `Theta AR`, `pAA`/`pRR`, `Y2 AR`,
      !! `T2 AR`, `t2ARAR` and `Theta 2 AR` (and the BS mirrors)
      type(sapt2_cache_t), intent(in) :: c2
      character(len=1), intent(in) :: side
      type(sapt2_amps_t), intent(out) :: amps

      real(dp), allocatable :: ovd(:, :, :)

      if (side == "A") then
         ovd = dressed_ar(c2)
         call amps_from_blocks(c2%b_aa, c2%b_ar, c2%b_rr, ovd, c2%eps_a, &
                               c2%nocc_a, c2%nvir_a, c2%ndf, amps)
      else
         ovd = dressed_bs(c2)
         call amps_from_blocks(c2%b_bb, c2%b_bs, c2%b_ss, ovd, c2%eps_b, &
                               c2%nocc_b, c2%nvir_b, c2%ndf, amps)
      end if
      deallocate (ovd)
   end subroutine build_sapt2_amps

   subroutine amps_from_blocks(b_oo, b_ov, b_vv, ovd, eps, no, nv, ndf, amps)
      !! The amplitude ladder for one monomer, `amplitudes.cc:42-105`
      real(dp), intent(in) :: b_oo(:, :, :), b_ov(:, :, :), b_vv(:, :, :)
      real(dp), intent(in) :: ovd(:, :, :)   !! The dressed occ-vir block
      real(dp), intent(in) :: eps(:)
      integer, intent(in) :: no, nv, ndf
      type(sapt2_amps_t), intent(out) :: amps

      real(dp), allocatable :: v_oovv(:, :, :, :), v_oooo(:, :, :, :)
      real(dp), allocatable :: v_vvvv(:, :, :, :), n(:, :, :, :), th2(:, :, :, :)
      real(dp), allocatable :: xk(:)
      integer :: a, r, c, q, x, y, k
      real(dp) :: denom, acc

      ! tOVOV (amplitudes.cc:106-137): the denominator is INCLUDED and
      ! nothing is antisymmetrized. The contraction runs over the exact
      ! factors only -- psi4's GEMM there is over ndf_, not ndf_+3.
      allocate (amps%t(no, nv, no, nv))
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  denom = eps(a) + eps(c) - eps(no + r) - eps(no + q)
                  amps%t(a, r, c, q) = dot_product(b_ov(a, r, 1:ndf), &
                                                   b_ov(c, q, 1:ndf))/denom
               end do
            end do
         end do
      end do

      ! theta = 2t - t with the occupieds swapped (utils.cc:1326).
      allocate (amps%theta(no, nv, no, nv))
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  amps%theta(a, r, c, q) = 2.0_dp*amps%t(a, r, c, q) &
                                           - amps%t(c, r, a, q)
               end do
            end do
         end do
      end do

      ! Theta^P = theta . dressed OV (amplitudes.cc:167-200).
      allocate (amps%th_ov(no, nv, size(ovd, 3)))
      do k = 1, size(ovd, 3)
         do r = 1, nv
            do a = 1, no
               amps%th_ov(a, r, k) = sum(amps%theta(a, r, :, :)*ovd(:, :, k))
            end do
         end do
      end do

      ! The unrelaxed density blocks (amplitudes.cc:138-166), WITHOUT the
      ! -2/+2 -- those factors live in the consuming terms, as in psi4.
      allocate (amps%p_oo(no, no), amps%p_vv(nv, nv))
      do c = 1, no
         do a = 1, no
            amps%p_oo(a, c) = sum(amps%t(a, :, :, :)*amps%theta(c, :, :, :))
         end do
      end do
      do q = 1, nv
         do r = 1, nv
            amps%p_vv(r, q) = sum(amps%t(:, :, :, r)*amps%theta(:, :, :, q))
         end do
      end do

      ! Y2 (amplitudes.cc:201-317). The Y2_3 part first; the T2 singles are
      ! written from it ALONE, before Y2_1/Y2_2 are added (:216).
      allocate (amps%y2(no, nv), amps%t_singles(no, nv))
      do r = 1, nv
         do a = 1, no
            acc = 0.0_dp
            do c = 1, no
               acc = acc - dot_product(b_oo(a, c, :), amps%th_ov(c, r, :))
            end do
            do q = 1, nv
               acc = acc + dot_product(amps%th_ov(a, q, :), b_vv(r, q, :))
            end do
            amps%y2(a, r) = acc
            amps%t_singles(a, r) = acc/(eps(a) - eps(no + r))
         end do
      end do
      ! Y2_1 (:230-262)
      allocate (xk(size(b_vv, 3)))
      do k = 1, size(xk)
         xk(k) = sum(b_vv(:, :, k)*amps%p_vv)
      end do
      do r = 1, nv
         do a = 1, no
            acc = 2.0_dp*dot_product(b_ov(a, r, :), xk)
            do x = 1, nv
               do q = 1, nv
                  acc = acc - amps%p_vv(q, x) &
                        *dot_product(b_ov(a, q, :), b_vv(r, x, :))
               end do
            end do
            amps%y2(a, r) = amps%y2(a, r) + acc
         end do
      end do
      ! Y2_2 (:263-291)
      do k = 1, size(xk)
         xk(k) = sum(b_oo(:, :, k)*amps%p_oo)
      end do
      do r = 1, nv
         do a = 1, no
            acc = -2.0_dp*dot_product(b_ov(a, r, :), xk)
            do x = 1, no
               do c = 1, no
                  acc = acc + amps%p_oo(c, x) &
                        *dot_product(b_oo(x, a, :), b_ov(c, r, :))
               end do
            end do
            amps%y2(a, r) = amps%y2(a, r) + acc
         end do
      end do
      deallocate (xk)

      ! t2OVOV (amplitudes.cc:318-395): the second-order doubles.
      allocate (v_oovv(no, no, nv, nv))
      do q = 1, nv
         do r = 1, nv
            do c = 1, no
               do a = 1, no
                  v_oovv(a, c, r, q) = dot_product(b_oo(a, c, 1:ndf), &
                                                   b_vv(r, q, 1:ndf))
               end do
            end do
         end do
      end do
      allocate (v_oooo(no, no, no, no))
      do y = 1, no
         do c = 1, no
            do x = 1, no
               do a = 1, no
                  v_oooo(a, x, c, y) = dot_product(b_oo(a, x, 1:ndf), &
                                                   b_oo(c, y, 1:ndf))
               end do
            end do
         end do
      end do
      allocate (v_vvvv(nv, nv, nv, nv))
      do y = 1, nv
         do q = 1, nv
            do x = 1, nv
               do r = 1, nv
                  v_vvvv(r, x, q, y) = dot_product(b_vv(r, x, 1:ndf), &
                                                   b_vv(q, y, 1:ndf))
               end do
            end do
         end do
      end do

      allocate (n(no, nv, no, nv))
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  acc = dot_product(b_ov(a, r, :), amps%th_ov(c, q, :))
                  do y = 1, nv
                     do x = 1, no
                        acc = acc - amps%t(x, r, c, y)*v_oovv(a, x, q, y) &
                              - amps%t(a, r, x, y)*v_oovv(c, x, q, y)
                     end do
                  end do
                  do y = 1, no
                     do x = 1, no
                        acc = acc + 0.5_dp*v_oooo(a, x, c, y)*amps%t(x, r, y, q)
                     end do
                  end do
                  do y = 1, nv
                     do x = 1, nv
                        acc = acc + 0.5_dp*v_vvvv(r, x, q, y)*amps%t(a, x, c, y)
                     end do
                  end do
                  n(a, r, c, q) = acc
               end do
            end do
         end do
      end do
      deallocate (v_oovv, v_oooo, v_vvvv)

      allocate (amps%t2(no, nv, no, nv))
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  denom = eps(a) + eps(c) - eps(no + r) - eps(no + q)
                  amps%t2(a, r, c, q) = (n(a, r, c, q) + n(c, q, a, r))/denom
               end do
            end do
         end do
      end do
      deallocate (n)

      ! Theta 2: the same theta() applied to t2 (amplitudes.cc:100-105).
      allocate (th2(no, nv, no, nv))
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  th2(a, r, c, q) = 2.0_dp*amps%t2(a, r, c, q) &
                                    - amps%t2(c, r, a, q)
               end do
            end do
         end do
      end do
      allocate (amps%th2_ov(no, nv, size(ovd, 3)))
      do k = 1, size(ovd, 3)
         do r = 1, nv
            do a = 1, no
               amps%th2_ov(a, r, k) = sum(th2(a, r, :, :)*ovd(:, :, k))
            end do
         end do
      end do
      deallocate (th2)
   end subroutine amps_from_blocks

   subroutine sapt2_zero_amps(amps)
      !! Zero every stored amplitude, for the t2 -> 0 fallback test
      !!
      !! With these, `Elst12`, `Exch11`, `Exch12` and `Ind22`'s pieces 2-7
      !! vanish identically. `Ind22`'s first piece does NOT: its perturbed
      !! doubles are built from bare integrals times the uncoupled response,
      !! not from anything stored here -- see `sapt_ind22`.
      type(sapt2_amps_t), intent(inout) :: amps

      amps%t = 0.0_dp
      amps%theta = 0.0_dp
      amps%th_ov = 0.0_dp
      amps%p_oo = 0.0_dp
      amps%p_vv = 0.0_dp
      amps%y2 = 0.0_dp
      amps%t_singles = 0.0_dp
      amps%t2 = 0.0_dp
      amps%th2_ov = 0.0_dp
   end subroutine sapt2_zero_amps

   function sapt2_amps_mp2_energy(c2, amps, side) result(e2)
      !! The monomer MP2 correlation energy the amplitudes imply
      !!
      !! `E2 = sum t theta D`, since `v = t D` elementwise and the denominator is
      !! symmetric under the occupied swap. Must equal `run_czt_mp2` on the
      !! same ghosted monomer to machine precision: the sharpest check available
      !! on the amplitude conventions.
      type(sapt2_cache_t), intent(in) :: c2
      type(sapt2_amps_t), intent(in) :: amps
      character(len=1), intent(in) :: side
      real(dp) :: e2

      integer :: a, r, c, q, no, nv
      real(dp) :: denom

      if (side == "A") then
         no = c2%nocc_a
         nv = c2%nvir_a
      else
         no = c2%nocc_b
         nv = c2%nvir_b
      end if

      e2 = 0.0_dp
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  if (side == "A") then
                     denom = c2%eps_a(a) + c2%eps_a(c) &
                             - c2%eps_a(no + r) - c2%eps_a(no + q)
                  else
                     denom = c2%eps_b(a) + c2%eps_b(c) &
                             - c2%eps_b(no + r) - c2%eps_b(no + q)
                  end if
                  e2 = e2 + amps%t(a, r, c, q)*amps%theta(a, r, c, q)*denom
               end do
            end do
         end do
      end do
   end function sapt2_amps_mp2_energy

   subroutine sapt2_k2f(c2, u_ar, v_ar, u_bs, v_bs)
      !! The `uAR`/`uBS` ("Exch12 K2f") and `vAR`/`vBS` ("Exch-Ind")
      !! potentials -- `exch_ind20rA_B`/`exch_ind20rB_A`,
      !! `exch-ind20.cc:682-1145`. Thirteen shared pieces per side; the two
      !! outputs differ only in a few coefficients.
      !!
      !! `vAR` must equal minus the USAPT0-factorized exchange-induction
      !! potential of `exch_ind_potential`, which the tests assert.
      type(sapt2_cache_t), intent(in) :: c2
      real(dp), allocatable, intent(out) :: u_ar(:, :), v_ar(:, :)
      real(dp), allocatable, intent(out) :: u_bs(:, :), v_bs(:, :)

      real(dp), allocatable :: ab1(:, :, :), ab2(:, :, :), rb1(:, :, :)
      real(dp), allocatable :: as1(:, :, :), aa1(:, :, :), ar1(:, :, :)
      real(dp), allocatable :: bb1(:, :, :), bs1(:, :, :)
      real(dp), allocatable :: p1(:, :), p2(:, :), p3(:, :), p4(:, :)
      real(dp), allocatable :: p5(:, :), p6(:, :), p7(:, :), p8(:, :)
      real(dp), allocatable :: p9(:, :), p10(:, :), p11(:, :), p12(:, :)
      real(dp), allocatable :: p13(:, :), cp(:), x_oo(:, :), x_mn(:, :)
      integer :: no_a, no_b, nv_a, nv_b, a, b, c, d, r, s, e
      real(dp) :: acc

      no_a = c2%nocc_a
      no_b = c2%nocc_b
      nv_a = c2%nvir_a
      nv_b = c2%nvir_b
      ab1 = dressed_ab(c2, 1)
      ab2 = dressed_ab(c2, 2)
      rb1 = dressed_rb(c2, 1)
      as1 = dressed_as(c2, 1)
      aa1 = dressed_aa(c2)
      ar1 = dressed_ar(c2)
      bb1 = dressed_bb(c2)
      bs1 = dressed_bs(c2)

      associate (s_oo => c2%s_ab(1:no_a, 1:no_b), &
                 s_vo => c2%s_ab(no_a + 1:c2%nmo_a, 1:no_b), &
                 s_ov => c2%s_ab(1:no_a, no_b + 1:c2%nmo_b))

         ! -- the A side, exch-ind20.cc:682-911 --------------------------
         allocate (p1(no_a, nv_a), p2(no_a, nv_a), p3(no_a, nv_a), &
                   p4(no_a, nv_a), p5(no_a, nv_a), p6(no_a, nv_a), &
                   p7(no_a, nv_a), p8(no_a, nv_a), p9(no_a, nv_a), &
                   p10(no_a, nv_a), p11(no_a, nv_a), p12(no_a, nv_a), &
                   p13(no_a, nv_a), cp(c2%nd3))
         do r = 1, nv_a
            do a = 1, no_a
               acc = 0.0_dp
               do b = 1, no_b
                  acc = acc + dot_product(ab1(a, b, :), rb1(r, b, :))
               end do
               p1(a, r) = acc
               acc = 0.0_dp
               do b = 1, no_b
                  acc = acc + s_oo(a, b)*dot_product(rb1(r, b, :), c2%diag_aa)
               end do
               p2(a, r) = acc
               acc = 0.0_dp
               do b = 1, no_b
                  do c = 1, no_a
                     acc = acc + s_oo(c, b)*dot_product(aa1(a, c, :), &
                                                        rb1(r, b, :))
                  end do
               end do
               p3(a, r) = acc
            end do
         end do
         do e = 1, c2%nd3
            cp(e) = sum(s_oo*ab2(:, :, e))
         end do
         do r = 1, nv_a
            do a = 1, no_a
               p4(a, r) = dot_product(ar1(a, r, :), cp)
               acc = 0.0_dp
               do b = 1, no_b
                  do c = 1, no_a
                     ! psi4's t5: sum_a AB2(a,b) AR1(a,r) against S(this,b)
                     acc = acc + s_oo(a, b)*dot_product(ab2(c, b, :), &
                                                        ar1(c, r, :))
                  end do
               end do
               p5(a, r) = acc
               acc = 0.0_dp
               do b = 1, no_b
                  acc = acc + dot_product(ab1(a, b, :), c2%diag_bb)*s_vo(r, b)
               end do
               p6(a, r) = acc
               acc = 0.0_dp
               do b = 1, no_b
                  do d = 1, no_b
                     acc = acc + dot_product(ab1(a, d, :), bb1(b, d, :)) &
                           *s_vo(r, b)
                  end do
               end do
               p7(a, r) = acc
            end do
         end do
         allocate (x_oo(no_b, no_b))
         do d = 1, no_b
            do b = 1, no_b
               x_oo(b, d) = dot_product(bb1(b, d, :), c2%diag_aa)
            end do
         end do
         do r = 1, nv_a
            do a = 1, no_a
               acc = 0.0_dp
               do b = 1, no_b
                  do d = 1, no_b
                     acc = acc + s_oo(a, d)*x_oo(b, d)*s_vo(r, b)
                  end do
               end do
               p8(a, r) = acc
            end do
         end do
         deallocate (x_oo)
         allocate (x_mn(no_a, no_a))
         do c = 1, no_a
            do a = 1, no_a
               x_mn(a, c) = dot_product(aa1(a, c, :), c2%diag_bb)
            end do
         end do
         do r = 1, nv_a
            do a = 1, no_a
               acc = 0.0_dp
               do c = 1, no_a
                  do b = 1, no_b
                     acc = acc + x_mn(a, c)*s_oo(c, b)*s_vo(r, b)
                  end do
               end do
               p9(a, r) = acc
               acc = 0.0_dp
               do d = 1, no_b
                  do b = 1, no_b
                     do c = 1, no_a
                        acc = acc + s_oo(c, b)*s_vo(r, d) &
                              *dot_product(aa1(a, c, :), bb1(d, b, :))
                     end do
                  end do
               end do
               p10(a, r) = acc
            end do
         end do
         deallocate (x_mn)
         allocate (x_oo(no_b, no_b))
         do d = 1, no_b
            do b = 1, no_b
               x_oo(b, d) = dot_product(s_oo(:, b), s_oo(:, d))
            end do
         end do
         do e = 1, c2%nd3
            cp(e) = sum(bb1(:, :, e)*x_oo)
         end do
         do r = 1, nv_a
            do a = 1, no_a
               p11(a, r) = dot_product(ar1(a, r, :), cp)
            end do
         end do
         deallocate (x_oo)
         allocate (x_mn(no_a, nv_a))
         do r = 1, nv_a
            do a = 1, no_a
               x_mn(a, r) = dot_product(ar1(a, r, :), c2%diag_bb)
            end do
         end do
         do r = 1, nv_a
            do a = 1, no_a
               acc = 0.0_dp
               do c = 1, no_a
                  acc = acc + dot_product(s_oo(a, :), s_oo(c, :))*x_mn(c, r)
               end do
               p12(a, r) = acc
               acc = 0.0_dp
               do d = 1, no_b
                  do c = 1, no_a
                     do b = 1, no_b
                        acc = acc + s_oo(a, b)*s_oo(c, d) &
                              *dot_product(bb1(b, d, :), ar1(c, r, :))
                     end do
                  end do
               end do
               p13(a, r) = acc
            end do
         end do
         deallocate (x_mn)

         u_ar = p1 + 2.0_dp*p2 - p3 + 4.0_dp*p4 - p5 + 2.0_dp*p6 - p7 &
                - 4.0_dp*p8 - 4.0_dp*p9 + 2.0_dp*p10 - 4.0_dp*p11 &
                - 4.0_dp*p12 + 2.0_dp*p13
         v_ar = p1 + 2.0_dp*p2 - p3 + 2.0_dp*p4 - p5 + 2.0_dp*p6 - p7 &
                - 2.0_dp*p8 - 2.0_dp*p9 + p10 - 2.0_dp*p11 &
                - 2.0_dp*p12 + p13
         deallocate (p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13)

         ! -- the B side, exch-ind20.cc:913-1145 -------------------------
         allocate (p1(no_b, nv_b), p2(no_b, nv_b), p3(no_b, nv_b), &
                   p4(no_b, nv_b), p5(no_b, nv_b), p6(no_b, nv_b), &
                   p7(no_b, nv_b), p8(no_b, nv_b), p9(no_b, nv_b), &
                   p10(no_b, nv_b), p11(no_b, nv_b), p12(no_b, nv_b), &
                   p13(no_b, nv_b))
         do s = 1, nv_b
            do b = 1, no_b
               acc = 0.0_dp
               do a = 1, no_a
                  acc = acc + dot_product(ab2(a, b, :), as1(a, s, :))
               end do
               p1(b, s) = acc
               acc = 0.0_dp
               do a = 1, no_a
                  acc = acc + s_oo(a, b)*dot_product(as1(a, s, :), c2%diag_bb)
               end do
               p2(b, s) = acc
               acc = 0.0_dp
               do c = 1, no_b
                  do a = 1, no_a
                     acc = acc + s_oo(a, c)*dot_product(bb1(b, c, :), &
                                                        as1(a, s, :))
                  end do
               end do
               p3(b, s) = acc
            end do
         end do
         do e = 1, c2%nd3
            cp(e) = sum(s_oo*ab1(:, :, e))
         end do
         do s = 1, nv_b
            do b = 1, no_b
               p4(b, s) = dot_product(bs1(b, s, :), cp)
               acc = 0.0_dp
               do c = 1, no_b
                  do a = 1, no_a
                     acc = acc + s_oo(a, b)*dot_product(ab1(a, c, :), &
                                                        bs1(c, s, :))
                  end do
               end do
               p5(b, s) = acc
               acc = 0.0_dp
               do a = 1, no_a
                  acc = acc + dot_product(ab2(a, b, :), c2%diag_aa)*s_ov(a, s)
               end do
               p6(b, s) = acc
               ! s7 uses AB2 -- exch_ind20rB_A reassigns B_p_AB just above
               ! this term, where the A side's mirror uses dress 1.
               acc = 0.0_dp
               do c = 1, no_a
                  do a = 1, no_a
                     acc = acc + dot_product(aa1(a, c, :), ab2(a, b, :)) &
                           *s_ov(c, s)
                  end do
               end do
               p7(b, s) = acc
            end do
         end do
         allocate (x_mn(no_a, no_a))
         do c = 1, no_a
            do a = 1, no_a
               x_mn(a, c) = dot_product(aa1(a, c, :), c2%diag_bb)
            end do
         end do
         do s = 1, nv_b
            do b = 1, no_b
               acc = 0.0_dp
               do c = 1, no_a
                  do a = 1, no_a
                     acc = acc + x_mn(a, c)*s_oo(c, b)*s_ov(a, s)
                  end do
               end do
               p8(b, s) = acc
            end do
         end do
         deallocate (x_mn)
         allocate (x_oo(no_b, no_b))
         do c = 1, no_b
            do b = 1, no_b
               x_oo(b, c) = dot_product(bb1(b, c, :), c2%diag_aa)
            end do
         end do
         do s = 1, nv_b
            do b = 1, no_b
               acc = 0.0_dp
               do c = 1, no_b
                  acc = acc + x_oo(b, c)*dot_product(s_oo(:, c), s_ov(:, s))
               end do
               p9(b, s) = acc
               acc = 0.0_dp
               do c = 1, no_a
                  do d = 1, no_b
                     do e = 1, no_a
                        acc = acc + s_oo(c, d)*s_ov(e, s) &
                              *dot_product(bb1(b, d, :), aa1(e, c, :))
                     end do
                  end do
               end do
               p10(b, s) = acc
            end do
         end do
         deallocate (x_oo)
         allocate (x_mn(no_a, no_a))
         do c = 1, no_a
            do a = 1, no_a
               x_mn(a, c) = dot_product(s_oo(a, :), s_oo(c, :))
            end do
         end do
         do e = 1, c2%nd3
            cp(e) = sum(aa1(:, :, e)*x_mn)
         end do
         do s = 1, nv_b
            do b = 1, no_b
               p11(b, s) = dot_product(bs1(b, s, :), cp)
            end do
         end do
         deallocate (x_mn)
         do s = 1, nv_b
            do b = 1, no_b
               acc = 0.0_dp
               do c = 1, no_b
                  acc = acc + dot_product(s_oo(:, b), s_oo(:, c)) &
                        *dot_product(bs1(c, s, :), c2%diag_aa)
               end do
               p12(b, s) = acc
               ! s13: C_p_AB(e,d,P) = sum_a s(a,d) AA1(e,a,P); C_p_BB(d,b,P)
               ! = sum_e s(e,d) C_p_AB(e,b,P); then sum_d C_p_BB . BS1(d,s,P)
               acc = 0.0_dp
               do d = 1, no_b
                  do e = 1, no_a
                     do a = 1, no_a
                        acc = acc + s_oo(e, d)*s_oo(a, b) &
                              *dot_product(aa1(e, a, :), bs1(d, s, :))
                     end do
                  end do
               end do
               p13(b, s) = acc
            end do
         end do

         u_bs = p1 + 2.0_dp*p2 - p3 + 4.0_dp*p4 - p5 + 2.0_dp*p6 - p7 &
                - 4.0_dp*p8 - 4.0_dp*p9 + 2.0_dp*p10 - 4.0_dp*p11 &
                - 4.0_dp*p12 + 2.0_dp*p13
         v_bs = p1 + 2.0_dp*p2 - p3 + 2.0_dp*p4 - p5 + 2.0_dp*p6 - p7 &
                - 2.0_dp*p8 - 2.0_dp*p9 + p10 - 2.0_dp*p11 &
                - 2.0_dp*p12 + p13
         deallocate (p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13)
      end associate

      deallocate (ab1, ab2, rb1, as1, aa1, ar1, bb1, bs1, cp)
   end subroutine sapt2_k2f

   function sapt_elst12(c2, amps_a, amps_b, chf_a, chf_b) result(energy)
      !! `Elst12,r` -- `elst12.cc:38-95`
      !!
      !! Each half is the partner's electrostatic potential against the
      !! monomer's MP2 density: the oo and vv blocks directly, and the
      !! orbital-relaxation (ov) block through the CPHF coefficients that
      !! `Ind20,r` already solved for -- `Tr(Z w) = Tr(Y x_w)`, so no
      !! Z-vector is ever formed. Using the unrelaxed density alone would be
      !! wrong by exactly that third contraction, and nothing would look
      !! broken.
      type(sapt2_cache_t), intent(in) :: c2
      type(sapt2_amps_t), intent(in) :: amps_a, amps_b
      real(dp), intent(in) :: chf_a(:, :), chf_b(:, :)
      real(dp) :: energy

      energy = -2.0_dp*sum(amps_a%p_oo*c2%w_baa) &
               + 2.0_dp*sum(amps_a%p_vv*c2%w_brr) &
               + 4.0_dp*sum(amps_a%y2*chf_a) &
               - 2.0_dp*sum(amps_b%p_oo*c2%w_abb) &
               + 2.0_dp*sum(amps_b%p_vv*c2%w_ass) &
               + 4.0_dp*sum(amps_b%y2*chf_b)
   end function sapt_elst12

   function exch110(c2, th_ov) result(energy)
      !! First-order exchange against monomer A's correlated density --
      !! `exch110`, `exch11.cc:57-127`. Reused by `Exch12` with the Theta 2
      !! intermediates (`exch12.cc:44`), which is why the theta enters as an
      !! argument rather than through the amplitude set.
      type(sapt2_cache_t), intent(in) :: c2
      real(dp), intent(in) :: th_ov(:, :, :)
      real(dp) :: energy

      real(dp), allocatable :: ab2(:, :, :), bb1(:, :, :), rb1(:, :, :)
      real(dp), allocatable :: cab(:, :, :), cbb(:, :, :)
      real(dp) :: e1, e2, e3, e4, acc
      integer :: no_a, no_b, nv_a, a, b, c, r, k

      no_a = c2%nocc_a
      no_b = c2%nocc_b
      nv_a = c2%nvir_a
      ab2 = dressed_ab(c2, 2)
      bb1 = dressed_bb(c2)
      rb1 = dressed_rb(c2, 1)

      associate (s_oo => c2%s_ab(1:no_a, 1:no_b), &
                 s_vo => c2%s_ab(no_a + 1:c2%nmo_a, 1:no_b))
         allocate (cab(no_a, no_b, c2%nd3))
         do k = 1, c2%nd3
            do b = 1, no_b
               do a = 1, no_a
                  cab(a, b, k) = dot_product(s_vo(:, b), th_ov(a, :, k))
               end do
            end do
         end do
         e1 = -2.0_dp*sum(cab*ab2)

         allocate (cbb(no_b, no_b, c2%nd3))
         do k = 1, c2%nd3
            do c = 1, no_b
               do b = 1, no_b
                  cbb(b, c, k) = dot_product(s_oo(:, b), cab(:, c, k))
               end do
            end do
         end do
         e2 = 4.0_dp*sum(bb1*cbb)
         deallocate (cab, cbb)

         e3 = 0.0_dp
         do r = 1, nv_a
            do a = 1, no_a
               acc = 0.0_dp
               do b = 1, no_b
                  acc = acc + s_oo(a, b)*dot_product(th_ov(a, r, :), &
                                                     rb1(r, b, :))
               end do
               e3 = e3 - 2.0_dp*acc
            end do
         end do

         e4 = 0.0_dp
         do r = 1, nv_a
            do a = 1, no_a
               e4 = e4 - 8.0_dp*dot_product(s_oo(a, :), s_vo(r, :)) &
                    *dot_product(th_ov(a, r, :), c2%diag_bb)
            end do
         end do
      end associate

      energy = e1 + e2 + e3 + e4
      deallocate (ab2, bb1, rb1)
   end function exch110

   function exch101(c2, th_ov) result(energy)
      !! The B-side mirror -- `exch101`, `exch11.cc:129-194`
      type(sapt2_cache_t), intent(in) :: c2
      real(dp), intent(in) :: th_ov(:, :, :)
      real(dp) :: energy

      real(dp), allocatable :: ab1(:, :, :), aa1(:, :, :), as1(:, :, :)
      real(dp), allocatable :: cab(:, :, :), caa(:, :, :)
      real(dp) :: e1, e2, e3, e4, acc
      integer :: no_a, no_b, nv_b, a, b, c, s, k

      no_a = c2%nocc_a
      no_b = c2%nocc_b
      nv_b = c2%nvir_b
      ab1 = dressed_ab(c2, 1)
      aa1 = dressed_aa(c2)
      as1 = dressed_as(c2, 1)

      associate (s_oo => c2%s_ab(1:no_a, 1:no_b), &
                 s_ov => c2%s_ab(1:no_a, c2%nocc_b + 1:c2%nmo_b))
         allocate (cab(no_a, no_b, c2%nd3))
         do k = 1, c2%nd3
            do b = 1, no_b
               do a = 1, no_a
                  cab(a, b, k) = dot_product(s_ov(a, :), th_ov(b, :, k))
               end do
            end do
         end do
         e1 = -2.0_dp*sum(cab*ab1)

         allocate (caa(no_a, no_a, c2%nd3))
         do k = 1, c2%nd3
            do c = 1, no_a
               do a = 1, no_a
                  caa(a, c, k) = dot_product(s_oo(c, :), cab(a, :, k))
               end do
            end do
         end do
         e2 = 4.0_dp*sum(aa1*caa)
         deallocate (cab, caa)

         e3 = 0.0_dp
         do s = 1, nv_b
            do b = 1, no_b
               acc = 0.0_dp
               do a = 1, no_a
                  acc = acc + s_oo(a, b)*dot_product(th_ov(b, s, :), &
                                                     as1(a, s, :))
               end do
               e3 = e3 - 2.0_dp*acc
            end do
         end do

         e4 = 0.0_dp
         do s = 1, nv_b
            do b = 1, no_b
               e4 = e4 - 8.0_dp*dot_product(s_oo(:, b), s_ov(:, s)) &
                    *dot_product(th_ov(b, s, :), c2%diag_aa)
            end do
         end do
      end associate

      energy = e1 + e2 + e3 + e4
      deallocate (ab1, aa1, as1)
   end function exch101

   function sapt_exch11(c2, amps_a, amps_b) result(energy)
      !! `Exch11 = Exch110 + Exch101` -- `exch11.cc:37-55`
      type(sapt2_cache_t), intent(in) :: c2
      type(sapt2_amps_t), intent(in) :: amps_a, amps_b
      real(dp) :: energy

      energy = exch110(c2, amps_a%th_ov) + exch101(c2, amps_b%th_ov)
   end function sapt_exch11

   function exch111(c2, th_a, th_b) result(energy)
      !! Both monomers' theta densities exchanging -- `exch12.cc:99-155`
      type(sapt2_cache_t), intent(in) :: c2
      real(dp), intent(in) :: th_a(:, :, :), th_b(:, :, :)
      real(dp) :: energy

      real(dp) :: e1, e2, acc1, acc2
      integer :: no_a, no_b, a, b, s, k

      no_a = c2%nocc_a
      no_b = c2%nocc_b

      associate (s_oo => c2%s_ab(1:no_a, 1:no_b), &
                 s_vo => c2%s_ab(no_a + 1:c2%nmo_a, 1:no_b), &
                 s_ov => c2%s_ab(1:no_a, c2%nocc_b + 1:c2%nmo_b), &
                 s_vv => c2%s_ab(no_a + 1:c2%nmo_a, c2%nocc_b + 1:c2%nmo_b))
         e1 = 0.0_dp
         do k = 1, c2%nd3
            do b = 1, no_b
               do a = 1, no_a
                  acc1 = dot_product(s_vo(:, b), th_a(a, :, k))
                  acc2 = dot_product(s_ov(a, :), th_b(b, :, k))
                  e1 = e1 - 4.0_dp*acc1*acc2
               end do
            end do
         end do

         e2 = 0.0_dp
         do k = 1, c2%nd3
            do s = 1, size(th_b, 2)
               do b = 1, no_b
                  acc1 = 0.0_dp
                  do a = 1, no_a
                     acc1 = acc1 + s_oo(a, b)*dot_product(s_vv(:, s), &
                                                          th_a(a, :, k))
                  end do
                  e2 = e2 - 4.0_dp*acc1*th_b(b, s, k)
               end do
            end do
         end do
      end associate
      energy = e1 + e2
   end function exch111

   function exch120_k2f(c2, amps_a, u_ar) result(energy)
      !! The T2 singles against the K2f potential -- `exch12.cc:157-241`
      type(sapt2_cache_t), intent(in) :: c2
      type(sapt2_amps_t), intent(in) :: amps_a
      real(dp), intent(in) :: u_ar(:, :)
      real(dp) :: energy

      real(dp), allocatable :: rb2(:, :, :), ab2(:, :, :), bb1(:, :, :)
      real(dp), allocatable :: aa1(:, :, :), ar1(:, :, :)
      real(dp), allocatable :: cab(:, :, :), yab(:, :)
      real(dp) :: e(SAPT2_PIECES)
      real(dp) :: acc
      integer :: no_a, no_b, nv_a, a, b, c, r, x, k

      no_a = c2%nocc_a
      no_b = c2%nocc_b
      nv_a = c2%nvir_a
      rb2 = dressed_rb(c2, 2)
      ab2 = dressed_ab(c2, 2)
      bb1 = dressed_bb(c2)
      aa1 = dressed_aa(c2)
      ar1 = dressed_ar(c2)

      associate (t1 => amps_a%t_singles, s_oo => c2%s_ab(1:no_a, 1:no_b), &
                 s_vo => c2%s_ab(no_a + 1:c2%nmo_a, 1:no_b))
         e(1) = -2.0_dp*sum(t1*u_ar)

         allocate (cab(no_a, no_b, c2%nd3))
         do k = 1, c2%nd3
            do b = 1, no_b
               do a = 1, no_a
                  cab(a, b, k) = dot_product(t1(a, :), rb2(:, b, k))
               end do
            end do
         end do
         e(2) = -2.0_dp*sum(ab2*cab)

         e(3) = 0.0_dp
         do k = 1, c2%nd3
            do c = 1, no_b
               do b = 1, no_b
                  e(3) = e(3) + 2.0_dp*bb1(b, c, k) &
                         *dot_product(s_oo(:, b), cab(:, c, k))
               end do
            end do
         end do

         e(4) = 0.0_dp
         do b = 1, no_b
            do a = 1, no_a
               e(4) = e(4) - 4.0_dp*s_oo(a, b) &
                      *dot_product(cab(a, b, :), c2%diag_bb)
            end do
         end do
         deallocate (cab)

         allocate (yab(no_a, no_b))
         do b = 1, no_b
            do a = 1, no_a
               yab(a, b) = dot_product(t1(a, :), s_vo(:, b))
            end do
         end do
         e(5) = 0.0_dp
         do b = 1, no_b
            do a = 1, no_a
               e(5) = e(5) - 4.0_dp*dot_product(ab2(a, b, :), c2%diag_aa) &
                      *yab(a, b)
            end do
         end do

         ! e6: D_p_AB(c,b,P) = sum_x yAB(x,b) AA1(c,x,P)
         e(6) = 0.0_dp
         do k = 1, c2%nd3
            do b = 1, no_b
               do c = 1, no_a
                  acc = 0.0_dp
                  do x = 1, no_a
                     acc = acc + yab(x, b)*aa1(c, x, k)
                  end do
                  e(6) = e(6) + 2.0_dp*ab2(c, b, k)*acc
               end do
            end do
         end do
         deallocate (yab)

         ! e7: C_p_AA(c,a,P) = sum_r t1(c,r) AR1(a,r,P) against
         !     D_p_AA(c,a,P) = sum_b s(a,b) AB2(c,b,P)
         e(7) = 0.0_dp
         do a = 1, no_a
            do c = 1, no_a
               acc = 0.0_dp
               do k = 1, c2%nd3
                  acc = acc + dot_product(t1(c, :), ar1(a, :, k)) &
                        *dot_product(s_oo(a, :), ab2(c, :, k))
               end do
               e(7) = e(7) + 2.0_dp*acc
            end do
         end do
      end associate

      energy = sum(e)
      deallocate (rb2, ab2, bb1, aa1, ar1)
   end function exch120_k2f

   function exch102_k2f(c2, amps_b, u_bs) result(energy)
      !! The B-side mirror -- `exch12.cc:243-380`
      type(sapt2_cache_t), intent(in) :: c2
      type(sapt2_amps_t), intent(in) :: amps_b
      real(dp), intent(in) :: u_bs(:, :)
      real(dp) :: energy

      real(dp), allocatable :: as2(:, :, :), ab1(:, :, :), aa1(:, :, :)
      real(dp), allocatable :: bb1(:, :, :), bs1(:, :, :)
      real(dp), allocatable :: cab(:, :, :), yab(:, :)
      real(dp) :: e(SAPT2_PIECES)
      real(dp) :: acc
      integer :: no_a, no_b, a, b, c, x, k

      no_a = c2%nocc_a
      no_b = c2%nocc_b
      as2 = dressed_as(c2, 2)
      ab1 = dressed_ab(c2, 1)
      aa1 = dressed_aa(c2)
      bb1 = dressed_bb(c2)
      bs1 = dressed_bs(c2)

      associate (t1 => amps_b%t_singles, s_oo => c2%s_ab(1:no_a, 1:no_b), &
                 s_ov => c2%s_ab(1:no_a, c2%nocc_b + 1:c2%nmo_b))
         e(1) = -2.0_dp*sum(t1*u_bs)

         allocate (cab(no_a, no_b, c2%nd3))
         do k = 1, c2%nd3
            do b = 1, no_b
               do a = 1, no_a
                  cab(a, b, k) = dot_product(t1(b, :), as2(a, :, k))
               end do
            end do
         end do
         e(2) = -2.0_dp*sum(ab1*cab)

         e(3) = 0.0_dp
         do k = 1, c2%nd3
            do c = 1, no_a
               do a = 1, no_a
                  e(3) = e(3) + 2.0_dp*aa1(a, c, k) &
                         *dot_product(s_oo(c, :), cab(a, :, k))
               end do
            end do
         end do

         e(4) = 0.0_dp
         do b = 1, no_b
            do a = 1, no_a
               e(4) = e(4) - 4.0_dp*s_oo(a, b) &
                      *dot_product(cab(a, b, :), c2%diag_aa)
            end do
         end do
         deallocate (cab)

         allocate (yab(no_a, no_b))
         do b = 1, no_b
            do a = 1, no_a
               yab(a, b) = dot_product(s_ov(a, :), t1(b, :))
            end do
         end do
         e(5) = 0.0_dp
         do b = 1, no_b
            do a = 1, no_a
               e(5) = e(5) - 4.0_dp*dot_product(ab1(a, b, :), c2%diag_bb) &
                      *yab(a, b)
            end do
         end do

         ! e6: D_p_AB(a,b,P) = sum_x yAB(a,x) BB1(x,b,P)
         e(6) = 0.0_dp
         do k = 1, c2%nd3
            do b = 1, no_b
               do a = 1, no_a
                  acc = 0.0_dp
                  do x = 1, no_b
                     acc = acc + yab(a, x)*bb1(x, b, k)
                  end do
                  e(6) = e(6) + 2.0_dp*ab1(a, b, k)*acc
               end do
            end do
         end do
         deallocate (yab)

         ! e7: C_p_BB(c,b,P) = sum_s t1(c,s) BS1(b,s,P) against
         !     D_p_BB(c,b,P) = sum_a s(a,b) AB1(a,c,P)
         e(7) = 0.0_dp
         do b = 1, no_b
            do c = 1, no_b
               acc = 0.0_dp
               do k = 1, c2%nd3
                  acc = acc + dot_product(t1(c, :), bs1(b, :, k)) &
                        *dot_product(s_oo(:, b), ab1(:, c, k))
               end do
               e(7) = e(7) + 2.0_dp*acc
            end do
         end do
      end associate

      energy = sum(e)
      deallocate (as2, ab1, aa1, bb1, bs1)
   end function exch102_k2f

   function exch120_k11u(c2, amps_a) result(energy)
      !! The six `exch120_k11u` pieces -- `exch12.cc:382-1786`, A side
      type(sapt2_cache_t), intent(in) :: c2
      type(sapt2_amps_t), intent(in) :: amps_a
      real(dp) :: energy

      real(dp), allocatable :: rb1(:, :, :), rb2(:, :, :), ar1(:, :, :)
      real(dp), allocatable :: ab1(:, :, :), ab2(:, :, :), aa1(:, :, :)
      real(dp), allocatable :: bb1(:, :, :), rr1(:, :, :)
      real(dp), allocatable :: d_rb(:, :, :), e_rb(:, :, :), d_br(:, :, :)
      real(dp), allocatable :: z_rb(:, :), xv(:), tt(:, :, :, :)
      real(dp), allocatable :: ths(:, :, :, :), ts(:, :, :, :), tp(:, :, :)
      real(dp), allocatable :: tq(:, :, :), m4(:, :, :, :), cpb(:, :, :)
      real(dp), allocatable :: p_vv(:, :), p_oo(:, :), t(:, :, :, :)
      real(dp), allocatable :: th(:, :, :, :), s_oo(:, :), s_vo(:, :)
      real(dp) :: e1, e2, e3, e4, e5, e6, acc, acc2
      integer :: no_a, no_b, nv_a, a, b, c, d, g, h, q, r, x, k

      no_a = c2%nocc_a
      no_b = c2%nocc_b
      nv_a = c2%nvir_a
      rb1 = dressed_rb(c2, 1)
      rb2 = dressed_rb(c2, 2)
      ar1 = dressed_ar(c2)
      ab1 = dressed_ab(c2, 1)
      ab2 = dressed_ab(c2, 2)
      aa1 = dressed_aa(c2)
      bb1 = dressed_bb(c2)
      rr1 = dressed_rr(c2)

      ! Locals rather than an associate block: gfortran 13 miscompiles the
      ! parse of a large associate with sliced selectors (ICE), and copies of
      ! these small arrays cost nothing.
      p_vv = amps_a%p_vv
      p_oo = amps_a%p_oo
      t = amps_a%t
      th = amps_a%theta
      s_oo = c2%s_ab(1:no_a, 1:no_b)
      s_vo = c2%s_ab(no_a + 1:c2%nmo_a, 1:no_b)

      ! -- k11u_1 (exch12.cc:382-718): the pRR density ------------------
      e1 = 0.0_dp
      do q = 1, nv_a
         do r = 1, nv_a
            acc = 0.0_dp
            do b = 1, no_b
               acc = acc + dot_product(rb1(r, b, :), rb2(q, b, :))
            end do
            e1 = e1 + 2.0_dp*p_vv(r, q)*acc
         end do
      end do
      allocate (d_rb(nv_a, no_b, c2%nd3), e_rb(nv_a, no_b, c2%nd3), &
                d_br(no_b, nv_a, c2%nd3))
      do k = 1, c2%nd3
         do b = 1, no_b
            do r = 1, nv_a
               d_rb(r, b, k) = dot_product(p_vv(r, :), rb1(:, b, k))
               e_rb(r, b, k) = dot_product(p_vv(r, :), rb2(:, b, k))
               d_br(b, r, k) = dot_product(s_oo(:, b), ar1(:, r, k))
            end do
         end do
      end do
      do b = 1, no_b
         do r = 1, nv_a
            e1 = e1 - 2.0_dp*dot_product(d_br(b, r, :), d_rb(r, b, :))
            e1 = e1 + 4.0_dp*s_vo(r, b)*dot_product(d_rb(r, b, :), &
                                                    c2%diag_aa)
            e1 = e1 + 4.0_dp*s_vo(r, b)*dot_product(e_rb(r, b, :), &
                                                    c2%diag_bb)
         end do
      end do
      do k = 1, c2%nd3
         do c = 1, no_b
            do b = 1, no_b
               e1 = e1 - 2.0_dp*bb1(b, c, k) &
                    *dot_product(s_vo(:, b), e_rb(:, c, k))
            end do
         end do
      end do
      allocate (z_rb(nv_a, no_b))
      do b = 1, no_b
         do r = 1, nv_a
            z_rb(r, b) = dot_product(p_vv(r, :), s_vo(:, b))
         end do
      end do
      do b = 1, no_b
         do r = 1, nv_a
            acc = 0.0_dp
            do a = 1, no_a
               acc = acc + dot_product(ar1(a, r, :), ab2(a, b, :))
            end do
            e1 = e1 - 2.0_dp*acc*z_rb(r, b)
            e1 = e1 - 8.0_dp*dot_product(d_br(b, r, :), c2%diag_bb) &
                 *z_rb(r, b)
         end do
      end do
      do k = 1, c2%nd3
         do c = 1, no_b
            do b = 1, no_b
               e1 = e1 + 4.0_dp*bb1(b, c, k) &
                    *dot_product(z_rb(:, c), d_br(b, :, k))
            end do
         end do
      end do
      do c = 1, no_b
         do b = 1, no_b
            e1 = e1 - 4.0_dp*dot_product(bb1(b, c, :), c2%diag_aa) &
                 *dot_product(s_vo(:, b), z_rb(:, c))
         end do
      end do
      deallocate (d_rb, e_rb, d_br)
      allocate (xv(c2%nd3))
      do k = 1, c2%nd3
         xv(k) = sum(rr1(:, :, k)*p_vv)
      end do
      do b = 1, no_b
         do a = 1, no_a
            e1 = e1 + 4.0_dp*s_oo(a, b)*dot_product(ab2(a, b, :), xv)
         end do
      end do
      do c = 1, no_b
         do b = 1, no_b
            e1 = e1 - 4.0_dp*dot_product(bb1(b, c, :), xv) &
                 *dot_product(s_oo(:, b), s_oo(:, c))
         end do
      end do
      deallocate (xv, z_rb)
      e1 = -e1

      ! -- k11u_2 (exch12.cc:720-876): the pAA density ------------------
      e2 = 0.0_dp
      do c = 1, no_a
         do a = 1, no_a
            acc = 0.0_dp
            do b = 1, no_b
               ! x1 + x3-with-sign + x4 + x6 pieces sharing the b loop
               acc = acc + dot_product(ab1(a, b, :), ab2(c, b, :))
            end do
            e2 = e2 + 2.0_dp*p_oo(a, c)*acc
         end do
      end do
      allocate (xv(c2%nd3))
      do k = 1, c2%nd3
         xv(k) = sum(ab2(:, :, k)*s_oo)
      end do
      do c = 1, no_a
         do a = 1, no_a
            e2 = e2 + 4.0_dp*p_oo(a, c)*dot_product(aa1(a, c, :), xv)
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            acc = 0.0_dp
            do b = 1, no_b
               acc2 = 0.0_dp
               do k = 1, c2%nd3
                  acc2 = acc2 + dot_product(s_oo(:, b), aa1(a, :, k)) &
                         *ab2(c, b, k)
               end do
               acc = acc + acc2
            end do
            e2 = e2 - 2.0_dp*p_oo(a, c)*acc
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            acc = 0.0_dp
            do b = 1, no_b
               acc = acc + s_oo(a, b)*dot_product(ab2(c, b, :), c2%diag_aa)
            end do
            e2 = e2 + 4.0_dp*p_oo(a, c)*acc
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            acc = 0.0_dp
            do b = 1, no_b
               acc2 = 0.0_dp
               do d = 1, no_a
                  acc2 = acc2 + dot_product(aa1(d, c, :), ab2(d, b, :))
               end do
               acc = acc + s_oo(a, b)*acc2
            end do
            e2 = e2 - 2.0_dp*p_oo(a, c)*acc
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            acc = 0.0_dp
            do b = 1, no_b
               acc = acc + dot_product(ab1(a, b, :), c2%diag_bb)*s_oo(c, b)
            end do
            e2 = e2 + 4.0_dp*p_oo(a, c)*acc
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            acc = 0.0_dp
            do b = 1, no_b
               acc2 = 0.0_dp
               do d = 1, no_b
                  acc2 = acc2 + dot_product(ab1(a, d, :), bb1(b, d, :))
               end do
               acc = acc + acc2*s_oo(c, b)
            end do
            e2 = e2 - 2.0_dp*p_oo(a, c)*acc
         end do
      end do
      do k = 1, c2%nd3
         xv(k) = 0.0_dp
         do d = 1, no_b
            do b = 1, no_b
               xv(k) = xv(k) + bb1(b, d, k) &
                       *dot_product(s_oo(:, b), s_oo(:, d))
            end do
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            e2 = e2 - 4.0_dp*p_oo(a, c)*dot_product(aa1(a, c, :), xv)
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            acc = 0.0_dp
            do b = 1, no_b
               acc2 = 0.0_dp
               do d = 1, no_b
                  acc2 = acc2 + s_oo(a, d) &
                         *dot_product(bb1(d, b, :), c2%diag_aa)
               end do
               acc = acc + acc2*s_oo(c, b)
            end do
            e2 = e2 - 4.0_dp*p_oo(a, c)*acc
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            acc = 0.0_dp
            do d = 1, no_a
               acc = acc + dot_product(s_oo(a, :), s_oo(d, :)) &
                     *dot_product(aa1(c, d, :), c2%diag_bb)
            end do
            e2 = e2 - 8.0_dp*p_oo(a, c)*acc
         end do
      end do
      ! x11: C_p_aB(a,b,P) = sum_d s(a,d) BB1(d,b,P); C_p_aA(a,d,P) =
      ! sum_b s(d,b) C_p_aB(a,b,P); contract with AA1(c,d,P) and paa.
      allocate (cpb(no_a, no_b, c2%nd3))
      do k = 1, c2%nd3
         do b = 1, no_b
            do a = 1, no_a
               cpb(a, b, k) = dot_product(s_oo(a, :), bb1(:, b, k))
            end do
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            acc = 0.0_dp
            do d = 1, no_a
               do k = 1, c2%nd3
                  acc = acc + aa1(c, d, k) &
                        *dot_product(s_oo(d, :), cpb(a, :, k))
               end do
            end do
            e2 = e2 + 4.0_dp*p_oo(a, c)*acc
         end do
      end do
      deallocate (xv, cpb)

      ! -- k11u_3 (exch12.cc:1042-1170): the theta.t squares ------------
      allocate (ths(no_a, nv_a, no_a, no_b), ts(no_a, nv_a, no_a, no_b))
      do b = 1, no_b
         do g = 1, no_a
            do r = 1, nv_a
               do c = 1, no_a
                  ths(c, r, g, b) = dot_product(th(c, r, g, :), s_vo(:, b))
                  ts(c, r, g, b) = dot_product(t(c, r, g, :), s_vo(:, b))
               end do
            end do
         end do
      end do
      e3 = 0.0_dp
      do b = 1, no_b
         do q = 1, nv_a
            do x = 1, nv_a
               do r = 1, nv_a
                  acc = dot_product(rr1(r, x, 1:c2%nd3), rb1(q, b, 1:c2%nd3))
                  acc2 = 0.0_dp
                  do g = 1, no_a
                     do c = 1, no_a
                        acc2 = acc2 + t(c, x, g, q)*ths(c, r, g, b)
                     end do
                  end do
                  e3 = e3 + 2.0_dp*acc*acc2
               end do
            end do
         end do
      end do
      do q = 1, nv_a
         do r = 1, nv_a
            e3 = e3 + 4.0_dp*sum(ts(:, r, :, :)*ths(:, q, :, :)) &
                 *dot_product(rr1(r, q, :), c2%diag_bb)
         end do
      end do
      do d = 1, no_b
         do b = 1, no_b
            do x = 1, nv_a
               do r = 1, nv_a
                  acc = dot_product(rr1(r, x, :), bb1(b, d, :))
                  acc2 = 0.0_dp
                  do g = 1, no_a
                     do c = 1, no_a
                        acc2 = acc2 + ts(c, x, g, b)*ths(c, r, g, d)
                     end do
                  end do
                  e3 = e3 - 2.0_dp*acc*acc2
               end do
            end do
         end do
      end do
      e3 = -e3

      ! -- k11u_4 (exch12.cc:1332-1424): the occupied square ------------
      allocate (m4(no_a, no_a, no_a, no_a))
      do h = 1, no_a
         do c = 1, no_a
            do g = 1, no_a
               do a = 1, no_a
                  m4(a, g, c, h) = sum(th(a, :, g, :)*t(c, :, h, :))
               end do
            end do
         end do
      end do
      allocate (tp(no_a, no_a, c2%nd3))
      do k = 1, c2%nd3
         do c = 1, no_a
            do a = 1, no_a
               tp(a, c, k) = sum(m4(a, :, c, :)*aa1(:, :, k))
            end do
         end do
      end do
      deallocate (m4)
      e4 = 0.0_dp
      do c = 1, no_a
         do a = 1, no_a
            acc = 0.0_dp
            do k = 1, c2%nd3
               acc = acc + tp(a, c, k)*dot_product(s_oo(c, :), ab2(a, :, k))
            end do
            e4 = e4 + 2.0_dp*acc
            e4 = e4 + 4.0_dp*dot_product(tp(a, c, :), c2%diag_bb) &
                 *dot_product(s_oo(a, :), s_oo(c, :))
            acc = 0.0_dp
            do k = 1, c2%nd3
               acc2 = 0.0_dp
               do b = 1, no_b
                  acc2 = acc2 + s_oo(c, b) &
                         *dot_product(s_oo(a, :), bb1(:, b, k))
               end do
               acc = acc + tp(a, c, k)*acc2
            end do
            e4 = e4 - 2.0_dp*acc
         end do
      end do
      deallocate (tp)
      e4 = -e4

      ! -- k11u_5 (exch12.cc:1471-1570): theta against Theta^P ----------
      allocate (tq(no_a, nv_a, c2%nd3))
      do k = 1, c2%nd3
         do r = 1, nv_a
            do a = 1, no_a
               tq(a, r, k) = sum(th(a, r, :, :)*amps_a%th_ov(:, :, k))
            end do
         end do
      end do
      e5 = 0.0_dp
      do b = 1, no_b
         do r = 1, nv_a
            acc = 0.0_dp
            do k = 1, c2%nd3
               acc = acc + dot_product(s_oo(:, b), tq(:, r, k))*rb1(r, b, k)
            end do
            e5 = e5 + acc
         end do
      end do
      do k = 1, c2%nd3
         do b = 1, no_b
            do a = 1, no_a
               acc = dot_product(s_vo(:, b), tq(a, :, k))
               e5 = e5 + ab2(a, b, k)*acc
            end do
         end do
      end do
      do k = 1, c2%nd3
         do c = 1, no_b
            do b = 1, no_b
               acc = 0.0_dp
               do a = 1, no_a
                  acc = acc + s_oo(a, b)*dot_product(s_vo(:, c), &
                                                     tq(a, :, k))
               end do
               e5 = e5 - 2.0_dp*bb1(b, c, k)*acc
            end do
         end do
      end do
      do r = 1, nv_a
         do a = 1, no_a
            e5 = e5 + 4.0_dp*dot_product(s_oo(a, :), s_vo(r, :)) &
                 *dot_product(tq(a, r, :), c2%diag_bb)
         end do
      end do
      deallocate (tq)
      e5 = -2.0_dp*e5

      ! -- k11u_6 (exch12.cc:1662-1786): the 3 t.t + theta.theta square -
      allocate (tt(no_a, nv_a, no_a, nv_a))
      do h = 1, nv_a
         do g = 1, no_a
            do r = 1, nv_a
               do a = 1, no_a
                  tt(a, r, g, h) = 3.0_dp*sum(t(a, r, :, :)*t(g, h, :, :)) &
                                   + sum(th(:, r, a, :)*th(:, h, g, :))
               end do
            end do
         end do
      end do
      allocate (tp(no_a, no_a, c2%nd3), tq(nv_a, nv_a, c2%nd3))
      do k = 1, c2%nd3
         do g = 1, no_a
            do a = 1, no_a
               tp(a, g, k) = sum(tt(a, :, g, :)*rr1(:, :, k))
            end do
         end do
         do h = 1, nv_a
            do r = 1, nv_a
               tq(r, h, k) = sum(tt(:, r, :, h)*aa1(:, :, k))
            end do
         end do
      end do
      deallocate (tt)
      e6 = 0.0_dp
      do k = 1, c2%nd3
         do b = 1, no_b
            do g = 1, no_a
               acc = dot_product(s_oo(:, b), tp(:, g, k))
               e6 = e6 - acc*ab2(g, b, k)
            end do
         end do
         do d = 1, no_b
            do b = 1, no_b
               acc = 0.0_dp
               do g = 1, no_a
                  acc = acc + s_oo(g, d)*dot_product(s_oo(:, b), &
                                                     tp(:, g, k))
               end do
               e6 = e6 + bb1(b, d, k)*acc
            end do
         end do
         do r = 1, nv_a
            do b = 1, no_b
               acc = dot_product(s_vo(:, b), tq(:, r, k))
               e6 = e6 - rb1(r, b, k)*acc
            end do
         end do
         do d = 1, no_b
            do b = 1, no_b
               acc = 0.0_dp
               do r = 1, nv_a
                  acc = acc + s_vo(r, d)*dot_product(s_vo(:, b), &
                                                     tq(:, r, k))
               end do
               e6 = e6 + bb1(b, d, k)*acc
            end do
         end do
      end do
      do g = 1, no_a
         do a = 1, no_a
            e6 = e6 - 2.0_dp*dot_product(tp(a, g, :), c2%diag_bb) &
                 *dot_product(s_oo(a, :), s_oo(g, :))
         end do
      end do
      do h = 1, nv_a
         do r = 1, nv_a
            e6 = e6 - 2.0_dp*dot_product(tq(r, h, :), c2%diag_bb) &
                 *dot_product(s_vo(r, :), s_vo(h, :))
         end do
      end do
      deallocate (tp, tq, ths, ts)
      e6 = -e6

      energy = e1 + e2 + e3 + e4 + e5 + e6
      deallocate (rb1, rb2, ar1, ab1, ab2, aa1, bb1, rr1)
   end function exch120_k11u

   function exch102_k11u(c2, amps_b) result(energy)
      !! The six `exch102_k11u` pieces -- `exch12.cc:562-1900`, B side
      type(sapt2_cache_t), intent(in) :: c2
      type(sapt2_amps_t), intent(in) :: amps_b
      real(dp) :: energy

      real(dp), allocatable :: as1(:, :, :), as2(:, :, :), bs1(:, :, :)
      real(dp), allocatable :: ab1(:, :, :), ab2(:, :, :), aa1(:, :, :)
      real(dp), allocatable :: bb1(:, :, :), ss1(:, :, :)
      real(dp), allocatable :: d_as(:, :, :), e_as(:, :, :), f_as(:, :, :)
      real(dp), allocatable :: z_as(:, :), xv(:), tt(:, :, :, :)
      real(dp), allocatable :: ths(:, :, :, :), ts(:, :, :, :), tp(:, :, :)
      real(dp), allocatable :: tq(:, :, :), m4(:, :, :, :), cpb(:, :, :)
      real(dp), allocatable :: p_ss(:, :), p_bb(:, :), t(:, :, :, :)
      real(dp), allocatable :: th(:, :, :, :), s_oo(:, :), s_ov(:, :)
      real(dp) :: e1, e2, e3, e4, e5, e6, acc, acc2
      integer :: no_a, no_b, nv_b, a, b, c, d, g, h, s, t_, x, k

      no_a = c2%nocc_a
      no_b = c2%nocc_b
      nv_b = c2%nvir_b
      as1 = dressed_as(c2, 1)
      as2 = dressed_as(c2, 2)
      bs1 = dressed_bs(c2)
      ab1 = dressed_ab(c2, 1)
      ab2 = dressed_ab(c2, 2)
      aa1 = dressed_aa(c2)
      bb1 = dressed_bb(c2)
      ss1 = dressed_ss(c2)

      ! Locals rather than an associate block, as in exch120_k11u.
      p_ss = amps_b%p_vv
      p_bb = amps_b%p_oo
      t = amps_b%t
      th = amps_b%theta
      s_oo = c2%s_ab(1:no_a, 1:no_b)
      s_ov = c2%s_ab(1:no_a, c2%nocc_b + 1:c2%nmo_b)

      ! -- k11u_1 (exch12.cc:562-718): the pSS density ------------------
      e1 = 0.0_dp
      do t_ = 1, nv_b
         do s = 1, nv_b
            acc = 0.0_dp
            do a = 1, no_a
               acc = acc + dot_product(as1(a, s, :), as2(a, t_, :))
            end do
            e1 = e1 + 2.0_dp*p_ss(s, t_)*acc
         end do
      end do
      allocate (d_as(no_a, nv_b, c2%nd3), e_as(no_a, nv_b, c2%nd3), &
                f_as(no_a, nv_b, c2%nd3))
      do k = 1, c2%nd3
         do s = 1, nv_b
            do a = 1, no_a
               d_as(a, s, k) = dot_product(p_ss(s, :), as1(a, :, k))
               e_as(a, s, k) = dot_product(p_ss(s, :), as2(a, :, k))
               f_as(a, s, k) = dot_product(s_oo(a, :), bs1(:, s, k))
            end do
         end do
      end do
      do s = 1, nv_b
         do a = 1, no_a
            e1 = e1 - 2.0_dp*dot_product(d_as(a, s, :), f_as(a, s, :))
            e1 = e1 + 4.0_dp*s_ov(a, s)*dot_product(d_as(a, s, :), &
                                                    c2%diag_bb)
            e1 = e1 + 4.0_dp*s_ov(a, s)*dot_product(e_as(a, s, :), &
                                                    c2%diag_aa)
         end do
      end do
      do k = 1, c2%nd3
         do c = 1, no_a
            do a = 1, no_a
               e1 = e1 - 2.0_dp*aa1(a, c, k) &
                    *dot_product(s_ov(c, :), e_as(a, :, k))
            end do
         end do
      end do
      allocate (z_as(no_a, nv_b))
      do s = 1, nv_b
         do a = 1, no_a
            z_as(a, s) = dot_product(s_ov(a, :), p_ss(:, s))
         end do
      end do
      do s = 1, nv_b
         do a = 1, no_a
            acc = 0.0_dp
            do b = 1, no_b
               acc = acc + dot_product(ab1(a, b, :), bs1(b, s, :))
            end do
            e1 = e1 - 2.0_dp*acc*z_as(a, s)
            e1 = e1 - 8.0_dp*dot_product(f_as(a, s, :), c2%diag_aa) &
                 *z_as(a, s)
         end do
      end do
      do k = 1, c2%nd3
         do c = 1, no_a
            do a = 1, no_a
               e1 = e1 + 4.0_dp*aa1(a, c, k) &
                    *dot_product(z_as(c, :), f_as(a, :, k))
            end do
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            e1 = e1 - 4.0_dp*dot_product(aa1(a, c, :), c2%diag_bb) &
                 *dot_product(z_as(a, :), s_ov(c, :))
         end do
      end do
      deallocate (d_as, e_as, f_as)
      allocate (xv(c2%nd3))
      do k = 1, c2%nd3
         xv(k) = sum(ss1(:, :, k)*p_ss)
      end do
      do b = 1, no_b
         do a = 1, no_a
            e1 = e1 + 4.0_dp*s_oo(a, b)*dot_product(ab1(a, b, :), xv)
         end do
      end do
      do c = 1, no_a
         do a = 1, no_a
            e1 = e1 - 4.0_dp*dot_product(aa1(a, c, :), xv) &
                 *dot_product(s_oo(a, :), s_oo(c, :))
         end do
      end do
      deallocate (xv, z_as)
      e1 = -e1

      ! -- k11u_2 (exch12.cc:878-1040): the pBB density -----------------
      e2 = 0.0_dp
      do d = 1, no_b
         do b = 1, no_b
            acc = 0.0_dp
            do a = 1, no_a
               acc = acc + dot_product(ab2(a, b, :), ab1(a, d, :))
            end do
            e2 = e2 + 2.0_dp*p_bb(b, d)*acc
         end do
      end do
      allocate (xv(c2%nd3))
      do k = 1, c2%nd3
         xv(k) = sum(ab1(:, :, k)*s_oo)
      end do
      do d = 1, no_b
         do b = 1, no_b
            e2 = e2 + 4.0_dp*p_bb(b, d)*dot_product(bb1(b, d, :), xv)
         end do
      end do
      do d = 1, no_b
         do b = 1, no_b
            acc = 0.0_dp
            do a = 1, no_a
               acc2 = 0.0_dp
               do k = 1, c2%nd3
                  acc2 = acc2 + dot_product(s_oo(a, :), bb1(b, :, k)) &
                         *ab1(a, d, k)
               end do
               acc = acc + acc2
            end do
            e2 = e2 - 2.0_dp*p_bb(b, d)*acc
         end do
      end do
      do d = 1, no_b
         do b = 1, no_b
            acc = 0.0_dp
            do a = 1, no_a
               acc = acc + s_oo(a, b)*dot_product(ab1(a, d, :), c2%diag_bb)
            end do
            e2 = e2 + 4.0_dp*p_bb(b, d)*acc
         end do
      end do
      do d = 1, no_b
         do b = 1, no_b
            acc = 0.0_dp
            do a = 1, no_a
               acc2 = 0.0_dp
               do c = 1, no_b
                  acc2 = acc2 + dot_product(ab1(a, c, :), bb1(c, d, :))
               end do
               acc = acc + s_oo(a, b)*acc2
            end do
            e2 = e2 - 2.0_dp*p_bb(b, d)*acc
         end do
      end do
      do d = 1, no_b
         do b = 1, no_b
            acc = 0.0_dp
            do a = 1, no_a
               acc = acc + dot_product(ab2(a, b, :), c2%diag_aa)*s_oo(a, d)
            end do
            e2 = e2 + 4.0_dp*p_bb(b, d)*acc
         end do
      end do
      do d = 1, no_b
         do b = 1, no_b
            acc = 0.0_dp
            do a = 1, no_a
               acc2 = 0.0_dp
               do c = 1, no_a
                  acc2 = acc2 + dot_product(aa1(c, a, :), ab2(c, b, :))
               end do
               acc = acc + acc2*s_oo(a, d)
            end do
            e2 = e2 - 2.0_dp*p_bb(b, d)*acc
         end do
      end do
      do k = 1, c2%nd3
         xv(k) = 0.0_dp
         do c = 1, no_a
            do a = 1, no_a
               xv(k) = xv(k) + aa1(a, c, k) &
                       *dot_product(s_oo(a, :), s_oo(c, :))
            end do
         end do
      end do
      do d = 1, no_b
         do b = 1, no_b
            e2 = e2 - 4.0_dp*p_bb(b, d)*dot_product(bb1(b, d, :), xv)
         end do
      end do
      do d = 1, no_b
         do b = 1, no_b
            acc = 0.0_dp
            do a = 1, no_a
               acc2 = 0.0_dp
               do c = 1, no_a
                  acc2 = acc2 + dot_product(aa1(a, c, :), c2%diag_bb) &
                         *s_oo(c, b)
               end do
               acc = acc + acc2*s_oo(a, d)
            end do
            e2 = e2 - 4.0_dp*p_bb(b, d)*acc
         end do
      end do
      do d = 1, no_b
         do b = 1, no_b
            acc = 0.0_dp
            do c = 1, no_b
               acc = acc + dot_product(s_oo(:, b), s_oo(:, c)) &
                     *dot_product(bb1(d, c, :), c2%diag_aa)
            end do
            e2 = e2 - 8.0_dp*p_bb(b, d)*acc
         end do
      end do
      ! x11: C_p_Ab(b,a,P) = sum_c s(c,b) AA1(c,a,P); C_p_bB(b,d,P) =
      ! sum_a s(a,d) C_p_Ab(b,a,P); contract with BB1(c,d,P) and pbb.
      allocate (cpb(no_b, no_a, c2%nd3))
      do k = 1, c2%nd3
         do a = 1, no_a
            do b = 1, no_b
               cpb(b, a, k) = dot_product(s_oo(:, b), aa1(:, a, k))
            end do
         end do
      end do
      do c = 1, no_b
         do b = 1, no_b
            acc = 0.0_dp
            do d = 1, no_b
               do k = 1, c2%nd3
                  acc = acc + bb1(c, d, k) &
                        *dot_product(s_oo(:, d), cpb(b, :, k))
               end do
            end do
            e2 = e2 + 4.0_dp*p_bb(b, c)*acc
         end do
      end do
      deallocate (xv, cpb)

      ! -- k11u_3 (exch12.cc:1172-1330) ---------------------------------
      allocate (ths(no_b, nv_b, no_b, no_a), ts(no_b, nv_b, no_b, no_a))
      do a = 1, no_a
         do g = 1, no_b
            do s = 1, nv_b
               do c = 1, no_b
                  ths(c, s, g, a) = dot_product(th(c, s, g, :), s_ov(a, :))
                  ts(c, s, g, a) = dot_product(t(c, s, g, :), s_ov(a, :))
               end do
            end do
         end do
      end do
      e3 = 0.0_dp
      do a = 1, no_a
         do t_ = 1, nv_b
            do x = 1, nv_b
               do s = 1, nv_b
                  acc = dot_product(ss1(s, x, :), as1(a, t_, :))
                  acc2 = 0.0_dp
                  do g = 1, no_b
                     do c = 1, no_b
                        acc2 = acc2 + t(c, x, g, t_)*ths(c, s, g, a)
                     end do
                  end do
                  e3 = e3 + 2.0_dp*acc*acc2
               end do
            end do
         end do
      end do
      do t_ = 1, nv_b
         do s = 1, nv_b
            e3 = e3 + 4.0_dp*sum(ts(:, s, :, :)*ths(:, t_, :, :)) &
                 *dot_product(ss1(s, t_, :), c2%diag_aa)
         end do
      end do
      do d = 1, no_a
         do a = 1, no_a
            do x = 1, nv_b
               do s = 1, nv_b
                  acc = dot_product(ss1(s, x, :), aa1(a, d, :))
                  acc2 = 0.0_dp
                  do g = 1, no_b
                     do c = 1, no_b
                        acc2 = acc2 + ts(c, x, g, a)*ths(c, s, g, d)
                     end do
                  end do
                  e3 = e3 - 2.0_dp*acc*acc2
               end do
            end do
         end do
      end do
      e3 = -e3

      ! -- k11u_4 (exch12.cc:1426-1500) ---------------------------------
      allocate (m4(no_b, no_b, no_b, no_b))
      do h = 1, no_b
         do c = 1, no_b
            do g = 1, no_b
               do b = 1, no_b
                  m4(b, g, c, h) = sum(th(b, :, g, :)*t(c, :, h, :))
               end do
            end do
         end do
      end do
      allocate (tp(no_b, no_b, c2%nd3))
      do k = 1, c2%nd3
         do c = 1, no_b
            do b = 1, no_b
               tp(b, c, k) = sum(m4(b, :, c, :)*bb1(:, :, k))
            end do
         end do
      end do
      deallocate (m4)
      e4 = 0.0_dp
      do c = 1, no_b
         do b = 1, no_b
            acc = 0.0_dp
            do k = 1, c2%nd3
               acc = acc + tp(b, c, k)*dot_product(s_oo(:, b), ab1(:, c, k))
            end do
            e4 = e4 + 2.0_dp*acc
            e4 = e4 + 4.0_dp*dot_product(tp(b, c, :), c2%diag_aa) &
                 *dot_product(s_oo(:, b), s_oo(:, c))
            acc = 0.0_dp
            do k = 1, c2%nd3
               acc2 = 0.0_dp
               do a = 1, no_a
                  acc2 = acc2 + s_oo(a, c) &
                         *dot_product(s_oo(:, b), aa1(a, :, k))
               end do
               acc = acc + tp(b, c, k)*acc2
            end do
            e4 = e4 - 2.0_dp*acc
         end do
      end do
      deallocate (tp)
      e4 = -e4

      ! -- k11u_5 (exch12.cc:1572-1660) ---------------------------------
      allocate (tq(no_b, nv_b, c2%nd3))
      do k = 1, c2%nd3
         do s = 1, nv_b
            do b = 1, no_b
               tq(b, s, k) = sum(th(b, s, :, :)*amps_b%th_ov(:, :, k))
            end do
         end do
      end do
      e5 = 0.0_dp
      do k = 1, c2%nd3
         do s = 1, nv_b
            do a = 1, no_a
               acc = dot_product(s_oo(a, :), tq(:, s, k))
               e5 = e5 + as1(a, s, k)*acc
            end do
         end do
         do b = 1, no_b
            do a = 1, no_a
               acc = dot_product(s_ov(a, :), tq(b, :, k))
               e5 = e5 + ab1(a, b, k)*acc
            end do
         end do
         do c = 1, no_a
            do a = 1, no_a
               acc = 0.0_dp
               do b = 1, no_b
                  acc = acc + s_oo(a, b)*dot_product(s_ov(c, :), &
                                                     tq(b, :, k))
               end do
               e5 = e5 - 2.0_dp*aa1(a, c, k)*acc
            end do
         end do
      end do
      do s = 1, nv_b
         do b = 1, no_b
            e5 = e5 + 4.0_dp*dot_product(s_oo(:, b), s_ov(:, s)) &
                 *dot_product(tq(b, s, :), c2%diag_aa)
         end do
      end do
      deallocate (tq)
      e5 = -2.0_dp*e5

      ! -- k11u_6 (exch12.cc:1788-1900) ---------------------------------
      allocate (tt(no_b, nv_b, no_b, nv_b))
      do h = 1, nv_b
         do g = 1, no_b
            do s = 1, nv_b
               do b = 1, no_b
                  tt(b, s, g, h) = 3.0_dp*sum(t(b, s, :, :)*t(g, h, :, :)) &
                                   + sum(th(:, s, b, :)*th(:, h, g, :))
               end do
            end do
         end do
      end do
      allocate (tp(no_b, no_b, c2%nd3), tq(nv_b, nv_b, c2%nd3))
      do k = 1, c2%nd3
         do g = 1, no_b
            do b = 1, no_b
               tp(b, g, k) = sum(tt(b, :, g, :)*ss1(:, :, k))
            end do
         end do
         do h = 1, nv_b
            do s = 1, nv_b
               tq(s, h, k) = sum(tt(:, s, :, h)*bb1(:, :, k))
            end do
         end do
      end do
      deallocate (tt)
      e6 = 0.0_dp
      do k = 1, c2%nd3
         do g = 1, no_b
            do a = 1, no_a
               acc = dot_product(s_oo(a, :), tp(:, g, k))
               e6 = e6 - acc*ab1(a, g, k)
            end do
         end do
         do c = 1, no_a
            do a = 1, no_a
               acc = 0.0_dp
               do g = 1, no_b
                  acc = acc + s_oo(c, g)*dot_product(s_oo(a, :), &
                                                     tp(:, g, k))
               end do
               e6 = e6 + aa1(a, c, k)*acc
            end do
         end do
         do s = 1, nv_b
            do a = 1, no_a
               acc = dot_product(s_ov(a, :), tq(:, s, k))
               e6 = e6 - as1(a, s, k)*acc
            end do
         end do
         do c = 1, no_a
            do a = 1, no_a
               acc = 0.0_dp
               do s = 1, nv_b
                  acc = acc + s_ov(c, s)*dot_product(s_ov(a, :), &
                                                     tq(:, s, k))
               end do
               e6 = e6 + aa1(a, c, k)*acc
            end do
         end do
      end do
      do g = 1, no_b
         do b = 1, no_b
            e6 = e6 - 2.0_dp*dot_product(tp(b, g, :), c2%diag_aa) &
                 *dot_product(s_oo(:, b), s_oo(:, g))
         end do
      end do
      do h = 1, nv_b
         do s = 1, nv_b
            e6 = e6 - 2.0_dp*dot_product(tq(s, h, :), c2%diag_aa) &
                 *dot_product(s_ov(:, s), s_ov(:, h))
         end do
      end do
      deallocate (tp, tq, ths, ts)
      e6 = -e6

      energy = e1 + e2 + e3 + e4 + e5 + e6
      deallocate (as1, as2, bs1, ab1, ab2, aa1, bb1, ss1)
   end function exch102_k11u

   function sapt_exch12(c2, amps_a, amps_b, u_ar, u_bs) result(energy)
      !! `Exch12` -- `exch12.cc:37-97`: `exch111`, the K2u pieces (the
      !! `Exch11` machinery driven by the Theta 2 intermediates), K2f, and
      !! the six K11u pieces per side
      type(sapt2_cache_t), intent(in) :: c2
      type(sapt2_amps_t), intent(in) :: amps_a, amps_b
      real(dp), intent(in) :: u_ar(:, :), u_bs(:, :)
      real(dp) :: energy

      energy = exch111(c2, amps_a%th_ov, amps_b%th_ov) &
               + exch110(c2, amps_a%th2_ov) + exch101(c2, amps_b%th2_ov) &
               + exch120_k2f(c2, amps_a, u_ar) + exch102_k2f(c2, amps_b, u_bs) &
               + exch120_k11u(c2, amps_a) + exch102_k11u(c2, amps_b)
   end function sapt_exch12

   subroutine sapt_ind22(c2, amps_a, amps_b, energy, pieces)
      !! `Ind22 = Ind220 + Ind202` -- `ind22.cc:37-148`
      !!
      !! The response is UNCOUPLED throughout -- `i = w/(e_o - e_v)` -- and
      !! the seven pieces per side are MBPT contractions of it with the
      !! monomer's amplitudes; no CPHF solve appears in psi4's Ind22.
      !!
      !! Two things the plan did not anticipate, found against psi4:
      !! piece 5 consumes the SECOND-order doubles (`ind22.cc:88` passes
      !! "t2ARAR Amplitudes"), the one place `t2` enters Ind22; and piece 1
      !! builds its perturbed doubles from bare integrals times the response,
      !! so it does NOT vanish when the stored amplitudes are zeroed --
      !! the `t2 -> 0` identity covers pieces 2-7 only.
      type(sapt2_cache_t), intent(in) :: c2
      type(sapt2_amps_t), intent(in) :: amps_a, amps_b
      real(dp), intent(out) :: energy
      real(dp), intent(out), optional :: pieces(SAPT2_PIECES, 2)

      real(dp) :: pa(SAPT2_PIECES), pb(SAPT2_PIECES)

      call ind220_one(c2, amps_a, "A", pa)
      call ind220_one(c2, amps_b, "B", pb)
      energy = sum(pa) + sum(pb)
      if (present(pieces)) then
         pieces(:, 1) = pa
         pieces(:, 2) = pb
      end if
   end subroutine sapt_ind22

   subroutine ind220_one(c2, amps, side, p)
      !! One direction of `Ind22` -- `ind220`, `ind22.cc:58-460`
      type(sapt2_cache_t), intent(in) :: c2
      type(sapt2_amps_t), intent(in) :: amps
      character(len=1), intent(in) :: side
      real(dp), intent(out) :: p(SAPT2_PIECES)

      real(dp), allocatable :: i_ov(:, :), i_other(:, :)
      real(dp), allocatable :: w_oo(:, :), w_ov(:, :), w_vv(:, :)
      real(dp), allocatable :: b_oo(:, :, :), b_ov(:, :, :), b_vv(:, :, :)
      real(dp), allocatable :: b_ov_o(:, :, :), eps(:)
      real(dp), allocatable :: x1(:, :, :, :), y1(:, :, :, :), g(:, :, :, :)
      real(dp), allocatable :: xw(:), yw(:), zw(:), ww(:)
      real(dp), allocatable :: xar(:, :), yar(:, :)
      integer :: no, nv, no_o, nv_o, a, r, c, q, x, k
      real(dp) :: denom, acc

      if (side == "A") then
         no = c2%nocc_a
         nv = c2%nvir_a
         w_oo = c2%w_baa
         w_ov = c2%w_bar
         w_vv = c2%w_brr
         b_oo = c2%b_aa
         b_ov = c2%b_ar
         b_vv = c2%b_rr
         eps = c2%eps_a
         no_o = c2%nocc_b
         nv_o = c2%nvir_b
         b_ov_o = c2%b_bs
         allocate (i_other(no_o, nv_o))
         do q = 1, nv_o
            do c = 1, no_o
               i_other(c, q) = c2%w_abs(c, q)/(c2%eps_b(c) - c2%eps_b(no_o + q))
            end do
         end do
      else
         no = c2%nocc_b
         nv = c2%nvir_b
         w_oo = c2%w_abb
         w_ov = c2%w_abs
         w_vv = c2%w_ass
         b_oo = c2%b_bb
         b_ov = c2%b_bs
         b_vv = c2%b_ss
         eps = c2%eps_b
         no_o = c2%nocc_a
         nv_o = c2%nvir_a
         b_ov_o = c2%b_ar
         allocate (i_other(no_o, nv_o))
         do q = 1, nv_o
            do c = 1, no_o
               i_other(c, q) = c2%w_bar(c, q)/(c2%eps_a(c) - c2%eps_a(no_o + q))
            end do
         end do
      end if

      allocate (i_ov(no, nv))
      do r = 1, nv
         do a = 1, no
            i_ov(a, r) = w_ov(a, r)/(eps(a) - eps(no + r))
         end do
      end do

      ! ind220_1 (ind22.cc:150-229): the wB-perturbed first-order doubles.
      ! Their RHS mixes bare integrals times the response with w times the
      ! stored t -- which is why this piece alone survives t -> 0.
      allocate (x1(no, nv, no, nv))
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  acc = 0.0_dp
                  do x = 1, nv
                     acc = acc + i_ov(a, x)*dot_product(b_vv(x, r, :), &
                                                        b_ov(c, q, :))
                  end do
                  do x = 1, no
                     acc = acc - i_ov(x, r)*dot_product(b_oo(a, x, :), &
                                                        b_ov(c, q, :))
                     acc = acc - w_oo(a, x)*amps%t(x, r, c, q)
                  end do
                  do x = 1, nv
                     acc = acc + amps%t(a, r, c, x)*w_vv(q, x)
                  end do
                  x1(a, r, c, q) = acc
               end do
            end do
         end do
      end do
      allocate (y1(no, nv, no, nv))
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  y1(a, r, c, q) = x1(a, r, c, q) + x1(c, q, a, r)
               end do
            end do
         end do
      end do
      x1 = y1
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  y1(a, r, c, q) = 2.0_dp*x1(a, r, c, q) - x1(c, r, a, q)
               end do
            end do
         end do
      end do
      p(1) = 0.0_dp
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  denom = eps(a) + eps(c) - eps(no + r) - eps(no + q)
                  p(1) = p(1) + x1(a, r, c, q)*y1(a, r, c, q)/denom
               end do
            end do
         end do
      end do
      deallocate (x1, y1)

      ! ind220_2 (:231-255)
      p(2) = 0.0_dp
      do r = 1, nv
         do a = 1, no
            acc = dot_product(i_ov(a, :), w_vv(r, :)) &
                  - dot_product(w_oo(a, :), i_ov(:, r))
            p(2) = p(2) + 4.0_dp*amps%t_singles(a, r)*acc
         end do
      end do

      ! ind220_3 (:257-296)
      p(3) = 0.0_dp
      do c = 1, no
         do a = 1, no
            p(3) = p(3) - 2.0_dp*amps%p_oo(a, c) &
                   *dot_product(i_ov(a, :), w_ov(c, :))
         end do
      end do
      do q = 1, nv
         do r = 1, nv
            p(3) = p(3) - 2.0_dp*amps%p_vv(r, q) &
                   *dot_product(i_ov(:, r), w_ov(:, q))
         end do
      end do

      ! ind220_4 (:298-334)
      p(4) = 0.0_dp
      do k = 1, size(b_ov, 3)
         do r = 1, nv
            do a = 1, no
               acc = 0.0_dp
               do x = 1, no
                  acc = acc + dot_product(i_ov(a, :), i_ov(x, :))*b_ov(x, r, k)
               end do
               do x = 1, nv
                  acc = acc + dot_product(i_ov(:, r), i_ov(:, x))*b_ov(a, x, k)
               end do
               p(4) = p(4) - 2.0_dp*acc*amps%th_ov(a, r, k)
            end do
         end do
      end do

      ! ind220_5 (:336-366): antisym(t2) with the denominator restored --
      ! the SECOND-order doubles, the one place Ind22 consumes t2.
      p(5) = 0.0_dp
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  denom = eps(a) + eps(c) - eps(no + r) - eps(no + q)
                  p(5) = p(5) + 2.0_dp*(2.0_dp*amps%t2(a, r, c, q) &
                                        - amps%t2(c, r, a, q))*denom &
                         *i_ov(a, r)*i_ov(c, q)
               end do
            end do
         end do
      end do

      ! ind220_6 (:368-421): g = 2(ar|cq) - (ac|rq)
      allocate (g(no, nv, no, nv))
      do q = 1, nv
         do c = 1, no
            do r = 1, nv
               do a = 1, no
                  g(a, r, c, q) = 2.0_dp*dot_product(b_ov(a, r, :), &
                                                     b_ov(c, q, :)) &
                                  - dot_product(b_oo(a, c, :), b_vv(r, q, :))
               end do
            end do
         end do
      end do
      allocate (xar(no, nv), yar(no, nv))
      do r = 1, nv
         do a = 1, no
            xar(a, r) = sum(g(a, r, :, :)*i_ov)
            yar(a, r) = sum(amps%theta(a, r, :, :)*i_ov)
         end do
      end do
      p(6) = -4.0_dp*sum(xar*yar)
      deallocate (g, xar, yar)

      ! ind220_7 (:423-460): this monomer's density against the OTHER
      ! monomer's uncoupled response, raw blocks -- pure Coulomb coupling.
      allocate (ww(size(b_oo, 3)), xw(size(b_oo, 3)), yw(size(b_oo, 3)), &
                zw(size(b_oo, 3)))
      do k = 1, size(b_oo, 3)
         ww(k) = sum(b_oo(:, :, k)*amps%p_oo)
         xw(k) = sum(b_vv(:, :, k)*amps%p_vv)
         yw(k) = sum(b_ov(:, :, k)*amps%t_singles)
         zw(k) = sum(b_ov_o(:, :, k)*i_other)
      end do
      p(7) = -8.0_dp*dot_product(ww, zw) + 8.0_dp*dot_product(xw, zw) &
             + 16.0_dp*dot_product(yw, zw)
      deallocate (ww, xw, yw, zw, i_ov, i_other)
   end subroutine ind220_one

   subroutine run_sapt2(mols, terms, error)
      !! Every SAPT0 and SAPT2 term, and both totals
      !!
      !!     E(SAPT2) = E(SAPT0) + Elst12 + Exch11 + Exch12
      !!                + Ind22 + Exch-Ind22
      !!
      !! from `sapt2.cc:262`. The SAPT0 part goes through exactly the code
      !! `run_sapt0` runs, so zeroed amplitudes reproduce it identically.
      !! `Exch-Ind22` is `Ind22` scaled by `Exch-Ind20,r / Ind20,r`
      !! (`ind22.cc:52`); psi4 does not guard that division, and for a
      !! weakly-polarizable pair at long range the ratio is noise over noise, so
      !! below `IND20_FLOOR` the scaled term is set to zero instead.
      type(sapt_molecules_t), intent(inout) :: mols
      type(sapt_terms_t), intent(out) :: terms
      type(error_t), intent(inout) :: error

      real(dp), parameter :: IND20_FLOOR = 1.0e-12_dp
      type(sapt_cache_t) :: c
      type(sapt2_cache_t) :: c2
      type(sapt2_amps_t) :: amps_a, amps_b
      real(dp), allocatable :: chf_a(:, :), chf_b(:, :)
      real(dp), allocatable :: u_ar(:, :), v_ar(:, :), u_bs(:, :), v_bs(:, :)

      call build_sapt_cache(mols, c, error)
      if (error%has_error()) return
      call sapt0_from_cache(mols, c, terms, error, chf_a, chf_b)
      if (error%has_error()) return
      call build_sapt2_cache(c, c2, error)
      call c%destroy()
      if (error%has_error()) return

      call build_sapt2_amps(c2, "A", amps_a)
      call build_sapt2_amps(c2, "B", amps_b)
      call sapt2_k2f(c2, u_ar, v_ar, u_bs, v_bs)

      terms%elst12 = sapt_elst12(c2, amps_a, amps_b, chf_a, chf_b)
      terms%exch11 = sapt_exch11(c2, amps_a, amps_b)
      terms%exch12 = sapt_exch12(c2, amps_a, amps_b, u_ar, u_bs)
      call sapt_ind22(c2, amps_a, amps_b, terms%ind22)

      if (abs(terms%ind20_r) > IND20_FLOOR) then
         terms%exch_ind22 = terms%ind22*(terms%exch_ind20_r/terms%ind20_r)
      else
         terms%exch_ind22 = 0.0_dp
      end if

      terms%total_sapt2 = terms%total + terms%elst12 + terms%exch11 &
                          + terms%exch12 + terms%ind22 + terms%exch_ind22

      call amps_a%destroy()
      call amps_b%destroy()
      call c2%destroy()
   end subroutine run_sapt2

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

   subroutine sapt2_cache_destroy(self)
      class(sapt2_cache_t), intent(inout) :: self

      if (allocated(self%eps_a)) deallocate (self%eps_a)
      if (allocated(self%eps_b)) deallocate (self%eps_b)
      if (allocated(self%s_ab)) deallocate (self%s_ab)
      if (allocated(self%v_baa)) deallocate (self%v_baa)
      if (allocated(self%v_abb)) deallocate (self%v_abb)
      if (allocated(self%v_aab)) deallocate (self%v_aab)
      if (allocated(self%v_bab)) deallocate (self%v_bab)
      if (allocated(self%w_baa)) deallocate (self%w_baa)
      if (allocated(self%w_bar)) deallocate (self%w_bar)
      if (allocated(self%w_brr)) deallocate (self%w_brr)
      if (allocated(self%w_abb)) deallocate (self%w_abb)
      if (allocated(self%w_abs)) deallocate (self%w_abs)
      if (allocated(self%w_ass)) deallocate (self%w_ass)
      if (allocated(self%b_aa)) deallocate (self%b_aa)
      if (allocated(self%b_ar)) deallocate (self%b_ar)
      if (allocated(self%b_rr)) deallocate (self%b_rr)
      if (allocated(self%b_bb)) deallocate (self%b_bb)
      if (allocated(self%b_bs)) deallocate (self%b_bs)
      if (allocated(self%b_ss)) deallocate (self%b_ss)
      if (allocated(self%b_ab)) deallocate (self%b_ab)
      if (allocated(self%b_as)) deallocate (self%b_as)
      if (allocated(self%b_rb)) deallocate (self%b_rb)
      if (allocated(self%diag_aa)) deallocate (self%diag_aa)
      if (allocated(self%diag_bb)) deallocate (self%diag_bb)
      self%ndf = 0
      self%nd3 = 0
   end subroutine sapt2_cache_destroy

   subroutine sapt2_amps_destroy(self)
      class(sapt2_amps_t), intent(inout) :: self

      if (allocated(self%t)) deallocate (self%t)
      if (allocated(self%theta)) deallocate (self%theta)
      if (allocated(self%th_ov)) deallocate (self%th_ov)
      if (allocated(self%p_oo)) deallocate (self%p_oo)
      if (allocated(self%p_vv)) deallocate (self%p_vv)
      if (allocated(self%y2)) deallocate (self%y2)
      if (allocated(self%t_singles)) deallocate (self%t_singles)
      if (allocated(self%t2)) deallocate (self%t2)
      if (allocated(self%th2_ov)) deallocate (self%th2_ov)
   end subroutine sapt2_amps_destroy

   subroutine sapt_molecules_destroy(self)
      class(sapt_molecules_t), intent(inout) :: self

      call self%dimer%destroy()
      call self%mono_a%destroy()
      call self%mono_b%destroy()
      self%n_atoms_a = 0
      self%n_atoms_b = 0
   end subroutine sapt_molecules_destroy

end module mqc_sapt
