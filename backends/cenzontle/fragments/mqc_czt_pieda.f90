module mqc_czt_pieda
   !! Pair interaction energy decomposition (PIEDA), GAMESS's `IPIEDA=1`
   !!
   !! The numerics only PIEDA needs: occupied orbitals recovered from a
   !! converged closed-shell density, the union ("higher-level", HL) state of
   !! two monomers, and its internal energy. Which pairs are decomposed, and
   !! `Ees` itself (`es_dimer_energy`), live in [[mqc_czt_fmo]].
   !!
   !!     Ees      = es_dimer_energy(I, J)
   !!     Eex      = E'^HL - E'_I - E'_J - Ees
   !!     Ect+mix  = dE_IJ - Ees - Eex - Edi  (a residual, not charge transfer alone)
   !!
   !! with `E'` an internal energy, `E - Tr(D u)`, and
   !! `D_HL = 2 C (C^T S C)^-1 C^T` for `C = [C_I, C_J]`. `E'^HL` and the two
   !! `E'` are those of the reference: Hartree-Fock's or the Kohn-Sham
   !! functional's, without the correlation of an MP2-family method.
   !!
   !! `Edi` is the pair's correlation interaction `Ec(IJ) - Ec(I) - Ec(J)` for
   !! an MP2-family method, which is inside `dE_IJ`. It is zero here otherwise:
   !! the empirical `Edi` (`fmo_options_t%pieda_dispersion`) is in no
   !! `dE_IJ` and is computed in [[mqc_czt_fmo]] by `pieda_dispersion_term`.
   use pic_types, only: dp, default_int
   use pic_io, only: to_char
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_direct, only: build_fock_direct, direct_stats_t
   use mqc_czt_xc, only: xc_context_t
   use mqc_czt_rhf, only: density_energy
   implicit none
   private

   public :: pieda_pair_terms_t
   public :: cholesky_occupied_orbitals
   public :: hl_density_from_orbitals
   public :: hl_prime_energy
   public :: combine_pieda_terms

   type :: pieda_pair_terms_t
      !! One pair's decomposition, valid only where `decomposed` is true
      logical :: decomposed = .false.
      real(dp) :: ees = 0.0_dp       !! Electrostatics, `es_dimer_energy`
      real(dp) :: eex = 0.0_dp       !! Exact HL exchange
      real(dp) :: ect_mix = 0.0_dp   !! Residual: charge transfer, mixing, response
      real(dp) :: edi = 0.0_dp
         !! `Ec(IJ) - Ec(I) - Ec(J)` for an MP2-family method, inside `dE_IJ`;
         !! zero otherwise
   end type pieda_pair_terms_t

   interface
      subroutine dpstrf(uplo, n, a, lda, piv, rank, tol, work, info)
         !! LAPACK's pivoted Cholesky with rank detection, not wrapped by
         !! pic-blas: no caller elsewhere in this project has needed the rank
         !! it recovers, only the factorisation `pic_potrf` already gives an
         !! unpivoted form of.
         import :: dp, default_int
         implicit none
         character, intent(in) :: uplo
         integer(default_int), intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda, n)
         integer(default_int), intent(out) :: piv(n)
         integer(default_int), intent(out) :: rank
         real(dp), intent(in) :: tol
         real(dp), intent(out) :: work(2*n)
         integer(default_int), intent(out) :: info
      end subroutine dpstrf
   end interface

contains

   subroutine cholesky_occupied_orbitals(density, n_occ, label, c_occ, error)
      !! `C_occ` from a converged closed-shell density, with its rank checked
      !!
      !! `density = 2 C_occ C_occ^T` for an idempotent closed-shell density, so
      !! a pivoted Cholesky of `density/2` returns `C_occ` up to an orthogonal
      !! rotation within the occupied space. The factorisation's own rank is a
      !! free electron-count check: a mismatch against `n_occ` means `density`
      !! is not the converged projector this depends on, and is refused as an
      !! error rather than a warning.
      real(dp), intent(in) :: density(:, :)
      integer, intent(in) :: n_occ
      character(len=*), intent(in) :: label
         !! Named in any error, so a rank mismatch points at which monomer
      real(dp), allocatable, intent(out) :: c_occ(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: a(:, :), work(:)
      integer(default_int), allocatable :: piv(:)
      integer(default_int) :: n, rank_found, info
      integer :: i, n_ao

      n_ao = size(density, 1)
      n = int(n_ao, default_int)
      allocate (a(n_ao, n_ao), source=density/2.0_dp)
      allocate (piv(n_ao), source=0_default_int)
      allocate (work(2*n_ao), source=0.0_dp)

      ! A negative tolerance takes LAPACK's own default,
      ! N*U*MAX(A(K,K)), which is what "no tolerance known a priori" means.
      call dpstrf("L", n, a, n, piv, rank_found, -1.0_dp, work, info)
      if (info < 0) then
         call error%set(ERROR_VALIDATION, "pieda: dpstrf on "//trim(label)// &
                        " was called with an illegal argument (position "// &
                        to_char(int(-info))//"); this is an internal error")
         return
      end if
      if (int(rank_found) /= n_occ) then
         call error%set(ERROR_VALIDATION, "pieda: "//trim(label)//"'s density has "// &
                        "Cholesky rank "//to_char(int(rank_found))//" where "// &
                        to_char(n_occ)//" occupied orbitals were expected -- it is "// &
                        "not a converged closed-shell projector")
         return
      end if

      ! LAPACK writes the lower factor into `a`'s lower triangle and is silent
      ! about the rest: the strict upper triangle still holds `density/2`'s
      ! own elements, not zero, unless cleared by hand.
      do i = 1, n_ao - 1
         a(i, i + 1:n_ao) = 0.0_dp
      end do

      allocate (c_occ(n_ao, n_occ), source=0.0_dp)
      do i = 1, n_ao
         c_occ(int(piv(i)), 1:n_occ) = a(i, 1:n_occ)
      end do
   end subroutine cholesky_occupied_orbitals

   subroutine hl_density_from_orbitals(c, s, d_hl, error)
      !! `D_HL = 2 C (C^T S C)^-1 C^T`, the union's own projector density
      !!
      !! Depends only on the span of `c`'s columns, so `c`'s two blocks -- one
      !! monomer's occupied orbitals each -- need no Gram-Schmidt against each
      !! other first.
      real(dp), intent(in) :: c(:, :)
      real(dp), intent(in) :: s(:, :)
      real(dp), allocatable, intent(out) :: d_hl(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: m(:, :), x(:, :)
      integer(default_int), allocatable :: ipiv(:)
      integer(default_int) :: info
      integer :: k, n_ao

      k = size(c, 2)
      n_ao = size(c, 1)
      allocate (m(k, k))
      m = matmul(transpose(c), matmul(s, c))
      allocate (x(k, n_ao), source=transpose(c))
      allocate (ipiv(k))

      call pic_getrf(m, ipiv, info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "pieda: the union of the two monomers' "// &
                        "occupied spaces is linearly dependent (dgetrf info = "// &
                        to_char(int(info))//")")
         return
      end if
      call pic_getrs(m, ipiv, x, info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "pieda: solving for the HL density failed "// &
                        "(dgetrs info = "//to_char(int(info))//")")
         return
      end if

      allocate (d_hl(n_ao, n_ao))
      d_hl = 2.0_dp*matmul(c, x)
   end subroutine hl_density_from_orbitals

   subroutine hl_prime_energy(mol, bounds, d_hl, e_hl_prime, error, xc)
      !! `E'^HL`, the union state's internal energy, from one Fock build
      !!
      !! Hartree, nuclear repulsion of `mol` included. `bounds` are `mol`'s
      !! Schwarz bounds. No frozen-orbital projector enters: this is GAMESS's
      !! iteration-1 PIEDA energy after its `EPROJ` correction, which removes
      !! exactly the projector's contribution.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: bounds(:, :)
      real(dp), intent(in) :: d_hl(:, :)
      real(dp), intent(out) :: e_hl_prime
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional :: xc
         !! A Kohn-Sham functional's context, built on `mol`: the energy is
         !! then that functional's (`E_xc[D_HL]`, its exact-exchange fraction
         !! and range separation), one exchange-correlation quadrature more.
         !! Absent is Hartree-Fock.

      real(dp), allocatable :: h(:, :), fock(:, :)
      type(direct_stats_t) :: stats

      ! `E[D; h + u] - Tr(D u)` is `E[D; h]` for any one-electron `u`, whether
      ! `E` is Hartree-Fock's or a Kohn-Sham functional's (`E_xc` does not see
      ! `u`), so the field the pair was solved in never needs building here.
      call mol%core_hamiltonian(h)
      if (present(xc)) then
         call density_energy(mol, h, d_hl, bounds, xc, e_hl_prime, error)
         return
      end if
      allocate (fock(size(h, 1), size(h, 2)))
      call build_fock_direct(mol, h, d_hl, bounds, fock, stats, error)
      if (error%has_error()) return

      e_hl_prime = 0.5_dp*sum(d_hl*(h + fock)) + mol%nuclear_repulsion()
   end subroutine hl_prime_energy

   pure subroutine combine_pieda_terms(delta_e, ees, eex, terms, edi)
      !! `Ect+mix = dE_IJ - Ees - Eex - Edi`, the residual that closes the sum
      !!
      !! Taken after `dE_IJ` is final -- reduced across ranks and reduced by
      !! `subtract_subsets` -- which is why this is a separate step from
      !! `hl_prime_energy` rather than folded into it.
      real(dp), intent(in) :: delta_e, ees, eex
      type(pieda_pair_terms_t), intent(inout) :: terms
      real(dp), intent(in), optional :: edi
         !! The correlation interaction already inside `delta_e`. Absent is
         !! zero.

      terms%ees = ees
      terms%eex = eex
      terms%ect_mix = delta_e - ees - eex
      terms%edi = 0.0_dp
      if (present(edi)) then
         terms%ect_mix = terms%ect_mix - edi
         terms%edi = edi
      end if
      terms%decomposed = .true.
   end subroutine combine_pieda_terms

end module mqc_czt_pieda
