!! Indirect nuclear spin-spin couplings: the Fermi-contact term
module mqc_czt_spin_coupling
   !! What splits an NMR line into a multiplet. Two nuclear moments are
   !! coupled through the electrons: one polarises the spin density around
   !! it, and the other feels that polarisation. Of Ramsey's four terms this
   !! module computes the Fermi-contact one, which is the interaction of a
   !! nuclear moment with the electron spin density *at* the nucleus and is
   !! most of a coupling between hydrogens or between a hydrogen and a
   !! carbon. The spin-dipole and the two spin-orbit terms are not built.
   !!
   !! The contact operator of nucleus `A` acts on the two spins with opposite
   !! signs, so the response to it is a triplet one: the coupled-perturbed
   !! equations are solved against the triplet operator, one right-hand side
   !! per nucleus, and
   !!
   !!     K_AB = 4 alpha^4 sum_ai h^A_ai U^B_ai,
   !!     h^A_ai = (4 pi / 3) phi_a(R_A) phi_i(R_A)
   !!
   !! term for term `pyscf.prop.ssc`'s `make_fc`. Hartree-Fock or Kohn-Sham,
   !! closed shell.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t
   use mqc_physical_constants, only: PI
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_ao, only: eval_ao_block
   use mqc_czt_cphf, only: cphf_solve
   use mqc_czt_xc, only: xc_context_t
   use mqc_nmr_isotopes, only: nmr_isotope, reduced_to_hertz
   implicit none
   private

   public :: czt_fermi_contact
   public :: print_coupling_report

   real(dp), parameter :: FINE_STRUCTURE = 7.2973525693e-3_dp
   real(dp), parameter :: RESPONSE_TOL = 1.0e-10_dp
   integer, parameter :: RESPONSE_MAX_ITER = 200
   integer, parameter :: LINE_LEN = 160
   integer, parameter :: ISOTOPE_LEN = 8

contains

   subroutine czt_fermi_contact(mol, orbitals, orbital_energies, n_occ, coupling, error, &
                                xc, density)
      !! The Fermi-contact reduced coupling between every pair of nuclei
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)          !! Canonical, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: coupling(:, :)
         !! (n_atoms, n_atoms), symmetric, zero on the diagonal. The reduced
         !! coupling in atomic units: the second derivative of the energy
         !! with respect to the two nuclear moments, each in units of the
         !! nuclear magneton times its g factor. `reduced_to_hertz` turns an
         !! element into `J` for a pair of isotopes.
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional :: xc
         !! The reference's functional. Absent or inactive is Hartree-Fock.
      real(dp), intent(in), optional :: density(:, :)
         !! The reference density. Required with an active `xc`.

      real(dp), allocatable :: ao(:, :), mo(:, :), rhs(:, :, :), u(:, :, :)
      integer :: natm, n_mo, n_vir, a, i, ia, ib
      logical :: kohn_sham

      if (error%has_error()) return
      natm = mol%natm
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ
      kohn_sham = .false.
      if (present(xc)) kohn_sham = xc%active

      ! Every orbital at every nucleus.
      call eval_ao_block(mol, mol%coords, ao, error)
      if (error%has_error()) return
      allocate (mo(natm, n_mo))
      call pic_gemm(ao, orbitals, mo)

      allocate (rhs(n_vir, n_occ, natm))
      do ia = 1, natm
         do i = 1, n_occ
            do a = 1, n_vir
               rhs(a, i, ia) = 4.0_dp*PI/3.0_dp*mo(ia, n_occ + a)*mo(ia, i)
            end do
         end do
      end do

      if (kohn_sham) then
         call cphf_solve(mol, orbitals, orbital_energies, n_occ, response=u, error=error, &
                         max_iter=RESPONSE_MAX_ITER, tol=RESPONSE_TOL, mo_rhs=rhs, &
                         xc=xc, density=density, triplet=.true.)
      else
         call cphf_solve(mol, orbitals, orbital_energies, n_occ, response=u, error=error, &
                         max_iter=RESPONSE_MAX_ITER, tol=RESPONSE_TOL, mo_rhs=rhs, &
                         triplet=.true.)
      end if
      if (error%has_error()) then
         call error%add_context("spin-spin coupling: the triplet response to a "// &
                                "nucleus's contact field")
         return
      end if

      allocate (coupling(natm, natm))
      coupling = 0.0_dp
      do ib = 1, natm
         do ia = 1, natm
            if (ia == ib) cycle
            coupling(ia, ib) = 4.0_dp*FINE_STRUCTURE**4*sum(rhs(:, :, ia)*u(:, :, ib))
         end do
      end do
      ! Symmetric by the operator's symmetry; the two triangles differ by
      ! the solves' residual and are averaged.
      coupling = 0.5_dp*(coupling + transpose(coupling))
   end subroutine czt_fermi_contact

   subroutine print_coupling_report(coupling, element_numbers, symbols)
      !! The couplings in Hz, one line per pair, for the observed isotopes
      real(dp), intent(in) :: coupling(:, :)        !! Reduced, atomic units
      integer, intent(in) :: element_numbers(:)
      character(len=*), intent(in) :: symbols(:)

      character(len=LINE_LEN) :: line
      character(len=ISOTOPE_LEN) :: iso_a, iso_b
      real(dp) :: g_a, g_b
      integer :: ia, ib

      call logger%info("")
      call logger%info("  Spin-spin couplings, Fermi-contact term (Hz)")
      call logger%info("  ----------------------------------------------------")
      do ia = 1, size(coupling, 1)
         call nmr_isotope(element_numbers(ia), iso_a, g_a)
         do ib = ia + 1, size(coupling, 1)
            call nmr_isotope(element_numbers(ib), iso_b, g_b)
            if (len_trim(iso_a) == 0 .or. len_trim(iso_b) == 0) cycle
            write (line, "(4x,a,1x,i0,t16,a,1x,i0,t30,a,a,a,t46,f12.4)") &
               trim(symbols(ia)), ia, trim(symbols(ib)), ib, trim(iso_a), "-", trim(iso_b), &
               reduced_to_hertz(coupling(ia, ib), g_a, g_b)
            call logger%info(trim(line))
         end do
      end do
      call logger%info("")
   end subroutine print_coupling_report

end module mqc_czt_spin_coupling
