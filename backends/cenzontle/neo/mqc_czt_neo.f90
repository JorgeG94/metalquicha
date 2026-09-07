module mqc_czt_neo
   !! Nuclear-electronic orbital Hartree-Fock: chosen nuclei get orbitals too
   !!
   !! The nuclear-electronic orbital (NEO) method of Hammes-Schiffer treats
   !! selected nuclei -- protons, here -- quantum mechanically, in a Gaussian
   !! basis of their own centred where the classical nucleus stood, and solves
   !! them self-consistently with the electrons. What comes out is a wavefunction
   !! that carries the proton's zero-point motion and delocalisation, rather than
   !! a harmonic correction bolted on afterwards.
   !!
   !! **Each quantum nucleus is one particle in its own component**, so inside a
   !! component there is no exchange and no self-Coulomb; a component's Fock
   !! matrix is its one-body Hamiltonian plus the Coulomb field of every *other*
   !! component. The nuclear charge of a quantised nucleus leaves the electronic
   !! Hamiltonian entirely: the electrons see it only through its density.
   !!
   !!     E = E_e[D_e; classical nuclei]  + sum_p <D_p| T/m_p + V_classical |D_p>
   !!       - sum_p J(D_e, D_p)           + sum_{p<q} J(D_p, D_q)
   !!
   !!     F_e = h_e + G_e[D_e] - sum_p J_ee[D_p]                  the last term is `h_extra`
   !!     F_p = T/m_p + V_classical - J_pp[D_e] + sum_{q/=p} J_pp[D_q]
   !!
   !! **The cross Coulomb terms need no new integrals.** `J_ee[D_p]` is
   !! `sum_pq (mu nu | p q) D_pq` over two different bases, and a *combined*
   !! molecule carrying the electronic basis on every atom and the proton basis
   !! on a second copy of each quantum centre makes those ordinary quartets: a
   !! Coulomb-only batched build with a block density `0 (+) D_p` returns
   !! `J_ee[D_p]` in its electronic block, and with `D_e (+) 0` returns
   !! `J_pp[D_e]` in the proton block. Density screening skips the quartets a
   !! block density never touches, and the proton basis is a few dozen
   !! functions, so the cross terms cost a fraction of one Fock build.
   !!
   !! **Coupling is by macro-iteration.** The electrons converge in the field of
   !! the current proton densities, then every proton is re-solved in the field
   !! of the new electrons; the loop stops when the total energy and the proton
   !! densities stop moving. NEO-HF converges in a handful of cycles this way.
   !! PySCF-NEO runs one DIIS over all components at once, which is the upgrade
   !! if a case ever needs it.
   !!
   !! **NEO-DFT** runs the electrons as Kohn-Sham through the same SCF, and adds
   !! the electron-proton correlation functional of Yang, Brorsen, Culpitt, Pak
   !! and Hammes-Schiffer (epc17), a local functional of the two densities on
   !! the electronic grid, restricted to the points a proton basis reaches:
   !!
   !!     E_epc = -int rho_e rho_p / (a - b sqrt(rho_e rho_p) + c rho_e rho_p)
   !!
   !! Its electronic potential rides in `h_extra` with the proton Coulomb term,
   !! frozen for one electronic SCF and refreshed every macro-iteration -- the
   !! fixed point is the same as a full Kohn-Sham treatment, only reached from
   !! outside -- and the energy replaces the frozen `Tr(D_e V_epc)` the SCF
   !! counted by the functional itself.
   !!
   !! The reference implementation is the `pyscf/neo` module of Yang Yang's
   !! PySCF fork (github.com/theorychemyang/pyscf); its mass convention, its
   !! proton basis sets (`basis_sets/neo/`) and its HCN test energies are what
   !! this is validated against.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type
   use mqc_basis_utils, only: find_basis_file
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_direct, only: build_fock_direct_many, schwarz_bounds, direct_stats_t
   use mqc_czt_ao, only: eval_ao_block, eval_rho
   use mqc_czt_xc, only: xc_context_t, xc_context_create
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid
   use mqc_program_limits, only: MAX_LINE_LENGTH
   implicit none
   private

   public :: neo_result_t
   public :: run_czt_neo_hf
   public :: PROTON_MASS

   real(dp), parameter :: PROTON_MASS = 1836.1526473649712_dp
      !! In electron masses. PySCF-NEO's convention: the most common isotope's
      !! atomic mass (1.00782503223 u) minus one electron, over the electron
      !! mass. CODATA's m_p/m_e is 1836.15267343; the difference moves a
      !! proton's kinetic energy in the eighth decimal.

   integer, parameter :: MAX_MACRO = 200
      !! Macro-iterations before giving up; NEO-HF takes ten or twenty, NEO-DFT
      !! with epc several times that
   integer, parameter :: MAX_INNER = 500
      !! Proton self-consistency steps per macro-iteration, with a functional
   real(dp), parameter :: INNER_MIX = 0.2_dp
      !! The proton's own map under epc17 is far from contractive -- half and
      !! half oscillates -- so it is mixed this gently and takes a hundred cheap
      !! steps
   real(dp), parameter :: MIX_START = 1.0_dp
      !! How much of a proton's new density enters between cycles, with a
      !! functional; halved whenever the residual grows
   real(dp), parameter :: MIX_FLOOR = 0.05_dp

   type :: neo_result_t
      real(dp) :: energy = 0.0_dp
         !! Total: electrons, classical nuclei, quantum nuclei and every coupling
      real(dp) :: electronic = 0.0_dp
         !! `E_e` with the classical nuclear repulsion and the electron-proton
         !! attraction, i.e. the RHF energy in the field of the proton densities
      real(dp) :: nuclear_one_body = 0.0_dp
         !! `sum_p <D_p| T/m_p + V_classical |D_p>`
      real(dp) :: electron_nucleus = 0.0_dp
         !! `-sum_p J(D_e, D_p)`, already inside `electronic`; reported apart
      real(dp) :: nucleus_nucleus = 0.0_dp
         !! `sum_{p<q} J(D_p, D_q)`, zero with one quantum nucleus
      real(dp) :: epc_energy = 0.0_dp
         !! The electron-proton correlation energy, zero without a functional
      logical :: kohn_sham = .false.       !! Whether the electrons ran as DFT
      integer :: n_quantum = 0
      integer :: iterations = 0                !! Macro-iterations run
      logical :: converged = .false.
      logical :: cartesian = .false.           !! The angular form every basis was built in
      type(rhf_result_t) :: electrons          !! The last electronic SCF
      real(dp), allocatable :: nuclear_orbital_energies(:, :)  !! (nao_p, n_quantum)
      real(dp), allocatable :: nuclear_orbitals(:, :, :)       !! (nao_p, nao_p, n_quantum)
      real(dp), allocatable :: nuclear_densities(:, :, :)      !! (nao_p, nao_p, n_quantum), `c c^T`
      real(dp), allocatable :: nuclear_overlaps(:, :, :)       !! (nao_p, nao_p, n_quantum)
   end type neo_result_t

contains

   subroutine run_czt_neo_hf(atomic_numbers, element_symbols, coordinates, basis_name, &
                             nuclear_basis, quantum, nelec, max_iter, energy_tol, &
                             density_tol, verbose, result, error, force_cartesian, in_core, &
                             functional, grid_level, epc)
      !! NEO-HF, or NEO-DFT, for a closed-shell electronic structure and any number of quantum protons
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)     !! (3, natm), Bohr
      character(len=*), intent(in) :: basis_name    !! The electronic basis
      character(len=*), intent(in) :: nuclear_basis  !! The proton basis, e.g. `pb4-d`
      logical, intent(in) :: quantum(:)             !! Which atoms are quantised
      integer, intent(in) :: nelec
      integer, intent(in) :: max_iter               !! Per electronic SCF
      real(dp), intent(in) :: energy_tol, density_tol
         !! The macro-iteration stops when the total energy moves by less than
         !! the first and every proton density by less than the second; the
         !! electronic SCFs inside use the same two.
      logical, intent(in) :: verbose
      type(neo_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: force_cartesian
      logical, intent(in), optional :: in_core
      character(len=*), intent(in), optional :: functional
         !! An exchange-correlation functional for the electrons; absent or
         !! empty is Hartree-Fock
      integer, intent(in), optional :: grid_level
         !! The DFT grid level, for the functional and for the epc integration
      character(len=*), intent(in), optional :: epc
         !! Electron-proton correlation: "17-1", "17-2", or empty for none.
         !! Needs a functional: a Hartree-Fock electron with a correlation
         !! functional bolted on is refused.

      type(czt_molecule_t) :: mol_e, mol_c
      type(czt_molecule_t), allocatable :: mol_p(:)
      type(molecular_basis_type) :: e_basis, n_basis, one_basis, c_basis
      character(len=:), allocatable :: e_path, n_path
      integer, allocatable :: which(:), z_c(:), offset(:)
      logical, allocatable :: ghost_c(:)
      real(dp), allocatable :: xyz_c(:, :)
      real(dp), allocatable :: s_p(:, :, :), h_p(:, :, :), x_p(:, :, :), d_p(:, :, :)
      real(dp), allocatable :: c_p(:, :, :), eps_p(:, :)
      real(dp), allocatable :: work(:, :), h_extra(:, :), zero_h(:, :), bounds(:, :)
      real(dp), allocatable :: dens(:, :, :), j_p(:, :, :), j_e(:, :, :), f_p(:, :)
      real(dp), allocatable :: d_e_prev(:, :), t_p(:, :), hc_p(:, :)
      type(direct_stats_t) :: stats
      type(error_t) :: read_error
      type(xc_context_t), target :: xc
      type(xc_context_t), pointer :: xc_arg
      type(dft_grid_t) :: grid
      logical :: kohn_sham, use_epc
      real(dp) :: epc_a, epc_b, epc_c, alpha_min, r_cut, e_epc
      integer, allocatable :: sel(:)
      real(dp), allocatable :: pts(:, :), w_sel(:), ao_e(:, :), ao_p(:, :, :)
      real(dp), allocatable :: rho_e(:), rho_p(:), v_grid(:), v_epc_e(:, :), v_epc_p(:, :, :)
      real(dp), allocatable :: ao_one(:, :), f_fixed(:, :), d_mix(:, :)
      integer :: inner
      real(dp) :: e_dummy, mix, dd_prev
      integer :: natm, nq, nao_e, nao_p, nao_c, p, q, it, i, lo, hi, qlo, qhi, n_sel, level, g
      real(dp) :: e_total, e_prev, dd, dd_max, e_p1, e_pp, e_ep
      character(len=MAX_LINE_LENGTH) :: line

      natm = size(atomic_numbers)
      if (size(quantum) /= natm) then
         call error%set(ERROR_VALIDATION, "NEO: the quantum-nucleus mask does not "// &
                        "match the atom count")
         return
      end if
      nq = count(quantum)
      if (nq == 0) then
         call error%set(ERROR_VALIDATION, "NEO: no nucleus is marked quantum")
         return
      end if
      allocate (which(nq))
      q = 0
      do i = 1, natm
         if (quantum(i)) then
            q = q + 1
            which(q) = i
            if (atomic_numbers(i) /= 1) then
               call error%set(ERROR_VALIDATION, "NEO: only hydrogen can be quantised "// &
                              "for now; atom "//trim(itoa(i))//" is "// &
                              trim(element_symbols(i)))
               return
            end if
         end if
      end do

      kohn_sham = .false.
      if (present(functional)) kohn_sham = len_trim(functional) > 0
      use_epc = .false.
      if (present(epc)) use_epc = len_trim(epc) > 0
      if (use_epc .and. .not. kohn_sham) then
         call error%set(ERROR_VALIDATION, "NEO: the electron-proton correlation functional "// &
                        "needs a Kohn-Sham electron; name a functional or drop epc")
         return
      end if
      if (use_epc) then
         select case (trim(adjustl(epc)))
         case ("17-1")
            epc_a = 2.35_dp
            epc_b = 2.4_dp
            epc_c = 3.2_dp
         case ("17-2")
            epc_a = 2.35_dp
            epc_b = 2.4_dp
            epc_c = 6.6_dp
         case default
            call error%set(ERROR_VALIDATION, "NEO: unknown electron-proton correlation "// &
                           "functional '"//trim(epc)//"'. Accepted: 17-1, 17-2")
            return
         end select
      end if
      level = 3
      if (present(grid_level)) level = grid_level
      result%kohn_sham = kohn_sham

      ! --- the molecules ------------------------------------------------------
      ! Electrons: every atom with its basis, the quantum ones ghosted, which
      ! takes their charge out of the nuclear attraction and the nuclear
      ! repulsion and leaves everything else as it was.
      call build_czt_molecule(atomic_numbers, element_symbols, coordinates, basis_name, &
                              mol_e, error, force_cartesian=force_cartesian, ghost=quantum)
      if (error%has_error()) return
      result%cartesian = mol_e%cartesian

      call find_basis_file(basis_name, e_path, read_error)
      if (read_error%has_error()) then
         call error%set(ERROR_VALIDATION, "NEO: "//read_error%get_message())
         return
      end if
      call build_molecular_basis_json(e_path, element_symbols, e_basis, read_error)
      if (read_error%has_error()) then
         call error%set(ERROR_VALIDATION, "NEO: could not read "//trim(basis_name)//": "// &
                        read_error%get_message())
         return
      end if
      call find_basis_file(nuclear_basis, n_path, read_error)
      if (read_error%has_error()) then
         call error%set(ERROR_VALIDATION, "NEO: no nuclear basis '"//trim(nuclear_basis)// &
                        "': "//read_error%get_message())
         return
      end if
      call build_molecular_basis_json(n_path, ["H "], n_basis, read_error)
      if (read_error%has_error()) then
         call error%set(ERROR_VALIDATION, "NEO: could not read "//trim(nuclear_basis)// &
                        ": "//read_error%get_message())
         return
      end if
      ! The proton basis says nothing about its angular form; it takes the
      ! electronic one so that every molecule below is built the same way.
      n_basis%angular_form = e_basis%angular_form
      ! Its most diffuse exponent bounds how far a proton density reaches,
      ! which is how far the epc integration has to look from the centre.
      alpha_min = huge(1.0_dp)
      do i = 1, n_basis%elements(1)%nshells
         alpha_min = min(alpha_min, minval(n_basis%elements(1)%shells(i)%exponents))
      end do

      ! One molecule per proton for its one-body terms: every atom present,
      ! for the classical potential, and the proton shells on the quantum
      ! centre alone. All quantum nuclei ghosted, this one included -- a
      ! proton does not feel its own charge, and the others reach it through
      ! their densities.
      allocate (mol_p(nq))
      do p = 1, nq
         call basis_on_one_atom(element_symbols, which(p), n_basis, one_basis)
         call mol_p(p)%build(atomic_numbers, coordinates, one_basis, error, &
                             force_cartesian=force_cartesian, ghost=quantum)
         call one_basis%destroy()
         if (error%has_error()) return
      end do

      ! The combined molecule: the electronic atoms first, in their order, then
      ! a second copy of each quantum centre carrying the proton shells.
      allocate (z_c(natm + nq), xyz_c(3, natm + nq), ghost_c(natm + nq))
      z_c(1:natm) = atomic_numbers
      xyz_c(:, 1:natm) = coordinates
      ghost_c(1:natm) = quantum
      call c_basis%allocate_elements(natm + nq)
      do i = 1, natm
         c_basis%elements(i) = e_basis%elements(i)
      end do
      do p = 1, nq
         z_c(natm + p) = atomic_numbers(which(p))
         xyz_c(:, natm + p) = coordinates(:, which(p))
         ghost_c(natm + p) = .true.
         c_basis%elements(natm + p) = n_basis%elements(1)
      end do
      c_basis%angular_form = e_basis%angular_form
      call mol_c%build(z_c, xyz_c, c_basis, error, force_cartesian=force_cartesian, &
                       ghost=ghost_c)
      call c_basis%destroy()
      call e_basis%destroy()
      call n_basis%destroy()
      if (error%has_error()) return

      nao_e = mol_e%nao
      nao_p = mol_p(1)%nao
      nao_c = mol_c%nao
      if (nao_c /= nao_e + nq*nao_p) then
         call error%set(ERROR_VALIDATION, "NEO: the combined basis does not decompose "// &
                        "into the electronic and proton blocks")
         return
      end if
      allocate (offset(nq))
      do p = 1, nq
         offset(p) = nao_e + (p - 1)*nao_p
      end do

      ! --- proton one-body terms and orthogonalisers ---------------------------
      ! `core_hamiltonian` is `T + V` with the electron's sign on `V`; the proton
      ! has the opposite charge and a mass, so `h_p = T/m - V = T (1/m + 1) - H`.
      allocate (s_p(nao_p, nao_p, nq), h_p(nao_p, nao_p, nq), x_p(nao_p, nao_p, nq))
      allocate (d_p(nao_p, nao_p, nq), c_p(nao_p, nao_p, nq), eps_p(nao_p, nq))
      allocate (t_p(nao_p, nao_p), hc_p(nao_p, nao_p), f_p(nao_p, nao_p))
      allocate (f_fixed(nao_p, nao_p), d_mix(nao_p, nao_p))
      do p = 1, nq
         call mol_p(p)%overlap(work)
         s_p(:, :, p) = work
         call mol_p(p)%kinetic(t_p)
         call mol_p(p)%core_hamiltonian(hc_p)
         h_p(:, :, p) = t_p*(1.0_dp/PROTON_MASS + 1.0_dp) - hc_p
         call inverse_sqrt(s_p(:, :, p), x_p(:, :, p), error)
         if (error%has_error()) return
         ! The starting proton: the lowest state of its one-body Hamiltonian.
         call lowest_state(h_p(:, :, p), x_p(:, :, p), c_p(:, :, p), eps_p(:, p), &
                           d_p(:, :, p), error)
         if (error%has_error()) return
      end do

      call schwarz_bounds(mol_c, bounds, error)
      if (error%has_error()) return
      allocate (zero_h(nao_c, nao_c), h_extra(nao_e, nao_e), dens(nao_c, nao_c, nq))
      zero_h = 0.0_dp

      ! --- Kohn-Sham, and the epc grid ----------------------------------------
      xc_arg => null()
      if (kohn_sham) then
         call xc_context_create(mol_e, trim(functional), xc, error, level=level)
         if (error%has_error()) return
         xc_arg => xc
      end if
      allocate (v_epc_e(nao_e, nao_e), v_epc_p(nao_p, nao_p, nq))
      v_epc_e = 0.0_dp
      v_epc_p = 0.0_dp
      e_epc = 0.0_dp
      n_sel = 0
      if (use_epc) then
         ! The electronic grid, kept only where a proton density is not
         ! negligible: exp(-2 alpha r^2) at `r_cut` is below 1e-11 for the most
         ! diffuse proton function, and everything else decays faster.
         call build_dft_grid(coordinates, atomic_numbers, grid, error, level=level)
         if (error%has_error()) return
         r_cut = sqrt(25.0_dp/alpha_min)
         allocate (sel(grid%n_points))
         do g = 1, grid%n_points
            do p = 1, nq
               if (norm2(grid%coords(:, g) - coordinates(:, which(p))) <= r_cut) then
                  n_sel = n_sel + 1
                  sel(n_sel) = g
                  exit
               end if
            end do
         end do
         if (n_sel == 0) then
            call error%set(ERROR_VALIDATION, "NEO: the epc grid has no points near a proton")
            return
         end if
         allocate (pts(3, n_sel), w_sel(n_sel))
         pts = grid%coords(:, sel(1:n_sel))
         w_sel = grid%weights(sel(1:n_sel))
         call eval_ao_block(mol_e, pts, ao_e, error)
         if (error%has_error()) return
         allocate (ao_p(n_sel, nao_p, nq), rho_p(n_sel), v_grid(n_sel))
         do p = 1, nq
            call eval_ao_block(mol_p(p), pts, ao_one, error)
            if (error%has_error()) return
            ao_p(:, :, p) = ao_one
         end do
         call grid%destroy()
      end if

      if (verbose) then
         write (line, "(A,I0,A,A,I0,A,I0,A)") "  NEO-HF: ", nq, " quantum nucle", &
            trim(merge("us", "i ", nq == 1))//"; ", nao_e, " electronic and ", nao_p, &
            " proton functions each"
         call logger%info(trim(line))
         call logger%info("  ----------------------------------------------------------------")
         call logger%info("      macro                 energy          dE         dD_p   SCF")
         call logger%info("  ----------------------------------------------------------------")
      end if

      e_prev = 0.0_dp
      e_total = 0.0_dp
      mix = MIX_START
      dd_prev = huge(1.0_dp)
      result%converged = .false.
      do it = 1, MAX_MACRO
         ! Pass A: the Coulomb field of every proton, on the combined basis.
         dens = 0.0_dp
         do p = 1, nq
            lo = offset(p) + 1
            hi = offset(p) + nao_p
            dens(lo:hi, lo:hi, p) = d_p(:, :, p)
         end do
         call build_fock_direct_many(mol_c, zero_h, dens, bounds, j_p, stats, error, &
                                     k_scale=0.0_dp, density_screen=.true.)
         if (error%has_error()) return
         h_extra = 0.0_dp
         do p = 1, nq
            h_extra = h_extra - j_p(1:nao_e, 1:nao_e, p)
         end do
         ! The epc potential on the electrons, from the last electronic density;
         ! there is none to build it from on the first cycle.
         v_epc_e = 0.0_dp
         if (use_epc .and. it > 1) then
            call eval_rho(ao_e, d_e_prev, rho_e)
            rho_e = max(rho_e, 0.0_dp)
            v_grid = 0.0_dp
            do p = 1, nq
               call eval_rho(ao_p(:, :, p), d_p(:, :, p), rho_p)
               rho_p = max(rho_p, 0.0_dp)
               call epc17_electron_potential(epc_a, epc_b, epc_c, rho_e, rho_p, v_grid)
            end do
            call grid_matrix(ao_e, w_sel*v_grid, v_epc_e)
            h_extra = h_extra + v_epc_e
         end if
         e_p1 = 0.0_dp
         e_pp = 0.0_dp
         do p = 1, nq
            e_p1 = e_p1 + sum(d_p(:, :, p)*h_p(:, :, p))
            lo = offset(p) + 1
            hi = offset(p) + nao_p
            do q = p + 1, nq
               e_pp = e_pp + sum(d_p(:, :, p)*j_p(lo:hi, lo:hi, q))
            end do
         end do

         ! The electrons, in that field. Warm-started after the first cycle:
         ! `d_e_prev` is unallocated on the first and so absent.
         call run_czt_rhf(mol_e, nelec, max_iter, energy_tol, density_tol, .false., &
                          result%electrons, error, h_extra=h_extra, in_core=in_core, &
                          guess_density=d_e_prev, xc=xc_arg)
         if (error%has_error()) return
         if (.not. result%electrons%converged) then
            call error%set(ERROR_VALIDATION, "NEO: the electronic SCF did not converge "// &
                           "in macro-iteration "//trim(itoa(it)))
            return
         end if
         d_e_prev = result%electrons%density
         e_ep = sum(result%electrons%density*(h_extra - v_epc_e))
         ! The functional's energy replaces the frozen potential's `Tr(D V)`
         ! that the SCF counted, and the proton potentials come off the new
         ! electronic density.
         e_epc = 0.0_dp
         if (use_epc) then
            call eval_rho(ao_e, result%electrons%density, rho_e)
            rho_e = max(rho_e, 0.0_dp)
            do p = 1, nq
               call eval_rho(ao_p(:, :, p), d_p(:, :, p), rho_p)
               rho_p = max(rho_p, 0.0_dp)
               call epc17_proton_terms(epc_a, epc_b, epc_c, rho_e, rho_p, w_sel, v_grid, e_epc)
               call grid_matrix(ao_p(:, :, p), w_sel*v_grid, v_epc_p(:, :, p))
            end do
         end if
         e_total = result%electrons%energy - sum(result%electrons%density*v_epc_e) &
                   + e_epc + e_p1 + e_pp

         ! Pass B: the electrons' Coulomb field on the proton blocks.
         deallocate (dens)
         allocate (dens(nao_c, nao_c, 1))
         dens = 0.0_dp
         dens(1:nao_e, 1:nao_e, 1) = result%electrons%density
         call build_fock_direct_many(mol_c, zero_h, dens, bounds, j_e, stats, error, &
                                     k_scale=0.0_dp, density_screen=.true.)
         if (error%has_error()) return
         deallocate (dens)
         allocate (dens(nao_c, nao_c, nq))

         ! Every proton, in the field of the new electrons and the other protons.
         ! With a correlation functional the proton's own potential depends on
         ! its density, and one diagonalisation per cycle overshoots -- the
         ! macro-iteration then flips between two states forever -- so each
         ! proton is taken to self-consistency at fixed electrons first, with
         ! the density mixed half and half between steps.
         dd_max = 0.0_dp
         do p = 1, nq
            lo = offset(p) + 1
            hi = offset(p) + nao_p
            f_fixed = h_p(:, :, p) - j_e(lo:hi, lo:hi, 1)
            do q = 1, nq
               if (q == p) cycle
               qlo = offset(p) + 1
               qhi = offset(p) + nao_p
               f_fixed = f_fixed + j_p(qlo:qhi, qlo:qhi, q)
            end do
            d_mix = d_p(:, :, p)
            do inner = 1, MAX_INNER
               f_p = f_fixed + v_epc_p(:, :, p)
               call lowest_state(f_p, x_p(:, :, p), c_p(:, :, p), eps_p(:, p), work, error)
               if (error%has_error()) return
               if (.not. use_epc) exit
               dd = maxval(abs(work - d_mix))
               if (dd < density_tol) exit
               d_mix = (1.0_dp - INNER_MIX)*d_mix + INNER_MIX*work
               call eval_rho(ao_p(:, :, p), d_mix, rho_p)
               rho_p = max(rho_p, 0.0_dp)
               e_dummy = 0.0_dp
               call epc17_proton_terms(epc_a, epc_b, epc_c, rho_e, rho_p, w_sel, v_grid, e_dummy)
               call grid_matrix(ao_p(:, :, p), w_sel*v_grid, v_epc_p(:, :, p))
            end do
            dd = maxval(abs(work - d_p(:, :, p)))
            dd_max = max(dd_max, dd)
            ! Damped between cycles with a functional when the residual grows:
            ! the electrons answer the proton's move and the proton answers
            ! back, and a step that is too long sends the two around a limit
            ! cycle at a few tenths in the density.
            if (use_epc) then
               d_p(:, :, p) = (1.0_dp - mix)*d_p(:, :, p) + mix*work
            else
               d_p(:, :, p) = work
            end if
         end do
         if (use_epc .and. it > 1) then
            if (dd_max > dd_prev) mix = max(0.5_dp*mix, MIX_FLOOR)
         end if
         dd_prev = dd_max

         if (verbose) then
            write (line, "(I10,F24.12,2ES12.3,I6)") it, e_total, e_total - e_prev, dd_max, &
               result%electrons%iterations
            call logger%info(trim(line))
         end if
         result%iterations = it
         if (it > 1) then
            if (abs(e_total - e_prev) < energy_tol .and. dd_max < density_tol) then
               result%converged = .true.
               exit
            end if
         end if
         e_prev = e_total
      end do
      if (verbose) then
         call logger%info("  ----------------------------------------------------------------")
      end if
      if (.not. result%converged) then
         call error%set(ERROR_VALIDATION, "NEO: the macro-iteration did not converge in "// &
                        trim(itoa(MAX_MACRO))//" cycles")
         return
      end if

      result%energy = e_total
      result%electronic = result%electrons%energy
      result%nuclear_one_body = e_p1
      result%electron_nucleus = e_ep
      result%nucleus_nucleus = e_pp
      result%epc_energy = e_epc
      result%n_quantum = nq
      call move_alloc(eps_p, result%nuclear_orbital_energies)
      call move_alloc(c_p, result%nuclear_orbitals)
      call move_alloc(d_p, result%nuclear_densities)
      call move_alloc(s_p, result%nuclear_overlaps)

      if (verbose) then
         write (line, "(A,F20.10)") "  NEO-HF energy      ", result%energy
         call logger%info(trim(line))
         write (line, "(A,F20.10)") "    electrons + classical nuclei + e-p ", result%electronic
         call logger%info(trim(line))
         write (line, "(A,F20.10)") "    of which electron-proton           ", result%electron_nucleus
         call logger%info(trim(line))
         write (line, "(A,F20.10)") "    proton one-body (T/m + V)          ", result%nuclear_one_body
         call logger%info(trim(line))
         if (nq > 1) then
            write (line, "(A,F20.10)") "    proton-proton                      ", result%nucleus_nucleus
            call logger%info(trim(line))
         end if
         if (use_epc) then
            write (line, "(A,F20.10)") "    electron-proton correlation (epc)  ", result%epc_energy
            call logger%info(trim(line))
         end if
         do p = 1, nq
            write (line, "(A,I0,A,F14.8)") "    proton ", p, " orbital energy ", &
               result%nuclear_orbital_energies(1, p)
            call logger%info(trim(line))
         end do
      end if

      do p = 1, nq
         call mol_p(p)%destroy()
      end do
      call mol_c%destroy()
      call mol_e%destroy()
      if (kohn_sham) call xc%destroy()
   end subroutine run_czt_neo_hf

   subroutine epc17_electron_potential(a, b, c, rho_e, rho_p, v)
      !! `dE_epc/drho_e` at every point, **added** to `v`, one proton at a time
      !!
      !!     E_epc = -int rho_e rho_p / (a - b sqrt(rho_e rho_p) + c rho_e rho_p)
      !!
      !! epc17 of Yang, Brorsen, Culpitt, Pak and Hammes-Schiffer, J. Chem.
      !! Phys. 147, 114113 (2017), in the form PySCF-NEO evaluates it.
      real(dp), intent(in) :: a, b, c
      real(dp), intent(in) :: rho_e(:), rho_p(:)
      real(dp), intent(inout) :: v(:)
      real(dp) :: prod, root, denom
      integer :: g
      do g = 1, size(rho_e)
         prod = rho_e(g)*rho_p(g)
         root = sqrt(prod)
         denom = a - b*root + c*prod
         v(g) = v(g) + (-a*rho_p(g) + 0.5_dp*b*rho_p(g)*root)/denom**2
      end do
   end subroutine epc17_electron_potential

   subroutine epc17_proton_terms(a, b, c, rho_e, rho_p, w, v, energy)
      !! `dE_epc/drho_p` at every point, and the energy **added** to `energy`
      real(dp), intent(in) :: a, b, c
      real(dp), intent(in) :: rho_e(:), rho_p(:), w(:)
      real(dp), intent(out) :: v(:)
      real(dp), intent(inout) :: energy
      real(dp) :: prod, root, denom
      integer :: g
      do g = 1, size(rho_e)
         prod = rho_e(g)*rho_p(g)
         root = sqrt(prod)
         denom = a - b*root + c*prod
         v(g) = (-a*rho_e(g) + 0.5_dp*b*rho_e(g)*root)/denom**2
         energy = energy - w(g)*prod/denom
      end do
   end subroutine epc17_proton_terms

   subroutine grid_matrix(ao, wv, v)
      !! `V_uv = sum_g chi_u(g) wv(g) chi_v(g)`, a potential on the grid as a matrix
      real(dp), intent(in) :: ao(:, :)     !! (n_points, n_ao)
      real(dp), intent(in) :: wv(:)        !! weight times potential per point
      real(dp), intent(out) :: v(:, :)
      real(dp), allocatable :: scaled(:, :)
      integer :: mu
      allocate (scaled(size(ao, 1), size(ao, 2)))
      do mu = 1, size(ao, 2)
         scaled(:, mu) = ao(:, mu)*wv
      end do
      call pic_gemm(ao, scaled, v, transa="T")
   end subroutine grid_matrix

   subroutine basis_on_one_atom(element_symbols, atom, shells, basis)
      !! A molecular basis with `shells` on `atom` and nothing anywhere else
      character(len=*), intent(in) :: element_symbols(:)
      integer, intent(in) :: atom
      type(molecular_basis_type), intent(in) :: shells   !! One element, the proton set
      type(molecular_basis_type), intent(out) :: basis
      integer :: i

      call basis%allocate_elements(size(element_symbols))
      do i = 1, size(element_symbols)
         if (i == atom) then
            basis%elements(i) = shells%elements(1)
         else
            basis%elements(i)%element = trim(element_symbols(i))
            call basis%elements(i)%allocate_shells(0)
         end if
      end do
      basis%angular_form = shells%angular_form
   end subroutine basis_on_one_atom

   subroutine inverse_sqrt(s, x, error)
      !! `X = S^(-1/2)`, symmetric, for a small well-conditioned overlap
      real(dp), intent(in) :: s(:, :)
      real(dp), intent(out) :: x(:, :)
      type(error_t), intent(inout) :: error
      real(dp), allocatable :: u(:, :), w(:), scaled(:, :)
      integer :: n, i, info

      n = size(s, 1)
      allocate (u(n, n), w(n), scaled(n, n))
      u = s
      call pic_syev(u, w, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "NEO: the proton overlap could not be diagonalised")
         return
      end if
      if (minval(w) <= 1.0e-10_dp) then
         call error%set(ERROR_VALIDATION, "NEO: the proton basis is linearly dependent")
         return
      end if
      do i = 1, n
         scaled(:, i) = u(:, i)/sqrt(w(i))
      end do
      call pic_gemm(scaled, u, x, transb="T")
   end subroutine inverse_sqrt

   subroutine lowest_state(f, x, c, eps, d, error)
      !! Diagonalise `F` in the orthogonal basis; the lowest orbital is the proton
      real(dp), intent(in) :: f(:, :), x(:, :)
      real(dp), intent(out) :: c(:, :)       !! All orbitals, back in the AO basis
      real(dp), intent(out) :: eps(:)
      real(dp), intent(out) :: d(:, :)       !! `c_1 c_1^T`
      type(error_t), intent(inout) :: error
      real(dp), allocatable :: work(:, :), f_ortho(:, :)
      integer :: n, info

      n = size(f, 1)
      allocate (work(n, n), f_ortho(n, n))
      call pic_gemm(f, x, work)
      call pic_gemm(x, work, f_ortho, transa="T")
      call pic_syev(f_ortho, eps, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "NEO: the proton Fock matrix could not be diagonalised")
         return
      end if
      call pic_gemm(x, f_ortho, c)
      call pic_gemm(c(:, 1:1), c(:, 1:1), d, transb="T")
   end subroutine lowest_state

   function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=16) :: s
      write (s, "(I0)") i
   end function itoa

end module mqc_czt_neo
