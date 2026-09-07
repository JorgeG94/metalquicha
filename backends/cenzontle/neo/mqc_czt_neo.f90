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

   integer, parameter :: MAX_MACRO = 100
      !! Macro-iterations before giving up; NEO-HF takes ten or twenty.

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
                             density_tol, verbose, result, error, force_cartesian, in_core)
      !! NEO-HF for a closed-shell electronic structure and any number of quantum protons
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
      integer :: natm, nq, nao_e, nao_p, nao_c, p, q, it, i, lo, hi, qlo, qhi
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

         ! The electrons, in that field. Warm-started after the first cycle.
         if (it == 1) then
            call run_czt_rhf(mol_e, nelec, max_iter, energy_tol, density_tol, .false., &
                             result%electrons, error, h_extra=h_extra, in_core=in_core)
         else
            call run_czt_rhf(mol_e, nelec, max_iter, energy_tol, density_tol, .false., &
                             result%electrons, error, h_extra=h_extra, in_core=in_core, &
                             guess_density=d_e_prev)
         end if
         if (error%has_error()) return
         if (.not. result%electrons%converged) then
            call error%set(ERROR_VALIDATION, "NEO: the electronic SCF did not converge "// &
                           "in macro-iteration "//trim(itoa(it)))
            return
         end if
         d_e_prev = result%electrons%density
         e_ep = sum(result%electrons%density*h_extra)
         e_total = result%electrons%energy + e_p1 + e_pp

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
         dd_max = 0.0_dp
         do p = 1, nq
            lo = offset(p) + 1
            hi = offset(p) + nao_p
            f_p = h_p(:, :, p) - j_e(lo:hi, lo:hi, 1)
            do q = 1, nq
               if (q == p) cycle
               qlo = offset(p) + 1
               qhi = offset(p) + nao_p
               f_p = f_p + j_p(qlo:qhi, qlo:qhi, q)
            end do
            call lowest_state(f_p, x_p(:, :, p), c_p(:, :, p), eps_p(:, p), work, error)
            if (error%has_error()) return
            dd = maxval(abs(work - d_p(:, :, p)))
            dd_max = max(dd_max, dd)
            d_p(:, :, p) = work
         end do

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
   end subroutine run_czt_neo_hf

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
