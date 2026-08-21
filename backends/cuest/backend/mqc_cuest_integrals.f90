!! Per-molecule cuEST integral objects and the matrices built from them
module mqc_cuest_integrals
   !! Holds everything cuEST needs for one molecule -- AO basis, auxiliary
   !! basis, pair list and integral plans -- and computes the matrices an SCF
   !! consumes: S, T, V, and the density-fitted J and K.
   !!
   !! Lifetime: one `cuest_system_t` per fragment. The cuEST *handle* is shared
   !! across fragments (see `mqc_cuest_context`), but every object below is
   !! geometry- and basis-specific and must be rebuilt per fragment.
   !!
   !! Two conventions that are invisible in the signatures and cost an
   !! afternoon each if got wrong:
   !!
   !!   * **Host vs device pointers.** Shell exponents/coefficients and the
   !!     pair list's coordinates are HOST arrays; every matrix in and out is a
   !!     DEVICE buffer. Both are `type(c_ptr)`, and confusing them segfaults
   !!     rather than returning a status.
   !!   * **Nuclear charges are stored as -Z.** cuEST does not assume the sign
   !!     of the electron charge, so the potential routine wants negated
   !!     nuclear charges.
   !!
   !! Matrix layout: cuEST matrices are row-major, Fortran arrays are
   !! column-major, so a Fortran `(n_ao, n_ao)` array is passed as its own
   !! transpose. Every matrix crossing this boundary (S, T, V, D, J, K) is
   !! symmetric, so that is a no-op. The one non-symmetric array, the occupied
   !! MO coefficients, is stored `(n_ao, n_occ)` here, whose column-major
   !! layout is byte-identical to the row-major `(n_occ, n_ao)` cuEST wants.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_int32_t, c_int64_t, &
                                                                             c_size_t, c_double, c_loc, c_associated
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_method_config, only: pcm_config_t
   use mqc_pcm_radii, only: cavity_radius
   use mqc_cgto, only: molecular_basis_type
   use mqc_cuest_basis, only: cuest_shell_set_t, build_cuest_shells
   use mqc_cuest_grid, only: atom_grid_set_t, build_atom_grids
   use mqc_cuest_context, only: cuest_context_t
   use mqc_cuest_runtime, only: cuest_status_check, cublas_status_check, &
                                workspace_alloc, workspace_free, &
                                copy_to_device, copy_to_host, device_sync, &
                                device_offset, &
                                device_alloc, device_free
   use cublas, only: cublasDaxpy, cublasDcopy, cublasDdot
   use cuest, only: cuestWorkspace_t, cuestWorkspaceDescriptor_t, &
                    cuestParametersCreate, cuestParametersDestroy, &
                    cuestPCMIntPlanCreate, cuestPCMIntPlanCreateWorkspaceQuery, &
                    cuestPCMIntPlanDestroy, &
                    cuestPCMPotentialCompute, cuestPCMPotentialComputeWorkspaceQuery, &
                    cuestPCMDerivativeCompute, cuestPCMDerivativeComputeWorkspaceQuery, &
                    cuestResultsCreate, cuestResultsDestroy, &
                    CUEST_PCMINTPLAN, CUEST_PCMINTPLAN_NUM_POINT, &
                    CUEST_PCMINTPLAN_NUM_ACTIVE_POINT, &
                    CUEST_PCMINTPLAN_PARAMETERS, &
                    CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS, &
                    CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS_CONVERGENCE_THRESHOLD, &
                    CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS_MAX_ITERATIONS, &
                    CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS, &
                    CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS_CONVERGENCE_THRESHOLD, &
                    CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS_MAX_ITERATIONS, &
                    CUEST_PCM_RESULTS, CUEST_PCMRESULT_PCM_DIELECTRIC_ENERGY, &
                    CUEST_PCMRESULT_CONVERGED_RESIDUAL, &
                    CUEST_PCMRESULT_NUM_ITERATIONS_TAKEN, CUEST_PCMRESULT_CONVERGED, &
                    cuestAOBasisCreate, cuestAOBasisCreateWorkspaceQuery, cuestAOBasisDestroy, &
                    cuestAOPairListCreate, cuestAOPairListCreateWorkspaceQuery, &
                    cuestAOPairListDestroy, &
                    cuestOEIntPlanCreate, cuestOEIntPlanCreateWorkspaceQuery, &
                    cuestOEIntPlanDestroy, &
                    cuestDFIntPlanCreate, cuestDFIntPlanCreateWorkspaceQuery, &
                    cuestDFIntPlanDestroy, &
                    cuestOverlapCompute, cuestOverlapComputeWorkspaceQuery, &
                    cuestKineticCompute, cuestKineticComputeWorkspaceQuery, &
                    cuestPotentialCompute, cuestPotentialComputeWorkspaceQuery, &
                    cuestMultipoleCompute, cuestMultipoleComputeWorkspaceQuery, &
                    CUEST_MULTIPOLECOMPUTE_PARAMETERS, &
                    cuestDFCoulombCompute, cuestDFCoulombComputeWorkspaceQuery, &
                    cuestDFSymmetricExchangeCompute, &
                    cuestDFSymmetricExchangeComputeWorkspaceQuery, &
                    CUEST_AOBASIS, CUEST_AOBASIS_NUM_AO, &
                    CUEST_AOBASIS_PARAMETERS, CUEST_AOPAIRLIST_PARAMETERS, &
                    CUEST_OEINTPLAN_PARAMETERS, CUEST_DFINTPLAN_PARAMETERS, &
                    CUEST_OVERLAPCOMPUTE_PARAMETERS, CUEST_KINETICCOMPUTE_PARAMETERS, &
                    CUEST_POTENTIALCOMPUTE_PARAMETERS, CUEST_DFCOULOMBCOMPUTE_PARAMETERS, &
                    CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS, &
                    CUEST_DFINTPLAN_PARAMETERS_EXCHANGE_FRACTION, &
                    CUEST_DFINTPLAN_PARAMETERS_LRC_EXCHANGE_FRACTION, &
                    CUEST_DFINTPLAN_PARAMETERS_LRC_EXCHANGE_OMEGA, &
                    cuestMolecularGridCreate, cuestMolecularGridCreateWorkspaceQuery, &
                    cuestMolecularGridDestroy, &
                    cuestXCIntPlanCreate, cuestXCIntPlanCreateWorkspaceQuery, &
                    cuestXCIntPlanDestroy, &
                    cuestXCPotentialRKSCompute, cuestXCPotentialRKSComputeWorkspaceQuery, &
                    cuestXCPotentialUKSCompute, cuestXCPotentialUKSComputeWorkspaceQuery, &
                    cuestXCDerivativeUKSCompute, cuestXCDerivativeUKSComputeWorkspaceQuery, &
                    CUEST_XCPOTENTIALUKSCOMPUTE_PARAMETERS, &
                    CUEST_XCDERIVATIVEUKSCOMPUTE_PARAMETERS, &
                    cuestOverlapDerivativeCompute, cuestOverlapDerivativeComputeWorkspaceQuery, &
                    cuestKineticDerivativeCompute, cuestKineticDerivativeComputeWorkspaceQuery, &
                    cuestPotentialDerivativeCompute, cuestPotentialDerivativeComputeWorkspaceQuery, &
                    cuestDFSymmetricDerivativeCompute, &
                    cuestDFSymmetricDerivativeComputeWorkspaceQuery, &
                    cuestXCDerivativeRKSCompute, cuestXCDerivativeRKSComputeWorkspaceQuery, &
                    CUEST_OVERLAPDERIVATIVECOMPUTE_PARAMETERS, &
                    CUEST_KINETICDERIVATIVECOMPUTE_PARAMETERS, &
                    CUEST_POTENTIALDERIVATIVECOMPUTE_PARAMETERS, &
                    CUEST_DFSYMMETRICDERIVATIVECOMPUTE_PARAMETERS, &
                    CUEST_XCDERIVATIVERKSCOMPUTE_PARAMETERS, &
                    CUEST_MOLECULARGRID_PARAMETERS, CUEST_XCINTPLAN_PARAMETERS, &
                    CUEST_XCPOTENTIALRKSCOMPUTE_PARAMETERS, &
                    CUEST_XCINTPLAN, CUEST_XCINTPLAN_EXCHANGE_SCALE, &
                    CUEST_XCINTPLAN_LRC_EXCHANGE_SCALE, CUEST_XCINTPLAN_LRC_OMEGA
   use cuest_helpers, only: cuest_query_i64, cuest_query_f64, cuest_param_set_f64, &
                            cuest_param_set_i64, cuest_results_query_f64, &
                            cuest_results_query_i64
   implicit none
   private

   public :: cuest_system_t  !! All cuEST objects for one molecule

   real(dp), parameter :: PAIR_LIST_THRESHOLD = 1.0e-14_dp
      !! Pair screening tolerance; the value cuEST's own samples use

   integer(c_size_t), parameter :: DF_EXCHANGE_BUFFER_BYTES = 2000000000_c_size_t
      !! Cap on DF-K intermediates. The algorithm exploits as much scratch as
      !! it is given; 2 GB is the reference samples' default.

   type :: cuest_system_t
      !! cuEST objects and device buffers for a single molecule
      type(c_ptr) :: handle = c_null_ptr  !! Borrowed handle, not owned here
      type(c_ptr) :: cublas = c_null_ptr  !! Borrowed cuBLAS handle, not owned here
      type(c_ptr) :: cusolver = c_null_ptr  !! Borrowed cuSOLVER handle, not owned here
      integer(c_int) :: solver_lwork = 0
         !! Doubles in `d_solver`, from a cuSOLVER bufferSize query

      ! cuEST objects, each with the persistent workspace that must outlive it
      type(c_ptr) :: basis = c_null_ptr       !! Primary (orbital) AO basis
      type(c_ptr) :: aux_basis = c_null_ptr   !! Auxiliary (fitting) AO basis
      type(c_ptr) :: pair_list = c_null_ptr   !! AO pair list, primary basis
      type(c_ptr) :: oe_plan = c_null_ptr     !! One-electron integral plan
      type(c_ptr) :: df_plan = c_null_ptr     !! Density-fitted integral plan
      type(c_ptr) :: molecular_grid = c_null_ptr  !! XC quadrature grid (DFT only)
      type(c_ptr) :: xc_plan = c_null_ptr     !! XC integral plan (DFT only)
      type(c_ptr) :: pcm_plan = c_null_ptr    !! Continuum plan (solvated only)
      type(cuestWorkspace_t) :: ws_basis, ws_aux_basis, ws_pair_list
      type(cuestWorkspace_t) :: ws_oe_plan, ws_df_plan
      type(cuestWorkspace_t) :: ws_grid, ws_xc_plan, ws_pcm_plan

      ! ---- polarizable continuum -------------------------------------------
      logical :: has_pcm = .false.
         !! Whether a continuum plan exists, and so whether `d_pcm` holds a
         !! potential this iteration. Read by `assemble_fock`, which skips the
         !! term rather than adding whatever the last fragment left behind.
      type(c_ptr) :: d_pcm = c_null_ptr
         !! The continuum's contribution to the Fock matrix, (n_ao, n_ao).
         !! Spin-independent -- the surface charges come from the total density --
         !! so an unrestricted iteration adds the same matrix to both channels,
         !! which is why this lives on the system rather than being passed in.
      type(c_ptr) :: d_q_in = c_null_ptr    !! Surface charges, previous iteration
      type(c_ptr) :: d_q_out = c_null_ptr   !! Surface charges, this iteration
         !! Two buffers rather than one because the charge solve is iterative and
         !! is handed the previous solution as its starting point. Aliasing input
         !! to output would be a bet on cuEST reading before it writes.
      integer(c_int64_t) :: n_pcm_points = 0      !! Cavity surface points
      integer(c_int64_t) :: n_pcm_active = 0      !! Of those, the ones not buried
      real(dp) :: pcm_tolerance = 1.0e-8_dp
      integer :: pcm_max_iter = 100
      ! Kept for reporting. A solvation energy is meaningless without the
      ! dielectric and the cavity that produced it, and none of these is
      ! recoverable from the output afterwards.
      real(dp) :: pcm_dielectric = 0.0_dp
      real(dp) :: pcm_zeta = 0.0_dp
      real(dp) :: pcm_radii_scale = 0.0_dp
      integer :: pcm_angular_points = 0
      real(dp) :: pcm_residual = 0.0_dp           !! Last solve's residual
      integer :: pcm_iterations = 0               !! Last solve's iteration count
      logical :: pcm_solved = .true.
         !! Whether the last charge solve converged. False is not fatal on its own
         !! -- an early SCF iteration need not solve tightly -- but the driver
         !! refuses a *converged* SCF whose final continuum did not solve.

      ! Dimensions
      integer(c_int64_t) :: n_atoms = 0  !! Atoms in this molecule
      integer(c_int64_t) :: n_ao = 0     !! AO basis functions
      integer(c_int64_t) :: n_occ = 0    !! Doubly occupied orbitals (RKS)
      integer(c_int64_t) :: n_occ_beta = 0
         !! Beta occupied orbitals. Unrestricted when this is >= 0 and
         !! `unrestricted` is set; n_occ then means the alpha count.
      logical :: unrestricted = .false.
         !! Whether the device buffers for a second spin channel exist

      real(dp) :: exchange_fraction = 1.0_dp
         !! Fraction of full-range exact exchange baked into the DF plan.
         !! 1.0 for Hartree-Fock, 0.25 for PBE0, 0 for a pure functional.
         !! For DFT this is queried from the XC plan, never hardcoded.
      real(dp) :: lrc_exchange_fraction = 0.0_dp
         !! Long-range exact exchange fraction, for range-separated hybrids
      real(dp) :: lrc_omega = 0.0_dp
         !! Range-separation parameter

      logical :: has_xc = .false.
         !! Whether an XC plan exists, i.e. this is DFT rather than HF
      logical :: needs_exchange = .true.
         !! False for a pure functional, where K would be computed and then
         !! scaled by zero. cuEST's own header advises skipping the call.

      ! Geometry, mirrored to the device for the potential integrals
      real(dp), allocatable :: xyz_host(:)      !! 3*n_atoms, Bohr, atom-major
      real(dp), allocatable :: charges_host(:)  !! n_atoms, stored as -Z

      ! Device buffers BORROWED from the context's pools, never owned here:
      ! `destroy` nulls them rather than freeing, so the next fragment reuses
      ! the same allocation instead of paying for another cudaMalloc.
      type(c_ptr) :: xyz_device = c_null_ptr
      type(c_ptr) :: charges_device = c_null_ptr
      type(c_ptr) :: d_matrix = c_null_ptr  !! Density
      type(c_ptr) :: d_c_occ = c_null_ptr   !! Occupied MO coefficients (alpha)
      type(c_ptr) :: d_c_occ_beta = c_null_ptr   !! Beta occupied MOs (UKS)
      type(c_ptr) :: d_result_beta = c_null_ptr  !! Second output matrix (UKS)
      type(c_ptr) :: d_result = c_null_ptr  !! Whichever matrix is being built
      type(c_ptr) :: d_gradient = c_null_ptr        !! natom x 3 gradient output
      type(c_ptr) :: d_charge_gradient = c_null_ptr  !! Hellmann-Feynman half

      ! One buffer per Fock term, so the three can coexist and be combined on
      ! the device instead of each being fetched before the next overwrites it.
      type(c_ptr) :: d_j = c_null_ptr     !! Coulomb output
      type(c_ptr) :: d_k = c_null_ptr     !! Exchange output
      type(c_ptr) :: d_xc = c_null_ptr    !! XC potential output
      type(c_ptr) :: d_fock = c_null_ptr  !! Assembled Fock
      type(c_ptr) :: d_core = c_null_ptr  !! Core Hamiltonian, iteration-invariant
      type(c_ptr) :: d_ovlp = c_null_ptr  !! Overlap, iteration-invariant
      type(c_ptr) :: d_error = c_null_ptr  !! DIIS error vector, (n_mo, n_mo)
      type(c_ptr) :: d_transform = c_null_ptr   !! Orthogonalizer X, (n_ao, n_mo)
      type(c_ptr) :: d_commutator = c_null_ptr  !! FDS - SDF in the AO basis
      type(c_ptr) :: d_work = c_null_ptr
         !! General scratch for the SCF's own linear algebra, n_ao^2. Every
         !! user must treat it as clobbered on entry and must not expect it to
         !! survive a call into anything else.

      ! Fock diagonalization
      type(c_ptr) :: d_fock_ortho = c_null_ptr    !! F' = X^T F X, then C'
      type(c_ptr) :: d_orbitals = c_null_ptr      !! C = X C', (n_ao, n_mo)
      type(c_ptr) :: d_eigenvalues = c_null_ptr   !! Orbital energies, (n_mo)
      type(c_ptr) :: d_solver = c_null_ptr        !! cuSOLVER workspace
      type(c_ptr) :: d_devinfo = c_null_ptr       !! cuSOLVER convergence flag

      ! Beta channel, null unless `unrestricted`. `d_fock_beta` is an OFFSET
      ! into the same allocation as `d_fock`, not a buffer of its own, so the
      ! two Fock matrices form one contiguous vector for the DIIS history to
      ! extrapolate over in a single call.
      type(c_ptr) :: d_density_alpha = c_null_ptr   !! D^a
      type(c_ptr) :: d_density_beta = c_null_ptr    !! D^b
      type(c_ptr) :: d_k_beta = c_null_ptr          !! K[C_b]
      type(c_ptr) :: d_xc_beta = c_null_ptr         !! Vxc_b
      type(c_ptr) :: d_fock_beta = c_null_ptr       !! F^b, at offset n_ao^2 of d_fock
      type(c_ptr) :: d_eigenvalues_beta = c_null_ptr  !! Beta orbital energies
   contains
      procedure :: create => system_create
      procedure :: destroy => system_destroy
      procedure :: compute_overlap => system_compute_overlap
      procedure :: compute_kinetic => system_compute_kinetic
      procedure :: compute_potential => system_compute_potential
      procedure :: compute_coulomb => system_compute_coulomb
      procedure :: compute_exchange => system_compute_exchange
      procedure :: compute_xc => system_compute_xc
      procedure :: compute_dipole => system_compute_dipole
      procedure :: compute_xc_uks => system_compute_xc_uks

      ! Device-resident SCF path: stage inputs, compute into a chosen output
      ! buffer without fetching, combine, and fetch once at the end.
      procedure :: stage_density => system_stage_density
      procedure :: stage_occupied => system_stage_occupied
      procedure :: stage_to => system_stage_to
      procedure :: coulomb_device => system_coulomb_device
      procedure :: exchange_device => system_exchange_device
      procedure :: xc_device => system_xc_device
      procedure :: pcm_device => system_pcm_device
      procedure :: build_pcm => build_pcm_plan
      procedure :: xc_uks_device => system_xc_uks_device
      procedure :: assemble_fock => system_assemble_fock
      procedure :: add_into => system_add_into
      procedure :: matrix_dot => system_matrix_dot
      procedure :: fetch => system_fetch
      procedure :: gradient_overlap => system_gradient_overlap
      procedure :: gradient_kinetic => system_gradient_kinetic
      procedure :: gradient_potential => system_gradient_potential
      procedure :: gradient_two_electron => system_gradient_two_electron
      procedure :: gradient_xc => system_gradient_xc
      procedure :: gradient_two_electron_uks => system_gradient_two_electron_uks
      procedure :: gradient_xc_uks => system_gradient_xc_uks
      procedure :: gradient_pcm => system_gradient_pcm
   end type cuest_system_t

contains

   subroutine system_create(this, context, atomic_numbers, coordinates, mol_basis, &
                            aux_mol_basis, use_spherical, n_occ, functional_id, &
                            n_radial, n_angular, error, n_occ_beta, n_guess_columns, pcm)
      !! Build every cuEST object needed for one molecule
      !!
      !! Coordinates are in Bohr, matching metalquicha's internal units and
      !! cuEST's expectation, so no conversion happens here.
      !!
      !! The handle and all plain device scratch come from `context` and are
      !! shared with every other fragment this rank evaluates.
      class(cuest_system_t), intent(out) :: this
      type(cuest_context_t), intent(inout) :: context          !! Per-rank handle and scratch pools
      integer, intent(in) :: atomic_numbers(:)                 !! Z per atom
      real(dp), intent(in) :: coordinates(:, :)                !! (3, n_atoms), Bohr
      type(molecular_basis_type), intent(in) :: mol_basis      !! Orbital basis, per atom
      type(molecular_basis_type), intent(in) :: aux_mol_basis  !! JKFIT basis, per atom
      logical, intent(in) :: use_spherical                     !! Pure vs Cartesian
      integer, intent(in) :: n_occ                             !! Doubly occupied orbitals
      integer, intent(in) :: functional_id                     !! cuEST XC functional, or <0 for pure HF
      type(pcm_config_t), intent(in), optional :: pcm
         !! A polarizable continuum, when one was asked for. Absent, or present
         !! with `enabled` false, leaves the calculation in the gas phase.
      integer, intent(in) :: n_radial                          !! Radial grid points (DFT only)
      integer, intent(in) :: n_angular                         !! Angular grid points (DFT only)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_occ_beta
         !! Beta occupied count. Present means unrestricted, and `n_occ` is
         !! then the alpha count rather than the doubly occupied count.
      integer, intent(in), optional :: n_guess_columns
         !! Widest coefficient matrix that will be passed in. A superposition
         !! guess carries the summed atomic occupations, which usually exceeds
         !! the molecule's own count. Sizing here matters: growing a pool later
         !! would reallocate it and invalidate the pointers already borrowed.

      integer :: iatom, n_atoms
      integer(c_int64_t) :: widest, n_spin

      n_atoms = size(atomic_numbers)
      this%handle = context%handle
      this%cublas = context%cublas_handle
      this%cusolver = context%cusolver_handle
      this%n_atoms = int(n_atoms, c_int64_t)
      this%n_occ = int(n_occ, c_int64_t)
      this%unrestricted = present(n_occ_beta)
      if (this%unrestricted) this%n_occ_beta = int(n_occ_beta, c_int64_t)
      this%has_xc = (functional_id >= 0)

      if (size(coordinates, 2) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "cuEST system: coordinate/atom count mismatch")
         return
      end if

      ! ---- geometry, host then device --------------------------------------
      allocate (this%xyz_host(3*n_atoms), this%charges_host(n_atoms))
      do iatom = 1, n_atoms
         this%xyz_host(3*(iatom - 1) + 1:3*iatom) = coordinates(1:3, iatom)
         ! Negated: cuEST does not assume the sign of the electron charge.
         this%charges_host(iatom) = -real(atomic_numbers(iatom), dp)
      end do

      call context%scratch_xyz%ensure(int(3*n_atoms, c_int64_t), "atom coordinates", error)
      call context%scratch_charges%ensure(int(n_atoms, c_int64_t), "nuclear charges", error)
      if (error%has_error()) return
      this%xyz_device = context%scratch_xyz%ptr
      this%charges_device = context%scratch_charges%ptr

      call copy_to_device(this%xyz_device, this%xyz_host, "atom coordinates", error)
      call copy_to_device(this%charges_device, this%charges_host, "nuclear charges", error)
      if (error%has_error()) return

      ! ---- bases, pair list, plans -----------------------------------------
      call build_basis(this, mol_basis, use_spherical, .false., error)
      if (error%has_error()) return

      call build_basis(this, aux_mol_basis, use_spherical, .true., error)
      if (error%has_error()) return

      call cuest_status_check(cuest_query_i64(this%handle, CUEST_AOBASIS, this%basis, &
                                              CUEST_AOBASIS_NUM_AO, this%n_ao), &
                              "query AOBASIS_NUM_AO", error)
      if (error%has_error()) return

      call build_pair_list(this, error)
      if (error%has_error()) return

      call build_oe_plan(this, error)
      if (error%has_error()) return

      ! The XC plan comes first for DFT: it is the authority on how much exact
      ! exchange the functional wants, and the DF plan has to be built with a
      ! matching operator. Querying rather than tabulating means a hybrid can
      ! never end up with the Coulomb and XC sides disagreeing.
      if (this%has_xc) then
         call build_xc_plan(this, atomic_numbers, functional_id, n_radial, n_angular, error)
         if (error%has_error()) return
      else
         this%exchange_fraction = 1.0_dp
         this%lrc_exchange_fraction = 0.0_dp
         this%lrc_omega = 0.0_dp
      end if

      this%needs_exchange = (this%exchange_fraction /= 0.0_dp) .or. &
                            (this%lrc_exchange_fraction /= 0.0_dp)

      ! The continuum, last: it is built on the one-electron plan and is
      ! independent of everything above it.
      if (present(pcm)) then
         if (pcm%enabled) then
            call this%build_pcm(atomic_numbers, pcm, error)
            if (error%has_error()) return
         end if
      end if

      call build_df_plan(this, error)
      if (error%has_error()) return

      ! ---- SCF scratch, borrowed from the rank's pools ---------------------
      call context%scratch_density%ensure(this%n_ao*this%n_ao, "density matrix", error)
      call context%scratch_result%ensure(this%n_ao*this%n_ao, "result matrix", error)
      if (error%has_error()) then
         call error%add_context("mqc_cuest_integrals:create")
         return
      end if
      this%d_matrix = context%scratch_density%ptr
      this%d_result = context%scratch_result%ptr

      ! Everything the device-resident SCF works on: the Fock terms, the
      ! orthogonalizer, and the commutator and its scratch. Ten n_ao^2 buffers
      ! is 80 MB at n_ao = 1000, against the 40-80 GB of an A100 or H100 --
      ! cheap next to the six transfers and three synchronises per iteration it
      ! removes, and pooled, so a fragmented run pays for it once per rank.
      !
      ! Several are sized n_ao^2 when they only need n_ao*n_mo or n_mo^2, since
      ! n_mo <= n_ao always. That costs a little memory and no correctness.
      !
      ! The same goes for the orthogonalizer and the commutator scratch: X is
      ! (n_ao, n_mo) and the projected error (n_mo, n_mo), but n_mo comes out
      ! of the overlap diagonalization and is not known this early.
      call context%scratch_j%ensure(this%n_ao*this%n_ao, "Coulomb matrix", error)
      call context%scratch_k%ensure(this%n_ao*this%n_ao, "exchange matrix", error)
      call context%scratch_xc%ensure(this%n_ao*this%n_ao, "XC potential", error)
      ! Two spin channels mean two Fock matrices and two error vectors, and the
      ! DIIS drives both from one stacked vector, so these two pools carry both
      ! back to back rather than being paired with buffers of their own.
      n_spin = merge(2_c_int64_t, 1_c_int64_t, this%unrestricted)
      call context%scratch_fock%ensure(n_spin*this%n_ao*this%n_ao, "Fock matrix", error)
      call context%scratch_core%ensure(this%n_ao*this%n_ao, "core Hamiltonian", error)
      call context%scratch_ovlp%ensure(this%n_ao*this%n_ao, "overlap matrix", error)
      call context%scratch_error%ensure(n_spin*this%n_ao*this%n_ao, "DIIS error vector", error)
      call context%scratch_transform%ensure(this%n_ao*this%n_ao, "orthogonalizer", error)
      call context%scratch_commutator%ensure(this%n_ao*this%n_ao, "commutator", error)
      call context%scratch_work%ensure(this%n_ao*this%n_ao, "commutator scratch", error)
      call context%scratch_fock_ortho%ensure(this%n_ao*this%n_ao, "orthogonal Fock", error)
      call context%scratch_orbitals%ensure(this%n_ao*this%n_ao, "molecular orbitals", error)
      call context%scratch_eigenvalues%ensure(this%n_ao, "orbital energies", error)
      call context%scratch_devinfo%ensure(1_c_int64_t, "eigensolver info", error)
      if (error%has_error()) then
         call error%add_context("mqc_cuest_integrals:create")
         return
      end if
      this%d_j = context%scratch_j%ptr
      this%d_k = context%scratch_k%ptr
      this%d_xc = context%scratch_xc%ptr
      this%d_fock = context%scratch_fock%ptr
      this%d_core = context%scratch_core%ptr
      this%d_ovlp = context%scratch_ovlp%ptr
      this%d_error = context%scratch_error%ptr
      this%d_transform = context%scratch_transform%ptr
      this%d_commutator = context%scratch_commutator%ptr
      this%d_work = context%scratch_work%ptr
      this%d_fock_ortho = context%scratch_fock_ortho%ptr
      this%d_orbitals = context%scratch_orbitals%ptr
      this%d_eigenvalues = context%scratch_eigenvalues%ptr
      this%d_devinfo = context%scratch_devinfo%ptr
      ! `d_solver` is sized by a cuSOLVER bufferSize query, which needs n_mo --
      ! not known until the overlap has been diagonalized. The SCF fills it in.
      this%d_solver = c_null_ptr

      call context%scratch_gradient%ensure(3*this%n_atoms, "gradient", error)
      call context%scratch_charge_gradient%ensure(3*this%n_atoms, "charge gradient", error)
      if (error%has_error()) then
         call error%add_context("mqc_cuest_integrals:create")
         return
      end if
      this%d_gradient = context%scratch_gradient%ptr
      this%d_charge_gradient = context%scratch_charge_gradient%ptr

      if (this%n_occ > 0) then
         widest = this%n_occ
         if (present(n_guess_columns)) widest = max(widest, int(n_guess_columns, c_int64_t))
         call context%scratch_c_occ%ensure(this%n_ao*widest, &
                                           "occupied MO coefficients", error)
         if (error%has_error()) then
            call error%add_context("mqc_cuest_integrals:create")
            return
         end if
         this%d_c_occ = context%scratch_c_occ%ptr
      end if

      if (this%unrestricted) then
         ! The DF derivative wants the two coefficient matrices concatenated,
         ! so the beta buffer is sized to hold both back to back.
         widest = this%n_occ + this%n_occ_beta
         if (present(n_guess_columns)) widest = max(widest, 2_c_int64_t*int(n_guess_columns, c_int64_t))
         call context%scratch_c_occ_beta%ensure(this%n_ao*widest, &
                                                "beta occupied MO coefficients", error)
         call context%scratch_result_beta%ensure(this%n_ao*this%n_ao, "beta result matrix", error)
         if (error%has_error()) then
            call error%add_context("mqc_cuest_integrals:create")
            return
         end if
         this%d_c_occ_beta = context%scratch_c_occ_beta%ptr
         this%d_result_beta = context%scratch_result_beta%ptr

         ! Beta twins for the device-resident SCF. The alpha channel reuses the
         ! buffers above, so a restricted run never allocates any of these.
         call context%scratch_density_alpha%ensure(this%n_ao*this%n_ao, "alpha density", error)
         call context%scratch_density_beta%ensure(this%n_ao*this%n_ao, "beta density", error)
         call context%scratch_k_beta%ensure(this%n_ao*this%n_ao, "beta exchange", error)
         call context%scratch_xc_beta%ensure(this%n_ao*this%n_ao, "beta XC potential", error)
         call context%scratch_eigenvalues_beta%ensure(this%n_ao, "beta orbital energies", error)
         if (error%has_error()) then
            call error%add_context("mqc_cuest_integrals:create")
            return
         end if
         this%d_density_alpha = context%scratch_density_alpha%ptr
         this%d_density_beta = context%scratch_density_beta%ptr
         this%d_k_beta = context%scratch_k_beta%ptr
         this%d_xc_beta = context%scratch_xc_beta%ptr
         this%d_eigenvalues_beta = context%scratch_eigenvalues_beta%ptr

         ! Not a buffer of its own: the second half of the Fock allocation, so
         ! the two spins are one contiguous vector for the DIIS.
         this%d_fock_beta = device_offset(this%d_fock, this%n_ao*this%n_ao)
      end if
   end subroutine system_create

   subroutine build_basis(this, mol_basis, use_spherical, is_auxiliary, error)
      !! Create one AO basis (primary or auxiliary) from a molecular basis
      class(cuest_system_t), intent(inout) :: this
      type(molecular_basis_type), intent(in) :: mol_basis
      logical, intent(in) :: use_spherical
      logical, intent(in) :: is_auxiliary  !! Which slot to fill
      type(error_t), intent(inout) :: error

      type(cuest_shell_set_t) :: shell_set
      type(cuestWorkspaceDescriptor_t) :: persistent_desc, temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: basis_params, basis_handle
      integer(c_int) :: status
      character(len=16) :: label

      label = merge("auxiliary       ", "primary         ", is_auxiliary)

      call build_cuest_shells(this%handle, mol_basis, use_spherical, shell_set, error)
      if (error%has_error()) return

      basis_params = c_null_ptr
      basis_handle = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_AOBASIS_PARAMETERS, basis_params), &
                              "cuestParametersCreate(AO basis)", error)
      if (error%has_error()) then
         call shell_set%destroy()
         return
      end if

      call cuest_status_check(cuestAOBasisCreateWorkspaceQuery(this%handle, shell_set%n_atoms, &
                                                               shell_set%n_shells_per_atom, &
                                                               shell_set%shells, basis_params, &
                                                               persistent_desc, temporary_desc, &
                                                               basis_handle), &
                              "cuestAOBasisCreateWorkspaceQuery("//trim(label)//")", error)
      if (.not. error%has_error()) then
         if (is_auxiliary) then
            call workspace_alloc(this%ws_aux_basis, persistent_desc, error)
         else
            call workspace_alloc(this%ws_basis, persistent_desc, error)
         end if
         call workspace_alloc(temporary_ws, temporary_desc, error)
      end if

      if (.not. error%has_error()) then
         if (is_auxiliary) then
            call cuest_status_check(cuestAOBasisCreate(this%handle, shell_set%n_atoms, &
                                                       shell_set%n_shells_per_atom, &
                                                       shell_set%shells, basis_params, &
                                                       this%ws_aux_basis, temporary_ws, &
                                                       this%aux_basis), &
                                    "cuestAOBasisCreate(auxiliary)", error)
         else
            call cuest_status_check(cuestAOBasisCreate(this%handle, shell_set%n_atoms, &
                                                       shell_set%n_shells_per_atom, &
                                                       shell_set%shells, basis_params, &
                                                       this%ws_basis, temporary_ws, &
                                                       this%basis), &
                                    "cuestAOBasisCreate(primary)", error)
         end if
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_AOBASIS_PARAMETERS, basis_params)
      call cuest_status_check(status, "cuestParametersDestroy(AO basis)", error)

      ! Shells are copied into the basis, so they can go now.
      call shell_set%destroy()
   end subroutine build_basis

   subroutine build_pair_list(this, error)
      !! Create the AO pair list over the primary basis
      !!
      !! Takes HOST coordinates, hence the C_LOC on a TARGET local copy.
      class(cuest_system_t), intent(inout) :: this
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: persistent_desc, temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: pair_list_params
      real(dp), allocatable, target :: xyz(:)
      integer(c_int) :: status

      allocate (xyz(size(this%xyz_host)))
      xyz = this%xyz_host

      pair_list_params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_AOPAIRLIST_PARAMETERS, pair_list_params), &
                              "cuestParametersCreate(pair list)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestAOPairListCreateWorkspaceQuery(this%handle, this%basis, &
                                                                  this%n_atoms, c_loc(xyz), &
                                                                  PAIR_LIST_THRESHOLD, &
                                                                  pair_list_params, &
                                                                  persistent_desc, temporary_desc, &
                                                                  this%pair_list), &
                              "cuestAOPairListCreateWorkspaceQuery", error)
      if (.not. error%has_error()) then
         call workspace_alloc(this%ws_pair_list, persistent_desc, error)
         call workspace_alloc(temporary_ws, temporary_desc, error)
      end if

      if (.not. error%has_error()) then
         call cuest_status_check(cuestAOPairListCreate(this%handle, this%basis, this%n_atoms, &
                                                       c_loc(xyz), PAIR_LIST_THRESHOLD, &
                                                       pair_list_params, this%ws_pair_list, &
                                                       temporary_ws, this%pair_list), &
                                 "cuestAOPairListCreate", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_AOPAIRLIST_PARAMETERS, pair_list_params)
      call cuest_status_check(status, "cuestParametersDestroy(pair list)", error)
   end subroutine build_pair_list

   subroutine build_oe_plan(this, error)
      !! Create the one-electron integral plan
      class(cuest_system_t), intent(inout) :: this
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: persistent_desc, temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: plan_params
      integer(c_int) :: status

      plan_params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_OEINTPLAN_PARAMETERS, plan_params), &
                              "cuestParametersCreate(OE plan)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestOEIntPlanCreateWorkspaceQuery(this%handle, this%basis, &
                                                                 this%pair_list, plan_params, &
                                                                 persistent_desc, temporary_desc, &
                                                                 this%oe_plan), &
                              "cuestOEIntPlanCreateWorkspaceQuery", error)
      if (.not. error%has_error()) then
         call workspace_alloc(this%ws_oe_plan, persistent_desc, error)
         call workspace_alloc(temporary_ws, temporary_desc, error)
      end if

      if (.not. error%has_error()) then
         call cuest_status_check(cuestOEIntPlanCreate(this%handle, this%basis, this%pair_list, &
                                                      plan_params, this%ws_oe_plan, temporary_ws, &
                                                      this%oe_plan), &
                                 "cuestOEIntPlanCreate", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_OEINTPLAN_PARAMETERS, plan_params)
      call cuest_status_check(status, "cuestParametersDestroy(OE plan)", error)
   end subroutine build_oe_plan

   subroutine build_df_plan(this, error)
      !! Create the density-fitted integral plan over both bases
      class(cuest_system_t), intent(inout) :: this
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: persistent_desc, temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: plan_params
      integer(c_int) :: status

      plan_params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_DFINTPLAN_PARAMETERS, plan_params), &
                              "cuestParametersCreate(DF plan)", error)
      if (error%has_error()) return

      ! The exchange operator is baked into the plan, not into the K call, so
      ! it must be set before the plan is built. Hartree-Fock wants the full
      ! 1.0; a hybrid functional would set its own fraction here.
      call cuest_status_check(cuest_param_set_f64(CUEST_DFINTPLAN_PARAMETERS, plan_params, &
                                                  CUEST_DFINTPLAN_PARAMETERS_EXCHANGE_FRACTION, &
                                                  this%exchange_fraction), &
                              "configure DF plan exchange fraction", error)
      call cuest_status_check(cuest_param_set_f64(CUEST_DFINTPLAN_PARAMETERS, plan_params, &
                                                  CUEST_DFINTPLAN_PARAMETERS_LRC_EXCHANGE_FRACTION, &
                                                  this%lrc_exchange_fraction), &
                              "configure DF plan long-range exchange fraction", error)
      call cuest_status_check(cuest_param_set_f64(CUEST_DFINTPLAN_PARAMETERS, plan_params, &
                                                  CUEST_DFINTPLAN_PARAMETERS_LRC_EXCHANGE_OMEGA, &
                                                  this%lrc_omega), &
                              "configure DF plan range-separation parameter", error)
      if (error%has_error()) return

      call cuest_status_check(cuestDFIntPlanCreateWorkspaceQuery(this%handle, this%basis, &
                                                                 this%aux_basis, this%pair_list, &
                                                                 plan_params, persistent_desc, &
                                                                 temporary_desc, this%df_plan), &
                              "cuestDFIntPlanCreateWorkspaceQuery", error)
      if (.not. error%has_error()) then
         call workspace_alloc(this%ws_df_plan, persistent_desc, error)
         call workspace_alloc(temporary_ws, temporary_desc, error)
      end if

      if (.not. error%has_error()) then
         call cuest_status_check(cuestDFIntPlanCreate(this%handle, this%basis, this%aux_basis, &
                                                      this%pair_list, plan_params, this%ws_df_plan, &
                                                      temporary_ws, this%df_plan), &
                                 "cuestDFIntPlanCreate", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_DFINTPLAN_PARAMETERS, plan_params)
      call cuest_status_check(status, "cuestParametersDestroy(DF plan)", error)
   end subroutine build_df_plan

   subroutine system_xc_uks_device(this, d_c_alpha, d_c_beta, n_occ_beta, &
                                   d_out_alpha, d_out_beta, xc_energy, error)
      !! Spin-resolved XC potentials, device in and device out
      !!
      !! `n_occ_beta` is passed rather than read from the system because the
      !! atomic guess drives this with coefficient blocks whose width is not
      !! the molecule's own occupation.
      !!
      !! cuEST requires numOccupiedBeta > 0, but a system can legitimately have
      !! no beta electrons -- a hydrogen atom, which is exactly what an atomic
      !! guess has to solve. The caller passes 1 and guarantees that column is
      !! zeroed, which gives a beta density of zero: the correct physics through
      !! an argument the library will accept.
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: d_c_alpha, d_c_beta   !! Occupied MOs on device
      integer, intent(in) :: n_occ_beta                !! Beta columns, at least 1
      type(c_ptr), intent(in) :: d_out_alpha, d_out_beta  !! Vxc per spin, on device
      real(dp), intent(out) :: xc_energy               !! Exc, Hartree
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      integer(c_int64_t) :: beta_occupancy
      real(dp), target :: energy

      xc_energy = 0.0_dp
      if (error%has_error() .or. .not. this%has_xc) return

      beta_occupancy = max(int(n_occ_beta, c_int64_t), 1_c_int64_t)

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_XCPOTENTIALUKSCOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(UKS XC potential)", error)
      if (error%has_error()) return

      variable_buffer%hostBufferSizeInBytes = 0_c_size_t
      variable_buffer%deviceBufferSizeInBytes = DF_EXCHANGE_BUFFER_BYTES
      energy = 0.0_dp

      call cuest_status_check(cuestXCPotentialUKSComputeWorkspaceQuery(this%handle, this%xc_plan, &
                                                                       params, variable_buffer, &
                                                                       temporary_desc, this%n_occ, &
                                                                       beta_occupancy, d_c_alpha, &
                                                                       d_c_beta, c_loc(energy), &
                                                                       d_out_alpha, d_out_beta), &
                              "cuestXCPotentialUKSComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestXCPotentialUKSCompute(this%handle, this%xc_plan, params, &
                                                            variable_buffer, temporary_ws, &
                                                            this%n_occ, beta_occupancy, &
                                                            d_c_alpha, d_c_beta, &
                                                            c_loc(energy), d_out_alpha, &
                                                            d_out_beta), &
                                 "cuestXCPotentialUKSCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_XCPOTENTIALUKSCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(UKS XC potential)", error)

      if (.not. error%has_error()) xc_energy = energy
   end subroutine system_xc_uks_device

   subroutine system_compute_xc_uks(this, c_occ_alpha, c_occ_beta, xc_energy, &
                                    vxc_alpha, vxc_beta, error)
      !! Spin-resolved exchange-correlation energy and potentials
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: c_occ_alpha(:, :)   !! (n_ao, n_occ_alpha)
      real(dp), intent(in) :: c_occ_beta(:, :)    !! (n_ao, n_occ_beta)
      real(dp), intent(out) :: xc_energy          !! Exc, Hartree
      real(dp), intent(out) :: vxc_alpha(:, :)    !! (n_ao, n_ao)
      real(dp), intent(out) :: vxc_beta(:, :)     !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error

      integer(c_int64_t) :: beta_occupancy
      real(dp), allocatable :: flat(:)

      xc_energy = 0.0_dp
      vxc_alpha = 0.0_dp
      vxc_beta = 0.0_dp
      if (error%has_error() .or. .not. this%has_xc) return

      allocate (flat(size(c_occ_alpha)))
      flat = reshape(c_occ_alpha, [size(c_occ_alpha)])
      call copy_to_device(this%d_c_occ, flat, "alpha occupied MOs (XC)", error)
      deallocate (flat)

      ! The zero column an empty beta channel needs is materialized here, on
      ! the host side of the boundary; the device routine only promises to pass
      ! the count on.
      beta_occupancy = max(this%n_occ_beta, 1_c_int64_t)
      allocate (flat(this%n_ao*beta_occupancy))
      flat = 0.0_dp
      if (this%n_occ_beta > 0) flat = reshape(c_occ_beta, [size(c_occ_beta)])
      call copy_to_device(this%d_c_occ_beta, flat, "beta occupied MOs (XC)", error)
      deallocate (flat)
      if (error%has_error()) return

      call this%xc_uks_device(this%d_c_occ, this%d_c_occ_beta, int(beta_occupancy), &
                              this%d_result, this%d_result_beta, xc_energy, error)

      call fetch_matrix(this, vxc_alpha, "Vxc alpha", error)
      call this%fetch(this%d_result_beta, vxc_beta, "Vxc beta", error)
      if (error%has_error()) xc_energy = 0.0_dp
   end subroutine system_compute_xc_uks

   subroutine system_fetch(this, device_ptr, matrix, label, error)
      !! Synchronize, then copy an n_ao x n_ao device buffer into a host matrix
      !!
      !! The synchronize is what makes this safe to call on a buffer some cuEST
      !! or cuBLAS call has only *queued* work into. Every fetch in this module
      !! goes through here for that reason.
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: device_ptr
      real(dp), intent(out) :: matrix(:, :)
      character(len=*), intent(in) :: label
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat(:)

      matrix = 0.0_dp
      if (error%has_error()) return

      call device_sync(label, error)
      if (error%has_error()) return

      allocate (flat(this%n_ao*this%n_ao))
      call copy_to_host(flat, device_ptr, label, error)
      if (.not. error%has_error()) matrix = reshape(flat, [int(this%n_ao), int(this%n_ao)])
      deallocate (flat)
   end subroutine system_fetch

   subroutine build_xc_plan(this, atomic_numbers, functional_id, n_radial, n_angular, error)
      !! Build the molecular quadrature grid and the XC integral plan
      !!
      !! Afterwards the plan is interrogated for the exact-exchange content of
      !! the functional, which the caller feeds into the DF plan.
      class(cuest_system_t), intent(inout) :: this
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: functional_id
      integer, intent(in) :: n_radial, n_angular
      type(error_t), intent(inout) :: error

      type(atom_grid_set_t) :: grid_set
      type(cuestWorkspaceDescriptor_t) :: persistent_desc, temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: grid_params, plan_params
      real(dp), allocatable, target :: xyz(:)
      integer(c_int) :: status

      ! ---- per-atom grids ---------------------------------------------------
      call build_atom_grids(this%handle, atomic_numbers, n_radial, n_angular, &
                            grid_set, error)
      if (error%has_error()) return

      ! ---- molecular grid (HOST coordinates and HOST atom-grid array) -------
      allocate (xyz(size(this%xyz_host)))
      xyz = this%xyz_host

      grid_params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_MOLECULARGRID_PARAMETERS, grid_params), &
                              "cuestParametersCreate(molecular grid)", error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestMolecularGridCreateWorkspaceQuery(this%handle, &
                                                                        this%n_atoms, &
                                                                        grid_set%grids, c_loc(xyz), &
                                                                        grid_params, persistent_desc, &
                                                                        temporary_desc, &
                                                                        this%molecular_grid), &
                                 "cuestMolecularGridCreateWorkspaceQuery", error)
      end if
      if (.not. error%has_error()) then
         call workspace_alloc(this%ws_grid, persistent_desc, error)
         call workspace_alloc(temporary_ws, temporary_desc, error)
      end if
      if (.not. error%has_error()) then
         call cuest_status_check(cuestMolecularGridCreate(this%handle, this%n_atoms, &
                                                          grid_set%grids, c_loc(xyz), &
                                                          grid_params, this%ws_grid, &
                                                          temporary_ws, this%molecular_grid), &
                                 "cuestMolecularGridCreate", error)
      end if
      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_MOLECULARGRID_PARAMETERS, grid_params)
      call cuest_status_check(status, "cuestParametersDestroy(molecular grid)", error)

      ! The atom grids are copied into the molecular grid and are dead weight now.
      call grid_set%destroy()
      deallocate (xyz)
      if (error%has_error()) return

      ! ---- XC integral plan -------------------------------------------------
      plan_params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_XCINTPLAN_PARAMETERS, plan_params), &
                              "cuestParametersCreate(XC plan)", error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestXCIntPlanCreateWorkspaceQuery(this%handle, this%basis, &
                                                                    this%molecular_grid, &
                                                                    functional_id, plan_params, &
                                                                    persistent_desc, temporary_desc, &
                                                                    this%xc_plan), &
                                 "cuestXCIntPlanCreateWorkspaceQuery", error)
      end if
      if (.not. error%has_error()) then
         call workspace_alloc(this%ws_xc_plan, persistent_desc, error)
         call workspace_alloc(temporary_ws, temporary_desc, error)
      end if
      if (.not. error%has_error()) then
         call cuest_status_check(cuestXCIntPlanCreate(this%handle, this%basis, &
                                                      this%molecular_grid, functional_id, &
                                                      plan_params, this%ws_xc_plan, &
                                                      temporary_ws, this%xc_plan), &
                                 "cuestXCIntPlanCreate", error)
      end if
      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_XCINTPLAN_PARAMETERS, plan_params)
      call cuest_status_check(status, "cuestParametersDestroy(XC plan)", error)
      if (error%has_error()) return

      ! ---- ask the functional how much exact exchange it wants --------------
      call cuest_status_check(cuest_query_f64(this%handle, CUEST_XCINTPLAN, this%xc_plan, &
                                              CUEST_XCINTPLAN_EXCHANGE_SCALE, &
                                              this%exchange_fraction), &
                              "query XC plan exchange scale", error)
      call cuest_status_check(cuest_query_f64(this%handle, CUEST_XCINTPLAN, this%xc_plan, &
                                              CUEST_XCINTPLAN_LRC_EXCHANGE_SCALE, &
                                              this%lrc_exchange_fraction), &
                              "query XC plan long-range exchange scale", error)
      call cuest_status_check(cuest_query_f64(this%handle, CUEST_XCINTPLAN, this%xc_plan, &
                                              CUEST_XCINTPLAN_LRC_OMEGA, this%lrc_omega), &
                              "query XC plan range-separation parameter", error)
   end subroutine build_xc_plan

   subroutine build_pcm_plan(this, atomic_numbers, pcm, error)
      !! Build the continuum cavity and the PCM integral plan
      !!
      !! cuEST owns the hard part. Given the angular point count per atom, the
      !! cavity radii and the Gaussian switching exponents it builds the surface
      !! itself, so there is no tesserae construction here and no switching
      !! function to get wrong -- which is the half of a PCM implementation that
      !! is fiddly on the CPU. What is ours is the *input data*: the radii, which
      !! come from `mqc_pcm_radii` with a citation, and the exponents below.
      !!
      !! Done once per geometry. The plan depends on the cavity and the
      !! dielectric, neither of which moves during an SCF, so the iteration only
      !! ever calls `pcm_device`.
      class(cuest_system_t), intent(inout) :: this
      integer, intent(in) :: atomic_numbers(:)
      type(pcm_config_t), intent(in) :: pcm
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: persistent_desc, temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: plan_params
      integer(c_int64_t), allocatable :: angular_per_atom(:)
      real(dp), allocatable, target :: radii(:), zetas(:), charges(:)
         !! HOST arrays, and `target` so they can be passed by address. cuEST
         !! reads plan-creation data on the host -- as `cuestMolecularGridCreate`
         !! and `cuestAOBasisCreate` do with their coordinates -- so these must
         !! not be device buffers. They were, and cuEST segfaulted inside
         !! `cuestPCMIntPlanCreate` dereferencing a device address on the host.
      integer(c_int) :: status
      integer :: iatom, natm
      real(dp) :: radius

      if (error%has_error()) return

      ! A dielectric is required rather than defaulted: every solvent has a
      ! different one, and a default would silently pick a solvent.
      if (pcm%dielectric <= 1.0_dp) then
         call error%set(ERROR_VALIDATION, "the polarizable continuum needs a solvent "// &
                        "dielectric greater than one; set keywords.pcm.dielectric. "// &
                        "There is no solvent-name table on this path, so nothing can "// &
                        "be assumed from a name.")
         return
      end if
      if (pcm%angular_points <= 0) then
         call error%set(ERROR_VALIDATION, "keywords.pcm.angular_points must be positive")
         return
      end if
      ! The model is cuEST's own: the plan takes a cavity and a dielectric and
      ! solves its fixed formulation. "cpcm" -- the default -- is the only value
      ! accepted here; "iefpcm" exists for the CPU backend, which solves either,
      ! and substituting one model's charges for the other's would change the
      ! energy without saying so.
      if (trim(pcm%method) /= "cpcm") then
         call error%set(ERROR_VALIDATION, "keywords.pcm.method = '"//trim(pcm%method)// &
                        "' is not available on the cuEST backend, whose continuum "// &
                        "solver is fixed. Run this model on the CPU backend, or drop "// &
                        "the method keyword.")
         return
      end if

      natm = size(atomic_numbers)
      allocate (angular_per_atom(natm), radii(natm), zetas(natm), charges(natm))

      do iatom = 1, natm
         call cavity_radius(atomic_numbers(iatom), pcm%radii_scale, radius, error)
         if (error%has_error()) then
            deallocate (angular_per_atom, radii, zetas, charges)
            return
         end if
         radii(iatom) = radius
         angular_per_atom(iatom) = int(pcm%angular_points, c_int64_t)

         ! The Gaussian that replaces each surface point on a smooth cavity.
         ! Sharper for a small sphere and for a denser grid, since the exponent
         ! has to scale with the point spacing -- which on a sphere of radius R
         ! carrying n points goes as R/sqrt(n).
         !
         ! **The convention here is unverified against cuEST's.** It is the
         ! Lange-Herbert form, and cuEST takes exponents rather than the
         ! prefactor, so a different definition on its side smooths the cavity by
         ! the wrong amount and changes the solvation energy without failing.
         ! `keywords.pcm.zeta` exists so that can be tested on hardware.
         zetas(iatom) = pcm%zeta*sqrt(real(pcm%angular_points, dp))/radii(iatom)

         ! Effective nuclear charge. Equal to Z because no ECP is wired up on
         ! this backend; with one, the core electrons it replaces would have to
         ! come off here or the cavity would see the wrong nuclear potential.
         charges(iatom) = real(atomic_numbers(iatom), dp)
      end do

      plan_params = c_null_ptr
      if (.not. error%has_error()) then
         call cuest_status_check(cuestParametersCreate(CUEST_PCMINTPLAN_PARAMETERS, plan_params), &
                                 "cuestParametersCreate(PCM plan)", error)
      end if
      if (.not. error%has_error()) then
         call cuest_status_check(cuestPCMIntPlanCreateWorkspaceQuery(this%handle, this%oe_plan, &
                                                                     plan_params, persistent_desc, &
                                                                     temporary_desc, &
                                                                     angular_per_atom, &
                                                                     pcm%dielectric, &
                                                                     c_loc(zetas), c_loc(radii), &
                                                                     c_loc(charges), &
                                                                     this%pcm_plan), &
                                 "cuestPCMIntPlanCreateWorkspaceQuery", error)
      end if
      if (.not. error%has_error()) then
         call workspace_alloc(this%ws_pcm_plan, persistent_desc, error)
         call workspace_alloc(temporary_ws, temporary_desc, error)
      end if
      if (.not. error%has_error()) then
         call cuest_status_check(cuestPCMIntPlanCreate(this%handle, this%oe_plan, plan_params, &
                                                       this%ws_pcm_plan, temporary_ws, &
                                                       angular_per_atom, pcm%dielectric, &
                                                       c_loc(zetas), c_loc(radii), &
                                                       c_loc(charges), &
                                                       this%pcm_plan), &
                                 "cuestPCMIntPlanCreate", error)
      end if
      call workspace_free(temporary_ws)
      if (c_associated(plan_params)) then
         status = cuestParametersDestroy(CUEST_PCMINTPLAN_PARAMETERS, plan_params)
         call cuest_status_check(status, "cuestParametersDestroy(PCM plan)", error)
      end if

      ! Copied into the plan by now, as the molecular grid copies its coordinates.
      deallocate (angular_per_atom, radii, zetas, charges)
      if (error%has_error()) return

      ! How many surface points the cavity ended up with, which sizes the charge
      ! buffers. Asked rather than computed: the buried points are dropped by
      ! cuEST, so only it knows the answer.
      call cuest_status_check(cuest_query_i64(this%handle, CUEST_PCMINTPLAN, this%pcm_plan, &
                                              CUEST_PCMINTPLAN_NUM_POINT, this%n_pcm_points), &
                              "query PCM plan point count", error)
      call cuest_status_check(cuest_query_i64(this%handle, CUEST_PCMINTPLAN, this%pcm_plan, &
                                              CUEST_PCMINTPLAN_NUM_ACTIVE_POINT, &
                                              this%n_pcm_active), &
                              "query PCM plan active point count", error)
      if (error%has_error()) return

      call device_alloc(this%d_pcm, this%n_ao*this%n_ao, "PCM potential matrix", error)
      call device_alloc(this%d_q_in, this%n_pcm_points, "PCM charges (in)", error)
      call device_alloc(this%d_q_out, this%n_pcm_points, "PCM charges (out)", error)
      if (error%has_error()) return

      ! Zero, so the first solve starts from an uncharged surface rather than from
      ! whatever the allocation happened to contain. Uploaded from the host rather
      ! than memset on the device: there is no zeroing helper here, and this runs
      ! once per geometry.
      block
         real(dp), allocatable :: zeros(:)
         allocate (zeros(this%n_pcm_points))
         zeros = 0.0_dp
         call copy_to_device(this%d_q_in, zeros, "PCM charges (in)", error)
         deallocate (zeros)
      end block
      if (error%has_error()) return

      this%pcm_tolerance = pcm%tolerance
      this%pcm_max_iter = pcm%max_iter
      this%pcm_dielectric = pcm%dielectric
      this%pcm_zeta = pcm%zeta
      this%pcm_radii_scale = pcm%radii_scale
      this%pcm_angular_points = pcm%angular_points
      this%has_pcm = .true.
   end subroutine build_pcm_plan

   subroutine system_pcm_device(this, d_density, pcm_energy, error)
      !! Solve for the surface charges and build the continuum's Fock term
      !!
      !! One call per SCF iteration, with the same shape as `xc_device`: density
      !! in, a potential matrix left in a device buffer, and an energy scalar
      !! back on the host. The energy is the dielectric (polarization) energy and
      !! already carries its factor of one half, so the caller adds it to the
      !! total directly -- exactly as it does E_xc, and for the same reason. The
      !! Fock matrix meanwhile takes the *full* potential, because the surface
      !! charges are determined variationally and the half does not appear in the
      !! derivative.
      !!
      !! The previous iteration's charges are handed in as the starting point.
      !! That is not just an optimisation: near convergence the solve then starts
      !! from an almost-correct surface and takes very few iterations, which is
      !! what keeps a continuum from doubling the cost of an SCF.
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: d_density   !! Total density, (n_ao, n_ao) on device
      real(dp), intent(out) :: pcm_energy    !! Dielectric energy, Hartree
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params, results
      integer(c_int) :: status
      integer(c_int64_t) :: max_iter, iterations_taken
      real(dp) :: residual

      pcm_energy = 0.0_dp
      if (error%has_error() .or. .not. this%has_pcm) return

      params = c_null_ptr
      results = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(PCM potential)", error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuest_param_set_f64(CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS, params, &
                                                     CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS_CONVERGENCE_THRESHOLD, &
                                                     this%pcm_tolerance), &
                                 "PCM potential convergence threshold", error)
      end if
      if (.not. error%has_error()) then
         max_iter = int(this%pcm_max_iter, c_int64_t)
         call cuest_status_check(cuest_param_set_i64(CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS, params, &
                                                     CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS_MAX_ITERATIONS, &
                                                     max_iter), &
                                 "PCM potential iteration limit", error)
      end if
      if (.not. error%has_error()) then
         call cuest_status_check(cuestResultsCreate(CUEST_PCM_RESULTS, results), &
                                 "cuestResultsCreate(PCM)", error)
      end if

      if (.not. error%has_error()) then
         call cuest_status_check(cuestPCMPotentialComputeWorkspaceQuery(this%handle, this%pcm_plan, &
                                                                        params, temporary_desc, &
                                                                        d_density, this%d_q_in, &
                                                                        this%d_q_out, results, &
                                                                        this%d_pcm), &
                                 "cuestPCMPotentialComputeWorkspaceQuery", error)
      end if
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestPCMPotentialCompute(this%handle, this%pcm_plan, params, &
                                                          temporary_ws, d_density, this%d_q_in, &
                                                          this%d_q_out, results, this%d_pcm), &
                                 "cuestPCMPotentialCompute", error)
      end if

      if (.not. error%has_error()) then
         call cuest_status_check(cuest_results_query_f64(CUEST_PCM_RESULTS, results, &
                                                         CUEST_PCMRESULT_PCM_DIELECTRIC_ENERGY, &
                                                         pcm_energy), &
                                 "query PCM dielectric energy", error)
      end if
      if (.not. error%has_error()) then
         call cuest_status_check(cuest_results_query_f64(CUEST_PCM_RESULTS, results, &
                                                         CUEST_PCMRESULT_CONVERGED_RESIDUAL, &
                                                         residual), &
                                 "query PCM residual", error)
         if (.not. error%has_error()) this%pcm_residual = residual
      end if
      if (.not. error%has_error()) then
         call cuest_status_check(cuest_results_query_i64(CUEST_PCM_RESULTS, results, &
                                                         CUEST_PCMRESULT_NUM_ITERATIONS_TAKEN, &
                                                         iterations_taken), &
                                 "query PCM iteration count", error)
         if (.not. error%has_error()) this%pcm_iterations = int(iterations_taken)
      end if
      ! Whether it solved, derived from the residual rather than read from
      ! `CUEST_PCMRESULT_CONVERGED`.
      !
      ! That attribute is not eight bytes wide -- querying it as an int64 returns
      ! CUEST_STATUS_INVALID_SIZE, while the residual and the iteration count
      ! read back fine as f64 and i64 -- and nothing in the bindings gives its
      ! type, so reading it means guessing a width. The residual against the
      ! threshold we asked for is the same statement without the guess: it is
      ! what "converged" means for an iterative solve, and it is a number this
      ! call is already known to return.
      if (.not. error%has_error()) then
         this%pcm_solved = (this%pcm_residual <= this%pcm_tolerance)
      end if

      ! This iteration's charges become the next one's starting point.
      if (.not. error%has_error()) then
         call cublas_status_check(cublasDcopy(this%cublas, int(this%n_pcm_points, c_int), &
                                              this%d_q_out, 1, this%d_q_in, 1), &
                                  "cublasDcopy(PCM charges out -> in)", error)
      end if

      call workspace_free(temporary_ws)
      if (c_associated(results)) then
         status = cuestResultsDestroy(CUEST_PCM_RESULTS, results)
         call cuest_status_check(status, "cuestResultsDestroy(PCM)", error)
      end if
      if (c_associated(params)) then
         status = cuestParametersDestroy(CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS, params)
         call cuest_status_check(status, "cuestParametersDestroy(PCM potential)", error)
      end if

      if (error%has_error()) pcm_energy = 0.0_dp
   end subroutine system_pcm_device

   subroutine system_gradient_pcm(this, density, gradient, error)
      !! The continuum's contribution to the nuclear gradient
      !!
      !! Everything in the dielectric energy depends on where the nuclei are:
      !! the surface points sit on atom-centred spheres, the switching weights
      !! are functions of interatomic distances, and the surface charges
      !! interact with the nuclei directly. cuEST differentiates all of that
      !! behind one call, in the same way it builds the cavity itself.
      !!
      !! The surface charges do *not* need a response term. They are determined
      !! variationally, so the energy is stationary with respect to them and
      !! their derivative contributes nothing at first order -- the same reason
      !! `pcm_device` hands the Fock matrix the full potential rather than half
      !! of it. What would otherwise be the expensive part of a continuum
      !! gradient is therefore simply absent.
      !!
      !! `d_q_in` holds the converged charges on entry: the last SCF iteration
      !! copies `q_out` into it, so a gradient taken after a converged SCF
      !! starts from the right surface rather than from zero.
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: density(:, :)    !! Total density, (n_ao, n_ao)
      real(dp), intent(out) :: gradient(:, :)  !! (3, n_atoms), Hartree/Bohr
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params, results
      integer(c_int) :: status
      integer(c_int64_t) :: max_iter
      real(dp) :: residual

      gradient = 0.0_dp
      if (error%has_error() .or. .not. this%has_pcm) return

      call stage_matrix(this, this%d_matrix, density, "density (PCM gradient)", error)
      if (error%has_error()) return

      params = c_null_ptr
      results = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(PCM derivative)", error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuest_param_set_f64(CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS, params, &
                                                     CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS_CONVERGENCE_THRESHOLD, &
                                                     this%pcm_tolerance), &
                                 "PCM derivative convergence threshold", error)
      end if
      if (.not. error%has_error()) then
         max_iter = int(this%pcm_max_iter, c_int64_t)
         call cuest_status_check(cuest_param_set_i64(CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS, params, &
                                                     CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS_MAX_ITERATIONS, &
                                                     max_iter), &
                                 "PCM derivative iteration limit", error)
      end if
      if (.not. error%has_error()) then
         call cuest_status_check(cuestResultsCreate(CUEST_PCM_RESULTS, results), &
                                 "cuestResultsCreate(PCM derivative)", error)
      end if

      if (.not. error%has_error()) then
         call cuest_status_check(cuestPCMDerivativeComputeWorkspaceQuery(this%handle, this%pcm_plan, &
                                                                         params, temporary_desc, &
                                                                         this%d_matrix, this%d_q_in, &
                                                                         this%d_q_out, results, &
                                                                         this%d_gradient), &
                                 "cuestPCMDerivativeComputeWorkspaceQuery", error)
      end if
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestPCMDerivativeCompute(this%handle, this%pcm_plan, params, &
                                                           temporary_ws, this%d_matrix, this%d_q_in, &
                                                           this%d_q_out, results, this%d_gradient), &
                                 "cuestPCMDerivativeCompute", error)
      end if

      ! The derivative runs its own iterative solve, so it converges or it does
      ! not, exactly as the energy's does. Read the residual rather than
      ! `CUEST_PCMRESULT_CONVERGED` for the reason given in `pcm_device`: that
      ! attribute's width is not recoverable from the bindings.
      if (.not. error%has_error()) then
         call cuest_status_check(cuest_results_query_f64(CUEST_PCM_RESULTS, results, &
                                                         CUEST_PCMRESULT_CONVERGED_RESIDUAL, &
                                                         residual), &
                                 "query PCM derivative residual", error)
         if (.not. error%has_error()) then
            this%pcm_residual = residual
            this%pcm_solved = (residual <= this%pcm_tolerance)
         end if
      end if

      call workspace_free(temporary_ws)
      if (c_associated(results)) then
         status = cuestResultsDestroy(CUEST_PCM_RESULTS, results)
         call cuest_status_check(status, "cuestResultsDestroy(PCM derivative)", error)
      end if
      if (c_associated(params)) then
         status = cuestParametersDestroy(CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS, params)
         call cuest_status_check(status, "cuestParametersDestroy(PCM derivative)", error)
      end if

      call fetch_gradient(this, gradient, "dE_diel/dR", error)
   end subroutine system_gradient_pcm

   subroutine system_xc_device(this, d_c_occ, d_out, xc_energy, error)
      !! RKS exchange-correlation potential, device in and device out
      !!
      !! `xc_energy` is a host scalar even here -- cuEST writes Exc through a
      !! host pointer, which is why this one call does synchronize internally
      !! whatever the surrounding code does. On a pure Hartree-Fock system
      !! there is no XC plan, so `d_out` is left untouched and the energy is
      !! zero: skip the term rather than adding a stale buffer.
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: d_c_occ  !! Occupied MOs, (n_ao, n_occ) on device
      type(c_ptr), intent(in) :: d_out    !! Vxc, (n_ao, n_ao) on device
      real(dp), intent(out) :: xc_energy  !! Exc, Hartree
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      real(dp), target :: energy

      xc_energy = 0.0_dp
      if (error%has_error() .or. .not. this%has_xc) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_XCPOTENTIALRKSCOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(RKS XC potential)", error)
      if (error%has_error()) return

      variable_buffer%hostBufferSizeInBytes = 0_c_size_t
      variable_buffer%deviceBufferSizeInBytes = DF_EXCHANGE_BUFFER_BYTES

      ! Unlike every other output here, the energy is a HOST scalar.
      energy = 0.0_dp
      call cuest_status_check(cuestXCPotentialRKSComputeWorkspaceQuery(this%handle, this%xc_plan, &
                                                                       params, variable_buffer, &
                                                                       temporary_desc, this%n_occ, &
                                                                       d_c_occ, c_loc(energy), &
                                                                       d_out), &
                              "cuestXCPotentialRKSComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestXCPotentialRKSCompute(this%handle, this%xc_plan, params, &
                                                            variable_buffer, temporary_ws, &
                                                            this%n_occ, d_c_occ, &
                                                            c_loc(energy), d_out), &
                                 "cuestXCPotentialRKSCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_XCPOTENTIALRKSCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(RKS XC potential)", error)

      if (.not. error%has_error()) xc_energy = energy
   end subroutine system_xc_device

   subroutine system_compute_xc(this, c_occ, xc_energy, xc_potential, error)
      !! Exchange-correlation energy and potential from host MO coefficients
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: c_occ(:, :)           !! (n_ao, n_occ)
      real(dp), intent(out) :: xc_energy            !! Exc, Hartree
      real(dp), intent(out) :: xc_potential(:, :)   !! Vxc, (n_ao, n_ao)
      type(error_t), intent(inout) :: error

      xc_energy = 0.0_dp
      xc_potential = 0.0_dp
      if (error%has_error() .or. .not. this%has_xc) return

      call this%stage_occupied(c_occ, error)
      if (error%has_error()) return

      call this%xc_device(this%d_c_occ, this%d_result, xc_energy, error)
      call fetch_matrix(this, xc_potential, "Vxc", error)
      if (error%has_error()) xc_energy = 0.0_dp
   end subroutine system_compute_xc

   subroutine system_compute_dipole(this, density, dipole, error)
      !! Electric dipole moment, in atomic units (electron-Bohr)
      !!
      !!   mu = sum_A Z_A (R_A - O)  -  sum_uv D_uv <u| r - O |v>
      !!
      !! The electronic term is negative because the electron charge is. The
      !! origin O is the centre of nuclear charge, which makes the result
      !! origin-independent for a neutral system and at least well defined for
      !! a charged one. Note that dmu/dR -- what IR intensities actually need
      !! -- is origin-independent only for a neutral system.
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: density(:, :)   !! Total density, (n_ao, n_ao)
      real(dp), intent(out) :: dipole(3)      !! mu_x, mu_y, mu_z
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      integer(c_int32_t), target :: multipole_order(3)
      real(dp), target :: origin(3)
      real(dp), allocatable :: component(:, :)
      real(dp) :: total_charge
      integer :: icomp, iatom

      dipole = 0.0_dp
      if (error%has_error()) return

      ! Centre of nuclear charge. charges_host holds -Z, hence the negation.
      total_charge = -sum(this%charges_host)
      origin = 0.0_dp
      if (total_charge > 0.0_dp) then
         do iatom = 1, int(this%n_atoms)
            origin = origin - this%charges_host(iatom) &
                     *this%xyz_host(3*(iatom - 1) + 1:3*iatom)
         end do
         origin = origin/total_charge
      end if

      ! Nuclear term.
      do iatom = 1, int(this%n_atoms)
         dipole = dipole - this%charges_host(iatom) &
                  *(this%xyz_host(3*(iatom - 1) + 1:3*iatom) - origin)
      end do

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_MULTIPOLECOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(multipole)", error)
      if (error%has_error()) return

      allocate (component(int(this%n_ao), int(this%n_ao)))
      do icomp = 1, 3
         ! Exponents of x^l y^m z^n: (1,0,0), (0,1,0), (0,0,1).
         multipole_order = 0_c_int32_t
         multipole_order(icomp) = 1_c_int32_t

         call cuest_status_check(cuestMultipoleComputeWorkspaceQuery(this%handle, this%oe_plan, &
                                                                     params, temporary_desc, &
                                                                     multipole_order, c_loc(origin), &
                                                                     this%d_result), &
                                 "cuestMultipoleComputeWorkspaceQuery", error)
         if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
         if (.not. error%has_error()) then
            call cuest_status_check(cuestMultipoleCompute(this%handle, this%oe_plan, params, &
                                                          temporary_ws, multipole_order, &
                                                          c_loc(origin), this%d_result), &
                                    "cuestMultipoleCompute", error)
         end if
         call workspace_free(temporary_ws)
         if (error%has_error()) exit

         call fetch_matrix(this, component, "multipole", error)
         if (error%has_error()) exit
         dipole(icomp) = dipole(icomp) - sum(density*component)
      end do
      deallocate (component)

      status = cuestParametersDestroy(CUEST_MULTIPOLECOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(multipole)", error)
      if (error%has_error()) dipole = 0.0_dp
   end subroutine system_compute_dipole

   subroutine system_compute_overlap(this, overlap, error)
      !! Overlap matrix S
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(out) :: overlap(:, :)  !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status

      if (error%has_error()) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_OVERLAPCOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(overlap)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestOverlapComputeWorkspaceQuery(this%handle, this%oe_plan, &
                                                                params, temporary_desc, &
                                                                this%d_result), &
                              "cuestOverlapComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestOverlapCompute(this%handle, this%oe_plan, params, &
                                                     temporary_ws, this%d_result), &
                                 "cuestOverlapCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_OVERLAPCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(overlap)", error)

      call fetch_matrix(this, overlap, "S", error)
   end subroutine system_compute_overlap

   subroutine system_compute_kinetic(this, kinetic, error)
      !! Kinetic energy matrix T
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(out) :: kinetic(:, :)  !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status

      if (error%has_error()) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_KINETICCOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(kinetic)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestKineticComputeWorkspaceQuery(this%handle, this%oe_plan, &
                                                                params, temporary_desc, &
                                                                this%d_result), &
                              "cuestKineticComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestKineticCompute(this%handle, this%oe_plan, params, &
                                                     temporary_ws, this%d_result), &
                                 "cuestKineticCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_KINETICCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(kinetic)", error)

      call fetch_matrix(this, kinetic, "T", error)
   end subroutine system_compute_kinetic

   subroutine system_compute_potential(this, potential, error)
      !! Nuclear attraction matrix V
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(out) :: potential(:, :)  !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status

      if (error%has_error()) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_POTENTIALCOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(potential)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestPotentialComputeWorkspaceQuery(this%handle, this%oe_plan, &
                                                                  params, temporary_desc, &
                                                                  this%n_atoms, this%xyz_device, &
                                                                  this%charges_device, this%d_result), &
                              "cuestPotentialComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestPotentialCompute(this%handle, this%oe_plan, params, &
                                                       temporary_ws, this%n_atoms, &
                                                       this%xyz_device, this%charges_device, &
                                                       this%d_result), &
                                 "cuestPotentialCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_POTENTIALCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(potential)", error)

      call fetch_matrix(this, potential, "V", error)
   end subroutine system_compute_potential

   subroutine system_coulomb_device(this, d_density, d_out, error)
      !! Density-fitted Coulomb matrix J, device in and device out
      !!
      !! Neither uploads nor fetches: the density is already on the device and
      !! J is left there for the caller to combine. The work is queued, not
      !! completed -- synchronize before reading `d_out` from the host.
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: d_density  !! Density, (n_ao, n_ao) on device
      type(c_ptr), intent(in) :: d_out      !! J, (n_ao, n_ao) on device
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status

      if (error%has_error()) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_DFCOULOMBCOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(DF Coulomb)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestDFCoulombComputeWorkspaceQuery(this%handle, this%df_plan, &
                                                                  params, temporary_desc, &
                                                                  d_density, d_out), &
                              "cuestDFCoulombComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestDFCoulombCompute(this%handle, this%df_plan, params, &
                                                       temporary_ws, d_density, d_out), &
                                 "cuestDFCoulombCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_DFCOULOMBCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(DF Coulomb)", error)
   end subroutine system_coulomb_device

   subroutine system_compute_coulomb(this, density, coulomb, error)
      !! Density-fitted Coulomb matrix J from a host density matrix
      !!
      !! The host-to-host form, kept for the gradient path and the atomic
      !! guess. The SCF drives `coulomb_device` instead.
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: density(:, :)   !! (n_ao, n_ao), symmetric
      real(dp), intent(out) :: coulomb(:, :)  !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return

      call stage_matrix(this, this%d_matrix, density, "density matrix", error)
      if (error%has_error()) return

      call this%coulomb_device(this%d_matrix, this%d_result, error)
      call fetch_matrix(this, coulomb, "J", error)
   end subroutine system_compute_coulomb

   subroutine system_exchange_device(this, d_c_occ, d_out, error, n_occupied)
      !! Density-fitted exchange matrix K, device in and device out
      !!
      !! `n_occupied` overrides the stored count, which is what the two spin
      !! channels of an unrestricted calculation need. An empty spin channel is
      !! the caller's problem here: with nothing to contract, `d_out` is left
      !! untouched rather than being zeroed, so skip the term instead of
      !! adding it.
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: d_c_occ  !! Occupied MOs, (n_ao, n_occ) on device
      type(c_ptr), intent(in) :: d_out    !! K, (n_ao, n_ao) on device
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_occupied

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      integer(c_int64_t) :: occupancy

      if (error%has_error()) return

      occupancy = this%n_occ
      if (present(n_occupied)) occupancy = int(n_occupied, c_int64_t)
      if (occupancy <= 0) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS, &
                                                    params), &
                              "cuestParametersCreate(DF exchange)", error)
      if (error%has_error()) return

      variable_buffer%hostBufferSizeInBytes = 0_c_size_t
      variable_buffer%deviceBufferSizeInBytes = DF_EXCHANGE_BUFFER_BYTES

      call cuest_status_check(cuestDFSymmetricExchangeComputeWorkspaceQuery(this%handle, &
                                                                            this%df_plan, params, &
                                                                            variable_buffer, &
                                                                            temporary_desc, &
                                                                            occupancy, d_c_occ, &
                                                                            d_out), &
                              "cuestDFSymmetricExchangeComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestDFSymmetricExchangeCompute(this%handle, this%df_plan, &
                                                                 params, variable_buffer, &
                                                                 temporary_ws, occupancy, &
                                                                 d_c_occ, d_out), &
                                 "cuestDFSymmetricExchangeCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(DF exchange)", error)
   end subroutine system_exchange_device

   subroutine system_compute_exchange(this, c_occ, exchange, error, n_occupied)
      !! Density-fitted exchange matrix K from host occupied MO coefficients
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: c_occ(:, :)      !! (n_ao, n_occ)
      real(dp), intent(out) :: exchange(:, :)  !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_occupied

      integer(c_int64_t) :: occupancy

      if (error%has_error()) return

      occupancy = this%n_occ
      if (present(n_occupied)) occupancy = int(n_occupied, c_int64_t)
      if (occupancy <= 0) then
         exchange = 0.0_dp
         return
      end if

      call this%stage_occupied(c_occ, error)
      if (error%has_error()) return

      call this%exchange_device(this%d_c_occ, this%d_result, error, n_occupied=n_occupied)
      call fetch_matrix(this, exchange, "K", error)
   end subroutine system_compute_exchange

   subroutine fetch_matrix(this, matrix, label, error)
      !! Synchronize, then copy the shared result buffer into a host matrix
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(out) :: matrix(:, :)
      character(len=*), intent(in) :: label
      type(error_t), intent(inout) :: error

      call this%fetch(this%d_result, matrix, label, error)
   end subroutine fetch_matrix

   ! ==========================================================================
   !  Device-resident SCF support
   !
   !  Staging puts a host matrix into a named device buffer; the compute_*
   !  _device routines above leave their results there; assembly and traces
   !  then run on the device through cuBLAS. Only the assembled Fock and two
   !  scalars cross back per SCF iteration.
   ! ==========================================================================

   subroutine system_stage_to(this, device_ptr, matrix, label, error)
      !! Copy a host matrix into a caller-chosen device buffer
      !!
      !! The size comes from the array, so this serves the n_ao x n_ao matrices
      !! and the n_mo x n_mo DIIS error alike; the caller owns the sizing.
      !!
      !! Matrices crossing this boundary are symmetric, so the row-major /
      !! column-major difference is a no-op -- see the note at the top of this
      !! module before staging anything that is not.
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: device_ptr
      real(dp), intent(in) :: matrix(:, :)
      character(len=*), intent(in) :: label
      type(error_t), intent(inout) :: error

      call stage_matrix(this, device_ptr, matrix, label, error)
   end subroutine system_stage_to

   subroutine system_stage_density(this, density, error)
      !! Copy the host density into the shared density buffer
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: density(:, :)  !! (n_ao, n_ao), symmetric
      type(error_t), intent(inout) :: error

      call stage_matrix(this, this%d_matrix, density, "density matrix", error)
   end subroutine system_stage_density

   subroutine system_stage_occupied(this, c_occ, error)
      !! Copy the host occupied MO coefficients into the shared alpha buffer
      !!
      !! K and Vxc take the same coefficients, so the SCF stages once per
      !! iteration and both device routines read the result.
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: c_occ(:, :)  !! (n_ao, n_occ)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat(:)

      if (error%has_error()) return

      allocate (flat(size(c_occ)))
      flat = reshape(c_occ, [size(c_occ)])
      call copy_to_device(this%d_c_occ, flat, "occupied MO coefficients", error)
      deallocate (flat)
   end subroutine system_stage_occupied

   subroutine system_assemble_fock(this, d_out, d_exchange, d_xc_term, &
                                   add_exchange, add_xc, error)
      !! F := H + J - K + Vxc, entirely on the device
      !!
      !! H and J are spin-independent -- the Coulomb term is built from the
      !! total density -- so they are read from the system's own buffers, while
      !! the output and the two spin-dependent terms are named by the caller.
      !! An unrestricted iteration therefore calls this twice, once per channel.
      !!
      !! `add_exchange` and `add_xc` say whether those buffers hold anything.
      !! A pure functional never runs the exchange call, Hartree-Fock never runs
      !! the XC call, and an empty beta channel runs neither; in every case the
      !! buffer holds whatever the last fragment left there, so adding it would
      !! be the classic silent stale-memory failure. The term is skipped rather
      !! than zeroed.
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: d_out        !! F for this spin, (n_ao, n_ao) on device
      type(c_ptr), intent(in) :: d_exchange   !! K for this spin
      type(c_ptr), intent(in) :: d_xc_term    !! Vxc for this spin
      logical, intent(in) :: add_exchange     !! Whether d_exchange was written this iteration
      logical, intent(in) :: add_xc           !! Whether d_xc_term was written this iteration
      type(error_t), intent(inout) :: error

      integer(c_int) :: n

      if (error%has_error()) return
      n = int(this%n_ao*this%n_ao, c_int)

      call cublas_status_check(cublasDcopy(this%cublas, n, this%d_core, 1, d_out, 1), &
                               "cublasDcopy(H -> F)", error)
      call cublas_status_check(cublasDaxpy(this%cublas, n, 1.0_dp, this%d_j, 1, d_out, 1), &
                               "cublasDaxpy(+J)", error)
      if (add_exchange) then
         call cublas_status_check(cublasDaxpy(this%cublas, n, -1.0_dp, d_exchange, 1, &
                                              d_out, 1), &
                                  "cublasDaxpy(-K)", error)
      end if
      if (add_xc) then
         call cublas_status_check(cublasDaxpy(this%cublas, n, 1.0_dp, d_xc_term, 1, &
                                              d_out, 1), &
                                  "cublasDaxpy(+Vxc)", error)
      end if
      ! The continuum, read from the system's own buffer like H and J rather than
      ! passed in: the surface charges come from the *total* density, so both spin
      ! channels see the same potential and there is nothing per-spin to name.
      if (this%has_pcm) then
         call cublas_status_check(cublasDaxpy(this%cublas, n, 1.0_dp, this%d_pcm, 1, &
                                              d_out, 1), &
                                  "cublasDaxpy(+Vpcm)", error)
      end if
   end subroutine system_assemble_fock

   subroutine system_add_into(this, d_a, d_b, d_out, error)
      !! d_out := d_a + d_b for n_ao x n_ao device matrices
      !!
      !! The total density an unrestricted Coulomb build needs, formed where
      !! its two halves already are.
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: d_a, d_b
      type(c_ptr), intent(in) :: d_out
      type(error_t), intent(inout) :: error

      integer(c_int) :: n

      if (error%has_error()) return
      n = int(this%n_ao*this%n_ao, c_int)

      call cublas_status_check(cublasDcopy(this%cublas, n, d_a, 1, d_out, 1), &
                               "cublasDcopy(D^a -> D^t)", error)
      call cublas_status_check(cublasDaxpy(this%cublas, n, 1.0_dp, d_b, 1, d_out, 1), &
                               "cublasDaxpy(+D^b)", error)
   end subroutine system_add_into

   subroutine system_matrix_dot(this, d_a, d_b, dot, label, error)
      !! sum_uv A_uv B_uv for two device-resident n_ao x n_ao matrices
      !!
      !! The Frobenius inner product is the dot product of the flattenings, so
      !! the row-major/column-major difference does not enter -- and for the
      !! symmetric matrices the SCF feeds it, neither does which one is which.
      !!
      !! The result comes back through a HOST pointer, so this call blocks
      !! until the stream drains. Two per iteration is fine; a loop over one
      !! per DIIS vector would be a synchronise per vector, which is what
      !! cublasDgemv exists to avoid.
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: d_a, d_b     !! (n_ao, n_ao) on device
      real(dp), intent(out) :: dot
      character(len=*), intent(in) :: label   !! What is being contracted, for errors
      type(error_t), intent(inout) :: error

      dot = 0.0_dp
      if (error%has_error()) return

      call cublas_status_check(cublasDdot(this%cublas, int(this%n_ao*this%n_ao, c_int), &
                                          d_a, 1, d_b, 1, dot), &
                               "cublasDdot("//trim(label)//")", error)
      if (error%has_error()) dot = 0.0_dp
   end subroutine system_matrix_dot

   ! ==========================================================================
   !  Nuclear derivatives
   !
   !  cuEST never exposes raw derivative integrals: the contraction with a
   !  density (or with occupied orbitals) happens inside the call and the
   !  result is a natom x 3 gradient. Every one of these OVERWRITES its output
   !  buffer rather than accumulating, so the caller sums the contributions.
   ! ==========================================================================

   subroutine fetch_gradient(this, gradient, label, error)
      !! Synchronize and pull an natom x 3 gradient back to the host
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(out) :: gradient(:, :)   !! (3, n_atoms)
      character(len=*), intent(in) :: label
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat(:)

      if (error%has_error()) then
         gradient = 0.0_dp
         return
      end if

      call device_sync(label, error)
      if (error%has_error()) return

      allocate (flat(3*this%n_atoms))
      call copy_to_host(flat, this%d_gradient, label, error)
      if (.not. error%has_error()) then
         ! cuEST lays the gradient out atom-major (natom x 3 row-major), which
         ! is byte-identical to a Fortran (3, natom) array.
         gradient = reshape(flat, [3, int(this%n_atoms)])
      else
         gradient = 0.0_dp
      end if
      deallocate (flat)
   end subroutine fetch_gradient

   subroutine system_gradient_overlap(this, weighted_density, gradient, error)
      !! Overlap (Pulay) derivative contracted with a matrix
      !!
      !! Pass the energy-weighted density; the caller subtracts the result,
      !! since the Pulay term enters the gradient with a minus sign.
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: weighted_density(:, :)  !! (n_ao, n_ao)
      real(dp), intent(out) :: gradient(:, :)         !! (3, n_atoms)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status

      if (error%has_error()) then
         gradient = 0.0_dp
         return
      end if

      call stage_matrix(this, this%d_matrix, weighted_density, "energy-weighted density", error)
      if (error%has_error()) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_OVERLAPDERIVATIVECOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(overlap derivative)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestOverlapDerivativeComputeWorkspaceQuery(this%handle, this%oe_plan, &
                                                                          params, temporary_desc, &
                                                                          this%d_matrix, this%d_gradient), &
                              "cuestOverlapDerivativeComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestOverlapDerivativeCompute(this%handle, this%oe_plan, params, &
                                                               temporary_ws, this%d_matrix, &
                                                               this%d_gradient), &
                                 "cuestOverlapDerivativeCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_OVERLAPDERIVATIVECOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(overlap derivative)", error)

      call fetch_gradient(this, gradient, "dS/dR", error)
   end subroutine system_gradient_overlap

   subroutine system_gradient_kinetic(this, density, gradient, error)
      !! Kinetic energy derivative contracted with the total density
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: density(:, :)    !! (n_ao, n_ao)
      real(dp), intent(out) :: gradient(:, :)  !! (3, n_atoms)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status

      if (error%has_error()) then
         gradient = 0.0_dp
         return
      end if

      call stage_matrix(this, this%d_matrix, density, "density (kinetic gradient)", error)
      if (error%has_error()) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_KINETICDERIVATIVECOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(kinetic derivative)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestKineticDerivativeComputeWorkspaceQuery(this%handle, this%oe_plan, &
                                                                          params, temporary_desc, &
                                                                          this%d_matrix, this%d_gradient), &
                              "cuestKineticDerivativeComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestKineticDerivativeCompute(this%handle, this%oe_plan, params, &
                                                               temporary_ws, this%d_matrix, &
                                                               this%d_gradient), &
                                 "cuestKineticDerivativeCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_KINETICDERIVATIVECOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(kinetic derivative)", error)

      call fetch_gradient(this, gradient, "dT/dR", error)
   end subroutine system_gradient_kinetic

   subroutine system_gradient_potential(this, density, gradient, error)
      !! Nuclear attraction derivative, both basis-centre and Hellmann-Feynman
      !!
      !! The two pieces are returned separately -- the derivative with respect
      !! to the AO centres and the derivative with respect to the nuclear
      !! positions themselves -- and both belong in the total gradient.
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: density(:, :)    !! (n_ao, n_ao)
      real(dp), intent(out) :: gradient(:, :)  !! (3, n_atoms)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      real(dp), allocatable :: charge_gradient(:, :)

      if (error%has_error()) then
         gradient = 0.0_dp
         return
      end if

      call stage_matrix(this, this%d_matrix, density, "density (potential gradient)", error)
      if (error%has_error()) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_POTENTIALDERIVATIVECOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(potential derivative)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestPotentialDerivativeComputeWorkspaceQuery(this%handle, this%oe_plan, &
                                                                            params, temporary_desc, &
                                                                            this%n_atoms, this%xyz_device, &
                                                                            this%charges_device, this%d_matrix, &
                                                                            this%d_gradient, &
                                                                            this%d_charge_gradient), &
                              "cuestPotentialDerivativeComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestPotentialDerivativeCompute(this%handle, this%oe_plan, params, &
                                                                 temporary_ws, this%n_atoms, &
                                                                 this%xyz_device, this%charges_device, &
                                                                 this%d_matrix, this%d_gradient, &
                                                                 this%d_charge_gradient), &
                                 "cuestPotentialDerivativeCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_POTENTIALDERIVATIVECOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(potential derivative)", error)

      call fetch_gradient(this, gradient, "dV/dR (basis)", error)
      if (error%has_error()) return

      ! Pull the Hellmann-Feynman half from its own buffer and add it in.
      allocate (charge_gradient(3, this%n_atoms))
      call fetch_named_gradient(this, this%d_charge_gradient, charge_gradient, &
                                "dV/dR (charges)", error)
      if (.not. error%has_error()) gradient = gradient + charge_gradient
      deallocate (charge_gradient)
   end subroutine system_gradient_potential

   subroutine system_gradient_two_electron(this, half_density, c_occ, gradient, error)
      !! Density-fitted Coulomb + exchange derivative
      !!
      !! cuEST defines E_JK = s_D E_J[D] + s_C E_K[C] and differentiates that.
      !! For closed-shell RKS the intended convention -- straight from the
      !! header -- is D = sum_i C_i C_i (the HALF density, not the total one
      !! the SCF carries), s_D = 2 and s_C = -1.
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: half_density(:, :)  !! (n_ao, n_ao), = D_total/2
      real(dp), intent(in) :: c_occ(:, :)         !! (n_ao, n_occ)
      real(dp), intent(out) :: gradient(:, :)     !! (3, n_atoms)
      type(error_t), intent(inout) :: error

      real(dp), parameter :: DENSITY_SCALE = 2.0_dp
      real(dp), parameter :: RKS_EXCHANGE_SCALE = -1.0_dp

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params, coefficient_ptr
      integer(c_int) :: status
      integer(c_int64_t) :: occupancies(1)
      integer(c_int64_t) :: n_matrices
      real(dp) :: coefficient_scale
      real(dp), allocatable :: flat(:)

      if (error%has_error()) then
         gradient = 0.0_dp
         return
      end if

      call stage_matrix(this, this%d_matrix, half_density, "half density (JK gradient)", error)
      if (error%has_error()) return

      ! A pure functional has no exchange term; the header allows passing zero
      ! coefficient matrices in that case.
      if (this%needs_exchange) then
         allocate (flat(size(c_occ)))
         flat = reshape(c_occ, [size(c_occ)])
         call copy_to_device(this%d_c_occ, flat, "occupied MOs (JK gradient)", error)
         deallocate (flat)
         if (error%has_error()) return
         n_matrices = 1_c_int64_t
         occupancies(1) = this%n_occ
         coefficient_ptr = this%d_c_occ
         coefficient_scale = RKS_EXCHANGE_SCALE
      else
         n_matrices = 0_c_int64_t
         occupancies(1) = 0_c_int64_t
         coefficient_ptr = c_null_ptr
         coefficient_scale = 0.0_dp
      end if

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_DFSYMMETRICDERIVATIVECOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(DF derivative)", error)
      if (error%has_error()) return

      variable_buffer%hostBufferSizeInBytes = 0_c_size_t
      variable_buffer%deviceBufferSizeInBytes = DF_EXCHANGE_BUFFER_BYTES

      call cuest_status_check(cuestDFSymmetricDerivativeComputeWorkspaceQuery(this%handle, this%df_plan, &
                                                                              params, variable_buffer, &
                                                                              temporary_desc, DENSITY_SCALE, &
                                                                              this%d_matrix, coefficient_scale, &
                                                                              n_matrices, occupancies, &
                                                                              coefficient_ptr, this%d_gradient), &
                              "cuestDFSymmetricDerivativeComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestDFSymmetricDerivativeCompute(this%handle, this%df_plan, params, &
                                                                   variable_buffer, temporary_ws, &
                                                                   DENSITY_SCALE, this%d_matrix, &
                                                                   coefficient_scale, n_matrices, &
                                                                   occupancies, coefficient_ptr, &
                                                                   this%d_gradient), &
                                 "cuestDFSymmetricDerivativeCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_DFSYMMETRICDERIVATIVECOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(DF derivative)", error)

      call fetch_gradient(this, gradient, "dE_JK/dR", error)
   end subroutine system_gradient_two_electron

   subroutine system_gradient_xc(this, c_occ, gradient, error)
      !! Exchange-correlation derivative at fixed grid
      !!
      !! Complete for cuEST's built-in functionals: the grid-weight (Becke
      !! partition) terms are included here. The separate
      !! cuestXCGridDerivativeCompute serves the user-supplied-functional
      !! path, where the caller evaluates the XC energy density itself.
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: c_occ(:, :)      !! (n_ao, n_occ)
      real(dp), intent(out) :: gradient(:, :)  !! (3, n_atoms)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      real(dp), allocatable :: flat(:)

      gradient = 0.0_dp
      if (error%has_error() .or. .not. this%has_xc) return

      allocate (flat(size(c_occ)))
      flat = reshape(c_occ, [size(c_occ)])
      call copy_to_device(this%d_c_occ, flat, "occupied MOs (XC gradient)", error)
      deallocate (flat)
      if (error%has_error()) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_XCDERIVATIVERKSCOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(RKS XC derivative)", error)
      if (error%has_error()) return

      variable_buffer%hostBufferSizeInBytes = 0_c_size_t
      variable_buffer%deviceBufferSizeInBytes = DF_EXCHANGE_BUFFER_BYTES

      call cuest_status_check(cuestXCDerivativeRKSComputeWorkspaceQuery(this%handle, this%xc_plan, &
                                                                        params, variable_buffer, &
                                                                        temporary_desc, this%n_occ, &
                                                                        this%d_c_occ, this%d_gradient), &
                              "cuestXCDerivativeRKSComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestXCDerivativeRKSCompute(this%handle, this%xc_plan, params, &
                                                             variable_buffer, temporary_ws, this%n_occ, &
                                                             this%d_c_occ, this%d_gradient), &
                                 "cuestXCDerivativeRKSCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_XCDERIVATIVERKSCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(RKS XC derivative)", error)

      call fetch_gradient(this, gradient, "dExc/dR", error)
   end subroutine system_gradient_xc

   subroutine system_gradient_two_electron_uks(this, total_density, c_occ_alpha, c_occ_beta, &
                                               gradient, error)
      !! Density-fitted Coulomb + exchange derivative, unrestricted
      !!
      !! cuEST wants the two coefficient matrices concatenated, an occupancy
      !! per matrix, and -- from the header's own UKS worked example --
      !! densityScale = 0.5 with the TOTAL density, coefficientScale = -0.5.
      !! That is the same E_JK as the restricted case written for two spin
      !! channels, which is the check that the factors are right.
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: total_density(:, :)  !! D^a + D^b, no factor 2
      real(dp), intent(in) :: c_occ_alpha(:, :), c_occ_beta(:, :)
      real(dp), intent(out) :: gradient(:, :)      !! (3, n_atoms)
      type(error_t), intent(inout) :: error

      real(dp), parameter :: UKS_DENSITY_SCALE = 0.5_dp
      real(dp), parameter :: UKS_EXCHANGE_SCALE = -0.5_dp

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params, coefficient_ptr
      integer(c_int) :: status
      integer(c_int64_t) ::  n_matrices
      integer(c_int64_t) :: occupancies(2)
      real(dp) :: coefficient_scale
      real(dp), allocatable :: flat(:)
      integer :: n_a, n_b

      if (error%has_error()) then
         gradient = 0.0_dp
         return
      end if

      n_a = int(this%n_occ)
      n_b = int(this%n_occ_beta)

      call stage_matrix(this, this%d_matrix, total_density, "total density (UKS JK gradient)", error)
      if (error%has_error()) return

      if (this%needs_exchange) then
         ! Concatenated: alpha columns then beta columns, in one buffer.
         allocate (flat(this%n_ao*(n_a + n_b)))
         flat(1:this%n_ao*n_a) = reshape(c_occ_alpha(:, 1:n_a), [this%n_ao*n_a])
         if (n_b > 0) then
            flat(this%n_ao*n_a + 1:) = reshape(c_occ_beta(:, 1:n_b), [this%n_ao*n_b])
         end if
         call copy_to_device(this%d_c_occ_beta, flat, "concatenated occupied MOs", error)
         deallocate (flat)
         if (error%has_error()) return

         n_matrices = 2_c_int64_t
         occupancies(1) = this%n_occ
         occupancies(2) = this%n_occ_beta
         coefficient_ptr = this%d_c_occ_beta
         coefficient_scale = UKS_EXCHANGE_SCALE
      else
         n_matrices = 0_c_int64_t
         occupancies = 0_c_int64_t
         coefficient_ptr = c_null_ptr
         coefficient_scale = 0.0_dp
      end if

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_DFSYMMETRICDERIVATIVECOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(DF derivative, UKS)", error)
      if (error%has_error()) return

      variable_buffer%hostBufferSizeInBytes = 0_c_size_t
      variable_buffer%deviceBufferSizeInBytes = DF_EXCHANGE_BUFFER_BYTES

      call cuest_status_check(cuestDFSymmetricDerivativeComputeWorkspaceQuery(this%handle, this%df_plan, &
                                                                              params, variable_buffer, &
                                                                              temporary_desc, UKS_DENSITY_SCALE, &
                                                                              this%d_matrix, coefficient_scale, &
                                                                              n_matrices, occupancies, &
                                                                              coefficient_ptr, this%d_gradient), &
                              "cuestDFSymmetricDerivativeComputeWorkspaceQuery(UKS)", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestDFSymmetricDerivativeCompute(this%handle, this%df_plan, params, &
                                                                   variable_buffer, temporary_ws, &
                                                                   UKS_DENSITY_SCALE, this%d_matrix, &
                                                                   coefficient_scale, n_matrices, &
                                                                   occupancies, coefficient_ptr, &
                                                                   this%d_gradient), &
                                 "cuestDFSymmetricDerivativeCompute(UKS)", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_DFSYMMETRICDERIVATIVECOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(DF derivative, UKS)", error)

      call fetch_gradient(this, gradient, "dE_JK/dR (UKS)", error)
   end subroutine system_gradient_two_electron_uks

   subroutine system_gradient_xc_uks(this, c_occ_alpha, c_occ_beta, gradient, error)
      !! Spin-resolved exchange-correlation derivative
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: c_occ_alpha(:, :), c_occ_beta(:, :)
      real(dp), intent(out) :: gradient(:, :)  !! (3, n_atoms)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      integer(c_int64_t) :: beta_occupancy
      real(dp), allocatable :: flat(:)

      gradient = 0.0_dp
      if (error%has_error() .or. .not. this%has_xc) return

      allocate (flat(size(c_occ_alpha)))
      flat = reshape(c_occ_alpha, [size(c_occ_alpha)])
      call copy_to_device(this%d_c_occ, flat, "alpha occupied MOs (XC gradient)", error)
      deallocate (flat)

      ! See compute_xc_uks: an empty beta channel goes in as one zero column.
      beta_occupancy = max(this%n_occ_beta, 1_c_int64_t)
      allocate (flat(this%n_ao*beta_occupancy))
      flat = 0.0_dp
      if (this%n_occ_beta > 0) flat = reshape(c_occ_beta, [size(c_occ_beta)])
      call copy_to_device(this%d_c_occ_beta, flat, "beta occupied MOs (XC gradient)", error)
      deallocate (flat)
      if (error%has_error()) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_XCDERIVATIVEUKSCOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(UKS XC derivative)", error)
      if (error%has_error()) return

      variable_buffer%hostBufferSizeInBytes = 0_c_size_t
      variable_buffer%deviceBufferSizeInBytes = DF_EXCHANGE_BUFFER_BYTES

      call cuest_status_check(cuestXCDerivativeUKSComputeWorkspaceQuery(this%handle, this%xc_plan, &
                                                                        params, variable_buffer, &
                                                                        temporary_desc, this%n_occ, &
                                                                        beta_occupancy, this%d_c_occ, &
                                                                        this%d_c_occ_beta, this%d_gradient), &
                              "cuestXCDerivativeUKSComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestXCDerivativeUKSCompute(this%handle, this%xc_plan, params, &
                                                             variable_buffer, temporary_ws, &
                                                             this%n_occ, beta_occupancy, &
                                                             this%d_c_occ, this%d_c_occ_beta, &
                                                             this%d_gradient), &
                                 "cuestXCDerivativeUKSCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_XCDERIVATIVEUKSCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(UKS XC derivative)", error)

      call fetch_gradient(this, gradient, "dExc/dR (UKS)", error)
   end subroutine system_gradient_xc_uks

   subroutine stage_matrix(this, device_ptr, matrix, label, error)
      !! Copy an n_ao x n_ao host matrix into a named device buffer
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: device_ptr
      real(dp), intent(in) :: matrix(:, :)
      character(len=*), intent(in) :: label
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat(:)

      if (error%has_error()) return

      allocate (flat(size(matrix)))
      flat = reshape(matrix, [size(matrix)])
      call copy_to_device(device_ptr, flat, label, error)
      deallocate (flat)
   end subroutine stage_matrix

   subroutine fetch_named_gradient(this, device_ptr, gradient, label, error)
      !! Pull an natom x 3 gradient from an explicitly named device buffer
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: device_ptr
      real(dp), intent(out) :: gradient(:, :)
      character(len=*), intent(in) :: label
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat(:)

      gradient = 0.0_dp
      if (error%has_error()) return

      allocate (flat(3*this%n_atoms))
      call copy_to_host(flat, device_ptr, label, error)
      if (.not. error%has_error()) gradient = reshape(flat, [3, int(this%n_atoms)])
      deallocate (flat)
   end subroutine fetch_named_gradient

   subroutine system_destroy(this)
      !! Release every cuEST object and device buffer, in reverse creation order
      !!
      !! Each persistent workspace is freed only after the object it backs has
      !! been destroyed -- cuEST keeps using that memory until then.
      class(cuest_system_t), intent(inout) :: this
      integer(c_int) :: status

      ! Borrowed from the context's pools -- drop the references, do not free.
      this%d_matrix = c_null_ptr
      this%d_c_occ = c_null_ptr
      this%d_c_occ_beta = c_null_ptr
      this%d_result_beta = c_null_ptr
      this%d_result = c_null_ptr
      this%d_gradient = c_null_ptr
      this%d_charge_gradient = c_null_ptr
      this%xyz_device = c_null_ptr
      this%charges_device = c_null_ptr
      this%d_j = c_null_ptr
      this%d_k = c_null_ptr
      this%d_xc = c_null_ptr
      this%d_fock = c_null_ptr
      this%d_core = c_null_ptr
      this%d_ovlp = c_null_ptr
      this%d_error = c_null_ptr
      this%d_transform = c_null_ptr
      this%d_commutator = c_null_ptr
      this%d_work = c_null_ptr
      this%d_fock_ortho = c_null_ptr
      this%d_orbitals = c_null_ptr
      this%d_eigenvalues = c_null_ptr
      this%d_solver = c_null_ptr
      this%d_devinfo = c_null_ptr
      this%d_density_alpha = c_null_ptr
      this%d_density_beta = c_null_ptr
      this%d_k_beta = c_null_ptr
      this%d_xc_beta = c_null_ptr
      this%d_fock_beta = c_null_ptr
      this%d_eigenvalues_beta = c_null_ptr

      ! The continuum, before the one-electron plan it was built on.
      call device_free(this%d_pcm)
      call device_free(this%d_q_in)
      call device_free(this%d_q_out)
      this%d_pcm = c_null_ptr
      this%d_q_in = c_null_ptr
      this%d_q_out = c_null_ptr
      if (c_associated(this%pcm_plan)) status = cuestPCMIntPlanDestroy(this%pcm_plan)
      this%pcm_plan = c_null_ptr
      call workspace_free(this%ws_pcm_plan)
      this%has_pcm = .false.
      this%n_pcm_points = 0
      this%n_pcm_active = 0

      if (c_associated(this%xc_plan)) status = cuestXCIntPlanDestroy(this%xc_plan)
      this%xc_plan = c_null_ptr
      call workspace_free(this%ws_xc_plan)

      if (c_associated(this%molecular_grid)) status = cuestMolecularGridDestroy(this%molecular_grid)
      this%molecular_grid = c_null_ptr
      call workspace_free(this%ws_grid)

      if (c_associated(this%df_plan)) status = cuestDFIntPlanDestroy(this%df_plan)
      this%df_plan = c_null_ptr
      call workspace_free(this%ws_df_plan)

      if (c_associated(this%oe_plan)) status = cuestOEIntPlanDestroy(this%oe_plan)
      this%oe_plan = c_null_ptr
      call workspace_free(this%ws_oe_plan)

      if (c_associated(this%pair_list)) status = cuestAOPairListDestroy(this%pair_list)
      this%pair_list = c_null_ptr
      call workspace_free(this%ws_pair_list)

      if (c_associated(this%aux_basis)) status = cuestAOBasisDestroy(this%aux_basis)
      this%aux_basis = c_null_ptr
      call workspace_free(this%ws_aux_basis)

      if (c_associated(this%basis)) status = cuestAOBasisDestroy(this%basis)
      this%basis = c_null_ptr
      call workspace_free(this%ws_basis)

      if (allocated(this%xyz_host)) deallocate (this%xyz_host)
      if (allocated(this%charges_host)) deallocate (this%charges_host)

      this%handle = c_null_ptr
      this%cublas = c_null_ptr
      this%cusolver = c_null_ptr
      this%solver_lwork = 0
      this%n_atoms = 0
      this%n_ao = 0
      this%n_occ = 0
   end subroutine system_destroy

end module mqc_cuest_integrals
