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
   use mqc_cgto, only: molecular_basis_type
   use mqc_cuest_basis, only: cuest_shell_set_t, build_cuest_shells
   use mqc_cuest_grid, only: atom_grid_set_t, build_atom_grids
   use mqc_cuest_context, only: cuest_context_t
   use mqc_cuest_runtime, only: cuest_status_check, workspace_alloc, workspace_free, &
                                copy_to_device, copy_to_host, device_sync
   use cuest, only: cuestWorkspace_t, cuestWorkspaceDescriptor_t, &
                    cuestParametersCreate, cuestParametersDestroy, &
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
   use cuest_helpers, only: cuest_query_i64, cuest_query_f64, cuest_param_set_f64
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

      ! cuEST objects, each with the persistent workspace that must outlive it
      type(c_ptr) :: basis = c_null_ptr       !! Primary (orbital) AO basis
      type(c_ptr) :: aux_basis = c_null_ptr   !! Auxiliary (fitting) AO basis
      type(c_ptr) :: pair_list = c_null_ptr   !! AO pair list, primary basis
      type(c_ptr) :: oe_plan = c_null_ptr     !! One-electron integral plan
      type(c_ptr) :: df_plan = c_null_ptr     !! Density-fitted integral plan
      type(c_ptr) :: molecular_grid = c_null_ptr  !! XC quadrature grid (DFT only)
      type(c_ptr) :: xc_plan = c_null_ptr     !! XC integral plan (DFT only)
      type(cuestWorkspace_t) :: ws_basis, ws_aux_basis, ws_pair_list
      type(cuestWorkspace_t) :: ws_oe_plan, ws_df_plan
      type(cuestWorkspace_t) :: ws_grid, ws_xc_plan

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
      type(c_ptr) :: d_charge_gradient = c_null_ptr !! Hellmann-Feynman half
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
      procedure :: gradient_overlap => system_gradient_overlap
      procedure :: gradient_kinetic => system_gradient_kinetic
      procedure :: gradient_potential => system_gradient_potential
      procedure :: gradient_two_electron => system_gradient_two_electron
      procedure :: gradient_xc => system_gradient_xc
      procedure :: gradient_two_electron_uks => system_gradient_two_electron_uks
      procedure :: gradient_xc_uks => system_gradient_xc_uks
   end type cuest_system_t

contains

   subroutine system_create(this, context, atomic_numbers, coordinates, mol_basis, &
                            aux_mol_basis, use_spherical, n_occ, functional_id, &
                            n_radial, n_angular, error, n_occ_beta, n_guess_columns)
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
      integer(c_int64_t) :: widest

      n_atoms = size(atomic_numbers)
      this%handle = context%handle
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

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      integer(c_int64_t) :: beta_occupancy
      real(dp), allocatable :: flat(:)
      real(dp), target :: energy

      xc_energy = 0.0_dp
      vxc_alpha = 0.0_dp
      vxc_beta = 0.0_dp
      if (error%has_error() .or. .not. this%has_xc) return

      allocate (flat(size(c_occ_alpha)))
      flat = reshape(c_occ_alpha, [size(c_occ_alpha)])
      call copy_to_device(this%d_c_occ, flat, "alpha occupied MOs (XC)", error)
      deallocate (flat)

      ! cuEST requires numOccupiedBeta > 0, but a system can legitimately have
      ! no beta electrons -- a hydrogen atom, for one, which is exactly what an
      ! atomic guess has to solve. Passing a single all-zero column gives a
      ! beta density of zero, which is the correct physics, through an argument
      ! the library will accept.
      beta_occupancy = max(this%n_occ_beta, 1_c_int64_t)
      allocate (flat(this%n_ao*beta_occupancy))
      flat = 0.0_dp
      if (this%n_occ_beta > 0) flat = reshape(c_occ_beta, [size(c_occ_beta)])
      call copy_to_device(this%d_c_occ_beta, flat, "beta occupied MOs (XC)", error)
      deallocate (flat)
      if (error%has_error()) return

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
                                                                       beta_occupancy, this%d_c_occ, &
                                                                       this%d_c_occ_beta, c_loc(energy), &
                                                                       this%d_result, this%d_result_beta), &
                              "cuestXCPotentialUKSComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestXCPotentialUKSCompute(this%handle, this%xc_plan, params, &
                                                            variable_buffer, temporary_ws, &
                                                            this%n_occ, beta_occupancy, &
                                                            this%d_c_occ, this%d_c_occ_beta, &
                                                            c_loc(energy), this%d_result, &
                                                            this%d_result_beta), &
                                 "cuestXCPotentialUKSCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_XCPOTENTIALUKSCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(UKS XC potential)", error)

      call fetch_matrix(this, vxc_alpha, "Vxc alpha", error)
      call fetch_named_matrix(this, this%d_result_beta, vxc_beta, "Vxc beta", error)
      if (.not. error%has_error()) xc_energy = energy
   end subroutine system_compute_xc_uks

   subroutine fetch_named_matrix(this, device_ptr, matrix, label, error)
      !! Copy an n_ao x n_ao matrix back from an explicitly named device buffer
      class(cuest_system_t), intent(inout) :: this
      type(c_ptr), intent(in) :: device_ptr
      real(dp), intent(out) :: matrix(:, :)
      character(len=*), intent(in) :: label
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat(:)

      matrix = 0.0_dp
      if (error%has_error()) return

      allocate (flat(this%n_ao*this%n_ao))
      call copy_to_host(flat, device_ptr, label, error)
      if (.not. error%has_error()) matrix = reshape(flat, [int(this%n_ao), int(this%n_ao)])
      deallocate (flat)
   end subroutine fetch_named_matrix

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

   subroutine system_compute_xc(this, c_occ, xc_energy, xc_potential, error)
      !! Exchange-correlation energy and potential from occupied MO coefficients
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: c_occ(:, :)           !! (n_ao, n_occ)
      real(dp), intent(out) :: xc_energy            !! Exc, Hartree
      real(dp), intent(out) :: xc_potential(:, :)   !! Vxc, (n_ao, n_ao)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      real(dp), allocatable :: flat(:)
      real(dp), target :: energy

      xc_energy = 0.0_dp
      if (error%has_error()) then
         xc_potential = 0.0_dp
         return
      end if
      if (.not. this%has_xc) then
         xc_potential = 0.0_dp
         return
      end if

      allocate (flat(size(c_occ)))
      flat = reshape(c_occ, [size(c_occ)])
      call copy_to_device(this%d_c_occ, flat, "occupied MO coefficients (XC)", error)
      deallocate (flat)
      if (error%has_error()) return

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
                                                                       this%d_c_occ, c_loc(energy), &
                                                                       this%d_result), &
                              "cuestXCPotentialRKSComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestXCPotentialRKSCompute(this%handle, this%xc_plan, params, &
                                                            variable_buffer, temporary_ws, &
                                                            this%n_occ, this%d_c_occ, &
                                                            c_loc(energy), this%d_result), &
                                 "cuestXCPotentialRKSCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_XCPOTENTIALRKSCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(RKS XC potential)", error)

      call fetch_matrix(this, xc_potential, "Vxc", error)
      if (.not. error%has_error()) xc_energy = energy
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

   subroutine system_compute_coulomb(this, density, coulomb, error)
      !! Density-fitted Coulomb matrix J from a density matrix
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: density(:, :)   !! (n_ao, n_ao), symmetric
      real(dp), intent(out) :: coulomb(:, :)  !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      real(dp), allocatable :: flat(:)

      if (error%has_error()) return

      allocate (flat(size(density)))
      flat = reshape(density, [size(density)])
      call copy_to_device(this%d_matrix, flat, "density matrix", error)
      deallocate (flat)
      if (error%has_error()) return

      params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_DFCOULOMBCOMPUTE_PARAMETERS, params), &
                              "cuestParametersCreate(DF Coulomb)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestDFCoulombComputeWorkspaceQuery(this%handle, this%df_plan, &
                                                                  params, temporary_desc, &
                                                                  this%d_matrix, this%d_result), &
                              "cuestDFCoulombComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestDFCoulombCompute(this%handle, this%df_plan, params, &
                                                       temporary_ws, this%d_matrix, this%d_result), &
                                 "cuestDFCoulombCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_DFCOULOMBCOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(DF Coulomb)", error)

      call fetch_matrix(this, coulomb, "J", error)
   end subroutine system_compute_coulomb

   subroutine system_compute_exchange(this, c_occ, exchange, error, n_occupied)
      !! Density-fitted exchange matrix K from occupied MO coefficients
      !!
      !! `n_occupied` overrides the stored count, which is what the two spin
      !! channels of an unrestricted calculation need.
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: c_occ(:, :)      !! (n_ao, n_occ)
      real(dp), intent(out) :: exchange(:, :)  !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_occupied
      integer(c_int64_t) :: occupancy

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      real(dp), allocatable :: flat(:)

      if (error%has_error()) return

      occupancy = this%n_occ
      if (present(n_occupied)) occupancy = int(n_occupied, c_int64_t)

      if (occupancy <= 0) then
         exchange = 0.0_dp
         return
      end if

      allocate (flat(size(c_occ)))
      flat = reshape(c_occ, [size(c_occ)])
      call copy_to_device(this%d_c_occ, flat, "occupied MO coefficients", error)
      deallocate (flat)
      if (error%has_error()) return

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
                                                                            occupancy, this%d_c_occ, &
                                                                            this%d_result), &
                              "cuestDFSymmetricExchangeComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestDFSymmetricExchangeCompute(this%handle, this%df_plan, &
                                                                 params, variable_buffer, &
                                                                 temporary_ws, occupancy, &
                                                                 this%d_c_occ, this%d_result), &
                                 "cuestDFSymmetricExchangeCompute", error)
      end if

      call workspace_free(temporary_ws)
      status = cuestParametersDestroy(CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS, params)
      call cuest_status_check(status, "cuestParametersDestroy(DF exchange)", error)

      call fetch_matrix(this, exchange, "K", error)
   end subroutine system_compute_exchange

   subroutine fetch_matrix(this, matrix, label, error)
      !! Synchronize, then copy the result buffer back into a host matrix
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(out) :: matrix(:, :)
      character(len=*), intent(in) :: label
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat(:)

      if (error%has_error()) then
         matrix = 0.0_dp
         return
      end if

      call device_sync(label, error)
      if (error%has_error()) return

      allocate (flat(this%n_ao*this%n_ao))
      call copy_to_host(flat, this%d_result, label, error)
      if (.not. error%has_error()) then
         matrix = reshape(flat, [this%n_ao, this%n_ao])
      else
         matrix = 0.0_dp
      end if
      deallocate (flat)
   end subroutine fetch_matrix

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

      call stage_matrix(this, weighted_density, "energy-weighted density", error)
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

      call stage_matrix(this, density, "density (kinetic gradient)", error)
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

      call stage_matrix(this, density, "density (potential gradient)", error)
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

      call stage_matrix(this, half_density, "half density (JK gradient)", error)
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
      integer(c_int64_t) :: occupancies(2), n_matrices
      real(dp) :: coefficient_scale
      real(dp), allocatable :: flat(:)
      integer :: n_a, n_b

      if (error%has_error()) then
         gradient = 0.0_dp
         return
      end if

      n_a = int(this%n_occ)
      n_b = int(this%n_occ_beta)

      call stage_matrix(this, total_density, "total density (UKS JK gradient)", error)
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

   subroutine stage_matrix(this, matrix, label, error)
      !! Copy an n_ao x n_ao host matrix into the shared device input buffer
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: matrix(:, :)
      character(len=*), intent(in) :: label
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat(:)

      allocate (flat(size(matrix)))
      flat = reshape(matrix, [size(matrix)])
      call copy_to_device(this%d_matrix, flat, label, error)
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
      this%n_atoms = 0
      this%n_ao = 0
      this%n_occ = 0
   end subroutine system_destroy

end module mqc_cuest_integrals
