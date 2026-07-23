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
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_int64_t, &
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
      integer(c_int64_t) :: n_occ = 0    !! Doubly occupied orbitals

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
      type(c_ptr) :: d_c_occ = c_null_ptr   !! Occupied MO coefficients
      type(c_ptr) :: d_result = c_null_ptr  !! Whichever matrix is being built
   contains
      procedure :: create => system_create
      procedure :: destroy => system_destroy
      procedure :: compute_overlap => system_compute_overlap
      procedure :: compute_kinetic => system_compute_kinetic
      procedure :: compute_potential => system_compute_potential
      procedure :: compute_coulomb => system_compute_coulomb
      procedure :: compute_exchange => system_compute_exchange
      procedure :: compute_xc => system_compute_xc
   end type cuest_system_t

contains

   subroutine system_create(this, context, atomic_numbers, coordinates, mol_basis, &
                            aux_mol_basis, use_spherical, n_occ, functional_id, &
                            n_radial, n_angular, error)
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

      integer :: iatom, n_atoms

      n_atoms = size(atomic_numbers)
      this%handle = context%handle
      this%n_atoms = int(n_atoms, c_int64_t)
      this%n_occ = int(n_occ, c_int64_t)
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

      if (this%n_occ > 0) then
         call context%scratch_c_occ%ensure(this%n_ao*this%n_occ, &
                                           "occupied MO coefficients", error)
         if (error%has_error()) then
            call error%add_context("mqc_cuest_integrals:create")
            return
         end if
         this%d_c_occ = context%scratch_c_occ%ptr
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

   subroutine system_compute_exchange(this, c_occ, exchange, error)
      !! Density-fitted exchange matrix K from occupied MO coefficients
      class(cuest_system_t), intent(inout) :: this
      real(dp), intent(in) :: c_occ(:, :)      !! (n_ao, n_occ)
      real(dp), intent(out) :: exchange(:, :)  !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error

      type(cuestWorkspaceDescriptor_t) :: temporary_desc, variable_buffer
      type(cuestWorkspace_t) :: temporary_ws
      type(c_ptr) :: params
      integer(c_int) :: status
      real(dp), allocatable :: flat(:)

      if (error%has_error()) return

      if (this%n_occ <= 0) then
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
                                                                           this%n_occ, this%d_c_occ, &
                                                                           this%d_result), &
                              "cuestDFSymmetricExchangeComputeWorkspaceQuery", error)
      if (.not. error%has_error()) call workspace_alloc(temporary_ws, temporary_desc, error)
      if (.not. error%has_error()) then
         call cuest_status_check(cuestDFSymmetricExchangeCompute(this%handle, this%df_plan, &
                                                                 params, variable_buffer, &
                                                                 temporary_ws, this%n_occ, &
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
      this%d_result = c_null_ptr
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
