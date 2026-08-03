!! Conversion of a metalquicha molecular basis into cuEST AO shell handles
module mqc_cuest_basis
   !! Bridges `mqc_cgto`'s `molecular_basis_type` to the array of
   !! `cuestAOShell_t` handles cuEST needs to build an AO basis.
   !!
   !! `molecular_basis_type%elements` is already indexed per *atom* in geometry
   !! order, which is exactly the layout cuEST's `numShellsPerAtom` expects, so
   !! the mapping is direct. Contraction coefficients are normalized on the way
   !! through -- basis files store unnormalized ones.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_int32_t, &
                                                                             c_int64_t, c_double, c_loc, c_associated
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type
   use mqc_basis_normalization, only: normalized_coefficients
   use mqc_cuest_runtime, only: cuest_status_check
   use cuest, only: cuestAOShellCreate, cuestAOShellDestroy, &
                    cuestParametersCreate, cuestParametersDestroy, &
                    CUEST_AOSHELL_PARAMETERS
   implicit none
   private

   public :: cuest_shell_set_t   !! Owned collection of cuEST AO shell handles
   public :: build_cuest_shells  !! molecular_basis_type -> AO shell handles

   type :: cuest_shell_set_t
      !! The AO shell handles for one molecule, plus the per-atom counts
      !!
      !! Handles are owned by this type: `destroy` releases every one of them.
      !! The set must stay alive at least until the AO basis is created.
      integer(c_int64_t) :: n_atoms = 0        !! Number of atoms
      integer(c_int64_t) :: n_shells_total = 0  !! Total shells over all atoms
      integer(c_int64_t), allocatable :: n_shells_per_atom(:)  !! Shells on each atom
      type(c_ptr), allocatable :: shells(:)    !! cuestAOShell_t handles
   contains
      procedure :: destroy => shell_set_destroy
   end type cuest_shell_set_t

contains

   subroutine build_cuest_shells(handle, mol_basis, use_spherical, shell_set, error)
      !! Create the cuEST AO shell handles for a molecular basis
      type(c_ptr), intent(in) :: handle                   !! cuEST handle
      type(molecular_basis_type), intent(in) :: mol_basis  !! Basis, one entry per atom
      logical, intent(in) :: use_spherical                 !! Pure (spherical) vs Cartesian
      type(cuest_shell_set_t), intent(out) :: shell_set
      type(error_t), intent(inout) :: error

      type(c_ptr) :: shell_params
      integer(c_int32_t) :: is_pure
      integer(c_int) :: status
      integer :: iatom, ishell, shell_index, n_prim
      real(dp), allocatable, target :: exponents(:), coefficients(:)

      is_pure = merge(1_c_int32_t, 0_c_int32_t, use_spherical)

      if (mol_basis%nelements <= 0) then
         call error%set(ERROR_VALIDATION, "cuEST basis build: molecular basis is empty")
         return
      end if

      ! ---- per-atom shell counts -------------------------------------------
      shell_set%n_atoms = int(mol_basis%nelements, c_int64_t)
      allocate (shell_set%n_shells_per_atom(mol_basis%nelements))
      shell_set%n_shells_total = 0
      do iatom = 1, mol_basis%nelements
         shell_set%n_shells_per_atom(iatom) = int(mol_basis%elements(iatom)%nshells, c_int64_t)
         shell_set%n_shells_total = shell_set%n_shells_total + shell_set%n_shells_per_atom(iatom)
      end do

      allocate (shell_set%shells(shell_set%n_shells_total))
      shell_set%shells = c_null_ptr

      ! Even object types with no configurable parameters need a live
      ! parameters handle -- passing a null pointer returns NULL_POINTER.
      shell_params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_AOSHELL_PARAMETERS, shell_params), &
                              "cuestParametersCreate(AO shell)", error)
      if (error%has_error()) return

      ! ---- create one handle per shell, in atom order ----------------------
      shell_index = 0
      atom_loop: do iatom = 1, mol_basis%nelements
         do ishell = 1, mol_basis%elements(iatom)%nshells
            n_prim = mol_basis%elements(iatom)%shells(ishell)%nfunc

            allocate (exponents(n_prim), coefficients(n_prim))
            exponents = mol_basis%elements(iatom)%shells(ishell)%exponents(1:n_prim)

            call normalized_coefficients(mol_basis%elements(iatom)%shells(ishell)%ang_mom, &
                                         exponents, &
                                         mol_basis%elements(iatom)%shells(ishell)%coefficients(1:n_prim), &
                                         1.0_dp, coefficients, error)
            if (error%has_error()) then
               deallocate (exponents, coefficients)
               exit atom_loop
            end if

            shell_index = shell_index + 1
            call create_shell(handle, is_pure, &
                              int(mol_basis%elements(iatom)%shells(ishell)%ang_mom, c_int64_t), &
                              int(n_prim, c_int64_t), exponents, coefficients, &
                              shell_params, shell_set%shells(shell_index), error)

            deallocate (exponents, coefficients)
            if (error%has_error()) exit atom_loop
         end do
      end do atom_loop

      status = cuestParametersDestroy(CUEST_AOSHELL_PARAMETERS, shell_params)
      call cuest_status_check(status, "cuestParametersDestroy(AO shell)", error)

      if (error%has_error()) then
         call error%add_context("mqc_cuest_basis:build_cuest_shells")
         call shell_set%destroy()
      end if
   end subroutine build_cuest_shells

   subroutine create_shell(handle, is_pure, ang_mom, n_prim, exponents, coefficients, &
                           shell_params, shell, error)
      !! Create one cuEST AO shell from host exponent/coefficient arrays
      !!
      !! Wrapped separately because C_LOC needs an actual argument with the
      !! TARGET attribute; these are HOST arrays, not device buffers.
      type(c_ptr), intent(in) :: handle, shell_params
      integer(c_int32_t), intent(in) :: is_pure
      integer(c_int64_t), intent(in) :: ang_mom, n_prim
      real(dp), intent(in), target :: exponents(:), coefficients(:)
      type(c_ptr), intent(inout) :: shell
      type(error_t), intent(inout) :: error

      call cuest_status_check(cuestAOShellCreate(handle, is_pure, ang_mom, n_prim, &
                                                 c_loc(exponents), c_loc(coefficients), &
                                                 shell_params, shell), &
                              "cuestAOShellCreate", error)
   end subroutine create_shell

   subroutine shell_set_destroy(this)
      !! Release every AO shell handle held by the set
      class(cuest_shell_set_t), intent(inout) :: this
      integer(c_int) :: status
      integer :: i

      if (allocated(this%shells)) then
         do i = 1, size(this%shells)
            if (c_associated(this%shells(i))) status = cuestAOShellDestroy(this%shells(i))
         end do
         deallocate (this%shells)
      end if
      if (allocated(this%n_shells_per_atom)) deallocate (this%n_shells_per_atom)
      this%n_atoms = 0
      this%n_shells_total = 0
   end subroutine shell_set_destroy

end module mqc_cuest_basis
