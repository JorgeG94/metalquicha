!! Stand-in for the CPU SCF when the build has no libcint
module mqc_libcint_bridge
   !! Same name and same entry points as the real bridge, declining.
   !!
   !! A stub rather than a preprocessor guard at the call site, matching
   !! mqc_cuest_bridge_stub next door: the method layer reads the same either
   !! way, and `libcint_backend_available` is how it asks in advance so the
   !! refusal can name the build option.
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_error, only: ERROR_VALIDATION
   use mqc_cuest_iface, only: cuest_scf_settings_t
   implicit none
   private

   public :: run_libcint_hf
   public :: run_libcint_mcscf
   public :: run_libcint_fmo
   public :: run_libcint_makefp
   public :: run_libcint_charges
   public :: run_libcint_efp
   public :: run_libcint_sapt0
   public :: libcint_backend_available
   public :: xc_available

contains

   pure function xc_available() result(available)
      !! No libcint means no `mqc_libcint_xc`, and therefore no functionals --
      !! whatever libxc itself was configured to do.
      logical :: available
      available = .false.
   end function xc_available

   pure function libcint_backend_available() result(available)
      !! Whether this build can run an SCF on the CPU
      logical :: available
      available = .false.
   end function libcint_backend_available

   subroutine run_libcint_sapt0(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, basis_name, &
                                terms, error)
      !! No-op stand-in: SAPT0 needs the CPU integral backend
      use pic_types, only: dp
      use mqc_program_limits, only: N_SAPT_TERMS
      use mqc_error, only: error_t
      integer, intent(in) :: z_a(:), z_b(:)
      character(len=*), intent(in) :: sym_a(:), sym_b(:)
      real(dp), intent(in) :: xyz_a(:, :), xyz_b(:, :)
      character(len=*), intent(in) :: basis_name
      real(dp), intent(out) :: terms(N_SAPT_TERMS)
      type(error_t), intent(inout) :: error

      terms = 0.0_dp
      call error%set(ERROR_VALIDATION, &
                     "SAPT needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_LIBCINT=ON")
      if (size(z_a) < 0 .or. size(z_b) < 0) return
      if (len_trim(sym_a(1))*len_trim(sym_b(1))*len_trim(basis_name) < 0) return
      if (size(xyz_a) < 0 .or. size(xyz_b) < 0) return
   end subroutine run_libcint_sapt0

   subroutine run_libcint_charges(atomic_numbers, element_symbols, coordinates, &
                                  basis_name, scheme, total_charge, charges, error)
      !! No-op stand-in: atomic charges need the CPU integral backend
      !!
      !! The molecule builder, the SCF and both partition schemes all live behind
      !! `MQC_ENABLE_LIBCINT`.
      use pic_types, only: dp
      use mqc_error, only: error_t
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      character(len=*), intent(in) :: basis_name
      character(len=*), intent(in) :: scheme
      integer, intent(in) :: total_charge
      real(dp), allocatable, intent(out) :: charges(:)
      type(error_t), intent(inout) :: error

      allocate (charges(0))
      call error%set(ERROR_VALIDATION, &
                     "atomic charges need the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_LIBCINT=ON")
      if (size(atomic_numbers) < 0 .or. size(coordinates) < 0) return
      if (len_trim(element_symbols(1))*len_trim(basis_name)*len_trim(scheme) < 0) return
      if (total_charge < -huge(1)) return
   end subroutine run_libcint_charges

   subroutine run_libcint_efp(potentials, fragment_sizes, fragment_atoms, &
                              coordinates, terms, error)
      !! No-op stand-in: an EFP interaction energy needs the CPU backend
      !!
      !! Every piece of it does -- the multipole and polarizability machinery, the
      !! integrals two fragments' basis sets share, and the potential reader itself
      !! all live behind `MQC_ENABLE_LIBCINT`.
      use pic_types, only: dp
      use mqc_program_limits, only: N_EFP_TERMS
      use mqc_error, only: error_t
      character(len=*), intent(in) :: potentials(:)
      integer, intent(in) :: fragment_sizes(:)
      integer, intent(in) :: fragment_atoms(:, :)
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), intent(out) :: terms(N_EFP_TERMS)
      type(error_t), intent(inout) :: error

      terms = 0.0_dp
      call error%set(ERROR_VALIDATION, &
                     "EFP needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_LIBCINT=ON")
      if (len_trim(potentials(1)) < 0) return
      if (size(fragment_sizes) < 0 .or. size(fragment_atoms) < 0) return
      if (size(coordinates) < 0) return
   end subroutine run_libcint_efp
   subroutine run_libcint_fmo(atomic_numbers, element_symbols, coordinates, owner, &
                              basis_name, esp, expansion, far_field, resppc, &
                              level, max_outer, outer_tol, scf_max_iter, &
                              scf_energy_tol, scf_density_tol, energy, error, comm)
      !! Run FMO2 (or EE-MBE) over a partitioned system
      !!
      !! Options arrive as plain scalars rather than the backend's own options
      !! type, so the layer above never has to see a type it cannot compile
      !! without the backend. Coordinates are Bohr; `owner(i)` is atom i's
      !! fragment, numbered from one with no gaps.
      use pic_types, only: dp
      use mqc_error, only: error_t
      use pic_mpi_lib, only: comm_t
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      integer, intent(in) :: owner(:)
      character(len=*), intent(in) :: basis_name, esp, expansion, far_field
      real(dp), intent(in) :: resppc
      integer, intent(in) :: level
      integer, intent(in) :: max_outer
      real(dp), intent(in) :: outer_tol
      integer, intent(in) :: scf_max_iter
      real(dp), intent(in) :: scf_energy_tol, scf_density_tol
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm
         !! Present means distribute the fragment work over this communicator.
         !! Absent means one rank does all of it.

      energy = 0.0_dp
      call error%set(ERROR_VALIDATION, &
                     "FMO needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_LIBCINT=ON")
      if (size(atomic_numbers) < 0 .or. size(coordinates) < 0 .or. size(owner) < 0) return
      if (len_trim(element_symbols(1)) < 0) return
      if (len_trim(basis_name)*len_trim(esp)*len_trim(expansion)*len_trim(far_field) < 0) return
      if (resppc < -huge(1.0_dp) .or. level < 0 .or. max_outer < 0) return
      if (outer_tol < 0.0_dp .or. scf_max_iter < 0 .or. scf_energy_tol < 0.0_dp) return
      if (scf_density_tol < 0.0_dp) return
      if (present(comm)) return
   end subroutine run_libcint_fmo

   subroutine run_libcint_makefp(atomic_numbers, element_symbols, coordinates, &
                                 basis_name, name, path, error, charge, verbose, &
                                 aux_basis, guess)
      !! No-op stand-in: an effective fragment potential needs the CPU backend
      use pic_types, only: dp
      use mqc_error, only: error_t
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      character(len=*), intent(in) :: basis_name, name, path
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: charge   !! Net charge; a fragment may be an ion
      logical, intent(in), optional :: verbose
      character(len=*), intent(in), optional :: aux_basis
      character(len=*), intent(in), optional :: guess

      call error%set(ERROR_VALIDATION, &
                     "MAKEFP needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_LIBCINT=ON")
      if (size(atomic_numbers) < 0 .or. size(coordinates) < 0) return
      if (len_trim(element_symbols(1)) < 0) return
      if (len_trim(basis_name)*len_trim(name)*len_trim(path) < 0) return
      if (present(charge) .or. present(verbose)) return
      if (present(aux_basis) .or. present(guess)) return
   end subroutine run_libcint_makefp

   subroutine run_libcint_hf(settings, fragment, result, want_gradient, want_hessian)
      !! No-op stand-in: report the missing backend, compute nothing
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient
      logical, intent(in), optional :: want_hessian

      call result%error%set(ERROR_VALIDATION, &
                            "This calculation needs an integral backend; build with "// &
                            "-DMQC_ENABLE_LIBCINT=ON for the CPU one, or "// &
                            "-DMQC_ENABLE_CUEST=ON for the GPU one")
      result%has_error = .true.
      result%has_energy = .false.
      if (len_trim(settings%basis_set) == 0 .or. fragment%n_atoms < 0) return
      if (present(want_gradient)) return
      if (present(want_hessian)) return
   end subroutine run_libcint_hf

   subroutine run_libcint_mcscf(settings, fragment, result)
      !! No-op stand-in: CASSCF and CASCI need the CPU integral backend
      !!
      !! All of it does -- the reference SCF, the active-space transform, the
      !! determinant machinery and the orbital optimiser are all behind
      !! `MQC_ENABLE_LIBCINT`, and there is no GPU path to fall through to.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result

      call result%error%set(ERROR_VALIDATION, &
                            "a multiconfigurational calculation needs the CPU integral "// &
                            "backend; build with -DMQC_ENABLE_LIBCINT=ON")
      result%has_error = .true.
      result%has_energy = .false.
      if (len_trim(settings%basis_set) == 0 .or. fragment%n_atoms < 0) return
   end subroutine run_libcint_mcscf

end module mqc_libcint_bridge
