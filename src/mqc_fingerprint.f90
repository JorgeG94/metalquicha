!! Identifies the calculation that produced a number
module mqc_fingerprint
   !! A short hash of everything that decides a fragment energy, so a stored
   !! energy can be matched against the calculation asking to reuse it.
   !!
   !! **This exists so that restart cannot lie.** Reusing energies is the one
   !! feature whose failure looks like success: a resume that quietly mixes
   !! terms from a different geometry, basis or functional finishes early and
   !! reports a converged total, and there is nothing in the number to see. A
   !! checkpoint carrying this hash can be refused instead.
   !!
   !! Two rules decide what goes in.
   !!
   !!   * **In, if it changes the number.** Coordinates, elements, charge and
   !!     multiplicity; the monomer partition, because a term names monomers by
   !!     index and a repartition silently renames every one of them; the
   !!     bonds, because they decide where caps go; the method, basis,
   !!     functional and the thresholds the SCF converges to.
   !!   * **Out, if it does not.** Log level, output paths, rank counts,
   !!     verbosity, GPU binding. These change between the run that wrote a
   !!     checkpoint and the run that resumes it as a matter of course, and a
   !!     hash that rejected those would be a hash nobody could use.
   !!
   !! When in doubt a field goes in. A false mismatch costs a recomputation; a
   !! false match costs a wrong answer that nothing will catch.
   !!
   !! Reals are hashed as their IEEE bits rather than as text, so no formatting
   !! choice can round two different geometries onto one hash. The single
   !! exception is negative zero, normalised to positive, because it compares
   !! equal to zero everywhere else in the program and would otherwise make an
   !! identical geometry hash differently depending on how it was built.
   use pic_types, only: dp, int32, int64
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_method_config, only: method_config_t
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, &
                               METHOD_TYPE_HF, METHOD_TYPE_DFT
   implicit none
   private

   public :: hasher_t
   public :: calculation_fingerprint
   public :: FINGERPRINT_LEN

   integer, parameter :: FINGERPRINT_LEN = 16
      !! Hex digits in a fingerprint: 64 bits

   integer(int64), parameter :: FNV_OFFSET = -3750763034362895579_int64
      !! 14695981039346656037 as a signed 64-bit value
   integer(int64), parameter :: FNV_PRIME = 1099511628211_int64

   type :: hasher_t
      !! FNV-1a, accumulated field by field
      !!
      !! Order matters and is part of the identity: the same values fed in a
      !! different order give a different hash, which is what makes two arrays
      !! of the same numbers in different roles distinguishable.
      integer(int64), private :: state = FNV_OFFSET
   contains
      procedure :: byte => hash_byte
      procedure :: int => hash_int
      procedure :: int_array => hash_int_array
      procedure :: real => hash_real
      procedure :: real_array => hash_real_array
      procedure :: text => hash_text
      procedure :: flag => hash_flag
      procedure :: hex => hash_hex
   end type hasher_t

contains

   function calculation_fingerprint(sys_geom, config, calc_type) result(hex)
      !! The identity of a calculation: what it computes, and on what
      !!
      !! Two runs agreeing here may share energies. Two runs disagreeing may
      !! not, and the difference need not be visible -- swapping a functional
      !! leaves every array shape identical.
      type(system_geometry_t), intent(in) :: sys_geom
      type(method_config_t), intent(in) :: config
      integer(int32), intent(in) :: calc_type
      character(len=FINGERPRINT_LEN) :: hex

      type(hasher_t) :: h

      call h%text("mqc-fingerprint-v1")
      call add_system(h, sys_geom)
      call add_method(h, config)
      call h%int(int(calc_type))
      hex = h%hex()
   end function calculation_fingerprint

   subroutine add_system(h, sys)
      !! Everything about the molecule and how it is cut up
      type(hasher_t), intent(inout) :: h
      type(system_geometry_t), intent(in) :: sys

      integer :: imon, islot

      call h%text("system")
      call h%int(sys%total_atoms)
      call h%int(sys%charge)
      call h%int(sys%multiplicity)
      if (allocated(sys%element_numbers)) call h%int_array(sys%element_numbers)
      if (allocated(sys%coordinates)) call h%real_array(reshape(sys%coordinates, &
                                                                [size(sys%coordinates)]))

      ! The partition. A term list names monomers by index, so two runs that
      ! agree on the atoms and disagree on the grouping agree on nothing that
      ! matters -- term (1,2) is a different dimer in each.
      call h%text("monomers")
      call h%int(sys%n_monomers)
      if (allocated(sys%fragment_sizes)) then
         call h%int_array(sys%fragment_sizes)
         do imon = 1, sys%n_monomers
            do islot = 1, sys%fragment_sizes(imon)
               call h%int(sys%fragment_atoms(islot, imon))
            end do
         end do
      end if
      if (allocated(sys%fragment_charges)) call h%int_array(sys%fragment_charges)
      if (allocated(sys%fragment_multiplicities)) then
         call h%int_array(sys%fragment_multiplicities)
      end if

      ! Bonds decide where caps go, and a cap is an atom that was not there
      ! before. `is_broken` is derived from the partition, so it is not hashed
      ! separately -- the partition above already covers it.
      call h%text("bonds")
      if (allocated(sys%bonds)) then
         call h%int(size(sys%bonds))
         do imon = 1, size(sys%bonds)
            call h%int(sys%bonds(imon)%atom_i)
            call h%int(sys%bonds(imon)%atom_j)
            call h%int(sys%bonds(imon)%order)
         end do
      else
         call h%int(0)
      end if
   end subroutine add_system

   subroutine add_method(h, config)
      !! Everything that decides what number a fragment gets
      !!
      !! Only the branch that will actually run is hashed. A deck carrying
      !! stale `dft` settings while running xTB would otherwise fingerprint
      !! differently from the same xTB calculation written cleanly, and the two
      !! compute the same thing.
      type(hasher_t), intent(inout) :: h
      type(method_config_t), intent(in) :: config

      call h%text("method")
      call h%int(int(config%method_type))

      select case (config%method_type)
      case (METHOD_TYPE_GFN1, METHOD_TYPE_GFN2)
         call h%text("xtb")
         call h%real(config%xtb%accuracy)
         call h%real(config%xtb%electronic_temp)
         call h%text(trim(config%xtb%solvent))
         call h%text(trim(config%xtb%solvation_model))
         call h%flag(config%xtb%use_cds)
         call h%flag(config%xtb%use_shift)
         call h%real(config%xtb%dielectric)
         call h%int(config%xtb%cpcm_nang)
         call h%real(config%xtb%cpcm_rscale)
         ! maxiter is in because it decides whether the SCF converged, and an
         ! energy from a run that stopped early is not the same number.
         call h%int(config%scf%max_iter)

      case default
         ! HF, DFT and anything else built on the SCF.
         call h%text(trim(config%basis_set))
         call h%flag(config%use_spherical)
         call h%text(trim(config%scf%aux_basis_set))
         ! Whether it was asked for, not only what it is. A double hybrid fits
         ! its perturbative term when a deck names an auxiliary basis and not
         ! when one merely defaulted, so two configs carrying the same string
         ! for different reasons are two different energies.
         call h%flag(config%scf%aux_basis_named)
         call h%text(trim(config%scf%guess))
         call h%flag(config%scf%unrestricted)
         call h%int(config%scf%max_iter)
         call h%real(config%scf%energy_convergence)
         call h%real(config%scf%density_convergence)
         if (config%method_type == METHOD_TYPE_DFT) then
            call h%text("dft")
            call h%text(trim(config%dft%functional))
            call h%text(trim(config%dft%grid_type))
            call h%int(config%dft%radial_points)
            call h%int(config%dft%angular_points)
         end if
      end select
   end subroutine add_method

   ! ==========================================================================
   !  The hasher
   ! ==========================================================================

   subroutine hash_byte(this, value)
      !! Fold one byte in
      class(hasher_t), intent(inout) :: this
      integer(int64), intent(in) :: value

      this%state = ieor(this%state, iand(value, 255_int64))
      this%state = this%state*FNV_PRIME
   end subroutine hash_byte

   subroutine hash_int(this, value)
      class(hasher_t), intent(inout) :: this
      integer, intent(in) :: value

      integer :: i
      integer(int64) :: wide

      wide = int(value, int64)
      do i = 0, 7
         call this%byte(ishft(wide, -8*i))
      end do
   end subroutine hash_int

   subroutine hash_int_array(this, values)
      class(hasher_t), intent(inout) :: this
      integer, intent(in) :: values(:)

      integer :: i

      call this%int(size(values))
      do i = 1, size(values)
         call this%int(values(i))
      end do
   end subroutine hash_int_array

   subroutine hash_real(this, value)
      !! Fold a double in by its bits, not its digits
      class(hasher_t), intent(inout) :: this
      real(dp), intent(in) :: value

      integer :: i
      integer(int64) :: bits
      real(dp) :: normalised

      ! -0.0 and 0.0 are equal everywhere else in this program, so they must
      ! not produce different fingerprints. Adding zero collapses the two.
      normalised = value
      if (normalised == 0.0_dp) normalised = 0.0_dp
      bits = transfer(normalised, bits)
      do i = 0, 7
         call this%byte(ishft(bits, -8*i))
      end do
   end subroutine hash_real

   subroutine hash_real_array(this, values)
      class(hasher_t), intent(inout) :: this
      real(dp), intent(in) :: values(:)

      integer :: i

      call this%int(size(values))
      do i = 1, size(values)
         call this%real(values(i))
      end do
   end subroutine hash_real_array

   subroutine hash_text(this, value)
      class(hasher_t), intent(inout) :: this
      character(len=*), intent(in) :: value

      integer :: i

      call this%int(len(value))
      do i = 1, len(value)
         call this%byte(int(iachar(value(i:i)), int64))
      end do
   end subroutine hash_text

   subroutine hash_flag(this, value)
      class(hasher_t), intent(inout) :: this
      logical, intent(in) :: value

      call this%int(merge(1, 0, value))
   end subroutine hash_flag

   function hash_hex(this) result(hex)
      !! The accumulated hash, as 16 lowercase hex digits
      class(hasher_t), intent(in) :: this
      character(len=FINGERPRINT_LEN) :: hex

      character(len=*), parameter :: DIGITS = "0123456789abcdef"
      integer :: i, nibble

      do i = 1, FINGERPRINT_LEN
         nibble = int(iand(ishft(this%state, -4*(FINGERPRINT_LEN - i)), 15_int64))
         hex(i:i) = DIGITS(nibble + 1:nibble + 1)
      end do
   end function hash_hex

end module mqc_fingerprint
