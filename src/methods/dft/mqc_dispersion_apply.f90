!! Which dispersion library serves which correction, in one place
module mqc_dispersion_apply
   !! The only module that knows the mapping from a `keywords.dft.dispersion`
   !! value to a library.
   !!
   !! There are two of those libraries now and they are not interchangeable:
   !! s-dftd3 computes "d3bj" from the nuclei alone, dftd4 computes "d4" from
   !! the nuclei *and the total molecular charge*, each is behind its own CMake
   !! option, and a build may have either, both or neither. Three places would
   !! otherwise need that table -- the reader, the method, and the version
   !! banner -- and three copies of a routing table is how a correction comes to
   !! be accepted by a deck and then not applied.
   !!
   !! Always compiled. The wrappers behind it are the pieces that come in two
   !! forms; this one calls whichever pair was built and reports what it found.
   use pic_types, only: dp, default_int
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_dispersion_names, only: dispersion_kind_is_d3, dispersion_kind_is_d4, DISPERSION_KINDS
   use mqc_dispersion, only: dispersion_available, dispersion_correction
   use mqc_dispersion_d4, only: dispersion_d4_available, dispersion_d4_correction
   implicit none
   private

   public :: dispersion_kind_available
   public :: dispersion_kind_option
   public :: dispersion_apply

contains

   pure function dispersion_kind_available(kind) result(available)
      !! Whether this build linked the library that serves `kind`
      !!
      !! Spelling is `dispersion_kind_is_known`'s question and is not asked
      !! again here: an unknown name is served by nothing, so it answers false
      !! and the caller's message about spelling is the one that should be seen.
      character(len=*), intent(in) :: kind
      logical :: available

      if (dispersion_kind_is_d3(kind)) then
         available = dispersion_available()
      else if (dispersion_kind_is_d4(kind)) then
         available = dispersion_d4_available()
      else
         available = .false.
      end if
   end function dispersion_kind_available

   pure function dispersion_kind_option(kind) result(option)
      !! The CMake option that would have supplied `kind`
      !!
      !! So that a refusal names the flag that fixes it, and the right one: the
      !! two options are independent, and telling someone who asked for "d4" to
      !! turn D3 on would send them round the loop twice.
      character(len=*), intent(in) :: kind
      character(len=16) :: option

      if (dispersion_kind_is_d4(kind)) then
         option = "MQC_ENABLE_DFTD4"
      else
         option = "MQC_ENABLE_DFTD3"
      end if
   end function dispersion_kind_option

   subroutine dispersion_apply(kind, functional, charge, atomic_numbers, coordinates, &
                               energy, gradient, error)
      !! The dispersion energy for `kind`, and its gradient when one is asked for
      !!
      !! `coordinates` is (3, natoms) in Bohr and `gradient` comes back in the
      !! same shape in Hartree/Bohr. `charge` is the total molecular charge and
      !! is required of every caller even though only D4 reads it -- an
      !! optional argument here would make "this correction does not need a
      !! charge" and "nobody passed one" the same call.
      character(len=*), intent(in) :: kind
      character(len=*), intent(in) :: functional
      real(dp), intent(in) :: charge
      integer(default_int), intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), intent(out) :: energy
      real(dp), intent(out), optional :: gradient(:, :)
      type(error_t), intent(out) :: error

      energy = 0.0_dp
      if (present(gradient)) gradient = 0.0_dp

      if (dispersion_kind_is_d3(kind)) then
         ! D3 has no charge dependence at all -- s-dftd3's structure
         ! constructor takes none -- so `charge` is dropped here rather than
         ! carried into a wrapper that would have to ignore it.
         if (present(gradient)) then
            call dispersion_correction(kind, functional, atomic_numbers, coordinates, &
                                       energy, gradient, error)
         else
            call dispersion_correction(kind, functional, atomic_numbers, coordinates, &
                                       energy, error=error)
         end if
      else if (dispersion_kind_is_d4(kind)) then
         if (present(gradient)) then
            call dispersion_d4_correction(kind, functional, charge, atomic_numbers, coordinates, &
                                          energy, gradient, error)
         else
            call dispersion_d4_correction(kind, functional, charge, atomic_numbers, coordinates, &
                                          energy, error=error)
         end if
      else
         call error%set(ERROR_VALIDATION, "unknown dispersion correction '"// &
                        trim(adjustl(kind))//"'. Known: "//DISPERSION_KINDS//".")
      end if
   end subroutine dispersion_apply

end module mqc_dispersion_apply
