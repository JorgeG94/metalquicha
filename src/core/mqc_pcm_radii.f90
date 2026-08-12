!! Cavity radii for a polarizable continuum
module mqc_pcm_radii
   !! Van der Waals radii, and the scaling that turns them into a solute cavity.
   !!
   !! **This is our data, which is the whole risk.** A continuum model's answer
   !! depends on the cavity more strongly than on almost anything else in it -- the
   !! solvation energy goes roughly as 1/R -- so a radius that is wrong by a few
   !! percent produces a solvation energy that is wrong by a few percent while
   !! looking entirely reasonable. There is no internal identity that catches it,
   !! the way the electron count catches a wrong density. So: one table, in one
   !! place, with the paper it came from written beside it, and a test that pins
   !! the values a reader is most likely to recognise.
   !!
   !! The values are Bondi's [J. Phys. Chem. 68, 441 (1964)], filled in for the
   !! elements Bondi did not cover from Mantina et al. [J. Phys. Chem. A 113, 5806
   !! (2009)], which is the same combination PySCF, ORCA and Q-Chem use for their
   !! default PCM cavities. They are *van der Waals* radii and not the
   !! Bragg-Slater ones in `mqc_dft_radial`: those size an integration grid, are
   !! systematically smaller, and using them here would shrink every cavity.
   !!
   !! Radii are returned in **Bohr**, like every other length in this program,
   !! though the table is written in Angstrom because that is how both papers
   !! tabulate them and a converted table cannot be checked against its source.
   use pic_types, only: dp
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: vdw_radius
   public :: cavity_radius
   public :: DEFAULT_RADII_SCALE
   public :: MAX_PCM_ELEMENT

   !> Highest atomic number with a radius here.
   !>
   !> Argon, matching the range the CPU validation suite covers. A heavier atom
   !> is refused rather than given a guessed radius: an invented cavity radius is
   !> exactly the failure this module is arranged to prevent.
   integer, parameter :: MAX_PCM_ELEMENT = 18

   !> Scaling from a van der Waals radius to a cavity radius.
   !>
   !> 1.2 is the long-standing convention, from the observation that the solute
   !> cavity that reproduces experimental solvation energies is somewhat larger
   !> than the bare van der Waals surface. It is a fitted number, not a derived
   !> one, and codes differ: Gaussian's UFF cavity and ORCA's default differ from
   !> this and from each other. So it is a default rather than a constant, and a
   !> deck can override it.
   real(dp), parameter :: DEFAULT_RADII_SCALE = 1.2_dp

   !> Van der Waals radii in Angstrom, H through Ar, indexed by atomic number.
   !>
   !> Bondi (1964) except Li, Be, B, Na, Mg and Al, which are Mantina (2009).
   real(dp), parameter :: VDW_ANGSTROM(MAX_PCM_ELEMENT) = [ &
                          1.20_dp, &   ! H
                          1.40_dp, &   ! He
                          1.82_dp, &   ! Li
                          1.53_dp, &   ! Be   (Mantina)
                          1.92_dp, &   ! B    (Mantina)
                          1.70_dp, &   ! C
                          1.55_dp, &   ! N
                          1.52_dp, &   ! O
                          1.47_dp, &   ! F
                          1.54_dp, &   ! Ne
                          2.27_dp, &   ! Na   (Mantina)
                          1.73_dp, &   ! Mg   (Mantina)
                          1.84_dp, &   ! Al   (Mantina)
                          2.10_dp, &   ! Si
                          1.80_dp, &   ! P
                          1.80_dp, &   ! S
                          1.75_dp, &   ! Cl
                          1.88_dp]     ! Ar

contains

   subroutine vdw_radius(atomic_number, radius, error)
      !! The van der Waals radius of an element, in Bohr
      !!
      !! Refused rather than extrapolated outside the table. A continuum model
      !! given a plausible-looking wrong radius returns a plausible-looking wrong
      !! solvation energy, and nothing downstream would question it.
      integer, intent(in) :: atomic_number
      real(dp), intent(out) :: radius
      type(error_t), intent(inout) :: error

      character(len=8) :: number_text

      radius = 0.0_dp
      if (atomic_number < 1 .or. atomic_number > MAX_PCM_ELEMENT) then
         write (number_text, "(i0)") atomic_number
         call error%set(ERROR_VALIDATION, "no cavity radius for element "// &
                        trim(number_text)//": the continuum model here carries van "// &
                        "der Waals radii for hydrogen through argon only. Refused "// &
                        "rather than guessed, since a wrong cavity radius gives a "// &
                        "wrong solvation energy and nothing would say so.")
         return
      end if

      radius = ANGSTROM_TO_BOHR*VDW_ANGSTROM(atomic_number)
   end subroutine vdw_radius

   subroutine cavity_radius(atomic_number, scale, radius, error)
      !! A scaled cavity radius, in Bohr
      integer, intent(in) :: atomic_number
      real(dp), intent(in) :: scale        !! Typically `DEFAULT_RADII_SCALE`
      real(dp), intent(out) :: radius
      type(error_t), intent(inout) :: error

      call vdw_radius(atomic_number, radius, error)
      if (error%has_error()) return
      if (scale <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "the cavity radius scale must be positive")
         return
      end if
      radius = scale*radius
   end subroutine cavity_radius

end module mqc_pcm_radii
