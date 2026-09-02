!! The EFP2 interaction energy, all five terms
module mqc_efp_energy
   !! What the rest of the EFP code is for: given placed fragments, one number.
   !!
   !! The five terms live in two modules because they need different things.
   !! `mqc_efp_interaction` works on a flattened point set and covers the terms that
   !! are geometry and stored tensors -- electrostatics, polarization, and the
   !! undamped `E6`. `mqc_efp_pair` works on fragment *pairs* and covers everything
   !! needing integrals over two fragments' basis sets at once: exchange repulsion,
   !! charge transfer, and the damped dispersion, whose damping is an overlap between
   !! localized orbitals on different fragments. This module is where the two meet,
   !! and it exists so no caller has to know which term came from where.
   !!
   !! **Every term here is separately validated against GAMESS.** The references are
   !! in `EFP_PLAN.md`, all for the same water dimer, and the point of keeping the
   !! breakdown in the result rather than returning a bare total is that a regression
   !! then names the term it broke.
   !!
   !! **Dispersion is the damped sum `E6 + E7 + E8`**, which is what GAMESS totals
   !! (`efdrvr.src:1923`) when the potential carries the tensors for all three. The
   !! undamped `dispersion_energy_e6` in `mqc_efp_interaction` is not used here; it
   !! exists because it was the rung the ladder started on.
   !!
   !! **Fragments arrive already turned.** `place_fragment` gives the rotation a deck
   !! implies and `mqc_efp_rotate` applies it, so by the time a fragment reaches this
   !! module its own stored frame *is* the working frame and all that is left is the
   !! translation each term takes as an offset.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_efp_read, only: efp_fragment_t
   use mqc_efp_interaction, only: efp_system_t, build_efp_system, &
                                  electrostatic_energy, polarization_energy
   use mqc_efp_pair, only: exchange_repulsion, charge_transfer, &
                           dispersion_e6_damped, dispersion_e7_damped, &
                           dispersion_e8_damped
   use mqc_efp_rotate, only: superpose
   implicit none
   private

   public :: efp_energy_t
   public :: efp_interaction_energy
   public :: place_fragment

   !! How far a deck atom may sit from where the potential's own geometry puts it,
   !! after the rigid shift, before the placement is refused. A fragment is rigid, so
   !! this is not a fitting tolerance -- it is the width of the round trip through
   !! the file, which carries ten decimals.
   real(dp), parameter :: PLACEMENT_TOL = 1.0e-6_dp

   !! Every multipole rank the electrostatics carries: charges, dipoles,
   !! quadrupoles and octupoles.
   integer, parameter :: MAX_RANK = 3

   type :: efp_energy_t
      !! One interaction energy, kept broken down by term
      real(dp) :: electrostatics = 0.0_dp
      real(dp) :: polarization = 0.0_dp
      real(dp) :: exchange_repulsion = 0.0_dp
      real(dp) :: dispersion = 0.0_dp     !! Damped `E6 + E7 + E8`
      real(dp) :: dispersion_e6 = 0.0_dp
      real(dp) :: dispersion_e7 = 0.0_dp
      real(dp) :: dispersion_e8 = 0.0_dp
      real(dp) :: charge_transfer = 0.0_dp
      real(dp) :: total = 0.0_dp
   end type efp_energy_t

contains

   function efp_interaction_energy(fragments, translations, error) result(energy)
      !! The interaction energy of a set of placed fragments
      !!
      !! `translations(:, k)` places `fragments(k)`, in Bohr, as a rigid shift of the
      !! geometry the potential itself carries.
      !!
      !! Electrostatics and polarization are computed over the whole system at once --
      !! polarization *must* be, since the induced dipoles are solved together and are
      !! not a sum of pair terms. The other three are pairwise and are summed over
      !! `a < b`. Nothing here assumes two fragments.
      !!
      !! A term whose data the potential does not carry is left at zero rather than
      !! erroring: a potential written without the dynamic polarizabilities has no
      !! dispersion to contribute, and refusing the whole calculation over it would
      !! make a legitimate potential unusable. What is *not* tolerated is a term
      !! failing for any other reason, which propagates.
      type(efp_fragment_t), intent(in) :: fragments(:)
      real(dp), intent(in) :: translations(:, :)   !! (3, n_fragments), Bohr
      type(error_t), intent(inout) :: error
      type(efp_energy_t) :: energy

      type(efp_system_t) :: system
      integer :: n, a, b
      logical :: have_dynamic, have_lmo, have_ct

      n = size(fragments)
      if (size(translations, 1) /= 3 .or. size(translations, 2) /= n) then
         call error%set(ERROR_VALIDATION, "efp: one translation per fragment is "// &
                        "needed, as (3, n_fragments)")
         return
      end if
      if (n < 2) return      ! a single fragment interacts with nothing

      call build_efp_system(fragments, translations, system, error)
      if (error%has_error()) return

      energy%electrostatics = electrostatic_energy(system, MAX_RANK, screen=.true.)
      energy%polarization = polarization_energy(system, fragments, error)
      if (error%has_error()) return

      do a = 1, n
         do b = a + 1, n
            have_dynamic = fragments(a)%has_dynamic .and. fragments(b)%has_dynamic
            have_lmo = fragments(a)%has_lmo .and. fragments(b)%has_lmo
            have_ct = fragments(a)%has_ctvec .and. fragments(b)%has_ctvec &
                      .and. fragments(a)%has_ctfok .and. fragments(b)%has_ctfok

            if (have_lmo) then
               energy%exchange_repulsion = energy%exchange_repulsion &
                                           + exchange_repulsion(fragments(a), fragments(b), &
                                                                translations(:, a), translations(:, b), error)
               if (error%has_error()) return
            end if

            if (have_dynamic .and. have_lmo) then
               energy%dispersion_e6 = energy%dispersion_e6 &
                                      + dispersion_e6_damped(fragments(a), fragments(b), &
                                                             translations(:, a), translations(:, b), error)
               if (error%has_error()) return
               ! E7 and E8 need the higher tensor blocks, which a potential may omit
               ! -- and GAMESS treats their absence as switching the term off rather
               ! than as an error (`efinp.src:5559, 5577`), so this does too.
               if (fragments(a)%has_dipquad .and. fragments(b)%has_dipquad) then
                  energy%dispersion_e7 = energy%dispersion_e7 &
                                         + dispersion_e7_damped(fragments(a), fragments(b), &
                                                                translations(:, a), translations(:, b), error)
                  if (error%has_error()) return
               end if
               if (fragments(a)%has_quadquad .and. fragments(b)%has_quadquad) then
                  energy%dispersion_e8 = energy%dispersion_e8 &
                                         + dispersion_e8_damped(fragments(a), fragments(b), &
                                                                translations(:, a), translations(:, b), error)
                  if (error%has_error()) return
               end if
            end if

            if (have_ct) then
               energy%charge_transfer = energy%charge_transfer &
                                        + charge_transfer(fragments(a), fragments(b), &
                                                          translations(:, a), translations(:, b), error)
               if (error%has_error()) return
            end if
         end do
      end do

      energy%dispersion = energy%dispersion_e6 + energy%dispersion_e7 &
                          + energy%dispersion_e8
      energy%total = energy%electrostatics + energy%polarization &
                     + energy%exchange_repulsion + energy%dispersion &
                     + energy%charge_transfer

      call system%destroy()
   end function efp_interaction_energy

   subroutine place_fragment(frag, coords, rot, translation, error)
      !! Where a deck's atoms put a fragment: a rotation and a shift
      !!
      !! A potential carries the geometry it was made at; a deck names atoms in its
      !! own coordinates. Placing the fragment is finding the rigid transform between
      !! the two, `deck = rot . own + translation`, which `superpose` builds from the
      !! frame three atoms define.
      !!
      !! **The rotation is returned, not applied.** Turning the fragment is
      !! `rotate_fragment`, and it is a separate step because it rewrites every tensor
      !! the potential carries and the caller may want to know the transform without
      !! paying for that.
      !!
      !! Every atom is then required to land where the potential says, which is a test
      !! rather than a fit: an EFP fragment is rigid, so a residual means the deck and
      !! the potential disagree about something. It catches a potential paired with the
      !! wrong species, atoms listed in a different order, and a geometry that has been
      !! relaxed since the potential was made -- each of which would otherwise place
      !! the fragment somewhere plausible and return an energy.
      type(efp_fragment_t), intent(in) :: frag
      real(dp), intent(in) :: coords(:, :)     !! (3, n_atoms) from the deck, Bohr
      real(dp), intent(out) :: rot(3, 3)
      real(dp), intent(out) :: translation(3)
      type(error_t), intent(inout) :: error

      real(dp) :: rmsd
      integer :: n

      n = size(coords, 2)
      if (size(coords, 1) /= 3) then
         call error%set(ERROR_VALIDATION, "efp: fragment coordinates must be (3, n)")
         return
      end if
      if (n /= frag%n_atoms) then
         call error%set(ERROR_VALIDATION, "efp: this fragment has "//count_text(n)// &
                        " atoms in the deck but its potential describes "// &
                        count_text(frag%n_atoms))
         return
      end if

      call superpose(frag%points(:, 1:n), coords, rot, translation, rmsd, error)
      if (error%has_error()) return

      if (rmsd > PLACEMENT_TOL) then
         call error%set(ERROR_VALIDATION, "efp: this fragment is not a rigid "// &
                        "placement of its potential's own geometry -- the atoms miss "// &
                        "by "//real_text(rmsd)//" Bohr after the best rigid "// &
                        "transform. Check that the potential describes this species, "// &
                        "that the atoms are listed in the potential's own order, and "// &
                        "that the geometry is the one the potential was built at.")
         return
      end if
   end subroutine place_fragment

   pure function real_text(x) result(text)
      !! A number as text, for a message
      real(dp), intent(in) :: x
      character(len=:), allocatable :: text

      character(len=24) :: buffer

      write (buffer, "(ES10.3)") x
      text = trim(adjustl(buffer))
   end function real_text

   pure function count_text(n) result(text)
      !! A small integer as text, for a message
      integer, intent(in) :: n
      character(len=:), allocatable :: text

      character(len=16) :: buffer

      write (buffer, "(I0)") n
      text = trim(buffer)
   end function count_text

end module mqc_efp_energy
