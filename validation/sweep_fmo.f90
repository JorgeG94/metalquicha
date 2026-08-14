program sweep_fmo
   !! How should a fragment's neighbours be represented to it?
   !!
   !! A measurement, not a test: it prints a table and asserts nothing, because
   !! what it measures is a property of the method and would be a brittle thing
   !! to pin an inequality on. `check_fmo` holds the assertions -- including the
   !! long-separation one, which is the sharpest conclusion drawn here.
   !!
   !! Stacked water clusters of 3, 4 and 5 monomers at six separations, each run
   !! four ways against a supermolecular RHF in the same basis. What it found:
   !!
   !!   * **The exact ESP is exact where it can be.** At 6 to 9 Angstrom, where
   !!     exchange and charge transfer between fragments have died but the
   !!     electrostatic field has not, FMO2 with the density-based embedding
   !!     reproduces the supermolecule to 1e-13. Only a correct embedding
   !!     operator does that, which makes those rows the real verification of
   !!     the method and not merely a comparison between variants.
   !!   * **Point charges plateau.** ESP-PTC stops improving at about 1e-08
   !!     (Mulliken) or 1e-10 (CHELPG) no matter how far apart the fragments are
   !!     moved. That floor is the charge approximation itself, and it is what
   !!     separates electrostatically embedded MBE from FMO.
   !!   * **Do not read the short-separation rows as favouring point charges.**
   !!     The exact ESP's residual near contact is the genuine three-body term,
   !!     which FMO2 does not contain; it is negative throughout and shrinks
   !!     monotonically as that term dies. The point-charge error changes sign
   !!     across the range, so near contact it cancels against the three-body
   !!     term and looks better than it is. Two unrelated errors happening to
   !!     have opposite signs is not accuracy, and nothing keeps them matched
   !!     for a different system.
   !!   * Embedding of either kind beats none by one to three orders of
   !!     magnitude everywhere, which is the argument for FMO over an
   !!     unembedded expansion.
   !!
   !! Rerun after touching the embedding, the charges, or the ESP integrals: the
   !! pattern above is what would break.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger, warning_level
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   implicit none

   real(dp), parameter :: A2B = 1.8897261254578281_dp
   integer, parameter :: N_SEP = 6
   integer, parameter :: N_SIZE = 3
   real(dp) :: sep(N_SEP) = [2.7_dp, 2.9_dp, 3.2_dp, 4.0_dp, 6.0_dp, 9.0_dp]
   integer :: nw(N_SIZE) = [3, 4, 5]
   integer :: i, k, wins_exact, total
   real(dp) :: e_ex, e_exact, e_mul, e_che, e_none

   call logger%configure(warning_level)
   wins_exact = 0
   total = 0

   print "(a)", "  waters  sep(A)   exact-ESP err   ptc-mulliken   ptc-chelpg      no embed   closest"
   do k = 1, size(nw)
      do i = 1, size(sep)
         call one(nw(k), sep(i), e_ex, e_exact, e_mul, e_che, e_none)
         if (e_ex == 0.0_dp) cycle
         total = total + 1
         if (abs(e_exact - e_ex) <= min(abs(e_mul - e_ex), abs(e_che - e_ex))) then
            wins_exact = wins_exact + 1
         end if
         print "(i6,f9.2,4es15.2,a)", nw(k), sep(i), &
            e_exact - e_ex, e_mul - e_ex, e_che - e_ex, e_none - e_ex, &
            closest(abs(e_exact - e_ex), abs(e_mul - e_ex), abs(e_che - e_ex))
      end do
   end do
   print "(a,i0,a,i0)", "exact ESP closest in ", wins_exact, " of ", total

contains

   function closest(a, b, c) result(name)
      real(dp), intent(in) :: a, b, c
      character(len=13) :: name

      name = "   ptc-chelpg"
      if (b <= a .and. b <= c) name = " ptc-mulliken"
      if (a <= b .and. a <= c) name = "    exact-ESP"
   end function closest

   subroutine one(n, d, e_reference, e_exact, e_mul, e_che, e_none)
      integer, intent(in) :: n
      real(dp), intent(in) :: d
      real(dp), intent(out) :: e_reference, e_exact, e_mul, e_che, e_none

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: r
      type(error_t) :: error
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      integer, allocatable :: z(:), owner(:)
      character(len=2), allocatable :: sym(:)
      real(dp), allocatable :: xyz(:, :)
      real(dp) :: mono(3, 3)
      integer :: w, j, at

      e_reference = 0.0_dp
      e_exact = 0.0_dp
      e_mul = 0.0_dp
      e_che = 0.0_dp
      e_none = 0.0_dp
      mono = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                      0.0_dp, -0.7572_dp, 0.5865_dp, &
                      0.0_dp, 0.7572_dp, 0.5865_dp], [3, 3])
      allocate (z(3*n), sym(3*n), owner(3*n), xyz(3, 3*n))
      at = 0
      do w = 1, n
         do j = 1, 3
            at = at + 1
            z(at) = merge(8, 1, j == 1)
            sym(at) = merge("O ", "H ", j == 1)
            xyz(:, at) = mono(:, j)
            xyz(3, at) = xyz(3, at) + real(w - 1, dp)*d
            owner(at) = w
         end do
      end do
      xyz = xyz*A2B

      call build_libcint_molecule(z, sym, xyz, "6-31g", mol, error)
      if (error%has_error()) return
      call run_libcint_rhf(mol, sum(z), 300, 1.0e-10_dp, 1.0e-8_dp, .false., scf, error)
      if (error%has_error() .or. .not. scf%converged) return
      e_reference = scf%energy

      opts%basis = "6-31g"
      opts%esp = "exact"
      call run_fmo2(z, sym, xyz, owner, opts, r, error)
      if (error%has_error()) then
         call error%clear()
      else
         e_exact = r%energy
      end if

      opts%esp = "ptc"
      opts%ptc_charges = "mulliken"
      call run_fmo2(z, sym, xyz, owner, opts, r, error)
      if (error%has_error()) then
         call error%clear()
      else
         e_mul = r%energy
      end if

      opts%esp = "ptc"
      opts%ptc_charges = "chelpg"
      call run_fmo2(z, sym, xyz, owner, opts, r, error)
      if (error%has_error()) then
         call error%clear()
      else
         e_che = r%energy
      end if

      opts%esp = "none"
      call run_fmo2(z, sym, xyz, owner, opts, r, error)
      if (error%has_error()) then
         call error%clear()
      else
         e_none = r%energy
      end if
   end subroutine one

end program sweep_fmo
