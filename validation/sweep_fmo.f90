program sweep_fmo
   !! Which point charges should represent a fragment to its neighbours?
   !!
   !! A measurement, not a test -- it prints a table and asserts nothing, because
   !! what it is measuring is a property of the *method* and would be a brittle
   !! thing to pin an inequality on. `check_fmo` holds the assertions.
   !!
   !! Stacked water clusters of 3, 4 and 5 monomers at four separations, each
   !! run three ways and compared against a supermolecular RHF in the same
   !! basis. What it found, and what `fmo_options_t%embedding` documents:
   !!
   !!   * Embedding of either kind beats none by one to two orders of
   !!     magnitude, in every case. That is the argument for FMO over an
   !!     unembedded many-body expansion, and it is not close.
   !!   * CHELPG is closer than Mulliken in 9 of 12 cases, and all three
   !!     exceptions are at the tightest separation, 2.70 A. That is where a
   !!     neighbouring fragment sits inside the van der Waals surface CHELPG
   !!     excludes from its fit -- so it is being represented by charges that
   !!     carry no information about where it actually is. Mulliken has no
   !!     excluded region and degrades more gently.
   !!   * Past 2.9 A, where hydrogen-bonded monomers actually sit, CHELPG wins
   !!     by roughly an order of magnitude, and the margin grows with both
   !!     separation and cluster size.
   !!
   !! Rerun it after touching the embedding, the charges, or the ESP integrals:
   !! the pattern above is the thing that would break.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger, warning_level
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   implicit none

   real(dp), parameter :: A2B = 1.8897261254578281_dp
   real(dp) :: sep(4) = [2.7_dp, 2.9_dp, 3.2_dp, 4.0_dp]
   integer :: nw(3) = [3, 4, 5]
   integer :: i, k, wins_chelpg, total
   real(dp) :: e_ex, e_mul, e_che, e_none

   call logger%configure(warning_level)
   wins_chelpg = 0
   total = 0

   print "(a)", "  waters   sep(A)      exact          mulliken err    chelpg err     none err    better"
   do k = 1, size(nw)
      do i = 1, size(sep)
         call one(nw(k), sep(i), e_ex, e_mul, e_che, e_none)
         if (e_ex == 0.0_dp) cycle
         total = total + 1
         if (abs(e_che - e_ex) < abs(e_mul - e_ex)) wins_chelpg = wins_chelpg + 1
         print "(i6,f10.2,f16.8,3es15.2,a)", nw(k), sep(i), e_ex, &
            e_mul - e_ex, e_che - e_ex, e_none - e_ex, &
            merge("  chelpg", "mulliken", abs(e_che - e_ex) < abs(e_mul - e_ex))
      end do
   end do
   print "(a,i0,a,i0)", "chelpg closer in ", wins_chelpg, " of ", total

contains

   subroutine one(n, d, e_exact, e_mul, e_che, e_none)
      integer, intent(in) :: n
      real(dp), intent(in) :: d
      real(dp), intent(out) :: e_exact, e_mul, e_che, e_none

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
      e_exact = scf%energy

      opts%basis = "6-31g"
      opts%embedding = "mulliken"
      call run_fmo2(z, sym, xyz, owner, opts, r, error)
      if (error%has_error()) then
         call error%clear()
      else
         e_mul = r%energy
      end if

      opts%embedding = "chelpg"
      call run_fmo2(z, sym, xyz, owner, opts, r, error)
      if (error%has_error()) then
         call error%clear()
      else
         e_che = r%energy
      end if

      opts%embedding = "none"
      call run_fmo2(z, sym, xyz, owner, opts, r, error)
      if (error%has_error()) then
         call error%clear()
      else
         e_none = r%energy
      end if
   end subroutine one

end program sweep_fmo
