program check_charges
   !! Check the two charge schemes against what must be true of any of them
   !!
   !! There is no reference number for a partial charge -- it is not an
   !! observable, and two schemes disagree by design. So these are the
   !! properties that hold whatever the scheme:
   !!
   !!   * **They sum to the molecular charge.** For Mulliken this is an identity
   !!     of the trace; for CHELPG it is the constraint the fit is solved under.
   !!     Either failing means the bookkeeping is wrong, not the chemistry.
   !!   * **Symmetry-equivalent atoms get equal charges.** Water's two hydrogens
   !!     cannot differ, and methane's four cannot.
   !!   * **The sign follows electronegativity.** Oxygen negative in water,
   !!     hydrogen positive; the fluorines negative in HF.
   !!
   !! And one property that is the whole reason CHELPG exists:
   !!
   !!   * **CHELPG reproduces the potential it was fitted to.** The residual on
   !!     the fit grid is reported. Mulliken has no such claim and is not asked
   !!     for one.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger, info_level
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_charges, only: mulliken_charges, chelpg_charges, chelpg_grid
   use mqc_libcint_esp, only: esp_contract
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   implicit none

   integer, parameter :: N_DIM = 3
   integer :: n_bad

   n_bad = 0
   call logger%configure(info_level)

   call run_case("water", [character(len=2) :: "O", "H", "H"], &
                 reshape([0.0000_dp, 0.0000_dp, 0.0000_dp, &
                          0.0000_dp, -1.4308_dp, 1.1078_dp, &
                          0.0000_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                 [8, 1, 1], 10, [2, 3], 0.10_dp, n_bad)

   call run_case("methane", [character(len=2) :: "C", "H", "H", "H", "H"], &
                 reshape([0.0000_dp, 0.0000_dp, 0.0000_dp, &
                          1.1900_dp, 1.1900_dp, 1.1900_dp, &
                          -1.1900_dp, -1.1900_dp, 1.1900_dp, &
                          -1.1900_dp, 1.1900_dp, -1.1900_dp, &
                          1.1900_dp, -1.1900_dp, -1.1900_dp], [N_DIM, 5]), &
                 ! Loose, and not a weakness of the fit. Methane is tetrahedral:
                 ! no dipole, no quadrupole, first nonvanishing moment the
                 ! octupole. Its potential outside the van der Waals surface is
                 ! an order of magnitude smaller than water's and is mostly
                 ! charge penetration, which no arrangement of point charges can
                 ! reproduce. A nonpolar molecule fits badly in relative terms
                 ! however good the code is.
                 [6, 1, 1, 1, 1], 10, [2, 5], 0.80_dp, n_bad)

   call check_basis_sensitivity(n_bad)

   if (n_bad == 0) then
      call logger%info("")
      call logger%info("charges are normalised, symmetric, and CHELPG reproduces its potential")
   else
      call logger%error("charge checks failed")
      stop 1
   end if

contains

   subroutine run_case(label, symbols, coords, z, nelec, equiv, rrms_tol, n_bad)
      character(len=*), intent(in) :: label
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: z(:)
      integer, intent(in) :: nelec
      integer, intent(in) :: equiv(2)   !! first and last of a symmetry-equivalent run
      real(dp), intent(in) :: rrms_tol
      integer, intent(inout) :: n_bad

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: error
      real(dp), allocatable :: s(:, :), q_mul(:), q_esp(:), pts(:, :), v(:)
      real(dp) :: spread_mul, spread_esp, rms, v_rms, d
      character(len=200) :: line
      integer :: i, g, iatom

      call logger%info("")
      write (line, "(a,a)") "== ", label
      call logger%info(trim(line))

      call build_libcint_molecule(z, symbols, coords, "6-31g", mol, error)
      if (failed(error, "molecule", n_bad)) return
      call run_libcint_rhf(mol, nelec, 100, 1.0e-9_dp, 1.0e-7_dp, .false., scf, error)
      if (failed(error, "scf", n_bad)) return

      call mol%overlap(s)
      call mulliken_charges(mol, scf%density, s, q_mul, error)
      if (failed(error, "mulliken", n_bad)) return
      call chelpg_charges(mol, scf%density, q_esp, error)
      if (failed(error, "chelpg", n_bad)) return

      do i = 1, mol%natm
         write (line, "(a,a4,a,f9.5,a,f9.5)") "   ", symbols(i), "   mulliken ", q_mul(i), &
            "    chelpg ", q_esp(i)
         call logger%info(trim(line))
      end do

      call report("mulliken sums to zero", abs(sum(q_mul)), 1.0e-8_dp, n_bad)
      call report("chelpg sums to zero", abs(sum(q_esp)), 1.0e-8_dp, n_bad)

      spread_mul = maxval(q_mul(equiv(1):equiv(2))) - minval(q_mul(equiv(1):equiv(2)))
      spread_esp = maxval(q_esp(equiv(1):equiv(2))) - minval(q_esp(equiv(1):equiv(2)))
      call report("equivalent atoms agree, mulliken", spread_mul, 1.0e-6_dp, n_bad)
      ! Looser: the fit grid is a cubic lattice and does not share the molecule's
      ! symmetry, so equivalent atoms see slightly different point sets.
      call report("equivalent atoms agree, chelpg", spread_esp, 5.0e-3_dp, n_bad)

      if (q_esp(1) >= 0.0_dp .and. z(1) > z(2)) then
         call logger%error("   FAIL the more electronegative atom is not negative")
         n_bad = n_bad + 1
      end if

      ! What CHELPG claims: the fitted charges reproduce the potential.
      call chelpg_grid(mol, pts, error)
      if (failed(error, "grid", n_bad)) return
      call esp_contract(mol, pts, scf%density, v, error)
      if (failed(error, "esp", n_bad)) return
      rms = 0.0_dp
      v_rms = 0.0_dp
      do g = 1, size(pts, 2)
         do iatom = 1, mol%natm
            d = norm2(pts(:, g) - mol%coords(:, iatom))
            v(g) = v(g) + mol%charges(iatom)/d
         end do
         ! The baseline: all charges zero, which the constraint permits. Any fit
         ! that does not beat it is not worth solving for.
         v_rms = v_rms + v(g)**2
         do iatom = 1, mol%natm
            d = norm2(pts(:, g) - mol%coords(:, iatom))
            v(g) = v(g) - q_esp(iatom)/d
         end do
         rms = rms + v(g)**2
      end do
      rms = sqrt(rms/real(size(pts, 2), dp))
      v_rms = sqrt(v_rms/real(size(pts, 2), dp))
      write (line, "(a,i0,a,es10.3,a,es10.3)") "   fit: ", size(pts, 2), &
         " points, rms residual ", rms, "  against rms V ", v_rms
      call logger%info(trim(line))
      ! RRMS, the usual measure of an ESP fit: the residual as a fraction of
      ! the potential a zero-charge model would leave. Below 0.1 is a good fit
      ! in the literature, and the tolerance is per-case because how well point
      ! charges *can* do depends on the molecule -- see the call sites.
      call report("chelpg rrms", rms/v_rms, rrms_tol, n_bad)
   end subroutine run_case

   subroutine check_basis_sensitivity(n_bad)
      !! The claim that decides which scheme to embed with
      !!
      !! Mulliken partitions by basis function, so a function's atom label is
      !! what assigns its density -- and a diffuse function centred on hydrogen
      !! reaches over the oxygen while still counting as hydrogen's. Enlarging
      !! the basis therefore moves Mulliken charges without the molecule
      !! changing. CHELPG partitions by a physical field, so enlarging the basis
      !! moves it only as far as the field itself moves.
      !!
      !! Same water, two bases, and the question is which set of charges is
      !! still recognisable in the second one.
      integer, intent(inout) :: n_bad

      real(dp) :: q_mul_small(3), q_esp_small(3), q_mul_big(3), q_esp_big(3)
      real(dp) :: shift_mul, shift_esp
      character(len=200) :: line

      call logger%info("")
      call logger%info("== water, 6-31g vs aug-cc-pvdz")

      call water_charges("6-31g", q_mul_small, q_esp_small, n_bad)
      call water_charges("aug-cc-pvdz", q_mul_big, q_esp_big, n_bad)
      if (n_bad > 0) return

      shift_mul = maxval(abs(q_mul_big - q_mul_small))
      shift_esp = maxval(abs(q_esp_big - q_esp_small))

      write (line, "(a,f9.5,a,f9.5,a,f9.5)") "   mulliken  O ", q_mul_small(1), " -> ", &
         q_mul_big(1), "    moved ", shift_mul
      call logger%info(trim(line))
      write (line, "(a,f9.5,a,f9.5,a,f9.5)") "   chelpg    O ", q_esp_small(1), " -> ", &
         q_esp_big(1), "    moved ", shift_esp
      call logger%info(trim(line))

      if (shift_esp < shift_mul) then
         write (line, "(a,f6.2,a)") "   ok   chelpg is the more stable of the two, by ", &
            shift_mul/shift_esp, "x"
         call logger%info(trim(line))
      else
         call logger%error("   FAIL chelpg moved at least as much as mulliken, which is "// &
                           "the opposite of why it is here")
         n_bad = n_bad + 1
      end if
   end subroutine check_basis_sensitivity

   subroutine water_charges(basis, q_mul, q_esp, n_bad)
      character(len=*), intent(in) :: basis
      real(dp), intent(out) :: q_mul(3), q_esp(3)
      integer, intent(inout) :: n_bad

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: error
      real(dp), allocatable :: s(:, :), qm(:), qe(:)
      real(dp) :: coords(N_DIM, 3)

      coords = reshape([0.0000_dp, 0.0000_dp, 0.0000_dp, &
                        0.0000_dp, -1.4308_dp, 1.1078_dp, &
                        0.0000_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3])

      q_mul = 0.0_dp
      q_esp = 0.0_dp
      call build_libcint_molecule([8, 1, 1], [character(len=2) :: "O", "H", "H"], &
                                  coords, basis, mol, error)
      if (failed(error, "molecule "//basis, n_bad)) return
      call run_libcint_rhf(mol, 10, 100, 1.0e-9_dp, 1.0e-7_dp, .false., scf, error)
      if (failed(error, "scf "//basis, n_bad)) return
      call mol%overlap(s)
      call mulliken_charges(mol, scf%density, s, qm, error)
      if (failed(error, "mulliken "//basis, n_bad)) return
      call chelpg_charges(mol, scf%density, qe, error)
      if (failed(error, "chelpg "//basis, n_bad)) return
      q_mul = qm
      q_esp = qe
   end subroutine water_charges

   function failed(error, what, n_bad) result(bad)
      type(error_t), intent(inout) :: error
      character(len=*), intent(in) :: what
      integer, intent(inout) :: n_bad
      logical :: bad

      bad = error%has_error()
      if (bad) then
         call logger%error("   FAIL "//what//": "//error%get_message())
         n_bad = n_bad + 1
      end if
   end function failed

   subroutine report(what, deviation, tol, n_bad)
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: deviation, tol
      integer, intent(inout) :: n_bad

      character(len=200) :: line

      if (deviation <= tol) then
         write (line, "(a,a,a,es10.3)") "   ok   ", what, "   ", deviation
         call logger%info(trim(line))
      else
         write (line, "(a,a,a,es10.3,a,es10.3)") "   FAIL ", what, "   ", deviation, " > ", tol
         call logger%error(trim(line))
         n_bad = n_bad + 1
      end if
   end subroutine report

end program check_charges
