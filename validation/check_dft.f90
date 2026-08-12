!! Manual check that an exchange-correlation energy comes out right
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON -DMQC_ENABLE_LIBXC=ON
!!     ./build/check_dft && python3 validation/check_dft.py
!!
!! Milestone 1 of the DFT plan: a converged RHF density, its value on the grid, and
!! an LDA exchange-correlation energy from libxc. No self-consistency yet -- the
!! density is Hartree-Fock's and stays fixed, so this isolates the grid, the
!! density assembly and the functional call from anything the SCF might do.
!!
!! Two checks, and only the second needs PySCF:
!!
!!   * the integrated density must be the electron count. Exact, and the one that
!!     catches a factor of two in the restricted density. Also asserted as a unit
!!     test, where it belongs; repeated here because a wrong electron count makes
!!     the energy below meaningless and it is better to say so than to print it.
!!   * E_xc against PySCF's, on the *same grid* and the *same density*. Same-grid
!!     is what makes it sharp: against a different grid it would only say the two
!!     grids are of similar quality.
!!
!! The grid and the density are dumped for the Python side, so the comparison is of
!! the functional evaluation alone and cannot be confounded by either.
program check_dft
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_ao, only: eval_ao_block, eval_rho
   use mqc_libcint_xc, only: xc_context_t, xc_context_create
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_ri_mp2, run_libcint_mp2
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid
   use mqc_error, only: error_t
   use xc_f03_lib_m, only: xc_f03_func_t, xc_f03_func_init, xc_f03_func_end, &
                           xc_f03_lda_exc, XC_UNPOLARIZED
   use xc_f03_funcs_m, only: XC_LDA_X, XC_LDA_C_VWN
   implicit none

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer :: failures

   failures = 0
   call scf_case("cc-pvdz", "lda_x", 3)
   call scf_case("cc-pvdz", "svwn", 3)
   call scf_case("cc-pvdz", "gga_x_pbe", 3)
   call scf_case("cc-pvdz", "pbe", 3)
   call scf_case("cc-pvdz", "b3lyp", 3)
   call scf_case("cc-pvdz", "pbe0", 3)
   call dh_case("cc-pvdz", "cc-pvdz-rifit", "b2plyp", 3)
   call dh_case("cc-pvdz", "cc-pvdz-rifit", "b2gp-plyp", 3)
   call dh_case("cc-pvdz", "cc-pvdz-rifit", "mpw2plyp", 3)
   ! The same functional without fitting the perturbative part, so the
   ! fitting error is measured against our own conventional MP2.
   call dh_case("cc-pvdz", "cc-pvdz-rifit", "b2plyp", 3, fitted=.false.)
   call one_case("cc-pvdz", 3)
   call one_case("cc-pvdz", 5)
   call one_case("sto-3g", 4)

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[dft] wrote every case; run validation/check_dft.py to compare"
   else
      write (*, "(A,I0,A)") "[dft] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine dh_case(basis, auxbasis, functional, level, fitted)
      !! A double hybrid: the Kohn-Sham part, then its perturbative part
      !!
      !! Assembled rather than requested, because libxc carries no double hybrid.
      !! Two steps, both of which already existed:
      !!
      !!   1. a hybrid Kohn-Sham SCF over the functional's semilocal components and
      !!      its exact-exchange fraction -- which `xc_spec_t` already describes and
      !!      `run_libcint_rhf` already runs;
      !!   2. RI-MP2 over the *Kohn-Sham* orbitals, scaled by the perturbative
      !!      fraction. That falls out for free: `run_libcint_ri_mp2` takes
      !!      coefficients and orbital energies and never asks where they came from.
      !!
      !! Density-fitted deliberately. The reference is pyscf-forge's `DFDH`, whose
      !! MP2 part is fitted, so a conventional transform here would sit ~1e-4 away
      !! and look exactly like a wrong coefficient.
      character(len=*), intent(in) :: basis, auxbasis, functional
      integer, intent(in) :: level
      logical, intent(in), optional :: fitted
         !! Fit the perturbative part (default) or transform it exactly. Both are
         !! wanted: the fitted one is comparable to the external reference, and the
         !! pair together measure the fitting error against our own conventional
         !! MP2 -- which is the only check of that error that does not import
         !! somebody else's auxiliary basis along with it.

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2
      type(xc_context_t) :: xc
      type(error_t) :: err
      real(dp) :: c(3, 3)
      real(dp) :: e_total, e_pt2
      integer :: unit
      logical :: use_df
      character(len=4) :: tag

      c = reshape([0.0_dp, 0.00000000009155_dp*ANG, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])

      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[dft] basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, auxbasis, aux, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[dft] aux basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call xc_context_create(mol, functional, xc, err, level=level)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[dft] ", functional, " context failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call run_libcint_rhf(mol, 10, 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err, &
                           in_core=.true., xc=xc)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A,A,A)") "[dft] ", functional, " SCF failed"
         failures = failures + 1
         return
      end if

      ! The perturbative part, over the orbitals the functional just produced.
      ! Frozen core off: DFDH correlates every electron.
      use_df = .true.
      if (present(fitted)) use_df = fitted
      if (use_df) then
         call run_libcint_ri_mp2(mol, aux, scf%orbitals, scf%orbital_energies, 5, &
                                 scf%energy, mp2, err, n_frozen=0)
      else
         call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, 5, &
                              scf%energy, mp2, err, n_frozen=0)
      end if
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[dft] ", functional, " RI-MP2 failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      e_pt2 = xc%pt2_fraction*mp2%correlation
      e_total = scf%energy + e_pt2

      tag = "_ri"
      if (.not. use_df) tag = "_conv"
      open (newunit=unit, file="/tmp/mqc_dh_"//functional//"_"//basis//trim(tag)//".txt", &
            status="replace", action="write")
      write (unit, "(3(es25.16e3,1x))") e_total, scf%energy, e_pt2
      close (unit)

      write (*, "(A,A,A,A,A,F18.10,A,F18.10,A,F14.10)") "  DH ", functional, &
         trim(tag), " ", ": E = ", e_total, "  KS = ", scf%energy, "  PT2 = ", e_pt2
      call xc%destroy()
      call aux%destroy()
      call mol%destroy()
   end subroutine dh_case

   subroutine scf_case(basis, functional, level)
      !! A self-consistent Kohn-Sham energy, for the Python side to compare
      !!
      !! The point of milestone 2: the potential reaches the Fock matrix and the
      !! SCF converges on it. Nothing here is new machinery -- it is
      !! `run_libcint_rhf` with one extra argument, which is the whole design.
      character(len=*), intent(in) :: basis, functional
      integer, intent(in) :: level

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc
      type(error_t) :: err
      real(dp) :: c(3, 3)
      integer :: unit
      character(len=8) :: lvl

      c = reshape([0.0_dp, 0.00000000009155_dp*ANG, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])

      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[dft] basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call xc_context_create(mol, functional, xc, err, level=level)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[dft] ", functional, " context failed: ", err%get_message()
         failures = failures + 1
         call mol%destroy()
         return
      end if

      call run_libcint_rhf(mol, 10, 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err, &
                           in_core=.true., xc=xc)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A,A,A,A)") "[dft] ", functional, " SCF failed: ", err%get_message()
         failures = failures + 1
         call xc%destroy()
         call mol%destroy()
         return
      end if

      write (lvl, "(I0)") level
      open (newunit=unit, file="/tmp/mqc_ks_"//functional//"_"//basis//"_L"//trim(lvl)// &
            ".txt", status="replace", action="write")
      write (unit, "(es25.16e3)") scf%energy
      write (unit, "(I0)") scf%iterations
      close (unit)

      write (*, "(A,A,A,A,A,I0,A,F18.10,A,I0)") "  KS ", functional, "/", basis, &
         " level ", level, ": E = ", scf%energy, "  iterations ", scf%iterations
      call xc%destroy()
      call mol%destroy()
   end subroutine scf_case

   subroutine one_case(basis, level)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: level

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(dft_grid_t) :: grid
      type(error_t) :: err
      type(xc_f03_func_t) :: fx, fc
      real(dp), allocatable :: ao(:, :), rho(:), ex(:), ec(:)
      real(dp) :: c(3, 3)
      real(dp) :: n_elec, e_x, e_c
      integer :: unit, ig
      character(len=8) :: lvl

      c = reshape([0.0_dp, 0.00000000009155_dp*ANG, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])

      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[dft] basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call run_libcint_rhf(mol, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                           in_core=.true.)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A)") "[dft] SCF failed"
         failures = failures + 1
         return
      end if

      call build_dft_grid(mol%coords, [8, 1, 1], grid, err, level=level)
      if (err%has_error()) then
         write (*, "(A,A)") "[dft] grid failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call eval_ao_block(mol, grid%coords, ao, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[dft] AO evaluation failed: ", err%get_message()
         failures = failures + 1
         return
      end if
      call eval_rho(ao, scf%density, rho)

      ! The exact check first: an energy from a density with the wrong number of
      ! electrons in it is not worth printing.
      n_elec = sum(grid%weights*rho)
      if (abs(n_elec - 10.0_dp) > 1.0e-4_dp) then
         write (*, "(A,F14.8)") "[dft] integrated density is not ten electrons: ", n_elec
         failures = failures + 1
         return
      end if

      ! LDA: Slater exchange plus VWN correlation, each as an energy density per
      ! electron, so the quadrature carries a factor of rho.
      allocate (ex(grid%n_points), ec(grid%n_points))
      call xc_f03_func_init(fx, XC_LDA_X, XC_UNPOLARIZED)
      call xc_f03_lda_exc(fx, int(grid%n_points, 8), rho, ex)
      call xc_f03_func_end(fx)
      call xc_f03_func_init(fc, XC_LDA_C_VWN, XC_UNPOLARIZED)
      call xc_f03_lda_exc(fc, int(grid%n_points, 8), rho, ec)
      call xc_f03_func_end(fc)

      e_x = sum(grid%weights*rho*ex)
      e_c = sum(grid%weights*rho*ec)

      write (lvl, "(I0)") level
      open (newunit=unit, file="/tmp/mqc_dft_"//basis//"_L"//trim(lvl)//".txt", &
            status="replace", action="write")
      write (unit, "(I0,1x,I0)") grid%n_points, mol%nao
      write (unit, "(3(es25.16e3,1x))") n_elec, e_x, e_c
      do ig = 1, grid%n_points
         write (unit, "(4(es25.16e3,1x))") grid%coords(:, ig), grid%weights(ig)
      end do
      do ig = 1, grid%n_points
         write (unit, "(es25.16e3)") rho(ig)
      end do
      close (unit)

      write (*, "(A,A,A,I0,A,I0,A,F16.10,A,F14.10,A,F14.10)") &
         "  ", basis, " level ", level, ": points ", grid%n_points, &
         "  N = ", n_elec, "  Ex = ", e_x, "  Ec = ", e_c

      call grid%destroy()
      call mol%destroy()
   end subroutine one_case

end program check_dft
