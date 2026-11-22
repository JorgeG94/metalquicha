!! Nuclear magnetic shielding tensors, against PySCF and against themselves
module test_mqc_czt_nmr
   !! The pinned tensors are `pyscf.prop.nmr.RHF` on water in 6-31G at this
   !! geometry, fed this repository's basis: all nine elements of each
   !! nucleus's tensor, which fixes every sign and index convention at once.
   !! The tensor is not symmetric, so its two off-diagonal elements in the
   !! molecular plane are separate checks.
   !!
   !! Beside them, what needs no other program:
   !!
   !!   * moving the molecule leaves every tensor unchanged, which is the
   !!     property gauge-including orbitals exist to give and which a
   !!     common-origin calculation in this basis does not have;
   !!   * the diamagnetic and paramagnetic parts add up to the total;
   !!   * the summary's isotropic value is a third of the trace and its
   !!     anisotropy that of a tensor whose principal values are known.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_nmr, only: czt_nmr_shielding
   use mqc_czt_giao, only: giao_coulomb_exchange
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_integrals, only: shell_dim
   use libcint_fortran, only: libcint_2e_ig1_sph
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_czt_pcm, only: pcm_context_t
   use mqc_method_config, only: pcm_config_t
   use mqc_nmr_summary, only: shielding_summary
   implicit none
   private

   public :: collect_mqc_czt_nmr_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   ! PySCF, water/6-31G, ppm. Rows are the field's direction.
   real(dp), parameter :: O_REF(3, 3) = reshape( &
                          [331.646944_dp, 0.0_dp, -17.748779_dp, &
                           0.0_dp, 307.327544_dp, 0.0_dp, &
                           -17.701335_dp, 0.0_dp, 341.403792_dp], [3, 3])
   real(dp), parameter :: H1_REF(3, 3) = reshape( &
                          [27.797327_dp, 0.0_dp, -2.421815_dp, &
                           0.0_dp, 25.909842_dp, 0.0_dp, &
                           -0.918204_dp, 0.0_dp, 44.357334_dp], [3, 3])
   real(dp), parameter :: H2_REF(3, 3) = reshape( &
                          [42.551405_dp, 0.0_dp, -4.740960_dp, &
                           0.0_dp, 25.938788_dp, 0.0_dp, &
                           -6.247660_dp, 0.0_dp, 29.677511_dp], [3, 3])
   real(dp), parameter :: PPM_TOL = 5.0e-5_dp
      !! The references are printed to six decimals.

   ! The same molecule under two functionals, written as PySCF prints them,
   ! a row per field direction. PBE is a pure GGA, so the field induces no
   ! potential and the response is uncoupled; B3LYP keeps a fifth of the
   ! exchange and with it a coupled response.
   real(dp), parameter :: PBE_REF(3, 3, 3) = reshape( &
                          [324.108712_dp, 0.0_dp, -12.598363_dp, &
                           0.0_dp, 296.121611_dp, 0.0_dp, &
                           -12.546327_dp, 0.0_dp, 331.274399_dp, &
                           28.925377_dp, 0.0_dp, -3.116130_dp, &
                           0.0_dp, 26.651564_dp, 0.0_dp, &
                           -0.644932_dp, 0.0_dp, 43.899715_dp, &
                           42.082937_dp, 0.0_dp, -4.056759_dp, &
                           0.0_dp, 26.674670_dp, 0.0_dp, &
                           -6.532487_dp, 0.0_dp, 30.802582_dp], [3, 3, 3])
   real(dp), parameter :: B3LYP_REF(3, 3, 3) = reshape( &
                          [324.452054_dp, 0.0_dp, -14.467035_dp, &
                           0.0_dp, 297.645243_dp, 0.0_dp, &
                           -14.414577_dp, 0.0_dp, 332.581913_dp, &
                           28.748962_dp, 0.0_dp, -2.947949_dp, &
                           0.0_dp, 26.515985_dp, 0.0_dp, &
                           -0.722904_dp, 0.0_dp, 43.913947_dp, &
                           42.108277_dp, 0.0_dp, -4.186792_dp, &
                           0.0_dp, 26.539775_dp, 0.0_dp, &
                           -6.415966_dp, 0.0_dp, 30.616770_dp], [3, 3, 3])
      !! Stored column by column, so each block above is one atom's tensor
      !! transposed: the first three numbers are its first column.

contains

   subroutine collect_mqc_czt_nmr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_tensors_match_pyscf", test_water), &
                  new_unittest("translation_leaves_the_tensors_unchanged", test_translation), &
                  new_unittest("kohn_sham_tensors_match_pyscf", test_kohn_sham), &
                  new_unittest("kohn_sham_translation_and_refusals", test_kohn_sham_gauge), &
                  new_unittest("continuum_shieldings_do_not_move", test_continuum), &
                  new_unittest("folded_two_electron_term_is_the_full_sum", test_folding), &
                  new_unittest("stored_and_direct_response_agree", test_response_routes), &
                  new_unittest("summary_of_a_known_tensor", test_summary) &
                  ]
   end subroutine collect_mqc_czt_nmr_tests

   subroutine shield_water(shift, sigma, dia, para, err, functional, dielectric)
      !! Water in 6-31G, moved by `shift`, and its shielding tensors
      real(dp), intent(in) :: shift(3)
      real(dp), allocatable, intent(out) :: sigma(:, :, :), dia(:, :, :), para(:, :, :)
      type(error_t), intent(inout) :: err
      character(len=*), intent(in), optional :: functional
         !! Kohn-Sham with this functional; absent is Hartree-Fock
      real(dp), intent(in), optional :: dielectric
         !! In a conductor-like continuum of this dielectric constant;
         !! Hartree-Fock only here. Absent is the gas phase.

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc
      type(pcm_context_t) :: pcm
      type(pcm_config_t) :: config
      real(dp) :: c(3, 3)
      integer :: i

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      do i = 1, 3
         c(:, i) = c(:, i) + shift
      end do
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
      if (err%has_error()) return
      if (present(functional)) then
         call xc_context_create(mol, functional, xc, err)
         if (err%has_error()) return
         call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                          grad_tol=1.0e-10_dp, xc=xc)
         if (err%has_error()) return
         call czt_nmr_shielding(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                                sigma, err, diamagnetic=dia, paramagnetic=para, xc=xc)
         call mol%destroy()
         return
      end if
      if (present(dielectric)) then
         config%enabled = .true.
         config%method = "cpcm"
         config%dielectric = dielectric
         config%angular_points = 302
         call pcm%build(mol, [8, 1, 1], config, err)
         if (err%has_error()) return
         call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                          grad_tol=1.0e-10_dp, pcm=pcm)
         if (err%has_error()) return
         call czt_nmr_shielding(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                                sigma, err, diamagnetic=dia, paramagnetic=para, pcm=pcm)
         call pcm%destroy()
         call mol%destroy()
         return
      end if
      call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                       grad_tol=1.0e-10_dp)
      if (err%has_error()) return
      call czt_nmr_shielding(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                             sigma, err, diamagnetic=dia, paramagnetic=para)
      call mol%destroy()
   end subroutine shield_water

   subroutine test_folding(error)
      !! The two-electron field derivative against its formulas as written.
      !!
      !! The routine computes a quarter of the shell quartets at most: only
      !! pairs on different atoms, each pair in one order, and the second
      !! pair in one order, with every value put in the four places the full
      !! sum would have put it, and it skips what its screen calls
      !! negligible. Here the full sum is done with none of that -- every
      !! ordered quartet, contracted exactly as
      !!
      !!     vj_ij = sum_kl (g i j|k l) D_lk,      vk_il = sum_jk (g i j|k l) D_jk
      !!
      !! reads -- and the two have to agree. That is the sign of the
      !! antisymmetry, the factor of two on the Coulomb term, the four
      !! exchange placements and the screen at once. The molecule is moved
      !! off the origin because the screen's estimate of the gauge factor
      !! grows with the distance from it, and a molecule at the origin would
      !! not notice an estimate that was too small.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: vj(:, :, :), vk(:, :, :), vj_full(:, :, :), vk_full(:, :, :)
      real(dp), allocatable :: vj_plain(:, :, :), vk_plain(:, :, :), bounds(:, :), buf(:)
      real(dp) :: c(3, 3)
      integer :: shls(4)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, i, j, k, l, comp, idx, n

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      do i = 1, 3
         c(:, i) = c(:, i) + [3.0_dp, -5.0_dp, 7.0_dp]
      end do
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
      call run_czt_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "water must converge: "//err%get_full_trace())
      if (allocated(error)) return
      n = mol%nao

      allocate (vj_full(n, n, 3), vk_full(n, n, 3), buf(3*6**4))
      vj_full = 0.0_dp
      vk_full = 0.0_dp
      do ish = 1, mol%nbas
         do jsh = 1, mol%nbas
            do ksh = 1, mol%nbas
               do lsh = 1, mol%nbas
                  di = shell_dim(mol%cartesian, ish - 1, mol%bas)
                  dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
                  dk = shell_dim(mol%cartesian, ksh - 1, mol%bas)
                  dl = shell_dim(mol%cartesian, lsh - 1, mol%bas)
                  io = mol%shell_offset(ish)
                  jo = mol%shell_offset(jsh)
                  ko = mol%shell_offset(ksh)
                  lo = mol%shell_offset(lsh)
                  shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
                  if (libcint_2e_ig1_sph(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, &
                                         mol%env) == 0) cycle
                  idx = 0
                  do comp = 1, 3
                     do l = lo + 1, lo + dl
                        do k = ko + 1, ko + dk
                           do j = jo + 1, jo + dj
                              do i = io + 1, io + di
                                 idx = idx + 1
                                 vj_full(i, j, comp) = vj_full(i, j, comp) &
                                                       + buf(idx)*scf%density(l, k)
                                 vk_full(i, l, comp) = vk_full(i, l, comp) &
                                                       + buf(idx)*scf%density(j, k)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      call schwarz_bounds(mol, bounds, err)
      call giao_coulomb_exchange(mol, scf%density, vj, vk, err, bounds=bounds)
      call giao_coulomb_exchange(mol, scf%density, vj_plain, vk_plain, err)
      call check(error,.not. err%has_error(), "field derivative: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(4x,a,2es10.2)") "folded and screened against the full sum, J and K ", &
         maxval(abs(vj - vj_full)), maxval(abs(vk - vk_full))
      write (*, "(4x,a,2es10.2)") "folded, unscreened, against the full sum, J and K ", &
         maxval(abs(vj_plain - vj_full)), maxval(abs(vk_plain - vk_full))
      call check(error, maxval(abs(vj_full)) > 1.0e-2_dp .and. maxval(abs(vk_full)) > 1.0e-2_dp, &
                 "the reference sum is empty: nothing was compared")
      if (allocated(error)) return
      call check(error, maxval(abs(vj_plain - vj_full)) < 1.0e-11_dp &
                 .and. maxval(abs(vk_plain - vk_full)) < 1.0e-11_dp, &
                 "the folded two-electron field derivative is not the full sum")
      if (allocated(error)) return
      call check(error, maxval(abs(vj - vj_full)) < 1.0e-10_dp &
                 .and. maxval(abs(vk - vk_full)) < 1.0e-10_dp, &
                 "the screen dropped quartets of the field derivative that matter")
      if (allocated(error)) return

      ! With no exact exchange to feed, the Coulomb matrix is the same and
      ! the exchange matrix is left at zero.
      call giao_coulomb_exchange(mol, scf%density, vj, vk, err, bounds=bounds, &
                                 coulomb_only=.true.)
      call check(error, maxval(abs(vj - vj_full)) < 1.0e-10_dp, &
                 "the Coulomb-only field derivative differs from the full one's")
      if (allocated(error)) return
      call check(error, maxval(abs(vk)) <= 0.0_dp, &
                 "the Coulomb-only field derivative filled the exchange matrix")
      call mol%destroy()
   end subroutine test_folding

   subroutine test_response_routes(error)
      !! The field response gives the same tensors whether its coupling
      !! matrix is built once by transforming the integrals or the integrals
      !! are recomputed on every pass.
      !!
      !! Water is small enough that the stored route is always taken, so
      !! without this the integral-direct route -- the one every molecule too
      !! large for the matrix takes -- would be exercised by nothing. The
      !! bound is the solver's own: both converge the same equations to the
      !! same residual.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: stored(:, :, :), direct(:, :, :)
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
      call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                       grad_tol=1.0e-10_dp)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "water must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call czt_nmr_shielding(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                             stored, err)
      call czt_nmr_shielding(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                             direct, err, integral_direct=.true.)
      call check(error,.not. err%has_error(), "shielding: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(4x,a,es10.2)") "stored against integral-direct response, ppm ", &
         maxval(abs(stored - direct))
      call check(error, maxval(abs(stored - direct)) < 1.0e-6_dp, &
                 "the stored and integral-direct field responses disagree")
      call mol%destroy()
   end subroutine test_response_routes

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: sigma(:, :, :), dia(:, :, :), para(:, :, :)
      integer :: iatom, x

      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], sigma, dia, para, err)
      call check(error,.not. err%has_error(), "shielding: "//err%get_full_trace())
      if (allocated(error)) return
      do iatom = 1, 3
         write (*, "(4x,a,i0)") "atom ", iatom
         do x = 1, 3
            write (*, "(6x,3f14.6)") sigma(x, :, iatom)
         end do
      end do
      call check(error, maxval(abs(sigma(:, :, 1) - O_REF)) < PPM_TOL, &
                 "oxygen's shielding tensor does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(sigma(:, :, 2) - H1_REF)) < PPM_TOL, &
                 "the first hydrogen's shielding tensor does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(sigma(:, :, 3) - H2_REF)) < PPM_TOL, &
                 "the second hydrogen's shielding tensor does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(sigma - dia - para)) < 1.0e-10_dp, &
                 "the two parts do not add up to the total")
      if (allocated(error)) return
      ! The tensor is not symmetric, and the reference says by how much; a
      ! transposed result would pass a symmetric check and fail this one.
      call check(error, abs(sigma(1, 3, 2) - sigma(3, 1, 2)) > 1.0_dp, &
                 "the hydrogen tensor should be visibly asymmetric")
   end subroutine test_water

   subroutine test_translation(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: here(:, :, :), there(:, :, :), dia(:, :, :), para(:, :, :)
      real(dp), allocatable :: dia_moved(:, :, :)

      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], here, dia, para, err)
      call shield_water([3.0_dp, -5.0_dp, 7.0_dp], there, dia_moved, para, err)
      call check(error,.not. err%has_error(), "shielding: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(4x,a,es10.2)") "largest change on moving the molecule ", &
         maxval(abs(there - here))
      call check(error, maxval(abs(there - here)) < 1.0e-6_dp, &
                 "the shielding changed when the molecule was moved")
   end subroutine test_translation

   subroutine test_kohn_sham(error)
      !! The exchange-correlation term over gauge-including orbitals, and the
      !! scaling of exact exchange in a hybrid, against `pyscf.prop.nmr.rks`.
      !!
      !! That module's quadrature loop fails under the current PySCF, so the
      !! numbers come from the same formulas run in one pass, and that copy
      !! reproduces the values in the module's own self-test -- hydrogen
      !! fluoride at LDA and at the B3LYP of its day -- to 1e-4 ppm.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: sigma(:, :, :), dia(:, :, :), para(:, :, :)

      if (.not. xc_available()) return
      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], sigma, dia, para, err, functional="pbe")
      call check(error,.not. err%has_error(), "PBE shielding: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(4x,a,es10.2)") "PBE, worst element against PySCF   ", &
         maxval(abs(sigma - PBE_REF))
      call check(error, maxval(abs(sigma - PBE_REF)) < PPM_TOL, &
                 "the PBE shielding tensors do not match PySCF's")
      if (allocated(error)) return

      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], sigma, dia, para, err, functional="b3lyp")
      call check(error,.not. err%has_error(), "B3LYP shielding: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(4x,a,es10.2)") "B3LYP, worst element against PySCF ", &
         maxval(abs(sigma - B3LYP_REF))
      call check(error, maxval(abs(sigma - B3LYP_REF)) < PPM_TOL, &
                 "the B3LYP shielding tensors do not match PySCF's")
      if (allocated(error)) return
      ! The functional matters: correlation deshields oxygen by ten ppm.
      call check(error, O_REF(1, 1) - sigma(1, 1, 1) > 5.0_dp, &
                 "the Kohn-Sham tensor is Hartree-Fock's")
   end subroutine test_kohn_sham

   subroutine test_kohn_sham_gauge(error)
      !! Moving the molecule changes nothing under a functional either,
      !! which is the exchange-correlation term's gauge factor being right:
      !! without it the shieldings move by ppm. And what is not built says so.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: here(:, :, :), there(:, :, :), dia(:, :, :), para(:, :, :)

      if (.not. xc_available()) return
      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], here, dia, para, err, functional="pbe")
      call shield_water([3.0_dp, -5.0_dp, 7.0_dp], there, dia, para, err, functional="pbe")
      call check(error,.not. err%has_error(), "shielding: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(4x,a,es10.2)") "PBE, largest change on moving the molecule ", &
         maxval(abs(there - here))
      call check(error, maxval(abs(there - here)) < 1.0e-5_dp, &
                 "the PBE shielding changed when the molecule was moved")
      if (allocated(error)) return

      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], here, dia, para, err, functional="tpss")
      call check(error, err%has_error(), "a meta-GGA shielding was not refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "meta-GGA") > 0, &
                 "the meta-GGA refusal does not name the reason: "//err%get_message())
      if (allocated(error)) return
      call err%clear()
      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], here, dia, para, err, &
                        functional="cam-b3lyp")
      call check(error, err%has_error() .and. &
                 index(err%get_message(), "range-separated") > 0, &
                 "a range-separated shielding was not refused by name")
   end subroutine test_kohn_sham_gauge

   subroutine test_continuum(error)
      !! Shieldings in a continuum, which no other program here computes.
      !!
      !! What can be checked without one is the thing the continuum's term
      !! exists for. The surface charges' potential sits in the Fock matrix,
      !! so its derivative through the gauge factor has to sit in the
      !! first-order one, or the shielding depends on where the molecule is:
      !! with that term left out, moving water by the vector below changes
      !! oxygen's shielding by 87 ppm. With it the change is the response
      !! solve's residual. Any error of sign or size in the term shows here.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: EPS_WATER = 78.3553_dp
      type(error_t) :: err
      real(dp), allocatable :: gas(:, :, :), here(:, :, :), there(:, :, :)
      real(dp), allocatable :: dia(:, :, :), para(:, :, :)
      real(dp) :: iso_gas, iso_here

      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], gas, dia, para, err)
      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], here, dia, para, err, dielectric=EPS_WATER)
      call shield_water([3.0_dp, -5.0_dp, 7.0_dp], there, dia, para, err, &
                        dielectric=EPS_WATER)
      call check(error,.not. err%has_error(), "shielding: "//err%get_full_trace())
      if (allocated(error)) return
      iso_gas = (gas(1, 1, 1) + gas(2, 2, 1) + gas(3, 3, 1))/3.0_dp
      iso_here = (here(1, 1, 1) + here(2, 2, 1) + here(3, 3, 1))/3.0_dp
      write (*, "(4x,a,es10.2)") "in the continuum, largest change on moving ", &
         maxval(abs(there - here))
      write (*, "(4x,a,2f12.4)") "oxygen, gas phase and continuum            ", &
         iso_gas, iso_here
      call check(error, maxval(abs(there - here)) < 1.0e-5_dp, &
                 "the shielding in a continuum changed when the molecule was moved")
      if (allocated(error)) return
      ! And the continuum did something: oxygen's shielding rises by twelve
      ! ppm in water's dielectric.
      call check(error, abs(iso_here - iso_gas - 11.68_dp) < 0.05_dp, &
                 "the continuum's effect on oxygen's shielding is not the one measured")
   end subroutine test_continuum

   subroutine test_summary(error)
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: tensor(3, 3), iso, aniso

      ! Principal values 10, 20 and 60: the mean is 30 and 60 is furthest
      ! from it, so the anisotropy is 60 - 15. The antisymmetric part must
      ! not enter either number.
      tensor = 0.0_dp
      tensor(1, 1) = 10.0_dp
      tensor(2, 2) = 20.0_dp
      tensor(3, 3) = 60.0_dp
      tensor(1, 2) = 4.0_dp
      tensor(2, 1) = -4.0_dp
      call shielding_summary(tensor, iso, aniso)
      call check(error, abs(iso - 30.0_dp) < 1.0e-12_dp, "isotropic value")
      if (allocated(error)) return
      call check(error, abs(aniso - 45.0_dp) < 1.0e-10_dp, "anisotropy")
      if (allocated(error)) return

      ! And when the outlier is the smallest principal value: -30, 20, 40
      ! average 10, -30 is furthest, so the anisotropy is -30 - 30.
      tensor = 0.0_dp
      tensor(1, 1) = 40.0_dp
      tensor(2, 2) = -30.0_dp
      tensor(3, 3) = 20.0_dp
      call shielding_summary(tensor, iso, aniso)
      call check(error, abs(aniso + 60.0_dp) < 1.0e-10_dp, &
                 "anisotropy when the outlying principal value is the lowest")
   end subroutine test_summary

end module test_mqc_czt_nmr

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_nmr, only: collect_mqc_czt_nmr_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_nmr", collect_mqc_czt_nmr_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
