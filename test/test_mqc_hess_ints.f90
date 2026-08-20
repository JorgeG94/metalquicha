!! Second derivatives of the one-electron integrals
module test_mqc_hess_ints
   !! These are checked against identities rather than against stored numbers,
   !! because the identities say something the numbers cannot: that the two
   !! derivative orderings belong to the same integral and are consistent with
   !! each other.
   !!
   !! They were also checked against PySCF once, while being written -- all six
   !! blocks agreed to between 1e-16 and 1e-13. That is worth doing and not
   !! worth keeping: it needs Python, and it would pin the layout to one
   !! external library's conventions rather than to anything true.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_hess_ints, only: hess_1e_block, HESS_OVLP_II, HESS_OVLP_IJ, &
                                    HESS_KIN_II, HESS_KIN_IJ, HESS_NUC_II, HESS_NUC_IJ, &
                                    hess_2e_block, HESS_ERI_II, HESS_ERI_IJ, HESS_ERI_IK, &
                                    hess_rinv_block, HESS_RINV_II, HESS_RINV_IJ
   use mqc_libcint_hessian, only: hcore_deriv_atom, make_h1_atom, overlap_deriv_atom, &
                                  solve_mo1_atom, nuclear_repulsion_hessian, partial_hessian, &
                                  response_hessian, rhf_hessian, hessian_to_matrix
   use mqc_vibrational_analysis, only: compute_vibrational_analysis
   use mqc_libcint_hess_ints, only: eri_ip1_block
   use pic_blas_interfaces, only: pic_gemm
   use mqc_libcint_hess_ints, only: eri_ip1_block
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   implicit none
   private

   public :: collect_mqc_hess_ints_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])

contains

   subroutine collect_mqc_hess_ints_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("translations_leave_the_integrals_alone", translations_leave_them_alone), &
                  new_unittest("derivatives_on_one_centre_commute", derivatives_commute), &
                  new_unittest("every_block_is_populated", every_block_is_populated), &
                  new_unittest("components_are_where_they_belong", components_are_where_they_belong), &
                  new_unittest("two_electron_blocks", two_electron_blocks), &
           new_unittest("hcore_derivative_is_translationally_invariant", hcore_derivative_is_translationally_invariant), &
                  new_unittest("the_perturbation_sums_to_nothing", the_perturbation_sums_to_nothing), &
                  new_unittest("overlap_derivative_moves_only_one_atom", overlap_derivative_moves_only_one_atom), &
                  new_unittest("first_order_density_against_finite_difference", first_order_density_fd), &
                  new_unittest("an_unknown_selector_is_refused", an_unknown_selector_is_refused), &
                  new_unittest("nuclear_repulsion_hessian_is_right", nuclear_repulsion_hessian_ok), &
                  new_unittest("rinv_blocks_sum_to_the_nuclear_one", rinv_blocks_sum_to_nuc), &
                  new_unittest("partial_hessian_against_finite_difference", partial_hessian_fd), &
                  new_unittest("partial_hessian_is_translationally_invariant", partial_hessian_translates), &
                  new_unittest("rhf_hessian_against_finite_difference", rhf_hessian_fd) &
                  ]
   end subroutine collect_mqc_hess_ints_tests

   subroutine build_water(mol, err, ok)
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      ok = .not. err%has_error()
   end subroutine build_water

   subroutine frozen_pieces(geo, dens, weight, energy, err)
      !! The energy expression at one geometry, with the density held fixed
      !!
      !! This is the function whose second derivative the partial Hessian is --
      !! the SCF Lagrangian with `D` and `W` frozen at the reference geometry
      !! and only the integrals allowed to follow the nuclei. Differencing it
      !! needs no response theory and no converged SCF at the displaced point,
      !! which is what makes it a check on the partial Hessian *alone* rather
      !! than on the whole assembly.
      real(dp), intent(in) :: geo(3, 3)
      real(dp), intent(in) :: dens(:, :), weight(:, :)
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      real(dp), allocatable :: h(:, :), s(:, :), eri(:, :, :, :)
      integer :: i, j, k, l, n
      real(dp) :: two

      call build_libcint_molecule(WATER_Z, WATER_SYM, geo, "sto-3g", mol, err)
      if (err%has_error()) return

      call mol%core_hamiltonian(h)
      call mol%overlap(s)
      call mol%eris(eri)
      n = size(h, 1)

      two = 0.0_dp
      do l = 1, n
         do k = 1, n
            do j = 1, n
               do i = 1, n
                  two = two + eri(i, j, k, l)*(0.5_dp*dens(i, j)*dens(k, l) &
                                               - 0.25_dp*dens(i, l)*dens(k, j))
               end do
            end do
         end do
      end do

      energy = sum(dens*h) + two - sum(weight*s) + mol%nuclear_repulsion()

      deallocate (h, s, eri)
      call mol%destroy()
   end subroutine frozen_pieces

   subroutine reference_scf(dens, weight, mol, err)
      !! Water in STO-3G, and the two matrices the partial Hessian contracts
      real(dp), allocatable, intent(out) :: dens(:, :), weight(:, :)
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      type(rhf_result_t) :: scf
      integer :: i, n, nocc

      nocc = 5
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err)
      if (err%has_error()) return

      dens = scf%density
      n = size(dens, 1)
      allocate (weight(n, n))
      weight = 0.0_dp
      do i = 1, nocc
         weight = weight + 2.0_dp*scf%orbital_energies(i) &
                  *matmul(reshape(scf%orbitals(:, i), [n, 1]), &
                          reshape(scf%orbitals(:, i), [1, n]))
      end do
   end subroutine reference_scf

   subroutine partial_hessian_fd(error)
      !! Every element of the frozen-density Hessian, against central differences
      !!
      !! **The check that decides whether the assembly is right.** The
      !! translation and permutation identities below and above are necessary
      !! and nowhere near sufficient: a term dropped from all sixteen
      !! two-electron contractions at once satisfies every one of them, and so
      !! does a factor of two on a piece that is itself symmetric. Differencing
      !! the energy expression pins the magnitudes, one atom pair at a time.
      !!
      !! Four displaced evaluations per element, so the mixed second derivative
      !! is `(f(+,+) - f(+,-) - f(-,+) + f(-,-)) / 4h^2` and the error is
      !! second order in the step. The step is where that error and the
      !! cancellation of nearly equal energies are both tolerable.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), weight(:, :), hess(:, :, :, :), nn(:, :, :, :)
      real(dp) :: geo(3, 3), pp, pm, mp_, mm, fd, worst, scale
      real(dp), parameter :: H = 2.0e-3_dp
      integer :: ia, ja, a, b

      call reference_scf(dens, weight, mol, err)
      call check(error, .not. err%has_error(), "the reference did not converge")
      if (allocated(error)) return

      call partial_hessian(mol, dens, weight, hess, err)
      call nuclear_repulsion_hessian(WATER_Z, WATER, nn, err)
      call check(error, .not. err%has_error(), "the partial Hessian did not build")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      hess = hess + nn

      worst = 0.0_dp
      scale = 0.0_dp
      do ja = 1, 3
         do ia = 1, 3
            do b = 1, 3
               do a = 1, 3
                  geo = WATER
                  geo(a, ia) = geo(a, ia) + H
                  geo(b, ja) = geo(b, ja) + H
                  call frozen_pieces(geo, dens, weight, pp, err)
                  geo = WATER
                  geo(a, ia) = geo(a, ia) + H
                  geo(b, ja) = geo(b, ja) - H
                  call frozen_pieces(geo, dens, weight, pm, err)
                  geo = WATER
                  geo(a, ia) = geo(a, ia) - H
                  geo(b, ja) = geo(b, ja) + H
                  call frozen_pieces(geo, dens, weight, mp_, err)
                  geo = WATER
                  geo(a, ia) = geo(a, ia) - H
                  geo(b, ja) = geo(b, ja) - H
                  call frozen_pieces(geo, dens, weight, mm, err)
                  if (err%has_error()) exit

                  fd = (pp - pm - mp_ + mm)/(4.0_dp*H*H)
                  worst = max(worst, abs(fd - hess(a, b, ia, ja)))
                  scale = max(scale, abs(fd))
               end do
            end do
         end do
      end do

      call check(error, .not. err%has_error(), "a displaced evaluation failed")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      ! The largest element is about `0.35` Hartree per Bohr squared here, so
      ! this only refuses a comparison against nothing.
      call check(error, scale > 0.1_dp, &
                 "the frozen energy barely curves, so this comparison proves nothing")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      ! One part in `1e4` of the largest element. Central differences at this
      ! step carry about `1.4e-5` of relative error on this molecule, which is
      ! the `h^2` term and not the Hessian's; the margin is five times that.
      call check(error, worst < 1.0e-4_dp*scale, &
                 "the partial Hessian disagrees with finite differences")
      call mol%destroy()
   end subroutine partial_hessian_fd

   subroutine partial_hessian_translates(error)
      !! Moving the whole molecule costs nothing, at second order too
      !!
      !! Weaker than the finite differences above and worth keeping anyway,
      !! because it is the identity that fails when a term is deposited on the
      !! wrong atom rather than computed wrongly -- and the Hellmann-Feynman
      !! terms, which reach atom pairs that carry no basis functions between
      !! them, are exactly where that can happen.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), weight(:, :), hess(:, :, :, :), nn(:, :, :, :)
      real(dp) :: row(3, 3), worst, scale
      integer :: ja, ia

      call reference_scf(dens, weight, mol, err)
      call check(error, .not. err%has_error(), "the reference did not converge")
      if (allocated(error)) return

      call partial_hessian(mol, dens, weight, hess, err)
      call nuclear_repulsion_hessian(WATER_Z, WATER, nn, err)
      call check(error, .not. err%has_error(), "the partial Hessian did not build")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      hess = hess + nn

      scale = maxval(abs(hess))
      worst = 0.0_dp
      do ja = 1, 3
         row = 0.0_dp
         do ia = 1, 3
            row = row + hess(:, :, ia, ja)
         end do
         worst = max(worst, maxval(abs(row)))
      end do

      call check(error, scale > 0.1_dp, "the Hessian is empty")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call check(error, worst < 1.0e-8_dp*scale, &
                 "a rigid translation does not leave the partial Hessian alone")

      ! And it is symmetric under swapping the two atoms with the two
      ! components together, which the deposit rules do not enforce.
      worst = 0.0_dp
      do ja = 1, 3
         do ia = 1, 3
            worst = max(worst, maxval(abs(hess(:, :, ia, ja) - transpose(hess(:, :, ja, ia)))))
         end do
      end do
      call check(error, worst < 1.0e-9_dp*scale, "the partial Hessian is not symmetric")
      call mol%destroy()
   end subroutine partial_hessian_translates

   subroutine scf_energy_at(geo, energy, err)
      !! A converged restricted Hartree-Fock energy at one geometry
      real(dp), intent(in) :: geo(3, 3)
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      energy = 0.0_dp
      call build_libcint_molecule(WATER_Z, WATER_SYM, geo, "sto-3g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err)
      if (err%has_error()) return
      energy = scf%energy
      call mol%destroy()
   end subroutine scf_energy_at

   subroutine rhf_hessian_fd(error)
      !! The whole thing, against differencing converged energies
      !!
      !! **This is the test the rest of the file exists to make debuggable.**
      !! Everything else checks one piece against an identity or against a
      !! frozen-density energy; this checks the assembled Hessian against the
      !! quantity it claims to be -- the curvature of the Hartree-Fock energy
      !! surface, obtained by converging an SCF at each of four displaced
      !! geometries per element.
      !!
      !! It is also the only check that says anything about the response at
      !! all, because the response has no meaning on its own: it is the
      !! difference between the energy of a molecule whose density is allowed
      !! to follow the nuclei and one whose density is not, and only the first
      !! of those is an observable.
      !!
      !! Note what the response is worth here. It is not a small correction --
      !! the frozen-density Hessian is off by tens of percent on the diagonal
      !! and gets some off-diagonal blocks qualitatively wrong.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: scf
      real(dp), allocatable :: hess(:, :, :, :)
      real(dp), allocatable :: mine(:, :), theirs(:, :)
      real(dp), allocatable :: f_mine(:), f_theirs(:), rm(:), fc(:), disp(:, :)
      real(dp) :: geo(3, 3), pp, pm, mp_, mm, fd, worst, scale
      real(dp), parameter :: H = 2.0e-3_dp
      integer :: ia, ja, a, b, nocc, k

      nocc = 5
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call run_libcint_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err)
      call check(error, .not. err%has_error(), "the reference did not converge")
      if (allocated(error)) return

      call rhf_hessian(mol, WATER_Z, scf%density, scf%orbitals, &
                       scf%orbital_energies, nocc, hess, err)
      call check(error, .not. err%has_error(), "the Hessian did not build")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      allocate (theirs(9, 9))
      worst = 0.0_dp
      scale = 0.0_dp
      do ja = 1, 3
         do ia = 1, 3
            do b = 1, 3
               do a = 1, 3
                  geo = WATER
                  geo(a, ia) = geo(a, ia) + H
                  geo(b, ja) = geo(b, ja) + H
                  call scf_energy_at(geo, pp, err)
                  geo = WATER
                  geo(a, ia) = geo(a, ia) + H
                  geo(b, ja) = geo(b, ja) - H
                  call scf_energy_at(geo, pm, err)
                  geo = WATER
                  geo(a, ia) = geo(a, ia) - H
                  geo(b, ja) = geo(b, ja) + H
                  call scf_energy_at(geo, mp_, err)
                  geo = WATER
                  geo(a, ia) = geo(a, ia) - H
                  geo(b, ja) = geo(b, ja) - H
                  call scf_energy_at(geo, mm, err)
                  if (err%has_error()) exit

                  fd = (pp - pm - mp_ + mm)/(4.0_dp*H*H)
                  theirs(3*(ia - 1) + a, 3*(ja - 1) + b) = fd
                  worst = max(worst, abs(fd - hess(a, b, ia, ja)))
                  scale = max(scale, abs(fd))
               end do
            end do
         end do
      end do

      call check(error, .not. err%has_error(), "a displaced SCF did not converge")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call check(error, scale > 0.1_dp, "the surface barely curves")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      ! Agreement is 4.5e-6 relative, which is the `h^2` term of the
      ! difference formula rather than anything in the Hessian; the gate is at
      ! twenty times that.
      call check(error, worst < 1.0e-4_dp*scale, &
                 "the analytic Hessian disagrees with finite differences")

      ! Both hold for the assembled Hessian and neither is implied by the
      ! finite differences, which compare element by element and would pass on
      ! a Hessian that agreed everywhere and was symmetric nowhere.
      worst = 0.0_dp
      do ja = 1, 3
         do ia = 1, 3
            worst = max(worst, maxval(abs(hess(:, :, ia, ja) - transpose(hess(:, :, ja, ia)))))
         end do
      end do
      call check(error, worst < 1.0e-8_dp*scale, "the Hessian is not symmetric")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      worst = 0.0_dp
      do ja = 1, 3
         worst = max(worst, maxval(abs(hess(:, :, 1, ja) + hess(:, :, 2, ja) + hess(:, :, 3, ja))))
      end do
      call check(error, worst < 1.0e-8_dp*scale, &
                 "a rigid translation does not leave the Hessian alone")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! ---- and the same comparison in the units anyone actually reads -------
      !
      ! The matrix agreeing element by element does not by itself mean the
      ! frequencies do, because between the two lies `hessian_to_matrix` and
      ! the mass weighting it feeds. Interleaving the atom and Cartesian
      ! indices the wrong way round leaves a matrix that is still symmetric,
      ! still translationally invariant and still agrees with nothing -- it
      ! only shows up once masses are attached, which is here.
      !
      ! Both Hessians go through the same analysis, so what is compared is the
      ! Hessians and not the analysis.
      call hessian_to_matrix(hess, mine)
      call compute_vibrational_analysis(mine, WATER_Z, f_mine, rm, fc, disp, &
                                        coordinates=WATER, project_trans_rot=.true.)
      deallocate (rm, fc, disp)
      call compute_vibrational_analysis(theirs, WATER_Z, f_theirs, rm, fc, disp, &
                                        coordinates=WATER, project_trans_rot=.true.)

      ! The three real modes, which are the last three after projection.
      worst = 0.0_dp
      scale = 0.0_dp
      do k = 7, 9
         worst = max(worst, abs(f_mine(k) - f_theirs(k)))
         scale = max(scale, abs(f_theirs(k)))
      end do
      call check(error, scale > 1000.0_dp, "no vibrations came back")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      ! A tenth of a wavenumber on modes of a few thousand. The residual is the
      ! finite-difference Hessian's, not this one's.
      call check(error, worst < 0.5_dp, &
                 "the analytic and finite-difference frequencies disagree")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! Six modes projected out, and they should be gone rather than small.
      worst = 0.0_dp
      do k = 1, 6
         worst = max(worst, abs(f_mine(k)))
      end do
      call check(error, worst < 1.0_dp, &
                 "the translations and rotations did not come out at zero")
      call mol%destroy()
   end subroutine rhf_hessian_fd

   subroutine rinv_blocks_sum_to_nuc(error)
      !! `-sum_C Z_C rinv(C)` is the nuclear attraction operator, at any order
      !!
      !! The whole reason the `rinv` blocks exist is that `ipipnuc` has already
      !! summed over the nuclei and so cannot say which one moved. This checks
      !! the two really are the same operator taken apart and put back
      !! together -- if the origin were being placed on the wrong atom, or the
      !! nine components came back permuted relative to `ipipnuc`, the sum
      !! would not close.
      !!
      !! It is a strong check on placement and layout and a weak one on value:
      !! whatever `ipiprinv` gets wrong, `ipipnuc` almost certainly gets wrong
      !! the same way, since one is built from the other. What pins the value
      !! is the finite-difference Hessian downstream.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      logical :: ok
      real(dp), allocatable :: nuc(:, :, :), rinv(:, :, :), total(:, :, :)
      integer :: c, which
      real(dp) :: worst

      call build_water(mol, err, ok)
      call check(error, ok, "water did not build")
      if (allocated(error)) return

      do which = 1, 2
         if (which == 1) then
            call hess_1e_block(mol, HESS_NUC_II, nuc, err)
         else
            call hess_1e_block(mol, HESS_NUC_IJ, nuc, err)
         end if
         call check(error, .not. err%has_error(), "nuclear block failed")
         if (allocated(error)) return

         allocate (total(size(nuc, 1), size(nuc, 2), size(nuc, 3)))
         total = 0.0_dp
         do c = 1, mol%natm
            if (which == 1) then
               call hess_rinv_block(mol, c, HESS_RINV_II, rinv, err)
            else
               call hess_rinv_block(mol, c, HESS_RINV_IJ, rinv, err)
            end if
            call check(error, .not. err%has_error(), "rinv block failed")
            if (allocated(error)) return
            total = total - mol%charges(c)*rinv
            deallocate (rinv)
         end do

         worst = maxval(abs(total - nuc))
         call check(error, worst < 1.0e-11_dp, "rinv blocks do not sum to the nuclear one")
         if (allocated(error)) return
         deallocate (total, nuc)
      end do
   end subroutine rinv_blocks_sum_to_nuc

   subroutine translations_leave_them_alone(error)
      !! Moving both centres together cannot change a two-centre integral
      !!
      !! `S` depends on where the bra and ket sit, so under a rigid translation
      !!
      !!     (d/dA + d/dB)^2 S = 0   ->   S_AA + 2 S_AB + S_BB = 0
      !!
      !! and `S_BB` is `S_AA` with the two indices exchanged. **This is the
      !! check that matters most here**, because it is exactly what fails when
      !! only one derivative ordering is used, or when the two are mixed up: a
      !! Hessian built from the wrong combination is not slightly wrong, it is
      !! one where a rigid translation costs energy, and the symptom is
      !! translational modes that no longer come out at zero frequency.
      !!
      !! Nuclear attraction is excluded on purpose. Its operator depends on the
      !! nuclei too, so translating the two basis centres alone does not leave
      !! it alone and this identity does not apply to it.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: aa(:, :, :), ab(:, :, :)
      logical :: ok
      integer :: n, c
      real(dp) :: worst

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_1e_block(mol, HESS_OVLP_II, aa, err)
      call hess_1e_block(mol, HESS_OVLP_IJ, ab, err)
      call check(error,.not. err%has_error(), "overlap blocks failed")
      if (.not. allocated(error)) then
         n = size(aa, 1)
         worst = 0.0_dp
         do c = 1, 9
            worst = max(worst, maxval(abs(aa(:, :, c) + transpose(aa(:, :, c)) &
                                          + 2.0_dp*ab(:, :, c))))
         end do
         call check(error, worst < 1.0e-10_dp, &
                    "the overlap second derivatives are not translationally invariant")
      end if
      if (allocated(aa)) deallocate (aa)
      if (allocated(ab)) deallocate (ab)
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call hess_1e_block(mol, HESS_KIN_II, aa, err)
      call hess_1e_block(mol, HESS_KIN_IJ, ab, err)
      call check(error,.not. err%has_error(), "kinetic blocks failed")
      if (.not. allocated(error)) then
         worst = 0.0_dp
         do c = 1, 9
            worst = max(worst, maxval(abs(aa(:, :, c) + transpose(aa(:, :, c)) &
                                          + 2.0_dp*ab(:, :, c))))
         end do
         call check(error, worst < 1.0e-10_dp, &
                    "the kinetic second derivatives are not translationally invariant")
      end if
      call mol%destroy()
   end subroutine translations_leave_them_alone

   subroutine derivatives_commute(error)
      !! Two derivatives on the same centre commute, so xy equals yx
      !!
      !! True of the `ipip` blocks and **not** of the `ipXip` ones, where the
      !! two derivatives act on different centres and the component index means
      !! something asymmetric. Asserting it of both would be wrong; asserting it
      !! of neither would miss a transposed component block.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: aa(:, :, :)
      logical :: ok

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_1e_block(mol, HESS_OVLP_II, aa, err)
      call check(error,.not. err%has_error(), "overlap block failed")
      if (.not. allocated(error)) then
         ! xy against yx, xz against zx, yz against zy
         call check(error, maxval(abs(aa(:, :, 2) - aa(:, :, 4))) < 1.0e-12_dp, &
                    "xy and yx differ on one centre")
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(aa(:, :, 3) - aa(:, :, 7))) < 1.0e-12_dp, &
                    "xz and zx differ on one centre")
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(aa(:, :, 6) - aa(:, :, 8))) < 1.0e-12_dp, &
                    "yz and zy differ on one centre")
      end if
      call mol%destroy()
   end subroutine derivatives_commute

   subroutine every_block_is_populated(error)
      !! Each selector returns a differently-shaped, non-zero array
      !!
      !! A dispatch that fell through to the same entry point for two selectors
      !! would satisfy every identity above, since both copies would be equally
      !! valid integrals. What catches it is that the six are genuinely
      !! different numbers.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: a(:, :, :), b(:, :, :)
      logical :: ok
      integer :: k, sel(6)

      sel = [HESS_OVLP_II, HESS_OVLP_IJ, HESS_KIN_II, HESS_KIN_IJ, &
             HESS_NUC_II, HESS_NUC_IJ]

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      do k = 1, 6
         call hess_1e_block(mol, sel(k), a, err)
         call check(error,.not. err%has_error(), "a block failed to build")
         if (allocated(error)) exit
         call check(error, maxval(abs(a)) > 1.0e-8_dp, "a block came back empty")
         if (allocated(error)) exit
         if (k > 1) then
            call check(error, maxval(abs(a - b)) > 1.0e-8_dp, &
                       "two selectors returned the same integral")
            if (allocated(error)) exit
         end if
         if (allocated(b)) deallocate (b)
         b = a
         deallocate (a)
      end do
      call mol%destroy()
   end subroutine every_block_is_populated

   subroutine components_are_where_they_belong(error)
      !! Which of the nine components is which
      !!
      !! **The identities above cannot see this.** A permutation of the
      !! component index applied consistently to every block satisfies
      !! translational invariance and the commuting test alike, because each of
      !! those holds component by component -- so a wrong stride in the
      !! unpacking reproduces every structural property and still hands back xy
      !! where xz belongs.
      !!
      !! Pinned as a norm per component rather than as individual elements.
      !! Elements are a bad choice here and the first version of this test found
      !! out the hard way: water in this orientation is full of symmetry zeros,
      !! so two elements picked by eye can sit where a scrambled stride happens
      !! to read an equal value, and the test passes while the layout is wrong.
      !! A norm over the whole block cannot be fooled that way -- moving values
      !! between components changes what each one sums to.
      !!
      !! The numbers came from PySCF, which agreed with this implementation to
      !! between 1e-16 and 1e-13 across all six blocks when it was written.
      !! Recorded here because that comparison needs Python and cannot run in
      !! this suite.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: a(:, :, :)
      logical :: ok
      integer :: c
      real(dp), parameter :: TOL = 1.0e-5_dp
      !> `int1e_ipipovlp`, sum of absolute values per component
      real(dp), parameter :: OVLP_NORM(9) = [ &
                             28.464886_dp, 2.724628_dp, 2.566234_dp, &
                             2.724628_dp, 27.896234_dp, 3.212813_dp, &
                             2.566234_dp, 3.212813_dp, 27.811801_dp]
      !> `int1e_ipnucip`, which is a different shape entirely and pins the
      !> dispatch as well as the layout
      real(dp), parameter :: NUC_NORM(9) = [ &
                             1209.762119_dp, 57.194237_dp, 55.848487_dp, &
                             57.194237_dp, 1217.340172_dp, 67.229615_dp, &
                             55.848487_dp, 67.229615_dp, 1210.321515_dp]

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_1e_block(mol, HESS_OVLP_II, a, err)
      call check(error,.not. err%has_error(), "overlap block failed")
      if (.not. allocated(error)) then
         do c = 1, 9
            call check(error, sum(abs(a(:, :, c))), OVLP_NORM(c), thr=TOL, &
                       more="an overlap component is not where it belongs")
            if (allocated(error)) exit
         end do
      end if
      if (allocated(a)) deallocate (a)
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call hess_1e_block(mol, HESS_NUC_IJ, a, err)
      call check(error,.not. err%has_error(), "nuclear block failed")
      if (.not. allocated(error)) then
         do c = 1, 9
            call check(error, sum(abs(a(:, :, c))), NUC_NORM(c), thr=1.0e-4_dp, &
                       more="a nuclear component is not where it belongs")
            if (allocated(error)) exit
         end do
      end if
      call mol%destroy()
   end subroutine components_are_where_they_belong

   subroutine two_electron_blocks(error)
      !! The three two-electron second derivatives
      !!
      !! Same two things the one-electron blocks are held to, for the same
      !! reasons. Two derivatives on the same centre commute, so `ipip1` is
      !! symmetric in its component index and `ip1ip2` -- one derivative on the
      !! bra and one on the ket -- is not, which is why only the first is
      !! checked that way.
      !!
      !! The per-component norms pin the layout, which the symmetry cannot: a
      !! stride error that permutes components consistently leaves every
      !! symmetry intact. They came from PySCF, which agreed with this
      !! implementation to between 3e-15 and 3e-14 on all three blocks when it
      !! was written.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: e(:, :, :, :, :)
      logical :: ok
      integer :: c
      real(dp), parameter :: TOL = 1.0e-4_dp
      real(dp), parameter :: IPIP1_NORM(9) = [ &
                             486.37821_dp, 45.44354_dp, 43.45087_dp, &
                             45.44354_dp, 480.76078_dp, 54.79232_dp, &
                             43.45087_dp, 54.79232_dp, 481.26138_dp]
      real(dp), parameter :: IP1IP2_NORM(9) = [ &
                             70.54127_dp, 46.69994_dp, 43.90662_dp, &
                             46.69994_dp, 89.39699_dp, 56.08919_dp, &
                             43.90662_dp, 56.08919_dp, 81.08118_dp]

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_2e_block(mol, HESS_ERI_II, e, err)
      call check(error,.not. err%has_error(), "ipip1 failed")
      if (.not. allocated(error)) then
         ! Both derivatives on centre one, so the component index is symmetric.
         call check(error, maxval(abs(e(:, :, :, :, 2) - e(:, :, :, :, 4))) < 1.0e-10_dp, &
                    "xy and yx differ with both derivatives on one centre")
      end if
      if (.not. allocated(error)) then
         do c = 1, 9
            call check(error, sum(abs(e(:, :, :, :, c))), IPIP1_NORM(c), thr=TOL, &
                       more="an ipip1 component is not where it belongs")
            if (allocated(error)) exit
         end do
      end if
      if (allocated(e)) deallocate (e)
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! A different integral entirely, so this pins the dispatch as well.
      call hess_2e_block(mol, HESS_ERI_IK, e, err)
      call check(error,.not. err%has_error(), "ip1ip2 failed")
      if (.not. allocated(error)) then
         do c = 1, 9
            call check(error, sum(abs(e(:, :, :, :, c))), IP1IP2_NORM(c), thr=TOL, &
                       more="an ip1ip2 component is not where it belongs")
            if (allocated(error)) exit
         end do
      end if
      if (.not. allocated(error)) then
         ! **A norm over the whole array cannot see a permutation inside it**,
         ! and the quartet unpacking has three indices that could be transposed
         ! against each other. A single element is no better: the first one
         ! tried sat in a quartet where two shell dimensions were both one, so
         ! the correct and transposed strides evaluated to the same index and
         ! the break was invisible.
         !
         ! Slice norms are pinned below and they do catch a component-level
         ! scramble. **They do not catch every index permutation**, and neither
         ! does anything else here: transposing the two ket indices in the
         ! unpacking was tried deliberately and survives all five tests. What
         ! rules that out is the comparison this file cannot run -- the whole
         ! `n_ao^4 x 9` array was checked against PySCF element by element when
         ! this was written, agreeing to 3e-15 on all three blocks. That is the
         ! verification of record for the quartet layout; what is below is
         ! regression cover, and the distinction is worth keeping honest.
         call check(error, sum(abs(e(:, :, 1, :, 2))), 11.128496_dp, thr=1.0e-4_dp, &
                    more="the third ket index is not where it belongs")
      end if
      if (.not. allocated(error)) then
         call check(error, sum(abs(e(:, :, :, 1, 2))), 11.553514_dp, thr=1.0e-4_dp, &
                    more="the fourth ket index is not where it belongs")
      end if
      call mol%destroy()
   end subroutine two_electron_blocks

   subroutine hcore_derivative_is_translationally_invariant(error)
      !! Moving every atom together cannot change the core Hamiltonian
      !!
      !!     sum_A dH/dR_A = 0
      !!
      !! exactly, because a rigid translation is not a change of the molecule.
      !! This is worth more than it looks: `dH/dR_A` has a basis-function part
      !! that touches only atom `A`'s block and a Hellmann-Feynman part that
      !! touches everything, and the two enter with opposite signs and different
      !! index ranges. Getting either wrong -- the sign the library's nabla
      !! implies, the nuclear charge, which block the basis term lands in --
      !! leaves a residue here.
      !!
      !! Checked against PySCF's `hcore_generator` while being written, matching
      !! to 4e-16 on every atom. That comparison needs Python; this one does
      !! not.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: h(:, :, :), total(:, :, :)
      logical :: ok
      integer :: a

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      do a = 1, 3
         call hcore_deriv_atom(mol, a, h, err)
         call check(error,.not. err%has_error(), "a core Hamiltonian derivative failed")
         if (allocated(error)) exit
         if (a == 1) then
            allocate (total(size(h, 1), size(h, 2), 3))
            total = 0.0_dp
         end if
         total = total + h
         deallocate (h)
      end do

      if (.not. allocated(error)) then
         call check(error, maxval(abs(total)) < 1.0e-10_dp, &
                    "the core Hamiltonian derivatives do not sum to zero over atoms")
      end if
      if (.not. allocated(error)) then
         ! **The sum rule does not see a missing kinetic term.** Its own
         ! contribution sums to zero over atoms independently, so dropping it
         ! altogether leaves the identity above satisfied -- which is what
         ! happened when that was tried. So the magnitude is pinned as well,
         ! per component and per atom, against the values PySCF's
         ! `hcore_generator` gave when this agreed with it to 4e-16.
         call hcore_deriv_atom(mol, 1, h, err)
         call check(error, sum(abs(h(:, :, 1))), 6.542609_dp, thr=1.0e-5_dp, &
                    more="the oxygen x component has the wrong magnitude")
      end if
      if (.not. allocated(error)) then
         call check(error, sum(abs(h(:, :, 2))), 28.624459_dp, thr=1.0e-5_dp, &
                    more="the oxygen y component has the wrong magnitude")
      end if
      if (.not. allocated(error)) then
         deallocate (h)
         call hcore_deriv_atom(mol, 2, h, err)
         call check(error, sum(abs(h(:, :, 3))), 13.350037_dp, thr=1.0e-5_dp, &
                    more="the hydrogen z component has the wrong magnitude")
      end if
      call mol%destroy()
   end subroutine hcore_derivative_is_translationally_invariant

   subroutine the_perturbation_sums_to_nothing(error)
      !! The per-atom perturbation, summed over atoms, vanishes
      !!
      !! `h1_A` is the core Hamiltonian and the mean field differentiated with
      !! respect to atom `A`. Translating every atom together changes neither,
      !! so the sum over atoms is zero -- and it has to be zero for the whole
      !! Hessian, since this is what drives the response.
      !!
      !! **What this catches and what it misses.** It catches a mishandled
      !! index in the two-electron assembly, where a quartet contributes
      !! through each of its four positions and only the ones on this atom
      !! should count -- get the ownership test wrong and the sum stops
      !! cancelling. It does not catch a term that is translationally invariant
      !! on its own, which is why the magnitudes are pinned beside it, exactly
      !! as for the core Hamiltonian derivative next door.
      !!
      !! Checked against PySCF's own `make_h1` while being written: agreement
      !! to 4e-11 on a quantity of order one, which is the level the two SCF
      !! densities differ at rather than anything structural.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: scf
      real(dp), allocatable :: h(:, :, :), total(:, :, :), ip1(:, :, :, :, :)
      logical :: ok
      integer :: a
      real(dp), parameter :: PIN = 1.0e-4_dp

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return
      call run_libcint_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error(), "the reference did not converge")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call eri_ip1_block(mol, ip1, err)

      do a = 1, 3
         call make_h1_atom(mol, scf%density, ip1, a, h, err)
         call check(error,.not. err%has_error(), "a perturbation failed to build")
         if (allocated(error)) exit
         if (a == 1) then
            allocate (total(size(h, 1), size(h, 2), 3))
            total = 0.0_dp
         end if
         total = total + h
         if (a < 3) deallocate (h)
      end do

      if (.not. allocated(error)) then
         call check(error, maxval(abs(total)) < 1.0e-9_dp, &
                    "the perturbations do not sum to zero over atoms")
      end if
      if (.not. allocated(error)) then
         ! The last one built was atom 3, whose norms match atom 2 by symmetry.
         call check(error, sum(abs(h(:, :, 2))), 6.984370_dp, thr=PIN, &
                    more="the hydrogen y component has the wrong magnitude")
      end if
      if (.not. allocated(error)) then
         call check(error, sum(abs(h(:, :, 3))), 5.044607_dp, thr=PIN, &
                    more="the hydrogen z component has the wrong magnitude")
      end if
      call mol%destroy()
   end subroutine the_perturbation_sums_to_nothing

   subroutine overlap_derivative_moves_only_one_atom(error)
      !! `dS/dR_A` touches only the rows and columns of atom `A`
      !!
      !! Displacing a nucleus moves the functions centred on it and nothing
      !! else, so every element of `dS/dR_A` outside atom `A`'s rows and columns
      !! is exactly zero. And summed over atoms the whole thing vanishes,
      !! because translating the molecule does not change how its functions
      !! overlap.
      !!
      !! **This matrix is why a nuclear Hessian needs a different response
      !! solve from anything else here.** A field perturbation leaves the
      !! overlap alone, so every existing caller of the coupled-perturbed
      !! machinery has `dS = 0` and an orbital response with no
      !! occupied-occupied block. Here orthonormality has to be maintained while
      !! the functions move, and that block is fixed at minus half of this.
      !!
      !! Matches PySCF's `s1ao` to 1e-16 on every atom.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: s1(:, :, :), total(:, :, :)
      logical :: ok
      integer :: a, i, j, c
      real(dp) :: outside

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      do a = 1, 3
         call overlap_deriv_atom(mol, a, s1, err)
         call check(error,.not. err%has_error(), "an overlap derivative failed")
         if (allocated(error)) exit
         if (a == 1) then
            allocate (total(size(s1, 1), size(s1, 2), 3))
            total = 0.0_dp
         end if
         total = total + s1
         if (a == 2) then
            ! The oxygen owns basis functions 1 to 5 in this basis, so a
            ! hydrogen's derivative must vanish on that whole block.
            outside = 0.0_dp
            do c = 1, 3
               do j = 1, 5
                  do i = 1, 5
                     outside = max(outside, abs(s1(i, j, c)))
                  end do
               end do
            end do
            call check(error, outside < 1.0e-14_dp, &
                       "a hydrogen's overlap derivative reaches the oxygen block")
         end if
         if (allocated(error)) exit
         deallocate (s1)
      end do

      if (.not. allocated(error)) then
         call check(error, maxval(abs(total)) < 1.0e-12_dp, &
                    "the overlap derivatives do not sum to zero over atoms")
      end if
      if (.not. allocated(error)) then
         call overlap_deriv_atom(mol, 1, s1, err)
         call check(error, maxval(abs(s1)) > 1.0e-3_dp, &
                    "the derivative came back empty, so the sum proves nothing")
      end if
      call mol%destroy()
   end subroutine overlap_derivative_moves_only_one_atom

   subroutine density_at(geo, dens, orbitals, energies, mol, err)
      !! A converged density at one geometry, and the orbitals with it
      real(dp), intent(in) :: geo(3, 3)
      real(dp), allocatable, intent(out) :: dens(:, :), orbitals(:, :), energies(:)
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      type(rhf_result_t) :: scf

      call build_libcint_molecule(WATER_Z, WATER_SYM, geo, "sto-3g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err)
      if (err%has_error()) return
      dens = scf%density
      orbitals = scf%orbitals
      energies = scf%orbital_energies
   end subroutine density_at

   subroutine first_order_density_fd(error)
      !! The analytic first-order density, against differencing two SCFs
      !!
      !! **The check that needs nothing external.** Displacing a nucleus and
      !! differencing the converged densities gives `dD/dR` directly, in the
      !! atomic-orbital basis, with no reference implementation and no
      !! molecular-orbital phase convention to agree about -- and phases are
      !! exactly what makes a coefficient-by-coefficient comparison against
      !! another program meaningless.
      !!
      !! Central differences, so the error is second order in the step, and the
      !! step is chosen where that error and the SCF's own convergence noise
      !! are both small: too large and the quadratic term shows, too small and
      !! differencing two nearly equal densities loses the digits that carry
      !! the answer.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol, mol_p, mol_m
      type(error_t) :: err
      real(dp), allocatable :: d0(:, :), dp_(:, :), dm(:, :), c0(:, :), e0(:)
      real(dp), allocatable :: cp(:, :), ep(:), cm(:, :), em(:)
      real(dp), allocatable :: h1(:, :, :), s1(:, :, :), mo1(:, :, :)
      real(dp), allocatable :: ip1(:, :, :, :, :), half(:, :), analytic(:, :), fd(:, :)
      real(dp) :: geo(3, 3)
      real(dp), parameter :: H = 1.0e-3_dp
      integer :: n, nocc
      real(dp) :: worst, scale

      nocc = 5
      call density_at(WATER, d0, c0, e0, mol, err)
      call check(error,.not. err%has_error(), "the reference did not converge")
      if (allocated(error)) return
      n = size(d0, 1)

      ! Atom 1 along z, chosen because water lies in the yz plane so the
      ! displacement is not along a symmetry axis that makes the answer zero.
      geo = WATER
      geo(3, 1) = geo(3, 1) + H
      call density_at(geo, dp_, cp, ep, mol_p, err)
      geo = WATER
      geo(3, 1) = geo(3, 1) - H
      call density_at(geo, dm, cm, em, mol_m, err)
      call check(error,.not. err%has_error(), "a displaced point did not converge")
      if (allocated(error)) return

      allocate (fd(n, n))
      fd = (dp_ - dm)/(2.0_dp*H)

      call eri_ip1_block(mol, ip1, err)
      call make_h1_atom(mol, d0, ip1, 1, h1, err)
      call overlap_deriv_atom(mol, 1, s1, err)
      call solve_mo1_atom(mol, c0, e0, nocc, h1, s1, mo1, err)
      call check(error,.not. err%has_error(), "the response did not solve")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! dD/dx = 2 (C mo1 Cocc^T + transpose), the z component
      allocate (half(n, nocc), analytic(n, n))
      call pic_gemm(c0, mo1(:, :, 3), half, alpha=2.0_dp, beta=0.0_dp)
      call pic_gemm(half, c0(:, 1:nocc), analytic, transb="T")
      analytic = analytic + transpose(analytic)

      worst = maxval(abs(analytic - fd))
      scale = maxval(abs(fd))

      ! Not vacuous: the density really does move when the oxygen does, so
      ! agreeing to nothing is not the same as agreeing.
      call check(error, scale > 1.0e-2_dp, &
                 "the density barely moved, so this comparison proves nothing")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call check(error, worst < 1.0e-5_dp*scale + 1.0e-6_dp, &
                 "the analytic first-order density disagrees with finite differences")
      call mol%destroy()
      call mol_p%destroy()
      call mol_m%destroy()
   end subroutine first_order_density_fd

   subroutine an_unknown_selector_is_refused(error)
      !! A selector outside the named set is an error, not a default
      !!
      !! Both dispatches end in a `case default`, so without this an
      !! unrecognised selector returns whichever integral that branch happens to
      !! name -- a real integral of exactly the right shape, and the wrong one.
      !! Nothing downstream could tell.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: m(:, :, :)
      real(dp), allocatable :: e(:, :, :, :, :)
      logical :: ok

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_1e_block(mol, 99, m, err)
      call check(error, err%has_error(), "an unknown one-electron selector was accepted")
      call err%clear()
      if (.not. allocated(error)) then
         call hess_2e_block(mol, 99, e, err)
         call check(error, err%has_error(), &
                    "an unknown two-electron selector was accepted")
         call err%clear()
      end if
      call mol%destroy()
   end subroutine an_unknown_selector_is_refused

   pure function repulsion(z, geo) result(e)
      !! `sum_{A<B} Z_A Z_B / R_AB`, written out so the test differences it
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: geo(:, :)
      real(dp) :: e
      integer :: a, b

      e = 0.0_dp
      do a = 1, size(z)
         do b = a + 1, size(z)
            e = e + real(z(a), dp)*real(z(b), dp) &
                /norm2(geo(:, a) - geo(:, b))
         end do
      end do
   end function repulsion

   subroutine nuclear_repulsion_hessian_ok(error)
      !! The analytic nuclear Hessian against differencing the energy twice
      !!
      !! No electrons and no integrals, so this can be checked against the
      !! closed-form energy directly -- which makes it the one block of the
      !! Hessian that needs nothing but arithmetic to trust.
      !!
      !! Two structural properties come with it. Each row sums to zero over
      !! atoms, because translating the molecule does not change how its nuclei
      !! repel; and the whole thing is symmetric under exchanging the two atoms
      !! together with the two Cartesian directions, because it is a second
      !! derivative.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: hess(:, :, :, :)
      real(dp) :: geo(3, 3), fd, worst, plus_plus, plus_minus, minus_plus, minus_minus
      real(dp), parameter :: H = 1.0e-4_dp
      integer :: a, b, i, j

      call nuclear_repulsion_hessian(WATER_Z, WATER, hess, err)
      call check(error,.not. err%has_error(), "the nuclear Hessian failed")
      if (allocated(error)) return

      ! One mixed second derivative by central differences in both variables:
      ! the oxygen z against a hydrogen y.
      geo = WATER
      geo(3, 1) = geo(3, 1) + H; geo(2, 2) = geo(2, 2) + H
      plus_plus = repulsion(WATER_Z, geo)
      geo = WATER
      geo(3, 1) = geo(3, 1) + H; geo(2, 2) = geo(2, 2) - H
      plus_minus = repulsion(WATER_Z, geo)
      geo = WATER
      geo(3, 1) = geo(3, 1) - H; geo(2, 2) = geo(2, 2) + H
      minus_plus = repulsion(WATER_Z, geo)
      geo = WATER
      geo(3, 1) = geo(3, 1) - H; geo(2, 2) = geo(2, 2) - H
      minus_minus = repulsion(WATER_Z, geo)
      fd = (plus_plus - plus_minus - minus_plus + minus_minus)/(4.0_dp*H*H)

      call check(error, hess(3, 2, 1, 2), fd, thr=1.0e-6_dp, &
                 more="the mixed second derivative disagrees with finite differences")
      if (allocated(error)) return
      call check(error, abs(fd) > 1.0e-3_dp, &
                 "the reference derivative is zero, so that comparison proves nothing")
      if (allocated(error)) return

      ! **And one where the two Cartesian directions agree.** The isotropic
      ! term only enters when they do, so a mixed component alone cannot see it
      ! -- a wrong power there survives the check above, the translation sum and
      ! the symmetry, because all three treat that term consistently however
      ! wrong it is.
      geo = WATER
      geo(3, 1) = geo(3, 1) + H; geo(3, 2) = geo(3, 2) + H
      plus_plus = repulsion(WATER_Z, geo)
      geo = WATER
      geo(3, 1) = geo(3, 1) + H; geo(3, 2) = geo(3, 2) - H
      plus_minus = repulsion(WATER_Z, geo)
      geo = WATER
      geo(3, 1) = geo(3, 1) - H; geo(3, 2) = geo(3, 2) + H
      minus_plus = repulsion(WATER_Z, geo)
      geo = WATER
      geo(3, 1) = geo(3, 1) - H; geo(3, 2) = geo(3, 2) - H
      minus_minus = repulsion(WATER_Z, geo)
      fd = (plus_plus - plus_minus - minus_plus + minus_minus)/(4.0_dp*H*H)

      call check(error, hess(3, 3, 1, 2), fd, thr=1.0e-6_dp, &
                 more="the diagonal second derivative disagrees with finite differences")
      if (allocated(error)) return

      ! Translations are free: every row sums to zero over the second atom.
      worst = 0.0_dp
      do j = 1, 3
         do i = 1, 3
            do a = 1, 3
               worst = max(worst, abs(sum(hess(i, j, a, :))))
            end do
         end do
      end do
      call check(error, worst < 1.0e-10_dp, &
                 "the nuclear Hessian is not translationally invariant")
      if (allocated(error)) return

      ! And it is a second derivative, so exchanging both indices together
      ! leaves it alone.
      worst = 0.0_dp
      do b = 1, 3
         do a = 1, 3
            do j = 1, 3
               do i = 1, 3
                  worst = max(worst, abs(hess(i, j, a, b) - hess(j, i, b, a)))
               end do
            end do
         end do
      end do
      call check(error, worst < 1.0e-10_dp, "the nuclear Hessian is not symmetric")
   end subroutine nuclear_repulsion_hessian_ok

end module test_mqc_hess_ints

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_hess_ints, only: collect_mqc_hess_ints_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_hess_ints", collect_mqc_hess_ints_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
