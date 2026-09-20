!! Tamm-Dancoff singlet excitation energies, against PySCF
module test_mqc_czt_tddft
   !! The Layer 2 gates of `TDDFT_PLAN.md`: the TDA operator `A` itself, and
   !! the excitation energies a Davidson finds in it.
   !!
   !! ## Why the matrix is checked before the spectrum
   !!
   !! An excitation energy is one number out of an eigensolve, and almost every
   !! way of getting the operator wrong -- a factor of two on the Coulomb term,
   !! a missing exchange pass, the kernel at half weight -- moves it by
   !! something that still looks like an excitation energy. So the first two
   !! cases here reproduce the whole `10x10` matrix of H2O/STO-3G element by
   !! element, where a dropped term has nowhere to hide, and only then do the
   !! remaining cases ask the solver for roots.
   !!
   !! `A` is built by applying the operator to the ten unit vectors. That is
   !! the same construction PySCF's side of the reference used, so the two are
   !! compared as matrices and not as two different summaries of one.
   !!
   !! ## Phases, and what is compared elementwise
   !!
   !! The occupied-virtual index runs `idx = (i-1)*n_vir + a`, virtual fastest,
   !! which is the layout `response_product` and `mqc_davidson` already share.
   !! Two codes converging the same closed shell agree on the orbitals up to a
   !! sign per orbital, and `A_{ia,jb}` carries one factor of each of four MO
   !! phases, so its **off-diagonal** signs are not code-independent. The
   !! diagonal is (every phase appears squared), and so is the spectrum. Hence
   !! `abs` off the diagonal, and nothing weaker anywhere else.
   !!
   !! ## The references
   !!
   !! PySCF 2.14 fed this repository's own basis JSON through `bse_to_pyscf`,
   !! `conv_tol = 1e-12`, geometry and grid as recorded in `TDDFT_PLAN.md`;
   !! the matrices are transcribed from `TDDFT_PLAN_matrices.txt` beside it.
   !! Both files are working documents that are not committed, so the numbers
   !! live here.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_response_product, only: response_product
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_czt_bridge, only: run_czt_hf
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t, STATE_SPIN_SINGLET
   implicit none
   private

   public :: collect_mqc_czt_tddft_tests

   integer, parameter :: N_OV = 10
      !! Occupied-virtual rotations of H2O in STO-3G: five doubly occupied
      !! orbitals and two virtuals.

   !! The geometry every reference in this file was taken at, **in Bohr**.
   !!
   !! `TDDFT_PLAN.md` states it in Angstrom, and converting it here would put
   !! `BOHR_TO_ANGSTROM` inside the comparison: PySCF 2.14 carries the CODATA
   !! 2010 value and this program the 2018 one, four parts in 1e10 apart, and
   !! that alone moved the orbital energies by 5e-10 and every element of the
   !! matrices below with them. Worse, `-DMQC_CODATA_YEAR=2014` would then
   !! silently move the geometry a test hard-codes references for. The plan's
   !! Angstrom numbers times 1.8897261254578281, written out, and given to
   !! PySCF as Bohr.
   !!
   !! No symmetry is imposed on either side, so the C2v structure of the
   !! matrices below is something the codes reproduce rather than something
   !! either was told.
   real(dp), parameter :: WATER_BOHR(3, 3) = reshape([ &
                                                     0.0_dp, 0.0_dp, 0.22259084031767759_dp, &
                                                     0.0_dp, 1.4275992706554927_dp, -0.89036525099683572_dp, &
                                                     0.0_dp, -1.4275992706554927_dp, -0.89036525099683572_dp], [3, 3])

   real(dp), parameter :: RHF_STO3G_ENERGY = -74.963146800103_dp
   real(dp), parameter :: PBE_STO3G_ENERGY = -75.225769660158_dp

   !! The RHF singlet TDA matrix of H2O/STO-3G, column major -- which for this
   !! matrix is also row major: PySCF's own is symmetric to 1.7e-16 and what is
   !! written out is its symmetric part, so a transposed transcription would
   !! not show up as a failure and is not something to guard against here.
   real(dp), parameter :: RHF_A(N_OV, N_OV) = reshape([ &
                                           20.107272486240_dp, -0.000000000000_dp, 0.020734288381_dp, 0.000000000000_dp, &
                                           -0.000000000000_dp, -0.023185679792_dp, 0.021878446764_dp, 0.000000000000_dp, &
                                         -0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, 20.157156130688_dp, &
                                            0.000000000000_dp, 0.045924795840_dp, -0.001436020409_dp, 0.000000000000_dp, &
                                            0.000000000000_dp, 0.030592731572_dp, 0.000000000000_dp, -0.000000000000_dp, &
                                            0.020734288381_dp, 0.000000000000_dp, 1.463576630280_dp, -0.000000000000_dp, &
                                            0.000000000000_dp, 0.054921127790_dp, -0.063303237277_dp, 0.000000000000_dp, &
                                             0.000000000000_dp, 0.000000000000_dp, 0.000000000000_dp, 0.045924795840_dp, &
                                           -0.000000000000_dp, 1.508587745587_dp, -0.027721006971_dp, 0.000000000000_dp, &
                                           0.000000000000_dp, 0.019550756739_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                          -0.000000000000_dp, -0.001436020409_dp, 0.000000000000_dp, -0.027721006971_dp, &
                                           0.792350538909_dp, 0.000000000000_dp, -0.000000000000_dp, -0.041151660792_dp, &
                                           -0.000000000000_dp, 0.000000000000_dp, -0.023185679792_dp, 0.000000000000_dp, &
                                             0.054921127790_dp, 0.000000000000_dp, 0.000000000000_dp, 1.051801661240_dp, &
                                            -0.109774457818_dp, 0.000000000000_dp, 0.000000000000_dp, 0.000000000000_dp, &
                                            0.021878446764_dp, 0.000000000000_dp, -0.063303237277_dp, 0.000000000000_dp, &
                                          -0.000000000000_dp, -0.109774457818_dp, 0.648258819693_dp, -0.000000000000_dp, &
                                           -0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, 0.030592731572_dp, &
                                            0.000000000000_dp, 0.019550756739_dp, -0.041151660792_dp, 0.000000000000_dp, &
                                          -0.000000000000_dp, 0.724583875320_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                           -0.000000000000_dp, 0.000000000000_dp, 0.000000000000_dp, -0.000000000000_dp, &
                                          -0.000000000000_dp, 0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                          0.485080278078_dp, -0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                            0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, 0.000000000000_dp, &
                             -0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, 0.555807696384_dp], [N_OV, N_OV])
   real(dp), parameter :: RHF_A_EIG(N_OV) = [ &
                          0.485080278078_dp, 0.555807696384_dp, 0.617580730515_dp, &
                          0.705081200245_dp, 0.810167907032_dp, 1.067858555914_dp, &
                          1.478121961140_dp, 1.510111516358_dp, 20.107348349883_dp, &
                          20.157317666869_dp]
   !! PySCF's own converged orbital energies, so the orbital-energy diagonal
   !! can be taken back out of both matrices and the coupling compared on its
   !! own -- which is the part Layer 2 builds.
   real(dp), parameter :: RHF_A_EPS(7) = [ &
                          -20.242376966048_dp, -1.268534558650_dp, -0.616911145520_dp, &
                          -0.453874586492_dp, -0.391502303404_dp, 0.605693789318_dp, &
                          0.740404043489_dp]

   !! The same matrix for PBE at `grids.level = 5` -- 90064 points on this
   !! water. Ours is `xc_context_create(..., level=5)`: the same prescription,
   !! not the same points, which is what `TOL_GRID` allows for.
   real(dp), parameter :: PBE_A(N_OV, N_OV) = reshape([ &
                                           18.797981919748_dp, 0.000000000000_dp, 0.008233889644_dp, -0.000000000000_dp, &
                                           0.000000000000_dp, -0.018007273607_dp, 0.007008313487_dp, -0.000000000000_dp, &
                                           0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, 18.886699937829_dp, &
                                           -0.000000000000_dp, 0.027746929597_dp, 0.004204543167_dp, -0.000000000000_dp, &
                                           -0.000000000000_dp, 0.026303971585_dp, 0.000000000000_dp, -0.000000000000_dp, &
                                            0.008233889644_dp, -0.000000000000_dp, 1.318181772063_dp, 0.000000000000_dp, &
                                          -0.000000000000_dp, 0.126105419044_dp, -0.080262805272_dp, -0.000000000000_dp, &
                                            0.000000000000_dp, 0.000000000000_dp, -0.000000000000_dp, 0.027746929597_dp, &
                                             0.000000000000_dp, 1.332734083532_dp, 0.037529056894_dp, 0.000000000000_dp, &
                                           -0.000000000000_dp, 0.031095656494_dp, 0.000000000000_dp, -0.000000000000_dp, &
                                            0.000000000000_dp, 0.004204543167_dp, -0.000000000000_dp, 0.037529056894_dp, &
                                           0.781107194550_dp, -0.000000000000_dp, 0.000000000000_dp, -0.056737441515_dp, &
                                          -0.000000000000_dp, 0.000000000000_dp, -0.018007273607_dp, -0.000000000000_dp, &
                                            0.126105419044_dp, 0.000000000000_dp, -0.000000000000_dp, 1.030301623731_dp, &
                                           -0.120791940908_dp, -0.000000000000_dp, 0.000000000000_dp, 0.000000000000_dp, &
                                          0.007008313487_dp, -0.000000000000_dp, -0.080262805272_dp, -0.000000000000_dp, &
                                           0.000000000000_dp, -0.120791940908_dp, 0.558579997408_dp, -0.000000000000_dp, &
                                          -0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, 0.026303971585_dp, &
                                          -0.000000000000_dp, 0.031095656494_dp, -0.056737441515_dp, -0.000000000000_dp, &
                                           -0.000000000000_dp, 0.683095550042_dp, 0.000000000000_dp, -0.000000000000_dp, &
                                             0.000000000000_dp, 0.000000000000_dp, 0.000000000000_dp, 0.000000000000_dp, &
                                           -0.000000000000_dp, 0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, &
                                           0.419805809301_dp, 0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                            0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, 0.000000000000_dp, &
                              -0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, 0.504469845475_dp], [N_OV, N_OV])
   real(dp), parameter :: PBE_A_EIG(N_OV) = [ &
                          0.419805809301_dp, 0.504469845475_dp, 0.526334486314_dp, &
                          0.654247084959_dp, 0.806236716878_dp, 0.997927938532_dp, &
                          1.336370051014_dp, 1.382776200364_dp, 18.798006687739_dp, &
                          18.886782913103_dp]

   integer, parameter :: N_CCPVDZ_STATES = 5
      !! Roots asked for in the cc-pVDZ cases, which is how many the plan's
      !! tables carry.

   !! H2O/cc-pVDZ, the singlet TDA roots, in Hartree.
   !!
   !! Regenerated rather than transcribed, because the STO-3G matrices above
   !! could not be: `TDDFT_PLAN_matrices.txt` was taken with PySCF's *internal*
   !! basis tables, which differ from this repository's own JSON in the eighth
   !! decimal on a Pople set -- 2.4e-8 in the STO-3G energy, 2e-7 in the
   !! orbital energies, and more than the matrix bound in every element. The
   !! Dunning sets are unaffected: fed our JSON, PySCF reproduces the plan's
   !! cc-pVDZ table to 1.5e-9 in the worst root, which is where a
   !! ten-decimal table stops saying anything. These are the same numbers to
   !! twelve.
   !!
   !! Each is an eigenvalue of the explicit `A`, probed out of PySCF's
   !! `TDA.gen_vind` on the 95 unit vectors, so no iterative solver's
   !! tolerance is in them.
   real(dp), parameter :: RHF_CCPVDZ_ENERGY = -76.026767997355_dp
   real(dp), parameter :: RHF_CCPVDZ_TDA(N_CCPVDZ_STATES) = [ &
                          0.338692376910_dp, 0.403909409291_dp, 0.435472490718_dp, &
                          0.501267843887_dp, 0.552877377503_dp]
   real(dp), parameter :: PBE_CCPVDZ_ENERGY = -76.333481668156_dp
   real(dp), parameter :: PBE_CCPVDZ_TDA(N_CCPVDZ_STATES) = [ &
                          0.270707248060_dp, 0.339473916603_dp, 0.356815758949_dp, &
                          0.430817321107_dp, 0.511144167228_dp]
   real(dp), parameter :: B3LYP_CCPVDZ_ENERGY = -76.420393534532_dp
   real(dp), parameter :: B3LYP_CCPVDZ_TDA(N_CCPVDZ_STATES) = [ &
                          0.280680694354_dp, 0.348407154560_dp, 0.368070529290_dp, &
                          0.440162767131_dp, 0.516279354331_dp]
   real(dp), parameter :: CAM_CCPVDZ_ENERGY = -76.391813183377_dp
   real(dp), parameter :: CAM_CCPVDZ_TDA(N_CCPVDZ_STATES) = [ &
                          0.283543256194_dp, 0.353278202364_dp, 0.371105546022_dp, &
                          0.445649836132_dp, 0.517104469420_dp]

   !! Hartree-Fock carries no quadrature, so nothing but the integrals and the
   !! two SCF thresholds sits between the codes, and both are converged far
   !! below this. Measured on the STO-3G matrix: 7.7e-13 on the worst diagonal
   !! element, 6.0e-13 off it, 1.3e-12 on the orbital energies and 1.2e-12 on
   !! the spectrum. The bound is two orders over that, which leaves room for
   !! the compiler-to-compiler spread and still sits eight orders under the
   !! smallest fault it has to catch: halving the kernel, dropping one
   !! exchange term or symmetrising the wrong way each move elements by 1e-2
   !! and up.
   real(dp), parameter :: TOL_EXACT = 1.0e-10_dp

   !! A Kohn-Sham comparison is limited by the grid, not by either solver.
   !! `level = 5` is the same prescription on both sides and not the same
   !! points -- 90058 here against PySCF's 90064 on this water -- so what this
   !! bound measures is the difference between two level-5 quadratures of one
   !! functional. Measured on the STO-3G PBE matrix: 3.4e-10 on the worst
   !! diagonal element, 2.0e-11 off it, 3.4e-10 on the spectrum, and 1.7e-13
   !! on the total energy. `test_mqc_czt_cphf`'s CAM-B3LYP polarizability, the
   !! nearest existing comparison, lands at 5.4e-9 against the same reference
   !! under the same grid; 1e-7 clears the worse of the two by more than an
   !! order and leaves the compiler spread room, which the double-hybrid
   !! Hessian shows can reach 2e-9 on a quadrature of this kind.
   real(dp), parameter :: TOL_GRID = 1.0e-7_dp

   !! The Davidson's own floor, against a dense diagonalisation of the same
   !! operator. Roots are accepted on a residual of 1e-8 here, and the
   !! eigenvalue error near an eigenvector is second order in the vector
   !! error, so 1e-9 on the eigenvalue is a bound the solver clears by
   !! construction rather than one fitted to it.
   real(dp), parameter :: TOL_DAVIDSON = 1.0e-9_dp

   !! The cc-pVDZ Hartree-Fock roots. Looser than `TOL_EXACT` because these
   !! come out of the Davidson rather than off a dense diagonalisation, so the
   !! solver's own floor is in them as well as the integrals'.
   real(dp), parameter :: TOL_CCPVDZ_HF = 1.0e-8_dp

contains

   subroutine collect_mqc_czt_tddft_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("the_rhf_tda_matrix_matches_pyscf", test_rhf_matrix), &
                  new_unittest("the_rhf_tda_spectrum_matches_pyscf", test_rhf_eigenvalues), &
                  new_unittest("the_pbe_tda_matrix_matches_pyscf", test_pbe_matrix), &
                  new_unittest("davidson_finds_the_dense_sto3g_roots", test_davidson), &
                  new_unittest("cc_pvdz_rhf_singlets_match_the_table", test_ccpvdz_rhf), &
                  new_unittest("cc_pvdz_pbe_singlets_match_the_table", test_ccpvdz_pbe), &
                  new_unittest("cc_pvdz_b3lyp_singlets_match_the_table", test_ccpvdz_b3lyp), &
                  new_unittest("cc_pvdz_cam_b3lyp_singlets_match_the_table", test_ccpvdz_cam), &
                  new_unittest("no_states_asked_for_means_no_spectrum", test_no_states), &
                  new_unittest("rpa_and_triplets_are_refused_by_name", test_later_layers) &
                  ]
   end subroutine collect_mqc_czt_tddft_tests

   subroutine water_sto3g(mol, scf, ctx, err, functional)
      !! Converge H2O/STO-3G at the plan geometry, Hartree-Fock or Kohn-Sham
      !!
      !! Tighter than the suite's usual SCF because what is compared is an
      !! operator built from the orbitals, not an energy: the orbital error
      !! goes as the commutator rather than its square, so the energy
      !! threshold has to be carried well past where the energy itself has
      !! stopped moving.
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(xc_context_t), intent(out) :: ctx
      type(error_t), intent(inout) :: err
      character(len=*), intent(in), optional :: functional
         !! Absent is Hartree-Fock, and leaves `ctx` untouched.

      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], WATER_BOHR, "sto-3g", &
                              mol, err)
      if (err%has_error()) return

      if (present(functional)) then
         call xc_context_create(mol, functional, ctx, err, level=5)
         if (err%has_error()) return
         call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-10_dp, .false., scf, err, &
                          xc=ctx, grad_tol=1.0e-9_dp)
      else
         call run_czt_rhf(mol, 10, 200, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err, &
                          in_core=.true., grad_tol=1.0e-12_dp)
      end if
   end subroutine water_sto3g

   subroutine dense_tda(mol, scf, ctx, kohn_sham, a, err)
      !! The explicit TDA matrix, from the operator applied to unit vectors
      !!
      !! `A = (1/2)[(A+B) + (A-B)]`, both halves from `response_product`. The
      !! identity is exact rather than an approximation of one: `(A+B)` is
      !! `dEps + 4(ai|bj) - c_x[(ab|ij)+(aj|ib)] + 4 f_xc` and `(A-B)` is
      !! `dEps - c_x[(ab|ij)-(aj|ib)]`, whose half sum is
      !! `dEps + 2(ai|bj) - c_x(ab|ij) + 2 f_xc`, which is `A`.
      type(czt_molecule_t), intent(in) :: mol
      type(rhf_result_t), intent(in) :: scf
      type(xc_context_t), intent(inout) :: ctx
      logical, intent(in) :: kohn_sham
      real(dp), allocatable, intent(out) :: a(:, :)
      type(error_t), intent(inout) :: err

      real(dp), allocatable :: c_occ(:, :), c_vir(:, :), gaps(:, :), zero_h(:, :)
      real(dp), allocatable :: bounds(:, :), u(:, :, :), ap(:, :, :), am(:, :, :)
      integer, allocatable :: idx(:)
      integer :: n_ao, n_mo, n_occ, n_vir, n_ov, i, jj, aa

      if (err%has_error()) return

      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      n_occ = scf%n_occupied
      n_vir = n_mo - n_occ
      n_ov = n_vir*n_occ

      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir), gaps(n_vir, n_occ))
      allocate (zero_h(n_ao, n_ao))
      c_occ = scf%orbitals(:, 1:n_occ)
      c_vir = scf%orbitals(:, n_occ + 1:n_mo)
      zero_h = 0.0_dp
      do i = 1, n_occ
         do aa = 1, n_vir
            gaps(aa, i) = scf%orbital_energies(n_occ + aa) - scf%orbital_energies(i)
         end do
      end do

      call schwarz_bounds(mol, bounds, err)
      if (err%has_error()) return

      allocate (u(n_vir, n_occ, n_ov), ap(n_vir, n_occ, n_ov), am(n_vir, n_occ, n_ov))
      allocate (idx(n_ov), a(n_ov, n_ov))
      u = 0.0_dp
      do jj = 1, n_ov
         aa = mod(jj - 1, n_vir) + 1
         i = (jj - 1)/n_vir + 1
         u(aa, i, jj) = 1.0_dp
         idx(jj) = jj
      end do
      ap = 0.0_dp
      am = 0.0_dp

      if (kohn_sham) then
         call response_product(mol, c_occ, c_vir, gaps, zero_h, u, idx, n_ov, .false., &
                               ap, err, bounds=bounds, k_scale=ctx%exx_fraction, &
                               xc=ctx, reference=scf%density)
         if (.not. err%has_error()) &
            call response_product(mol, c_occ, c_vir, gaps, zero_h, u, idx, n_ov, .true., &
                                  am, err, bounds=bounds, k_scale=ctx%exx_fraction, &
                                  xc=ctx, reference=scf%density)
      else
         call response_product(mol, c_occ, c_vir, gaps, zero_h, u, idx, n_ov, .false., &
                               ap, err, bounds=bounds, k_scale=1.0_dp)
         if (.not. err%has_error()) &
            call response_product(mol, c_occ, c_vir, gaps, zero_h, u, idx, n_ov, .true., &
                                  am, err, bounds=bounds, k_scale=1.0_dp)
      end if
      if (err%has_error()) return

      do jj = 1, n_ov
         a(:, jj) = 0.5_dp*(reshape(ap(:, :, jj), [n_ov]) + reshape(am(:, :, jj), [n_ov]))
      end do
   end subroutine dense_tda

   subroutine compare_matrix(error, a, reference, tol, what)
      !! Diagonal with its sign, off-diagonal in absolute value
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: a(:, :), reference(:, :), tol
      character(len=*), intent(in) :: what

      real(dp) :: worst_diag, worst_off
      integer :: i, j, n

      n = size(a, 1)
      worst_diag = 0.0_dp
      worst_off = 0.0_dp
      do j = 1, n
         do i = 1, n
            if (i == j) then
               worst_diag = max(worst_diag, abs(a(i, j) - reference(i, j)))
            else
               worst_off = max(worst_off, abs(abs(a(i, j)) - abs(reference(i, j))))
            end if
         end do
      end do

      call check(error, worst_diag < tol, "the "//what//" TDA diagonal disagrees "// &
                 "with PySCF")
      if (allocated(error)) return
      call check(error, worst_off < tol, "the "//what//" TDA off-diagonal magnitudes "// &
                 "disagree with PySCF")
   end subroutine compare_matrix

   function eigenvalues_of(a, ok) result(values)
      !! Every eigenvalue of a small symmetric matrix, ascending
      real(dp), intent(in) :: a(:, :)
      logical, intent(out) :: ok
      real(dp), allocatable :: values(:)

      real(dp), allocatable :: work(:, :)
      integer :: info

      work = a
      allocate (values(size(a, 1)))
      call pic_syev(work, values, jobz="N", uplo="U", info=info)
      ok = info == 0
   end function eigenvalues_of

   subroutine test_rhf_matrix(error)
      !! The ten-by-ten singlet TDA matrix of H2O/STO-3G, element by element
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: ctx
      type(error_t) :: err
      real(dp), allocatable :: a(:, :)

      call water_sto3g(mol, scf, ctx, err)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "the Hartree-Fock reference failed: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! The energy and the orbital energies first, and separately. Every
      ! element below is built out of the orbitals, so a geometry, a basis or
      ! a convergence threshold that does not match the reference fails ten
      ! thousand comparisons at once and says nothing about the operator.
      call check(error, abs(scf%energy - RHF_STO3G_ENERGY) < TOL_EXACT, &
                 "the Hartree-Fock energy is not the one the reference matrix "// &
                 "was taken at, so the orbitals are not either")
      if (.not. allocated(error)) &
         call check(error, maxval(abs(scf%orbital_energies(1:7) - RHF_A_EPS)) < &
                    TOL_EXACT, "the Hartree-Fock orbital energies are not PySCF's")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call dense_tda(mol, scf, ctx, .false., a, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the TDA operator failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, size(a, 1) == N_OV, "H2O/STO-3G did not give ten "// &
                 "occupied-virtual rotations")
      if (allocated(error)) return

      call compare_matrix(error, a, RHF_A, TOL_EXACT, "Hartree-Fock")
   end subroutine test_rhf_matrix

   subroutine test_rhf_eigenvalues(error)
      !! Every root of that matrix, not only the ones a solver would look for
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: ctx
      type(error_t) :: err
      real(dp), allocatable :: a(:, :), values(:)
      logical :: ok

      call water_sto3g(mol, scf, ctx, err)
      if (.not. err%has_error()) call dense_tda(mol, scf, ctx, .false., a, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the Hartree-Fock TDA matrix failed: "// &
                 err%get_message())
      if (allocated(error)) return

      values = eigenvalues_of(a, ok)
      call check(error, ok, "the dense diagonalisation of the TDA matrix failed")
      if (allocated(error)) return
      call check(error, maxval(abs(values - RHF_A_EIG)) < TOL_EXACT, &
                 "the Hartree-Fock TDA spectrum disagrees with PySCF")
   end subroutine test_rhf_eigenvalues

   subroutine test_pbe_matrix(error)
      !! The same matrix for a pure functional, where the kernel is the whole
      !! two-electron difference
      !!
      !! PBE carries no exact exchange, so `(A-B)` reduces to the orbital
      !! energy differences and everything that is not the diagonal comes
      !! from the Coulomb term and the kernel. A kernel at the wrong weight
      !! has nothing to hide behind here.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: ctx
      type(error_t) :: err
      real(dp), allocatable :: a(:, :), values(:)
      logical :: ok

      if (.not. xc_available()) return

      call water_sto3g(mol, scf, ctx, err, functional="pbe")
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "the PBE reference failed: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! Looser than the matrix bound: the energy is a quadrature of the
      ! functional itself rather than of its kernel, and the two grids differ.
      call check(error, abs(scf%energy - PBE_STO3G_ENERGY) < TOL_GRID, &
                 "the PBE energy is not the one the reference matrix was taken at")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call dense_tda(mol, scf, ctx, .true., a, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the PBE TDA operator failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call compare_matrix(error, a, PBE_A, TOL_GRID, "PBE")
      if (allocated(error)) return

      values = eigenvalues_of(a, ok)
      call check(error, ok, "the dense diagonalisation of the PBE TDA matrix failed")
      if (allocated(error)) return
      call check(error, maxval(abs(values - PBE_A_EIG)) < TOL_GRID, &
                 "the PBE TDA spectrum disagrees with PySCF")
   end subroutine test_pbe_matrix

   subroutine water_fragment(fragment)
      !! The plan geometry as the bridge wants it: element numbers and Bohr
      type(physical_fragment_t), intent(out) :: fragment

      fragment%n_atoms = 3
      fragment%charge = 0
      fragment%multiplicity = 1
      fragment%nelec = 10
      fragment%n_caps = 0
      allocate (fragment%element_numbers(3), fragment%coordinates(3, 3))
      fragment%element_numbers = [8, 1, 1]
      fragment%coordinates = WATER_BOHR
   end subroutine water_fragment

   subroutine excited_run(basis, functional, n_states, method, spin, result)
      !! One whole calculation through the bridge, with an excited-state block
      !!
      !! Through `run_czt_hf` rather than the solver directly, because what
      !! Layer 2 adds is as much the wiring as the operator: a spectrum that
      !! never reaches `calculation_result_t` is not a feature.
      character(len=*), intent(in) :: basis, functional, method, spin
      integer, intent(in) :: n_states
      type(calculation_result_t), intent(out) :: result

      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment

      call water_fragment(fragment)
      settings%basis_set = basis
      settings%functional = functional
      settings%grid_level = 5
      ! The reference spectra come from a dense diagonalisation and carry no
      ! solver floor of their own, so the SCF and the Davidson both have to be
      ! converged past the tolerance the roots are compared at.
      settings%energy_tol = 1.0e-12_dp
      settings%grad_tol = 1.0e-9_dp
      settings%density_tol = 1.0e-9_dp
      settings%max_iter = 200
      settings%excited%enabled = n_states > 0
      settings%excited%n_states = n_states
      settings%excited%method = method
      settings%excited%spin = spin
      settings%excited%tolerance = 1.0e-8_dp
      settings%excited%max_iter = 200

      call run_czt_hf(settings, fragment, result)
   end subroutine excited_run

   subroutine compare_roots(error, result, reference, tol, what)
      !! Every root of the table, against what the run reported
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), intent(in) :: result
      real(dp), intent(in) :: reference(:), tol
      character(len=*), intent(in) :: what

      integer :: i

      call check(error,.not. result%has_error, "the "//what//" run failed: "// &
                 result%error%get_message())
      if (allocated(error)) return
      call check(error, result%has_excited_states, "the "//what//" run reported no "// &
                 "excited states")
      if (allocated(error)) return
      call check(error, allocated(result%excitation_energies), "the "//what// &
                 " run set has_excited_states without allocating any energies")
      if (allocated(error)) return
      call check(error, size(result%excitation_energies) == size(reference), &
                 "the "//what//" run converged a different number of roots than "// &
                 "were asked for")
      if (allocated(error)) return

      do i = 1, size(reference)
         call check(error, abs(result%excitation_energies(i) - reference(i)) < tol, &
                    "a "//what//" TDA excitation energy disagrees with the table")
         if (allocated(error)) return
      end do

      call check(error, allocated(result%state_spin), "the "//what//" run labelled "// &
                 "no spins")
      if (allocated(error)) return
      call check(error, all(result%state_spin == STATE_SPIN_SINGLET), &
                 "a singlet solve reported a root that is not labelled a singlet")
   end subroutine compare_roots

   subroutine test_davidson(error)
      !! The solver's roots are the dense matrix's, on the same system
      !!
      !! Both sides are this program's: the matrix from applying the operator
      !! to unit vectors, the roots from the Davidson the bridge runs. What
      !! this measures is the solver and the guess, with the operator held
      !! fixed -- the PySCF comparison above is what says the operator is
      !! right.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: ctx
      type(error_t) :: err
      type(calculation_result_t) :: result
      real(dp), allocatable :: a(:, :), values(:)
      logical :: ok
      integer :: i

      call water_sto3g(mol, scf, ctx, err)
      if (.not. err%has_error()) call dense_tda(mol, scf, ctx, .false., a, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the dense TDA matrix failed: "// &
                 err%get_message())
      if (allocated(error)) return
      values = eigenvalues_of(a, ok)
      call check(error, ok, "the dense diagonalisation failed")
      if (allocated(error)) return

      call excited_run("sto-3g", "", 5, "tda", "singlet", result)
      call compare_roots(error, result, values(1:5), TOL_DAVIDSON, "STO-3G")
      if (allocated(error)) return

      ! And against PySCF, which the dense matrix already matched: a solver
      ! agreeing with a wrong matrix would pass the comparison above.
      do i = 1, 5
         call check(error, abs(result%excitation_energies(i) - RHF_A_EIG(i)) < &
                    TOL_DAVIDSON, "a converged STO-3G root disagrees with PySCF")
         if (allocated(error)) return
      end do
   end subroutine test_davidson

   subroutine test_ccpvdz_rhf(error)
      !! The five Hartree-Fock singlets of the plan's table
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      call excited_run("cc-pvdz", "", N_CCPVDZ_STATES, "tda", "singlet", result)
      call check(error, abs(result%energy%scf - RHF_CCPVDZ_ENERGY) < 1.0e-8_dp, &
                 "the cc-pVDZ Hartree-Fock energy is not the table's")
      if (allocated(error)) return
      call compare_roots(error, result, RHF_CCPVDZ_TDA, TOL_CCPVDZ_HF, "cc-pVDZ RHF")
   end subroutine test_ccpvdz_rhf

   subroutine test_ccpvdz_pbe(error)
      !! The five PBE singlets
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "pbe", N_CCPVDZ_STATES, "tda", "singlet", result)
      call check(error, abs(result%energy%scf - PBE_CCPVDZ_ENERGY) < TOL_GRID, &
                 "the cc-pVDZ PBE energy is not the table's")
      if (allocated(error)) return
      call compare_roots(error, result, PBE_CCPVDZ_TDA, TOL_GRID, "cc-pVDZ PBE")
   end subroutine test_ccpvdz_pbe

   subroutine test_ccpvdz_b3lyp(error)
      !! The five B3LYP singlets, libxc 402 on both sides
      !!
      !! `hyb_gga_xc_b3lyp` is the VWN-RPA flavour, which is what PySCF's
      !! bare `B3LYP` resolves to as well. The other spelling in circulation
      !! moves every root here by a few times 1e-4, which is four orders above
      !! the bound, so a functional mix-up fails rather than passes loosely.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "b3lyp", N_CCPVDZ_STATES, "tda", "singlet", result)
      call check(error, abs(result%energy%scf - B3LYP_CCPVDZ_ENERGY) < TOL_GRID, &
                 "the cc-pVDZ B3LYP energy is not the table's")
      if (allocated(error)) return
      call compare_roots(error, result, B3LYP_CCPVDZ_TDA, TOL_GRID, "cc-pVDZ B3LYP")
   end subroutine test_ccpvdz_b3lyp

   subroutine test_ccpvdz_cam(error)
      !! The five CAM-B3LYP singlets, which need the attenuated exchange pass
      !!
      !! The one case here with two exchange builds per product: `exx_fraction`
      !! of the full-range matrix and `rs_k_lr` of the one against
      !! `erf(omega r)/r`. Dropping the second leaves a converged spectrum of a
      !! plain 0.19-hybrid, which is the fault `test_mqc_czt_cphf` found in the
      !! static polarizability.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "cam-b3lyp", N_CCPVDZ_STATES, "tda", "singlet", result)
      call check(error, abs(result%energy%scf - CAM_CCPVDZ_ENERGY) < TOL_GRID, &
                 "the cc-pVDZ CAM-B3LYP energy is not the table's")
      if (allocated(error)) return
      call compare_roots(error, result, CAM_CCPVDZ_TDA, TOL_GRID, "cc-pVDZ CAM-B3LYP")
   end subroutine test_ccpvdz_cam

   subroutine test_no_states(error)
      !! `n_states = 0` leaves the calculation exactly as it was
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      call excited_run("sto-3g", "", 0, "tda", "singlet", result)
      call check(error,.not. result%has_error, "a deck asking for no excited states "// &
                 "failed: "//result%error%get_message())
      if (allocated(error)) return
      call check(error, result%has_energy, "a ground-state run returned no energy")
      if (allocated(error)) return
      call check(error, abs(result%energy%scf - RHF_STO3G_ENERGY) < TOL_EXACT, &
                 "the ground-state energy moved when the excited-state block was off")
      if (allocated(error)) return
      call check(error,.not. result%has_excited_states, "a run that asked for no "// &
                 "roots reported excited states")
      if (allocated(error)) return
      call check(error,.not. allocated(result%excitation_energies), &
                 "a run that asked for no roots allocated excitation energies")
   end subroutine test_no_states

   subroutine test_later_layers(error)
      !! What Layer 2 does not do is refused by name, not approximated
      !!
      !! A TDA number handed back for an RPA deck is the failure worth
      !! guarding against: it converges, it is of the right magnitude, and
      !! nothing in the output says which problem was solved.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      call excited_run("sto-3g", "", 2, "rpa", "singlet", result)
      call check(error, result%has_error, "an RPA deck was answered rather than "// &
                 "refused")
      if (allocated(error)) return
      call check(error,.not. result%has_excited_states, "a refused RPA deck still "// &
                 "reported excited states")
      if (allocated(error)) return

      call excited_run("sto-3g", "", 2, "tda", "triplet", result)
      call check(error, result%has_error, "a triplet deck was answered rather than "// &
                 "refused")
      if (allocated(error)) return
      call check(error,.not. result%has_excited_states, "a refused triplet deck "// &
                 "still reported excited states")
   end subroutine test_later_layers

end module test_mqc_czt_tddft

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_tddft, only: collect_mqc_czt_tddft_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_tddft", collect_mqc_czt_tddft_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
