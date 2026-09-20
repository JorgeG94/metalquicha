!! Singlet excitation energies, Tamm-Dancoff and full RPA, against PySCF
module test_mqc_czt_tddft
   !! The Layer 2 and Layer 4 gates of `TDDFT_PLAN.md`: the TDA operator `A`
   !! itself, the paired `(A+B)`/`(A-B)` operator behind it, and the
   !! excitation energies each solver finds.
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
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_tddft, only: tda_operator_t, build_tda_operator, tda_dense_matrix, &
                            rpa_operator_t, build_rpa_operator, rpa_dense_matrices, &
                            singlet_excitations
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

   integer, parameter :: N_F2_STATES = 3
      !! Roots asked for in the F2 case, which is the request that used to
      !! come back with the wrong third one.

   !! F2 at its equilibrium bond length, **in Bohr**, along z.
   !!
   !! 1.4119 Angstrom is the experimental `r_e`; 2.668 Bohr is that to the
   !! four figures the reference was taken at, and PySCF was given the same
   !! number as Bohr, so no conversion constant enters the comparison.
   real(dp), parameter :: F2_BOHR(3, 2) = reshape([ &
                                                  0.0_dp, 0.0_dp, 0.0_dp, &
                                                  0.0_dp, 0.0_dp, 2.668_dp], [3, 2])

   !! F2/cc-pVDZ, restricted Hartree-Fock, the three lowest singlet TDA roots.
   !!
   !! The spectrum this molecule is here for:
   !!
   !!     0.183610964  0.183610964  0.332986263  0.332986263  0.556013580
   !!
   !! two degenerate pairs and then a single root. What makes it the case to
   !! pin is where those roots sit in the *gaps* the Davidson guess is picked
   !! on: the three lowest are 0.759, 0.759 and 0.841 hartree, and the 0.841
   !! one carries the fifth root, while the pair carrying the third and
   !! fourth is at 0.903 -- outside a three-vector guess, and 0.06 hartree
   !! past anything a degeneracy window closes. A guess of one unit vector
   !! per root therefore converges roots one, two and *five* and reports the
   !! fifth as the third: every one of them a true eigenpair, none of them
   !! flagged. Asked for five roots the same solver finds all five, which is
   !! what says the starting space was too narrow rather than the solver
   !! broken.
   !!
   !! Taken the same way as the water tables: PySCF 2.14 fed this
   !! repository's own cc-pVDZ JSON through `bse_to_pyscf`, `conv_tol =
   !! 1e-15`, the whole 171 by 171 `A` probed out of `TDA.gen_vind` on the
   !! unit vectors and diagonalised densely, so no iterative tolerance is in
   !! the reference.
   real(dp), parameter :: F2_CCPVDZ_ENERGY = -198.685678500661_dp
   real(dp), parameter :: F2_CCPVDZ_TDA(N_F2_STATES) = [ &
                          0.183610964095_dp, 0.183610964095_dp, 0.332986263273_dp]

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
   !! under the same grid. On cc-pVDZ, where the roots come through the
   !! solver rather than off a dense matrix, the worst of the five is 1.9e-10
   !! for PBE, 1.6e-9 for B3LYP and 6.4e-10 for CAM-B3LYP -- B3LYP being the
   !! one to watch, at a sixtieth of this bound. The RPA column lands in the
   !! same places, 1.8e-10, 1.5e-9 and 6.3e-10, so the paired solver adds
   !! nothing to what the quadrature already costs. 1e-7 clears the worst of them
   !! by that margin and leaves the compiler spread room, which the
   !! double-hybrid Hessian shows can reach 2e-9 on a quadrature of this kind.
   !! Tightening it to the measured numbers would be re-recording a pin to
   !! make it pass on one compiler.
   real(dp), parameter :: TOL_GRID = 1.0e-7_dp

   !! The Davidson's own floor, against a dense diagonalisation of the same
   !! operator. Roots are accepted on a residual of 1e-8 here, and the
   !! eigenvalue error near an eigenvector is second order in the vector
   !! error, so this is a bound the solver clears by construction rather than
   !! one fitted to it -- measured 1.5e-11 on the worst of the five STO-3G
   !! roots, against both the dense spectrum and PySCF's.
   real(dp), parameter :: TOL_DAVIDSON = 1.0e-9_dp

   !! The cc-pVDZ Hartree-Fock roots. Looser than `TOL_EXACT` because these
   !! come out of the Davidson rather than off a dense diagonalisation, so the
   !! solver's own floor is in them as well as the integrals'. Measured
   !! 3.1e-12 on the worst of the five for TDA and 9.2e-11 for RPA, whose
   !! reduction amplifies the same disagreement the way it does at STO-3G.
   real(dp), parameter :: TOL_CCPVDZ_HF = 1.0e-8_dp

   !! ------------------------------------------------------------------
   !! Layer 4: the full RPA
   !! ------------------------------------------------------------------

   !! Every root of the H2O/STO-3G Hartree-Fock TDHF problem, in Hartree.
   !!
   !! From `pyscf.tdscf.rhf.get_ab` on the same reference the matrices above
   !! were taken from, reduced as `(A-B)^{1/2}(A+B)(A-B)^{1/2}` and
   !! diagonalised whole, so there is no solver tolerance in them. The
   !! non-Hermitian `2n` problem was diagonalised as well and agrees to
   !! 4.6e-14, and PySCF's own iterative `TDHF` to 1.0e-13 -- three routes to
   !! the same ten numbers, which is what makes them a reference rather than
   !! one code's output.
   !!
   !! Every one of them sits below its Tamm-Dancoff partner in `RHF_A_EIG`,
   !! which is the whole physical content of the approximation and is worth
   !! seeing in the table: 0.4851 against 0.4835, 0.5558 against 0.5553.
   real(dp), parameter :: RHF_STO3G_RPA(N_OV) = [ &
                          0.483544026053_dp, 0.555275192922_dp, 0.613540106713_dp, &
                          0.702292095334_dp, 0.806409059989_dp, 1.045220958387_dp, &
                          1.462314417289_dp, 1.508625576304_dp, 20.107306379835_dp, &
                          20.157271711632_dp]

   !! H2O/cc-pVDZ, the five singlet RPA roots, in Hartree, for the four
   !! functionals. Regenerated to twelve decimals the same way -- our own
   !! basis JSON through `bse_to_pyscf`, the geometry handed over in Bohr --
   !! and they reproduce the plan's ten-decimal table to 1.3e-9 at worst,
   !! which is where a table written to ten decimals stops saying anything.
   real(dp), parameter :: RHF_CCPVDZ_RPA(N_CCPVDZ_STATES) = [ &
                          0.336535689992_dp, 0.401350367586_dp, 0.432987560701_dp, &
                          0.497799802347_dp, 0.551224557096_dp]
   real(dp), parameter :: PBE_CCPVDZ_RPA(N_CCPVDZ_STATES) = [ &
                          0.269685669408_dp, 0.339275987762_dp, 0.354437859012_dp, &
                          0.428784470826_dp, 0.509472350901_dp]
   real(dp), parameter :: B3LYP_CCPVDZ_RPA(N_CCPVDZ_STATES) = [ &
                          0.279639965109_dp, 0.348189899392_dp, 0.365848830007_dp, &
                          0.438303579652_dp, 0.514775708061_dp]
   real(dp), parameter :: CAM_CCPVDZ_RPA(N_CCPVDZ_STATES) = [ &
                          0.282337416507_dp, 0.353035871122_dp, 0.368974575576_dp, &
                          0.443809832697_dp, 0.515659298412_dp]

   !! Hydrogen at 3.0 Angstrom in 6-31G, written out in Bohr for the reason
   !! the water geometry is.
   !!
   !! A stretched closed-shell H2 is the standard place to look for a broken
   !! reference, and the singlet channel is not where it breaks: `(A-B)` here
   !! has eigenvalues 0.0305, 0.9295 and 1.0355, all positive, and PySCF
   !! converges the spectrum without complaint. So this case gates the
   !! agreement, not the refusal -- the instability path is exercised in
   !! `test_mqc_czt_rpa_solver`, where an indefinite difference can be
   !! constructed rather than hoped for. What the case is still worth: the
   !! smallest of those eigenvalues is thirty times below the others, which
   !! is the near-singular `(A-B)` its square root has to survive.
   real(dp), parameter :: H2_BOHR(3, 2) = reshape([ &
                                                  0.0_dp, 0.0_dp, 0.0_dp, &
                                                  0.0_dp, 0.0_dp, 5.66917837637348399_dp], [3, 2])
   real(dp), parameter :: H2_631G_ENERGY = -0.815591795493_dp
   real(dp), parameter :: H2_631G_RPA(3) = [ &
                          0.112116666867_dp, 0.999987235658_dp, 1.099106035658_dp]

   !! The dense reduction against PySCF's, and why this is not `TOL_EXACT`.
   !!
   !! Hartree-Fock carries no quadrature, so the two codes' `(A+B)` and
   !! `(A-B)` agree to 7.7e-13 element by element -- and their *reduction*
   !! agrees three decimals worse than that. `(A-B)^{1/2}(A+B)(A-B)^{1/2}`
   !! has a norm near 800 on this system, so a part in 1e12 of the matrices
   !! is a part in 1e11 of `w^2`, and `w = sqrt(w^2)` halves nothing that
   !! matters at these magnitudes. The measured per-root disagreement is
   !! 2.3e-11, 2.5e-11, 1.6e-11, 1.9e-11, 3.8e-12, 5.0e-12, 3.1e-12,
   !! 1.2e-11, 7.5e-11 and 6.7e-11 -- the two worst being the 20-hartree core
   !! excitations, which carry the norm. A bound of 1e-10 would pass with a
   !! margin of 1.3, which is a pin fitted to one compiler rather than a
   !! gate; this one clears the worst by thirteen and still sits eight orders
   !! under the smallest fault it has to catch, a sign error or a dropped
   !! term in either half moving roots by 1e-2 and up.
   real(dp), parameter :: TOL_RPA_DENSE = 1.0e-9_dp

   !! The paired solver against a dense reduction of the operator it was
   !! given. Measured 7.0e-14 on the worst of the five STO-3G roots against
   !! the dense spectrum and 2.5e-11 against PySCF's, and 1.9e-12 on
   !! stretched H2. The bound is the Davidson's, for the same reason -- roots
   !! are accepted on a residual and the eigenvalue error near an
   !! eigenvector is second order in the vector error.
   real(dp), parameter :: TOL_RPA_SOLVER = 1.0e-9_dp

   !! `sum(X^2) - sum(Y^2)` against the half the routine documents. An
   !! algebraic identity the solver imposes by division rather than a
   !! converged quantity, so it holds to round-off however far the roots got:
   !! measured 2.2e-16.
   real(dp), parameter :: TOL_PAIRED_NORM = 1.0e-10_dp

   !! The Casida reduction against the paired solver, on the same PBE
   !! reference. Two solvers, two operators and one spectrum; measured
   !! 9.0e-14 on the worst of five.
   real(dp), parameter :: TOL_CASIDA = 1.0e-9_dp

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
                  new_unittest("three_roots_of_f2_are_the_lowest_three", test_f2_lowest_three), &
                  new_unittest("no_states_asked_for_means_no_spectrum", test_no_states), &
                  new_unittest("triplets_are_refused_by_name", test_later_layers), &
                  new_unittest("the_rhf_rpa_spectrum_matches_pyscf", test_rhf_rpa_matrix), &
                  new_unittest("the_paired_solver_matches_the_dense_reduction", &
                               test_rpa_solver), &
                  new_unittest("cc_pvdz_rhf_rpa_matches_the_table", test_ccpvdz_rhf_rpa), &
                  new_unittest("cc_pvdz_pbe_rpa_matches_the_table", test_ccpvdz_pbe_rpa), &
                  new_unittest("cc_pvdz_b3lyp_rpa_matches_the_table", &
                               test_ccpvdz_b3lyp_rpa), &
                  new_unittest("cc_pvdz_cam_b3lyp_rpa_matches_the_table", &
                               test_ccpvdz_cam_rpa), &
                  new_unittest("the_casida_reduction_matches_the_paired_solver", &
                               test_casida), &
                  new_unittest("stretched_h2_rpa_matches_pyscf", test_stretched_h2), &
                  new_unittest("an_unreachable_tolerance_stops_and_says_so", &
                               test_unreachable_tolerance) &
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
      !! The explicit TDA matrix, through the operator the solver uses
      !!
      !! `tda_dense_matrix` applies the shipped operator to the ten unit
      !! vectors, which is the construction PySCF's side of the reference
      !! used as well -- so the two are compared as matrices and not as two
      !! summaries of one. Going through the real operator rather than
      !! assembling `(1/2)[(A+B)+(A-B)]` here is the point: a test that redid
      !! the half sum itself would agree with a shipped operator that had
      !! dropped a term, since both halves come from the same routine.
      !!
      !! `mol` and `ctx` are `target` because the operator keeps pointers to
      !! them; they stay valid for as long as this call runs, which is longer
      !! than the operator is used for.
      type(czt_molecule_t), intent(in), target :: mol
      type(rhf_result_t), intent(in) :: scf
      type(xc_context_t), intent(inout), target :: ctx
      logical, intent(in) :: kohn_sham
      real(dp), allocatable, intent(out) :: a(:, :)
      type(error_t), intent(inout) :: err

      type(tda_operator_t) :: operator

      if (err%has_error()) return

      if (kohn_sham) then
         call build_tda_operator(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, operator, err, xc=ctx, &
                                 reference=scf%density)
      else
         call build_tda_operator(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, operator, err)
      end if
      if (err%has_error()) return

      call tda_dense_matrix(operator, a, err)
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

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
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

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
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

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
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

   subroutine difluorine_fragment(fragment)
      !! F2 at `F2_BOHR`, closed shell, eighteen electrons
      type(physical_fragment_t), intent(out) :: fragment

      fragment%n_atoms = 2
      fragment%charge = 0
      fragment%multiplicity = 1
      fragment%nelec = 18
      fragment%n_caps = 0
      allocate (fragment%element_numbers(2), fragment%coordinates(3, 2))
      fragment%element_numbers = [9, 9]
      fragment%coordinates = F2_BOHR
   end subroutine difluorine_fragment

   subroutine excited_run(basis, functional, n_states, method, spin, result, &
                          tolerance, molecule)
      !! One whole calculation through the bridge, with an excited-state block
      !!
      !! Through `run_czt_hf` rather than the solver directly, because what
      !! Layer 2 adds is as much the wiring as the operator: a spectrum that
      !! never reaches `calculation_result_t` is not a feature.
      character(len=*), intent(in) :: basis, functional, method, spin
      integer, intent(in) :: n_states
      type(calculation_result_t), intent(out) :: result
      real(dp), intent(in), optional :: tolerance
         !! What a root is accepted at. Absent is 1e-8, which every case here
         !! but the deliberately unreachable one is compared at.
      type(physical_fragment_t), intent(in), optional :: molecule
         !! What to run it on. Absent is the water every other case here uses.

      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment

      if (present(molecule)) then
         fragment = molecule
      else
         call water_fragment(fragment)
      end if
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
      if (present(tolerance)) settings%excited%tolerance = tolerance
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
                    "a "//what//" excitation energy disagrees with the table")
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

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
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

   subroutine test_f2_lowest_three(error)
      !! Three roots means the lowest three, not three true roots of any rank
      !!
      !! The failure this pins is not a wrong number: every root the old
      !! guess returned was an eigenvalue of the right matrix, converged to
      !! the tolerance asked for, and the third one was the fifth of the
      !! spectrum. Nothing in the run said so, which is why it is gated here
      !! and not left to a user to notice.
      !!
      !! Both halves are checked: that the third root is 0.3330 and not
      !! 0.5560, and that the first two are still the degenerate pair -- a
      !! guess that lost the *pair* would move those instead and is a
      !! different fault with the same cause.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result
      type(physical_fragment_t) :: fragment

      call difluorine_fragment(fragment)
      call excited_run("cc-pvdz", "", N_F2_STATES, "tda", "singlet", result, &
                       molecule=fragment)
      call check(error, abs(result%energy%scf - F2_CCPVDZ_ENERGY) < 1.0e-8_dp, &
                 "the F2/cc-pVDZ Hartree-Fock energy is not the one the reference "// &
                 "spectrum was taken at")
      if (allocated(error)) return
      call compare_roots(error, result, F2_CCPVDZ_TDA, TOL_CCPVDZ_HF, "F2/cc-pVDZ")
      if (allocated(error)) return

      ! Said separately, because the bound above would also be cleared by a
      ! run that reported two roots and stopped.
      call check(error, size(result%excitation_energies) == N_F2_STATES, &
                 "the F2 run did not report three roots")
      if (allocated(error)) return
      call check(error, result%excitation_energies(3) < 0.4_dp, &
                 "the third F2 root is the fifth of the spectrum: the Davidson "// &
                 "guess does not span the degenerate pair above the first two")
   end subroutine test_f2_lowest_three

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
      !! What is still not implemented is refused by name, not approximated
      !!
      !! A singlet number handed back for a triplet deck is the failure worth
      !! guarding against: it converges, it is of the right magnitude, and
      !! nothing in the output says which problem was solved. The RPA half of
      !! this case is gone because Layer 4 answers an RPA deck; what it
      !! answers with is `test_ccpvdz_rhf_rpa` and the three beside it.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      call excited_run("sto-3g", "", 2, "tda", "triplet", result)
      call check(error, result%has_error, "a triplet deck was answered rather than "// &
                 "refused")
      if (allocated(error)) return
      call check(error,.not. result%has_excited_states, "a refused triplet deck "// &
                 "still reported excited states")
      if (allocated(error)) return

      ! A method string that is neither is still refused, and named.
      call excited_run("sto-3g", "", 2, "cis", "singlet", result)
      call check(error, result%has_error, "an unknown excited-state method was "// &
                 "resolved to something rather than refused")
   end subroutine test_later_layers

   subroutine dense_rpa(mol, scf, ctx, kohn_sham, aplus, aminus, err)
      !! The explicit `(A+B)` and `(A-B)`, through the operator the solver uses
      !!
      !! The same argument as `dense_tda`: probing the shipped operator with
      !! unit vectors is what makes the comparison a comparison. Assembling
      !! the two halves here out of `build_hessian` would check one
      !! construction of the physics against another construction of the same
      !! physics, and both could be wrong together.
      type(czt_molecule_t), intent(in), target :: mol
      type(rhf_result_t), intent(in) :: scf
      type(xc_context_t), intent(inout), target :: ctx
      logical, intent(in) :: kohn_sham
      real(dp), allocatable, intent(out) :: aplus(:, :), aminus(:, :)
      type(error_t), intent(inout) :: err

      type(rpa_operator_t) :: operator

      if (err%has_error()) return

      if (kohn_sham) then
         call build_rpa_operator(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, operator, err, xc=ctx, &
                                 reference=scf%density)
      else
         call build_rpa_operator(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, operator, err)
      end if
      if (err%has_error()) return

      call rpa_dense_matrices(operator, aplus, aminus, err)
   end subroutine dense_rpa

   function paired_spectrum(aplus, aminus, ok) result(values)
      !! Every `w` of `[A B; B A]`, from the explicit reduction
      !!
      !! `(A-B)^{1/2}(A+B)(A-B)^{1/2}`, built with LAPACK and diagonalised
      !! whole -- no iteration, so this is the answer the iterative solver is
      !! measured against rather than a second approximation to it. It is
      !! also the arithmetic PySCF's side of `RHF_STO3G_RPA` performed, which
      !! is why the two can be compared to 1e-10.
      real(dp), intent(in) :: aplus(:, :), aminus(:, :)
      logical, intent(out) :: ok
      real(dp), allocatable :: values(:)

      real(dp), allocatable :: vectors(:, :), w(:), half(:, :), scaled(:, :)
      real(dp), allocatable :: work(:, :), reduced(:, :)
      integer :: n, info, k

      n = size(aplus, 1)
      ok = .false.
      allocate (values(n), w(n))
      vectors = 0.5_dp*(aminus + transpose(aminus))
      call pic_syev(vectors, w, jobz="V", uplo="U", info=info)
      if (info /= 0) return
      if (minval(w) <= 0.0_dp) return

      allocate (scaled(n, n), half(n, n), work(n, n), reduced(n, n))
      do k = 1, n
         scaled(:, k) = vectors(:, k)*sqrt(w(k))
      end do
      call pic_gemm(scaled, vectors, half, transb="T")
      call pic_gemm(0.5_dp*(aplus + transpose(aplus)), half, work)
      call pic_gemm(half, work, reduced)
      reduced = 0.5_dp*(reduced + transpose(reduced))
      call pic_syev(reduced, values, jobz="N", uplo="U", info=info)
      if (info /= 0) return
      values = sqrt(values)
      ok = .true.
   end function paired_spectrum

   subroutine test_rhf_rpa_matrix(error)
      !! Every RPA root of H2O/STO-3G, from the two explicit halves
      !!
      !! The matrix-level gate of Layer 4. `RHF_A` already says `(A+B)` and
      !! `(A-B)` are right in the combination the Tamm-Dancoff operator reads
      !! them in -- their half sum -- and that combination cannot see an
      !! error that moves the two halves oppositely. The reduction below
      !! reads them separately, so it can.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
      type(error_t) :: err
      real(dp), allocatable :: aplus(:, :), aminus(:, :), values(:)
      logical :: ok

      call water_sto3g(mol, scf, ctx, err)
      if (.not. err%has_error()) call dense_rpa(mol, scf, ctx, .false., aplus, aminus, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the Hartree-Fock RPA matrices failed: "// &
                 err%get_message())
      if (allocated(error)) return

      ! The two halves average to the matrix Layer 2 already pinned, which is
      ! what says the split itself is right before its spectrum is read.
      call compare_matrix(error, 0.5_dp*(aplus + aminus), RHF_A, TOL_EXACT, &
                          "Hartree-Fock half-sum")
      if (allocated(error)) return

      values = paired_spectrum(aplus, aminus, ok)
      call check(error, ok, "(A-B) of a converged closed shell was not positive "// &
                 "definite, so the reduction could not be formed")
      if (allocated(error)) return
      call check(error, maxval(abs(values - RHF_STO3G_RPA)) < TOL_RPA_DENSE, &
                 "the Hartree-Fock RPA spectrum disagrees with PySCF")
      if (allocated(error)) return

      ! Every RPA root below its Tamm-Dancoff partner: the physical content
      ! of dropping `B`, and a sanity check no tolerance would catch if the
      ! two spectra had been transcribed into each other's places.
      call check(error, all(RHF_STO3G_RPA < RHF_A_EIG), "an RPA root came out above "// &
                 "its Tamm-Dancoff partner, which the de-excitation coupling "// &
                 "cannot do")
   end subroutine test_rhf_rpa_matrix

   subroutine test_rpa_solver(error)
      !! The paired solver's roots are the dense reduction's, and its
      !! amplitudes carry the normalisation the routine documents
      !!
      !! Both sides are this program's: the matrices from applying the
      !! operator to unit vectors, the roots from the solver. What this
      !! measures is the solver, with the operator held fixed -- the PySCF
      !! comparison above is what says the operator is right.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
      type(error_t) :: err
      real(dp), allocatable :: aplus(:, :), aminus(:, :), dense(:, :)
      real(dp), allocatable :: values(:), omega(:), x(:, :), y(:, :)
      real(dp) :: weight
      logical :: ok
      integer :: k

      call water_sto3g(mol, scf, ctx, err)
      if (.not. err%has_error()) call dense_rpa(mol, scf, ctx, .false., aplus, aminus, err)
      if (err%has_error()) then
         call mol%destroy()
         call check(error, .false., "the Hartree-Fock RPA matrices failed: "// &
                    err%get_message())
         return
      end if
      values = paired_spectrum(aplus, aminus, ok)

      call singlet_excitations(mol, scf%orbitals, scf%orbital_energies, &
                               scf%n_occupied, 5, "rpa", omega, x, y, err, &
                               tolerance=1.0e-10_dp, max_iter=100)
      call mol%destroy()
      call check(error, ok .and. .not. err%has_error(), "the paired solve failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, size(omega) == 5, "the paired solve returned a different "// &
                 "number of roots than were asked for")
      if (allocated(error)) return

      call check(error, maxval(abs(omega - values(1:5))) < TOL_RPA_SOLVER, &
                 "a paired root disagrees with the dense reduction of the same "// &
                 "two matrices")
      if (allocated(error)) return
      call check(error, maxval(abs(omega - RHF_STO3G_RPA(1:5))) < TOL_RPA_SOLVER, &
                 "a converged paired root disagrees with PySCF")
      if (allocated(error)) return

      ! The normalisation, which is documented rather than derivable and is
      ! the one thing Layer 5 has to be able to rely on.
      do k = 1, 5
         weight = dot_product(x(:, k), x(:, k)) - dot_product(y(:, k), y(:, k))
         call check(error, abs(weight - 0.5_dp) < TOL_PAIRED_NORM, &
                    "a paired root's amplitudes are not at the |X|^2-|Y|^2 = 1/2 "// &
                    "the routine documents")
         if (allocated(error)) return
      end do

      ! `Y` is not zero here, which is the only difference between this and
      ! the Tamm-Dancoff answer: a solver that quietly dropped the
      ! de-excitation block would pass every energy check above on a system
      ! this small and fail here.
      dense = y
      call check(error, maxval(abs(dense)) > 1.0e-3_dp, "the paired solve returned "// &
                 "de-excitation amplitudes of zero, which is the Tamm-Dancoff "// &
                 "answer and not this one")
   end subroutine test_rpa_solver

   subroutine test_ccpvdz_rhf_rpa(error)
      !! The five Hartree-Fock RPA roots of the plan's table
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      call excited_run("cc-pvdz", "", N_CCPVDZ_STATES, "rpa", "singlet", result)
      call check(error, abs(result%energy%scf - RHF_CCPVDZ_ENERGY) < 1.0e-8_dp, &
                 "the cc-pVDZ Hartree-Fock energy is not the table's")
      if (allocated(error)) return
      call compare_roots(error, result, RHF_CCPVDZ_RPA, TOL_CCPVDZ_HF, "cc-pVDZ RHF RPA")
   end subroutine test_ccpvdz_rhf_rpa

   subroutine test_ccpvdz_pbe_rpa(error)
      !! The five PBE RPA roots, where `(A-B)` is the bare diagonal
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "pbe", N_CCPVDZ_STATES, "rpa", "singlet", result)
      call check(error, abs(result%energy%scf - PBE_CCPVDZ_ENERGY) < TOL_GRID, &
                 "the cc-pVDZ PBE energy is not the table's")
      if (allocated(error)) return
      call compare_roots(error, result, PBE_CCPVDZ_RPA, TOL_GRID, "cc-pVDZ PBE RPA")
   end subroutine test_ccpvdz_pbe_rpa

   subroutine test_ccpvdz_b3lyp_rpa(error)
      !! The five B3LYP RPA roots, libxc 402 on both sides
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "b3lyp", N_CCPVDZ_STATES, "rpa", "singlet", result)
      call check(error, abs(result%energy%scf - B3LYP_CCPVDZ_ENERGY) < TOL_GRID, &
                 "the cc-pVDZ B3LYP energy is not the table's")
      if (allocated(error)) return
      call compare_roots(error, result, B3LYP_CCPVDZ_RPA, TOL_GRID, "cc-pVDZ B3LYP RPA")
   end subroutine test_ccpvdz_b3lyp_rpa

   subroutine test_ccpvdz_cam_rpa(error)
      !! The five CAM-B3LYP RPA roots, which need the attenuated exchange
      !! pass in both halves
      !!
      !! `(A-B)` is exchange only, so this is the one case here where the
      !! long-range pass is the whole of a product rather than a correction
      !! to one. Dropping it moves the difference operator wholesale, which
      !! moves `w` in a direction no tolerance hides.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "cam-b3lyp", N_CCPVDZ_STATES, "rpa", "singlet", result)
      call check(error, abs(result%energy%scf - CAM_CCPVDZ_ENERGY) < TOL_GRID, &
                 "the cc-pVDZ CAM-B3LYP energy is not the table's")
      if (allocated(error)) return
      call compare_roots(error, result, CAM_CCPVDZ_RPA, TOL_GRID, "cc-pVDZ CAM-B3LYP RPA")
   end subroutine test_ccpvdz_cam_rpa

   subroutine test_casida(error)
      !! The Casida reduction and the paired solver agree on a pure functional
      !!
      !! Two solvers on two different operators -- a Hermitian Davidson on
      !! `dEps^{1/2}(A+B)dEps^{1/2}`, and the paired subspace method on
      !! `(A+B)` and `(A-B)` kept apart -- reaching one spectrum. They share
      !! the `(A+B)` product and nothing else, so an error in the paired
      !! solver's square root, its biorthonormalisation or its residuals
      !! shows up here and an error in the shared product does not. PBE, so
      !! that `(A-B)` really is the diagonal the reduction assumes.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
      type(error_t) :: err
      real(dp), allocatable :: paired(:), casida(:), x(:, :), y(:, :)

      if (.not. xc_available()) return

      call water_sto3g(mol, scf, ctx, err, functional="pbe")
      if (err%has_error()) then
         call mol%destroy()
         call check(error, .false., "the PBE reference failed: "//err%get_message())
         return
      end if

      call singlet_excitations(mol, scf%orbitals, scf%orbital_energies, &
                               scf%n_occupied, 5, "rpa", paired, x, y, err, &
                               xc=ctx, reference=scf%density, tolerance=1.0e-10_dp, &
                               max_iter=100)
      if (.not. err%has_error()) &
         call singlet_excitations(mol, scf%orbitals, scf%orbital_energies, &
                                  scf%n_occupied, 5, "casida", casida, x, y, err, &
                                  xc=ctx, reference=scf%density, &
                                  tolerance=1.0e-12_dp, max_iter=200)
      call mol%destroy()
      call check(error,.not. err%has_error(), "a PBE excitation solve failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, size(paired) == 5 .and. size(casida) == 5, &
                 "the two routes converged different numbers of roots")
      if (allocated(error)) return
      call check(error, maxval(abs(paired - casida)) < TOL_CASIDA, &
                 "the Casida reduction and the paired solver disagree on the PBE "// &
                 "spectrum")
   end subroutine test_casida

   subroutine test_stretched_h2(error)
      !! A near-singular `(A-B)`: stretched H2, against PySCF
      !!
      !! The lowest eigenvalue of `(A-B)` here is 0.0305 against 1.0 for the
      !! others, which is where a subspace square root is most likely to lose
      !! digits. PySCF converges this without complaint and so must this; the
      !! refusal path lives in `test_mqc_czt_rpa_solver`, where an indefinite
      !! difference can be built rather than waited for.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: omega(:), x(:, :), y(:, :)

      call build_czt_molecule([1, 1], ["H ", "H "], H2_BOHR, "6-31g", mol, err)
      if (.not. err%has_error()) &
         call run_czt_rhf(mol, 2, 200, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err, &
                          in_core=.true., grad_tol=1.0e-12_dp)
      if (err%has_error()) then
         call mol%destroy()
         call check(error, .false., "the stretched H2 reference failed: "// &
                    err%get_message())
         return
      end if
      call check(error, abs(scf%energy - H2_631G_ENERGY) < TOL_EXACT, &
                 "the stretched H2 energy is not the one the reference was taken at")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call singlet_excitations(mol, scf%orbitals, scf%orbital_energies, &
                               scf%n_occupied, 3, "rpa", omega, x, y, err, &
                               tolerance=1.0e-10_dp, max_iter=100)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the stretched H2 paired solve "// &
                 "failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, size(omega) == 3, "stretched H2 in 6-31G has three single "// &
                 "excitations and the solve did not return three")
      if (allocated(error)) return
      call check(error, maxval(abs(omega - H2_631G_RPA)) < TOL_RPA_SOLVER, &
                 "a stretched H2 RPA root disagrees with PySCF")
   end subroutine test_stretched_h2

   subroutine test_unreachable_tolerance(error)
      !! A tolerance nothing can reach comes back with a message, not a hang
      !!
      !! The failure this guards against is a solve that spends its whole
      !! iteration budget adding vectors that project to nothing, and then
      !! reports a spectrum with no indication of how good it is. Ten
      !! occupied-virtual rotations and a tolerance below round-off is the
      !! cheapest way to reach that state deliberately: the subspace becomes
      !! the whole space, every residual direction is already in it, and
      !! there is nothing left to add.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      call excited_run("sto-3g", "", 3, "rpa", "singlet", result, tolerance=1.0e-30_dp)
      call check(error, result%has_error, "a tolerance below round-off was reported "// &
                 "as met")
      if (allocated(error)) return
      call check(error,.not. result%has_excited_states, "a solve that could not "// &
                 "reach its tolerance still reported excited states")
      if (allocated(error)) return
      call check(error, index(result%error%get_message(), "residual") > 0, &
                 "a stalled solve did not report the residual it reached: "// &
                 result%error%get_message())
   end subroutine test_unreachable_tolerance

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
