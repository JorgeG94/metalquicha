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
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf
   use mqc_czt_tddft, only: tda_operator_t, build_tda_operator, tda_dense_matrix, &
                            rpa_operator_t, build_rpa_operator, rpa_dense_matrices, &
                            response_excitations, tda_operator_uhf_t, &
                            build_tda_operator_uhf, tda_dense_matrix_uhf, &
                            rpa_operator_uhf_t, build_rpa_operator_uhf, &
                            rpa_dense_matrices_uhf, response_excitations_uhf
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_czt_bridge, only: run_czt_hf
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t, STATE_SPIN_SINGLET, &
                               STATE_SPIN_TRIPLET, STATE_SPIN_UNRESTRICTED
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

   !! ------------------------------------------------------------------
   !! Layer 3: triplets
   !! ------------------------------------------------------------------

   !! The **triplet** TDA matrix of H2O/STO-3G, and the same for PBE.
   !!
   !! Regenerated from PySCF 2.14 through `bse_to_pyscf` exactly as the singlet
   !! matrices above were -- `td = tdscf.TDA(mf); td.singlet = False`, then
   !! `td.gen_vind` on the ten unit vectors. `TDDFT_PLAN.md`'s own triplet
   !! eigenvalues are 2.3e-8 away from these, which is the size of the
   !! internal-basis-table error Layer 2 found and not a disagreement about
   !! the physics; the plan's numbers are not usable at this bound. PySCF's
   !! own SCF is driven to `conv_tol = 1e-15` for these, not the 1e-12 the
   !! plan records: the two core-excited elements near 20 hartree move by
   !! 7e-11 between 1e-13 and 1e-15, which is most of `TOL_EXACT`.
   !!
   !! What separates these from the singlet matrices is the whole of Layer 3:
   !! no Coulomb term at all, and `(f_aa - f_ab)/2` in place of the singlet
   !! kernel. Either one left as it was moves the diagonal by tenths of a
   !! Hartree, so this comparison has no way to pass on half the change.
   real(dp), parameter :: RHF_TRIPLET_A(N_OV, N_OV) = reshape([ &
                                           20.044631057292_dp, 0.000000000000_dp, 0.006981609406_dp, -0.000000000000_dp, &
                                            0.000000000000_dp, -0.009192048725_dp, 0.021190548711_dp, 0.000000000000_dp, &
                                          -0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, 20.114461669827_dp, &
                                          -0.000000000000_dp, 0.009377471331_dp, -0.009192048725_dp, -0.000000000000_dp, &
                                           0.000000000000_dp, 0.004183310061_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                            0.006981609406_dp, -0.000000000000_dp, 1.259675258537_dp, 0.000000000000_dp, &
                                            0.000000000000_dp, -0.098427312376_dp, 0.058471503998_dp, 0.000000000000_dp, &
                                           -0.000000000000_dp, 0.000000000000_dp, -0.000000000000_dp, 0.009377471331_dp, &
                                           0.000000000000_dp, 1.384804093184_dp, -0.098427312376_dp, -0.000000000000_dp, &
                                            0.000000000000_dp, -0.013700945330_dp, 0.000000000000_dp, 0.000000000000_dp, &
                                           0.000000000000_dp, -0.009192048725_dp, 0.000000000000_dp, -0.098427312376_dp, &
                                            0.650988069850_dp, 0.000000000000_dp, -0.000000000000_dp, 0.047690906689_dp, &
                                         -0.000000000000_dp, -0.000000000000_dp, -0.009192048725_dp, -0.000000000000_dp, &
                                           -0.098427312376_dp, -0.000000000000_dp, 0.000000000000_dp, 0.746582265512_dp, &
                                           0.047690906689_dp, 0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                             0.021190548711_dp, 0.000000000000_dp, 0.058471503998_dp, 0.000000000000_dp, &
                                            -0.000000000000_dp, 0.047690906689_dp, 0.510245234739_dp, 0.000000000000_dp, &
                                           -0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, 0.004183310061_dp, &
                                            0.000000000000_dp, -0.013700945330_dp, 0.047690906689_dp, 0.000000000000_dp, &
                                           0.000000000000_dp, 0.586343429340_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                          -0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, &
                                         -0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                           0.407929748035_dp, 0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                           0.000000000000_dp, 0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, &
                              -0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, 0.507159756493_dp], [N_OV, N_OV])
   real(dp), parameter :: RHF_TRIPLET_A_EIG(N_OV) = [ &
                          0.407929748035_dp, 0.493107875364_dp, 0.507159756493_dp, &
                          0.559509762430_dp, 0.664359533721_dp, 0.742346594321_dp, &
                          1.281018295514_dp, 1.398256329686_dp, 20.044661050881_dp, &
                          20.114471636364_dp]
   real(dp), parameter :: PBE_TRIPLET_A(N_OV, N_OV) = reshape([ &
                                           18.735778954125_dp, 0.000000000000_dp, -0.005198341883_dp, 0.000000000000_dp, &
                                            0.000000000000_dp, -0.002468429366_dp, 0.002911017825_dp, 0.000000000000_dp, &
                                          -0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, 18.843978547526_dp, &
                                          0.000000000000_dp, -0.005849745439_dp, -0.002523027196_dp, -0.000000000000_dp, &
                                           0.000000000000_dp, -0.003704306619_dp, -0.000000000000_dp, 0.000000000000_dp, &
                                           -0.005198341883_dp, 0.000000000000_dp, 1.099749195038_dp, -0.000000000000_dp, &
                                            0.000000000000_dp, -0.031716282447_dp, 0.027507647096_dp, 0.000000000000_dp, &
                                          -0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, -0.005849745439_dp, &
                                          -0.000000000000_dp, 1.220043699155_dp, -0.034630603385_dp, -0.000000000000_dp, &
                                            0.000000000000_dp, 0.000256133849_dp, -0.000000000000_dp, 0.000000000000_dp, &
                                           0.000000000000_dp, -0.002523027196_dp, 0.000000000000_dp, -0.034630603385_dp, &
                                            0.646403184515_dp, 0.000000000000_dp, -0.000000000000_dp, 0.027066308972_dp, &
                                           0.000000000000_dp, 0.000000000000_dp, -0.002468429366_dp, -0.000000000000_dp, &
                                           -0.031716282447_dp, -0.000000000000_dp, 0.000000000000_dp, 0.732062077467_dp, &
                                            0.022740831929_dp, 0.000000000000_dp, 0.000000000000_dp, -0.000000000000_dp, &
                                             0.002911017825_dp, 0.000000000000_dp, 0.027507647096_dp, 0.000000000000_dp, &
                                            -0.000000000000_dp, 0.022740831929_dp, 0.438243244281_dp, 0.000000000000_dp, &
                                          -0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, -0.003704306619_dp, &
                                             0.000000000000_dp, 0.000256133849_dp, 0.027066308972_dp, 0.000000000000_dp, &
                                            0.000000000000_dp, 0.537318019557_dp, 0.000000000000_dp, -0.000000000000_dp, &
                                         -0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                            0.000000000000_dp, 0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, &
                                           0.346438094487_dp, 0.000000000000_dp, -0.000000000000_dp, -0.000000000000_dp, &
                                           -0.000000000000_dp, 0.000000000000_dp, 0.000000000000_dp, -0.000000000000_dp, &
                              -0.000000000000_dp, -0.000000000000_dp, 0.000000000000_dp, 0.455508915464_dp], [N_OV, N_OV])
   real(dp), parameter :: PBE_TRIPLET_A_EIG(N_OV) = [ &
                          0.346438094487_dp, 0.435145512444_dp, 0.455508915464_dp, &
                          0.530873240438_dp, 0.650760930701_dp, 0.731465596452_dp, &
                          1.103441080196_dp, 1.222127692703_dp, 18.735781281819_dp, &
                          18.843981586911_dp]

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

   integer, parameter :: N_CCPVDZ_TRIPLETS = 3
      !! Triplets asked for in the cc-pVDZ cases, which is how many the plan's
      !! tables carry.

   !! H2O/cc-pVDZ, the triplet roots, in Hartree, to twelve decimals, for the
   !! Tamm-Dancoff and the full problem.
   !!
   !! Regenerated rather than transcribed, as the singlets were. `get_ab`
   !! returns the singlet `A` and `B` whatever `td.singlet` says, so the RPA
   !! references are probed out of the paired `gen_vind`, which respects it,
   !! and reduced densely; for a pure functional PySCF routes that through
   !! `TDDFTNoHybrid`, whose `gen_vind` is the Casida matrix on one vector and
   !! whose eigenvalues are `w^2` outright. Reading the wrong one of those two
   !! shapes is silent and gives PBE roots four times too small.
   !!
   !! The plan's ten-digit triplet tables agree with these to 1.3e-9 except
   !! for B3LYP, whose whole column sits 2.3e-8 away in the plan -- singlets
   !! and triplets, TDA and RPA alike, so it is the plan's B3LYP reference and
   !! not anything about spin.
   real(dp), parameter :: RHF_CCPVDZ_TRIPLET(N_CCPVDZ_TRIPLETS) = [ &
                          0.304752995712_dp, 0.382448306917_dp, 0.383738345685_dp]
   real(dp), parameter :: PBE_CCPVDZ_TRIPLET(N_CCPVDZ_TRIPLETS) = [ &
                          0.245640585460_dp, 0.321365211218_dp, 0.322245501451_dp]
   real(dp), parameter :: B3LYP_CCPVDZ_TRIPLET(N_CCPVDZ_TRIPLETS) = [ &
                          0.254376470706_dp, 0.331305647032_dp, 0.331912037086_dp]
   real(dp), parameter :: CAM_CCPVDZ_TRIPLET(N_CCPVDZ_TRIPLETS) = [ &
                          0.256574750592_dp, 0.335035697057_dp, 0.336081588546_dp]
   real(dp), parameter :: RHF_CCPVDZ_TRIPLET_RPA(N_CCPVDZ_TRIPLETS) = [ &
                          0.299715458168_dp, 0.373956218363_dp, 0.376953572951_dp]
   real(dp), parameter :: PBE_CCPVDZ_TRIPLET_RPA(N_CCPVDZ_TRIPLETS) = [ &
                          0.244650838105_dp, 0.319741652798_dp, 0.321478892159_dp]
   real(dp), parameter :: B3LYP_CCPVDZ_TRIPLET_RPA(N_CCPVDZ_TRIPLETS) = [ &
                          0.253105574522_dp, 0.329819908019_dp, 0.330198597791_dp]
   real(dp), parameter :: CAM_CCPVDZ_TRIPLET_RPA(N_CCPVDZ_TRIPLETS) = [ &
                          0.255374444156_dp, 0.332919878823_dp, 0.335011282187_dp]

   !! Every root of the H2O/STO-3G Hartree-Fock **triplet** TDHF problem, the
   !! partner of `RHF_STO3G_RPA` below. Each sits under its Tamm-Dancoff
   !! partner in `RHF_TRIPLET_A_EIG`, as each singlet RPA root sits under its
   !! own, and the whole triplet spectrum sits under the singlet one.
   real(dp), parameter :: RHF_STO3G_TRIPLET_RPA(N_OV) = [ &
                          0.406101746192_dp, 0.474511202342_dp, 0.506576117767_dp, &
                          0.539432613518_dp, 0.659721221505_dp, 0.726963786914_dp, &
                          1.276436161205_dp, 1.394618237824_dp, 20.044633687193_dp, &
                          20.114433933272_dp]

   integer, parameter :: N_BOTH_STATES = 3
      !! Roots **per manifold** asked for by the `spin = "both"` case, so six
      !! come back.

   !! The union of the two H2O/STO-3G Tamm-Dancoff manifolds, as
   !! `spin = "both"` reports it: the three lowest of each, sorted together.
   !!
   !! The interleaving is the point. Concatenating the manifolds, or sorting
   !! within them and not across, gives a list in a different order that is
   !! made of the same numbers -- so the spins have to be checked alongside
   !! the energies, and both against a hand-merged reference rather than
   !! against whatever the code produced.
   real(dp), parameter :: BOTH_STO3G(2*N_BOTH_STATES) = [ &
                          0.407929748035_dp, 0.485080278078_dp, 0.493107875364_dp, &
                          0.507159756493_dp, 0.555807696384_dp, 0.617580730515_dp]
   integer, parameter :: BOTH_STO3G_SPIN(2*N_BOTH_STATES) = [ &
                         STATE_SPIN_TRIPLET, STATE_SPIN_SINGLET, STATE_SPIN_TRIPLET, &
                         STATE_SPIN_TRIPLET, STATE_SPIN_SINGLET, STATE_SPIN_SINGLET]

   !! Hartree-Fock carries no quadrature, so nothing but the integrals and the
   !! two SCF thresholds sits between the codes, and both are converged far
   !! below this. Measured on the STO-3G matrix: 7.7e-13 on the worst diagonal
   !! element, 6.0e-13 off it, 1.3e-12 on the orbital energies and 1.2e-12 on
   !! the spectrum. The bound is two orders over that, which leaves room for
   !! the compiler-to-compiler spread and still sits eight orders under the
   !! smallest fault it has to catch: halving the kernel, dropping one
   !! exchange term or symmetrising the wrong way each move elements by 1e-2
   !! and up. The triplet matrix lands in the same place, 6.1e-13 on the
   !! diagonal and 3.8e-13 off it.
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
   !!
   !! The triplet cases land in the same band: 3.5e-10 on the STO-3G PBE
   !! matrix and its spectrum, and on cc-pVDZ 1.8e-10 for PBE, 5.6e-10 for
   !! B3LYP and 1.9e-10 for CAM-B3LYP in the Tamm-Dancoff column, 1.9e-10,
   !! 6.6e-10 and 2.0e-10 in the RPA one. The polarised kernel is evaluated
   !! on the same points as the unpolarised one, so it adds no quadrature
   !! error of its own.
   real(dp), parameter :: TOL_GRID = 1.0e-7_dp

   !! The Davidson's own floor, against a dense diagonalisation of the same
   !! operator. Roots are accepted on a residual of 1e-8 here, and the
   !! eigenvalue error near an eigenvector is second order in the vector
   !! error, so this is a bound the solver clears by construction rather than
   !! one fitted to it -- measured 1.5e-11 on the worst of the five STO-3G
   !! roots, against both the dense spectrum and PySCF's, and 1.7e-11 on the
   !! worst of the six a `spin = "both"` solve merges.
   real(dp), parameter :: TOL_DAVIDSON = 1.0e-9_dp

   !! The cc-pVDZ Hartree-Fock roots. Looser than `TOL_EXACT` because these
   !! come out of the Davidson rather than off a dense diagonalisation, so the
   !! solver's own floor is in them as well as the integrals'. Measured
   !! 3.1e-12 on the worst of the five for TDA and 9.2e-11 for RPA, whose
   !! reduction amplifies the same disagreement the way it does at STO-3G.
   !! The triplets are 2.2e-12 and 2.4e-12 on the same two -- better than
   !! the singlet RPA because their reference was taken at `conv_tol = 1e-15`
   !! rather than 1e-13, which is what that 9.2e-11 mostly is.
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
   !!
   !! The triplet spectrum of the same system is 6.8e-13 against the same
   !! bound, for the reason the cc-pVDZ note below gives: its reference was
   !! taken at a tighter SCF.
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
   !! 9.0e-14 on the worst of five, and 3.3e-14 in the triplet manifold.
   real(dp), parameter :: TOL_CASIDA = 1.0e-9_dp

   ! --- Layer 6: the unrestricted gates -------------------------------------

   integer, parameter :: OH_N_OV = 130
      !! Spin-blocked rotations of the OH radical in cc-pVDZ: nineteen
      !! functions, five alpha and four beta electrons, so `5*14 + 4*15`.

   !! The OH radical, **in Bohr**: O at the origin, H at 0.9697 Angstrom on z.
   !!
   !! The plan's Angstrom distance times 1.8897261254578281, written out and
   !! handed to PySCF as Bohr, for the reason `WATER_BOHR` above is: the two
   !! codes carry different CODATA Bohr radii and converting on each side
   !! moves every orbital energy by 5e-10.
   real(dp), parameter :: OH_BOHR(3, 2) = reshape([ &
                                                  0.0_dp, 0.0_dp, 0.0_dp, &
                                                  0.0_dp, 0.0_dp, 1.832467423856456_dp], &
                                                  [3, 2])

   real(dp), parameter :: OH_UHF_ENERGY = -75.393846033464_dp

   !! `trace(A)` and `||A||_F` of the 130 by 130 unrestricted Tamm-Dancoff
   !! matrix.
   !!
   !! The matrix itself is too large to pin element by element, and its
   !! individual elements are not code-independent anyway -- OH is a 2-Pi
   !! radical, so its degenerate pi pair can be mixed arbitrarily between two
   !! codes. These two summaries are invariant under both the phase and that
   !! mixing, and between them every element contributes to one or the other.
   !! Their own run-to-run scatter from the threaded accumulation is 1e-11.
   real(dp), parameter :: OH_UHF_TRACE = 852.525539816128_dp
   real(dp), parameter :: OH_UHF_FROBENIUS = 119.939597984151_dp

   !! The five lowest unrestricted Tamm-Dancoff roots of OH / cc-pVDZ.
   !!
   !! **The first is not an excitation and is reported anyway.** 6.7e-3
   !! hartree is the rotation of the singly-occupied pi shell, which `A`
   !! alone is not singular along and the Tamm-Dancoff spectrum therefore
   !! keeps; the paired problem puts the same rotation at `w^2 = 0` and drops
   !! it, which is why the RPA list below starts one root higher. Neither is
   !! a fault in the solver, and `TDDFT_PLAN.md` says so.
   real(dp), parameter :: OH_UHF_TDA(5) = [ &
                          0.006697638062_dp, 0.173272241494_dp, 0.326244801914_dp, &
                          0.372883427182_dp, 0.431442526350_dp]

   !! The five lowest unrestricted RPA roots, the zero already dropped.
   real(dp), parameter :: OH_UHF_RPA(5) = [ &
                          0.169746047579_dp, 0.321231703197_dp, 0.370114042055_dp, &
                          0.415769391408_dp, 0.453280297821_dp]

   ! --- the Kohn-Sham gates, on a different doublet --------------------------
   !
   !! **Why the unrestricted Kohn-Sham cases are not the OH radical.**
   !!
   !! `TDDFT_PLAN.md` puts them there too, and they cannot go there. OH is a
   !! 2-Pi radical: the singly-occupied pi orbital is one of a degenerate
   !! pair, and a quadrature is not cylindrically symmetric, so the two
   !! orientations of that hole are **two distinct stationary points** of the
   !! Kohn-Sham energy, 6e-7 hartree apart. Both codes land on one or the
   !! other depending on the initial guess and, because the valley between
   !! them is nearly flat, on the order the threads finished the grid in:
   !! PySCF's `minao` and `atom` guesses reach one and its `1e` guess the
   !! other, and this program's own SCF was measured on both across two runs
   !! of the same test. The first root moves by 4.4e-5 between them, which is
   !! the rotation itself, and the next two by 3e-6.
   !!
   !! That is a property of the molecule, not of either code, and no
   !! tolerance makes it a gate. The water **cation** is the same
   !! unrestricted physics with a non-degenerate singly-occupied orbital, at
   !! the geometry this file already carries: PySCF converges it to
   !! |g| = 3e-10 from either guess, onto the same solution to twelve
   !! decimals. The Hartree-Fock gates stay on OH, where there is no
   !! quadrature to break the degeneracy and the two orientations are exactly
   !! degenerate.

   !! The water cation's five lowest roots, 175 spin-blocked rotations
   !! (`5*19 + 4*20`), from a dense diagonalisation of PySCF's own operator.
   !! `E(UKS PBE) = -75.881628961747`, `E(UKS B3LYP) = -75.967356416634`,
   !! both converged to `|g| < 5e-10` from either initial guess.
   real(dp), parameter :: CATION_PBE_TDA(5) = [ &
                          0.097595783031_dp, 0.240744153231_dp, 0.473433123899_dp, &
                          0.508215741547_dp, 0.515707923119_dp]

   real(dp), parameter :: CATION_B3LYP_TDA(5) = [ &
                          0.091902458748_dp, 0.237603694194_dp, 0.484648424585_dp, &
                          0.518092069836_dp, 0.529036148197_dp]

   real(dp), parameter :: CATION_B3LYP_RPA(5) = [ &
                          0.089312820588_dp, 0.236158714407_dp, 0.482329524466_dp, &
                          0.517089540365_dp, 0.527087149050_dp]

   !! Triplet H2 at 1.4 Bohr, a reference with two alpha electrons and **no
   !! beta electrons at all**: the beta spin contributes no rotations, so the
   !! trial vector is the alpha block alone and every beta half is empty.
   !! `E(UHF) = -0.766770390234`, `|g| = 1.8e-13`, 16 alpha excitations, and
   !! the Fock matrix is diagonal in both sets of orbitals to 1.2e-13, so
   !! these are canonical.
   !!
   !! Roots 3 and 4 are the two perpendicular pi components and are exactly
   !! degenerate. **Five roots, not three**: asked for three, the Davidson
   !! converges 1, 2 and 5 and reports the last as root 3 -- its guess is the
   !! three lowest diagonal gaps and `roots_to_solve` extends that over a
   !! degenerate *diagonal*, which this pair is not. That is the solver's
   !! guess, not the unrestricted operator, and it happens on a closed shell
   !! the same way; it is recorded here because a three-root gate on this
   !! molecule would pin the wrong spectrum.
   real(dp), parameter :: H2_TRIPLET_BOHR(3, 2) = reshape([ &
                                                          0.0_dp, 0.0_dp, 0.0_dp, &
                                                          0.0_dp, 0.0_dp, 1.4_dp], &
                                                          [3, 2])

   real(dp), parameter :: H2_TRIPLET_TDA(5) = [ &
                          0.251093280011_dp, 0.598938222312_dp, 0.870154715086_dp, &
                          0.870154715086_dp, 0.875564711106_dp]

   real(dp), parameter :: EXCITED_FLOOR = 1.0e-3_dp
      !! What the solver calls a rotation rather than an excitation, repeated
      !! here so the paired test can assert that the near-zero root fell below
      !! it rather than assume so.

   !! The unrestricted operator's spin sum and difference against the two
   !! restricted manifolds, for a Kohn-Sham reference.
   !!
   !! Looser than `TOL_EXACT` by an order because the two sides ask libxc for
   !! the same analytic quantity two ways: the restricted kernel out of the
   !! unpolarised functional, this one out of the polarised functional at
   !! `rho_a = rho_b`. Hartree-Fock has no such split and is held at
   !! `TOL_EXACT`.
   real(dp), parameter :: TOL_UKS_MANIFOLD = 1.0e-9_dp

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
                  new_unittest("an_unknown_method_or_spin_is_refused", &
                               test_later_layers), &
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
                  new_unittest("the_rhf_triplet_tda_matrix_matches_pyscf", &
                               test_rhf_triplet_matrix), &
                  new_unittest("the_pbe_triplet_tda_matrix_matches_pyscf", &
                               test_pbe_triplet_matrix), &
                  new_unittest("the_rhf_triplet_rpa_spectrum_matches_pyscf", &
                               test_rhf_triplet_rpa), &
                  new_unittest("cc_pvdz_rhf_triplets_match_pyscf", test_ccpvdz_rhf_t), &
                  new_unittest("cc_pvdz_pbe_triplets_match_pyscf", test_ccpvdz_pbe_t), &
                  new_unittest("cc_pvdz_b3lyp_triplets_match_pyscf", &
                               test_ccpvdz_b3lyp_t), &
                  new_unittest("cc_pvdz_cam_b3lyp_triplets_match_pyscf", &
                               test_ccpvdz_cam_t), &
                  new_unittest("cc_pvdz_rhf_triplet_rpa_matches_pyscf", &
                               test_ccpvdz_rhf_t_rpa), &
                  new_unittest("cc_pvdz_pbe_triplet_rpa_matches_pyscf", &
                               test_ccpvdz_pbe_t_rpa), &
                  new_unittest("cc_pvdz_b3lyp_triplet_rpa_matches_pyscf", &
                               test_ccpvdz_b3lyp_t_rpa), &
                  new_unittest("cc_pvdz_cam_b3lyp_triplet_rpa_matches_pyscf", &
                               test_ccpvdz_cam_t_rpa), &
                  new_unittest("the_casida_reduction_matches_the_triplet_solver", &
                               test_casida_triplet), &
                  new_unittest("both_manifolds_interleave_by_energy", test_both_spins), &
                  new_unittest("a_triplet_unstable_reference_is_named", test_instability), &
                  new_unittest("an_unreachable_tolerance_stops_and_says_so", &
                               test_unreachable_tolerance), &
                  new_unittest("the_oh_uhf_tda_matrix_matches_pyscf", &
                               test_oh_uhf_matrix), &
                  new_unittest("the_oh_uhf_rpa_spectrum_matches_pyscf", &
                               test_oh_uhf_rpa_matrix), &
                  new_unittest("oh_uhf_tda_roots_match_pyscf", &
                               test_oh_uhf_tda_solver), &
                  new_unittest("oh_uhf_rpa_roots_match_pyscf", &
                               test_oh_uhf_rpa_solver), &
                  new_unittest("cation_uks_pbe_tda_roots_match_pyscf", &
                               test_cation_uks_pbe), &
                  new_unittest("cation_uks_b3lyp_tda_roots_match_pyscf", &
                               test_cation_uks_b3lyp), &
                  new_unittest("cation_uks_b3lyp_rpa_roots_match_pyscf", &
                               test_cation_uks_b3lyp_rpa), &
                  new_unittest("a_reference_with_no_beta_electrons_has_a_spectrum", &
                               test_no_beta_electrons), &
                  new_unittest("unrestricted_amplitudes_carry_unit_norm", &
                               test_uhf_amplitude_norm), &
                  new_unittest("the_unrestricted_hf_operator_holds_both_manifolds", &
                               test_restricted_from_unrestricted_hf), &
                  new_unittest("the_unrestricted_pbe_operator_holds_both_manifolds", &
                               test_restricted_from_unrestricted_pbe) &
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

   subroutine dense_tda(mol, scf, ctx, kohn_sham, a, err, spin)
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
      character(len=*), intent(in), optional :: spin
         !! Which manifold, `singlet` when absent -- so every Layer 2 case
         !! below reaches exactly the operator it always did.

      type(tda_operator_t) :: operator
      character(len=16) :: manifold

      if (err%has_error()) return
      manifold = "singlet"
      if (present(spin)) manifold = spin

      if (kohn_sham) then
         call build_tda_operator(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, operator, err, xc=ctx, &
                                 reference=scf%density, spin=trim(manifold))
      else
         call build_tda_operator(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, operator, err, spin=trim(manifold))
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

   subroutine excited_run(basis, functional, n_states, method, spin, result, tolerance)
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
      if (present(tolerance)) settings%excited%tolerance = tolerance
      settings%excited%max_iter = 200

      call run_czt_hf(settings, fragment, result)
   end subroutine excited_run

   subroutine compare_roots(error, result, reference, tol, what, spin_code)
      !! Every root of the table, against what the run reported
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), intent(in) :: result
      real(dp), intent(in) :: reference(:), tol
      character(len=*), intent(in) :: what
      integer, intent(in), optional :: spin_code
         !! What every root should be labelled; `STATE_SPIN_SINGLET` absent.

      integer :: i, want_spin

      want_spin = STATE_SPIN_SINGLET
      if (present(spin_code)) want_spin = spin_code

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
      call check(error, all(result%state_spin == want_spin), &
                 "a single-manifold solve reported a root labelled with the other "// &
                 "spin")
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

   subroutine test_rhf_triplet_matrix(error)
      !! The ten-by-ten triplet TDA matrix of H2O/STO-3G, element by element
      !!
      !! Two independent changes separate this matrix from the singlet one:
      !! the Coulomb term is gone, and the kernel -- absent here, since this
      !! is Hartree-Fock -- would be the spin difference. So for Hartree-Fock
      !! this case is the `j_scale = 0` half of Layer 3 isolated from the
      !! kernel, and the PBE case below is the other half.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
      type(error_t) :: err
      real(dp), allocatable :: a(:, :), values(:)
      logical :: ok

      call water_sto3g(mol, scf, ctx, err)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "the Hartree-Fock reference failed: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call dense_tda(mol, scf, ctx, .false., a, err, spin="triplet")
      call mol%destroy()
      call check(error,.not. err%has_error(), "the triplet TDA operator failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call compare_matrix(error, a, RHF_TRIPLET_A, TOL_EXACT, "Hartree-Fock triplet")
      if (allocated(error)) return

      values = eigenvalues_of(a, ok)
      call check(error, ok, "the dense diagonalisation of the triplet matrix failed")
      if (allocated(error)) return
      call check(error, maxval(abs(values - RHF_TRIPLET_A_EIG)) < TOL_EXACT, &
                 "the Hartree-Fock triplet TDA spectrum disagrees with PySCF")
      if (allocated(error)) return

      ! Every triplet below its singlet partner, which is Hund's rule and a
      ! check no tolerance would catch if the two spectra had been
      ! transcribed into each other's places.
      call check(error, all(RHF_TRIPLET_A_EIG < RHF_A_EIG), "a triplet root came "// &
                 "out above its singlet partner")
   end subroutine test_rhf_triplet_matrix

   subroutine test_pbe_triplet_matrix(error)
      !! The same, for a pure functional, where the triplet kernel is the
      !! whole of the coupling
      !!
      !! PBE carries no exact exchange and a triplet has no Coulomb term, so
      !! every element off the orbital-energy diagonal here comes from
      !! `(f_aa - f_ab)/2` and nothing else. A kernel evaluated unpolarised --
      !! that is, the singlet one left in place -- moves the first root by
      !! 0.07 hartree, and the four GGA combinations mis-weighted against each
      !! other move it by less but never by less than this bound.
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

      call dense_tda(mol, scf, ctx, .true., a, err, spin="triplet")
      call mol%destroy()
      call check(error,.not. err%has_error(), "the PBE triplet TDA operator failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call compare_matrix(error, a, PBE_TRIPLET_A, TOL_GRID, "PBE triplet")
      if (allocated(error)) return

      values = eigenvalues_of(a, ok)
      call check(error, ok, "the dense diagonalisation of the PBE triplet matrix "// &
                 "failed")
      if (allocated(error)) return
      call check(error, maxval(abs(values - PBE_TRIPLET_A_EIG)) < TOL_GRID, &
                 "the PBE triplet TDA spectrum disagrees with PySCF")
   end subroutine test_pbe_triplet_matrix

   subroutine test_rhf_triplet_rpa(error)
      !! Every triplet RPA root of H2O/STO-3G, from the two explicit halves
      !!
      !! `(A-B)` is exchange only and so is the same operator for both spins;
      !! the whole of the difference is in `(A+B)`. This case reads the two
      !! apart, which the Tamm-Dancoff half sum above cannot, so an error
      !! that moved them oppositely would survive that gate and fail here.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
      type(error_t) :: err
      real(dp), allocatable :: aplus(:, :), aminus(:, :)
      real(dp), allocatable :: splus(:, :), sminus(:, :), values(:)
      logical :: ok

      call water_sto3g(mol, scf, ctx, err)
      if (.not. err%has_error()) &
         call dense_rpa(mol, scf, ctx, .false., aplus, aminus, err, spin="triplet")
      if (.not. err%has_error()) &
         call dense_rpa(mol, scf, ctx, .false., splus, sminus, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the triplet RPA matrices failed: "// &
                 err%get_message())
      if (allocated(error)) return

      ! The two halves average to the triplet `A` already pinned above.
      call compare_matrix(error, 0.5_dp*(aplus + aminus), RHF_TRIPLET_A, TOL_EXACT, &
                          "triplet half-sum")
      if (allocated(error)) return

      ! `(A-B)` is exchange only, so it does not know which manifold it is
      ! in. Not a tolerance on a physical quantity: the two builds are the
      ! same arithmetic on the same densities, so they agree exactly.
      call check(error, maxval(abs(aminus - sminus)) == 0.0_dp, &
                 "the triplet (A-B) differs from the singlet one, which an "// &
                 "exchange-only operator cannot")
      if (allocated(error)) return

      values = paired_spectrum(aplus, aminus, ok)
      call check(error, ok, "(A-B) of a converged closed shell was not positive "// &
                 "definite, so the triplet reduction could not be formed")
      if (allocated(error)) return
      call check(error, maxval(abs(values - RHF_STO3G_TRIPLET_RPA)) < TOL_RPA_DENSE, &
                 "the triplet RPA spectrum disagrees with PySCF")
      if (allocated(error)) return
      call check(error, all(RHF_STO3G_TRIPLET_RPA < RHF_TRIPLET_A_EIG), &
                 "a triplet RPA root came out above its Tamm-Dancoff partner")
   end subroutine test_rhf_triplet_rpa

   subroutine test_ccpvdz_rhf_t(error)
      !! The three Hartree-Fock triplets of H2O/cc-pVDZ
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      call excited_run("cc-pvdz", "", N_CCPVDZ_TRIPLETS, "tda", "triplet", result)
      call compare_roots(error, result, RHF_CCPVDZ_TRIPLET, TOL_CCPVDZ_HF, &
                         "cc-pVDZ RHF triplet", spin_code=STATE_SPIN_TRIPLET)
   end subroutine test_ccpvdz_rhf_t

   subroutine test_ccpvdz_pbe_t(error)
      !! The three PBE triplets
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "pbe", N_CCPVDZ_TRIPLETS, "tda", "triplet", result)
      call compare_roots(error, result, PBE_CCPVDZ_TRIPLET, TOL_GRID, &
                         "cc-pVDZ PBE triplet", spin_code=STATE_SPIN_TRIPLET)
   end subroutine test_ccpvdz_pbe_t

   subroutine test_ccpvdz_b3lyp_t(error)
      !! The three B3LYP triplets, libxc 402 on both sides
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "b3lyp", N_CCPVDZ_TRIPLETS, "tda", "triplet", result)
      call compare_roots(error, result, B3LYP_CCPVDZ_TRIPLET, TOL_GRID, &
                         "cc-pVDZ B3LYP triplet", spin_code=STATE_SPIN_TRIPLET)
   end subroutine test_ccpvdz_b3lyp_t

   subroutine test_ccpvdz_cam_t(error)
      !! The three CAM-B3LYP triplets, where both exchange passes and the
      !! triplet kernel have to be right at once
      !!
      !! The only Tamm-Dancoff case here exercising the attenuated exchange
      !! build and the polarised kernel together. Each is separately visible
      !! in an earlier case, so a failure that appears only here is the
      !! combination -- most likely the long-range pass being dropped when the
      !! Coulomb term is.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "cam-b3lyp", N_CCPVDZ_TRIPLETS, "tda", "triplet", &
                       result)
      call compare_roots(error, result, CAM_CCPVDZ_TRIPLET, TOL_GRID, &
                         "cc-pVDZ CAM-B3LYP triplet", spin_code=STATE_SPIN_TRIPLET)
   end subroutine test_ccpvdz_cam_t

   subroutine test_ccpvdz_rhf_t_rpa(error)
      !! The three Hartree-Fock triplets of the full problem
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      call excited_run("cc-pvdz", "", N_CCPVDZ_TRIPLETS, "rpa", "triplet", result)
      call compare_roots(error, result, RHF_CCPVDZ_TRIPLET_RPA, TOL_CCPVDZ_HF, &
                         "cc-pVDZ RHF triplet RPA", spin_code=STATE_SPIN_TRIPLET)
   end subroutine test_ccpvdz_rhf_t_rpa

   subroutine test_ccpvdz_pbe_t_rpa(error)
      !! The three PBE triplets of the full problem
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "pbe", N_CCPVDZ_TRIPLETS, "rpa", "triplet", result)
      call compare_roots(error, result, PBE_CCPVDZ_TRIPLET_RPA, TOL_GRID, &
                         "cc-pVDZ PBE triplet RPA", spin_code=STATE_SPIN_TRIPLET)
   end subroutine test_ccpvdz_pbe_t_rpa

   subroutine test_ccpvdz_b3lyp_t_rpa(error)
      !! The three B3LYP triplets of the full problem
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "b3lyp", N_CCPVDZ_TRIPLETS, "rpa", "triplet", result)
      call compare_roots(error, result, B3LYP_CCPVDZ_TRIPLET_RPA, TOL_GRID, &
                         "cc-pVDZ B3LYP triplet RPA", spin_code=STATE_SPIN_TRIPLET)
   end subroutine test_ccpvdz_b3lyp_t_rpa

   subroutine test_ccpvdz_cam_t_rpa(error)
      !! The three CAM-B3LYP triplets of the full problem
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      if (.not. xc_available()) return

      call excited_run("cc-pvdz", "cam-b3lyp", N_CCPVDZ_TRIPLETS, "rpa", "triplet", &
                       result)
      call compare_roots(error, result, CAM_CCPVDZ_TRIPLET_RPA, TOL_GRID, &
                         "cc-pVDZ CAM-B3LYP triplet RPA", spin_code=STATE_SPIN_TRIPLET)
   end subroutine test_ccpvdz_cam_t_rpa

   subroutine test_casida_triplet(error)
      !! The Casida reduction and the paired solver agree on a triplet too
      !!
      !! The reduction's assumption is about `(A-B)`, which for a triplet is
      !! the same exchange-only operator it is for a singlet -- so on a pure
      !! functional it is `diag(dEps)` in both manifolds and the cross-check
      !! carries over unchanged. Worth running rather than asserting: the two
      !! routes share only the `(A+B)` product, and only one of them applies
      !! the triplet kernel through a square root.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
      type(error_t) :: err
      real(dp), allocatable :: paired(:), casida(:), x(:, :), y(:, :)
      integer, allocatable :: spins(:)

      if (.not. xc_available()) return

      call water_sto3g(mol, scf, ctx, err, functional="pbe")
      if (err%has_error()) then
         call mol%destroy()
         call check(error, .false., "the PBE reference failed: "//err%get_message())
         return
      end if

      call response_excitations(mol, scf%orbitals, scf%orbital_energies, &
                                scf%n_occupied, 5, "rpa", "triplet", paired, spins, &
                                x, y, err, xc=ctx, reference=scf%density, &
                                tolerance=1.0e-10_dp, max_iter=100)
      if (.not. err%has_error()) &
         call response_excitations(mol, scf%orbitals, scf%orbital_energies, &
                                   scf%n_occupied, 5, "casida", "triplet", casida, &
                                   spins, x, y, err, xc=ctx, reference=scf%density, &
                                   tolerance=1.0e-12_dp, max_iter=200)
      call mol%destroy()
      call check(error,.not. err%has_error(), "a PBE triplet solve failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, size(paired) == 5 .and. size(casida) == 5, &
                 "the two routes converged different numbers of triplet roots")
      if (allocated(error)) return
      call check(error, maxval(abs(paired - casida)) < TOL_CASIDA, &
                 "the Casida reduction and the paired solver disagree on the PBE "// &
                 "triplet spectrum")
      if (allocated(error)) return
      call check(error, all(spins == STATE_SPIN_TRIPLET), &
                 "a triplet solve labelled a root a singlet")
   end subroutine test_casida_triplet

   subroutine test_both_spins(error)
      !! `spin = "both"` returns one spectrum, sorted across the manifolds
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result
      integer :: i

      call excited_run("sto-3g", "", N_BOTH_STATES, "tda", "both", result)
      call check(error,.not. result%has_error, "the both-manifolds run failed: "// &
                 result%error%get_message())
      if (allocated(error)) return
      call check(error, result%has_excited_states, "the both-manifolds run reported "// &
                 "no excited states")
      if (allocated(error)) return
      call check(error, size(result%excitation_energies) == 2*N_BOTH_STATES, &
                 "asking for both manifolds did not return both manifolds' roots")
      if (allocated(error)) return
      call check(error, allocated(result%state_spin), "the both-manifolds run "// &
                 "labelled no spins")
      if (allocated(error)) return
      call check(error, size(result%state_spin) == 2*N_BOTH_STATES, &
                 "the spin labels are not one per root")
      if (allocated(error)) return

      do i = 1, 2*N_BOTH_STATES
         call check(error, abs(result%excitation_energies(i) - BOTH_STO3G(i)) < &
                    TOL_DAVIDSON, "a root of the merged spectrum is not the one "// &
                    "that belongs at its position")
         if (allocated(error)) return
         call check(error, result%state_spin(i) == BOTH_STO3G_SPIN(i), &
                    "a root of the merged spectrum carries the wrong spin label")
         if (allocated(error)) return
      end do
   end subroutine test_both_spins

   subroutine h2_stretched_run(n_states, method, spin, result)
      !! H2 at 3.0 Angstrom through the bridge, restricted Hartree-Fock
      integer, intent(in) :: n_states
      character(len=*), intent(in) :: method, spin
      type(calculation_result_t), intent(out) :: result

      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment

      fragment%n_atoms = 2
      fragment%charge = 0
      fragment%multiplicity = 1
      fragment%nelec = 2
      fragment%n_caps = 0
      allocate (fragment%element_numbers(2), fragment%coordinates(3, 2))
      fragment%element_numbers = [1, 1]
      fragment%coordinates = H2_BOHR

      settings%basis_set = "6-31g"
      settings%functional = ""
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
   end subroutine h2_stretched_run

   subroutine test_instability(error)
      !! A triplet-unstable reference is named, not reported as a number
      !!
      !! Stretched H2 is where the restricted solution stops being a minimum
      !! against spin polarisation: the lowest triplet `A` eigenvalue is
      !! -0.172 hartree and one squared RPA frequency is -0.0118. Its
      !! *singlet* channel is perfectly well behaved -- `(A-B)` is positive
      !! definite in both manifolds, because it does not depend on the spin --
      !! which is the second half of this case: the diagnosis has to be about
      !! the triplet operator and not about the molecule being awkward.
      !!
      !! Two wrong answers are possible and both look plausible from outside.
      !! The Tamm-Dancoff route can report the negative root as an excitation
      !! or drop it under the floor and hand back the ones above it as though
      !! a state were merely missing. The paired route is worse: it *skips*
      !! an imaginary frequency and only complains when it runs out of real
      !! ones, so asking for one root succeeds and says nothing.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result
      character(len=:), allocatable :: message

      call h2_stretched_run(1, "tda", "triplet", result)
      call check(error, result%has_error, "a triplet-unstable reference produced a "// &
                 "Tamm-Dancoff spectrum rather than a diagnosis")
      if (allocated(error)) return
      message = result%error%get_message()
      call check(error, index(message, "triplet-unstable") > 0, &
                 "the failure of a triplet-unstable reference was reported as "// &
                 "something else: "//message)
      if (allocated(error)) return
      call check(error,.not. result%has_excited_states, "a refused triplet-unstable "// &
                 "run still reported excited states")
      if (allocated(error)) return

      ! One root, which the paired solver would otherwise answer happily out
      ! of the two real frequencies above the imaginary one.
      call h2_stretched_run(1, "rpa", "triplet", result)
      call check(error, result%has_error, "a triplet-unstable reference produced an "// &
                 "RPA spectrum rather than a diagnosis")
      if (allocated(error)) return
      message = result%error%get_message()
      call check(error, index(message, "triplet-unstable") > 0, &
                 "the failure of a triplet-unstable paired solve was reported as "// &
                 "something else: "//message)
      if (allocated(error)) return

      ! The same molecule, same reference, singlet manifold: an ordinary
      ! answer. Without this the case would pass just as well for a build
      ! that refused every excited-state run on stretched H2.
      call h2_stretched_run(1, "tda", "singlet", result)
      call check(error,.not. result%has_error, "the singlet manifold of the same "// &
                 "reference failed too: "//result%error%get_message())
      if (allocated(error)) return
      call check(error, abs(result%energy%scf - H2_631G_ENERGY) < TOL_CCPVDZ_HF, &
                 "the stretched-H2 Hartree-Fock energy is not PySCF's")
      if (allocated(error)) return
      call check(error, result%has_excited_states, "the singlet manifold of a "// &
                 "triplet-unstable reference reported no excited states")
   end subroutine test_instability

   subroutine test_later_layers(error)
      !! What is still not implemented is refused by name, not approximated
      !!
      !! Both of this case's original halves are gone: Layer 4 answers an RPA
      !! deck and Layer 3 answers a triplet one. What is left is the check
      !! that a method string naming neither is refused rather than resolved
      !! to whichever is nearer -- the two are different eigenproblems and
      !! both converge. The spin string has the same guard, one line below.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      call excited_run("sto-3g", "", 2, "cis", "singlet", result)
      call check(error, result%has_error, "an unknown excited-state method was "// &
                 "resolved to something rather than refused")
      if (allocated(error)) return
      call check(error,.not. result%has_excited_states, "a refused deck still "// &
                 "reported excited states")
      if (allocated(error)) return

      call excited_run("sto-3g", "", 2, "tda", "quintet", result)
      call check(error, result%has_error, "an unknown excited-state spin was "// &
                 "resolved to something rather than refused")
   end subroutine test_later_layers

   subroutine dense_rpa(mol, scf, ctx, kohn_sham, aplus, aminus, err, spin)
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
      character(len=*), intent(in), optional :: spin
         !! Which manifold, `singlet` when absent.

      type(rpa_operator_t) :: operator
      character(len=16) :: manifold

      if (err%has_error()) return
      manifold = "singlet"
      if (present(spin)) manifold = spin

      if (kohn_sham) then
         call build_rpa_operator(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, operator, err, xc=ctx, &
                                 reference=scf%density, spin=trim(manifold))
      else
         call build_rpa_operator(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, operator, err, spin=trim(manifold))
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
      integer, allocatable :: spins(:)
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

      call response_excitations(mol, scf%orbitals, scf%orbital_energies, &
                                scf%n_occupied, 5, "rpa", "singlet", omega, spins, &
                                x, y, err, tolerance=1.0e-10_dp, max_iter=100)
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
      integer, allocatable :: spins(:)

      if (.not. xc_available()) return

      call water_sto3g(mol, scf, ctx, err, functional="pbe")
      if (err%has_error()) then
         call mol%destroy()
         call check(error, .false., "the PBE reference failed: "//err%get_message())
         return
      end if

      call response_excitations(mol, scf%orbitals, scf%orbital_energies, &
                                scf%n_occupied, 5, "rpa", "singlet", paired, spins, &
                                x, y, err, xc=ctx, reference=scf%density, &
                                tolerance=1.0e-10_dp, max_iter=100)
      if (.not. err%has_error()) &
         call response_excitations(mol, scf%orbitals, scf%orbital_energies, &
                                   scf%n_occupied, 5, "casida", "singlet", casida, &
                                   spins, x, y, err, xc=ctx, reference=scf%density, &
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
      integer, allocatable :: spins(:)

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

      call response_excitations(mol, scf%orbitals, scf%orbital_energies, &
                                scf%n_occupied, 3, "rpa", "singlet", omega, spins, &
                                x, y, err, tolerance=1.0e-10_dp, max_iter=100)
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

   subroutine oh_reference(mol, scf, ctx, err, functional)
      !! Converge the OH radical in cc-pVDZ, unrestricted, at the plan geometry
      !!
      !! A doublet with a small beta gap -- 0.8 eV between the beta HOMO and
      !! LUMO -- so the SCF is driven hard: what is compared is an operator
      !! built from the orbitals, and the orbital error goes as the commutator
      !! rather than its square.
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(xc_context_t), intent(out) :: ctx
      type(error_t), intent(inout) :: err
      character(len=*), intent(in), optional :: functional
         !! Absent is Hartree-Fock, and leaves `ctx` untouched.

      call build_czt_molecule([8, 1], ["O ", "H "], OH_BOHR, "cc-pvdz", mol, err)
      if (err%has_error()) return

      if (present(functional)) then
         call xc_context_create(mol, functional, ctx, err, level=5, polarized=.true.)
         if (err%has_error()) return
         call run_czt_uhf(mol, 9, 2, 400, 1.0e-13_dp, 1.0e-10_dp, .false., scf, err, &
                          xc=ctx, grad_tol=1.0e-9_dp)
      else
         ! 1e-11 on the commutator was reachable on most runs and not on all:
         ! the doublet's beta gap is small and the last decade of the DIIS
         ! wanders with the thread schedule. 1e-10 is reached every time and
         ! is two decades below what the 1e-10 eigenvalue gate needs.
         call run_czt_uhf(mol, 9, 2, 500, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                          grad_tol=1.0e-10_dp)
      end if
   end subroutine oh_reference

   subroutine dense_tda_uhf(mol, scf, ctx, kohn_sham, a, err)
      !! The explicit unrestricted TDA matrix, through the shipped operator
      !!
      !! Both spin blocks at once, so `a` is the `(n_ov_a + n_ov_b)` square
      !! with the coupling blocks in it. Probed with unit vectors, which is
      !! the construction PySCF's `gen_vind` side of the reference used too.
      type(czt_molecule_t), intent(in), target :: mol
      type(rhf_result_t), intent(in) :: scf
      type(xc_context_t), intent(inout), target :: ctx
      logical, intent(in) :: kohn_sham
      real(dp), allocatable, intent(out) :: a(:, :)
      type(error_t), intent(inout) :: err

      type(tda_operator_uhf_t) :: operator

      if (err%has_error()) return
      if (kohn_sham) then
         call build_tda_operator_uhf(mol, scf%orbitals, scf%orbital_energies, &
                                     scf%n_occupied, scf%orbitals_beta, &
                                     scf%orbital_energies_beta, scf%n_occupied_beta, &
                                     operator, err, xc=ctx, ref_a=scf%density, &
                                     ref_b=scf%density_beta)
      else
         call build_tda_operator_uhf(mol, scf%orbitals, scf%orbital_energies, &
                                     scf%n_occupied, scf%orbitals_beta, &
                                     scf%orbital_energies_beta, scf%n_occupied_beta, &
                                     operator, err)
      end if
      if (err%has_error()) return

      call tda_dense_matrix_uhf(operator, a, err)
   end subroutine dense_tda_uhf

   subroutine dense_rpa_uhf(mol, scf, ctx, kohn_sham, aplus, aminus, err)
      !! The explicit unrestricted `(A+B)` and `(A-B)`, through the same operator
      type(czt_molecule_t), intent(in), target :: mol
      type(rhf_result_t), intent(in) :: scf
      type(xc_context_t), intent(inout), target :: ctx
      logical, intent(in) :: kohn_sham
      real(dp), allocatable, intent(out) :: aplus(:, :), aminus(:, :)
      type(error_t), intent(inout) :: err

      type(rpa_operator_uhf_t) :: operator

      if (err%has_error()) return
      if (kohn_sham) then
         call build_rpa_operator_uhf(mol, scf%orbitals, scf%orbital_energies, &
                                     scf%n_occupied, scf%orbitals_beta, &
                                     scf%orbital_energies_beta, scf%n_occupied_beta, &
                                     operator, err, xc=ctx, ref_a=scf%density, &
                                     ref_b=scf%density_beta)
      else
         call build_rpa_operator_uhf(mol, scf%orbitals, scf%orbital_energies, &
                                     scf%n_occupied, scf%orbitals_beta, &
                                     scf%orbital_energies_beta, scf%n_occupied_beta, &
                                     operator, err)
      end if
      if (err%has_error()) return

      call rpa_dense_matrices_uhf(operator, aplus, aminus, err)
   end subroutine dense_rpa_uhf

   subroutine oh_fragment(fragment)
      !! The OH radical as the bridge wants it: element numbers, Bohr, doublet
      type(physical_fragment_t), intent(out) :: fragment

      fragment%n_atoms = 2
      fragment%charge = 0
      fragment%multiplicity = 2
      fragment%nelec = 9
      fragment%n_caps = 0
      allocate (fragment%element_numbers(2), fragment%coordinates(3, 2))
      fragment%element_numbers = [8, 1]
      fragment%coordinates = OH_BOHR
   end subroutine oh_fragment

   subroutine cation_excited_run(functional, n_states, method, result)
      !! The water cation through the bridge, unrestricted, with a spectrum
      !!
      !! The same geometry as every restricted case in this file, one electron
      !! short: a doublet whose singly-occupied orbital is not degenerate, so
      !! the Kohn-Sham solution is unique and a 1e-7 comparison means
      !! something. See the note above `CATION_PBE_TDA`.
      character(len=*), intent(in) :: functional, method
      integer, intent(in) :: n_states
      type(calculation_result_t), intent(out) :: result

      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment

      call water_fragment(fragment)
      fragment%charge = 1
      fragment%multiplicity = 2
      fragment%nelec = 9
      settings%basis_set = "cc-pvdz"
      settings%functional = functional
      settings%grid_level = 5
      settings%energy_tol = 1.0e-12_dp
      settings%grad_tol = 1.0e-9_dp
      settings%density_tol = 1.0e-9_dp
      settings%max_iter = 300
      settings%excited%enabled = n_states > 0
      settings%excited%n_states = n_states
      settings%excited%method = method
      settings%excited%spin = "singlet"
      settings%excited%tolerance = 1.0e-9_dp
      settings%excited%max_iter = 200

      call run_czt_hf(settings, fragment, result)
   end subroutine cation_excited_run

   subroutine oh_excited_run(functional, n_states, method, result)
      !! One whole unrestricted calculation through the bridge
      character(len=*), intent(in) :: functional, method
      integer, intent(in) :: n_states
      type(calculation_result_t), intent(out) :: result

      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment

      call oh_fragment(fragment)
      settings%basis_set = "cc-pvdz"
      settings%functional = functional
      settings%grid_level = 5
      settings%energy_tol = 1.0e-12_dp
      settings%grad_tol = 1.0e-9_dp
      settings%density_tol = 1.0e-9_dp
      settings%max_iter = 400
      settings%excited%enabled = n_states > 0
      settings%excited%n_states = n_states
      settings%excited%method = method
      settings%excited%spin = "singlet"
      settings%excited%tolerance = 1.0e-9_dp
      settings%excited%max_iter = 200

      call run_czt_hf(settings, fragment, result)
   end subroutine oh_excited_run

   function invariants_of(a) result(pair)
      !! `trace(A)` and `||A||_F`, the two summaries a phase cannot move
      !!
      !! Two codes converging the same open shell agree on the orbitals only
      !! up to a sign per orbital, and on a degenerate block only up to an
      !! orthogonal mixing inside it. `A` is covariant under both -- it
      !! carries one occupied and one virtual index on each side -- so its
      !! trace and its Frobenius norm are the same numbers in either code
      !! while almost no individual element is. A 130 by 130 matrix is too
      !! large to pin element by element in a test file; these two are what
      !! can be pinned, and between them they see every element.
      real(dp), intent(in) :: a(:, :)
      real(dp) :: pair(2)

      integer :: i

      pair = 0.0_dp
      do i = 1, size(a, 1)
         pair(1) = pair(1) + a(i, i)
      end do
      pair(2) = sqrt(sum(a*a))
   end function invariants_of

   subroutine test_oh_uhf_matrix(error)
      !! The unrestricted TDA operator of the OH radical, against PySCF
      !!
      !! The whole 130 by 130 matrix is built by probing the shipped operator
      !! with every unit vector, which is what PySCF's `gen_vind` side of the
      !! reference did as well. What is compared is its symmetry, its two
      !! phase-independent invariants and its five lowest eigenvalues.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
      type(error_t) :: err
      real(dp), allocatable :: a(:, :), values(:)
      real(dp) :: pair(2)
      logical :: ok

      call oh_reference(mol, scf, ctx, err)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "the unrestricted Hartree-Fock reference failed: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call check(error, abs(scf%energy - OH_UHF_ENERGY) < TOL_EXACT, &
                 "the OH unrestricted Hartree-Fock energy is not PySCF's")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call dense_tda_uhf(mol, scf, ctx, .false., a, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the dense unrestricted TDA matrix "// &
                 "failed: "//err%get_message())
      if (allocated(error)) return

      call check(error, size(a, 1) == OH_N_OV, "the unrestricted operator is not the "// &
                 "length of the two occupied-virtual blocks together")
      if (allocated(error)) return
      ! `A` is symmetric for a real reference, and nothing in the build
      ! enforces it: the two spin blocks come from separate transforms and the
      ! coupling blocks from the Coulomb term of one against the other.
      call check(error, maxval(abs(a - transpose(a))) < TOL_EXACT, &
                 "the unrestricted TDA matrix is not symmetric")
      if (allocated(error)) return

      pair = invariants_of(a)
      call check(error, abs(pair(1) - OH_UHF_TRACE) < 1.0e-8_dp, &
                 "the trace of the unrestricted TDA matrix is not PySCF's")
      if (allocated(error)) return
      call check(error, abs(pair(2) - OH_UHF_FROBENIUS) < 1.0e-8_dp, &
                 "the Frobenius norm of the unrestricted TDA matrix is not PySCF's")
      if (allocated(error)) return

      values = eigenvalues_of(a, ok)
      call check(error, ok, "the dense diagonalisation failed")
      if (allocated(error)) return
      call check(error, maxval(abs(values(1:5) - OH_UHF_TDA)) < TOL_EXACT, &
                 "an OH unrestricted TDA root disagrees with PySCF")
   end subroutine test_oh_uhf_matrix

   subroutine test_oh_uhf_rpa_matrix(error)
      !! The paired unrestricted problem of OH, from the two explicit halves
      !!
      !! The lowest `w` here is the numerical zero the Tamm-Dancoff spectrum
      !! keeps as 6.7e-3: the half-filled shell's own rotation, which the
      !! paired problem puts at `w^2 = 0` because `(A+B)` is singular along
      !! it. It is skipped rather than compared, and skipping it is the
      !! statement that it is not an excitation.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
      type(error_t) :: err
      real(dp), allocatable :: aplus(:, :), aminus(:, :), values(:)
      logical :: ok

      call oh_reference(mol, scf, ctx, err)
      if (.not. err%has_error()) call dense_rpa_uhf(mol, scf, ctx, .false., aplus, &
                                                    aminus, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the dense unrestricted paired halves "// &
                 "failed: "//err%get_message())
      if (allocated(error)) return

      values = paired_spectrum(aplus, aminus, ok)
      call check(error, ok, "the dense paired reduction failed; (A-B) should be "// &
                 "positive definite on this doublet")
      if (allocated(error)) return
      ! `paired_spectrum` takes the square root of every `w^2`, and this one
      ! is a numerical zero that lands on either side of it -- PySCF's own
      ! comes out at +4e-15 on one run and -1e-15 on the next -- so what comes
      ! back here is either a number far below the floor or a NaN. Both say
      ! the same thing, and a comparison a NaN fails is how that is written.
      call check(error,.not. (values(1) >= EXCITED_FLOOR), "the rotation of the "// &
                 "half-filled shell did not come back at the numerical zero")
      if (allocated(error)) return
      call check(error, maxval(abs(values(2:6) - OH_UHF_RPA)) < TOL_RPA_DENSE, &
                 "an OH unrestricted RPA root disagrees with PySCF")
   end subroutine test_oh_uhf_rpa_matrix

   subroutine compare_uhf_roots(error, result, reference, tol, what)
      !! Every root of the reference, against what an unrestricted run reported
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
      call check(error, size(result%excitation_energies) == size(reference), &
                 "the "//what//" run converged a different number of roots than "// &
                 "were asked for")
      if (allocated(error)) return
      do i = 1, size(reference)
         call check(error, abs(result%excitation_energies(i) - reference(i)) < tol, &
                    "a "//what//" excitation energy disagrees with PySCF")
         if (allocated(error)) return
      end do
      call check(error, allocated(result%state_spin), "the "//what//" run labelled "// &
                 "no spins")
      if (allocated(error)) return
      call check(error, all(result%state_spin == STATE_SPIN_UNRESTRICTED), &
                 "an unrestricted root was labelled with a multiplicity it does "// &
                 "not have")
   end subroutine compare_uhf_roots

   subroutine test_oh_uhf_tda_solver(error)
      !! The five lowest unrestricted TDA roots of OH, through the bridge
      type(error_type), allocatable, intent(out) :: error

      type(calculation_result_t) :: result

      call oh_excited_run("", 5, "tda", result)
      call compare_uhf_roots(error, result, OH_UHF_TDA, TOL_CCPVDZ_HF, "OH UHF TDA")
   end subroutine test_oh_uhf_tda_solver

   subroutine test_oh_uhf_rpa_solver(error)
      !! The five lowest unrestricted RPA roots of OH, through the bridge
      !!
      !! The near-zero root the Tamm-Dancoff spectrum carries is not here:
      !! the paired solver puts it at `w^2` below its floor and skips it, so
      !! the five roots asked for are the five physical ones. That is the
      !! difference the plan records between the two columns, and it is
      !! gated here rather than worked around.
      type(error_type), allocatable, intent(out) :: error

      type(calculation_result_t) :: result

      call oh_excited_run("", 5, "rpa", result)
      call compare_uhf_roots(error, result, OH_UHF_RPA, TOL_CCPVDZ_HF, "OH UHF RPA")
   end subroutine test_oh_uhf_rpa_solver

   subroutine test_cation_uks_pbe(error)
      !! The five lowest UKS PBE Tamm-Dancoff roots of the water cation
      type(error_type), allocatable, intent(out) :: error

      type(calculation_result_t) :: result

      if (.not. xc_available()) then
         call check(error, .true.)
         return
      end if
      call cation_excited_run("pbe", 5, "tda", result)
      call compare_uhf_roots(error, result, CATION_PBE_TDA, TOL_GRID, "H2O+ UKS PBE TDA")
   end subroutine test_cation_uks_pbe

   subroutine test_cation_uks_b3lyp(error)
      !! The five lowest UKS B3LYP Tamm-Dancoff roots of the water cation
      type(error_type), allocatable, intent(out) :: error

      type(calculation_result_t) :: result

      if (.not. xc_available()) then
         call check(error, .true.)
         return
      end if
      call cation_excited_run("b3lyp", 5, "tda", result)
      call compare_uhf_roots(error, result, CATION_B3LYP_TDA, TOL_GRID, &
                             "H2O+ UKS B3LYP TDA")
   end subroutine test_cation_uks_b3lyp

   subroutine test_cation_uks_b3lyp_rpa(error)
      !! The five lowest UKS B3LYP RPA roots of the water cation
      !!
      !! The hybrid's paired route, which is the only case here where the
      !! attenuated-free exchange pass, the spin-resolved kernel and the
      !! antisymmetric unrestricted Fock build all run in one product.
      type(error_type), allocatable, intent(out) :: error

      type(calculation_result_t) :: result

      if (.not. xc_available()) then
         call check(error, .true.)
         return
      end if
      call cation_excited_run("b3lyp", 5, "rpa", result)
      call compare_uhf_roots(error, result, CATION_B3LYP_RPA, TOL_GRID, &
                             "H2O+ UKS B3LYP RPA")
   end subroutine test_cation_uks_b3lyp_rpa

   subroutine test_no_beta_electrons(error)
      !! Triplet H2: an unrestricted spectrum out of a reference with no beta
      !!
      !! A high-spin reference can have an empty beta spin, and its alpha
      !! excitations are as well defined as any other open shell's. The
      !! operator used to refuse this outright. What it exercises that nothing
      !! else does is the empty half of every unrestricted quantity: a
      !! zero-length beta block in the trial vector and the diagonal, a beta
      !! response density that is identically zero, a `(n_ao, 0)` orbital
      !! rectangle, and the guards that keep those out of BLAS rather than
      !! calling it with a vanishing inner dimension.
      !!
      !! Measured 2.7e-11 on the worst of the five roots.
      type(error_type), allocatable, intent(out) :: error

      type(calculation_result_t) :: result
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment

      fragment%n_atoms = 2
      fragment%charge = 0
      fragment%multiplicity = 3
      fragment%nelec = 2
      fragment%n_caps = 0
      allocate (fragment%element_numbers(2), fragment%coordinates(3, 2))
      fragment%element_numbers = [1, 1]
      fragment%coordinates = H2_TRIPLET_BOHR

      settings%basis_set = "cc-pvdz"
      settings%functional = ""
      settings%energy_tol = 1.0e-12_dp
      settings%grad_tol = 1.0e-10_dp
      settings%density_tol = 1.0e-10_dp
      settings%max_iter = 300
      settings%excited%enabled = .true.
      settings%excited%n_states = 5
      settings%excited%method = "tda"
      settings%excited%spin = "singlet"
      settings%excited%tolerance = 1.0e-9_dp
      settings%excited%max_iter = 200

      call run_czt_hf(settings, fragment, result)
      call compare_uhf_roots(error, result, H2_TRIPLET_TDA, TOL_CCPVDZ_HF, &
                             "triplet H2 UHF TDA")
   end subroutine test_no_beta_electrons

   subroutine test_uhf_amplitude_norm(error)
      !! `sum_spin(|X|^2 - |Y|^2) = 1` for every unrestricted root
      !!
      !! The unrestricted convention, and the one thing about the amplitudes a
      !! consumer cannot derive for itself. Checked on the paired route, where
      !! it is an identity the solver imposes rather than a property of a unit
      !! vector, and on the Tamm-Dancoff one, where it says the two routes
      !! agree -- which the restricted pair, at 1 and 1/2, do not.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx
      type(error_t) :: err
      real(dp), allocatable :: omega(:), x(:, :), y(:, :)
      integer, allocatable :: spins(:)
      real(dp) :: worst, norm
      integer :: k

      call oh_reference(mol, scf, ctx, err)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "the unrestricted reference failed: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      worst = 0.0_dp
      call response_excitations_uhf(mol, scf%orbitals, scf%orbital_energies, &
                                    scf%n_occupied, scf%orbitals_beta, &
                                    scf%orbital_energies_beta, scf%n_occupied_beta, &
                                    3, "rpa", omega, spins, x, y, err, &
                                    tolerance=1.0e-9_dp)
      if (.not. err%has_error()) then
         do k = 1, size(omega)
            norm = dot_product(x(:, k), x(:, k)) - dot_product(y(:, k), y(:, k))
            worst = max(worst, abs(norm - 1.0_dp))
         end do
      end if
      if (.not. err%has_error()) then
         deallocate (omega, spins, x, y)
         call response_excitations_uhf(mol, scf%orbitals, scf%orbital_energies, &
                                       scf%n_occupied, scf%orbitals_beta, &
                                       scf%orbital_energies_beta, scf%n_occupied_beta, &
                                       3, "tda", omega, spins, x, y, err, &
                                       tolerance=1.0e-9_dp)
         if (.not. err%has_error()) then
            do k = 1, size(omega)
               norm = dot_product(x(:, k), x(:, k)) - dot_product(y(:, k), y(:, k))
               worst = max(worst, abs(norm - 1.0_dp))
            end do
         end if
      end if
      call mol%destroy()
      call check(error,.not. err%has_error(), "an unrestricted solve failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, worst < TOL_PAIRED_NORM, "an unrestricted amplitude is not "// &
                 "normalised to sum_spin(|X|^2 - |Y|^2) = 1")
   end subroutine test_uhf_amplitude_norm

   subroutine restricted_manifolds_case(functional, worst_singlet, worst_triplet, &
                                        error, ok)
      !! The unrestricted operator on a closed shell, against the two restricted ones
      !!
      !! `A_aa + A_ab` is the singlet `A` and `A_aa - A_ab` the triplet one --
      !! Psi4's `test_RU_TDA_C1`, which is the strongest statement available
      !! about an unrestricted response operator without a second code, because
      !! it pins the cross-spin block that no closed-shell test can see.
      !!
      !! **One set of orbitals, not two SCFs.** The restricted and the
      !! unrestricted operator are built from the same converged orbitals, the
      !! second reading them as both spins and half the density as each. So
      !! the molecular-orbital phases are identical by construction and the
      !! two matrices are compared **element by element** rather than through
      !! their eigenvalues -- which is what makes this see the coupling block
      !! at all.
      character(len=*), intent(in) :: functional
         !! Empty is Hartree-Fock.
      real(dp), intent(out) :: worst_singlet, worst_triplet
      type(error_type), allocatable, intent(out) :: error
      logical, intent(out) :: ok

      type(czt_molecule_t), target :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t), target :: ctx, ctx_pol
      type(tda_operator_uhf_t) :: operator
      type(error_t) :: err
      real(dp), allocatable :: singlet(:, :), triplet(:, :), both(:, :)
      real(dp), allocatable :: half(:, :)
      integer :: n_ov
      logical :: kohn_sham

      ok = .false.
      worst_singlet = 0.0_dp
      worst_triplet = 0.0_dp
      kohn_sham = len_trim(functional) > 0

      if (kohn_sham) then
         call water_sto3g(mol, scf, ctx, err, functional=functional)
      else
         call water_sto3g(mol, scf, ctx, err)
      end if
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the closed-shell reference failed: "// &
                    err%get_message())
         call mol%destroy()
         return
      end if

      call dense_tda(mol, scf, ctx, kohn_sham, singlet, err, spin="singlet")
      if (.not. err%has_error()) then
         call dense_tda(mol, scf, ctx, kohn_sham, triplet, err, spin="triplet")
      end if
      if (err%has_error()) then
         call check(error, .false., "a restricted manifold failed: "//err%get_message())
         call mol%destroy()
         return
      end if

      ! The same orbitals as both spins, and half the density as each. A
      ! second, spin-polarised context because libxc fixes the spin channel
      ! when a functional is initialised; same functional, same grid level, so
      ! the quadrature is the same points.
      half = 0.5_dp*scf%density
      if (kohn_sham) then
         call xc_context_create(mol, functional, ctx_pol, err, level=5, &
                                polarized=.true.)
         if (err%has_error()) then
            call check(error, .false., "the polarised context failed: "// &
                       err%get_message())
            call mol%destroy()
            return
         end if
         call build_tda_operator_uhf(mol, scf%orbitals, scf%orbital_energies, &
                                     scf%n_occupied, scf%orbitals, &
                                     scf%orbital_energies, scf%n_occupied, operator, &
                                     err, xc=ctx_pol, ref_a=half, ref_b=half)
      else
         call build_tda_operator_uhf(mol, scf%orbitals, scf%orbital_energies, &
                                     scf%n_occupied, scf%orbitals, &
                                     scf%orbital_energies, scf%n_occupied, operator, err)
      end if
      if (.not. err%has_error()) call tda_dense_matrix_uhf(operator, both, err)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "the unrestricted operator failed: "// &
                    err%get_message())
         return
      end if

      n_ov = size(singlet, 1)
      worst_singlet = maxval(abs(both(1:n_ov, 1:n_ov) &
                                 + both(1:n_ov, n_ov + 1:2*n_ov) - singlet))
      worst_triplet = maxval(abs(both(1:n_ov, 1:n_ov) &
                                 - both(1:n_ov, n_ov + 1:2*n_ov) - triplet))
      ok = .true.
   end subroutine restricted_manifolds_case

   subroutine test_restricted_from_unrestricted_hf(error)
      !! Hartree-Fock: `A_aa +/- A_ab` is the singlet and triplet `A`
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst_singlet, worst_triplet
      logical :: ok

      call restricted_manifolds_case("", worst_singlet, worst_triplet, error, ok)
      if (allocated(error) .or. .not. ok) return
      call check(error, worst_singlet < TOL_EXACT, "the unrestricted operator's "// &
                 "spin sum is not the restricted singlet A")
      if (allocated(error)) return
      call check(error, worst_triplet < TOL_EXACT, "the unrestricted operator's "// &
                 "spin difference is not the restricted triplet A")
   end subroutine test_restricted_from_unrestricted_hf

   subroutine test_restricted_from_unrestricted_pbe(error)
      !! PBE: the same, and the statement that the polarised kernel is right
      !!
      !! Looser than the Hartree-Fock case by an order, and the reason is
      !! libxc: the restricted side evaluates the unpolarised functional's
      !! second derivative and this side evaluates the polarised one at
      !! `rho_a = rho_b`, which are the same number computed two ways.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst_singlet, worst_triplet
      logical :: ok

      if (.not. xc_available()) then
         call check(error, .true.)
         return
      end if
      call restricted_manifolds_case("pbe", worst_singlet, worst_triplet, error, ok)
      if (allocated(error) .or. .not. ok) return
      call check(error, worst_singlet < TOL_UKS_MANIFOLD, "the unrestricted PBE "// &
                 "operator's spin sum is not the restricted singlet A")
      if (allocated(error)) return
      call check(error, worst_triplet < TOL_UKS_MANIFOLD, "the unrestricted PBE "// &
                 "operator's spin difference is not the restricted triplet A")
   end subroutine test_restricted_from_unrestricted_pbe

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
