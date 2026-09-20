!! Tamm-Dancoff excitation energies, singlet and triplet, against PySCF
module test_mqc_czt_tddft
   !! The Layer 2 and Layer 3 gates of `TDDFT_PLAN.md`: the TDA operator `A`
   !! itself, singlet and triplet, and the excitation energies a Davidson
   !! finds in each.
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
   use mqc_czt_tddft, only: tda_operator_t, build_tda_operator, tda_dense_matrix
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_czt_bridge, only: run_czt_hf
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t, STATE_SPIN_SINGLET, &
                               STATE_SPIN_TRIPLET
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

   !! H2O/cc-pVDZ, the triplet TDA roots, in Hartree, to twelve decimals.
   !!
   !! Regenerated rather than transcribed, as the singlets above were; the
   !! plan's ten-digit triplet table agrees with these to 1.3e-9, which is
   !! where it stops saying anything. Same references, same geometry, same
   !! `grids.level = 5`.
   real(dp), parameter :: RHF_CCPVDZ_TRIPLET(N_CCPVDZ_TRIPLETS) = [ &
                          0.304752995712_dp, 0.382448306917_dp, 0.383738345685_dp]
   real(dp), parameter :: PBE_CCPVDZ_TRIPLET(N_CCPVDZ_TRIPLETS) = [ &
                          0.245640585460_dp, 0.321365211218_dp, 0.322245501451_dp]
   real(dp), parameter :: B3LYP_CCPVDZ_TRIPLET(N_CCPVDZ_TRIPLETS) = [ &
                          0.254376470706_dp, 0.331305647032_dp, 0.331912037086_dp]
   real(dp), parameter :: CAM_CCPVDZ_TRIPLET(N_CCPVDZ_TRIPLETS) = [ &
                          0.256574750592_dp, 0.335035697057_dp, 0.336081588546_dp]

   integer, parameter :: N_BOTH_STATES = 3
      !! Roots **per manifold** asked for by the `spin = "both"` case, so six
      !! come back.

   !! The union of the two H2O/STO-3G manifolds, as `spin = "both"` reports
   !! it: the three lowest of each, sorted together by energy.
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

   !! H2 stretched to 3.0 Angstrom, in Bohr, whose restricted Hartree-Fock
   !! solution is a saddle point against spin polarisation.
   !!
   !! At this separation the closed shell is the textbook wrong answer: the
   !! lowest triplet TDA root of PySCF's own operator is -0.172235168077
   !! hartree, and every root above it is an expansion about that saddle
   !! point. 6-31G rather than STO-3G so the occupied-virtual space is three
   !! dimensional and the Davidson has somewhere to iterate; the minimal basis
   !! gives a one-by-one problem, which tests the message and not the solver.
   real(dp), parameter :: H2_STRETCHED_BOHR(3, 2) = reshape([ &
                                                            0.0_dp, 0.0_dp, 0.0_dp, &
                                                            0.0_dp, 0.0_dp, 5.6691783763734843_dp], [3, 2])
   real(dp), parameter :: H2_STRETCHED_ENERGY = -0.815591795493_dp

   !! Hartree-Fock carries no quadrature, so nothing but the integrals and the
   !! two SCF thresholds sits between the codes, and both are converged far
   !! below this. Measured on the STO-3G matrix: 7.7e-13 on the worst diagonal
   !! element, 6.0e-13 off it, 1.3e-12 on the orbital energies and 1.2e-12 on
   !! the spectrum. The bound is two orders over that, which leaves room for
   !! the compiler-to-compiler spread and still sits eight orders under the
   !! smallest fault it has to catch: halving the kernel, dropping one
   !! exchange term or symmetrising the wrong way each move elements by 1e-2
   !! and up. The triplet matrix lands in the same place, 6.1e-13 on the
   !! diagonal and 3.8e-13 off it, but only once the reference's own SCF is
   !! driven to `conv_tol = 1e-15`: at PySCF's ordinary 1e-13 the core-excited
   !! diagonal elements near 20 hartree still carry 7e-11, which is most of
   !! this bound spent on the reference rather than on the operator.
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
   !! one to watch, at a sixtieth of this bound. 1e-7 clears the worst of them
   !! by that margin and leaves the compiler spread room, which the
   !! double-hybrid Hessian shows can reach 2e-9 on a quadrature of this kind.
   !! Tightening it to the measured numbers would be re-recording a pin to
   !! make it pass on one compiler.
   !!
   !! The triplet cases land in the same band: 3.5e-10 on the STO-3G PBE
   !! matrix and its spectrum, and on cc-pVDZ 1.8e-10 for PBE, 5.6e-10 for
   !! B3LYP and 1.9e-10 for CAM-B3LYP. The polarised kernel is evaluated on
   !! the same points as the unpolarised one, so it adds no quadrature error
   !! of its own.
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
   !! 3.1e-12 on the worst of the five singlets and 2.2e-12 on the worst of
   !! the three triplets, so the margin is four orders.
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
                  new_unittest("the_rhf_triplet_tda_matrix_matches_pyscf", &
                               test_rhf_triplet_matrix), &
                  new_unittest("the_pbe_triplet_tda_matrix_matches_pyscf", &
                               test_pbe_triplet_matrix), &
                  new_unittest("cc_pvdz_rhf_triplets_match_pyscf", test_ccpvdz_rhf_t), &
                  new_unittest("cc_pvdz_pbe_triplets_match_pyscf", test_ccpvdz_pbe_t), &
                  new_unittest("cc_pvdz_b3lyp_triplets_match_pyscf", test_ccpvdz_b3lyp_t), &
                  new_unittest("cc_pvdz_cam_b3lyp_triplets_match_pyscf", test_ccpvdz_cam_t), &
                  new_unittest("both_manifolds_interleave_by_energy", test_both_spins), &
                  new_unittest("a_triplet_unstable_reference_is_named", test_instability), &
                  new_unittest("no_states_asked_for_means_no_spectrum", test_no_states), &
                  new_unittest("rpa_is_refused_by_name", test_later_layers) &
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
                    "a "//what//" TDA excitation energy disagrees with the table")
         if (allocated(error)) return
      end do

      call check(error, allocated(result%state_spin), "the "//what//" run labelled "// &
                 "no spins")
      if (allocated(error)) return
      call check(error, all(result%state_spin == want_spin), &
                 "a single-manifold solve reported a root labelled with the other spin")
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

   subroutine test_rhf_triplet_matrix(error)
      !! The ten-by-ten triplet TDA matrix of H2O/STO-3G, element by element
      !!
      !! Two independent changes separate this matrix from the singlet one
      !! above, and each on its own is worth a tenth of a Hartree on the
      !! diagonal: the Coulomb term is gone, and the kernel -- absent here,
      !! since this is Hartree-Fock -- would be the spin difference. So for
      !! Hartree-Fock this case is exactly the `j_scale = 0` half of Layer 3,
      !! isolated from the kernel, and the PBE case below is the other half.
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
   end subroutine test_rhf_triplet_matrix

   subroutine test_pbe_triplet_matrix(error)
      !! The same, for a pure functional, where the triplet kernel is the
      !! whole of the coupling
      !!
      !! PBE has no exact exchange and a triplet has no Coulomb term, so every
      !! element off the orbital-energy diagonal here comes from
      !! `(f_aa - f_ab)/2` and nothing else. A triplet kernel evaluated
      !! unpolarised -- that is, the singlet one left in place -- moves the
      !! first root by 0.07 hartree, and the four GGA combinations mis-weighted
      !! against each other move it by less but never by less than this bound.
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
      call check(error, ok, "the dense diagonalisation of the PBE triplet matrix failed")
      if (allocated(error)) return
      call check(error, maxval(abs(values - PBE_TRIPLET_A_EIG)) < TOL_GRID, &
                 "the PBE triplet TDA spectrum disagrees with PySCF")
   end subroutine test_pbe_triplet_matrix

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
      !! The only case here that exercises the attenuated exchange build and
      !! the polarised kernel together. Each of them is separately visible in
      !! an earlier case, so a failure that appears only here is the
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

   subroutine h2_stretched_run(n_states, spin, result)
      !! H2 at 3.0 Angstrom through the bridge, restricted Hartree-Fock
      integer, intent(in) :: n_states
      character(len=*), intent(in) :: spin
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
      fragment%coordinates = H2_STRETCHED_BOHR

      settings%basis_set = "6-31g"
      settings%functional = ""
      settings%energy_tol = 1.0e-12_dp
      settings%grad_tol = 1.0e-9_dp
      settings%density_tol = 1.0e-9_dp
      settings%max_iter = 200
      settings%excited%enabled = n_states > 0
      settings%excited%n_states = n_states
      settings%excited%method = "tda"
      settings%excited%spin = spin
      settings%excited%tolerance = 1.0e-8_dp
      settings%excited%max_iter = 200

      call run_czt_hf(settings, fragment, result)
   end subroutine h2_stretched_run

   subroutine test_instability(error)
      !! A triplet-unstable reference is named, not reported as a number
      !!
      !! Stretched H2 is where the restricted solution stops being a minimum:
      !! its lowest triplet TDA root is -0.17 hartree. Two wrong answers are
      !! possible and both look plausible from outside -- reporting the
      !! negative root as an excitation, or silently dropping it under the
      !! floor and handing back the two above it as though a state were merely
      !! missing. The singlet manifold of the same reference is well behaved,
      !! which is what the second half of this checks: the diagnosis has to be
      !! about the triplet operator and not about the molecule being awkward.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result
      character(len=:), allocatable :: message

      call h2_stretched_run(1, "triplet", result)
      call check(error, result%has_error, "a triplet-unstable reference produced a "// &
                 "spectrum rather than a diagnosis")
      if (allocated(error)) return
      message = result%error%get_message()
      call check(error, index(message, "triplet-unstable") > 0, &
                 "the failure of a triplet-unstable reference was reported as "// &
                 "something else: "//message)
      if (allocated(error)) return
      call check(error,.not. result%has_excited_states, "a refused triplet-unstable "// &
                 "run still reported excited states")
      if (allocated(error)) return

      ! The same molecule, same reference, singlet manifold: an ordinary
      ! answer. Without this the case would pass just as well for a build
      ! that refused every excited-state run on stretched H2.
      call h2_stretched_run(1, "singlet", result)
      call check(error,.not. result%has_error, "the singlet manifold of the same "// &
                 "reference failed too: "//result%error%get_message())
      if (allocated(error)) return
      call check(error, abs(result%energy%scf - H2_STRETCHED_ENERGY) < TOL_CCPVDZ_HF, &
                 "the stretched-H2 Hartree-Fock energy is not PySCF's")
      if (allocated(error)) return
      call check(error, result%has_excited_states, "the singlet manifold of a "// &
                 "triplet-unstable reference reported no excited states")
   end subroutine test_instability

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
      !! What is still unimplemented is refused by name, not approximated
      !!
      !! A TDA number handed back for an RPA deck is the failure worth
      !! guarding against: it converges, it is of the right magnitude, and
      !! nothing in the output says which problem was solved.
      !!
      !! This case checked the triplet refusal too until Layer 3 removed it.
      !! What replaces that half is not a weaker assertion but a stronger one
      !! -- the four cc-pVDZ triplet cases and the two triplet matrices above
      !! -- so nothing is lost by the refusal going away.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      call excited_run("sto-3g", "", 2, "rpa", "singlet", result)
      call check(error, result%has_error, "an RPA deck was answered rather than "// &
                 "refused")
      if (allocated(error)) return
      call check(error,.not. result%has_excited_states, "a refused RPA deck still "// &
                 "reported excited states")
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
