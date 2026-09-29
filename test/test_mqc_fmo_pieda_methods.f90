!! PIEDA beyond Hartree-Fock: Kohn-Sham and the MP2 family
module test_mqc_fmo_pieda_methods
   !! What PIEDA is held to once the fragments are not plain Hartree-Fock:
   !!
   !! 1. the terms close the pair energy, for a GGA, a hybrid, a
   !!    range-separated hybrid, MP2, SCS-MP2 and RI-MP2;
   !! 2. an MP2 run's Ees and Eex are the Hartree-Fock run's, which is what
   !!    says Eex is built from the Hartree-Fock internal energies and not from
   !!    ones with the monomers' correlation in them;
   !! 3. an MP2 run's Edi is `Ec(IJ) - Ec(I) - Ec(J)` of separately computed
   !!    correlation energies, spin-component scaling included;
   !! 4. the Kohn-Sham energy of the union state is the SCF's own functional:
   !!    at a single fragment's density it is that fragment's SCF energy;
   !! 5. Kohn-Sham Eex vanishes with distance and is repulsive at contact;
   !! 6. `edi_in_energy`, and the empirical Edi under a functional;
   !! 7. a functional next to a detached bond, in both `pieda_hl` modes.
   !!
   !! Kohn-Sham cases return at once on a build without libxc, as
   !! `test_mqc_fmo_dft` does.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: to_bohr
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_rhf, only: run_czt_rhf, rhf_result_t
   use mqc_czt_mp2, only: mp2_result_t, run_czt_mp2, run_czt_ri_mp2
   use mqc_czt_xc, only: xc_context_t, xc_available
   use mqc_czt_pieda, only: cholesky_occupied_orbitals, hl_density_from_orbitals, &
                            hl_prime_energy
   use mqc_czt_fragment_solver, only: fragment_xc_context
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_method_config, only: correlation_config_t
   use mqc_elements, only: core_orbital_count
   use mqc_dispersion_apply, only: dispersion_apply, dispersion_kind_available
   implicit none
   private

   public :: collect_mqc_fmo_pieda_methods

   real(dp), parameter :: SUM_TOL = 1.0e-12_dp
      !! The closing term is a residual, so the sum is exact to round-off.
   real(dp), parameter :: HF_TOL = 1.0e-10_dp
      !! MP2 against Hartree-Fock: both take the same embedded densities and
      !! differ only in where the self-consistent field stopped.
   real(dp), parameter :: EC_TOL = 1.0e-10_dp
      !! Edi against correlation energies computed apart: the same
      !! Hartree-Fock reference, solved twice to the same tolerance.
   real(dp), parameter :: HL_TOL = 1.0e-9_dp
   real(dp), parameter :: SCF_GRAD_TOL = 1.0e-10_dp
      !! The commutator bound of the SCFs an MP2 energy is compared across: the
      !! correlation energy is not variational in the orbitals, so it inherits
      !! their error at first order.

   integer, parameter :: GRID_LEVEL = 3
   character(len=*), parameter :: AUX = "cc-pvdz-rifit"

contains

   subroutine collect_mqc_fmo_pieda_methods(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("kohn_sham_terms_close_every_pair_energy", test_ks_sum), &
                  new_unittest("mp2_family_terms_close_every_pair_energy", test_mp2_sum), &
                  new_unittest("mp2_pieda_has_the_hartree_fock_ees_and_eex", test_mp2_hf_terms), &
                  new_unittest("mp2_edi_is_the_correlation_interaction", test_mp2_edi), &
                  new_unittest("scs_mp2_edi_is_scaled", test_scs_edi), &
                  new_unittest("the_union_energy_of_one_fragment_is_its_kohn_sham_energy", &
                               test_hl_functional), &
                  new_unittest("kohn_sham_eex_vanishes_with_distance", test_ks_eex_far), &
                  new_unittest("kohn_sham_eex_is_repulsive_at_contact", test_ks_eex_contact), &
                  new_unittest("edi_in_energy_follows_the_method", test_edi_in_energy), &
                  new_unittest("empirical_edi_uses_the_functionals_parameters", test_functional_edi), &
                  new_unittest("mp2_with_empirical_dispersion_is_refused", test_mp2_dispersion), &
                  new_unittest("kohn_sham_pieda_next_to_a_cut", test_ks_cut) &
                  ]
   end subroutine collect_mqc_fmo_pieda_methods

   subroutine test_ks_sum(error)
      !! Kohn-Sham: a GGA, a global hybrid, a range-separated hybrid
      type(error_type), allocatable, intent(out) :: error

      character(len=9), parameter :: FUNCTIONALS(3) = ["pbe      ", "pbe0     ", "cam-b3lyp"]
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(error_t) :: err
      integer :: z(6), k
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water_dimer(z, sym, xyz, 0.0_dp)
      do k = 1, size(FUNCTIONALS)
         call pieda_options(opts, trim(FUNCTIONALS(k)))
         call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
         call check(error,.not. err%has_error(), trim(FUNCTIONALS(k))//" PIEDA run failed: "// &
                    err%get_message())
         if (allocated(error)) return
         call check_closure(error, res, trim(FUNCTIONALS(k)))
         if (allocated(error)) return
      end do
   end subroutine test_ks_sum

   subroutine test_mp2_sum(error)
      !! MP2 family: the four terms close a trimer's three pairs
      type(error_type), allocatable, intent(out) :: error

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(error_t) :: err
      integer :: z(9), k
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      character(len=6), parameter :: NAMES(3) = ["mp2   ", "scs   ", "ri-mp2"]

      call water_trimer(z, sym, xyz)
      do k = 1, 3
         call pieda_options(opts, "")
         call mp2_settings(opts, use_ri=k == 3, use_scs=k == 2)
         call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
         call check(error,.not. err%has_error(), trim(NAMES(k))//" PIEDA run failed: "// &
                    err%get_message())
         if (allocated(error)) return
         call check(error, size(res%pairs), 3, "the trimer should have three pairs")
         if (allocated(error)) return
         call check_closure(error, res, trim(NAMES(k)))
         if (allocated(error)) return
      end do
   end subroutine test_mp2_sum

   subroutine test_mp2_hf_terms(error)
      !! Ees and Eex do not see the correlation
      !!
      !! If Eex took the monomers' `energy` with their correlation in it, it
      !! would be off by the sum of their `Ec`, some 1e-2 Hartree.
      type(error_type), allocatable, intent(out) :: error

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: hf, mp2
      type(error_t) :: err
      integer :: z(9), p
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)

      call water_trimer(z, sym, xyz)
      call pieda_options(opts, "")
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, hf, err)
      call check(error,.not. err%has_error(), "the Hartree-Fock run failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call mp2_settings(opts, use_ri=.false., use_scs=.false.)
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, mp2, err)
      call check(error,.not. err%has_error(), "the MP2 run failed: "//err%get_message())
      if (allocated(error)) return

      do p = 1, size(hf%pairs)
         call check(error, abs(mp2%pairs(p)%ees - hf%pairs(p)%ees) < HF_TOL, &
                    "MP2 Ees is not the Hartree-Fock Ees")
         if (allocated(error)) return
         call check(error, abs(mp2%pairs(p)%eex - hf%pairs(p)%eex) < HF_TOL, &
                    "MP2 Eex is not the Hartree-Fock Eex")
         if (allocated(error)) then
            write (*, *) "   pair", p, " hf", hf%pairs(p)%eex, " mp2", mp2%pairs(p)%eex
            return
         end if
         ! The correlation part of the pair's energy is all that MP2 adds to it.
         call check(error, abs((mp2%pairs(p)%energy - hf%pairs(p)%energy) &
                               - mp2%pairs(p)%edi) < HF_TOL, &
                    "MP2 Edi is not the MP2 pair energy less the Hartree-Fock one")
         if (allocated(error)) return
      end do
   end subroutine test_mp2_hf_terms

   subroutine test_mp2_edi(error)
      type(error_type), allocatable, intent(out) :: error

      call edi_against_separate_correlation(error, use_ri=.false., use_scs=.false.)
      if (allocated(error)) return
      call edi_against_separate_correlation(error, use_ri=.true., use_scs=.false.)
   end subroutine test_mp2_edi

   subroutine test_scs_edi(error)
      type(error_type), allocatable, intent(out) :: error

      call edi_against_separate_correlation(error, use_ri=.false., use_scs=.true.)
   end subroutine test_scs_edi

   subroutine edi_against_separate_correlation(error, use_ri, use_scs)
      !! `Edi = Ec(IJ) - Ec(I) - Ec(J)`
      !!
      !! With no field (`esp = "none"`) each term of the expansion is an
      !! isolated molecule, so its correlation energy is what an ordinary
      !! Hartree-Fock-plus-MP2 calculation on those atoms reports.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: use_ri, use_scs

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: ec_ij, ec_i, ec_j

      call water_dimer(z, sym, xyz, 0.0_dp)
      call pieda_options(opts, "")
      opts%esp = "none"
      opts%scf%grad_tol = SCF_GRAD_TOL
      call mp2_settings(opts, use_ri=use_ri, use_scs=use_scs)
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the MP2 PIEDA run failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call correlation_energy(z, sym, xyz, use_ri, use_scs, ec_ij, err)
      if (.not. err%has_error()) call correlation_energy(z(1:3), sym(1:3), xyz(:, 1:3), &
                                                         use_ri, use_scs, ec_i, err)
      if (.not. err%has_error()) call correlation_energy(z(4:6), sym(4:6), xyz(:, 4:6), &
                                                         use_ri, use_scs, ec_j, err)
      call check(error,.not. err%has_error(), "the separate correlation energies failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call check(error, res%pairs(1)%pieda, "the pair was not decomposed")
      if (allocated(error)) return
      call check(error, abs(res%pairs(1)%edi - (ec_ij - ec_i - ec_j)) < EC_TOL, &
                 "Edi is not Ec(IJ) - Ec(I) - Ec(J)")
      if (allocated(error)) then
         write (*, *) "   edi   =", res%pairs(1)%edi
         write (*, *) "   ec sum=", ec_ij - ec_i - ec_j
         return
      end if
      call check(error, res%pairs(1)%edi < 0.0_dp, "the correlation interaction should attract")
   end subroutine edi_against_separate_correlation

   subroutine test_hl_functional(error)
      !! the union state's energy functional is the SCF's
      !!
      !! One fragment's `D_HL` is its own density, so evaluating the functional
      !! at it has to return the SCF energy -- exchange-correlation, the
      !! exact-exchange fraction, the range separation and a meta-GGA's kinetic
      !! energy density all included. The Hartree-Fock functional at the same
      !! density is far from it.
      type(error_type), allocatable, intent(out) :: error

      character(len=9), parameter :: FUNCTIONALS(4) = ["pbe      ", "b3lyp    ", "cam-b3lyp", "tpss     "]
      type(cuest_scf_settings_t) :: method
      type(czt_molecule_t) :: mol
      type(xc_context_t) :: xc
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: bounds(:, :), c(:, :), s(:, :), d_hl(:, :)
      real(dp) :: e_ks, e_hf
      integer :: z(3), k
      character(len=2) :: sym(3)
      real(dp) :: xyz(3, 3)

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water(z, sym, xyz)
      do k = 1, size(FUNCTIONALS)
         method%functional = trim(FUNCTIONALS(k))
         method%grid_level = GRID_LEVEL
         call build_czt_molecule(z, sym, xyz, "6-31g", mol, err)
         if (.not. err%has_error()) call fragment_xc_context(method, mol, xc, err)
         if (.not. err%has_error()) call run_czt_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, &
                                                     .false., scf, err, xc=xc)
         if (.not. err%has_error()) call schwarz_bounds(mol, bounds, err)
         if (.not. err%has_error()) call cholesky_occupied_orbitals(scf%density, 5, "water", &
                                                                    c, err)
         if (.not. err%has_error()) then
            call mol%overlap(s)
            call hl_density_from_orbitals(c, s, d_hl, err)
         end if
         if (.not. err%has_error()) call hl_prime_energy(mol, bounds, d_hl, e_ks, err, xc=xc)
         if (.not. err%has_error()) call hl_prime_energy(mol, bounds, d_hl, e_hf, err)
         call check(error,.not. err%has_error(), trim(FUNCTIONALS(k))//" failed: "// &
                    err%get_message())
         if (allocated(error)) return
         call xc%destroy()

         call check(error, abs(e_ks - scf%energy) < HL_TOL, &
                    "the union energy at one fragment's density is not its "// &
                    trim(FUNCTIONALS(k))//" energy")
         if (allocated(error)) then
            write (*, *) "   hl energy =", e_ks
            write (*, *) "   scf energy=", scf%energy
            return
         end if
         call check(error, abs(e_hf - scf%energy) > 1.0e-3_dp, &
                    "the Hartree-Fock functional should not agree with "// &
                    trim(FUNCTIONALS(k)))
         if (allocated(error)) return
      end do
   end subroutine test_hl_functional

   subroutine test_ks_eex_far(error)
      !! no overlap, no exchange
      type(error_type), allocatable, intent(out) :: error

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water_dimer(z, sym, xyz, 5.1_dp)
      call pieda_options(opts, "pbe")
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the separated PBE run failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, res%pairs(1)%distance > 6.0_dp, "the dimer is not pulled apart")
      if (allocated(error)) return
      call check(error, abs(res%pairs(1)%eex) < 1.0e-6_dp, &
                 "PBE Eex does not vanish for a dimer pulled far apart")
      if (allocated(error)) write (*, *) "   eex =", res%pairs(1)%eex
   end subroutine test_ks_eex_far

   subroutine test_ks_eex_contact(error)
      !! at the hydrogen-bonded geometry the exchange repels, as it
      !! does at Hartree-Fock
      type(error_type), allocatable, intent(out) :: error

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res, hf
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water_dimer(z, sym, xyz, 0.0_dp)
      call pieda_options(opts, "pbe")
      opts%basis = "6-31g"
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the PBE run failed: "//err%get_message())
      if (allocated(error)) return
      call pieda_options(opts, "")
      opts%basis = "6-31g"
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, hf, err)
      call check(error,.not. err%has_error(), "the Hartree-Fock run failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call check(error, res%pairs(1)%eex > 0.0_dp, "PBE Eex is not repulsive at contact")
      if (allocated(error)) then
         write (*, *) "   eex =", res%pairs(1)%eex
         return
      end if
      call check(error, res%pairs(1)%eex > 0.5_dp*hf%pairs(1)%eex .and. &
                 res%pairs(1)%eex < 2.0_dp*hf%pairs(1)%eex, &
                 "PBE Eex is not of the size of Hartree-Fock's")
      if (allocated(error)) then
         write (*, *) "   pbe eex =", res%pairs(1)%eex
         write (*, *) "   hf eex  =", hf%pairs(1)%eex
      end if
   end subroutine test_ks_eex_contact

   subroutine test_edi_in_energy(error)
      !! which Edi is inside the pair energy
      type(error_type), allocatable, intent(out) :: error

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      character(len=4) :: disp

      call water_dimer(z, sym, xyz, 0.0_dp)

      call pieda_options(opts, "")
      call mp2_settings(opts, use_ri=.false., use_scs=.false.)
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the MP2 run failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, res%edi_in_energy, "MP2's Edi should be inside the pair energy")
      if (allocated(error)) return

      call pieda_options(opts, "")
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the Hartree-Fock run failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error,.not. res%edi_in_energy, "Hartree-Fock has no Edi in its energy")
      if (allocated(error)) return

      disp = "d4"
      if (.not. dispersion_kind_available("d4")) disp = "d3bj"
      if (.not. dispersion_kind_available(trim(disp))) return

      opts%pieda_dispersion = disp
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the Hartree-Fock run with dispersion failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error,.not. res%edi_in_energy, &
                 "the empirical Edi is not inside the pair energy")
      if (allocated(error)) return
      call check(error, res%pairs(1)%edi < 0.0_dp, "the empirical Edi should attract")
      if (allocated(error)) return
      call check(error, abs(res%pairs(1)%ees + res%pairs(1)%eex + res%pairs(1)%ect_mix &
                            - res%pairs(1)%energy) < SUM_TOL, &
                 "the empirical Edi is in the sum that closes the pair energy")
   end subroutine test_edi_in_energy

   subroutine test_functional_edi(error)
      !! Empirical Edi under a functional takes that functional's damping
      type(error_type), allocatable, intent(out) :: error

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res, hf
      type(error_t) :: err, derr
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: e_ij, e_i, e_j
      character(len=4) :: disp

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check
      disp = "d4"
      if (.not. dispersion_kind_available("d4")) disp = "d3bj"
      if (.not. dispersion_kind_available(trim(disp))) return

      call water_dimer(z, sym, xyz, 0.0_dp)
      call pieda_options(opts, "pbe")
      opts%pieda_dispersion = disp
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the PBE run with dispersion failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call dispersion_apply(trim(disp), "pbe", 0.0_dp, z, xyz, e_ij, error=derr)
      if (.not. derr%has_error()) call dispersion_apply(trim(disp), "pbe", 0.0_dp, z(1:3), &
                                                        xyz(:, 1:3), e_i, error=derr)
      if (.not. derr%has_error()) call dispersion_apply(trim(disp), "pbe", 0.0_dp, z(4:6), &
                                                        xyz(:, 4:6), e_j, error=derr)
      call check(error,.not. derr%has_error(), "the PBE dispersion failed: "// &
                 derr%get_message())
      if (allocated(error)) return
      call check(error, abs(res%pairs(1)%edi - (e_ij - e_i - e_j)) < 1.0e-12_dp, &
                 "Edi under PBE is not the PBE-damped dispersion interaction")
      if (allocated(error)) return

      call pieda_options(opts, "")
      opts%pieda_dispersion = disp
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, hf, err)
      call check(error,.not. err%has_error(), "the Hartree-Fock run with dispersion failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, abs(res%pairs(1)%edi - hf%pairs(1)%edi) > 1.0e-6_dp, &
                 "Edi under PBE is the Hartree-Fock one")
   end subroutine test_functional_edi

   subroutine test_mp2_dispersion(error)
      !! The backend refuses what `fragment_refusal` refuses, by name
      type(error_type), allocatable, intent(out) :: error

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)

      call water_dimer(z, sym, xyz, 0.0_dp)
      call pieda_options(opts, "")
      call mp2_settings(opts, use_ri=.false., use_scs=.false.)
      opts%pieda_dispersion = "d4"
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error, err%has_error(), "MP2 with pieda_dispersion was not refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "pieda_dispersion") > 0, &
                 "the refusal does not name the key")
   end subroutine test_mp2_dispersion

   subroutine test_ks_cut(error)
      !! butane cut into two ethyls, with a water. The two pairs that
      !! are not themselves joined by the bond are decomposed; the union state
      !! only differs between the modes by the frozen virtuals it holds, and
      !! the functional does not see them.
      type(error_type), allocatable, intent(out) :: error

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: gamess, projected
      type(error_t) :: err
      integer :: z(17), p
      character(len=2) :: sym(17)
      real(dp) :: xyz(3, 17)

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call butane_water(z, sym, xyz)
      call pieda_options(opts, "pbe")
      opts%bond_breaking = "afo"
      opts%afo_localization = "er"
      opts%resppc = 2.0_dp
      opts%scf_energy_tol = 1.0e-10_dp
      opts%scf_density_tol = 1.0e-8_dp
      opts%outer_tol = 1.0e-9_dp
      call run_fmo2(z, sym, xyz, [1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3], &
                    opts, gamess, err)
      call check(error,.not. err%has_error(), "the PBE cut run failed: "//err%get_message())
      if (allocated(error)) return
      opts%pieda_hl = "projected"
      call run_fmo2(z, sym, xyz, [1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3], &
                    opts, projected, err)
      call check(error,.not. err%has_error(), "the projected PBE cut run failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call check(error, count(gamess%pairs%pieda), 2, &
                 "the two pairs the bond does not join should be decomposed")
      if (allocated(error)) return
      do p = 1, size(gamess%pairs)
         if (.not. gamess%pairs(p)%pieda) cycle
         call check(error, abs(gamess%pairs(p)%eex - projected%pairs(p)%eex) < 1.0e-4_dp, &
                    "the two pieda_hl modes disagree by more than the frozen virtuals allow")
         if (allocated(error)) then
            write (*, *) "   gamess   =", gamess%pairs(p)%eex
            write (*, *) "   projected=", projected%pairs(p)%eex
            return
         end if
      end do
      do p = 1, size(gamess%pairs)
         if (.not. gamess%pairs(p)%pieda) cycle
         call check(error, abs(gamess%pairs(p)%ees + gamess%pairs(p)%eex &
                               + gamess%pairs(p)%ect_mix - gamess%pairs(p)%energy) < SUM_TOL, &
                    "the terms do not close the energy of a pair next to a cut")
         if (allocated(error)) return
      end do
   end subroutine test_ks_cut

   subroutine check_closure(error, res, label)
      !! Every decomposed pair: `Ees + Eex + Ect+mix (+ Edi) = dE_IJ`
      type(error_type), allocatable, intent(out) :: error
      type(fmo_result_t), intent(in) :: res
      character(len=*), intent(in) :: label

      real(dp) :: total
      integer :: p

      do p = 1, size(res%pairs)
         call check(error, res%pairs(p)%pieda, label//": a pair was not decomposed")
         if (allocated(error)) return
         total = res%pairs(p)%ees + res%pairs(p)%eex + res%pairs(p)%ect_mix
         if (res%edi_in_energy) total = total + res%pairs(p)%edi
         call check(error, abs(total - res%pairs(p)%energy) < SUM_TOL, &
                    label//": the four terms do not close the pair energy")
         if (allocated(error)) then
            write (*, *) "   sum   =", total
            write (*, *) "   energy=", res%pairs(p)%energy
            return
         end if
      end do
   end subroutine check_closure

   subroutine correlation_energy(z, sym, xyz, use_ri, use_scs, ec, error)
      !! The MP2 correlation energy of a whole molecule, as an ordinary
      !! Hartree-Fock-plus-MP2 calculation gives it, scaled if asked
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz(:, :)
      logical, intent(in) :: use_ri, use_scs
      real(dp), intent(out) :: ec
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol, aux_mol
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2
      type(correlation_config_t) :: default_corr
      real(dp) :: ss, os
      integer :: nelec, frozen

      ec = 0.0_dp
      ss = 1.0_dp
      os = 1.0_dp
      if (use_scs) then
         ss = default_corr%scs_ss
         os = default_corr%scs_os
      end if
      nelec = sum(z)
      frozen = core_orbital_count(z)

      call build_czt_molecule(z, sym, xyz, "sto-3g", mol, error)
      if (error%has_error()) return
      call run_czt_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error, &
                       grad_tol=SCF_GRAD_TOL)
      if (error%has_error()) return
      if (use_ri) then
         call build_czt_molecule(z, sym, xyz, AUX, aux_mol, error)
         if (error%has_error()) return
         call run_czt_ri_mp2(mol, aux_mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                             scf%energy, mp2, error, n_frozen=frozen)
         call aux_mol%destroy()
      else
         call run_czt_mp2(mol, scf%orbitals, scf%orbital_energies, nelec/2, scf%energy, &
                          mp2, error, n_frozen=frozen)
      end if
      if (error%has_error()) return
      ec = ss*mp2%same_spin + os*mp2%opposite_spin
   end subroutine correlation_energy

   subroutine pieda_options(opts, functional)
      !! PIEDA on, STO-3G, the exact field, every pair solved, tight tolerances;
      !! Hartree-Fock when `functional` is empty
      type(fmo_options_t), intent(out) :: opts
      character(len=*), intent(in) :: functional

      opts%basis = "sto-3g"
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%resppc = -1.0_dp
      opts%resdim = 0.0_dp
      opts%outer_tol = 1.0e-12_dp
      opts%max_outer = 100
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-12_dp
      opts%scf_density_tol = 1.0e-10_dp
      opts%pieda = .true.
      if (len_trim(functional) > 0) then
         opts%method%functional = functional
         opts%method%grid_level = GRID_LEVEL
      end if
   end subroutine pieda_options

   subroutine mp2_settings(opts, use_ri, use_scs)
      !! MP2, or one of its variants, as the fragment method, frozen core
      type(fmo_options_t), intent(inout) :: opts
      logical, intent(in) :: use_ri, use_scs

      type(correlation_config_t) :: default_corr

      opts%method%run_mp2 = .true.
      opts%method%corr_density_fitting = use_ri
      opts%method%aux_basis_set = AUX
      opts%method%freeze_core = .true.
      if (use_scs) then
         opts%method%scs_ss = default_corr%scs_ss
         opts%method%scs_os = default_corr%scs_os
      end if
   end subroutine mp2_settings

   subroutine water(z, sym, xyz)
      !! One water, Angstrom converted to Bohr
      integer, intent(out) :: z(3)
      character(len=2), intent(out) :: sym(3)
      real(dp), intent(out) :: xyz(3, 3)

      z = [8, 1, 1]
      sym = ["O ", "H ", "H "]
      xyz = to_bohr(reshape([0.0_dp, 0.0_dp, 0.1173_dp, &
                             0.0_dp, 0.7572_dp, -0.4692_dp, &
                             0.0_dp, -0.7572_dp, -0.4692_dp], [3, 3]))
   end subroutine water

   subroutine water_dimer(z, sym, xyz, pull)
      !! The hydrogen-bonded water dimer (S22), O...O 2.9 A, with the second
      !! water moved `pull` Angstrom along x
      integer, intent(out) :: z(6)
      character(len=2), intent(out) :: sym(6)
      real(dp), intent(out) :: xyz(3, 6)
      real(dp), intent(in) :: pull

      real(dp) :: ang(3, 6)

      z = [8, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "O ", "H ", "H "]
      ang = reshape([-1.551007_dp, -0.114520_dp, 0.000000_dp, &
                     -1.934259_dp, 0.762503_dp, 0.000000_dp, &
                     -0.599677_dp, 0.040712_dp, 0.000000_dp, &
                     1.350625_dp, 0.111469_dp, 0.000000_dp, &
                     1.680398_dp, -0.373741_dp, -0.758561_dp, &
                     1.680398_dp, -0.373741_dp, 0.758561_dp], [3, 6])
      ang(1, 4:6) = ang(1, 4:6) + pull
      xyz = to_bohr(ang)
   end subroutine water_dimer

   subroutine butane_water(z, sym, xyz)
      !! Anti butane, cut into two ethyls, and a water 4.5 A beyond C4, as
      !! `test_mqc_fmo_pieda` and `test_mqc_afo_fmo` have it
      integer, intent(out) :: z(17)
      character(len=2), intent(out) :: sym(17)
      real(dp), intent(out) :: xyz(3, 17)

      z = [6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1]
      sym = ["C ", "C ", "C ", "C ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", &
             "H ", "H ", "O ", "H ", "H "]
      xyz = to_bohr(reshape([0.000000_dp, 0.000000_dp, 0.000000_dp, &
                             1.268400_dp, 0.855600_dp, 0.000000_dp, &
                             2.536900_dp, 0.000000_dp, 0.000000_dp, &
                             3.805300_dp, 0.855600_dp, 0.000000_dp, &
                             0.272900_dp, -1.055300_dp, 0.000000_dp, &
                             -0.588900_dp, 0.222400_dp, -0.889800_dp, &
                             -0.588900_dp, 0.222400_dp, 0.889800_dp, &
                             1.268400_dp, 1.484700_dp, 0.890100_dp, &
                             1.268400_dp, 1.484700_dp, -0.890100_dp, &
                             2.536900_dp, -0.629100_dp, -0.890100_dp, &
                             2.536900_dp, -0.629100_dp, 0.890100_dp, &
                             3.532400_dp, 1.910800_dp, 0.000000_dp, &
                             4.394200_dp, 0.633100_dp, -0.889800_dp, &
                             4.394200_dp, 0.633100_dp, 0.889800_dp, &
                             7.535896_dp, 3.372076_dp, 0.000000_dp, &
                             7.575283_dp, 4.330466_dp, 0.000000_dp, &
                             8.439273_dp, 3.049628_dp, 0.000000_dp], [3, 17]))
   end subroutine butane_water

   subroutine water_trimer(z, sym, xyz)
      !! Three waters stacked 2.9 A apart, as `test_mqc_fmo_dft` has them
      integer, intent(out) :: z(9)
      character(len=2), intent(out) :: sym(9)
      real(dp), intent(out) :: xyz(3, 9)

      z = [8, 1, 1, 8, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "O ", "H ", "H ", "O ", "H ", "H "]
      xyz = to_bohr(reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                             0.0_dp, -0.7572_dp, 0.5865_dp, &
                             0.0_dp, 0.7572_dp, 0.5865_dp, &
                             0.0_dp, 0.0_dp, 2.9_dp, &
                             0.0_dp, -0.7572_dp, 3.4865_dp, &
                             0.0_dp, 0.7572_dp, 3.4865_dp, &
                             0.0_dp, 0.0_dp, 5.8_dp, &
                             0.0_dp, -0.7572_dp, 6.3865_dp, &
                             0.0_dp, 0.7572_dp, 6.3865_dp], [3, 9]))
   end subroutine water_trimer

end module test_mqc_fmo_pieda_methods

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fmo_pieda_methods, only: collect_mqc_fmo_pieda_methods
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fmo_pieda_methods", collect_mqc_fmo_pieda_methods)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
