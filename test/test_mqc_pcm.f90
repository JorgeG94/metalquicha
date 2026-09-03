!! The CPU continuum: cavity, charges, and the solvated SCF
module test_mqc_pcm
   !! What can be asserted without a reference implementation is the physics:
   !! an isolated sphere's tessellation must integrate to its own surface area,
   !! points buried between overlapping spheres must vanish, and a central
   !! charge in a spherical cavity is Born's problem with a closed-form energy.
   !! C-PCM approaches it as -f q^2/2R with f = (eps-1)/eps; IEF-PCM's D-matrix
   !! terms make it exact for a sphere, so it must land on the *same* Born
   !! energy despite carrying f = (eps-1)/(eps+1) inside.
   !!
   !! The one pinned number is a solvated water SCF, validated against PySCF's
   !! `pyscf.solvent.pcm` (C-PCM, identical cavity radii and Lebedev order):
   !! PySCF gives -74.9684556809 on this geometry and this implementation
   !! -74.9684557057, 2.5e-8 apart. The pin catches drift, not correctness;
   !! the Born tests carry the correctness.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_constants, only: PI
   use mqc_pcm_radii, only: cavity_radius
   use mqc_method_config, only: pcm_config_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_pcm, only: pcm_context_t, build_pcm_surface
   implicit none
   private

   public :: collect_mqc_pcm_tests

   real(dp), parameter :: EPS_WATER = 78.3553_dp

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   !! Bohr, C2v
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])

contains

   subroutine collect_mqc_pcm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("one_sphere_area_is_exact", one_sphere_area_is_exact), &
                  new_unittest("overlapping_spheres_bury_points", overlapping_spheres_bury_points), &
                  new_unittest("born_ion_cpcm", born_ion_cpcm), &
                  new_unittest("born_ion_iefpcm_is_exact", born_ion_iefpcm_is_exact), &
                  new_unittest("solvated_water_has_not_drifted", solvated_water_has_not_drifted), &
                  new_unittest("refuses_an_unknown_model", refuses_an_unknown_model), &
                  new_unittest("refuses_a_gas_dielectric", refuses_a_gas_dielectric), &
                  new_unittest("refuses_an_untabulated_order", refuses_an_untabulated_order) &
                  ]
   end subroutine collect_mqc_pcm_tests

   subroutine born_energy(method, e_diel, q_total, err)
      !! A unit positive charge at the centre of one oxygen-sized cavity
      character(len=*), intent(in) :: method
      real(dp), intent(out) :: e_diel, q_total
      type(error_t), intent(inout) :: err

      type(pcm_context_t) :: ctx
      real(dp), allocatable :: v(:), q(:)
      real(dp) :: d
      integer :: k

      call build_pcm_surface([8], reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), &
                             302, 1.2_dp, ctx%surface, err)
      if (err%has_error()) return
      call ctx%init_model(method, EPS_WATER, err)
      if (err%has_error()) return
      call ctx%build_operators(err)
      if (err%has_error()) return

      ! The potential of the central charge at each smooth surface point --
      ! erf-attenuated exactly as the SCF's nuclear term is, though at these
      ! exponents the erf is 1 to machine precision.
      allocate (v(ctx%surface%n_points))
      do k = 1, ctx%surface%n_points
         d = norm2(ctx%surface%coords(:, k))
         v(k) = erf(ctx%surface%zeta(k)*d)/d
      end do

      call ctx%solve(v, q, err)
      if (err%has_error()) return
      e_diel = 0.5_dp*dot_product(q, v)
      q_total = sum(q)
   end subroutine born_energy

   subroutine one_sphere_area_is_exact(error)
      !! An isolated sphere keeps every point and integrates to 4 pi R^2
      !!
      !! The Lebedev weights sum to one and nothing is switched off, so the
      !! area sum is exact to roundoff -- this is the identity that catches a
      !! misplaced 4 pi, a radius applied twice, or a weight normalization
      !! change in the grid tables.
      type(error_type), allocatable, intent(out) :: error

      type(pcm_context_t) :: ctx
      type(error_t) :: err
      real(dp) :: radius

      call build_pcm_surface([8], reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), &
                             110, 1.2_dp, ctx%surface, err)
      call check(error,.not. err%has_error(), more="surface build failed")
      if (allocated(error)) return

      call check(error, ctx%surface%n_points, 110, more="an isolated sphere buried points")
      if (allocated(error)) return

      call cavity_radius(8, 1.2_dp, radius, err)
      call check(error, sum(ctx%surface%area), 4.0_dp*PI*radius**2, &
                 thr=1.0e-10_dp, more="the tessellated area is not the sphere's")
   end subroutine one_sphere_area_is_exact

   subroutine overlapping_spheres_bury_points(error)
      !! Two oxygens a bond apart lose the points between them
      type(error_type), allocatable, intent(out) :: error

      type(pcm_context_t) :: ctx
      type(error_t) :: err
      real(dp) :: radius

      call build_pcm_surface([8, 8], reshape( &
                             [0.0_dp, 0.0_dp, 0.0_dp, &
                              0.0_dp, 0.0_dp, 2.3_dp], [3, 2]), &
                             110, 1.2_dp, ctx%surface, err)
      call check(error,.not. err%has_error(), more="surface build failed")
      if (allocated(error)) return

      call check(error, ctx%surface%n_points < 220, &
                 more="no point was buried between overlapping spheres")
      if (allocated(error)) return
      call check(error, ctx%surface%n_discarded, 220 - ctx%surface%n_points, &
                 more="kept plus discarded is not the full count")
      if (allocated(error)) return

      call cavity_radius(8, 1.2_dp, radius, err)
      call check(error, sum(ctx%surface%area) < 2.0_dp*4.0_dp*PI*radius**2, &
                 more="the fused surface is not smaller than two full spheres")
   end subroutine overlapping_spheres_bury_points

   subroutine born_ion_cpcm(error)
      !! C-PCM reproduces its own Born expression, -f q^2/(2R), f = (eps-1)/eps
      !!
      !! For a conductor-like model on a sphere the discretization is the only
      !! error, and with 302 points it is parts in 1e6. The apparent charge
      !! must integrate to -f by Gauss's law.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: e_diel, q_total, radius, f, e_ref

      call born_energy("cpcm", e_diel, q_total, err)
      call check(error,.not. err%has_error(), more="the Born solve failed")
      if (allocated(error)) return

      call cavity_radius(8, 1.2_dp, radius, err)
      f = (EPS_WATER - 1.0_dp)/EPS_WATER
      e_ref = -0.5_dp*f/radius
      call check(error, e_diel, e_ref, thr=1.0e-5_dp*abs(e_ref), &
                 more="C-PCM misses the Born energy of its own model")
      if (allocated(error)) return
      call check(error, q_total, -f, thr=1.0e-4_dp, &
                 more="the apparent charge violates Gauss's law")
   end subroutine born_ion_cpcm

   subroutine born_ion_iefpcm_is_exact(error)
      !! IEF-PCM lands on the *true* Born energy, -q^2 (eps-1)/(2 R eps)
      !!
      !! Its internal scale factor is (eps-1)/(eps+1), so getting the physical
      !! (eps-1)/eps out is what exercises the D A S terms -- drop or misscale
      !! them and this is off by roughly 1/(2 eps), far outside the tolerance.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: e_diel, q_total, radius, e_ref

      call born_energy("iefpcm", e_diel, q_total, err)
      call check(error,.not. err%has_error(), more="the Born solve failed")
      if (allocated(error)) return

      call cavity_radius(8, 1.2_dp, radius, err)
      e_ref = -0.5_dp*(EPS_WATER - 1.0_dp)/EPS_WATER/radius
      call check(error, e_diel, e_ref, thr=1.0e-4_dp*abs(e_ref), &
                 more="IEF-PCM misses the exact Born energy for a sphere")
   end subroutine born_ion_iefpcm_is_exact

   subroutine solvated_water_has_not_drifted(error)
      !! The full path: cavity, charges, SCF -- pinned to the PySCF-validated run
      !!
      !! Water RHF/STO-3G, C-PCM, eps 78.3553, 302 points. The pin holds this
      !! implementation's own converged value; the tolerance is far looser than
      !! the 2.5e-8 gap to PySCF, so drifting past it means something moved.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(pcm_context_t) :: ctx
      type(pcm_config_t) :: config
      type(rhf_result_t) :: scf
      type(error_t) :: err

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), more="molecule build failed")
      if (allocated(error)) return

      config%enabled = .true.
      config%method = "cpcm"
      config%dielectric = EPS_WATER
      config%angular_points = 302
      call ctx%build(mol, WATER_Z, config, err)
      call check(error,.not. err%has_error(), more="context build failed")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call run_czt_rhf(mol, 10, 100, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err, &
                       pcm=ctx)
      call mol%destroy()
      call check(error,.not. err%has_error(), more="the solvated SCF failed")
      if (allocated(error)) return
      call check(error, scf%converged, more="the solvated SCF did not converge")
      if (allocated(error)) return

      call check(error, scf%energy, -74.96845570_dp, thr=1.0e-6_dp, &
                 more="the solvated total energy has drifted")
      if (allocated(error)) return
      call check(error, scf%pcm_energy, -0.00571472_dp, thr=1.0e-6_dp, &
                 more="the dielectric energy has drifted")
   end subroutine solvated_water_has_not_drifted

   subroutine refuses_an_unknown_model(error)
      !! "cosmo" is a real model this path does not solve; it must not be mapped
      type(error_type), allocatable, intent(out) :: error

      type(pcm_context_t) :: ctx
      type(error_t) :: err

      call ctx%init_model("cosmo", EPS_WATER, err)
      call check(error, err%has_error(), more="an unknown model was accepted")
   end subroutine refuses_an_unknown_model

   subroutine refuses_a_gas_dielectric(error)
      !! eps = 1 is the gas phase pretending to be a solvent
      type(error_type), allocatable, intent(out) :: error

      type(pcm_context_t) :: ctx
      type(error_t) :: err

      call ctx%init_model("cpcm", 1.0_dp, err)
      call check(error, err%has_error(), more="a dielectric of one was accepted")
   end subroutine refuses_a_gas_dielectric

   subroutine refuses_an_untabulated_order(error)
      !! Lebedev order 74 exists in mqc_lebedev but has no fitted SWIG exponent
      type(error_type), allocatable, intent(out) :: error

      type(pcm_context_t) :: ctx
      type(error_t) :: err

      call build_pcm_surface([8], reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), &
                             74, 1.2_dp, ctx%surface, err)
      call check(error, err%has_error(), more="an order with no fitted exponent "// &
                 "was accepted")
   end subroutine refuses_an_untabulated_order

end module test_mqc_pcm

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_pcm, only: collect_mqc_pcm_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_pcm", collect_mqc_pcm_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
