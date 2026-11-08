!! The analytic RI-MP2 correlation Hessian
module test_mqc_ri_mp2_hessian
   !! Two oracles, because each catches what the other cannot.
   !!
   !! **Differences of `czt_ri_mp2_gradient`.** A fixed-step stencil of the
   !! analytic gradient is a derivative of exactly the energy the Hessian claims
   !! to differentiate, and needs nothing external. Its tolerance has to absorb
   !! the stencil and two SCF convergences, which is loose enough to hide a sign
   !! error in a small term.
   !!
   !! **Numbers from the numpy prototype**, `tools/cpu_validation/
   !! ri_mp2_hessian.py --dump`, read from `tools/mp2_hessian_oracle/
   !! ri_mp2_hessian_reference.json`. Those agree with the Fortran to the
   !! convergence of the SCF, so the tolerance is tight, and the file also
   !! carries the four additive pieces so that a disagreement says where it is.
   !! Three cases: all-electron 6-31G, frozen-core 6-31G with a triple-zeta
   !! fitting set, and cc-pVDZ, whose d shells put the integral bookkeeping to
   !! work that s and p shells do not.
   !!
   !! The reference is the RI-MP2 energy of an **exact** RHF reference, so the
   !! gradient is called with `fit_response = .false.`.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use json_module, only: json_file
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, FIT_METRIC_THRESHOLD
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_ri_mp2_gradient, only: czt_ri_mp2_gradient
   use mqc_czt_ri_mp2_hessian, only: ri_mp2_correlation_hessian, ri_hessian_parts_t, &
                                     inverse_metric
   use mqc_czt_hessian, only: rhf_hessian
   implicit none
   private

   public :: collect_mqc_ri_mp2_hessian_tests

   character(len=*), parameter :: REFERENCE = &
                                  "../tools/mp2_hessian_oracle/ri_mp2_hessian_reference.json"
   real(dp), parameter :: TOL_PIECES = 1.0e-9_dp
      !! The four pieces and the total agree with the prototype to 1e-12 to
      !! 3e-12 in every case, threaded; this leaves three orders of magnitude for
      !! a platform's libm and BLAS while still failing on a term at the 1e-8
      !! level -- a missing cross term, say, is 1e-3.
   real(dp), parameter :: STEP = 2.0e-3_dp
   real(dp), parameter :: TOL_FD = 3.0e-7_dp
      !! Of the symmetrised stencil. The gradient's own Z-vector is iterative and
      !! its error, amplified by the step, makes the raw stencil asymmetric by
      !! 2.3e-7 (the same figure as its deviation from the analytic Hessian);
      !! symmetrised it sits 1.2e-7 away. The prototype's numbers are the tight
      !! oracle, at 1e-9; this one is for a sign or a factor, which are
      !! orders of magnitude larger.

contains

   subroutine collect_mqc_ri_mp2_hessian_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("reference_all_electron", reference_all_electron), &
                  new_unittest("reference_frozen_core", reference_frozen_core), &
                  new_unittest("reference_d_shells", reference_d_shells), &
                  new_unittest("against_gradient_differences", against_gradient_differences), &
                  new_unittest("frozen_core_against_gradient_differences", &
                               frozen_against_gradient_differences), &
                  new_unittest("batching_changes_nothing", batching_changes_nothing), &
                  new_unittest("translations_leave_it_alone", translations), &
                  new_unittest("a_metric_that_drops_a_mode_is_counted", dropped_modes) &
                  ]
   end subroutine collect_mqc_ri_mp2_hessian_tests

   function atomic_number(symbol) result(z)
      character(len=*), intent(in) :: symbol
      integer :: z

      select case (trim(adjustl(symbol)))
      case ("H")
         z = 1
      case ("C")
         z = 6
      case ("N")
         z = 7
      case ("O")
         z = 8
      case default
         z = 0
      end select
   end function atomic_number

   subroutine read_matrix(json, base, n, matrix)
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: base
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: matrix(:, :)

      real(dp), allocatable :: row(:)
      character(len=256) :: path
      integer :: i
      logical :: found

      allocate (matrix(n, n))
      do i = 1, n
         write (path, "(a,a,i0,a)") trim(base), "(", i, ")"
         call json%get(trim(path), row, found)
         if (.not. found) then
            matrix = huge(1.0_dp)
            return
         end if
         matrix(i, :) = row
      end do
   end subroutine read_matrix

   subroutine to_matrix(h4, matrix)
      !! `(3, 3, natm, natm)` to `(3N, 3N)`, row `3*(A-1)+a`, column `3*(B-1)+b`
      real(dp), intent(in) :: h4(:, :, :, :)
      real(dp), allocatable, intent(out) :: matrix(:, :)

      integer :: natm, ia, ib, a, b

      natm = size(h4, 3)
      allocate (matrix(3*natm, 3*natm))
      do ib = 1, natm
         do ia = 1, natm
            do b = 1, 3
               do a = 1, 3
                  matrix(3*(ia - 1) + a, 3*(ib - 1) + b) = h4(a, b, ia, ib)
               end do
            end do
         end do
      end do
   end subroutine to_matrix

   subroutine converge(geo, symbols, basis, aux_basis, mol, aux, scf, err)
      real(dp), intent(in) :: geo(:, :)
      character(len=2), intent(in) :: symbols(:)
      character(len=*), intent(in) :: basis, aux_basis
      type(czt_molecule_t), intent(out) :: mol, aux
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      integer, allocatable :: z(:)
      integer :: i

      allocate (z(size(symbols)))
      do i = 1, size(symbols)
         z(i) = atomic_number(symbols(i))
      end do
      call build_czt_molecule(z, symbols, geo, basis, mol, err)
      if (err%has_error()) return
      call build_czt_molecule(z, symbols, geo, aux_basis, aux, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, sum(z), 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err)
   end subroutine converge

   subroutine one_reference(error, icase)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: icase

      type(json_file) :: json
      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(ri_hessian_parts_t) :: parts
      type(error_t) :: err
      real(dp), allocatable :: hess4(:, :, :, :), hess(:, :)
      real(dp), allocatable :: ref_total(:, :), ref_p1r(:, :), ref_p1f(:, :), ref_orb(:, :)
      real(dp), allocatable :: ref_resp(:, :), geo(:, :), rows(:)
      character(len=2), allocatable :: symbols(:)
      character(len=:), allocatable :: sym_buf
      character(len=:), allocatable :: basis, aux_basis
      character(len=256) :: base, path
      integer :: n3, natm, nfrozen, i, n_occ
      logical :: found

      call json%initialize()
      call json%load(filename=REFERENCE)
      call check(error,.not. json%failed(), "the prototype's reference numbers must load")
      if (allocated(error)) return

      write (base, "(a,i0,a)") "cases(", icase, ")"
      call json%get(trim(base)//".basis", basis, found)
      call json%get(trim(base)//".aux_basis", aux_basis, found)
      call json%get(trim(base)//".n_frozen", nfrozen, found)
      call json%get(trim(base)//".n3", n3, found)
      call check(error, found, "the reference case must have a size")
      if (allocated(error)) return
      natm = n3/3
      allocate (symbols(natm), geo(3, natm))
      do i = 1, natm
         write (path, "(a,a,i0,a)") trim(base), ".symbols(", i, ")"
         call json%get(trim(path), sym_buf, found)
         symbols(i) = adjustl(sym_buf)
         write (path, "(a,a,i0,a)") trim(base), ".geometry_bohr(", i, ")"
         call json%get(trim(path), rows, found)
         geo(:, i) = rows
      end do
      call read_matrix(json, trim(base)//".hessian", n3, ref_total)
      call read_matrix(json, trim(base)//".pass1_reference", n3, ref_p1r)
      call read_matrix(json, trim(base)//".pass1_fitted", n3, ref_p1f)
      call read_matrix(json, trim(base)//".orbital", n3, ref_orb)
      call read_matrix(json, trim(base)//".response", n3, ref_resp)
      call json%destroy()

      call converge(geo, symbols, basis, aux_basis, mol, aux, scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge: "//err%get_full_trace())
      if (allocated(error)) return
      n_occ = sum([(atomic_number(symbols(i)), i=1, natm)])/2
      call ri_mp2_correlation_hessian(mol, aux, scf%orbitals, scf%orbital_energies, n_occ, &
                                      hess4, err, n_frozen=nfrozen, tol=1.0e-13_dp, parts=parts)
      call check(error,.not. err%has_error(), "the Hessian must run: "//err%get_full_trace())
      if (allocated(error)) return
      call to_matrix(hess4, hess)

      call check(error, maxval(abs(parts%pass1_reference - ref_p1r)) < TOL_PIECES, &
                 "pass 1, the reference operator's share, disagrees with the prototype")
      if (allocated(error)) return
      call check(error, maxval(abs(parts%orbital - ref_orb)) < TOL_PIECES, &
                 "the orbital response disagrees with the prototype")
      if (allocated(error)) return
      ! The fitted share of pass 1 splits differently: the cross terms ride with
      ! the pair sum here, in one sweep. Their sum is the same.
      call check(error, maxval(abs((parts%pass1_fitted + parts%response) &
                                   - (ref_p1f + ref_resp))) < TOL_PIECES, &
                 "the fitted pass 1 plus the density response disagrees with the prototype")
      if (allocated(error)) return
      call check(error, maxval(abs(hess - ref_total)) < TOL_PIECES, &
                 "the correlation Hessian disagrees with the prototype")
      if (allocated(error)) return
      call check(error, maxval(abs(hess - transpose(hess))) < 1.0e-8_dp, &
                 "the correlation Hessian is not symmetric")
      call mol%destroy()
      call aux%destroy()
   end subroutine one_reference

   subroutine reference_all_electron(error)
      type(error_type), allocatable, intent(out) :: error

      call one_reference(error, 1)
   end subroutine reference_all_electron

   subroutine reference_frozen_core(error)
      type(error_type), allocatable, intent(out) :: error

      call one_reference(error, 2)
   end subroutine reference_frozen_core

   subroutine reference_d_shells(error)
      type(error_type), allocatable, intent(out) :: error

      call one_reference(error, 3)
   end subroutine reference_d_shells

   subroutine total_gradient(geo, symbols, basis, aux_basis, nfrozen, grad, err)
      !! The RI-MP2 gradient (reference included) at one geometry, flattened
      real(dp), intent(in) :: geo(:, :)
      character(len=2), intent(in) :: symbols(:)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nfrozen
      real(dp), intent(out) :: grad(:)
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      real(dp), allocatable :: g(:, :)
      integer :: i, n_occ

      call converge(geo, symbols, basis, aux_basis, mol, aux, scf, err)
      if (err%has_error()) return
      n_occ = sum([(atomic_number(symbols(i)), i=1, size(symbols))])/2
      call czt_ri_mp2_gradient(mol, aux, scf%orbitals, scf%orbital_energies, n_occ, g, err, &
                               n_frozen=nfrozen, fit_response=.false.)
      if (err%has_error()) return
      grad = reshape(g, [size(grad)])
      call mol%destroy()
      call aux%destroy()
   end subroutine total_gradient

   subroutine differences_case(error, basis, aux_basis, nfrozen)
      !! RHF Hessian plus this block, against a five-point stencil of the gradient
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nfrozen

      real(dp), parameter :: WATER(3, 3) = reshape( &
                             [0.05_dp, 0.02_dp, -0.03_dp, &
                              0.0_dp, -1.4308_dp, 1.1078_dp, &
                              0.04_dp, 1.4308_dp, 1.1078_dp], [3, 3])
      character(len=2), parameter :: SYM(3) = ["O ", "H ", "H "]
      real(dp), parameter :: SHIFT(4) = [2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp]
      real(dp), parameter :: WEIGHT(4) = [-1.0_dp, 8.0_dp, -8.0_dp, 1.0_dp]
      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: hess4(:, :, :, :), corr(:, :), rhf4(:, :, :, :), rhf(:, :)
      real(dp) :: geo(3, 3), g(9, 4), fd(9, 9), dev
      integer :: y, k, atom, comp

      call converge(WATER, SYM, basis, aux_basis, mol, aux, scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call ri_mp2_correlation_hessian(mol, aux, scf%orbitals, scf%orbital_energies, 5, &
                                      hess4, err, n_frozen=nfrozen, tol=1.0e-13_dp)
      call check(error,.not. err%has_error(), "the Hessian must run: "//err%get_full_trace())
      if (allocated(error)) return
      call rhf_hessian(mol, [8, 1, 1], scf%density, scf%orbitals, scf%orbital_energies, 5, &
                       rhf4, err, tol=1.0e-13_dp)
      call check(error,.not. err%has_error(), "the RHF Hessian must run")
      if (allocated(error)) return
      call to_matrix(hess4, corr)
      call to_matrix(rhf4, rhf)

      do y = 1, 9
         atom = (y - 1)/3 + 1
         comp = y - 3*(atom - 1)
         do k = 1, 4
            geo = WATER
            geo(comp, atom) = geo(comp, atom) + SHIFT(k)*STEP
            call total_gradient(geo, SYM, basis, aux_basis, nfrozen, g(:, k), err)
            call check(error,.not. err%has_error(), "the displaced gradient failed")
            if (allocated(error)) return
         end do
         fd(:, y) = (WEIGHT(1)*g(:, 1) + WEIGHT(2)*g(:, 2) + WEIGHT(3)*g(:, 3) &
                     + WEIGHT(4)*g(:, 4))/(12.0_dp*STEP)
      end do
      fd = 0.5_dp*(fd + transpose(fd))
      dev = maxval(abs(corr + rhf - fd))
      call check(error, dev < TOL_FD, "RHF Hessian plus the correlation block is not the "// &
                 "derivative of the RI-MP2 gradient")
      if (allocated(error)) write (*, "(a,es10.3)") "largest deviation ", dev
      call mol%destroy()
      call aux%destroy()
   end subroutine differences_case

   subroutine against_gradient_differences(error)
      type(error_type), allocatable, intent(out) :: error

      call differences_case(error, "6-31g", "cc-pvdz-rifit", 0)
   end subroutine against_gradient_differences

   subroutine frozen_against_gradient_differences(error)
      type(error_type), allocatable, intent(out) :: error

      call differences_case(error, "6-31g", "cc-pvdz-rifit", 1)
   end subroutine frozen_against_gradient_differences

   subroutine batching_changes_nothing(error)
      !! Atoms per three-index batch and perturbations per Fock pass are
      !! performance knobs; the answer must not move with them
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: WATER(3, 3) = reshape( &
                             [0.05_dp, 0.02_dp, -0.03_dp, &
                              0.0_dp, -1.4308_dp, 1.1078_dp, &
                              0.04_dp, 1.4308_dp, 1.1078_dp], [3, 3])
      character(len=2), parameter :: SYM(3) = ["O ", "H ", "H "]
      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: whole(:, :, :, :), split(:, :, :, :)

      call converge(WATER, SYM, "6-31g", "cc-pvdz-rifit", mol, aux, scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge")
      if (allocated(error)) return
      call ri_mp2_correlation_hessian(mol, aux, scf%orbitals, scf%orbital_energies, 5, &
                                      whole, err, n_frozen=1, tol=1.0e-13_dp)
      call ri_mp2_correlation_hessian(mol, aux, scf%orbitals, scf%orbital_energies, 5, &
                                      split, err, n_frozen=1, tol=1.0e-13_dp, batch=2, &
                                      atom_batch=1)
      call check(error,.not. err%has_error(), "both runs must succeed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(whole - split)) < 1.0e-10_dp, &
                 "the batch sizes change the Hessian")
      call mol%destroy()
      call aux%destroy()
   end subroutine batching_changes_nothing

   subroutine translations(error)
      !! Moving every atom together moves nothing: each row sums to zero over
      !! the second atom. This is the correlation block alone, which has to be
      !! invariant by itself because the reference block is.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: WATER(3, 3) = reshape( &
                             [0.05_dp, 0.02_dp, -0.03_dp, &
                              0.0_dp, -1.4308_dp, 1.1078_dp, &
                              0.04_dp, 1.4308_dp, 1.1078_dp], [3, 3])
      character(len=2), parameter :: SYM(3) = ["O ", "H ", "H "]
      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: h4(:, :, :, :)
      real(dp) :: worst
      integer :: ia, a, b

      call converge(WATER, SYM, "6-31g", "cc-pvdz-rifit", mol, aux, scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge")
      if (allocated(error)) return
      call ri_mp2_correlation_hessian(mol, aux, scf%orbitals, scf%orbital_energies, 5, &
                                      h4, err, n_frozen=1, tol=1.0e-13_dp)
      call check(error,.not. err%has_error(), "the Hessian must run")
      if (allocated(error)) return
      worst = 0.0_dp
      do ia = 1, 3
         do b = 1, 3
            do a = 1, 3
               worst = max(worst, abs(sum(h4(a, b, ia, :))))
            end do
         end do
      end do
      call check(error, worst < 1.0e-8_dp, "translating all atoms changes the correlation energy")
      call mol%destroy()
      call aux%destroy()
   end subroutine translations

   subroutine dropped_modes(error)
      !! The pseudo-inverse of the metric, and the count the Hessian refuses on
      !!
      !! `-J^-1 dJ J^-1` is the derivative of `J^-1` only when no mode is dropped,
      !! so the count has to be right. A positive-definite matrix drops none and
      !! its inverse is an inverse; a matrix with an exact null direction drops
      !! one and its pseudo-inverse still inverts the rest.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: metric(3, 3), jinv(3, 3)
      real(dp), allocatable :: inv(:, :)
      type(error_t) :: err
      integer :: dropped

      metric = reshape([4.0_dp, 1.0_dp, 0.0_dp, &
                        1.0_dp, 3.0_dp, 0.5_dp, &
                        0.0_dp, 0.5_dp, 2.0_dp], [3, 3])
      call inverse_metric(metric, inv, dropped, err)
      call check(error, dropped, 0, "a positive-definite metric drops no mode")
      if (allocated(error)) return
      jinv = matmul(inv, metric)
      call check(error, maxval(abs(jinv - reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                                                   0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3]))) &
                 < 1.0e-13_dp, "J^-1 J must be the identity when nothing is dropped")
      if (allocated(error)) return

      ! rank two: the third column is the sum of the first two
      metric = reshape([2.0_dp, 1.0_dp, 3.0_dp, &
                        1.0_dp, 2.0_dp, 3.0_dp, &
                        3.0_dp, 3.0_dp, 6.0_dp], [3, 3])
      call inverse_metric(metric, inv, dropped, err)
      call check(error, dropped, 1, "an exact null direction is one dropped mode")
      if (allocated(error)) return
      call check(error, FIT_METRIC_THRESHOLD > 0.0_dp, "the threshold is a positive number")
   end subroutine dropped_modes

end module test_mqc_ri_mp2_hessian

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_ri_mp2_hessian, only: collect_mqc_ri_mp2_hessian_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_ri_mp2_hessian", collect_mqc_ri_mp2_hessian_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
