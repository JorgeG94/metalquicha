program dadr_probe2
   !! TEMPORARY: the analytic d(alpha)/dR so far, against the finite-difference
   !! oracle, so the size of the missing term is a number rather than a guess.
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_hessian, only: response_hessian, assemble_polarizability_derivative, &
                                  solve_mo1_batch, hcore_deriv_atom, overlap_deriv_atom, &
                                  make_h1_atom
   use mqc_libcint_hess_ints, only: eri_ip1_block
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   integer, parameter :: Z(3) = [8, 1, 1]
   character(len=2), parameter :: SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: W(3, 3) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                             0.0_dp, 0.0_dp, 1.814137_dp, &
                                             0.0_dp, 1.756_dp, -0.4543_dp], [3, 3])
   type(libcint_molecule_t) :: mol
   type(rhf_result_t) :: scf
   type(error_t) :: err
   real(dp), allocatable :: hess(:, :, :, :), dadr(:, :, :, :)
   real(dp), allocatable :: mo1(:, :, :, :), s1(:, :, :, :), h1(:, :, :, :)
   real(dp), allocatable :: hcore_a(:, :, :), s1a(:, :, :), h1a(:, :, :)
   real(dp), allocatable :: eri_ip1(:, :, :, :, :)
   integer :: ia, u, a, b, c, nao

   call build_libcint_molecule(Z, SYM, W, "6-31g", mol, err)
   call run_libcint_rhf(mol, 10, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err)
   if (err%has_error()) error stop "scf: "//err%get_message()
   nao = mol%nao

   ! The nuclear response, the same way `response_hessian` builds it.
   call eri_ip1_block(mol, eri_ip1, err)
   allocate (h1(nao, nao, 3, mol%natm), s1(nao, nao, 3, mol%natm))
   do ia = 1, mol%natm
      call make_h1_atom(mol, scf%density, eri_ip1, ia, h1a, err)
      call hcore_deriv_atom(mol, ia, hcore_a, err)
      call overlap_deriv_atom(mol, ia, s1a, err)
      if (err%has_error()) error stop "skeletons: "//err%get_message()
      h1(:, :, :, ia) = h1a + hcore_a
      s1(:, :, :, ia) = s1a
      deallocate (hcore_a, s1a, h1a)
   end do
   call solve_mo1_batch(mol, scf%orbitals, scf%orbital_energies, 5, h1, s1, mo1, err)
   if (err%has_error()) error stop "mo1: "//err%get_message()

   call assemble_polarizability_derivative(mol, scf%orbitals, scf%orbital_energies, 5, &
                                           mo1, s1, dadr, err)
   if (err%has_error()) error stop "dadr: "//err%get_message()

   open (newunit=u, file="dadr_analytic.txt", status="replace")
   do ia = 1, mol%natm
      do c = 1, 3
         do a = 1, 3
            write (u, "(3es24.15)") (dadr(a, b, c, ia), b=1, 3)
         end do
      end do
   end do
   close (u)
   write (*, "(a)") "wrote dadr_analytic.txt"
   call mol%destroy()
end program dadr_probe2
