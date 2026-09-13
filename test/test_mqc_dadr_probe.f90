program dadr_probe
   !! TEMPORARY: d(alpha)/dR for one nuclear coordinate, by differencing our own
   !! polarizability. The oracle the analytic assembly is built against.
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_polarizability, only: static_polarizability
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   integer, parameter :: Z(3) = [8, 1, 1]
   character(len=2), parameter :: SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: W(3, 3) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                             0.0_dp, 0.0_dp, 1.814137_dp, &
                                             0.0_dp, 1.756_dp, -0.4543_dp], [3, 3])
   real(dp), parameter :: WEIGHT(6) = [-1.0_dp, 9.0_dp, -45.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]
   integer, parameter :: OFFSET(6) = [-3, -2, -1, 1, 2, 3]
   real(dp), parameter :: STEP = 2.0e-3_dp

   real(dp) :: origin(3), coords(3, 3), alpha(3, 3), acc(3, 3)
   integer :: point, atom, cart, a, u
   type(error_t) :: err

   ! The origin is pinned across displacements: alpha is origin-independent, but
   ! pinning removes even that much variation from the difference.
   origin = (W(:, 1) + W(:, 2) + W(:, 3))/3.0_dp

   open (newunit=u, file="dadr_fd.txt", status="replace")
   do atom = 1, 3
      do cart = 1, 3
         acc = 0.0_dp
         do point = 1, 6
            coords = W
            coords(cart, atom) = coords(cart, atom) + real(OFFSET(point), dp)*STEP
            call alpha_at(coords, origin, alpha, err)
            if (err%has_error()) error stop "alpha: "//err%get_message()
            acc = acc + WEIGHT(point)*alpha
         end do
         acc = acc/(60.0_dp*STEP)
         do a = 1, 3
            write (u, "(3es24.15)") acc(a, :)
         end do
      end do
   end do
   close (u)
   write (*, "(a)") "wrote dadr_fd.txt (9 perturbations x 3x3)"

contains

   subroutine alpha_at(coords, origin, alpha, err)
      real(dp), intent(in) :: coords(3, 3), origin(3)
      real(dp), intent(out) :: alpha(3, 3)
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      call build_libcint_molecule(Z, SYM, coords, "6-31g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err)
      if (.not. err%has_error()) &
         call static_polarizability(mol, scf%orbitals, scf%orbital_energies, 5, &
                                    alpha, err, origin=origin)
      call mol%destroy()
   end subroutine alpha_at
end program dadr_probe
