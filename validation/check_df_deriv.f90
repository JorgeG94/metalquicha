!! The density-fitting derivative integrals, before anything is built on them
program check_df_deriv
   !! Prints a few elements of d(mu nu|P)/dR and d(P|Q)/dR for comparison
   !! against PySCF's int3c2e_ip1, int3c2e_ip2 and int2c2e_ip1.
   !!
   !! Checked at this level on purpose. An index or ordering mistake in a
   !! three-centre derivative does not announce itself in the assembled
   !! gradient -- it produces a number of the right magnitude that disagrees
   !! with finite differences for reasons that could be anywhere in the
   !! contraction. Comparing the integrals themselves puts the error where it
   !! happened.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_gradient, only: three_centre_deriv, two_centre_deriv
   implicit none

   integer, parameter :: N_DIM = 3
   type(czt_molecule_t) :: orb, aux
   type(error_t) :: error
   real(dp), allocatable :: ip1(:, :, :, :), ip2(:, :, :, :), ip2c(:, :, :)
   integer :: i, j, k, c

   call build_czt_molecule([1, 1], ["H", "H"], &
                           reshape([0.0_dp, 0.0_dp, -0.7_dp, &
                                    0.0_dp, 0.0_dp, 0.7_dp], [N_DIM, 2]), &
                           "sto-3g", orb, error)
   if (error%has_error()) then
      write (*, "(a,a)") "orbital basis failed: ", error%get_message()
      error stop 1
   end if
   call build_czt_molecule([1, 1], ["H", "H"], &
                           reshape([0.0_dp, 0.0_dp, -0.7_dp, &
                                    0.0_dp, 0.0_dp, 0.7_dp], [N_DIM, 2]), &
                           "cc-pvdz-rifit", aux, error)
   if (error%has_error()) then
      write (*, "(a,a)") "auxiliary basis failed: ", error%get_message()
      error stop 1
   end if

   write (*, "(a,i0,a,i0)") "nao = ", orb%nao, "   naux = ", aux%nao

   call three_centre_deriv(orb, aux, 1, ip1)
   call three_centre_deriv(orb, aux, 2, ip2)
   call two_centre_deriv(aux, ip2c)

   write (*, "(a)") "ip1 (mu,nu,P,comp) for mu,nu = 1,2 and P = 1..3"
   do k = 1, min(3, aux%nao)
      write (*, "(a,i0,a,3f18.12)") "   P=", k, "  ", (ip1(1, 2, k, c), c=1, 3)
   end do
   write (*, "(a)") "ip2 (mu,nu,P,comp) for mu,nu = 1,2 and P = 1..3"
   do k = 1, min(3, aux%nao)
      write (*, "(a,i0,a,3f18.12)") "   P=", k, "  ", (ip2(1, 2, k, c), c=1, 3)
   end do
   ! Q on the *other* atom. Both indices on one centre makes (P|Q)
   ! independent of position, so its derivative is zero and the comparison
   ! would pass without testing anything.
   write (*, "(a,i0,a)") "2c2e_ip1 (P,Q,comp) for P=1, Q = ", aux%nao/2 + 1, "..+2"
   do k = aux%nao/2 + 1, min(aux%nao/2 + 3, aux%nao)
      write (*, "(a,i0,a,3f18.12)") "   Q=", k, "  ", (ip2c(1, k, c), c=1, 3)
   end do

end program check_df_deriv
