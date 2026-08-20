!! Density matrices in the active space
module mqc_rdm
   !! The one- and two-particle density matrices of a CI wave function,
   !!
   !!     D_pq      = <Psi| E_pq |Psi>
   !!     d_pqrs    = <Psi| E_pq E_rs - delta_qr E_ps |Psi>
   !!
   !! with `E_pq` summed over both spins, so these are spin-traced. They are
   !! everything the rest of MCSCF needs from the CI: the orbital gradient, the
   !! generalised Fock matrix and the energy are all contractions of these two
   !! against integrals, and none of them ever sees a CI vector again.
   !!
   !! Both fall out of the same operation the sigma build uses. If
   !! `t_pq = E_pq |Psi>` -- which `apply_excitations` returns for every pair at
   !! once -- then `D_pq` is the overlap of `t_pq` with the wave function, and
   !! since `E_pq` is the adjoint of `E_qp`,
   !!
   !!     <Psi| E_pq E_rs |Psi> = <t_qp | t_rs>
   !!
   !! so the whole two-particle matrix is one matrix multiply of that
   !! intermediate against itself. The `delta_qr E_ps` correction is then
   !! subtracted, which is what makes `d` the quantity that contracts with
   !! `(pq|rs)` to give the two-electron energy.
   !!
   !! The convention is PySCF's `make_rdm12`, which is also Helgaker's and
   !! GAMESS's. It is worth being explicit about because a factor of two and an
   !! index transposition are the two ways a two-particle density matrix is
   !! usually wrong, and both give a plausible energy.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_determinants, only: link_table_t
   use mqc_ci, only: excitations_block, beta_strings_per_block
   implicit none
   private

   public :: active_space_rdms
   public :: rdm_energy

contains

   subroutine active_space_rdms(ci, alpha, beta, dm1, dm2, error)
      !! Spin-traced one- and two-particle density matrices
      real(dp), intent(in) :: ci(:, :)       !! (n_alpha_strings, n_beta_strings)
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), allocatable, intent(out) :: dm1(:, :)
      real(dp), allocatable, intent(out) :: dm2(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gathered(:, :), paired(:, :)
      real(dp), allocatable :: flat(:)
      integer :: norb, na, nb, npair, ndet, p, q, r, s, pq, qp, rs
      integer :: per_block, first, last, width

      if (error%has_error()) return
      norb = alpha%n_orbitals
      na = alpha%n_strings
      nb = beta%n_strings
      npair = norb*norb
      ndet = na*nb

      if (beta%n_orbitals /= norb) then
         call error%set(ERROR_VALIDATION, "the alpha and beta excitation tables "// &
                        "describe different active spaces: "//to_char(norb)//" and "// &
                        to_char(beta%n_orbitals)//" orbitals.")
         return
      end if
      if (size(ci, 1) /= na .or. size(ci, 2) /= nb) then
         call error%set(ERROR_VALIDATION, "the vector is "//to_char(size(ci, 1))// &
                        " by "//to_char(size(ci, 2))//" but the tables have "// &
                        to_char(na)//" alpha and "//to_char(nb)//" beta strings.")
         return
      end if

      ! Blocked over beta strings, for the reason `beta_strings_per_block`
      ! gives: the intermediate is `norb^2` by the determinant count, which for
      ! a large active space is tens of gigabytes. Both contractions below sum
      ! over determinants, so a block contributes a partial sum and the totals
      ! accumulate -- exactly, not approximately.
      per_block = beta_strings_per_block(npair, na, nb)
      allocate (gathered(npair, na*per_block))
      allocate (dm1(norb, norb), paired(npair, npair), flat(na*per_block))
      dm1 = 0.0_dp
      paired = 0.0_dp

      do first = 1, nb, per_block
         last = min(first + per_block - 1, nb)
         width = na*(last - first + 1)

         call excitations_block(ci, alpha, beta, first, last, gathered(:, 1:width))
         flat(1:width) = reshape(ci(:, first:last), [width])

         ! D_pq = <Psi| E_pq |Psi>
         do q = 1, norb
            do p = 1, norb
               pq = p + (q - 1)*norb
               dm1(p, q) = dm1(p, q) &
                           + dot_product(gathered(pq, 1:width), flat(1:width))
            end do
         end do

         ! <Psi| E_pq E_rs |Psi> = <E_qp Psi | E_rs Psi>, every pair against
         ! every other. `beta=1` so the blocks add rather than overwrite.
         call pic_gemm(gathered(:, 1:width), gathered(:, 1:width), paired, &
                       transb="T", alpha=1.0_dp, beta=1.0_dp)
      end do

      allocate (dm2(norb, norb, norb, norb))
      do s = 1, norb
         do r = 1, norb
            rs = r + (s - 1)*norb
            do q = 1, norb
               do p = 1, norb
                  qp = q + (p - 1)*norb
                  dm2(p, q, r, s) = paired(qp, rs)
               end do
            end do
         end do
      end do
      ! The delta_qr correction, which turns <E E> into the quantity that
      ! contracts with (pq|rs).
      do s = 1, norb
         do r = 1, norb
            do p = 1, norb
               dm2(p, r, r, s) = dm2(p, r, r, s) - dm1(p, s)
            end do
         end do
      end do

      deallocate (gathered, paired, flat)
   end subroutine active_space_rdms

   pure function rdm_energy(h1e, eri, dm1, dm2) result(energy)
      !! The active-space energy rebuilt from the density matrices
      !!
      !!     E = sum_pq h_pq D_pq + (1/2) sum_pqrs (pq|rs) d_pqrs
      !!
      !! Nothing needs this to run an MCSCF -- the CI already produced the
      !! energy as an eigenvalue. It exists because the two numbers agreeing is
      !! the strongest single check available on the density matrices: they are
      !! reached by completely different arithmetic, and a transposed index or a
      !! factor of two in `d` that leaves every trace identity intact will not
      !! survive it.
      real(dp), intent(in) :: h1e(:, :)
      real(dp), intent(in) :: eri(:, :, :, :)
      real(dp), intent(in) :: dm1(:, :)
      real(dp), intent(in) :: dm2(:, :, :, :)
      real(dp) :: energy

      integer :: norb, p, q, r, s

      norb = size(dm1, 1)
      energy = 0.0_dp
      do q = 1, norb
         do p = 1, norb
            energy = energy + h1e(p, q)*dm1(p, q)
         end do
      end do
      do s = 1, norb
         do r = 1, norb
            do q = 1, norb
               do p = 1, norb
                  energy = energy + 0.5_dp*eri(p, q, r, s)*dm2(p, q, r, s)
               end do
            end do
         end do
      end do
   end function rdm_energy

end module mqc_rdm
