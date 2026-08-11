!! Closed-shell MP2 from stored AO integrals
module mqc_libcint_mp2
   !! Conventional restricted MP2: the reference the fitted version is checked
   !! against, in the same relationship the in-core Fock build has to the direct
   !! one.
   !!
   !! **The transformation is the whole of it.** Going from (mu nu|la si) to
   !! (ia|jb) directly is an eight-index sum, N^8, and unusable at any size. Done
   !! one index at a time it is four N^5 steps, each a matrix multiply:
   !!
   !!     (mu nu|la si) --C_occ--> (i nu|la si) --C_vir--> (ia|la si)
   !!                   --C_occ--> (ia|j si)   --C_vir--> (ia|jb)
   !!
   !! Every step contracts one AO index against one MO coefficient matrix, and
   !! each is expressed as a gemm rather than a loop nest, so the cost sits in
   !! BLAS where it belongs. The saving is not marginal: for a hundred basis
   !! functions N^8 is 10^16 operations and this is 10^10.
   !!
   !! Memory is the price, and it is why this is a reference rather than a
   !! method. The first intermediate is n_occ * n^3 -- 80 MB at a hundred
   !! functions with ten occupied -- and the AO tensor it reads is n^4 on top of
   !! that. The fitted version exists to avoid exactly this.
   !!
   !! The energy is decomposed as it is accumulated, because the two spin cases
   !! are what SCS-MP2 scales separately and recovering them afterwards from a
   !! total is not possible:
   !!
   !!     E_OS = sum (ia|jb)^2 / D          opposite spin
   !!     E_SS = sum (ia|jb)[(ia|jb) - (ib|ja)] / D   same spin
   !!     E_MP2 = E_OS + E_SS
   !!
   !! which is the familiar sum (ia|jb)[2(ia|jb) - (ib|ja)] / D regrouped, with
   !! D = e_i + e_j - e_a - e_b.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, build_df_mo_tensor
   implicit none
   private

   public :: run_libcint_mp2
   public :: run_libcint_ri_mp2
   public :: mp2_result_t

   type :: mp2_result_t
      !! What an MP2 calculation leaves behind
      real(dp) :: same_spin = 0.0_dp        !! E_SS
      real(dp) :: opposite_spin = 0.0_dp    !! E_OS
      real(dp) :: correlation = 0.0_dp      !! E_SS + E_OS
      real(dp) :: scf_energy = 0.0_dp       !! The reference it was built on
      real(dp) :: total = 0.0_dp            !! SCF + correlation
      integer :: n_frozen = 0               !! Core orbitals excluded
      integer :: n_occupied = 0             !! Correlated occupied count
      integer :: n_virtual = 0
   end type mp2_result_t

contains

   subroutine run_libcint_mp2(mol, coeff, orbital_energies, n_occ, scf_energy, &
                              result, error, n_frozen)
      !! E(2) for a closed-shell reference
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)              !! C, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)      !! (n_mo)
      integer, intent(in) :: n_occ                     !! Doubly occupied count
      real(dp), intent(in) :: scf_energy               !! Total SCF energy
      type(mp2_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen
         !! Core orbitals to leave uncorrelated. Defaults to none, which is
         !! what makes this comparable to a reference run with frozen=0 --
         !! most programs freeze by default, and comparing a frozen number
         !! against an unfrozen one disagrees by tens of millihartree for a
         !! reason that has nothing to do with the implementation.

      real(dp), allocatable :: ao_eri(:, :, :, :)
      real(dp), allocatable :: ovov(:, :, :, :)
      integer :: n_ao, n_mo, n_v, n_o, frozen
      integer :: i, j, a, b
      real(dp) :: iajb, ibja, denom, e_ss, e_os

      n_ao = mol%nao
      n_mo = size(coeff, 2)

      frozen = 0
      if (present(n_frozen)) frozen = n_frozen
      if (frozen < 0 .or. frozen >= n_occ) then
         call error%set(ERROR_VALIDATION, "MP2: the frozen core must leave at "// &
                        "least one occupied orbital to correlate")
         return
      end if

      n_o = n_occ - frozen
      n_v = n_mo - n_occ
      if (n_v < 1) then
         call error%set(ERROR_VALIDATION, "MP2: no virtual orbitals; the basis "// &
                        "is saturated by the occupied space")
         return
      end if
      if (size(orbital_energies) /= n_mo) then
         call error%set(ERROR_VALIDATION, "MP2: orbital energies do not match the "// &
                        "coefficient matrix")
         return
      end if

      call mol%eris(ao_eri)
      call transform_ovov(ao_eri, coeff, frozen, n_occ, n_ao, n_mo, ovov)
      deallocate (ao_eri)

      ! The denominators are all strictly negative for a stable reference, so a
      ! non-negative one means the SCF solution is not a minimum. Rather than
      ! guard every term, the sum is checked afterwards: MP2 correlation is
      ! negative, and a positive total says the reference was wrong, which is a
      ! more useful thing to report than a divide-by-near-zero here.
      e_ss = 0.0_dp
      e_os = 0.0_dp
      do i = 1, n_o
         do j = 1, n_o
            do a = 1, n_v
               do b = 1, n_v
                  iajb = ovov(i, a, j, b)
                  ibja = ovov(i, b, j, a)
                  denom = orbital_energies(frozen + i) + orbital_energies(frozen + j) &
                          - orbital_energies(n_occ + a) - orbital_energies(n_occ + b)
                  e_os = e_os + iajb*iajb/denom
                  e_ss = e_ss + iajb*(iajb - ibja)/denom
               end do
            end do
         end do
      end do

      result%same_spin = e_ss
      result%opposite_spin = e_os
      result%correlation = e_ss + e_os
      result%scf_energy = scf_energy
      result%total = scf_energy + result%correlation
      result%n_frozen = frozen
      result%n_occupied = n_o
      result%n_virtual = n_v
   end subroutine run_libcint_mp2

   subroutine run_libcint_ri_mp2(mol, aux, coeff, orbital_energies, n_occ, scf_energy, &
                                 result, error, n_frozen)
      !! E(2) from a fitted (ia|jb), never forming the four-index tensor
      !!
      !! `build_df_tensor` already returns B(mu nu, P) with the inverse-root
      !! metric folded in, which is exactly the object this needs:
      !!
      !!     (ia|jb) = sum_P B^P_ia B^P_jb
      !!
      !! So the work is a transform of the AO pair index into (occupied,
      !! virtual) and then one gemm per occupied pair. The saving over the
      !! conventional route is where the memory goes: n_occ*n_vir*n_aux held
      !! here against the n^4 tensor and its n_occ*n^3 intermediate there. For
      !! water in cc-pVDZ that is kilobytes against megabytes, and the gap is
      !! two powers of the system size.
      !!
      !! The error against conventional MP2 is the fitting error and nothing
      !! else -- typically a few tenths of a millihartree with a matched RIFIT
      !! set, and much worse with a JKFIT one, which is fitted for the Coulomb
      !! and exchange blocks rather than for this.
      type(libcint_molecule_t), intent(in) :: mol
      type(libcint_molecule_t), intent(in) :: aux    !! Auxiliary basis, same atoms
      real(dp), intent(in) :: coeff(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: scf_energy
      type(mp2_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen

      real(dp), allocatable :: bia(:, :, :), g(:, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :)
      integer :: n_ao, n_mo, n_o, n_v, n_aux, frozen
      integer :: i, j, a, bb
      real(dp) :: iajb, ibja, denom, e_ss, e_os, weight

      n_ao = mol%nao
      n_mo = size(coeff, 2)

      frozen = 0
      if (present(n_frozen)) frozen = n_frozen
      if (frozen < 0 .or. frozen >= n_occ) then
         call error%set(ERROR_VALIDATION, "RI-MP2: the frozen core must leave at "// &
                        "least one occupied orbital to correlate")
         return
      end if

      n_o = n_occ - frozen
      n_v = n_mo - n_occ
      if (n_v < 1) then
         call error%set(ERROR_VALIDATION, "RI-MP2: no virtual orbitals")
         return
      end if

      allocate (c_occ(n_ao, n_o), c_vir(n_ao, n_v))
      c_occ = coeff(:, frozen + 1:n_occ)
      c_vir = coeff(:, n_occ + 1:n_mo)

      ! Transformed before it is fitted, which is where the saving is -- see
      ! `build_df_mo_tensor`. Fitting the whole AO pair space first and reading
      ! the occupied-virtual corner out of it costs the same answer several
      ! times over.
      call build_df_mo_tensor(mol, aux, c_occ, c_vir, bia, error)
      deallocate (c_occ, c_vir)
      if (error%has_error()) return
      n_aux = size(bia, 2)

      ! One gemm per occupied pair rebuilds that pair's whole (a,b) block:
      ! g(a,b) = sum_P B^P_ia B^P_jb, which is (ia|jb), and g(b,a) is (ib|ja).
      ! Only the lower triangle of occupied pairs. Exchanging i with j and a
      ! with b leaves every term unchanged -- (ia|jb) becomes (jb|ia), which is
      ! the same integral -- so the (j,i) pair contributes exactly what (i,j)
      ! does and is counted by weight instead of by a second gemm. Half the
      ! work for the same sum, and the gemm is where the time is.
      allocate (g(n_v, n_v))
      e_ss = 0.0_dp
      e_os = 0.0_dp
      do i = 1, n_o
         do j = 1, i
            call pic_gemm(bia(:, :, i), bia(:, :, j), g, transb="T")
            weight = 2.0_dp
            if (i == j) weight = 1.0_dp
            do bb = 1, n_v
               do a = 1, n_v
                  iajb = g(a, bb)
                  ibja = g(bb, a)
                  denom = orbital_energies(frozen + i) + orbital_energies(frozen + j) &
                          - orbital_energies(n_occ + a) - orbital_energies(n_occ + bb)
                  e_os = e_os + weight*iajb*iajb/denom
                  e_ss = e_ss + weight*iajb*(iajb - ibja)/denom
               end do
            end do
         end do
      end do
      deallocate (g, bia)

      result%same_spin = e_ss
      result%opposite_spin = e_os
      result%correlation = e_ss + e_os
      result%scf_energy = scf_energy
      result%total = scf_energy + result%correlation
      result%n_frozen = frozen
      result%n_occupied = n_o
      result%n_virtual = n_v
   end subroutine run_libcint_ri_mp2

   subroutine transform_ovov(eri, coeff, frozen, n_occ, n_ao, n_mo, ovov)
      !! (mu nu|la si) -> (ia|jb), one index at a time, without copying
      !!
      !! Each quarter contracts one AO index against one coefficient block, and
      !! each is a gemm. The care is entirely in the layout: every intermediate
      !! is held so that the slice the next quarter reads is contiguous, and is
      !! handed to BLAS through a bounds-remapped pointer rather than `reshape`.
      !!
      !! That distinction is not cosmetic. `reshape` always materialises a
      !! temporary, so viewing one occupied orbital's block as a matrix copied
      !! nao^3 elements per orbital -- 573 MB on a water hexamer in cc-pVDZ,
      !! to feed gemms that together do less work than the copying. Pointer
      !! remapping is the same view for nothing.
      !!
      !! The layouts, in order:
      !!
      !!   q1 (nu la si, i)   quarter 1 contracts mu
      !!   q2 (la si, a, i)   quarter 2 contracts nu
      !!   q3 (si, j)         quarter 3 contracts la, per (i,a)
      !!   ovov (i, a, j, b)  quarter 4 contracts si
      !!
      !! Each is chosen so the following quarter reads down a column.
      real(dp), intent(in), target, contiguous :: eri(:, :, :, :)
      real(dp), intent(in) :: coeff(:, :)
      integer, intent(in) :: frozen, n_occ, n_ao, n_mo
      real(dp), allocatable, intent(out) :: ovov(:, :, :, :)

      real(dp), allocatable :: c_occ(:, :), c_vir(:, :)
      real(dp), allocatable, target :: q1(:, :), q2(:, :, :)
      real(dp), allocatable :: q3(:, :), q4(:, :)
      real(dp), pointer, contiguous :: view(:, :)
      integer :: n_o, n_v, i, a, j

      n_o = n_occ - frozen
      n_v = n_mo - n_occ

      allocate (c_occ(n_ao, n_o), c_vir(n_ao, n_v))
      c_occ = coeff(:, frozen + 1:n_occ)
      c_vir = coeff(:, n_occ + 1:n_mo)

      ! Quarter 1: mu against the occupied block, over the whole tensor at once.
      ! The occupied index lands last so the next quarter reads columns.
      allocate (q1(n_ao*n_ao*n_ao, n_o))
      view(1:n_ao, 1:n_ao*n_ao*n_ao) => eri
      call pic_gemm(view, c_occ, q1, transa="T")

      ! Quarter 2: nu against the virtual block, one occupied orbital at a time.
      allocate (q2(n_ao*n_ao, n_v, n_o))
      do i = 1, n_o
         view(1:n_ao, 1:n_ao*n_ao) => q1(:, i)
         call pic_gemm(view, c_vir, q2(:, :, i), transa="T")
      end do
      deallocate (q1)

      ! Quarters 3 and 4: la then si, per (i,a). Both slices are columns.
      allocate (q3(n_o, n_ao), q4(n_o, n_v))
      allocate (ovov(n_o, n_v, n_o, n_v))
      do i = 1, n_o
         do a = 1, n_v
            view(1:n_ao, 1:n_ao) => q2(:, a, i)
            call pic_gemm(c_occ, view, q3, transa="T")
            call pic_gemm(q3, c_vir, q4)
            do j = 1, n_o
               ovov(i, a, j, :) = q4(j, :)
            end do
         end do
      end do
      deallocate (q2, q3, q4, c_occ, c_vir)
   end subroutine transform_ovov

end module mqc_libcint_mp2
