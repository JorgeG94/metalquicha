!! Closed-shell MP2 from stored AO integrals
module mqc_libcint_mp2
   !! Conventional restricted MP2: the reference the fitted version is checked
   !! against.
   !!
   !! **The transformation is the whole of it.** One index at a time it is four
   !! N^5 steps, each a gemm:
   !!
   !!     (mu nu|la si) --C_occ--> (i nu|la si) --C_vir--> (ia|la si)
   !!                   --C_occ--> (ia|j si)   --C_vir--> (ia|jb)
   !!
   !! Memory is the price, and it is why this is a reference rather than a
   !! method: the first intermediate is n_occ * n^3, and the AO tensor it reads
   !! is n^4 on top of that.
   !!
   !! The energy is decomposed as it is accumulated, since the two spin cases
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
   use mqc_timing, only: timing_report_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, build_df_mo_tensor, pair_index
   implicit none
   private

   public :: run_libcint_mp2
   public :: run_libcint_ri_mp2
   public :: run_libcint_ump2
   public :: run_libcint_uri_mp2
   public :: mp2_result_t
   ! Exported for coupled cluster, which needs every MO block rather than just
   ! (ia|jb).
   public :: transform_block
   public :: transform_ovov   !! The MP2 gradient builds its own amplitudes from this

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
         !! Core orbitals to leave uncorrelated. Defaults to none, where most
         !! programs freeze by default -- a frozen number against an unfrozen
         !! one disagrees by tens of millihartree.

      real(dp), allocatable :: ao_eri(:, :)
      real(dp), allocatable :: ovov(:, :, :, :)
      integer :: n_ao, n_mo, n_v, n_o, frozen
      integer :: i, j, a, b
      real(dp) :: iajb, ibja, denom, e_ss, e_os

      type(timing_report_t) :: clk

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

      call clk%start()
      call clk%begin("AO integrals")
      call mol%eris_packed(ao_eri)
      call clk%lap()
      call clk%begin("AO->MO transform")
      call transform_ovov(ao_eri, coeff, frozen, n_occ, n_ao, n_mo, ovov)
      deallocate (ao_eri)
      call clk%lap()

      ! The denominators are all strictly negative for a stable reference, so a
      ! non-negative one means the SCF solution is not a minimum.
      ! TODO(mqc): nothing guards the division and nothing checks the sum
      ! afterwards either -- no caller of `mp2_result_t` tests the sign of the
      ! correlation -- so a vanishing denominator returns an infinity, and a
      ! positive one a silently wrong energy.
      call clk%begin("energy denominators")
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

      call clk%lap()
      call clk%finish()
      call clk%report("MP2")

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
                                 result, error, n_frozen, b_ao_in)
      !! E(2) from a fitted (ia|jb), never forming the four-index tensor
      !!
      !! `build_df_tensor` returns B(mu nu, P) with the inverse-root metric
      !! folded in, which is exactly the object this needs:
      !!
      !!     (ia|jb) = sum_P B^P_ia B^P_jb
      !!
      !! So the work is a transform of the AO pair index into (occupied,
      !! virtual) and then one gemm per occupied pair, holding n_occ*n_vir*n_aux
      !! against the conventional route's n^4 tensor and n_occ*n^3 intermediate.
      !!
      !! The error against conventional MP2 is the fitting error and nothing
      !! else -- a few tenths of a millihartree with a matched RIFIT set, and
      !! much worse with a JKFIT one, which is fitted for the Coulomb and
      !! exchange blocks rather than for this.
      type(libcint_molecule_t), intent(in) :: mol
      type(libcint_molecule_t), intent(in) :: aux    !! Auxiliary basis, same atoms
      real(dp), intent(in) :: coeff(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: scf_energy
      type(mp2_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen
      real(dp), intent(in), optional :: b_ao_in(:, :)
         !! The AO-basis fitted tensor a fitted SCF already built with this same
         !! auxiliary basis. Given it, this reduces to the MO transform; see
         !! `build_df_mo_block`.

      real(dp), allocatable :: bia(:, :, :), g(:, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :)
      integer :: n_ao, n_mo, n_o, n_v, n_aux, frozen
      integer :: i, j, a, bb
      real(dp) :: iajb, ibja, denom, e_ss, e_os, weight

      type(timing_report_t) :: clk

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

      ! Transformed before it is fitted, which is where the saving is: fitting
      ! the whole AO pair space first and reading the occupied-virtual corner
      ! out of it costs the same answer several times over.
      call clk%start()
      call clk%begin("B tensor (3c/2c fit)")
      ! Energy only: `bia` is contracted with another `bia` and nothing
      ! else, so the cheap factor is safe here. The RI-MP2 *gradient*
      ! builds its own and must not.
      call build_df_mo_tensor(mol, aux, c_occ, c_vir, bia, error, fast_factor=.true., &
                              b_ao_in=b_ao_in)
      deallocate (c_occ, c_vir)
      if (error%has_error()) return
      call clk%lap()
      n_aux = size(bia, 2)

      ! One gemm per occupied pair rebuilds that pair's whole (a,b) block:
      ! g(a,b) = sum_P B^P_ia B^P_jb, which is (ia|jb), and g(b,a) is (ib|ja).
      ! Only the lower triangle of occupied pairs: exchanging i with j and a
      ! with b leaves every term unchanged, so the (j,i) pair is counted by
      ! weight instead of by a second gemm.
      !
      ! `g` is per-thread scratch, allocated inside the region, which is why the
      ! region is explicit instead of a bare `parallel do`. `schedule(dynamic)`
      ! because the inner loop runs to `i`, so the last occupied orbital does
      ! n_o times the work of the first.
      call clk%begin("gemm + denominators")
      e_ss = 0.0_dp
      e_os = 0.0_dp
      !$omp parallel default(none) &
      !$omp    shared(bia, orbital_energies, n_o, n_v, n_occ, frozen) &
      !$omp    private(i, j, a, bb, g, iajb, ibja, denom, weight) &
      !$omp    reduction(+:e_ss, e_os)
      allocate (g(n_v, n_v))
      !$omp do schedule(dynamic)
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
      !$omp end do
      deallocate (g)
      !$omp end parallel
      deallocate (bia)
      call clk%lap()
      call clk%finish()
      call clk%report("RI-MP2")

      result%same_spin = e_ss
      result%opposite_spin = e_os
      result%correlation = e_ss + e_os
      result%scf_energy = scf_energy
      result%total = scf_energy + result%correlation
      result%n_frozen = frozen
      result%n_occupied = n_o
      result%n_virtual = n_v
   end subroutine run_libcint_ri_mp2

   !
   ! THE UNRESTRICTED PAIR. Both routines below compute
   !
   !     E(2) = 1/4 sum_(ijab in a) |<ij||ab>|^2 / D
   !          + 1/4 sum_(ijab in b) |<ij||ab>|^2 / D
   !          +     sum_(ia in a, jb in b) (ia|jb)^2 / D
   !
   ! with D = e_i + e_j - e_a - e_b, and store the two like-spin terms together
   ! in `same_spin` and the mixed one in `opposite_spin`, which is what those
   ! fields already mean in the restricted result.
   !
   ! The like-spin terms are written as 1/2 sum (ia|jb)[(ia|jb) - (ib|ja)]/D
   ! rather than the equal 1/4 sum [(ia|jb) - (ib|ja)]^2/D: half the arithmetic
   ! and the same shape as the restricted loop above it.
   !
   ! Frozen orbitals are counted per spin, the same number from each. For a
   ! doublet that leaves one more correlated beta hole than alpha, which is
   ! what freezing a core rather than an electron count means.
   !
   subroutine run_libcint_ump2(mol, coeff_a, coeff_b, eps_a, eps_b, &
                               n_occ_a, n_occ_b, scf_energy, result, error, n_frozen)
      !! E(2) for an unrestricted reference, from the four-index transform
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff_a(:, :), coeff_b(:, :)   !! (n_ao, n_mo) per spin
      real(dp), intent(in) :: eps_a(:), eps_b(:)
      integer, intent(in) :: n_occ_a, n_occ_b
      real(dp), intent(in) :: scf_energy
      type(mp2_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen

      real(dp), allocatable :: ao_eri(:, :)
      real(dp), allocatable :: ovov_aa(:, :, :, :), ovov_bb(:, :, :, :), ovov_ab(:, :, :, :)
      real(dp), allocatable :: ca_o(:, :), ca_v(:, :), cb_o(:, :), cb_v(:, :)
      integer :: n_ao, n_mo, frozen, n_oa, n_ob, n_va, n_vb
      real(dp) :: e_aa, e_bb, e_ab
      type(timing_report_t) :: clk

      call ump2_dimensions(mol, coeff_a, coeff_b, eps_a, eps_b, n_occ_a, n_occ_b, &
                           n_frozen, n_ao, n_mo, frozen, n_oa, n_ob, n_va, n_vb, error)
      if (error%has_error()) return

      ca_o = coeff_a(:, frozen + 1:n_occ_a)
      ca_v = coeff_a(:, n_occ_a + 1:n_mo)
      cb_o = coeff_b(:, frozen + 1:n_occ_b)
      cb_v = coeff_b(:, n_occ_b + 1:n_mo)

      call clk%start()
      call clk%begin("AO integrals")
      call mol%eris_packed(ao_eri)
      call clk%lap()

      call clk%begin("AO->MO transform (3 blocks)")
      ! Three blocks, not one. The like-spin pair needs its own coefficients on
      ! both electrons; the mixed block needs alpha on one and beta on the
      ! other and is NOT antisymmetrised -- two electrons of different spin are
      ! already distinguishable, so there is no exchange term to subtract.
      call transform_block(ao_eri, ca_o, ca_v, ca_o, ca_v, ovov_aa)
      call transform_block(ao_eri, cb_o, cb_v, cb_o, cb_v, ovov_bb)
      call transform_block(ao_eri, ca_o, ca_v, cb_o, cb_v, ovov_ab)
      deallocate (ao_eri, ca_o, ca_v, cb_o, cb_v)
      call clk%lap()

      call clk%begin("energy denominators")
      e_aa = like_spin_energy(ovov_aa, eps_a, frozen, n_occ_a, n_oa, n_va)
      e_bb = like_spin_energy(ovov_bb, eps_b, frozen, n_occ_b, n_ob, n_vb)
      e_ab = mixed_spin_energy(ovov_ab, eps_a, eps_b, frozen, n_occ_a, n_occ_b, &
                               n_oa, n_ob, n_va, n_vb)
      deallocate (ovov_aa, ovov_bb, ovov_ab)
      call clk%lap()
      call clk%finish()
      call clk%report("UMP2")

      call fill_ump2_result(result, e_aa + e_bb, e_ab, scf_energy, frozen, &
                            n_oa + n_ob, n_va + n_vb)
   end subroutine run_libcint_ump2

   subroutine run_libcint_uri_mp2(mol, aux, coeff_a, coeff_b, eps_a, eps_b, &
                                  n_occ_a, n_occ_b, scf_energy, result, error, n_frozen)
      !! E(2) for an unrestricted reference, from a fitted (ia|jb)
      !!
      !! One B tensor per spin, both fitted against the same auxiliary basis so
      !! the mixed block is a gemm between them rather than anything new:
      !!
      !!     (ia|JB) = sum_P B(a)^P_ia B(b)^P_JB
      !!
      !! The like-spin blocks reuse the occupied-pair weighting of the
      !! restricted routine, so the lower triangle is enough. **The mixed block
      !! gets no such weighting**: i is alpha and j is beta, so (i,j) and (j,i)
      !! are different pairs of different spins and the full rectangle is the
      !! sum.
      type(libcint_molecule_t), intent(in) :: mol
      type(libcint_molecule_t), intent(in) :: aux
      real(dp), intent(in) :: coeff_a(:, :), coeff_b(:, :)
      real(dp), intent(in) :: eps_a(:), eps_b(:)
      integer, intent(in) :: n_occ_a, n_occ_b
      real(dp), intent(in) :: scf_energy
      type(mp2_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen

      real(dp), allocatable :: b_a(:, :, :), b_b(:, :, :), g(:, :)
      real(dp), allocatable :: ca_o(:, :), ca_v(:, :), cb_o(:, :), cb_v(:, :)
      integer :: n_ao, n_mo, frozen, n_oa, n_ob, n_va, n_vb
      integer :: i, j, a, b
      real(dp) :: e_aa, e_bb, e_ab, iajb, denom
      type(timing_report_t) :: clk

      call ump2_dimensions(mol, coeff_a, coeff_b, eps_a, eps_b, n_occ_a, n_occ_b, &
                           n_frozen, n_ao, n_mo, frozen, n_oa, n_ob, n_va, n_vb, error)
      if (error%has_error()) return

      ca_o = coeff_a(:, frozen + 1:n_occ_a)
      ca_v = coeff_a(:, n_occ_a + 1:n_mo)
      cb_o = coeff_b(:, frozen + 1:n_occ_b)
      cb_v = coeff_b(:, n_occ_b + 1:n_mo)

      call clk%start()
      call clk%begin("B tensors (3c/2c fit, both spins)")
      call build_df_mo_tensor(mol, aux, ca_o, ca_v, b_a, error, fast_factor=.true.)
      if (error%has_error()) return
      call build_df_mo_tensor(mol, aux, cb_o, cb_v, b_b, error, fast_factor=.true.)
      if (error%has_error()) return
      deallocate (ca_o, ca_v, cb_o, cb_v)
      call clk%lap()

      call clk%begin("gemm + denominators")
      e_aa = ri_like_spin_energy(b_a, eps_a, frozen, n_occ_a, n_oa, n_va)
      e_bb = ri_like_spin_energy(b_b, eps_b, frozen, n_occ_b, n_ob, n_vb)

      allocate (g(n_va, n_vb))
      e_ab = 0.0_dp
      do i = 1, n_oa
         do j = 1, n_ob
            call pic_gemm(b_a(:, :, i), b_b(:, :, j), g, transb="T")
            do b = 1, n_vb
               do a = 1, n_va
                  iajb = g(a, b)
                  denom = eps_a(frozen + i) + eps_b(frozen + j) &
                          - eps_a(n_occ_a + a) - eps_b(n_occ_b + b)
                  e_ab = e_ab + iajb*iajb/denom
               end do
            end do
         end do
      end do
      deallocate (g, b_a, b_b)
      call clk%lap()
      call clk%finish()
      call clk%report("URI-MP2")

      call fill_ump2_result(result, e_aa + e_bb, e_ab, scf_energy, frozen, &
                            n_oa + n_ob, n_va + n_vb)
   end subroutine run_libcint_uri_mp2

   subroutine ump2_dimensions(mol, coeff_a, coeff_b, eps_a, eps_b, n_occ_a, n_occ_b, &
                              n_frozen, n_ao, n_mo, frozen, n_oa, n_ob, n_va, n_vb, error)
      !! The shape checks both unrestricted routines need, in one place
      !!
      !! Checked per spin rather than once: a doublet has a different occupied
      !! count in each, and the beta channel is the one that runs out first.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff_a(:, :), coeff_b(:, :)
      real(dp), intent(in) :: eps_a(:), eps_b(:)
      integer, intent(in) :: n_occ_a, n_occ_b
      integer, intent(in), optional :: n_frozen
      integer, intent(out) :: n_ao, n_mo, frozen, n_oa, n_ob, n_va, n_vb
      type(error_t), intent(inout) :: error

      n_ao = mol%nao
      n_mo = size(coeff_a, 2)
      frozen = 0
      if (present(n_frozen)) frozen = n_frozen

      if (size(coeff_b, 2) /= n_mo) then
         call error%set(ERROR_VALIDATION, "UMP2: the alpha and beta coefficient "// &
                        "matrices have different numbers of orbitals")
         return
      end if
      if (size(eps_a) /= n_mo .or. size(eps_b) /= n_mo) then
         call error%set(ERROR_VALIDATION, "UMP2: orbital energies do not match the "// &
                        "coefficient matrices")
         return
      end if
      if (frozen < 0 .or. frozen >= n_occ_a .or. frozen >= n_occ_b) then
         call error%set(ERROR_VALIDATION, "UMP2: the frozen core must leave at least "// &
                        "one occupied orbital in each spin to correlate")
         return
      end if

      n_oa = n_occ_a - frozen
      n_ob = n_occ_b - frozen
      n_va = n_mo - n_occ_a
      n_vb = n_mo - n_occ_b
      if (n_va < 1 .or. n_vb < 1) then
         call error%set(ERROR_VALIDATION, "UMP2: no virtual orbitals in one spin; "// &
                        "the basis is saturated by the occupied space")
         return
      end if
   end subroutine ump2_dimensions

   subroutine fill_ump2_result(result, same, opposite, scf_energy, frozen, n_o, n_v)
      !! One place that decides what the result fields mean for an open shell
      type(mp2_result_t), intent(out) :: result
      real(dp), intent(in) :: same, opposite, scf_energy
      integer, intent(in) :: frozen, n_o, n_v
      result%same_spin = same
      result%opposite_spin = opposite
      result%correlation = same + opposite
      result%scf_energy = scf_energy
      result%total = scf_energy + result%correlation
      result%n_frozen = frozen
      ! Summed over spins, so these count spin orbitals where the restricted
      ! routines count spatial ones.
      result%n_occupied = n_o
      result%n_virtual = n_v
   end subroutine fill_ump2_result

   pure function like_spin_energy(ovov, eps, frozen, n_occ, n_o, n_v) result(e)
      !! 1/2 sum (ia|jb)[(ia|jb) - (ib|ja)] / D, one spin against itself
      real(dp), intent(in) :: ovov(:, :, :, :)
      real(dp), intent(in) :: eps(:)
      integer, intent(in) :: frozen, n_occ, n_o, n_v
      real(dp) :: e
      integer :: i, j, a, b
      real(dp) :: iajb, ibja, denom
      e = 0.0_dp
      do i = 1, n_o
         do j = 1, n_o
            do a = 1, n_v
               do b = 1, n_v
                  iajb = ovov(i, a, j, b)
                  ibja = ovov(i, b, j, a)
                  denom = eps(frozen + i) + eps(frozen + j) &
                          - eps(n_occ + a) - eps(n_occ + b)
                  e = e + 0.5_dp*iajb*(iajb - ibja)/denom
               end do
            end do
         end do
      end do
   end function like_spin_energy

   pure function mixed_spin_energy(ovov, eps_a, eps_b, frozen, n_occ_a, n_occ_b, &
                                   n_oa, n_ob, n_va, n_vb) result(e)
      !! sum (ia|JB)^2 / D, alpha on one electron and beta on the other
      real(dp), intent(in) :: ovov(:, :, :, :)
      real(dp), intent(in) :: eps_a(:), eps_b(:)
      integer, intent(in) :: frozen, n_occ_a, n_occ_b, n_oa, n_ob, n_va, n_vb
      real(dp) :: e
      integer :: i, j, a, b
      real(dp) :: iajb, denom
      e = 0.0_dp
      do i = 1, n_oa
         do j = 1, n_ob
            do a = 1, n_va
               do b = 1, n_vb
                  iajb = ovov(i, a, j, b)
                  denom = eps_a(frozen + i) + eps_b(frozen + j) &
                          - eps_a(n_occ_a + a) - eps_b(n_occ_b + b)
                  e = e + iajb*iajb/denom
               end do
            end do
         end do
      end do
   end function mixed_spin_energy

   function ri_like_spin_energy(bia, eps, frozen, n_occ, n_o, n_v) result(e)
      !! The like-spin term from a fitted B, lower occupied triangle only
      ! TODO(mqc): serial, where `run_libcint_ri_mp2` threads the identical
      ! occupied-pair loop over the same gemm. The unrestricted fitted path
      ! runs both spin blocks and the mixed one on a single thread.
      real(dp), intent(in) :: bia(:, :, :)
      real(dp), intent(in) :: eps(:)
      integer, intent(in) :: frozen, n_occ, n_o, n_v
      real(dp) :: e
      real(dp), allocatable :: g(:, :)
      integer :: i, j, a, b
      real(dp) :: iajb, ibja, denom, weight
      allocate (g(n_v, n_v))
      e = 0.0_dp
      do i = 1, n_o
         do j = 1, i
            call pic_gemm(bia(:, :, i), bia(:, :, j), g, transb="T")
            weight = 2.0_dp
            if (i == j) weight = 1.0_dp
            do b = 1, n_v
               do a = 1, n_v
                  iajb = g(a, b)
                  ibja = g(b, a)
                  denom = eps(frozen + i) + eps(frozen + j) &
                          - eps(n_occ + a) - eps(n_occ + b)
                  e = e + 0.5_dp*weight*iajb*(iajb - ibja)/denom
               end do
            end do
         end do
      end do
      deallocate (g)
   end function ri_like_spin_energy

   subroutine transform_ovov(eri, coeff, frozen, n_occ, n_ao, n_mo, ovov)
      !! (pq|rs) over packed AO pairs -> (ia|jb)
      real(dp), intent(in) :: eri(:, :)      !! (n_pair, n_pair), see `pair_index`
      real(dp), intent(in) :: coeff(:, :)
      integer, intent(in) :: frozen, n_occ, n_ao, n_mo
      real(dp), allocatable, intent(out) :: ovov(:, :, :, :)

      real(dp), allocatable :: c_occ(:, :), c_vir(:, :)

      allocate (c_occ(n_ao, n_occ - frozen), c_vir(n_ao, n_mo - n_occ))
      c_occ = coeff(:, frozen + 1:n_occ)
      c_vir = coeff(:, n_occ + 1:n_mo)
      call transform_block(eri, c_occ, c_vir, c_occ, c_vir, ovov)
      deallocate (c_occ, c_vir)
   end subroutine transform_ovov

   subroutine transform_block(eri, c1, c2, c3, c4, out)
      !! (pq|rs) over packed AO pairs -> (12|34) for any four coefficient blocks
      !!
      !! With the AO pairs packed, the four quarter transforms fall into two
      !! identical halves. One column of the packed tensor is a whole AO pair
      !! index held fixed while the other runs, so unpacking it into an nao by
      !! nao matrix and taking C_occ^T M C_vir transforms that pair to (ia) in
      !! two gemms. Doing it once down the columns and once down the rows is the
      !! whole transformation:
      !!
      !!     (pq|rs) --unpack pq--> (ia|rs) --unpack rs--> (ia|jb)
      !!
      !! Which blocks come out is entirely the caller's four coefficient
      !! matrices: C_occ/C_vir twice over gives (ia|jb) for MP2, and C_act four
      !! times over gives the whole active MO tensor, which is what coupled
      !! cluster needs.
      real(dp), intent(in) :: eri(:, :)      !! (n_pair, n_pair), see `pair_index`
      real(dp), intent(in) :: c1(:, :), c2(:, :)   !! Bra pair, left and right
      real(dp), intent(in) :: c3(:, :), c4(:, :)   !! Ket pair, left and right
      real(dp), allocatable, intent(out) :: out(:, :, :, :)

      real(dp), allocatable :: half(:, :)
      real(dp), allocatable :: m(:, :), tmp(:, :), blk(:, :)
      integer :: n_ao, n1, n2, n3, n4, n_pair, pq, l, mm, n, o, p, q

      n_ao = size(c1, 1)
      n1 = size(c1, 2)
      n2 = size(c2, 2)
      n3 = size(c3, 2)
      n4 = size(c4, 2)
      n_pair = n_ao*(n_ao + 1)/2

      ! First half: every packed bra pair transformed to (12), leaving (12|rs).
      ! Held with the packed index first so the second half reads columns.
      allocate (half(n_pair, n1*n2))
      !$omp parallel default(none) &
      !$omp    shared(eri, half, c1, c2, n_pair, n_ao, n1, n2) &
      !$omp    private(pq, p, q, l, mm, m, tmp, blk)
      allocate (m(n_ao, n_ao), tmp(n1, n_ao), blk(n1, n2))
      !$omp do schedule(static)
      do pq = 1, n_pair
         do q = 1, n_ao
            do p = 1, n_ao
               m(p, q) = eri(pair_index(p, q), pq)
            end do
         end do
         call pic_gemm(c1, m, tmp, transa="T")
         call pic_gemm(tmp, c2, blk)
         do mm = 1, n2
            do l = 1, n1
               half(pq, (mm - 1)*n1 + l) = blk(l, mm)
            end do
         end do
      end do
      !$omp end do
      deallocate (m, tmp, blk)
      !$omp end parallel

      ! Second half: the same operation on the ket pair, for each (12).
      allocate (out(n1, n2, n3, n4))
      !$omp parallel default(none) &
      !$omp    shared(half, out, c3, c4, n_pair, n_ao, n1, n2, n3, n4) &
      !$omp    private(pq, p, q, l, mm, n, o, m, tmp, blk)
      allocate (m(n_ao, n_ao), tmp(n3, n_ao), blk(n3, n4))
      !$omp do schedule(static) collapse(2)
      do mm = 1, n2
         do l = 1, n1
            do q = 1, n_ao
               do p = 1, n_ao
                  m(p, q) = half(pair_index(p, q), (mm - 1)*n1 + l)
               end do
            end do
            call pic_gemm(c3, m, tmp, transa="T")
            call pic_gemm(tmp, c4, blk)
            do o = 1, n4
               do n = 1, n3
                  out(l, mm, n, o) = blk(n, o)
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (m, tmp, blk)
      !$omp end parallel

      deallocate (half)
   end subroutine transform_block

end module mqc_libcint_mp2
