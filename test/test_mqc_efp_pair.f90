!! A molecule spanning two fragments, checked by its own block structure
module test_mqc_efp_pair
   !! The inter-fragment overlaps that exchange repulsion, charge transfer and the
   !! dispersion damping need come from a molecule covering both fragments, built out
   !! of the basis the potential itself carries. What makes that checkable without a
   !! reference is the block structure of the overlap it produces.
   !!
   !! With `a`'s atoms first, the leading `n_ao_a` rows and columns of the pair
   !! overlap must be *exactly* `a`'s own overlap, computed from `a` alone. That is a
   !! strong test: it fails if the basis was rebuilt with GAMESS's primitive
   !! normalization still folded in, if a shell's angular momentum was mistaken, if an
   !! `L` shell was split wrongly, or if the atoms and their bases were paired up out
   !! of order. All of those otherwise produce a plausible matrix.
   !!
   !! The off-diagonal block is what the physics wants, and there is nothing here to
   !! check it against yet -- so what is asserted about it is that it is not zero, and
   !! that it shrinks as the fragments separate. A block that came back empty would
   !! pass every diagonal check above.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_efp_potential, only: efp_potential_t, make_efp_potential, &
                                write_efp_potential
   use mqc_efp_read, only: efp_fragment_t, read_efp_potential
   use mqc_efp_pair, only: two_fragment_molecule, fragment_molecule, fragment_lmo, &
                           exchange_repulsion
   use pic_blas_interfaces, only: pic_gemm
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_efp_pair_tests

   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp

contains

   subroutine collect_mqc_efp_pair_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("efp_pair_diagonal_blocks", test_blocks), &
                  new_unittest("efp_pair_coupling_falls_off", test_falloff), &
                  new_unittest("efp_lmo_orthonormal", test_lmo), &
                  new_unittest("efp_exchange_repulsion", test_xr) &
                  ]
   end subroutine collect_mqc_efp_pair_tests

   subroutine water_fragment(frag, err)
      type(efp_fragment_t), intent(out) :: frag
      type(error_t), intent(inout) :: err

      type(efp_potential_t) :: pot
      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)
      character(len=*), parameter :: path = "test_efp_pair.efp"
      integer :: unit, stat

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.00000000000000_dp, 0.00000000009155_dp, 0.10077199490609_dp, &
                   0.00000000000000_dp, 0.77250895271063_dp, -0.46780199741728_dp, &
                   0.00000000000000_dp, -0.77250895280218_dp, -0.46780199748881_dp], &
                  [3, 3])*ANG
      call make_efp_potential(z, symbols, c, "6-31g*", "WATER", pot, err)
      if (err%has_error()) return
      call write_efp_potential(pot, path, err)
      if (err%has_error()) return
      call read_efp_potential(path, frag, err)
      call pot%destroy()
      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
   end subroutine water_fragment

   subroutine test_blocks(error)
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: a, b
      type(libcint_molecule_t) :: pair, single
      type(error_t) :: err
      real(dp), allocatable :: s_pair(:, :), s_single(:, :)
      real(dp) :: c(3, 3), offset(3)
      integer :: z(3), n_ao_a, n
      character(len=2) :: symbols(3)

      call water_fragment(a, err)
      call check(error,.not. err%has_error(), "building the fragment failed")
      if (allocated(error)) return
      call water_fragment(b, err)
      if (allocated(error)) return

      offset = [3.0_dp*ANG, 0.0_dp, 0.0_dp]
      call two_fragment_molecule(a, b, [0.0_dp, 0.0_dp, 0.0_dp], offset, pair, &
                                 n_ao_a, err)
      call check(error,.not. err%has_error(), "building the pair molecule failed")
      if (allocated(error)) return

      ! One water in 6-31G* is 19 Cartesian functions, so the pair must be 38 and the
      ! split must fall at 19. Hard-coded because both are properties of the basis
      ! rather than of this code, and getting either wrong means the block boundary
      ! is in the wrong place while every matrix still looks reasonable.
      call check(error, n_ao_a == 19, "the first fragment should span 19 functions")
      if (allocated(error)) return
      call check(error, pair%nao == 38, "the pair should span 38 functions")
      if (allocated(error)) return

      ! The same water on its own, from the basis file by name rather than from the
      ! potential -- so this compares the reconstructed basis against the original.
      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.00000000000000_dp, 0.00000000009155_dp, 0.10077199490609_dp, &
                   0.00000000000000_dp, 0.77250895271063_dp, -0.46780199741728_dp, &
                   0.00000000000000_dp, -0.77250895280218_dp, -0.46780199748881_dp], &
                  [3, 3])*ANG
      call build_libcint_molecule(z, symbols, c, "6-31g*", single, err)
      call check(error,.not. err%has_error(), "building the single molecule failed")
      if (allocated(error)) return

      call pair%overlap(s_pair)
      call single%overlap(s_single)
      n = n_ao_a
      ! Tolerance set by the file, not the arithmetic: a potential carries eight
      ! decimals on a contraction coefficient, so a basis recovered from it differs
      ! from the one it was written from at about 1e-8, and the overlap inherits that.
      call check(error, maxval(abs(s_pair(1:n, 1:n) - s_single)) < 1.0e-7_dp, &
                 "the pair's first diagonal block is not the fragment's own overlap")
      if (allocated(error)) return
      ! And the second fragment's block, which also tests that the offset was applied
      ! to its atoms and not to its basis.
      call check(error, maxval(abs(s_pair(n + 1:2*n, n + 1:2*n) - s_single)) < 1.0e-7_dp, &
                 "the pair's second diagonal block is not the fragment's own overlap")
      if (allocated(error)) return

      ! The block the physics wants: present, and small, at 3 Angstrom.
      call check(error, maxval(abs(s_pair(1:n, n + 1:2*n))) > 1.0e-8_dp, &
                 "the inter-fragment overlap block is empty")
      if (allocated(error)) return
      ! The measured peak here is 1.0e-01, from the diffuse hydrogen s functions at
      ! 5.67 Bohr. The bound is set loosely above it: it is a guard against a
      ! catastrophically wrong build -- two fragments placed on top of each other give
      ! order one -- and not a reference value. An earlier version guessed 0.1 and
      ! failed by a fraction of a percent, which is what guessing bounds earns.
      call check(error, maxval(abs(s_pair(1:n, n + 1:2*n))) < 0.5_dp, &
                 "the inter-fragment overlap is implausibly large at 3 Angstrom")

      call a%destroy()
      call b%destroy()
      call pair%destroy()
      call single%destroy()
   end subroutine test_blocks

   subroutine test_falloff(error)
      !! The coupling must fall as the fragments separate
      !!
      !! Overlap decays exponentially with distance, so this is a weak numerical claim
      !! and a strong structural one: it fails if the offset never reached the second
      !! fragment's coordinates, which would leave both copies on top of each other
      !! and the block identical at every separation.
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: a, b
      type(error_t) :: err
      real(dp) :: near, far

      call water_fragment(a, err)
      call check(error,.not. err%has_error(), "building the fragment failed")
      if (allocated(error)) return
      call water_fragment(b, err)
      if (allocated(error)) return

      near = coupling(a, b, 3.0_dp, err)
      call check(error,.not. err%has_error(), "the near pair failed")
      if (allocated(error)) return
      far = coupling(a, b, 6.0_dp, err)
      call check(error,.not. err%has_error(), "the far pair failed")
      if (allocated(error)) return

      call check(error, far < near, "the inter-fragment overlap did not fall off")
      if (allocated(error)) return
      ! An order of magnitude over three Angstrom is a modest ask of something
      ! exponential, and it fails outright if the offset was never applied.
      call check(error, far < 0.1_dp*near, &
                 "the inter-fragment overlap barely changed with separation")

      call a%destroy()
      call b%destroy()
   end subroutine test_falloff

   subroutine test_xr(error)
      !! Pauli exchange repulsion against GAMESS
      !!
      !! The term that pinned the valence-charge convention: the potential in the third
      !! term is built from each fragment's *valence* nuclear charges, not its full
      !! ones, because the localized orbitals a potential carries are valence only and
      !! the core has to appear as screening of its own nucleus. With full charges this
      !! comes out -0.002949404 against GAMESS's -0.001172851 -- a factor of 2.5, not a
      !! sign or a scale, which is what made it findable.
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: a, b
      type(error_t) :: err
      real(dp) :: e

      call water_fragment(a, err)
      call water_fragment(b, err)
      call check(error,.not. err%has_error(), "building the fragments failed")
      if (allocated(error)) return
      e = exchange_repulsion(a, b, [0.0_dp, 0.0_dp, 0.0_dp], &
                             [3.0_dp*ANG, 0.0_dp, 0.0_dp], err)
      call check(error,.not. err%has_error(), "exchange repulsion failed")
      if (allocated(error)) return
      ! GAMESS's own exchange repulsion for this dimer, every digit it prints. It
      ! computes this from data the potential supplies rather than from anything of its
      ! own, so an exact match is the right expectation and a near one would mean a
      ! wrong factor somewhere.
      call check(error, e, -0.001172851_dp, thr=1.0e-9_dp, &
                 message="exchange repulsion disagrees with GAMESS")
      call a%destroy()
      call b%destroy()
   end subroutine test_xr

   subroutine test_lmo(error)
      !! The localized orbitals come back orthonormal, which pins the AO conversion
      !!
      !! A potential stores its orbitals in GAMESS's AO order, so using them with our
      !! integrals means undoing a permutation and two normalizations. Applied in the
      !! wrong direction the result is still a plausible matrix, and every downstream
      !! energy would be quietly wrong -- so the test is a property no wrong
      !! permutation can satisfy: `C^T S C` must be the identity, because localized
      !! orbitals are orthonormal by construction.
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: a
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: lmo(:, :), s(:, :), sc(:, :), gram(:, :)
      integer :: i, j, n_lmo
      real(dp) :: worst

      call water_fragment(a, err)
      call check(error,.not. err%has_error(), "building the fragment failed")
      if (allocated(error)) return

      call fragment_molecule(a, [0.0_dp, 0.0_dp, 0.0_dp], mol, err)
      call check(error,.not. err%has_error(), "building the molecule failed")
      if (allocated(error)) return
      call fragment_lmo(a, mol, lmo, err)
      call check(error,.not. err%has_error(), "converting the orbitals failed")
      if (allocated(error)) return

      n_lmo = a%n_lmo_proj
      call mol%overlap(s)
      allocate (sc(mol%nao, n_lmo), gram(n_lmo, n_lmo))
      call pic_gemm(s, lmo, sc)
      call pic_gemm(lmo, sc, gram, transa="T")

      worst = 0.0_dp
      do j = 1, n_lmo
         do i = 1, n_lmo
            if (i == j) then
               worst = max(worst, abs(gram(i, j) - 1.0_dp))
            else
               worst = max(worst, abs(gram(i, j)))
            end if
         end do
      end do
      ! The file carries eight decimals on an orbital coefficient, so this is the
      ! format's precision rather than the arithmetic's.
      call check(error, worst < 1.0e-6_dp, &
                 "the localized orbitals are not orthonormal, so the AO order "// &
                 "conversion is wrong")

      call a%destroy()
      call mol%destroy()
   end subroutine test_lmo

   function coupling(a, b, separation, err) result(peak)
      type(efp_fragment_t), intent(in) :: a, b
      real(dp), intent(in) :: separation
      type(error_t), intent(inout) :: err
      real(dp) :: peak

      type(libcint_molecule_t) :: pair
      real(dp), allocatable :: s(:, :)
      integer :: n_ao_a

      peak = 0.0_dp
      call two_fragment_molecule(a, b, [0.0_dp, 0.0_dp, 0.0_dp], &
                                 [separation*ANG, 0.0_dp, 0.0_dp], pair, n_ao_a, err)
      if (err%has_error()) return
      call pair%overlap(s)
      peak = maxval(abs(s(1:n_ao_a, n_ao_a + 1:pair%nao)))
      call pair%destroy()
   end function coupling

end module test_mqc_efp_pair

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_efp_pair, only: collect_mqc_efp_pair_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_efp_pair", collect_mqc_efp_pair_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
