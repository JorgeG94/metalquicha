!! Occupation-restricted active spaces: classes, compatibility, and addressing
module test_mqc_ormas_space
   !! The stage an ORMAS CI is built on, and like the CAS strings before it,
   !! one that can be settled completely without any physics. The class lists
   !! are integer compositions, the counts are products of binomial
   !! coefficients, and the address is a bijection onto a contiguous range.
   !!
   !! Three kinds of check, and each catches what the others cannot.
   !!
   !! *Reference tables* are transcribed from a partition worked through by
   !! hand, entry by entry -- every class vector, the whole compatibility grid,
   !! every offset. They pin the enumeration order, which nothing else here
   !! does: an implementation that ordered the classes differently would still
   !! address a consistent bijection and still count correctly, and would be
   !! silently incompatible with every table built against it later.
   !!
   !! *Identities* -- the address is dense and injective, the grid is symmetric
   !! when the two spins are alike, tightening does not change what the space
   !! contains -- hold over every determinant of several partitions rather than
   !! being sampled, so they catch an implementation that is self-consistently
   !! wrong.
   !!
   !! *A brute-force oracle* counts determinants by generating every string of
   !! the full active space and testing the user's windows directly on the bit
   !! patterns. It shares no code at all with the module under test -- no
   !! classes, no compatibility, no offsets -- so it is the check that would
   !! survive the class machinery being wrong in a way that is internally
   !! consistent. It is affordable only because these spaces are tiny, which is
   !! exactly why the partitions here are small.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: int64
   use mqc_error, only: error_t
   use mqc_determinants, only: generate_strings
   use mqc_ormas_space, only: ormas_space_t, build_ormas_space, &
                              determinant_address, describes_a_cas, &
                              ormas_strings, ormas_string_address
   implicit none
   private

   public :: collect_mqc_ormas_space_tests

contains

   subroutine collect_mqc_ormas_space_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("two_subspace_tables", test_two_subspace), &
                  new_unittest("three_subspace_tables", test_three_subspace), &
                  new_unittest("counts_against_brute_force", test_brute_force), &
                  new_unittest("address_is_a_dense_bijection", test_bijection), &
                  new_unittest("compatibility_is_symmetric", test_symmetry), &
                  new_unittest("windows_are_tightened", test_tightening), &
                  new_unittest("open_windows_are_still_a_cas", test_cas_degeneracy), &
                  new_unittest("strings_match_the_worked_example", test_strings), &
                  new_unittest("string_addressing_round_trips", test_string_round_trip), &
                  new_unittest("dead_classes_pair_with_nothing", test_outside), &
                  new_unittest("unreachable_classes_are_dropped", test_pruning), &
                  new_unittest("the_worked_determinant_list", test_worked_determinants), &
                  new_unittest("refusals", test_refusals) &
                  ]
   end subroutine collect_mqc_ormas_space_tests

   subroutine test_two_subspace(error)
      !! Four orbitals cut in two, two electrons of each spin, at least two in
      !! the first subspace and at most two in the second
      !!
      !! Every table below was worked out by hand. The determinant count is 27,
      !! reached two ways: by summing `size(alpha class) * row length` over the
      !! classes, and by the closed form over compatible pairs.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      integer, parameter :: N_CLASSES = 3
      integer :: expected_class(2, N_CLASSES)
      integer :: ga, gb
      logical :: expected_compatible(3, 3)

      call build_ormas_space([1, 3], 4, 2, 2, [2, 0], [4, 2], space, err)
      call check(error,.not. err%has_error(), "building the partition failed")
      if (allocated(error)) return

      call check(error, space%n_alpha_classes, 3)
      if (allocated(error)) return
      call check(error, space%n_beta_classes, 3)
      if (allocated(error)) return

      expected_class = reshape([2, 0, 1, 1, 0, 2], [2, 3])
      do ga = 1, 3
         call check(error, all(space%alpha_class(:, ga) == expected_class(:, ga)), &
                    "alpha class vector")
         if (allocated(error)) return
      end do

      call check(error, all(space%alpha_class_size == [1_int64, 4_int64, 1_int64]), &
                 "strings per class")
      if (allocated(error)) return
      call check(error, all(space%alpha_offset == [0_int64, 1_int64, 5_int64, 6_int64]), &
                 "class string offsets")
      if (allocated(error)) return

      !  gb\ga   (2,0) (1,1) (0,2)
      !  (2,0)     T     T     T
      !  (1,1)     T     T     F
      !  (0,2)     T     F     F
      expected_compatible = reshape([.true., .true., .true., &
                                     .true., .true., .false., &
                                     .true., .false., .false.], [3, 3])
      do ga = 1, 3
         do gb = 1, 3
            call check(error, space%compatible(gb, ga) .eqv. expected_compatible(gb, ga), &
                       "compatibility grid")
            if (allocated(error)) return
         end do
      end do

      call check(error, all(space%row_length == [6_int64, 5_int64, 1_int64]), &
                 "determinants per alpha string")
      if (allocated(error)) return
      call check(error, all(space%class_base == [0_int64, 6_int64, 26_int64]), &
                 "determinants preceding each alpha class")
      if (allocated(error)) return

      call check(error, all(space%block_offset(:, 1) == [0_int64, 1_int64, 5_int64]), &
                 "beta block offsets within class (2,0)")
      if (allocated(error)) return
      call check(error, all(space%block_offset(1:2, 2) == [0_int64, 1_int64]), &
                 "beta block offsets within class (1,1)")
      if (allocated(error)) return
      call check(error, space%block_offset(1, 3) == 0_int64, &
                 "beta block offset within class (0,2)")
      if (allocated(error)) return

      call check(error, all(space%alpha_base == &
                            [0_int64, 6_int64, 11_int64, 16_int64, 21_int64, 26_int64]), &
                 "determinants preceding each alpha string")
      if (allocated(error)) return

      call check(error, space%n_determinants == 27_int64, "determinant count")
      if (allocated(error)) return
      call check(error,.not. describes_a_cas(space), "this partition restricts")
      if (allocated(error)) return

      call space%destroy()
   end subroutine test_two_subspace

   subroutine test_three_subspace(error)
      !! Six orbitals cut in three, three electrons of each spin
      !!
      !! The windows reduce to two live constraints -- at least two electrons
      !! in the first subspace, at most two in the third -- because each
      !! two-orbital subspace holds at most two of each spin anyway. 236
      !! determinants, checked here against the closed form and, in
      !! `test_brute_force`, against an enumeration that knows nothing of
      !! classes.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      integer, parameter :: N_CLASSES = 7
      integer :: expected_class(3, N_CLASSES)
      integer :: g
      integer(int64) :: closed_form
      integer :: ga, gb

      call build_ormas_space([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2], space, err)
      call check(error,.not. err%has_error(), "building the partition failed")
      if (allocated(error)) return

      call check(error, space%n_alpha_classes, 7)
      if (allocated(error)) return

      expected_class = reshape([2, 1, 0, &
                                2, 0, 1, &
                                1, 2, 0, &
                                1, 1, 1, &
                                1, 0, 2, &
                                0, 2, 1, &
                                0, 1, 2], [3, 7])
      do g = 1, 7
         call check(error, all(space%alpha_class(:, g) == expected_class(:, g)), &
                    "alpha class vector")
         if (allocated(error)) return
      end do

      call check(error, all(space%alpha_class_size == &
                            [2_int64, 2_int64, 2_int64, 8_int64, 2_int64, 2_int64, 2_int64]), &
                 "strings per class")
      if (allocated(error)) return

      ! Every string of this space is reachable -- no subspace can overflow with
      ! only three electrons of a spin and two orbitals per subspace.
      call check(error, space%alpha_offset(8) == 20_int64, "C(6,3) strings in total")
      if (allocated(error)) return

      call check(error, count(space%compatible), 26, "compatible class pairs")
      if (allocated(error)) return

      ! The boundary is inclusive: class 2 with class 4 puts exactly two
      ! electrons in the third subspace, which is the maximum and allowed.
      call check(error, space%compatible(4, 2), "an exactly-at-the-maximum pair")
      if (allocated(error)) return
      ! Three in the third subspace is one too many, though the first subspace
      ! is satisfied.
      call check(error,.not. space%compatible(2, 5), "one over the maximum")
      if (allocated(error)) return
      ! One electron in the first subspace is one too few.
      call check(error,.not. space%compatible(6, 3), "one under the minimum")
      if (allocated(error)) return

      closed_form = 0_int64
      do ga = 1, space%n_alpha_classes
         do gb = 1, space%n_beta_classes
            if (.not. space%compatible(gb, ga)) cycle
            closed_form = closed_form + space%alpha_class_size(ga)*space%beta_class_size(gb)
         end do
      end do
      call check(error, closed_form == 236_int64, "closed-form determinant count")
      if (allocated(error)) return
      call check(error, space%n_determinants == 236_int64, "determinant count")
      if (allocated(error)) return

      call space%destroy()
   end subroutine test_three_subspace

   pure function subspace_occupancy(string, first_orbital, n_orbitals) result(occupancy)
      !! How many electrons a string puts in one subspace
      integer(int64), intent(in) :: string
      integer, intent(in) :: first_orbital, n_orbitals
      integer :: occupancy

      integer :: orbital

      occupancy = 0
      do orbital = first_orbital, first_orbital + n_orbitals - 1
         if (btest(string, orbital - 1)) occupancy = occupancy + 1
      end do
   end function subspace_occupancy

   function brute_force_count(first_orbital, n_active, n_alpha, n_beta, &
                              min_electrons, max_electrons) result(count)
      !! Determinants, by testing every string pair against the windows
      !!
      !! Deliberately shares nothing with the module under test: it generates
      !! the full active space with `mqc_determinants` and counts bits. If this
      !! and `n_determinants` agree, they are not agreeing because of a common
      !! mistake.
      integer, intent(in) :: first_orbital(:), min_electrons(:), max_electrons(:)
      integer, intent(in) :: n_active, n_alpha, n_beta
      integer(int64) :: count

      integer(int64), allocatable :: alpha(:), beta(:)
      type(error_t) :: err
      integer :: n_subspaces, k, width, occupancy
      integer(int64) :: ia, ib
      logical :: allowed

      n_subspaces = size(first_orbital)
      call generate_strings(n_active, n_alpha, alpha, err)
      call generate_strings(n_active, n_beta, beta, err)

      count = 0_int64
      do ia = 1, size(alpha, kind=int64)
         do ib = 1, size(beta, kind=int64)
            allowed = .true.
            do k = 1, n_subspaces
               if (k < n_subspaces) then
                  width = first_orbital(k + 1) - first_orbital(k)
               else
                  width = n_active - first_orbital(k) + 1
               end if
               occupancy = subspace_occupancy(alpha(ia), first_orbital(k), width) &
                           + subspace_occupancy(beta(ib), first_orbital(k), width)
               if (occupancy < min_electrons(k)) allowed = .false.
               if (occupancy > max_electrons(k)) allowed = .false.
            end do
            if (allowed) count = count + 1_int64
         end do
      end do
   end function brute_force_count

   subroutine test_brute_force(error)
      !! Every count, twice over: against an enumeration that knows nothing of
      !! classes, and against a number computed outside this language
      !!
      !! The literals come from `tools/ormas_reference/ormas_reference.py`,
      !! which enumerates determinants by brute force over spin-orbital sets
      !! and has itself been checked against `pyscf.fci` on unrestricted
      !! windows. Two independent counts agreeing is worth more than either --
      !! the local oracle shares this file's understanding of what the windows
      !! mean, and the external one does not.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err

      call build_ormas_space([1, 3], 4, 2, 2, [2, 0], [4, 2], space, err)
      call check(error, space%n_determinants == &
                 brute_force_count([1, 3], 4, 2, 2, [2, 0], [4, 2]), "two subspaces")
      if (allocated(error)) return
      call check(error, space%n_determinants == 27_int64, "against the python reference")
      if (allocated(error)) return
      call space%destroy()

      call build_ormas_space([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2], space, err)
      call check(error, space%n_determinants == &
                 brute_force_count([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2]), &
                 "three subspaces")
      if (allocated(error)) return
      call check(error, space%n_determinants == 236_int64, "against the python reference")
      if (allocated(error)) return
      call space%destroy()

      ! A truncated CI: everything in the first subspace but at most two
      ! electrons promoted into the second.
      call build_ormas_space([1, 4], 7, 3, 3, [4, 0], [6, 2], space, err)
      call check(error, space%n_determinants == &
                 brute_force_count([1, 4], 7, 3, 3, [4, 0], [6, 2]), "a truncated CI")
      if (allocated(error)) return
      call check(error, space%n_determinants == 205_int64, "against the python reference")
      if (allocated(error)) return
      call space%destroy()

      ! Unequal spins, so the two class lists differ and the grid is not
      ! symmetric.
      call build_ormas_space([1, 3, 5], 6, 3, 1, [1, 0, 0], [4, 4, 2], space, err)
      call check(error, space%n_determinants == &
                 brute_force_count([1, 3, 5], 6, 3, 1, [1, 0, 0], [4, 4, 2]), &
                 "more alpha than beta")
      if (allocated(error)) return
      call check(error, space%n_determinants == 100_int64, "against the python reference")
      if (allocated(error)) return
      call space%destroy()

      ! Windows the tightening pass rewrites -- see `test_tightening`.
      call build_ormas_space([1, 4], 6, 2, 2, [0, 3], [6, 6], space, err)
      call check(error, space%n_determinants == &
                 brute_force_count([1, 4], 6, 2, 2, [0, 3], [6, 6]), &
                 "windows that get tightened")
      if (allocated(error)) return
      call check(error, space%n_determinants == 63_int64, "against the python reference")
      if (allocated(error)) return
      call space%destroy()
   end subroutine test_brute_force

   subroutine test_bijection(error)
      !! Every compatible string pair gets its own address, and together they
      !! cover `1 .. n_determinants` with no gaps
      !!
      !! This is the property the whole layout exists to have, and the one that
      !! fails quietly: an off-by-one in a single offset leaves the addressing
      !! self-consistent for most pairs and collides on a few.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      logical, allocatable :: seen(:)
      integer :: a, b, ga, gb
      integer(int64) :: address

      call build_ormas_space([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2], space, err)
      call check(error,.not. err%has_error(), "building the partition failed")
      if (allocated(error)) return

      allocate (seen(space%n_determinants))
      seen = .false.

      do a = 1, int(space%alpha_offset(space%n_alpha_classes + 1))
         do b = 1, int(space%beta_offset(space%n_beta_classes + 1))
            ga = space%alpha_string_class(a)
            gb = space%beta_string_class(b)
            if (.not. space%compatible(gb, ga)) cycle

            address = determinant_address(space, a, b)
            call check(error, address >= 1_int64 .and. address <= space%n_determinants, &
                       "address is inside the vector")
            if (allocated(error)) return
            call check(error,.not. seen(address), "two determinants share an address")
            if (allocated(error)) return
            seen(address) = .true.
         end do
      end do

      call check(error, all(seen), "the addresses leave a gap")
      if (allocated(error)) return

      deallocate (seen)
      call space%destroy()
   end subroutine test_bijection

   subroutine test_symmetry(error)
      !! With as many alpha electrons as beta, the grid is its own transpose
      !!
      !! The windows constrain the two spins only through their sum, so a pair
      !! of classes is allowed regardless of which spin is which. Nothing
      !! downstream would work without this -- it is what lets a determinant
      !! and its spin-transpose be addressed from the same tables -- and it is
      !! cheap enough to assert directly.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      integer :: g

      call build_ormas_space([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2], space, err)
      call check(error,.not. err%has_error(), "building the partition failed")
      if (allocated(error)) return

      call check(error, space%n_alpha_classes, space%n_beta_classes)
      if (allocated(error)) return

      do g = 1, space%n_alpha_classes
         call check(error, all(space%alpha_class(:, g) == space%beta_class(:, g)), &
                    "the two spins enumerate the same classes")
         if (allocated(error)) return
      end do

      call check(error, all(space%compatible .eqv. transpose(space%compatible)), &
                 "the compatibility grid is not symmetric")
      if (allocated(error)) return

      call space%destroy()
   end subroutine test_symmetry

   subroutine test_tightening(error)
      !! Windows are reported as what they reduce to, not as what was asked
      !!
      !! Six orbitals in two subspaces with four electrons: asking for up to
      !! six in each subspace is not a real allowance, because the second
      !! subspace must hold at least three and so the first can never exceed
      !! one. The pass that notices this runs before anything else reads the
      !! windows, which is why the numbers stored here differ from the
      !! arguments.
      !!
      !! And what it must *not* do is change the space. Tightening removes
      !! promises the partition could not have kept, so the determinant count
      !! is the same either way -- asserted against brute force, which is given
      !! the original windows.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err

      call build_ormas_space([1, 4], 6, 2, 2, [0, 3], [6, 6], space, err)
      call check(error,.not. err%has_error(), "building the partition failed")
      if (allocated(error)) return

      call check(error, all(space%max_electrons == [1, 4]), &
                 "the maxima were not tightened")
      if (allocated(error)) return
      call check(error, all(space%min_electrons == [0, 3]), &
                 "the minima should not have moved")
      if (allocated(error)) return

      call check(error, space%n_determinants == &
                 brute_force_count([1, 4], 6, 2, 2, [0, 3], [6, 6]), &
                 "tightening changed which determinants exist")
      if (allocated(error)) return

      call space%destroy()
   end subroutine test_tightening

   subroutine test_cas_degeneracy(error)
      !! A partition that restricts nothing is a complete active space, however
      !! many subspaces are drawn on it
      !!
      !! The second case here has three occupation classes per spin and a full
      !! compatibility grid, so it exercises all the machinery and must still
      !! come out as the plain CAS count. Recognising it matters later: in a
      !! genuine CAS, rotating one active orbital into another changes nothing,
      !! and an optimiser that believes that of a restricted space converges to
      !! the wrong answer.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err

      call build_ormas_space([1], 4, 2, 2, [4], [4], space, err)
      call check(error,.not. err%has_error(), "building the partition failed")
      if (allocated(error)) return
      call check(error, space%n_determinants == 36_int64, "CAS(4,4) determinant count")
      if (allocated(error)) return
      call check(error, describes_a_cas(space), "one open subspace is a CAS")
      if (allocated(error)) return
      call space%destroy()

      call build_ormas_space([1, 3], 4, 2, 2, [0, 0], [4, 4], space, err)
      call check(error,.not. err%has_error(), "building the partition failed")
      if (allocated(error)) return
      call check(error, space%n_alpha_classes, 3)
      if (allocated(error)) return
      call check(error, all(space%compatible), "open windows allow every pair")
      if (allocated(error)) return
      call check(error, space%n_determinants == 36_int64, &
                 "two subspaces with open windows are still CAS(4,4)")
      if (allocated(error)) return
      call check(error, describes_a_cas(space), &
                 "a partition that restricts nothing is a CAS")
      if (allocated(error)) return
      call space%destroy()
   end subroutine test_cas_degeneracy

   subroutine test_strings(error)
      !! The six strings of the two-subspace partition, as bit patterns
      !!
      !! Class (2,0) puts both electrons in orbitals 1-2; class (1,1) takes one
      !! from each subspace with the second varying fastest, giving (1,3),
      !! (1,4), (2,3), (2,4); class (0,2) puts both in orbitals 3-4. Bit `k` is
      !! orbital `k + 1`, so those are 3, 5, 9, 6, 10, 12 -- an order that is
      !! *not* ascending numerically, because the class grouping outranks it.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      integer(int64), allocatable :: alpha(:), beta(:)

      call build_ormas_space([1, 3], 4, 2, 2, [2, 0], [4, 2], space, err)
      call ormas_strings(space, alpha, beta, err)
      call check(error,.not. err%has_error(), "generating the strings failed")
      if (allocated(error)) return

      call check(error, size(alpha), 6, "how many alpha strings")
      if (allocated(error)) return
      call check(error, all(alpha == [3_int64, 5_int64, 9_int64, 6_int64, 10_int64, 12_int64]), &
                 "alpha strings")
      if (allocated(error)) return
      call check(error, all(beta == alpha), "the two spins have the same strings here")
      if (allocated(error)) return

      call space%destroy()
   end subroutine test_strings

   subroutine round_trips(space, error)
      !! Every string of a space addresses back to where it was found
      type(ormas_space_t), intent(in) :: space
      type(error_type), allocatable, intent(inout) :: error

      integer(int64), allocatable :: alpha(:), beta(:)
      type(error_t) :: err
      integer(int64) :: i

      call ormas_strings(space, alpha, beta, err)
      if (err%has_error()) return

      do i = 1, size(alpha, kind=int64)
         call check(error, ormas_string_address(space, alpha(i), .true.) == i, &
                    "an alpha string does not address back to itself")
         if (allocated(error)) return
         call check(error, int(popcnt(alpha(i))), space%n_alpha, "alpha electron count")
         if (allocated(error)) return

         ! Generation walks the classes in order, so a string's position must
         ! fall inside the range its class was allotted.
         call check(error, i > space%alpha_offset(space%alpha_string_class(int(i))) .and. &
                    i <= space%alpha_offset(space%alpha_string_class(int(i)) + 1), &
                    "a string sits outside its own class's range")
         if (allocated(error)) return
      end do

      do i = 1, size(beta, kind=int64)
         call check(error, ormas_string_address(space, beta(i), .false.) == i, &
                    "a beta string does not address back to itself")
         if (allocated(error)) return
         call check(error, int(popcnt(beta(i))), space%n_beta, "beta electron count")
         if (allocated(error)) return
      end do
   end subroutine round_trips

   subroutine test_string_round_trip(error)
      !! Generation and addressing are inverse, over every string of several
      !! partitions
      !!
      !! The two are written independently -- one builds a string from a class
      !! by taking a product over subspaces, the other takes a string apart and
      !! reconstructs a mixed-radix index -- so agreeing over the whole list is
      !! a real constraint on both, not a tautology.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err

      call build_ormas_space([1, 3], 4, 2, 2, [2, 0], [4, 2], space, err)
      call round_trips(space, error)
      if (allocated(error)) return
      call space%destroy()

      call build_ormas_space([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2], space, err)
      call round_trips(space, error)
      if (allocated(error)) return
      call space%destroy()

      ! Subspaces of unequal width, so the mixed-radix weights differ.
      call build_ormas_space([1, 4], 7, 3, 3, [4, 0], [6, 2], space, err)
      call round_trips(space, error)
      if (allocated(error)) return
      call space%destroy()

      ! Unequal spins, so the two lists are genuinely different.
      call build_ormas_space([1, 3, 5], 6, 3, 1, [1, 0, 0], [4, 4, 2], space, err)
      call round_trips(space, error)
      if (allocated(error)) return
      call space%destroy()
   end subroutine test_string_round_trip

   subroutine test_outside(error)
      !! Every string exists; some of them pair with nothing
      !!
      !! Six orbitals in two subspaces, the second needing at least three
      !! electrons. No determinant may put both of one spin's electrons in the
      !! first subspace -- but the *string* that does so still exists and still
      !! has an address, because a string is not a determinant and an excitation
      !! will need to name it. What it does not have is a partner: its
      !! occupation class is compatible with no class of the other spin.
      !!
      !! Keeping such strings is what lets one list serve both the space and
      !! anything that steps outside it, which is what an excitation does.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      integer(int64), allocatable :: alpha(:), beta(:)
      integer(int64) :: at
      integer :: dead

      call build_ormas_space([1, 4], 6, 2, 2, [0, 3], [6, 6], space, err)
      call ormas_strings(space, alpha, beta, err)
      call check(error,.not. err%has_error(), "generating the strings failed")
      if (allocated(error)) return

      call check(error, size(alpha), 15, "every string of the full space is kept")
      if (allocated(error)) return

      ! Orbitals 1 and 2, both in the first subspace: a real string, with an
      ! address, whose class pairs with nothing.
      at = ormas_string_address(space, 3_int64, .true.)
      call check(error, at > 0_int64, "a string was denied an address")
      if (allocated(error)) return
      dead = space%alpha_string_class(int(at))
      call check(error,.not. any(space%compatible(:, dead)), &
                 "a class the windows exclude found a partner")
      if (allocated(error)) return

      ! Orbitals 1 and 4 straddle the subspaces, so this class does pair.
      at = ormas_string_address(space, 9_int64, .true.)
      call check(error, any(space%compatible(:, space%alpha_string_class(int(at)))), &
                 "a usable class was left with no partner")
      if (allocated(error)) return

      ! The grid decides the determinant count, not which strings exist, so
      ! keeping the dead ones changes nothing.
      call check(error, space%n_determinants == 63_int64, "determinant count")
      if (allocated(error)) return

      call space%destroy()
   end subroutine test_outside

   subroutine test_worked_determinants(error)
      !! Rows of the determinant list worked out by hand, by their strings
      !!
      !! Now that a string index means something concrete, the addressing can
      !! be checked against determinants named the way a person would name
      !! them -- which orbitals each spin occupies -- rather than against
      !! indices that came out of the same machinery being tested.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      integer(int64), allocatable :: alpha(:), beta(:)
      integer :: a, b

      call build_ormas_space([1, 3], 4, 2, 2, [2, 0], [4, 2], space, err)
      call ormas_strings(space, alpha, beta, err)
      call check(error,.not. err%has_error(), "building the partition failed")
      if (allocated(error)) return

      ! alpha (1,3) with beta (2,4) is the eleventh determinant
      a = int(ormas_string_address(space, 5_int64, .true.))
      b = int(ormas_string_address(space, 10_int64, .false.))
      call check(error, determinant_address(space, a, b) == 11_int64, "det 11")
      if (allocated(error)) return

      ! alpha (1,4) with beta (1,3) is the thirteenth
      a = int(ormas_string_address(space, 9_int64, .true.))
      b = int(ormas_string_address(space, 5_int64, .false.))
      call check(error, determinant_address(space, a, b) == 13_int64, "det 13")
      if (allocated(error)) return

      ! alpha (2,3) with beta (2,3) is the twentieth
      a = int(ormas_string_address(space, 6_int64, .true.))
      call check(error, determinant_address(space, a, a) == 20_int64, "det 20")
      if (allocated(error)) return

      ! alpha (2,4) with beta (1,2) is the twenty-second
      a = int(ormas_string_address(space, 10_int64, .true.))
      b = int(ormas_string_address(space, 3_int64, .false.))
      call check(error, determinant_address(space, a, b) == 22_int64, "det 22")
      if (allocated(error)) return

      ! alpha (1,2) with beta (3,4) is the sixth -- the last of the first row,
      ! reached only because class (2,0) is compatible with every beta class
      a = int(ormas_string_address(space, 3_int64, .true.))
      b = int(ormas_string_address(space, 12_int64, .false.))
      call check(error, determinant_address(space, a, b) == 6_int64, "det 6")
      if (allocated(error)) return

      call space%destroy()
   end subroutine test_worked_determinants

   subroutine test_pruning(error)
      !! Classes further than one excitation from the space are not kept
      !!
      !! Twelve orbitals in three subspaces, the third barred from holding any
      !! electron at all. Of the ten ways to spread three electrons of one spin,
      !! only three appear in a determinant; four more sit one excitation away
      !! and have to be kept, because the sigma build's intermediate lands on
      !! them and needs to name where it landed. The remaining three -- two or
      !! three electrons in the barred subspace -- are two excitations out and
      !! can be reached by nothing.
      !!
      !! Small here, and the reason to do it at all is not. For singles and
      !! doubles from seven orbitals into twenty-nine this keeps 136,620 alpha
      !! strings out of 8,347,680, and every table indexed by a string shrinks
      !! by the same factor.
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      integer(int64), allocatable :: alpha(:), beta(:)
      integer :: g

      call build_ormas_space([1, 5, 9], 12, 3, 3, [4, 0, 0], [6, 2, 0], space, err)
      call check(error,.not. err%has_error(), "building the partition failed")
      if (allocated(error)) return

      call check(error, space%n_alpha_classes, 7, &
                 "seven of the ten classes are within one excitation of the space")
      if (allocated(error)) return

      ! Nothing kept may put two or more electrons in the barred subspace: that
      ! is two excitations from anything the space contains.
      do g = 1, space%n_alpha_classes
         call check(error, space%alpha_class(3, g) <= 1, &
                    "a class two excitations out was kept")
         if (allocated(error)) return
      end do

      call ormas_strings(space, alpha, beta, err)
      call check(error, size(alpha), 168, "and the strings shrink with the classes")
      if (allocated(error)) return

      ! Dropping them changes nothing about the space itself.
      call check(error, space%n_determinants == &
                 brute_force_count([1, 5, 9], 12, 3, 3, [4, 0, 0], [6, 2, 0]), &
                 "pruning changed which determinants exist")
      if (allocated(error)) return

      call space%destroy()
   end subroutine test_pruning

   subroutine test_refusals(error)
      !! The partitions that describe nothing, and are said so
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err

      ! Minima that cannot all be met at once.
      call build_ormas_space([1, 3], 4, 1, 1, [2, 2], [4, 4], space, err)
      call check(error, err%has_error(), "impossible minima were accepted")
      if (allocated(error)) return
      call space%destroy()

      ! Maxima that do not leave room for the electrons.
      call err%clear()
      call build_ormas_space([1, 3], 4, 2, 2, [0, 0], [1, 1], space, err)
      call check(error, err%has_error(), "maxima too small were accepted")
      if (allocated(error)) return
      call space%destroy()

      ! A subspace asked to hold more electrons than its orbitals can.
      call err%clear()
      call build_ormas_space([1, 3], 4, 2, 2, [0, 0], [6, 4], space, err)
      call check(error, err%has_error(), "an overfull subspace was accepted")
      if (allocated(error)) return
      call space%destroy()

      ! A minimum above its own maximum.
      call err%clear()
      call build_ormas_space([1, 3], 4, 2, 2, [3, 0], [2, 4], space, err)
      call check(error, err%has_error(), "an inverted window was accepted")
      if (allocated(error)) return
      call space%destroy()

      ! Windows given for the wrong number of subspaces.
      call err%clear()
      call build_ormas_space([1, 3], 4, 2, 2, [0], [4, 4], space, err)
      call check(error, err%has_error(), "a short window list was accepted")
      if (allocated(error)) return
      call space%destroy()

      ! Subspaces that do not start at the first active orbital.
      call err%clear()
      call build_ormas_space([2, 3], 4, 2, 2, [0, 0], [4, 4], space, err)
      call check(error, err%has_error(), "a partition with a gap was accepted")
      if (allocated(error)) return
      call space%destroy()
   end subroutine test_refusals

end module test_mqc_ormas_space

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_ormas_space, only: collect_mqc_ormas_space_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_ormas_space", collect_mqc_ormas_space_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
