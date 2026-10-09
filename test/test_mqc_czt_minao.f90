module test_mqc_czt_minao
   !! Pins the minao initial guess against PySCF's `init_guess_by_minao`.
   !!
   !! The reference numbers come from `tools/minao/minao_reference.py`, which
   !! reads the orbital basis from this repository's own JSON so the two codes
   !! share both the basis data and the AO order. What is checked:
   !!
   !!   * the density matches PySCF's, on its diagonal and first row, for water
   !!     in cc-pVDZ and methanol in def2-SVP;
   !!   * projected into the minimal basis itself the density is exactly the
   !!     diagonal occupation matrix, carrying the neutral electron count;
   !!   * the occupation table reproduces PySCF's `frac_occ`;
   !!   * the projection survives a linear-dependence threshold that drops
   !!     directions, and refuses an element outside H-Kr rather than guessing.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_minao, only: build_minao_guess, minao_occupation
   use mqc_czt_rhf, only: SCF_GUESS_MINAO
   use mqc_czt_atomic_guess, only: build_atomic_guess, parse_guess_name
   implicit none
   private
   public :: collect_mqc_czt_minao_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: PYSCF_TOL = 1.0e-9_dp
      !! The reference is printed to 1e-12; what is left is basis-data
      !! normalisation at the precision both codes evaluate it.

contains

   subroutine collect_mqc_czt_minao_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_matches_pyscf", test_water_pyscf), &
                  new_unittest("methanol_matches_pyscf", test_methanol_pyscf), &
                  new_unittest("self_projection_is_the_occupation", test_self_projection), &
                  new_unittest("occupations_follow_frac_occ", test_occupations), &
                  new_unittest("spins_are_halves", test_spin_halves), &
                  new_unittest("pruned_projection_stays_finite", test_pruned), &
                  new_unittest("heavy_element_is_refused", test_heavy_refused) &
                  ]
   end subroutine collect_mqc_czt_minao_tests

   subroutine water(mol, err, basis)
      !! The water of test_mqc_czt_guess, and of minao_reference.py
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
   end subroutine water

   subroutine methanol(mol, err, basis)
      !! Methanol, as in minao_reference.py
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis
      real(dp) :: c(3, 6)

      c = ANG*reshape([-0.0467_dp, 0.6635_dp, 0.0000_dp, &
                       -0.0467_dp, -0.7581_dp, 0.0000_dp, &
                       -1.0868_dp, 0.9790_dp, 0.0000_dp, &
                       0.4377_dp, 1.0739_dp, 0.8920_dp, &
                       0.4377_dp, 1.0739_dp, -0.8920_dp, &
                       0.8752_dp, -1.0659_dp, 0.0000_dp], [3, 6])
      call build_czt_molecule([6, 8, 1, 1, 1, 1], ["C ", "O ", "H ", "H ", "H ", "H "], &
                              c, basis, mol, err)
   end subroutine methanol

   function electrons(mol, density) result(n)
      !! Tr(D S)
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp) :: n
      real(dp), allocatable :: s(:, :)

      call mol%overlap(s)
      n = sum(density*s)
   end function electrons

   subroutine compare(error, label, got, want)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: got(:), want(:)
      character(len=160) :: msg
      integer :: i

      call check(error, size(got) == size(want), label//": length differs from the reference")
      if (allocated(error)) return
      do i = 1, size(want)
         if (abs(got(i) - want(i)) > PYSCF_TOL) then
            write (msg, "(a,a,i0,a,es12.4,a,es12.4)") label, " element ", i, ": ", got(i), &
               " where PySCF has ", want(i)
            call check(error, .false., trim(msg))
            return
         end if
      end do
   end subroutine compare

   subroutine test_water_pyscf(error)
      !! Water/cc-pVDZ, element by element against PySCF 2.14
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d(:, :)
      integer :: i
      real(dp), parameter :: DIAG(24) = [ &
                             1.999121660913_dp, 2.045448624773_dp, 0.013438847942_dp, 1.303452732810_dp, &
                             1.289113361100_dp, 1.308442373123_dp, 0.001835401070_dp, 0.000562967490_dp, &
                             0.002113412362_dp, 0.000000974913_dp, 0.000000574546_dp, 0.000048243751_dp, &
                             0.000005348451_dp, 0.000055266613_dp, 1.423353976779_dp, 0.047762202694_dp, &
                             0.000021553466_dp, 0.000010183876_dp, 0.000723795081_dp, 1.423619277045_dp, &
                             0.047762763257_dp, 0.000661472313_dp, 0.000010267070_dp, 0.000085584373_dp]
      real(dp), parameter :: ROW1(24) = [ &
                             1.999121660913_dp, 0.005538835341_dp, 0.001449798841_dp, 0.008146900334_dp, &
                             0.000000000000_dp, 0.006304338256_dp, -0.000826423747_dp, 0.000000000000_dp, &
                             -0.000639782220_dp, 0.000000000000_dp, 0.000000000000_dp, 0.000078876187_dp, &
                             -0.000101434907_dp, 0.000086307309_dp, 0.002565465528_dp, 0.001450217254_dp, &
                             0.000089044562_dp, 0.000000000000_dp, 0.000527792379_dp, 0.002543456623_dp, &
                             0.001458496631_dp, 0.000536727575_dp, 0.000000000000_dp, -0.000046660237_dp]

      call water(mol, err, "cc-pvdz")
      call build_minao_guess(mol, d, err)
      call check(error,.not. err%has_error(), "minao must build: "//err%get_full_trace())
      if (allocated(error)) return

      call compare(error, "water diagonal", [(d(i, i), i=1, size(d, 1))], DIAG)
      if (allocated(error)) return
      call compare(error, "water row 1", d(1, :), ROW1)
      if (allocated(error)) return
      call check(error, abs(electrons(mol, d) - 9.987286549344_dp) < PYSCF_TOL, &
                 "water Tr(D S) differs from PySCF's")
      if (allocated(error)) return
      call check(error, abs(norm2(d) - 4.201891741805_dp) < PYSCF_TOL, &
                 "water |D| differs from PySCF's")
      call mol%destroy()
   end subroutine test_water_pyscf

   subroutine test_methanol_pyscf(error)
      !! Methanol/def2-SVP, with d on the heavy atoms
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d(:, :)
      integer :: i
      real(dp), parameter :: DIAG(48) = [ &
                             2.160678948978_dp, 1.111139916932_dp, 0.251772453103_dp, 0.338769018709_dp, &
                             0.339745757473_dp, 0.338525473360_dp, 0.109958974175_dp, 0.112995570274_dp, &
                             0.110919558083_dp, 0.000037514272_dp, 0.000058452435_dp, 0.000056372239_dp, &
                             0.000033645821_dp, 0.000086273402_dp, 2.163648193602_dp, 0.909642497929_dp, &
                             0.330124861589_dp, 0.616769407147_dp, 0.617991732592_dp, 0.612571492633_dp, &
                             0.249047720214_dp, 0.256776296759_dp, 0.279999135251_dp, 0.000015733695_dp, &
                             0.000004784200_dp, 0.000024131267_dp, 0.000000244869_dp, 0.000051838252_dp, &
                             0.481321258989_dp, 0.188040793325_dp, 0.000212169850_dp, 0.000028401705_dp, &
                             0.000016221192_dp, 0.482323289337_dp, 0.186919726538_dp, 0.000061503654_dp, &
                             0.000039334583_dp, 0.000161243605_dp, 0.482323289337_dp, 0.186919726538_dp, &
                             0.000061503654_dp, 0.000039334583_dp, 0.000161243605_dp, 0.486955261580_dp, &
                             0.184436557515_dp, 0.000480694100_dp, 0.000065684192_dp, 0.000001351504_dp]
      real(dp), parameter :: ROW1(48) = [ &
                             2.160678948978_dp, -0.398509614273_dp, -0.248673889884_dp, -0.000081889193_dp, &
                             0.000218647144_dp, -0.000000000000_dp, 0.000437770903_dp, -0.002684533123_dp, &
                             0.000000000000_dp, 0.000030070597_dp, 0.000000000000_dp, 0.000022859814_dp, &
                             0.000000000000_dp, -0.000032923632_dp, -0.000268699102_dp, -0.000598666629_dp, &
                             -0.000422577698_dp, -0.000106116500_dp, 0.000196687643_dp, 0.000000000000_dp, &
                             -0.000262680600_dp, -0.002713162925_dp, -0.000000000000_dp, -0.000065242914_dp, &
                             -0.000000000000_dp, 0.000164080626_dp, -0.000000000000_dp, 0.000087892556_dp, &
                             -0.001259931310_dp, 0.004531022146_dp, -0.000441194251_dp, 0.000249917254_dp, &
                             -0.000000000000_dp, -0.001453475010_dp, 0.004490228807_dp, 0.000196088788_dp, &
                             0.000297012920_dp, 0.000414943365_dp, -0.001453475010_dp, 0.004490228807_dp, &
                             0.000196088788_dp, 0.000297012920_dp, -0.000414943365_dp, 0.000539385885_dp, &
                             -0.001860894502_dp, 0.000568424070_dp, 0.000193570161_dp, 0.000000000000_dp]

      call methanol(mol, err, "def2-svp")
      call build_minao_guess(mol, d, err)
      call check(error,.not. err%has_error(), "minao must build: "//err%get_full_trace())
      if (allocated(error)) return

      call compare(error, "methanol diagonal", [(d(i, i), i=1, size(d, 1))], DIAG)
      if (allocated(error)) return
      call compare(error, "methanol row 1", d(1, :), ROW1)
      if (allocated(error)) return
      call check(error, abs(electrons(mol, d) - 17.985552263370_dp) < PYSCF_TOL, &
                 "methanol Tr(D S) differs from PySCF's")
      if (allocated(error)) return
      call check(error, abs(norm2(d) - 4.280559076835_dp) < PYSCF_TOL, &
                 "methanol |D| differs from PySCF's")
      call mol%destroy()
   end subroutine test_methanol_pyscf

   subroutine test_self_projection(error)
      !! In its own basis the projection is the identity
      !!
      !! P = S^-1 S is one, so the guess must come back as the diagonal
      !! occupation matrix, and its trace against S as exactly the ten
      !! electrons of the neutral atoms.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d(:, :)
      real(dp) :: off
      integer :: i, j

      call water(mol, err, "minao/minao")
      call check(error,.not. err%has_error(), "water must build in the minimal basis: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call build_minao_guess(mol, d, err)
      call check(error,.not. err%has_error(), "minao must build: "//err%get_full_trace())
      if (allocated(error)) return

      ! O: 1s 2s 2p(x3) at 2, 2, 4/3 each; H: 1s at 1
      call check(error, size(d, 1) == 7, "water's minimal basis is seven functions")
      if (allocated(error)) return
      call compare(error, "self-projected diagonal", [(d(i, i), i=1, 7)], &
                   [2.0_dp, 2.0_dp, 4.0_dp/3.0_dp, 4.0_dp/3.0_dp, 4.0_dp/3.0_dp, 1.0_dp, 1.0_dp])
      if (allocated(error)) return
      off = 0.0_dp
      do j = 1, 7
         do i = 1, 7
            if (i /= j) off = max(off, abs(d(i, j)))
         end do
      end do
      call check(error, off < 1.0e-10_dp, "the self-projected density must be diagonal")
      if (allocated(error)) return
      call check(error, abs(electrons(mol, d) - 10.0_dp) < 1.0e-10_dp, &
                 "the self-projected density must carry ten electrons")
      call mol%destroy()
   end subroutine test_self_projection

   subroutine test_occupations(error)
      !! A few rows of PySCF's frac_occ, including the table's oddities
      type(error_type), allocatable, intent(out) :: error
      integer :: nd
      real(dp) :: f

      call minao_occupation(8, 1, nd, f)      ! O 2p4: the open shell, nothing doubly occupied
      call check(error, nd == 0 .and. abs(f - 4.0_dp/3.0_dp) < 1.0e-14_dp, "O p")
      if (allocated(error)) return
      call minao_occupation(10, 1, nd, f)     ! Ne: 2p closed, no open column
      call check(error, nd == 1 .and. f == 0.0_dp, "Ne p")
      if (allocated(error)) return
      call minao_occupation(1, 0, nd, f)      ! H 1s1
      call check(error, nd == 0 .and. abs(f - 1.0_dp) < 1.0e-14_dp, "H s")
      if (allocated(error)) return
      call minao_occupation(21, 1, nd, f)     ! Sc: PySCF puts 13 electrons in p
      call check(error, nd == 2 .and. abs(f - 1.0_dp/3.0_dp) < 1.0e-14_dp, "Sc p")
      if (allocated(error)) return
      call minao_occupation(26, 2, nd, f)     ! Fe 3d8 in this table
      call check(error, nd == 0 .and. abs(f - 1.6_dp) < 1.0e-14_dp, "Fe d")
      if (allocated(error)) return
      call minao_occupation(37, 0, nd, f)     ! Rb: outside the data
      call check(error, nd == 0 .and. f == 0.0_dp, "Rb is not covered")
   end subroutine test_occupations

   subroutine test_spin_halves(error)
      !! Through the atomic-guess entry point each spin gets half
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d(:, :), d_a(:, :), d_b(:, :)
      integer :: kind

      call parse_guess_name("minao", kind, err)
      call check(error,.not. err%has_error() .and. kind == SCF_GUESS_MINAO, &
                 "'minao' must parse to SCF_GUESS_MINAO")
      if (allocated(error)) return

      call water(mol, err, "cc-pvdz")
      call build_minao_guess(mol, d, err)
      call build_atomic_guess(mol, SCF_GUESS_MINAO, d_a, d_b, err)
      call check(error,.not. err%has_error(), "minao must build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(d_a - 0.5_dp*d)) < 1.0e-14_dp .and. &
                 maxval(abs(d_b - d_a)) == 0.0_dp, "alpha and beta must each be half")
      call mol%destroy()
   end subroutine test_spin_halves

   subroutine test_pruned(error)
      !! A threshold that drops directions still gives a finite guess
      !!
      !! aug-cc-pVDZ water has nothing near 1e-3 that matters to the minimal
      !! functions, so dropping those eigenvectors moves the electron count
      !! a little and must not blow anything up.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d(:, :)
      real(dp) :: n

      call water(mol, err, "aug-cc-pvdz")
      call build_minao_guess(mol, d, err, threshold=1.0e-3_dp)
      call check(error,.not. err%has_error(), "minao must build: "//err%get_full_trace())
      if (allocated(error)) return
      n = electrons(mol, d)
      call check(error, all(abs(d) < 10.0_dp) .and. abs(n - 10.0_dp) < 0.1_dp, &
                 "a pruned projection must stay a sensible density")
      call mol%destroy()
   end subroutine test_pruned

   subroutine test_heavy_refused(error)
      !! An element without data is an error, not a silent zero block
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: d(:, :)
      real(dp) :: c(3, 1)

      c = 0.0_dp
      call build_czt_molecule([54], ["Xe"], c, "def2-svp", mol, err)
      call check(error,.not. err%has_error(), "Xe must build: "//err%get_full_trace())
      if (allocated(error)) return
      call build_minao_guess(mol, d, err)
      call check(error, err%has_error(), "minao must refuse xenon")
      call mol%destroy()
   end subroutine test_heavy_refused

end module test_mqc_czt_minao

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_minao, only: collect_mqc_czt_minao_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_minao", collect_mqc_czt_minao_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
