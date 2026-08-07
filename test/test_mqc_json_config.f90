module test_mqc_json_config
   !! Guards that reading a deck as JSON gives the same config as reading the
   !! `.mqc` that `mqc_prep.py` generates from it.
   !!
   !! This is the whole safety argument for reading JSON directly. The `.mqc`
   !! path has been the only input path for the life of the code and every
   !! validation energy was produced through it, so it is the reference; the
   !! JSON reader is correct exactly insofar as it agrees with it. The repo
   !! already carries matched `.json`/`.mqc` pairs under `validation/inputs`,
   !! which makes that comparison free.
   !!
   !! The pairs below are named rather than globbed. Fortran has no portable
   !! directory listing, and naming them says what each one is *for* -- the
   !! list is a feature matrix, not a sample. A pair that stops existing fails
   !! loudly rather than silently reducing coverage.
   !!
   !! Note that `.mqc` inlines the geometry that the JSON only references by
   !! path, so agreement here also exercises the xyz resolution.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_config_types, only: mqc_config_t
   use mqc_config_parser, only: read_mqc_file
   use mqc_json_config_reader, only: read_json_config_file
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_json_config_tests

   character(len=*), parameter :: INPUT_DIR = "../validation/inputs/"
      !! Tests run with their source directory as the working directory

   real(dp), parameter :: COORD_TOL = 1.0e-9_dp
      !! The .mqc emitter writes coordinates as decimal text, so a round trip
      !! through it is not bit-exact against the .xyz the JSON points at.
      !! Anything above this is a real disagreement, not formatting.

contains

   subroutine collect_mqc_json_config_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("single_molecule_xyz_reference", test_h3o), &
                  new_unittest("fragmented_mbe_with_cutoffs", test_frag_water20), &
                  new_unittest("overlapping_fragments_and_connectivity", test_overlapping_gly3), &
                  new_unittest("multi_molecule_mode", test_multi_frag), &
                  new_unittest("dft_functional", test_dft_b3lyp), &
                  new_unittest("open_shell", test_uhf_oh), &
                  new_unittest("xtb_solvation", test_solvation), &
                  new_unittest("hessian_settings", test_hessian), &
                  new_unittest("charged_fragments", test_charged_cluster), &
                  new_unittest("missing_file_is_an_error", test_missing_file) &
                  ]
   end subroutine collect_mqc_json_config_tests

   subroutine test_h3o(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_pair("h3o", error)
   end subroutine test_h3o

   subroutine test_frag_water20(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_pair("frag_water20_hf", error)
   end subroutine test_frag_water20

   subroutine test_overlapping_gly3(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_pair("overlapping_gly3", error)
   end subroutine test_overlapping_gly3

   subroutine test_multi_frag(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_pair("multi_frag", error)
   end subroutine test_multi_frag

   subroutine test_dft_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_pair("dft_water_b3lyp", error)
   end subroutine test_dft_b3lyp

   subroutine test_uhf_oh(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_pair("uhf_oh", error)
   end subroutine test_uhf_oh

   subroutine test_solvation(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_pair("w1_water_cpcm", error)
   end subroutine test_solvation

   subroutine test_hessian(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_pair("hess_h3o", error)
   end subroutine test_hessian

   subroutine test_charged_cluster(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_pair("charged_cluster", error)
   end subroutine test_charged_cluster

   subroutine test_missing_file(error)
      !! A path that does not exist must be an error, not an empty config
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: read_error

      call read_json_config_file(INPUT_DIR//"no_such_deck.json", config, read_error)
      call check(error, read_error%has_error(), &
                 "reading a nonexistent JSON file should set an error")
   end subroutine test_missing_file

   subroutine compare_pair(stem, error)
      !! Read <stem>.json and <stem>.mqc and require the configs to agree
      character(len=*), intent(in) :: stem
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: from_json, from_mqc
      type(error_t) :: json_error, mqc_error

      call read_json_config_file(INPUT_DIR//stem//".json", from_json, json_error)
      call check(error, .not. json_error%has_error(), &
                 stem//".json failed to parse: "//json_error%get_message())
      if (allocated(error)) return

      call read_mqc_file(INPUT_DIR//stem//".mqc", from_mqc, mqc_error)
      call check(error, .not. mqc_error%has_error(), &
                 stem//".mqc failed to parse: "//mqc_error%get_message())
      if (allocated(error)) return

      call compare_configs(stem, from_json, from_mqc, error)

      call from_json%destroy()
      call from_mqc%destroy()
   end subroutine compare_pair

   subroutine compare_configs(stem, a, b, error)
      !! Field-by-field comparison of two configs
      character(len=*), intent(in) :: stem
      type(mqc_config_t), intent(in) :: a, b
      type(error_type), allocatable, intent(out) :: error

      integer :: imol

      call check_int(error, stem, "method", int(a%method), int(b%method))
      if (allocated(error)) return
      call check_int(error, stem, "calc_type", int(a%calc_type), int(b%calc_type))
      if (allocated(error)) return
      call check_int(error, stem, "index_base", a%index_base, b%index_base)
      if (allocated(error)) return

      call check_text(error, stem, "schema_name", a%schema_name, b%schema_name)
      if (allocated(error)) return
      call check_text(error, stem, "basis", a%basis, b%basis)
      if (allocated(error)) return
      call check_text(error, stem, "aux_basis", a%aux_basis, b%aux_basis)
      if (allocated(error)) return
      call check_text(error, stem, "functional", a%functional, b%functional)
      if (allocated(error)) return
      call check_text(error, stem, "log_level", a%log_level, b%log_level)
      if (allocated(error)) return
      call check_text(error, stem, "fragment_breakdown", a%fragment_breakdown, b%fragment_breakdown)
      if (allocated(error)) return
      call check_text(error, stem, "scf_guess", a%scf_guess, b%scf_guess)
      if (allocated(error)) return
      call check_text(error, stem, "solvent", a%solvent, b%solvent)
      if (allocated(error)) return
      call check_text(error, stem, "solvation_model", a%solvation_model, b%solvation_model)
      if (allocated(error)) return

      call check_int(error, stem, "scf_maxiter", a%scf_maxiter, b%scf_maxiter)
      if (allocated(error)) return
      call check_real(error, stem, "scf_tolerance", a%scf_tolerance, b%scf_tolerance)
      if (allocated(error)) return
      call check(error, a%scf_unrestricted .eqv. b%scf_unrestricted, &
                 stem//": scf_unrestricted differs")
      if (allocated(error)) return

      call check_real(error, stem, "hessian_displacement", a%hessian_displacement, &
                      b%hessian_displacement)
      if (allocated(error)) return
      call check_real(error, stem, "hessian_temperature", a%hessian_temperature, &
                      b%hessian_temperature)
      if (allocated(error)) return
      call check_real(error, stem, "hessian_pressure", a%hessian_pressure, b%hessian_pressure)
      if (allocated(error)) return

      call check_real(error, stem, "dielectric", a%dielectric, b%dielectric)
      if (allocated(error)) return
      call check_int(error, stem, "cpcm_nang", a%cpcm_nang, b%cpcm_nang)
      if (allocated(error)) return
      call check_real(error, stem, "cpcm_rscale", a%cpcm_rscale, b%cpcm_rscale)
      if (allocated(error)) return

      call check_text(error, stem, "frag_method", a%frag_method, b%frag_method)
      if (allocated(error)) return
      call check_int(error, stem, "frag_level", a%frag_level, b%frag_level)
      if (allocated(error)) return
      call check(error, a%allow_overlapping_fragments .eqv. b%allow_overlapping_fragments, &
                 stem//": allow_overlapping_fragments differs")
      if (allocated(error)) return
      call check_text(error, stem, "embedding", a%embedding, b%embedding)
      if (allocated(error)) return
      call check_text(error, stem, "cutoff_method", a%cutoff_method, b%cutoff_method)
      if (allocated(error)) return
      call check_text(error, stem, "distance_metric", a%distance_metric, b%distance_metric)
      if (allocated(error)) return
      call check_int(error, stem, "global_groups", a%global_groups, b%global_groups)
      if (allocated(error)) return
      call check_int(error, stem, "nodes_per_group", a%nodes_per_group, b%nodes_per_group)
      if (allocated(error)) return
      call check_cutoffs(error, stem, a, b)
      if (allocated(error)) return

      ! Molecules: either the single-molecule fields at the top of the config,
      ! or the molecules array. Which one is in use is itself part of what has
      ! to agree.
      call check_int(error, stem, "nmol", a%nmol, b%nmol)
      if (allocated(error)) return

      if (a%nmol == 0) then
         call check_molecule(error, stem//" (single)", &
                             a%charge, b%charge, a%multiplicity, b%multiplicity, &
                             a%geometry%natoms, b%geometry%natoms, &
                             a%geometry%elements, b%geometry%elements, &
                             a%geometry%coords, b%geometry%coords, &
                             a%nfrag, b%nfrag, a%nbonds, b%nbonds, a%nbroken, b%nbroken)
         if (allocated(error)) return
         call check_fragments(error, stem, a, b, 0)
         if (allocated(error)) return
         call check_bonds(error, stem, a, b, 0)
      else
         do imol = 1, a%nmol
            call check_molecule(error, stem//" (molecule "//int_text(imol)//")", &
                                a%molecules(imol)%charge, b%molecules(imol)%charge, &
                                a%molecules(imol)%multiplicity, b%molecules(imol)%multiplicity, &
                                a%molecules(imol)%geometry%natoms, &
                                b%molecules(imol)%geometry%natoms, &
                                a%molecules(imol)%geometry%elements, &
                                b%molecules(imol)%geometry%elements, &
                                a%molecules(imol)%geometry%coords, &
                                b%molecules(imol)%geometry%coords, &
                                a%molecules(imol)%nfrag, b%molecules(imol)%nfrag, &
                                a%molecules(imol)%nbonds, b%molecules(imol)%nbonds, &
                                a%molecules(imol)%nbroken, b%molecules(imol)%nbroken)
            if (allocated(error)) return
            call check_fragments(error, stem, a, b, imol)
            if (allocated(error)) return
            call check_bonds(error, stem, a, b, imol)
            if (allocated(error)) return
         end do
      end if
   end subroutine compare_configs

   subroutine check_molecule(error, label, charge_a, charge_b, mult_a, mult_b, &
                             natoms_a, natoms_b, elements_a, elements_b, &
                             coords_a, coords_b, nfrag_a, nfrag_b, &
                             nbonds_a, nbonds_b, nbroken_a, nbroken_b)
      !! Structure and geometry for one molecule
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      integer, intent(in) :: charge_a, charge_b, mult_a, mult_b, natoms_a, natoms_b
      character(len=*), intent(in) :: elements_a(:), elements_b(:)
      real(dp), intent(in) :: coords_a(:, :), coords_b(:, :)
      integer, intent(in) :: nfrag_a, nfrag_b, nbonds_a, nbonds_b, nbroken_a, nbroken_b

      integer :: iatom

      call check_int(error, label, "charge", charge_a, charge_b)
      if (allocated(error)) return
      call check_int(error, label, "multiplicity", mult_a, mult_b)
      if (allocated(error)) return
      call check_int(error, label, "natoms", natoms_a, natoms_b)
      if (allocated(error)) return
      call check_int(error, label, "nfrag", nfrag_a, nfrag_b)
      if (allocated(error)) return
      call check_int(error, label, "nbonds", nbonds_a, nbonds_b)
      if (allocated(error)) return
      call check_int(error, label, "nbroken", nbroken_a, nbroken_b)
      if (allocated(error)) return

      do iatom = 1, natoms_a
         call check(error, trim(elements_a(iatom)) == trim(elements_b(iatom)), &
                    label//": element "//int_text(iatom)//" differs: '"// &
                    trim(elements_a(iatom))//"' vs '"//trim(elements_b(iatom))//"'")
         if (allocated(error)) return
         call check(error, maxval(abs(coords_a(:, iatom) - coords_b(:, iatom))) < COORD_TOL, &
                    label//": coordinates of atom "//int_text(iatom)//" differ")
         if (allocated(error)) return
      end do
   end subroutine check_molecule

   subroutine check_fragments(error, stem, a, b, imol)
      !! Fragment definitions, charges and multiplicities
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: stem
      type(mqc_config_t), intent(in) :: a, b
      integer, intent(in) :: imol  !! 0 for single-molecule mode

      integer :: ifrag, nfrag

      if (imol == 0) then
         nfrag = a%nfrag
      else
         nfrag = a%molecules(imol)%nfrag
      end if
      if (nfrag <= 0) return

      do ifrag = 1, nfrag
         if (imol == 0) then
            call check_one_fragment(error, stem, ifrag, &
                                    a%fragments(ifrag)%charge, b%fragments(ifrag)%charge, &
                                    a%fragments(ifrag)%multiplicity, &
                                    b%fragments(ifrag)%multiplicity, &
                                    a%fragments(ifrag)%indices, b%fragments(ifrag)%indices)
         else
            call check_one_fragment(error, stem, ifrag, &
                                    a%molecules(imol)%fragments(ifrag)%charge, &
                                    b%molecules(imol)%fragments(ifrag)%charge, &
                                    a%molecules(imol)%fragments(ifrag)%multiplicity, &
                                    b%molecules(imol)%fragments(ifrag)%multiplicity, &
                                    a%molecules(imol)%fragments(ifrag)%indices, &
                                    b%molecules(imol)%fragments(ifrag)%indices)
         end if
         if (allocated(error)) return
      end do
   end subroutine check_fragments

   subroutine check_one_fragment(error, stem, ifrag, charge_a, charge_b, &
                                 mult_a, mult_b, indices_a, indices_b)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: stem
      integer, intent(in) :: ifrag, charge_a, charge_b, mult_a, mult_b
      integer, intent(in) :: indices_a(:), indices_b(:)

      character(len=:), allocatable :: label

      label = stem//" fragment "//int_text(ifrag)
      call check_int(error, label, "charge", charge_a, charge_b)
      if (allocated(error)) return
      call check_int(error, label, "multiplicity", mult_a, mult_b)
      if (allocated(error)) return
      call check_int(error, label, "atom count", size(indices_a), size(indices_b))
      if (allocated(error)) return
      call check(error, all(indices_a == indices_b), label//": atom indices differ")
   end subroutine check_one_fragment

   subroutine check_bonds(error, stem, a, b, imol)
      !! Connectivity, including the derived broken flag
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: stem
      type(mqc_config_t), intent(in) :: a, b
      integer, intent(in) :: imol  !! 0 for single-molecule mode

      integer :: ibond, nbonds

      if (imol == 0) then
         nbonds = a%nbonds
      else
         nbonds = a%molecules(imol)%nbonds
      end if
      if (nbonds <= 0) return

      do ibond = 1, nbonds
         if (imol == 0) then
            call check_one_bond(error, stem, ibond, a%bonds(ibond), b%bonds(ibond))
         else
            call check_one_bond(error, stem, ibond, a%molecules(imol)%bonds(ibond), &
                                b%molecules(imol)%bonds(ibond))
         end if
         if (allocated(error)) return
      end do
   end subroutine check_bonds

   subroutine check_one_bond(error, stem, ibond, bond_a, bond_b)
      use mqc_config_types, only: bond_t
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: stem
      integer, intent(in) :: ibond
      type(bond_t), intent(in) :: bond_a, bond_b

      character(len=:), allocatable :: label

      label = stem//" bond "//int_text(ibond)
      call check_int(error, label, "atom_i", bond_a%atom_i, bond_b%atom_i)
      if (allocated(error)) return
      call check_int(error, label, "atom_j", bond_a%atom_j, bond_b%atom_j)
      if (allocated(error)) return
      call check_int(error, label, "order", bond_a%order, bond_b%order)
      if (allocated(error)) return
      ! The one derived field, and the one most likely to be got wrong: it is
      ! computed from the fragment definitions, not read from the input.
      call check(error, bond_a%is_broken .eqv. bond_b%is_broken, &
                 label//": is_broken differs")
   end subroutine check_one_bond

   subroutine check_cutoffs(error, stem, a, b)
      !! Per-level distance cutoffs, allocated only when the input gives them
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: stem
      type(mqc_config_t), intent(in) :: a, b

      integer :: level

      call check(error, allocated(a%fragment_cutoffs) .eqv. allocated(b%fragment_cutoffs), &
                 stem//": one config has fragment_cutoffs and the other does not")
      if (allocated(error)) return
      if (.not. allocated(a%fragment_cutoffs)) return

      call check_int(error, stem, "fragment_cutoffs size", &
                     size(a%fragment_cutoffs), size(b%fragment_cutoffs))
      if (allocated(error)) return

      do level = 1, size(a%fragment_cutoffs)
         call check_real(error, stem, "fragment_cutoffs("//int_text(level)//")", &
                         a%fragment_cutoffs(level), b%fragment_cutoffs(level))
         if (allocated(error)) return
      end do
   end subroutine check_cutoffs

   ! ---- small comparison helpers -------------------------------------------

   subroutine check_int(error, label, field, a, b)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label, field
      integer, intent(in) :: a, b

      call check(error, a == b, label//": "//field//" differs: "// &
                 int_text(a)//" (json) vs "//int_text(b)//" (mqc)")
   end subroutine check_int

   subroutine check_real(error, label, field, a, b)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label, field
      real(dp), intent(in) :: a, b

      ! Relative where the values are large, absolute near zero -- the .mqc
      ! emitter's decimal formatting loses low bits either way.
      call check(error, abs(a - b) <= 1.0e-9_dp*max(1.0_dp, abs(a)), &
                 label//": "//field//" differs")
   end subroutine check_real

   subroutine check_text(error, label, field, a, b)
      !! Compare two optional strings, treating unallocated as absent
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label, field
      character(len=:), allocatable, intent(in) :: a, b

      call check(error, allocated(a) .eqv. allocated(b), &
                 label//": "//field//" is set in one config and not the other")
      if (allocated(error)) return
      if (.not. allocated(a)) return
      call check(error, trim(a) == trim(b), &
                 label//": "//field//" differs: '"//trim(a)//"' (json) vs '"// &
                 trim(b)//"' (mqc)")
   end subroutine check_text

   pure function int_text(value) result(text)
      integer, intent(in) :: value
      character(len=:), allocatable :: text
      character(len=16) :: buffer

      write (buffer, "(I0)") value
      text = trim(adjustl(buffer))
   end function int_text

end module test_mqc_json_config

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_json_config, only: collect_mqc_json_config_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_json_config", collect_mqc_json_config_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
