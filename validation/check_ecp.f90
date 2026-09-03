!! The ECP layout and the ECP matrix, against PySCF
!!
!! Iodine in def2-SVP with def2-ECP, and a hydrogen far enough away that the
!! overlap screen has to reject as well as accept. Both reference numbers came
!! from `mol.intor('ECPscalar_sph')` on the same system.
!!
!! **This reads mqc's own `def2-ecp.json`**, not PySCF's copy of the same
!! potential, and that is deliberate. The two agree here -- the parsed
!! exponents were checked against PySCF's channel by channel -- but a data
!! file that had drifted would otherwise surface as an integral error and be
!! chased in the wrong place, which has happened on this project before.
!!
!! Not a ctest case, like the rest of `validation/`: the reference is another
!! program's output rather than a number this one printed earlier.
program check_ecp
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_basis_utils, only: find_basis_file
   use mqc_cgto, only: molecular_basis_type
   use mqc_ecp, only: molecular_ecp_type
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_json_ecp_reader, only: build_molecular_ecp_json
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_ecp, only: ecp_matrix
   implicit none

   !! PySCF 2.14, mol.intor('ECPscalar_sph') on the system below
   real(dp), parameter :: REF_MAX = 2.54494167e+01_dp
   real(dp), parameter :: REF_TRACE = 9.05132744e+01_dp
   real(dp), parameter :: TOL = 1.0e-7_dp

   type(error_t) :: err
   type(molecular_basis_type) :: basis
   type(molecular_ecp_type) :: ecp
   type(libcint_molecule_t) :: mol
   real(dp), allocatable :: m(:, :)
   real(dp) :: xyz(3, 2)
   real(dp) :: got_max, got_trace, asym
   character(len=8) :: syms(2)
   character(len=:), allocatable :: bfile, efile
   integer :: z(2)
   integer :: i, bad

   bad = 0
   syms = ["I ", "H "]
   z = [53, 1]
   xyz = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                  0.0_dp, 0.0_dp, 25.0_dp], [3, 2])

   call find_basis_file("def2-svp", bfile, err)
   if (err%has_error()) then
      write (*, "(a,a)") "FAIL basis file: ", err%get_message()
      stop 1
   end if
   call find_basis_file("def2-ecp", efile, err)
   if (err%has_error()) then
      write (*, "(a,a)") "FAIL ecp file: ", err%get_message()
      stop 1
   end if

   call build_molecular_basis_json(bfile, syms, basis, err)
   if (err%has_error()) then
      write (*, "(a,a)") "FAIL basis: ", err%get_message()
      stop 1
   end if
   call build_molecular_ecp_json(efile, syms, ecp, err)
   if (err%has_error()) then
      write (*, "(a,a)") "FAIL ecp: ", err%get_message()
      stop 1
   end if

   ! The channels, before any integral. A projected channel that came back
   ! claiming l = -1 still integrates to *something*, and the something looks
   ! plausible -- which is exactly how a reader fault stayed hidden until an
   ! integral consumed it. Asserting the layout separates the two.
   if (ecp%atoms(1)%local%ang_mom /= 3) then
      write (*, "(a,i0)") "FAIL local channel l should be 3, got ", &
         ecp%atoms(1)%local%ang_mom
      bad = bad + 1
   end if
   if (ecp%atoms(1)%n_projected /= 3) then
      write (*, "(a,i0)") "FAIL iodine should carry 3 projected channels, got ", &
         ecp%atoms(1)%n_projected
      bad = bad + 1
   else
      do i = 1, 3
         if (ecp%atoms(1)%projected(i)%ang_mom /= i - 1) then
            write (*, "(a,i0,a,i0)") "FAIL projected channel ", i, " has l = ", &
               ecp%atoms(1)%projected(i)%ang_mom
            bad = bad + 1
         end if
      end do
   end if
   if (ecp%atoms(2)%has_ecp) then
      write (*, "(a)") "FAIL hydrogen should carry no ECP"
      bad = bad + 1
   end if

   call mol%build(z, xyz, basis, err, ecp=ecp)
   if (err%has_error()) then
      write (*, "(a,a)") "FAIL build: ", err%get_message()
      stop 1
   end if

   ! Z minus the electrons the potential replaced. The nuclear attraction and
   ! the nuclear repulsion both read this, so leaving the full 53 in place
   ! gives a calculation that converges and is wrong by hundreds of Hartree
   ! with nothing looking amiss.
   if (abs(mol%charges(1) - 25.0_dp) > TOL) then
      write (*, "(a,f8.3)") "FAIL iodine effective charge should be 25, got ", &
         mol%charges(1)
      bad = bad + 1
   end if
   if (mol%core_electrons(1) /= 28) then
      write (*, "(a,i0)") "FAIL iodine core electrons should be 28, got ", &
         mol%core_electrons(1)
      bad = bad + 1
   end if
   if (mol%necpbas <= 0) then
      write (*, "(a)") "FAIL no ECP rows were laid out"
      bad = bad + 1
   end if

   call ecp_matrix(mol%nao, mol%nbas, mol%natm, mol%cartesian, mol%atm, &
                   mol%bas_with_ecp, mol%env, mol%shell_offset, mol%necpbas, m)

   got_max = maxval(abs(m))
   got_trace = sum([(m(i, i), i=1, mol%nao)])
   asym = maxval(abs(m - transpose(m)))

   write (*, "(a,i0,a,i0,a,i0)") "check_ecp: nao ", mol%nao, "  nbas ", mol%nbas, &
      "  ecp rows ", mol%necpbas
   write (*, "(a,es16.8,a,es16.8)") "  max |ECP| ", got_max, "   PySCF ", REF_MAX
   write (*, "(a,es16.8,a,es16.8)") "  trace     ", got_trace, "   PySCF ", REF_TRACE
   write (*, "(a,es16.8)") "  asymmetry ", asym

   if (abs(got_max - REF_MAX) > TOL) then
      write (*, "(a)") "FAIL max |ECP| disagrees with PySCF"
      bad = bad + 1
   end if
   if (abs(got_trace - REF_TRACE) > TOL) then
      write (*, "(a)") "FAIL trace disagrees with PySCF"
      bad = bad + 1
   end if
   ! The ECP operator is symmetric, and the two scalars above can both match
   ! while a shell-pair block sits in the wrong place.
   if (asym > 1.0e-12_dp) then
      write (*, "(a)") "FAIL the matrix is not symmetric"
      bad = bad + 1
   end if

   if (bad > 0) then
      write (*, "(a,i0,a)") "check_ecp: FAILED (", bad, " problem(s))"
      stop 1
   end if
   write (*, "(a)") "check_ecp: OK"
end program check_ecp
