!! An effective fragment potential as a flat pair of buffers
module mqc_czt_efp_serialize
   !! `efp_fragment_t` flattened into one integer array and one real array, and
   !! rebuilt from them.
   !!
   !! **Why this exists.** EFMO distributes its monomers over MPI ranks because
   !! MAKEFP is the cost of the method, but *every* rank then needs *every*
   !! fragment's potential: the far pairs and the one induction over all
   !! fragments are not decomposable by owner. So each potential has to travel,
   !! and a derived type with thirty allocatable components of four different
   !! ranks does not travel through an MPI reduction on its own.
   !!
   !! **A sum-reduction over buffers that are zero where a rank computed
   !! nothing** is the shape used, the same one `mqc_czt_fmo` exchanges its
   !! monomer densities with. That needs the *layout* to be known everywhere
   !! before the contents are, which is why this is split in two: a fixed-size
   !! header of counts and flags, reduced first, and a body whose length every
   !! rank then computes from the reduced headers.
   !!
   !! **Exact, not nearly.** The alternative -- writing a `.efp` file and
   !! reading it back on the other ranks -- exists and is tested, but the format
   !! carries eight decimals on the orbital blocks, so two ranks would disagree
   !! in the eighth digit of the exchange repulsion. What goes through here is
   !! the bits.
   use pic_types, only: dp
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_efp_read, only: efp_fragment_t
   implicit none
   private

   public :: EFP_HEADER_INTS
   public :: fragment_header
   public :: fragment_buffer_sizes
   public :: fragment_pack
   public :: fragment_unpack

   integer, parameter :: EFP_HEADER_INTS = 54
      !! Counts, presence flags and one allocation flag per component array.
      !!
      !! Fixed size, because it is reduced before anything about the fragment is
      !! known on a rank that did not build it. The allocation flags are carried
      !! rather than derived from the counts: a potential written without, say,
      !! the octupoles has `n_points > 0` and no octupole array, and rebuilding
      !! one full of zeros would be a different fragment.

   integer, parameter :: N_REAL_ARRAYS = 23
   integer, parameter :: N_INT_ARRAYS = 4

   integer, parameter :: LABEL_LEN = 8
      !! Width of `labels`, which is what `efp_fragment_t` declares.

contains

   subroutine fragment_header(frag, header)
      !! The counts and flags that fix a fragment's layout
      type(efp_fragment_t), intent(in) :: frag
      integer, intent(out) :: header(EFP_HEADER_INTS)

      header = 0
      header(1) = frag%n_points
      header(2) = frag%n_atoms
      header(3) = frag%multiplicity
      header(4) = frag%n_lmo
      header(5) = frag%n_freq
      header(6) = frag%n_pol
      header(7) = frag%n_dipquad
      header(8) = frag%n_quadquad
      header(9) = frag%n_shells
      ! The primitive count has no field of its own on the type -- the exponent
      ! array's length *is* it -- so it is measured here.
      if (allocated(frag%prim_expo)) header(10) = size(frag%prim_expo)
      header(11) = frag%n_lmo_proj
      header(12) = frag%nao_proj
      header(13) = frag%n_occ_ct
      header(14) = frag%n_mo_ct
      header(15) = flag(frag%has_screen)
      header(16) = flag(frag%has_screen2)
      header(17) = flag(frag%has_dynamic)
      header(18) = flag(frag%has_static_pol)
      header(19) = flag(frag%has_dipquad)
      header(20) = flag(frag%has_quadquad)
      header(21) = flag(frag%has_basis)
      header(22) = flag(frag%has_lmo)
      header(23) = flag(frag%has_fock)
      header(24) = flag(frag%has_ctvec)
      header(25) = flag(frag%has_ctfok)
      if (allocated(frag%name)) header(26) = len(frag%name)

      header(27) = flag(allocated(frag%points))
      header(28) = flag(allocated(frag%mass))
      header(29) = flag(allocated(frag%charge))
      header(30) = flag(allocated(frag%q_elec))
      header(31) = flag(allocated(frag%q_nuc))
      header(32) = flag(allocated(frag%dipole))
      header(33) = flag(allocated(frag%quadrupole))
      header(34) = flag(allocated(frag%octopole))
      header(35) = flag(allocated(frag%screen))
      header(36) = flag(allocated(frag%screen2))
      header(37) = flag(allocated(frag%dyn_pol))
      header(38) = flag(allocated(frag%centroids))
      header(39) = flag(allocated(frag%frequencies))
      header(40) = flag(allocated(frag%static_pol))
      header(41) = flag(allocated(frag%pol_points))
      header(42) = flag(allocated(frag%dipquad))
      header(43) = flag(allocated(frag%quadquad))
      header(44) = flag(allocated(frag%prim_expo))
      header(45) = flag(allocated(frag%prim_coef))
      header(46) = flag(allocated(frag%lmo_gamess))
      header(47) = flag(allocated(frag%fock_lmo))
      header(48) = flag(allocated(frag%ctvec_gamess))
      header(49) = flag(allocated(frag%eps_occ))

      header(50) = flag(allocated(frag%shell_atom))
      header(51) = flag(allocated(frag%shell_l))
      header(52) = flag(allocated(frag%shell_first))
      header(53) = flag(allocated(frag%shell_nprim))
      header(54) = flag(allocated(frag%labels))
   end subroutine fragment_header

   pure function flag(present_) result(k)
      !! A logical as the 0 or 1 an integer reduction can carry
      logical, intent(in) :: present_
      integer :: k

      k = 0
      if (present_) k = 1
   end function flag

   subroutine fragment_buffer_sizes(header, n_ints, n_reals)
      !! How long the body buffers are for a fragment of this shape
      !!
      !! Derived from the header alone, so a rank that never saw the fragment
      !! computes the same two numbers as the rank that built it. Keep the
      !! order here identical to the walk in `pack_body`.
      integer, intent(in) :: header(EFP_HEADER_INTS)
      integer, intent(out) :: n_ints, n_reals

      integer :: np, nl, nf, npol, nsh, nprim, nlp, nao, nocc, nmo

      np = header(1)
      nl = header(4)
      nf = header(5)
      npol = header(6)
      nsh = header(9)
      nprim = header(10)
      nlp = header(11)
      nao = header(12)
      nocc = header(13)
      nmo = header(14)

      n_ints = header(26)                              ! name, one integer per character
      n_ints = n_ints + header(54)*np*LABEL_LEN        ! labels, likewise
      n_ints = n_ints + (header(50) + header(51) + header(52) + header(53))*nsh

      n_reals = header(27)*3*np &
                + header(28)*np + header(29)*np + header(30)*np + header(31)*np &
                + header(32)*3*np + header(33)*6*np + header(34)*10*np &
                + header(35)*np + header(36)*np &
                + header(37)*9*nl*nf + header(38)*3*nl + header(39)*nf &
                + header(40)*9*npol + header(41)*3*npol &
                + header(42)*header(7)*nl*nf + header(43)*header(8)*nl*nf &
                + header(44)*nprim + header(45)*nprim &
                + header(46)*nao*nlp + header(47)*nlp*nlp &
                + header(48)*nao*nmo + header(49)*nocc
   end subroutine fragment_buffer_sizes

   subroutine fragment_pack(frag, header, ibuf, rbuf, error)
      !! Write a fragment's contents into the two body buffers
      !!
      !! `ibuf` and `rbuf` are exactly the lengths `fragment_buffer_sizes`
      !! reported for `header`, which must be this fragment's own.
      type(efp_fragment_t), intent(in) :: frag
      integer, intent(in) :: header(EFP_HEADER_INTS)
      integer, intent(out) :: ibuf(:)
      real(dp), intent(out) :: rbuf(:)
      type(error_t), intent(inout) :: error

      integer :: n_ints, n_reals, ai, ar, i, k

      call fragment_buffer_sizes(header, n_ints, n_reals)
      if (size(ibuf) /= n_ints .or. size(rbuf) /= n_reals) then
         call error%set(ERROR_VALIDATION, "efp serialize: the buffers are "// &
                        to_char(size(ibuf))//" and "//to_char(size(rbuf))// &
                        " long where the header asks for "//to_char(n_ints)// &
                        " and "//to_char(n_reals))
         return
      end if

      ai = 1
      ar = 1
      if (header(26) > 0) then
         do i = 1, header(26)
            ibuf(ai) = iachar(frag%name(i:i))
            ai = ai + 1
         end do
      end if
      if (header(54) == 1) then
         do i = 1, header(1)
            do k = 1, LABEL_LEN
               ibuf(ai) = iachar(frag%labels(i) (k:k))
               ai = ai + 1
            end do
         end do
      end if
      call put_int(ibuf, ai, frag%shell_atom, header(50))
      call put_int(ibuf, ai, frag%shell_l, header(51))
      call put_int(ibuf, ai, frag%shell_first, header(52))
      call put_int(ibuf, ai, frag%shell_nprim, header(53))

      ! Every array flattened with `reshape`, in the order
      ! `fragment_buffer_sizes` counts them and `fragment_unpack` reads them.
      if (header(27) == 1) call put_real(rbuf, ar, reshape(frag%points, [size(frag%points)]))
      if (header(28) == 1) call put_real(rbuf, ar, frag%mass)
      if (header(29) == 1) call put_real(rbuf, ar, frag%charge)
      if (header(30) == 1) call put_real(rbuf, ar, frag%q_elec)
      if (header(31) == 1) call put_real(rbuf, ar, frag%q_nuc)
      if (header(32) == 1) call put_real(rbuf, ar, reshape(frag%dipole, [size(frag%dipole)]))
      if (header(33) == 1) call put_real(rbuf, ar, reshape(frag%quadrupole, [size(frag%quadrupole)]))
      if (header(34) == 1) call put_real(rbuf, ar, reshape(frag%octopole, [size(frag%octopole)]))
      if (header(35) == 1) call put_real(rbuf, ar, frag%screen)
      if (header(36) == 1) call put_real(rbuf, ar, frag%screen2)
      if (header(37) == 1) call put_real(rbuf, ar, reshape(frag%dyn_pol, [size(frag%dyn_pol)]))
      if (header(38) == 1) call put_real(rbuf, ar, reshape(frag%centroids, [size(frag%centroids)]))
      if (header(39) == 1) call put_real(rbuf, ar, frag%frequencies)
      if (header(40) == 1) call put_real(rbuf, ar, reshape(frag%static_pol, [size(frag%static_pol)]))
      if (header(41) == 1) call put_real(rbuf, ar, reshape(frag%pol_points, [size(frag%pol_points)]))
      if (header(42) == 1) call put_real(rbuf, ar, reshape(frag%dipquad, [size(frag%dipquad)]))
      if (header(43) == 1) call put_real(rbuf, ar, reshape(frag%quadquad, [size(frag%quadquad)]))
      if (header(44) == 1) call put_real(rbuf, ar, frag%prim_expo)
      if (header(45) == 1) call put_real(rbuf, ar, frag%prim_coef)
      if (header(46) == 1) call put_real(rbuf, ar, reshape(frag%lmo_gamess, [size(frag%lmo_gamess)]))
      if (header(47) == 1) call put_real(rbuf, ar, reshape(frag%fock_lmo, [size(frag%fock_lmo)]))
      if (header(48) == 1) call put_real(rbuf, ar, reshape(frag%ctvec_gamess, [size(frag%ctvec_gamess)]))
      if (header(49) == 1) call put_real(rbuf, ar, frag%eps_occ)
   end subroutine fragment_pack

   subroutine put_int(buf, at, values, present_)
      !! Append an integer array, if it is there, and advance the cursor
      integer, intent(inout) :: buf(:)
      integer, intent(inout) :: at
      integer, allocatable, intent(in) :: values(:)
      integer, intent(in) :: present_

      if (present_ /= 1) return
      buf(at:at + size(values) - 1) = values
      at = at + size(values)
   end subroutine put_int

   subroutine put_real(buf, at, values)
      !! Append a flattened real array and advance the cursor
      real(dp), intent(inout) :: buf(:)
      integer, intent(inout) :: at
      real(dp), intent(in) :: values(:)

      buf(at:at + size(values) - 1) = values
      at = at + size(values)
   end subroutine put_real

   subroutine fragment_unpack(header, ibuf, rbuf, frag, error)
      !! Rebuild a fragment from its header and body
      !!
      !! Field for field the inverse of `fragment_pack`, in the same order.
      !! Anything the header says was absent is left unallocated rather than
      !! allocated empty, so `allocated(...)` answers on the far rank what it
      !! answered on the near one.
      integer, intent(in) :: header(EFP_HEADER_INTS)
      integer, intent(in) :: ibuf(:)
      real(dp), intent(in) :: rbuf(:)
      type(efp_fragment_t), intent(out) :: frag
      type(error_t), intent(inout) :: error

      integer :: n_ints, n_reals, ai, ar, i, k
      integer :: np, nl, nf, npol, nsh, nprim, nlp, nao, nocc, nmo
      character(len=LABEL_LEN) :: label

      call fragment_buffer_sizes(header, n_ints, n_reals)
      if (size(ibuf) /= n_ints .or. size(rbuf) /= n_reals) then
         call error%set(ERROR_VALIDATION, "efp serialize: the buffers are "// &
                        to_char(size(ibuf))//" and "//to_char(size(rbuf))// &
                        " long where the header asks for "//to_char(n_ints)// &
                        " and "//to_char(n_reals))
         return
      end if

      np = header(1)
      nl = header(4)
      nf = header(5)
      npol = header(6)
      nsh = header(9)
      nprim = header(10)
      nlp = header(11)
      nao = header(12)
      nocc = header(13)
      nmo = header(14)

      frag%n_points = np
      frag%n_atoms = header(2)
      frag%multiplicity = header(3)
      frag%n_lmo = nl
      frag%n_freq = nf
      frag%n_pol = npol
      frag%n_dipquad = header(7)
      frag%n_quadquad = header(8)
      frag%n_shells = nsh
      frag%n_lmo_proj = nlp
      frag%nao_proj = nao
      frag%n_occ_ct = nocc
      frag%n_mo_ct = nmo
      frag%has_screen = header(15) == 1
      frag%has_screen2 = header(16) == 1
      frag%has_dynamic = header(17) == 1
      frag%has_static_pol = header(18) == 1
      frag%has_dipquad = header(19) == 1
      frag%has_quadquad = header(20) == 1
      frag%has_basis = header(21) == 1
      frag%has_lmo = header(22) == 1
      frag%has_fock = header(23) == 1
      frag%has_ctvec = header(24) == 1
      frag%has_ctfok = header(25) == 1

      ai = 1
      ar = 1
      if (header(26) > 0) then
         allocate (character(len=header(26)) :: frag%name)
         do i = 1, header(26)
            frag%name(i:i) = achar(ibuf(ai))
            ai = ai + 1
         end do
      end if
      if (header(54) == 1) then
         allocate (frag%labels(np))
         do i = 1, np
            do k = 1, LABEL_LEN
               label(k:k) = achar(ibuf(ai))
               ai = ai + 1
            end do
            frag%labels(i) = label
         end do
      end if
      call take_int(ibuf, ai, frag%shell_atom, nsh, header(50))
      call take_int(ibuf, ai, frag%shell_l, nsh, header(51))
      call take_int(ibuf, ai, frag%shell_first, nsh, header(52))
      call take_int(ibuf, ai, frag%shell_nprim, nsh, header(53))

      call take_2d(rbuf, ar, frag%points, 3, np, header(27))
      call take_1d(rbuf, ar, frag%mass, np, header(28))
      call take_1d(rbuf, ar, frag%charge, np, header(29))
      call take_1d(rbuf, ar, frag%q_elec, np, header(30))
      call take_1d(rbuf, ar, frag%q_nuc, np, header(31))
      call take_2d(rbuf, ar, frag%dipole, 3, np, header(32))
      call take_2d(rbuf, ar, frag%quadrupole, 6, np, header(33))
      call take_2d(rbuf, ar, frag%octopole, 10, np, header(34))
      call take_1d(rbuf, ar, frag%screen, np, header(35))
      call take_1d(rbuf, ar, frag%screen2, np, header(36))
      call take_4d(rbuf, ar, frag%dyn_pol, 3, 3, nl, nf, header(37))
      call take_2d(rbuf, ar, frag%centroids, 3, nl, header(38))
      call take_1d(rbuf, ar, frag%frequencies, nf, header(39))
      call take_3d(rbuf, ar, frag%static_pol, 3, 3, npol, header(40))
      call take_2d(rbuf, ar, frag%pol_points, 3, npol, header(41))
      call take_3d(rbuf, ar, frag%dipquad, header(7), nl, nf, header(42))
      call take_3d(rbuf, ar, frag%quadquad, header(8), nl, nf, header(43))
      call take_1d(rbuf, ar, frag%prim_expo, nprim, header(44))
      call take_1d(rbuf, ar, frag%prim_coef, nprim, header(45))
      call take_2d(rbuf, ar, frag%lmo_gamess, nao, nlp, header(46))
      call take_2d(rbuf, ar, frag%fock_lmo, nlp, nlp, header(47))
      call take_2d(rbuf, ar, frag%ctvec_gamess, nao, nmo, header(48))
      call take_1d(rbuf, ar, frag%eps_occ, nocc, header(49))
   end subroutine fragment_unpack

   subroutine take_int(buf, at, values, n, present_)
      !! Read back an integer array of `n` entries, or leave it unallocated
      integer, intent(in) :: buf(:)
      integer, intent(inout) :: at
      integer, allocatable, intent(out) :: values(:)
      integer, intent(in) :: n, present_

      if (present_ /= 1) return
      allocate (values(n))
      values = buf(at:at + n - 1)
      at = at + n
   end subroutine take_int

   subroutine take_1d(buf, at, values, n, present_)
      real(dp), intent(in) :: buf(:)
      integer, intent(inout) :: at
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in) :: n, present_

      if (present_ /= 1) return
      allocate (values(n))
      values = buf(at:at + n - 1)
      at = at + n
   end subroutine take_1d

   subroutine take_2d(buf, at, values, d1, d2, present_)
      real(dp), intent(in) :: buf(:)
      integer, intent(inout) :: at
      real(dp), allocatable, intent(out) :: values(:, :)
      integer, intent(in) :: d1, d2, present_

      integer :: n

      if (present_ /= 1) return
      n = d1*d2
      allocate (values(d1, d2))
      values = reshape(buf(at:at + n - 1), [d1, d2])
      at = at + n
   end subroutine take_2d

   subroutine take_3d(buf, at, values, d1, d2, d3, present_)
      real(dp), intent(in) :: buf(:)
      integer, intent(inout) :: at
      real(dp), allocatable, intent(out) :: values(:, :, :)
      integer, intent(in) :: d1, d2, d3, present_

      integer :: n

      if (present_ /= 1) return
      n = d1*d2*d3
      allocate (values(d1, d2, d3))
      values = reshape(buf(at:at + n - 1), [d1, d2, d3])
      at = at + n
   end subroutine take_3d

   subroutine take_4d(buf, at, values, d1, d2, d3, d4, present_)
      real(dp), intent(in) :: buf(:)
      integer, intent(inout) :: at
      real(dp), allocatable, intent(out) :: values(:, :, :, :)
      integer, intent(in) :: d1, d2, d3, d4, present_

      integer :: n

      if (present_ /= 1) return
      n = d1*d2*d3*d4
      allocate (values(d1, d2, d3, d4))
      values = reshape(buf(at:at + n - 1), [d1, d2, d3, d4])
      at = at + n
   end subroutine take_4d

end module mqc_czt_efp_serialize
