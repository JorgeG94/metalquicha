!! Natural orbitals of a one-particle density, and what their occupations say
module mqc_czt_natural_orbitals
   !! The orbitals that diagonalise a density, for any density.
   !!
   !! A single closed-shell determinant has occupations of exactly two and
   !! zero and its natural orbitals say nothing its canonical ones did not.
   !! Every other wave function here has a density that is not idempotent --
   !! an unrestricted reference summed over spin, a correlated relaxed
   !! density -- and for those the occupations are the compact answer to
   !! "how far is this from one determinant, and in which orbitals".
   !!
   !! `mqc_czt_mcscf`'s `natural_orbitals` does this for an active-space
   !! density given in the orbital basis. This module takes the density in the
   !! atomic-orbital basis instead, which is the form every other source
   !! arrives in.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: natural_occupation_summary_t
   public :: density_natural_orbitals
   public :: summarize_natural_occupations
   public :: print_natural_orbital_report
   public :: FRACTIONAL_OCCUPATION

   real(dp), parameter :: OVERLAP_FLOOR = 1.0e-8_dp
      !! Overlap eigenvalues at or below this are dropped from the
      !! orthogonalisation, relative to the largest
   real(dp), parameter :: FRACTIONAL_OCCUPATION = 0.02_dp
      !! An orbital whose occupation is further than this from both zero and
      !! two is counted as fractionally occupied. The conventional window for
      !! choosing an active space from natural occupations.
   integer, parameter :: LINE_LEN = 160

   type :: natural_occupation_summary_t
      !! What a set of natural occupations says about the wave function
      real(dp) :: n_electrons = 0.0_dp
         !! Their sum, which is the trace of the density with the overlap
      real(dp) :: unpaired_linear = 0.0_dp
         !! `sum_i min(n_i, 2 - n_i)`: Head-Gordon's count of effectively
         !! unpaired electrons. Zero for a closed-shell determinant, and the
         !! number of singly occupied orbitals for a high-spin one.
      real(dp) :: unpaired_quadratic = 0.0_dp
         !! `sum_i n_i^2 (2 - n_i)^2`, the same count weighted so that
         !! occupations near zero and two contribute at second order rather
         !! than first. The one to read when dynamic correlation is present.
      integer :: n_fractional = 0
         !! Orbitals with `FRACTIONAL_OCCUPATION < n < 2 - FRACTIONAL_OCCUPATION`
      real(dp) :: largest_deviation = 0.0_dp
         !! `max_i min(n_i, 2 - n_i)`
   end type natural_occupation_summary_t

contains

   subroutine density_natural_orbitals(density, overlap, natural, occupations, error)
      !! Natural orbitals and occupations of a density in the atomic-orbital basis
      !!
      !! Solves `S D S c = n S c`. The overlap is orthogonalised canonically,
      !! so near-null directions of the basis are dropped rather than divided
      !! by, and there are then fewer orbitals than basis functions.
      real(dp), intent(in) :: density(:, :)
         !! (n_ao, n_ao), symmetric. The total density, both spins summed.
      real(dp), intent(in) :: overlap(:, :)   !! (n_ao, n_ao)
      real(dp), allocatable, intent(out) :: natural(:, :)
         !! (n_ao, n_kept), orthonormal under `overlap`, most occupied first
      real(dp), allocatable, intent(out) :: occupations(:)
         !! (n_kept), descending
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: vectors(:, :), values(:), x(:, :), sx(:, :)
      real(dp), allocatable :: work(:, :), reduced(:, :), occ(:), sorted(:, :)
      integer :: n, kept, i, info

      if (error%has_error()) return
      n = size(density, 1)
      if (size(density, 2) /= n .or. size(overlap, 1) /= n .or. size(overlap, 2) /= n) then
         call error%set(ERROR_VALIDATION, "natural orbitals: the density and the "// &
                        "overlap are not the same square size")
         return
      end if

      allocate (vectors(n, n), values(n))
      vectors = overlap
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "natural orbitals: the overlap matrix "// &
                        "would not diagonalise")
         return
      end if
      kept = count(values > OVERLAP_FLOOR*maxval(values))
      if (kept < 1) then
         call error%set(ERROR_VALIDATION, "natural orbitals: the overlap matrix "// &
                        "has no positive eigenvalue")
         return
      end if
      ! Ascending eigenvalues, so the kept ones are the last `kept` columns.
      allocate (x(n, kept))
      do i = 1, kept
         x(:, i) = vectors(:, n - kept + i)/sqrt(values(n - kept + i))
      end do
      deallocate (vectors, values)

      ! X^T S D S X in the orthonormal basis X.
      allocate (sx(n, kept), work(n, kept), reduced(kept, kept))
      call pic_gemm(overlap, x, sx, beta=0.0_dp)
      call pic_gemm(density, sx, work, beta=0.0_dp)
      call pic_gemm(sx, work, reduced, transa="T", beta=0.0_dp)
      reduced = 0.5_dp*(reduced + transpose(reduced))

      allocate (occ(kept))
      call pic_syev(reduced, occ, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "natural orbitals: the density would not "// &
                        "diagonalise")
         return
      end if

      ! Ascending out of the solver, descending out of here.
      allocate (sorted(kept, kept), occupations(kept))
      do i = 1, kept
         sorted(:, i) = reduced(:, kept - i + 1)
         occupations(i) = occ(kept - i + 1)
      end do
      allocate (natural(n, kept))
      call pic_gemm(x, sorted, natural, beta=0.0_dp)
   end subroutine density_natural_orbitals

   pure subroutine summarize_natural_occupations(occupations, summary)
      !! The electron count, the unpaired-electron counts and the fractional set
      real(dp), intent(in) :: occupations(:)
      type(natural_occupation_summary_t), intent(out) :: summary

      real(dp) :: n, dev
      integer :: i

      do i = 1, size(occupations)
         ! Rounding leaves occupations a hair outside [0, 2]; inside the
         ! counts that is zero deviation, not a negative one.
         n = min(max(occupations(i), 0.0_dp), 2.0_dp)
         dev = min(n, 2.0_dp - n)
         summary%n_electrons = summary%n_electrons + occupations(i)
         summary%unpaired_linear = summary%unpaired_linear + dev
         summary%unpaired_quadratic = summary%unpaired_quadratic + (n*(2.0_dp - n))**2
         summary%largest_deviation = max(summary%largest_deviation, dev)
         if (dev > FRACTIONAL_OCCUPATION) summary%n_fractional = summary%n_fractional + 1
      end do
   end subroutine summarize_natural_occupations

   subroutine print_natural_orbital_report(occupations, summary, source)
      !! The occupations that are not two or zero, and the counts above
      real(dp), intent(in) :: occupations(:)
      type(natural_occupation_summary_t), intent(in) :: summary
      character(len=*), intent(in) :: source
         !! Which density these are the natural orbitals of, for the heading

      real(dp), parameter :: PRINT_FLOOR = 1.0e-3_dp
         !! Occupations within this of two or of zero are counted, not listed
      real(dp), parameter :: NEGATIVE_NOTICE = 1.0e-6_dp
         !! How far outside [0, 2] an occupation must be before it is mentioned
      character(len=LINE_LEN) :: line
      integer :: i, n_full, n_empty

      call logger%info("  natural orbitals of the "//trim(source)//" density")
      write (line, "(a,f14.8)") "    electrons                         ", summary%n_electrons
      call logger%info(trim(line))
      write (line, "(a,f14.8)") "    effectively unpaired (linear)     ", summary%unpaired_linear
      call logger%info(trim(line))
      write (line, "(a,f14.8)") "    effectively unpaired (quadratic)  ", &
         summary%unpaired_quadratic
      call logger%info(trim(line))
      write (line, "(a,i0,a,f4.2,a,f4.2)") "    ", summary%n_fractional, &
         " orbitals with occupation between ", FRACTIONAL_OCCUPATION, " and ", &
         2.0_dp - FRACTIONAL_OCCUPATION
      call logger%info(trim(line))
      ! A relaxed correlated density is not the density of any wave function,
      ! so its occupations are not bound to [0, 2]. Said when it happens: a
      ! reader comparing against a variational density should know.
      if (minval(occupations) < -NEGATIVE_NOTICE .or. &
          maxval(occupations) > 2.0_dp + NEGATIVE_NOTICE) then
         write (line, "(a,f12.8,a,f12.8)") "    occupations run outside [0, 2]: lowest ", &
            minval(occupations), ", highest ", maxval(occupations)
         call logger%info(trim(line))
      end if

      n_full = count(occupations > 2.0_dp - PRINT_FLOOR)
      n_empty = count(occupations < PRINT_FLOOR)
      write (line, "(a,i0,a,i0,a)") "    ", n_full, " within 0.001 of two, ", n_empty, &
         " within 0.001 of zero; the rest:"
      call logger%info(trim(line))
      call logger%info("      orbital    occupation")
      do i = 1, size(occupations)
         if (occupations(i) > 2.0_dp - PRINT_FLOOR .or. occupations(i) < PRINT_FLOOR) cycle
         write (line, "(i13,f14.8)") i, occupations(i)
         call logger%info(trim(line))
      end do
   end subroutine print_natural_orbital_report

end module mqc_czt_natural_orbitals
