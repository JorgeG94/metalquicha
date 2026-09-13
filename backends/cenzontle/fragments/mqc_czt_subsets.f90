!! Subsets of a fragment list, and the many-body difference over them
module mqc_czt_subsets
   !! The bookkeeping every many-body expansion in this backend shares: which
   !! groups of fragments to compute, and how to turn a quantity computed *on*
   !! each group into the part that group alone contributes.
   !!
   !! **The difference operator.** For any function `f` of a group of fragments,
   !!
   !!     df_S = f(S) - sum over proper non-empty subsets T of S of df_T
   !!
   !! so that `sum over S of df_S` over *every* subset of the system is exactly
   !! `f(whole system)`, and truncating the sum at `|S| <= n` is the n-body
   !! expansion of `f`. At `|S| = 1` it is `f(I)`, at `|S| = 2` it is
   !! `f(IJ) - f(I) - f(J)`, at `|S| = 3` the usual three-body form. Nothing
   !! here knows what `f` is: [[mqc_czt_fmo]] applies it to an n-mer's internal
   !! energy and [[mqc_czt_efmo]] applies it twice, to an n-mer's in-vacuo
   !! energy and to the induction energy of the same n-mer's potentials.
   !!
   !! **Order matters.** `enumerate_subsets` returns groups smallest first, and
   !! `subtract_subsets` relies on it: a group's own difference is only correct
   !! once every one of its subsets is final.
   !!
   !! **Cost is the binomial and nothing else.** There are `C(N, n)` groups of
   !! size `n`, so level three on twenty fragments is 1140 groups against 190
   !! for level two. No level is refused here; a caller that wants to warn has
   !! the count in hand before it computes anything.
   !!
   !! This lives apart from either method because both need it and neither
   !! should have to `use` the other: EFMO reaching into [[mqc_czt_fmo]] for the
   !! enumerator would drag in the adjusted frozen orbitals, the ESP matrices
   !! and the Fock projector, none of which it has any use for.
   use pic_types, only: dp
   implicit none
   private

   public :: enumerate_subsets
   public :: subtract_subsets
   public :: is_subset
   public :: n_choose

contains

   subroutine enumerate_subsets(n_frag, level, terms, term_size, n_terms)
      !! Every combination of fragments from one up to `level`, smallest first
      !!
      !! Size order matters: the difference for a group is reduced by its
      !! subsets, so each subset has to be final before anything containing it
      !! is touched.
      integer, intent(in) :: n_frag, level
      integer, allocatable, intent(out) :: terms(:, :), term_size(:)
         !! `terms(1:term_size(t), t)` are the fragments of group `t`, ascending
      integer, intent(out) :: n_terms

      integer, allocatable :: pick(:)
      integer :: m, total, k

      total = 0
      do m = 1, level
         total = total + n_choose(n_frag, m)
      end do
      allocate (terms(max(level, 1), max(total, 1)), source=0)
      allocate (term_size(max(total, 1)), source=0)

      n_terms = 0
      do m = 1, level
         allocate (pick(m))
         do k = 1, m
            pick(k) = k
         end do
         do
            n_terms = n_terms + 1
            terms(1:m, n_terms) = pick
            term_size(n_terms) = m
            if (.not. step_combination(pick, m, n_frag)) exit
         end do
         deallocate (pick)
      end do
   end subroutine enumerate_subsets

   function step_combination(pick, m, n) result(more)
      !! Advance a combination in lexicographic order; false when exhausted
      integer, intent(inout) :: pick(:)
      integer, intent(in) :: m, n
      logical :: more

      integer :: i, k

      more = .false.
      do i = m, 1, -1
         if (pick(i) < n - m + i) then
            pick(i) = pick(i) + 1
            do k = i + 1, m
               pick(k) = pick(k - 1) + 1
            end do
            more = .true.
            return
         end if
      end do
   end function step_combination

   function n_choose(n, k) result(c)
      !! Binomial coefficient, built up rather than from factorials
      integer, intent(in) :: n, k
      integer :: c

      integer :: i

      c = 1
      do i = 1, k
         c = c*(n - k + i)/i
      end do
   end function n_choose

   subroutine subtract_subsets(terms, term_size, n_terms, correction)
      !! Reduce each group's quantity by the differences its subsets already carry
      !!
      !! In on entry: `correction(t)` is `f` evaluated on group `t`. Out on
      !! exit: `correction(t)` is `df_t`, the part of `f` that group and no
      !! smaller group accounts for.
      !!
      !! **The group list has to be closed under taking subsets.** Every proper
      !! subset of a listed group must itself be listed, or its share is never
      !! removed and the difference is not a difference. `enumerate_subsets` is
      !! closed by construction; a caller that filters the list -- EFMO keeps
      !! only the groups whose every internal pair is near -- has to filter with
      !! a criterion that is itself inherited by subsets.
      integer, intent(in) :: terms(:, :)
      integer, intent(in) :: term_size(:)
      integer, intent(in) :: n_terms
      real(dp), intent(inout) :: correction(:)

      integer :: t, u

      do t = 1, n_terms
         if (term_size(t) < 2) cycle
         do u = 1, n_terms
            if (term_size(u) >= term_size(t)) cycle
            if (.not. is_subset(terms(1:term_size(u), u), terms(1:term_size(t), t))) cycle
            correction(t) = correction(t) - correction(u)
         end do
      end do
   end subroutine subtract_subsets

   function is_subset(small, big) result(inside)
      !! Whether every member of `small` appears in `big`
      integer, intent(in) :: small(:), big(:)
      logical :: inside

      integer :: i

      inside = .true.
      do i = 1, size(small)
         if (.not. any(big == small(i))) then
            inside = .false.
            return
         end if
      end do
   end function is_subset

end module mqc_czt_subsets
