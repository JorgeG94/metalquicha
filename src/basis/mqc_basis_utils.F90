!! Utilities for handling basis set names and files
#ifndef MQC_DEFAULT_BASIS_DIR
! The fpm build has no configure step to learn where the source tree is, so
! the compiled-in fallback is simply absent there; $MQC_BASIS_PATH and
! ./basis_sets still work.
#define MQC_DEFAULT_BASIS_DIR ""
#endif
module mqc_basis_utils
   !! Normalizes basis set names and locates the JSON file that defines them.
   !!
   !! A normalized name is exactly the name the Basis Set Exchange uses:
   !! lowercase, with `*` written `_st_` and everything else -- `+`,
   !! parentheses, commas -- left alone. So def2-SV(P) becomes def2-sv(p) and
   !! 6-31G** becomes 6-31g_st__st_.
   !!
   !! Deferring to BSE rather than inventing a scheme is what makes the whole
   !! catalogue addressable. An earlier version flattened harder -- parentheses
   !! and commas dropped, `*` to `s`, `+` to `p` -- which collapsed def2-SV(P)
   !! onto def2-SVP and def2-SV(P)-RIFIT onto def2-SVP-RIFIT: four real basis
   !! sets, two filenames. It also meant the extracted files had to be renamed
   !! out of their bundle names, so the mapping had to be kept in step in two
   !! languages. Now the file on disk carries its BSE name and there is nothing
   !! to keep in step.
   !!
   !! Basis sets are Basis Set Exchange JSON only. The Gaussian94 (`.gbs`) and
   !! GAMESS ($DATA `.txt`) readers were removed: BSE JSON keeps a combined
   !! shell's coefficient sets separate, so SP shells split cleanly with no
   !! `uncontract_spdf=1` download flag, it carries ECP data in the same file,
   !! and one format means one parser to trust.
   use mqc_error, only: error_t, ERROR_IO
   implicit none
   private

   public :: normalize_basis_name
   public :: find_basis_file
   public :: basis_search_path

   integer, parameter :: MAX_PATH = 4096  !! Longest directory list we accept

   character(len=*), parameter :: BASIS_PATH_VARIABLE = "MQC_BASIS_PATH"
      !! Environment variable holding a colon-separated list of directories

contains

   pure function normalize_basis_name(basis_name) result(normalized)
      !! Normalize a basis set name to its Basis Set Exchange spelling
      !!
      !! Rules:
      !!   lowercased
      !!   * -> _st_   (the only character BSE escapes, being awkward in a filename)
      !!   spaces dropped
      !!   everything else kept, `+` and parentheses included
      !!
      !! Examples:
      !!   6-31G*       -> 6-31g_st_
      !!   6-311++G**   -> 6-311++g_st__st_
      !!   6-311G(d,p)  -> 6-311g(d,p)
      !!   def2-SV(P)   -> def2-sv(p)
      !!   def2-SVP     -> def2-svp
      !!   cc-pVDZ      -> cc-pvdz
      character(len=*), intent(in) :: basis_name
      character(len=:), allocatable :: normalized
      integer :: i, out_pos, code
      character(len=256) :: buffer

      buffer = ""
      out_pos = 0

      do i = 1, len_trim(basis_name)
         select case (basis_name(i:i))
         case ("*")
            ! BSE writes the polarization star this way; a literal "*" in a
            ! filename is a glob waiting to go off.
            buffer(out_pos + 1:out_pos + 4) = "_st_"
            out_pos = out_pos + 4

         case (" ")
            continue

         case default
            out_pos = out_pos + 1
            code = iachar(basis_name(i:i))
            if (code >= iachar("A") .and. code <= iachar("Z")) then
               buffer(out_pos:out_pos) = achar(code + 32)
            else
               buffer(out_pos:out_pos) = basis_name(i:i)
            end if
         end select
      end do

      normalized = trim(buffer(1:out_pos))
   end function normalize_basis_name

   function basis_search_path() result(directories)
      !! Directories searched for basis set files, in order
      !!
      !! 1. every entry of $MQC_BASIS_PATH, colon-separated
      !! 2. ./basis_sets, so a run from the source tree works with no setup
      !! 3. the basis_sets/ of the tree this binary was configured from
      !!
      !! Three because each covers a case the others do not: the variable lets
      !! a user point at their own collection, the relative path is what the
      !! test suite and the validation scripts rely on, and the compiled-in
      !! default is what makes `mqc` work from an arbitrary directory.
      character(len=:), allocatable :: directories(:)

      character(len=MAX_PATH) :: env_value
      integer, parameter :: MAX_SEARCH_DIRECTORIES = 64
      character(len=MAX_PATH) :: collected(MAX_SEARCH_DIRECTORIES)
      integer :: env_length, status, n, start, colon

      n = 0

      call get_environment_variable(BASIS_PATH_VARIABLE, env_value, env_length, status)
      if (status == 0 .and. env_length > 0) then
         start = 1
         do while (start <= env_length .and. n < size(collected) - 2)
            colon = index(env_value(start:env_length), ":")
            if (colon == 0) then
               if (len_trim(env_value(start:env_length)) > 0) then
                  n = n + 1
                  collected(n) = adjustl(env_value(start:env_length))
               end if
               exit
            end if
            if (colon > 1) then
               if (len_trim(env_value(start:start + colon - 2)) > 0) then
                  n = n + 1
                  collected(n) = adjustl(env_value(start:start + colon - 2))
               end if
            end if
            start = start + colon
         end do
      end if

      n = n + 1
      collected(n) = "basis_sets"

      if (len(MQC_DEFAULT_BASIS_DIR) > 0) then
         n = n + 1
         collected(n) = MQC_DEFAULT_BASIS_DIR
      end if

      allocate (character(len=MAX_PATH) :: directories(n))
      directories = collected(1:n)
   end function basis_search_path

   subroutine find_basis_file(basis_name, filename, error)
      !! Locate the JSON file defining a basis set
      !!
      !! The name is normalized, then looked for as `<dir>/<normalized>.json`
      !! in each directory of `basis_search_path` until one exists.
      character(len=*), intent(in) :: basis_name
      character(len=:), allocatable, intent(out) :: filename
      type(error_t), intent(out) :: error

      character(len=:), allocatable :: normalized
      character(len=:), allocatable :: directories(:)
      character(len=:), allocatable :: candidate, tried
      logical :: exists
      integer :: idir

      normalized = normalize_basis_name(basis_name)
      if (len(normalized) == 0) then
         call error%set(ERROR_IO, "Empty basis set name")
         return
      end if

      directories = basis_search_path()
      tried = ""

      do idir = 1, size(directories)
         candidate = trim(directories(idir))//"/"//normalized//".json"
         inquire (file=candidate, exist=exists)
         if (exists) then
            filename = candidate
            return
         end if
         if (len(tried) > 0) tried = tried//", "
         tried = tried//candidate
      end do

      call error%set(ERROR_IO, "Basis set file not found for '"//trim(basis_name)// &
                     "'. Tried "//tried//". Add the basis to MQC_BASIS_SETS and "// &
                     "reconfigure, or set "//BASIS_PATH_VARIABLE//" to a directory "// &
                     "holding "//normalized//".json")
   end subroutine find_basis_file

end module mqc_basis_utils
