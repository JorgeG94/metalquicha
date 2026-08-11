!! Schema validation for MQC JSON input decks
module mqc_json_schema
   !! Checks a deck against the input schema before anything reads it.
   !!
   !! Why this is separate from the reader. The reader's `optional_*`
   !! accessors leave a component alone when a key is absent, which is what
   !! keeps the defaults declared in one place -- but it means a key the reader
   !! does not recognise is indistinguishable from a key nobody wrote.
   !! `"maxiter"` misspelled as `"max_iter"` would silently run 300 iterations
   !! instead of the 40 you asked for, and nothing would say so. Rejecting
   !! unknown keys is the only way to tell those two apart, and it cannot be
   !! done by the reader, which by construction only ever looks up names it
   !! already knows.
   !!
   !! So this module walks the document instead of querying it. That needs the
   !! `json_core`/`json_value` interface rather than `json_file`: the latter
   !! reaches a value by path and has no way to ask an object what keys it
   !! actually has.
   !!
   !! What is checked:
   !!
   !!   * every key, at every level, is one the schema defines
   !!   * required keys are present
   !!   * the parallel fragment lists have matching lengths
   !!   * fragment charges sum to the molecular charge
   !!
   !! What is not: types and ranges. The reader reports those where it reads
   !! them, with the value in hand, and duplicating the knowledge here would
   !! give two places to update and two chances to disagree.
   use mqc_error, only: error_t, ERROR_VALIDATION
   use json_module, only: json_file, json_core, json_value
   implicit none
   private

   public :: ensure_valid_json  !! Reject a deck the schema does not describe

   integer, parameter :: MAX_KEYS = 16
      !! Widest allow-list below; raising it is the only cost of a new key

   type :: key_set_t
      !! The keys one object may contain, and which of them it must
      character(len=32) :: allowed(MAX_KEYS) = ""
      integer :: n_allowed = 0
      character(len=32) :: required(MAX_KEYS) = ""
      integer :: n_required = 0
   end type key_set_t

contains

   subroutine ensure_valid_json(json, error, settings_only)
      !! Validate a loaded deck, setting `error` on the first problem found
      !!
      !! Reports one problem rather than all of them: the first unknown key is
      !! usually the whole story, and a reader that stops there gives a
      !! message short enough to act on.
      !!
      !! `settings_only` validates a document that describes what to compute
      !! but not what to compute it on -- the form that crosses to the workers
      !! when the system arrives by another route. It relaxes exactly one
      !! thing, the requirement that molecules be present; every other key is
      !! held to the same spelling and the same shape, which is the point of
      !! sending a document rather than a struct.
      type(json_file), intent(inout) :: json
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: settings_only

      type(json_core) :: core
      type(json_value), pointer :: root
      logical :: settings

      settings = .false.
      if (present(settings_only)) settings = settings_only

      call json%get_core(core)
      call json%get(root)
      if (.not. associated(root)) then
         call error%set(ERROR_VALIDATION, "Input file contains no JSON document")
         return
      end if

      call check_object(core, root, "", root_keys(settings), error)
      if (error%has_error()) return

      call check_child_object(core, root, "schema", schema_keys(), error)
      if (error%has_error()) return
      call check_child_object(core, root, "model", model_keys(), error)
      if (error%has_error()) return
      call check_child_object(core, root, "system", system_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "system", "logger", logger_keys(), error)
      if (error%has_error()) return

      call check_child_object(core, root, "keywords", keywords_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "scf", scf_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "correlation", &
                                   correlation_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "hessian", hessian_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "aimd", aimd_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "xtb", xtb_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "fragmentation", &
                                   fragmentation_keys(), error)
      if (error%has_error()) return
      call check_cutoffs(core, root, error)
      if (error%has_error()) return

      call check_molecules(core, root, error)
   end subroutine ensure_valid_json

   ! ==========================================================================
   !  The schema
   !
   !  One function per object, so adding a key is a one-line change next to
   !  the other keys of the same object rather than an edit to a table
   !  somewhere else.
   ! ==========================================================================

   function root_keys(settings_only) result(keys)
      logical, intent(in) :: settings_only
      type(key_set_t) :: keys
      call allow(keys, "schema")
      call allow(keys, "model")
      call allow(keys, "driver")
      call allow(keys, "molecules")
      call allow(keys, "keywords")
      call allow(keys, "system")
      call allow(keys, "title")
      call require(keys, "schema")
      call require(keys, "model")
      call require(keys, "driver")
      ! Still *allowed* when settings-only, so a document that carries one gets
      ! the reader's explanation of why it will not be used rather than a bare
      ! "unknown key: molecules", which would read as the opposite problem.
      if (.not. settings_only) call require(keys, "molecules")
   end function root_keys

   function schema_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "name")
      call allow(keys, "version")
      call allow(keys, "index_base")
      call allow(keys, "units")
      call require(keys, "name")
      call require(keys, "version")
   end function schema_keys

   function model_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "method")
      call allow(keys, "basis")
      call allow(keys, "aux_basis")
      call allow(keys, "functional")
      call require(keys, "method")
   end function model_keys

   function system_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "logger")
      call allow(keys, "skip_json_output")
      call allow(keys, "unchecked_input")
      call allow(keys, "checkpoint")
      call allow(keys, "fragment_breakdown")
   end function system_keys

   function logger_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "level")
   end function logger_keys

   function keywords_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "scf")
      call allow(keys, "hessian")
      call allow(keys, "aimd")
      call allow(keys, "fragmentation")
      call allow(keys, "xtb")
      call allow(keys, "correlation")
   end function keywords_keys

   function scf_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "maxiter")
      call allow(keys, "tolerance")
      call allow(keys, "unrestricted")
      call allow(keys, "guess")
      call allow(keys, "allow_crap_scf")
      call allow(keys, "density_fitting")
   end function scf_keys

   function correlation_keys() result(keys)
      !! Post-Hartree-Fock settings, deliberately not under "scf"
      !!
      !! The reference and the correlation treatment are configured apart
      !! because they can disagree: a density-fitted Hartree-Fock followed by a
      !! conventional MP2 is a combination someone will ask for, and one
      !! `density_fitting` shared between them could not express it. RI-MP2
      !! adds its own here rather than reading the SCF's.
      type(key_set_t) :: keys
      call allow(keys, "freeze_core")
      call allow(keys, "n_frozen_core")
      call allow(keys, "density_fitting")
      call allow(keys, "aux_basis")
      call allow(keys, "scs")
      call allow(keys, "scs_ss")
      call allow(keys, "scs_os")
   end function correlation_keys

   function hessian_keys() result(keys)
      type(key_set_t) :: keys
      ! Both spellings of the displacement are accepted, as the retired JSON
      ! generator accepted them.
      call allow(keys, "finite_difference_displacement")
      call allow(keys, "displacement")
      call allow(keys, "temperature")
      call allow(keys, "pressure")
   end function hessian_keys

   function aimd_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "dt")
      call allow(keys, "timestep")
      call allow(keys, "nsteps")
      call allow(keys, "steps")
      call allow(keys, "initial_temperature")
      call allow(keys, "temperature")
      call allow(keys, "output_frequency")
      call allow(keys, "output_freq")
   end function aimd_keys

   function xtb_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "solvent")
      call allow(keys, "solvation_model")
      call allow(keys, "dielectric")
      call allow(keys, "cpcm_nang")
      call allow(keys, "cpcm_rscale")
      call allow(keys, "use_cds")
      call allow(keys, "use_shift")
   end function xtb_keys

   function fragmentation_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "method")
      call allow(keys, "level")
      call allow(keys, "allow_overlapping_fragments")
      call allow(keys, "max_intersection_level")
      call allow(keys, "embedding")
      call allow(keys, "cutoff_method")
      call allow(keys, "distance_metric")
      call allow(keys, "cutoffs")
      call allow(keys, "global_groups")
      call allow(keys, "nodes_per_group")
      call require(keys, "method")
   end function fragmentation_keys

   function molecule_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "symbols")
      call allow(keys, "geometry")
      call allow(keys, "xyz")
      call allow(keys, "molecular_charge")
      call allow(keys, "molecular_multiplicity")
      call allow(keys, "connectivity")
      call allow(keys, "fragments")
      call allow(keys, "fragment_charges")
      call allow(keys, "fragment_multiplicities")
      call require(keys, "molecular_charge")
      call require(keys, "molecular_multiplicity")
   end function molecule_keys

   ! ==========================================================================
   !  Walking
   ! ==========================================================================

   subroutine check_object(core, object, path, keys, error)
      !! Every key of `object` must be in `keys`, and every required one present
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: object
      character(len=*), intent(in) :: path  !! For the message; "" at the root
      type(key_set_t), intent(in) :: keys
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: child
      character(len=:), allocatable :: name
      integer :: n_children, ichild, ikey
      logical :: found

      call core%info(object, n_children=n_children)

      do ichild = 1, n_children
         call core%get_child(object, ichild, child, found)
         if (.not. found) cycle
         call core%info(child, name=name)
         if (.not. allocated(name)) cycle
         if (.not. is_allowed(keys, name)) then
            call error%set(ERROR_VALIDATION, "Unknown key "//quoted(path, name)// &
                           ". Allowed here: "//allowed_list(keys))
            return
         end if
      end do

      do ikey = 1, keys%n_required
         call core%get(object, trim(keys%required(ikey)), child, found)
         if (.not. found .or. .not. associated(child)) then
            call error%set(ERROR_VALIDATION, "Missing required key "// &
                           quoted(path, trim(keys%required(ikey))))
            return
         end if
      end do
   end subroutine check_object

   subroutine check_child_object(core, parent, name, keys, error)
      !! Validate `parent.name` if it is there; absent is fine
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: parent
      character(len=*), intent(in) :: name
      type(key_set_t), intent(in) :: keys
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: child
      logical :: found

      call core%get(parent, name, child, found)
      if (.not. found .or. .not. associated(child)) return
      call check_object(core, child, name, keys, error)
   end subroutine check_child_object

   subroutine check_grandchild_object(core, root, parent_name, name, keys, error)
      !! Validate `parent_name.name` if both levels are there
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      character(len=*), intent(in) :: parent_name, name
      type(key_set_t), intent(in) :: keys
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: parent, child
      logical :: found

      call core%get(root, parent_name, parent, found)
      if (.not. found .or. .not. associated(parent)) return
      call core%get(parent, name, child, found)
      if (.not. found .or. .not. associated(child)) return
      call check_object(core, child, parent_name//"."//name, keys, error)
   end subroutine check_grandchild_object

   subroutine check_cutoffs(core, root, error)
      !! Cutoff keys are n-mer names or decimal levels, so they get their own rule
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: keywords, fragmentation, cutoffs, child
      character(len=:), allocatable :: name
      integer :: n_children, ichild
      logical :: found

      call core%get(root, "keywords", keywords, found)
      if (.not. found .or. .not. associated(keywords)) return
      call core%get(keywords, "fragmentation", fragmentation, found)
      if (.not. found .or. .not. associated(fragmentation)) return
      call core%get(fragmentation, "cutoffs", cutoffs, found)
      if (.not. found .or. .not. associated(cutoffs)) return

      call core%info(cutoffs, n_children=n_children)
      do ichild = 1, n_children
         call core%get_child(cutoffs, ichild, child, found)
         if (.not. found) cycle
         call core%info(child, name=name)
         if (.not. allocated(name)) cycle
         if (nmer_level(name) <= 0) then
            call error%set(ERROR_VALIDATION, "Unknown key "// &
                           quoted("keywords.fragmentation.cutoffs", name)// &
                           ". Use an n-mer name (dimer, trimer, tetramer, "// &
                           "pentamer, hexamer, heptamer, octamer) or a level "// &
                           "number of 2 or more")
            return
         end if
      end do
   end subroutine check_cutoffs

   subroutine check_molecules(core, root, error)
      !! Every molecule's keys, plus the cross-checks between its lists
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: molecules, molecule
      integer :: n_molecules, imol
      logical :: found

      call core%get(root, "molecules", molecules, found)
      if (.not. found .or. .not. associated(molecules)) return

      call core%info(molecules, n_children=n_molecules)
      if (n_molecules <= 0) then
         call error%set(ERROR_VALIDATION, "'molecules' must contain at least one entry")
         return
      end if

      do imol = 1, n_molecules
         call core%get_child(molecules, imol, molecule, found)
         if (.not. found) cycle
         call check_object(core, molecule, molecule_path(imol), molecule_keys(), error)
         if (error%has_error()) return
         call check_molecule_geometry(core, molecule, imol, error)
         if (error%has_error()) return
         call check_molecule_fragments(core, molecule, imol, error)
         if (error%has_error()) return
      end do
   end subroutine check_molecules

   subroutine check_molecule_geometry(core, molecule, imol, error)
      !! Exactly one of `xyz` or `symbols` + `geometry`
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: molecule
      integer, intent(in) :: imol
      type(error_t), intent(inout) :: error

      logical :: has_xyz, has_symbols, has_geometry

      has_xyz = has_key(core, molecule, "xyz")
      has_symbols = has_key(core, molecule, "symbols")
      has_geometry = has_key(core, molecule, "geometry")

      if (has_xyz .and. (has_symbols .or. has_geometry)) then
         call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                        ": give either 'xyz' or 'symbols' with 'geometry', not both")
         return
      end if
      if (.not. has_xyz) then
         if (.not. (has_symbols .and. has_geometry)) then
            call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                           ": must give either 'xyz' or 'symbols' with 'geometry'")
            return
         end if
      end if
   end subroutine check_molecule_geometry

   subroutine check_molecule_fragments(core, molecule, imol, error)
      !! Parallel list lengths, and charges that add up
      !!
      !! Fragment charges summing to something other than the molecular charge
      !! is the one of these that produces a plausible wrong answer rather than
      !! a crash, which is why it is worth catching here.
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: molecule
      integer, intent(in) :: imol
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: fragments, charges
      integer :: n_fragments, n_charges, n_mults, total, molecular_charge
      logical :: found

      call core%get(molecule, "fragments", fragments, found)
      if (.not. found .or. .not. associated(fragments)) return
      call core%info(fragments, n_children=n_fragments)
      if (n_fragments <= 0) then
         call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                        ": 'fragments' must not be empty if it is given")
         return
      end if

      n_charges = child_count(core, molecule, "fragment_charges")
      if (n_charges > 0 .and. n_charges /= n_fragments) then
         call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                        ": 'fragment_charges' has "//int_text(n_charges)// &
                        " entries for "//int_text(n_fragments)//" fragments")
         return
      end if

      n_mults = child_count(core, molecule, "fragment_multiplicities")
      if (n_mults > 0 .and. n_mults /= n_fragments) then
         call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                        ": 'fragment_multiplicities' has "//int_text(n_mults)// &
                        " entries for "//int_text(n_fragments)//" fragments")
         return
      end if

      if (n_charges <= 0) return
      call core%get(molecule, "fragment_charges", charges, found)
      if (.not. found) return
      total = sum_integer_array(core, charges, n_charges)

      call core%get(molecule, "molecular_charge", molecular_charge, found)
      if (.not. found) return
      if (total /= molecular_charge) then
         call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                        ": fragment charges sum to "//int_text(total)// &
                        " but the molecular charge is "//int_text(molecular_charge)// &
                        "; these must agree")
      end if
   end subroutine check_molecule_fragments

   ! ==========================================================================
   !  Helpers
   ! ==========================================================================

   pure subroutine allow(keys, name)
      !! Record a key an object may contain
      type(key_set_t), intent(inout) :: keys
      character(len=*), intent(in) :: name

      keys%n_allowed = keys%n_allowed + 1
      keys%allowed(keys%n_allowed) = name
   end subroutine allow

   pure subroutine require(keys, name)
      !! Record a key an object must contain
      !!
      !! Required keys are also allowed ones, so every call here is paired with
      !! an `allow` above rather than implying it -- the allow-list is what the
      !! error message offers the user, and it should read as the whole set.
      type(key_set_t), intent(inout) :: keys
      character(len=*), intent(in) :: name

      keys%n_required = keys%n_required + 1
      keys%required(keys%n_required) = name
   end subroutine require

   pure function is_allowed(keys, name) result(ok)
      type(key_set_t), intent(in) :: keys
      character(len=*), intent(in) :: name
      logical :: ok
      integer :: i

      ok = .false.
      do i = 1, keys%n_allowed
         if (trim(keys%allowed(i)) == trim(name)) then
            ok = .true.
            return
         end if
      end do
   end function is_allowed

   pure function allowed_list(keys) result(text)
      !! The allowed keys, for an error message that says what to write instead
      type(key_set_t), intent(in) :: keys
      character(len=:), allocatable :: text
      integer :: i

      text = ""
      do i = 1, keys%n_allowed
         if (i > 1) text = text//", "
         text = text//trim(keys%allowed(i))
      end do
   end function allowed_list

   pure function quoted(path, name) result(text)
      !! A key named as the user would write it, with its parent path
      character(len=*), intent(in) :: path, name
      character(len=:), allocatable :: text

      if (len_trim(path) == 0) then
         text = "'"//trim(name)//"'"
      else
         text = "'"//trim(path)//"."//trim(name)//"'"
      end if
   end function quoted

   pure function molecule_path(imol) result(text)
      character(len=*), parameter :: PREFIX = "molecules["
      integer, intent(in) :: imol
      character(len=:), allocatable :: text

      ! Written 0-based to match how a user counts entries in the array.
      text = PREFIX//int_text(imol - 1)//"]"
   end function molecule_path

   function has_key(core, object, name) result(present_key)
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: object
      character(len=*), intent(in) :: name
      logical :: present_key

      type(json_value), pointer :: child

      call core%get(object, name, child, present_key)
      present_key = present_key .and. associated(child)
   end function has_key

   function child_count(core, object, name) result(n)
      !! Entries in `object.name`, or 0 when it is absent
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: object
      character(len=*), intent(in) :: name
      integer :: n

      type(json_value), pointer :: child
      logical :: found

      n = 0
      call core%get(object, name, child, found)
      if (.not. found .or. .not. associated(child)) return
      call core%info(child, n_children=n)
   end function child_count

   function sum_integer_array(core, array, n) result(total)
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: array
      integer, intent(in) :: n
      integer :: total

      type(json_value), pointer :: element
      integer :: i, value
      logical :: found

      total = 0
      do i = 1, n
         call core%get_child(array, i, element, found)
         if (.not. found) cycle
         call core%get(element, value)
         total = total + value
      end do
   end function sum_integer_array

   pure function nmer_level(name) result(level)
      !! The n-mer level a cutoff key names, or 0 if it names none
      character(len=*), intent(in) :: name
      integer :: level

      integer :: status

      select case (trim(name))
      case ("dimer")
         level = 2
      case ("trimer")
         level = 3
      case ("tetramer")
         level = 4
      case ("pentamer")
         level = 5
      case ("hexamer")
         level = 6
      case ("heptamer")
         level = 7
      case ("octamer")
         level = 8
      case default
         read (name, *, iostat=status) level
         if (status /= 0 .or. level < 2) level = 0
      end select
   end function nmer_level

   pure function int_text(value) result(text)
      integer, intent(in) :: value
      character(len=:), allocatable :: text
      character(len=16) :: buffer

      write (buffer, "(I0)") value
      text = trim(adjustl(buffer))
   end function int_text

end module mqc_json_schema
