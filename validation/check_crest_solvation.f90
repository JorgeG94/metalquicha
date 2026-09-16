!! Manual check that a solvent named in an mqc deck reaches CREST's sampler
!!
!!     cmake -B build -DMQC_ENABLE_CREST=ON
!!     cmake --build build --target check_crest_solvation
!!     ./build/check_crest_solvation
!!
!! `check_crest_search` established that a whole sampling run can be configured
!! from an argument vector built in memory. This asks a narrower question about
!! that vector: whether a solvation flag placed in it actually reaches the
!! calculator CREST samples with.
!!
!! It is checked on its own because every way it fails is quiet. CREST copies a
!! solvation model onto its own xTB calculator only inside `if (env%gbsa)`
!! (`env2calc`, in `legacy_wrappers.f90`), and nothing but a flag sets
!! `env%gbsa`. A vector without one leaves the metadynamics, the optimisations
!! and the energy window that prunes the ensemble in gas phase, while mqc's
!! refinement level -- which reads the deck and knows about the solvent --
!! re-ranks the survivors in solution. Nothing reports an error: the run prints
!! a solvent from beginning to end and sampled in vacuum, and conformers that
!! are stable only in solution were pruned before refinement ever saw them.
!!
!! Parse only. `parseflags` and the calculator it builds are the whole subject;
!! no search is run, so this costs nothing and needs no compute node.
program check_crest_solvation
   use crest_data, only: systemdata
   implicit none

   !! Explicit interface for a bare external subroutine, as in
   !! `check_crest_search`. `crest_main.f90` calls it through an implicit
   !! interface; declaring it here costs nothing and lets the compiler check
   !! the call.
   interface
      subroutine parseflags(env, arg, nra)
         import :: systemdata
         implicit none
         type(systemdata), intent(inout) :: env
         integer, intent(in) :: nra
         !! No intent, deliberately: the definition in confparse.f90 declares
         !! it without one, and intent is part of a procedure's
         !! characteristics -- adding it here would make this interface
         !! disagree with the procedure it describes.
         ! allow(missing-intent)
         character(len=*) :: arg(nra)
      end subroutine parseflags
   end interface

   integer, parameter :: ARG_LEN = 64
   character(len=*), parameter :: PROBE_FILE = "crest_probe.xyz"

   type(systemdata) :: env_solv, env_gas
   character(len=:), allocatable :: argv(:)
   integer :: unit
   logical :: ok

   ok = .true.

   ! The geometry parseflags insists on reading. Which molecule it is does not
   ! matter here -- nothing is computed -- only that the calculator gets built.
   open (newunit=unit, file=PROBE_FILE, status="replace", action="write")
   write (unit, "(a)") "3"
   write (unit, "(a)") "water"
   write (unit, "(a)") "O   0.0000000   0.0000000   0.1173000"
   write (unit, "(a)") "H   0.0000000   0.7572000  -0.4692000"
   write (unit, "(a)") "H   0.0000000  -0.7572000  -0.4692000"
   close (unit)

   ! The vector mqc builds for a deck that named `keywords.xtb.solvent`, with
   ! the flag and its name in the last two places exactly as the driver puts
   ! them.
   allocate (character(len=ARG_LEN) :: argv(9))
   argv(1) = PROBE_FILE
   argv(2) = "-T"
   argv(3) = "1"
   argv(4) = "-chrg"
   argv(5) = "0"
   argv(6) = "-uhf"
   argv(7) = "0"
   argv(8) = "-alpb"
   argv(9) = "water"
   call parseflags(env_solv, argv, 9)
   deallocate (argv)

   ! `-chrg` and `-uhf` leave these behind, and CREST refuses to start in a
   ! directory holding either without the matching flag.
   call remove_file(".CHRG")
   call remove_file(".UHF")

   write (*, "(a)") ""
   write (*, "(a)") "  with -alpb water:"
   write (*, "(a,l1)") "    env%gbsa             ", env_solv%gbsa
   write (*, "(a,a)") "    env%solv             ", shown(env_solv%solv)
   write (*, "(a,a)") "    env%solvent          ", trim(env_solv%solvent)
   write (*, "(a,i0)") "    calculation levels   ", env_solv%calc%ncalculations

   if (.not. env_solv%gbsa) then
      write (*, "(a)") "  FAIL: the flag did not set env%gbsa, so env2calc will skip solvation"
      ok = .false.
   end if

   if (env_solv%calc%ncalculations < 1) then
      write (*, "(a)") "  FAIL: no calculation level was built to carry the solvation"
      ok = .false.
   else
      write (*, "(a,a)") "    level 1 solvmodel    ", shown(env_solv%calc%calcs(1)%solvmodel)
      write (*, "(a,a)") "    level 1 solvent      ", shown(env_solv%calc%calcs(1)%solvent)

      ! This is the component `api_engrad` hands to `tblite_add_solv`, so it is
      ! the one that decides whether a sampling gradient is solvated.
      if (.not. allocated(env_solv%calc%calcs(1)%solvmodel)) then
         write (*, "(a)") "  FAIL: the sampling calculator has no solvation model"
         ok = .false.
      else if (trim(env_solv%calc%calcs(1)%solvmodel) /= "alpb") then
         write (*, "(a)") "  FAIL: -alpb did not select ALPB on the sampling calculator"
         ok = .false.
      end if

      if (.not. allocated(env_solv%calc%calcs(1)%solvent)) then
         write (*, "(a)") "  FAIL: the sampling calculator has no solvent name"
         ok = .false.
      else if (trim(env_solv%calc%calcs(1)%solvent) /= "water") then
         write (*, "(a)") "  FAIL: the solvent name did not survive the parse"
         ok = .false.
      end if
   end if

   ! The control, and the reason this check exists: the same vector without the
   ! last two entries has to come back gas phase. If it did not, the flag would
   ! be proving nothing.
   allocate (character(len=ARG_LEN) :: argv(7))
   argv(1) = PROBE_FILE
   argv(2) = "-T"
   argv(3) = "1"
   argv(4) = "-chrg"
   argv(5) = "0"
   argv(6) = "-uhf"
   argv(7) = "0"
   call parseflags(env_gas, argv, 7)
   deallocate (argv)

   call remove_file(".CHRG")
   call remove_file(".UHF")

   write (*, "(a)") ""
   write (*, "(a)") "  without the flag:"
   write (*, "(a,l1)") "    env%gbsa             ", env_gas%gbsa
   write (*, "(a,i0)") "    calculation levels   ", env_gas%calc%ncalculations

   if (env_gas%gbsa) then
      write (*, "(a)") "  FAIL: env%gbsa was set without a flag asking for it"
      ok = .false.
   end if

   if (env_gas%calc%ncalculations >= 1) then
      write (*, "(a,a)") "    level 1 solvmodel    ", shown(env_gas%calc%calcs(1)%solvmodel)
      if (allocated(env_gas%calc%calcs(1)%solvmodel)) then
         write (*, "(a)") "  FAIL: the unflagged vector built a solvated calculator"
         ok = .false.
      end if
   end if

   write (*, "(a)") ""
   if (.not. ok) then
      write (*, "(a)") "  FAILED"
      error stop 1
   end if
   write (*, "(a)") "  OK: the solvation flag reaches the calculator CREST samples with,"
   write (*, "(a)") "      and its absence leaves that calculator in gas phase"

contains

   function shown(text) result(display)
      !! A deferred-length component's value, or a marker when it has none.
      !! Several of these are unallocated in a gas-phase parse, which is the
      !! result being reported rather than something to guard against.
      character(len=:), allocatable, intent(in) :: text
      character(len=:), allocatable :: display

      if (allocated(text)) then
         display = trim(text)
      else
         display = "(unallocated)"
      end if
   end function shown

   subroutine remove_file(path)
      !! Delete a file if it is there, and say nothing if it is not
      character(len=*), intent(in) :: path

      integer :: unit, io
      logical :: exists

      inquire (file=path, exist=exists)
      if (.not. exists) return
      ! `readwrite` rather than `read`: the unit exists only so that the close
      ! can remove the file, which is not a read.
      open (newunit=unit, file=path, status="old", action="readwrite", iostat=io)
      if (io == 0) close (unit, status="delete")
   end subroutine remove_file

end program check_crest_solvation
