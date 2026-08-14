!! Program limits and default parameter, publics
module mqc_program_limits
   !! Contains compile-time limits and default values for the metalquicha program.
   !! These are tunable parameter that control memory allocation and algorithm behavior.
   use pic_types, only: dp
   implicit none
   private

   !---------------------------------------------------------------------------
   ! Many-Body Expansion Limits
   !---------------------------------------------------------------------------

   !> Maximum MBE truncation order (1-body, 2-body, ..., N-body)
   !> Higher orders require factorial growth in fragment combinations
   integer, parameter, public :: MAX_MBE_LEVEL = 10

   !> Terms in an EFP2 interaction energy: electrostatics, polarization,
   !> exchange repulsion, dispersion, charge transfer, and their total. Here
   !> rather than beside the EFP code because the real CPU bridge and the stub
   !> that stands in for it must agree, and only one of those is compiled with
   !> the backend present.
   integer, parameter, public :: N_EFP_TERMS = 6
   !> Terms a SAPT0 interaction energy is reported as: electrostatics, exchange,
   !> induction, exchange-induction, dispersion, exchange-dispersion, the
   !> delta-HF correction, the supermolecular HF reference, and the total. Here
   !> rather than beside the SAPT code because the real CPU bridge and the stub
   !> standing in for it must agree, and only one is compiled with the backend.
   integer, parameter, public :: N_SAPT_TERMS = 12

   !> The slot each SAPT term occupies, and the name it is written out under.
   !> Both sides of the backend gate index this array, and so does the JSON
   !> writer, so the order is fixed in exactly one place. The names are the
   !> literature's rather than the console's prose, because their reason to
   !> exist is being compared against another code's output.
   !>
   !> `_u` is uncoupled and `_r` is response (coupled); `_s2` is the single
   !> exchange approximation against the full `S^inf`. A code that reports only
   !> one of each pair cannot be checked against one that reports the other.
   character(len=*), parameter, public :: SAPT_TERM_NAMES(N_SAPT_TERMS) = &
                                          [character(len=12) :: "elst10", "exch10_s2", "exch10", &
                                                                 "ind20_u", "ind20_r", "exch_ind20_u", "exch_ind20_r", &
                                                              "disp20", "exch_disp20", "delta_hf", "e_int_hf_cp", "total"]

   !> Group-global result batching size for MPI MBE (multi-global coordinator)
   integer, parameter, public :: GROUP_RESULT_BATCH_SIZE = 256

   !---------------------------------------------------------------------------
   ! Numerical Differentiation Defaults
   !---------------------------------------------------------------------------

   !> Default step size for finite difference calculations (Bohr)
   !> ~0.005 Bohr = ~0.0026 Angstrom, suitable for Hessian/gradient FD
   real(dp), parameter, public :: DEFAULT_FD_DISPLACEMENT = 0.005_dp

   !---------------------------------------------------------------------------
   ! I/O Limits
   !---------------------------------------------------------------------------

   !> Maximum length for input file lines
   integer, parameter, public :: MAX_LINE_LENGTH = 1024

   !> Maximum length for element symbols (e.g., "He", "Uue")
   integer, parameter, public :: MAX_ELEMENT_SYMBOL_LEN = 4

   !> JSON output format for real numbers (scientific notation)
   !> Valid values: 'G', 'E', 'EN', 'ES' (json-fortran uses machine precision)
   character(len=*), parameter, public :: JSON_REAL_FORMAT = "ES"

   !---------------------------------------------------------------------------
   ! Geometry/Structure Limits
   !---------------------------------------------------------------------------

   !> Minimum allowed distance between atoms (Bohr)
   !> Atoms closer than this are considered overlapping (error condition)
   real(dp), parameter, public :: MIN_ATOM_DISTANCE = 0.01_dp

end module mqc_program_limits
