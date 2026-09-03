!! A shared printer for iterative-solver convergence tables
module mqc_convergence_report
   !! Both the SCF (`mqc_czt_rhf`) and the FMO outer loop (`mqc_czt_fmo`)
   !! report convergence as the same table: a title, a ruled column header, one
   !! row per iteration, and a footer saying whether and in how many steps it
   !! converged. The rows differ -- the SCF carries DIIS depth and per-iteration
   !! timings, the FMO loop a single energy sum -- so each caller writes its own
   !! row; the frame around it is identical and lives here.
   !!
   !! Visibility is the caller's, passed as `show`. This module prints through
   !! `logger%info`; the caller sets `show` from whatever logger level and rank
   !! guard the table belongs behind, so the SCF table and the FMO outer table can
   !! sit at different levels without this module knowing about either.
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   implicit none
   private

   public :: convergence_header
   public :: convergence_footer

contains

   subroutine convergence_header(show, title, columns, width)
      !! A blank line, the title, a rule, the column headings and a closing rule
      logical, intent(in) :: show
      character(len=*), intent(in) :: title    !! e.g. "SCF iterations"
      character(len=*), intent(in) :: columns  !! the pre-spaced column headings
      integer, intent(in) :: width             !! rule width, in dashes

      if (.not. show) return
      call logger%info("")
      call logger%info("  "//title)
      call logger%info("  "//repeat("-", width))
      call logger%info("  "//columns)
      call logger%info("  "//repeat("-", width))
   end subroutine convergence_header

   subroutine convergence_footer(show, converged, iterations, noun, width)
      !! A closing rule, then whether it converged and in how many `noun`
      logical, intent(in) :: show, converged
      integer, intent(in) :: iterations
      character(len=*), intent(in) :: noun     !! "iterations", "outer iterations"
      integer, intent(in) :: width

      if (.not. show) return
      call logger%info("  "//repeat("-", width))
      if (converged) then
         call logger%info("  converged in "//to_char(iterations)//" "//noun)
      else
         call logger%info("  NOT converged after "//to_char(iterations)//" "//noun)
      end if
   end subroutine convergence_footer
end module mqc_convergence_report
