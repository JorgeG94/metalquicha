!! ASCII art logo display for metalquicha
module mqc_logo
   !! Provides the project branding sunflower logo and version information
   !! displayed at program startup.
   implicit none
   private

   public :: print_logo  !! Display ASCII sunflower logo and project info

contains

   subroutine print_logo()
      !! Print the PIC Chemistry ASCII sunflower logo

      write (*, '(A)') ' '
      write (*, '(A)') ' '
      write (*, '(A)') '                        __   __'
      write (*, '(A)') '                     .-(  ''.''  )-.'
      write (*, '(A)') '                    (   \  |  /   )'
      write (*, '(A)') '                   ( ''`-.;;;;;.-''` )'
      write (*, '(A)') '                  ( :-==;;;;;;;==-: )'
      write (*, '(A)') '                   (  .-'';;;;;''-.  )'
      write (*, '(A)') '                    (``  /  |  \  ``)'
      write (*, '(A)') '                     ''-(__.''.__)-'''
      write (*, '(A)') ' '
      write (*, '(A)') '                      (Art by jgs)'
      write (*, '(A)') ' '
      write (*, '(A)') '    ╔═══════════════════════════════════════════════╗'
      write (*, '(A)') '    ║              Met"al q"uicha                   ║'
      write (*, '(A)') '    ║                (Sunflower)                    ║'
      write (*, '(A)') '    ║   A hastily put together Fortran code for     ║'
      write (*, '(A)') '    ║     Fragmented Based Quantum Chemistry        ║'
      write (*, '(A)') '    ║                                               ║'
      write (*, '(A)') '    ║        Coded up by Jorge as a hobby           ║'
      write (*, '(A)') '    ╚═══════════════════════════════════════════════╝'
      write (*, '(A)') ' '

   end subroutine print_logo

end module mqc_logo
