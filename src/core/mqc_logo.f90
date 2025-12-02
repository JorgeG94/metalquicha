module mqc_logo
   !! ASCII art logo display for PIC Chemistry
   !!
   !! Provides the project branding sunflower logo and version information
   !! displayed at program startup.
   implicit none
   private

   public :: print_logo  !! Display ASCII sunflower logo and project info

contains

   subroutine print_logo()
      !! Print the PIC Chemistry ASCII sunflower logo

      write (*, '(A)') ''
      write (*, '(A)') ''
      write (*, '(A)') '                        __   __'
      write (*, '(A)') '                     .-(  ''.''  )-.'
      write (*, '(A)') '                    (   \  |  /   )'
      write (*, '(A)') '                   ( ''`-.;;;;;.-''` )'
      write (*, '(A)') '                  ( :-==;;;;;;;==-: )'
      write (*, '(A)') '                   (  .-'';;;;;''-.  )'
      write (*, '(A)') '                    (``  /  |  \  ``)'
      write (*, '(A)') '                     ''-(__.''.__)-'''
      write (*, '(A)') ''
      write (*, '(A)') '                      (Art by jgs)'
      write (*, '(A)') ''
      write (*, '(A)') '    ╔═══════════════════════════════════════════════╗'
      write (*, '(A)') '    ║              Met"al q"uicha                   ║'
      write (*, '(A)') '    ║                (Sunflower)                    ║'
      write (*, '(A)') '    ║     A hastily put together framework for      ║'
      write (*, '(A)') '    ║   Fortran Based High Performance Computing    ║'
      write (*, '(A)') '    ║                                               ║'
      write (*, '(A)') '    ║        Case study: quantum chemistry          ║'
      write (*, '(A)') '    ╚═══════════════════════════════════════════════╝'
      write (*, '(A)') ''

   end subroutine print_logo

end module mqc_logo
