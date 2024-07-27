!+
! Subroutine detab(str)
!
! Routine to replace tab characters with a single blank spaces.
!
! Input:
!   str         -- character(*): String with possible tabs.
!
! Output:
!   str         -- character(*): String with blank spaces substituted.
!-

subroutine detab(str)

implicit none

character(*) str
character(1), parameter :: tab = achar(9)
integer i

!

do i = 1, len(str)
  if (str(i:i) == tab) str(i:i) = ' '
enddo

end subroutine
