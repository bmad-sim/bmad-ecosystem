!+
! Subroutine BACK_TO_UNDER (OUTPUT, INPUT)
!
!   Subroutine to convert backslashes to underscores.  The purpose of this
! routine is to prevent the creation of filenames that contain backslashes.
! -- Created by Daniel Fromowitz, September 1999.
!
! Input:
!     INPUT  -- Character: Name that may contain backslashes
!
! Output:
!     OUTPUT -- Character: Name that has any backslashes changed to underscores
!-

subroutine back_to_under (output, input)
  implicit none

  character*(*) input, output
  integer i,l

  l = len(input)
  do i = 1,l
    if (input(i:i) == '\') then
      output(i:i) = '_'
    else
      output(i:i) = input(i:i)
    endif
  enddo

  return
  end
