!  WORD_LEN    SUBROUTINE  PROGRAMING  C.DCS.LIB   DCS         96.7
!+
! Function WORD_LEN (WORDING)
!
! Function to return the "length" of the first word in a string
! If string(1:1) = ' ' than WORD_LEN is set to 0
!
! Example:
!     ix = word_len ('to be or not    ')
! Result:
!     ix = 2
!-

!$Id$
!$Log$
!Revision 1.4  2004/08/03 15:40:36  dcs
!Fixed standards nonconforming code.
!
!Revision 1.3  2002/02/23 20:34:54  dcs
!Modified for Single/Double Real Toggle
!
!Revision 1.2  2001/09/27 17:47:13  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



  integer function word_len (wording)


  use precision_def

  implicit none

  integer i1, i
  character(*) wording

!

  i1 = len(wording)
  do i = 1, i1
    if (wording(i:i) == ' ' .or. wording(i:i) == char(9)) then
      word_len = i - 1
      return
    endif
  enddo

  word_len = i1
  return

  end
