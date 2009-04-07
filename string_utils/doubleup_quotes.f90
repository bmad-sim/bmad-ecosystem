!+
! Subroutine doubleup_quotes (str_in, str_out, quote)
!
! Subroutine to doubleup on quotation marks in a string
!
! Input:
!   str_in -- Character*(*): Input string
!   quote  -- Character*1: quotation mark to double up (" or ')
!
! Output:
!   str_out -- Character*(*): Output string
!
! Example:
!   call doubleup_quotes ('This:":', str_out, '"')
!   str_out -> 'This:"":'
!-                               

!$Id$
!$Log$
!Revision 1.5  2003/04/30 16:26:38  dcs
!F90 Standard conforming update.
!
!Revision 1.4  2002/02/23 20:34:41  dcs
!Modified for Single/Double Real Toggle
!
!Revision 1.3  2001/10/25 19:13:14  helms
!Added explicit variable declarations.
!mat_inv has been replaced by mat_inverse in BMAD (not DCSLIB)
!
!Revision 1.2  2001/09/27 17:47:03  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine doubleup_quotes (str_in, str_out, quote)

  use precision_def

  implicit none

  character*(*) str_in, str_out
  character quote*1, str*500

  integer ix_in, ix_out, ixq

!

  ix_in = 1
  ix_out = 1

  do
    ixq = index(str_in(ix_in:), quote)
    if (ixq == 0) then
      str(ix_out:) = str_in(ix_in:)
      if (len(str_out) .lt. len_trim(str)) then
        print *, 'WARNING FROM DOUBLEUP_QUOTES: OUTPUT STRING IS TOO SHORT!'
      endif
      str_out = str
      return
    else
      str(ix_out:) = str_in(ix_in:ix_in+ixq-1) // quote
      ix_in = ix_in + ixq
      ix_out = ix_out + ixq + 1
    endif
  enddo

end subroutine
