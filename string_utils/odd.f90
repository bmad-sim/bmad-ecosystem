!  ODD         FUNCTION    LOGICAL     C.DCS.LIB   DCS         98.3
!+
! Function ODD(NUM)
!
! Function to return true if NUM is odd
!
! Input:
!     NUM -- Integer. Number to test
!
! Output:
!     ODD -- Logical. True if NUM is odd.
!-

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:34:48  dcs
!Modified for Single/Double Real Toggle
!
!Revision 1.2  2001/09/27 17:47:09  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


  function odd (num)


  use precision_def

  implicit none

  integer num

  logical odd

!

  if (mod(num, 2) == 0) then
    odd = .false.
  else
    odd = .true.
  endif

  return

  end

