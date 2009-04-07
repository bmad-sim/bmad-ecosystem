!  INTEGER_RE  SUBROUTINE  PROGRAMING  C.DCS.LIB   DCS         96.7
!+
! function INTEGER_READ(error_message)
!
! function equivalent to TYPSCN function RELTYP with added automatic error
! message generated
!
! example:
!
!     integer = integer_read('ERROR: ')
!-

!$Id$
!$Log$
!Revision 1.4  2003/04/30 16:26:40  dcs
!F90 Standard conforming update.
!
!Revision 1.3  2002/02/23 20:34:44  dcs
!Modified for Single/Double Real Toggle
!
!Revision 1.2  2001/09/27 17:47:06  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


  integer function integer_read(error_message)


  use precision_def

  implicit none

  integer idelim, inttyp

  logical lwait

  character*(*) error_message

!

  lwait = .false.

  integer_read = inttyp(lwait, idelim)

  if (idelim > 0) then
    print *
    print *, error_message
    call zertyp(' ')
    stop
  endif

  return

  end
