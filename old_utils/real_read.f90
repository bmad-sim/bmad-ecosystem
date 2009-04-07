!  REAL_READ   SUBROUTINE  PROGRAMING  C.DCS.LIB   DCS         96.7
!+
! function REAL_READ(error_message)
!
! function equivalent to TYPSCN function RELTYP with added automatic error
! message generated
!
! example:
!
!     real(rp) = real_read('ERROR: ')
!-

!$Id$
!$Log$
!Revision 1.5  2003/07/09 01:29:31  dcs
!new bmad
!
!Revision 1.4  2003/04/30 16:26:42  dcs
!F90 Standard conforming update.
!
!Revision 1.3  2002/02/23 20:34:51  dcs
!Modified for Single/Double Real Toggle
!
!Revision 1.2  2001/09/27 17:47:11  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


  real(rp) function real_read(error_message)


  use precision_def

  implicit none

  integer idelim

  real(rp) reltyp

  logical lwait

  character*(*) error_message

!

  lwait = .false.

  real_read = reltyp(lwait, idelim)

  if (idelim > 0) then
    print *
    print *, error_message
    call zertyp(' ')
    stop
  endif

  return

  end
