!  IF_ERROR    SUBROUTINE  PROGRAMING  C.DCS.LIB   DCS         96.7
!+
! subroutine IF_ERROR (IDELIM, ICMD, ERROR_STRING, LINE_NUMBER, END_CHECK)
!
! Subroutine for checking input errors after using the TYPSCN function
! calls ICMTYP, INTTYP, RELTYP, or IOCTYP. On an input error IF_ERROR
! will generate an error message and exit the program. If there is no error
! IF_ERROR will do nothing.
!
! inputs:
!     integer  IDELIM = 2nd argument of ICMTYP, INTTYP, RELTYP, and IOCTYP
!                   functions
!
!     integer ICMD = result of ICMTYP function call. If IF_ERROR is being used
!               to check INTTYP, RELTYP, and IOCTYP this argument should be
!               set to some positive value
!
!     character ERROR_STRING = error string used in an error message.
!
!     integer  LINE_NUMBER = line number of line being parsed. This is uses in
!                   in the output error message. If LINE_NUMBER negative
!                   or zero then LINE_NUMBER is ignored.
!
!     integer END_CHECK = 0 => no error if this is last data on line
!                       = 1 => error if this is last data on line
!
! Example:
!
!     switch = icmtyp (lwait, idelim, switch_names)
!     call if_error (idelim, switch, 'SWITCH', 5, 1)
!
! An error message is generated and the program exited
! if idelim >= 0 or switch =< 0 .
!
! Example:
!
!     real_num = reltyp(lwait, idelim)
!     call if_error (idelim, 1, 'REAL_NUM', -1, 0)
!
!-

!$Id$
!$Log$
!Revision 1.4  2003/04/30 16:26:39  dcs
!F90 Standard conforming update.
!
!Revision 1.3  2002/02/23 20:34:44  dcs
!Modified for Single/Double Real Toggle
!
!Revision 1.2  2001/09/27 17:47:05  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine if_error (idelim, icmd, error_string, line_number, end_check)


  use precision_def

  implicit none

  integer idelim, line_number, end_check, icmd

  character*(*) error_string

! if we have an error...



  if (idelim >= 1 .or. icmd <= 0) then

    print *

    if (line_number > 0) then
      print *, 'ERROR READING: ' // error_string // &
                                             ' IN LINE ', line_number
    else
      print *, 'ERROR READING: ' // error_string
    endif

    call zertyp(' ')
    stop

  endif

! if end of line is an error....

  if (idelim < 0 .and. end_check == 1) then

    print *
    if (line_number > 0) then
      print *, 'ERROR: END OF LINE ENCOUNTERED WHERE NOT EXPECTED AFTER'
      print *, '       READING: ' // error_string // ' IN LINE ', line_number
    else
      print *, 'ERROR: END OF LINE ENCOUNTERED WHERE NOT EXPECTED AFTER'
      print *, '       READING:' // error_string
    endif

    call zertyp (' ')

    stop

  endif

  return

  end

