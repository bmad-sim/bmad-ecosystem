#include "CESR_platform.inc"

module command_line_mod

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function cesr_iargc ()
!
! Platform independent function to return the number of command
! line arguments. Use this with cesr_getarg.
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   none
!
! Output:
!   cesr_iargc -- Integer: Number of arguments.
!-

function cesr_iargc()

  implicit none

  integer cesr_iargc

#ifndef CESR_VMS
  integer iargc
  cesr_iargc = iargc()
#endif

#if defined (CESR_VMS)
  character(200) :: raw_line
  integer :: n_args, i, pos
  call lib$get_foreign(raw_line)
  cesr_iargc = 0
  do while (len_trim(raw_line) > 0)
     pos = index(raw_line, ' ')
     raw_line = raw_line(pos+1:)
     cesr_iargc = cesr_iargc + 1
  end do
#endif

end function cesr_iargc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine cesr_getarg (i_arg, arg)
!
! Platform independent function to return the i'th command
! line argument. Use this with cesr_iargc.
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   i_arg -- Integer: Index of argument to return.
!              i_arg = 0 => Entire line.
!              i_arg = 1 => First argument.
!
! Output:
!   arg   -- Character(*): i'th command line argument.
!              If i_arg > number_of_args then arg is a blank string. 
!-
  
subroutine cesr_getarg(i_arg, arg)

  implicit none

  integer, intent(in) :: i_arg
  character(*), intent(out) :: arg
  character(200) :: raw
  integer i, pos

!

#ifndef CESR_VMS
  if (i_arg == 0) then
    call getarg (1, arg)
    do i = 2, cesr_iargc()
      call getarg (i, raw)
      arg = trim(arg) // ' ' // trim(raw)
    enddo
  else
    call getarg(i_arg, arg)
  endif
#endif

#if defined (CESR_VMS)
  call lib$get_foreign(raw)
  if (i_arg == 0) then
    arg = raw
  else
    do i = 1, i_arg
       pos = index(raw, ' ')
       arg = raw(:pos-1)
       raw = raw(pos+1:)
    enddo
  endif
#endif

end subroutine cesr_getarg

end module command_line_mod
