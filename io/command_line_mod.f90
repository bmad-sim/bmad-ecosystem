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
!   use sim_utils
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

integer iargc
cesr_iargc = iargc()

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
!   use sim_utils
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

if (i_arg == 0) then
  call getarg (1, arg)
  do i = 2, cesr_iargc()
    call getarg (i, raw)
    arg = trim(arg) // ' ' // trim(raw)
  enddo
else
  call getarg(i_arg, arg)
endif

end subroutine cesr_getarg

end module command_line_mod
