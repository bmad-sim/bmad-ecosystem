module command_line_mod

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function cesr_iargc ()
!
! Note: Use the Fortran intrinsic command_argument_count instead
!
! Platform independent function to return the number of command line arguments. 
! Use this with cesr_getarg.
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

cesr_iargc = command_argument_count()

end function cesr_iargc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine cesr_getarg (i_arg, arg)
!
! Platform independent function to return the i'th command line argument. 
! Use this with cesr_iargc.
!
! Note: The difference between this routine and the Fortran instrinsic
! get_command_argument is that for i_arg = 0, this routine returns the 
! command line with the name of the executable removed from the beginning of the line. 
! get_command_argument, on the other hand returns the name of the
! executable when the argument is 0.
!
! Input:
!   i_arg -- Integer: Index of argument to return.
!              i_arg = 0 => Entire line minus the executable string.
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
  integer i

!

if (i_arg == 0) then
  call get_command_argument (1, arg)
  do i = 2, cesr_iargc()
    call get_command_argument (i, raw)
    arg = trim(arg) // ' ' // trim(raw)
  enddo
else
  call get_command_argument(i_arg, arg)
endif

end subroutine cesr_getarg

end module command_line_mod
