!+
! Subroutine tao_set_global_cmd (who, set_value)
!
! Routine to set global variables
! 
! Input:
!   who       -- Character(*): which global variable to set
!   set_value -- Character(*): Value to set to.
!
! Output:
!    %global  -- Global variables structure.
!-

subroutine tao_set_global_cmd (who, set_value)

use tao_mod
use quick_plot

implicit none

type (tao_global_struct) global

character(*) who, set_value
character(20) :: r_name = 'tao_set_global_cmd'

integer iu, ios

namelist / params / global

! open a scratch file for a namelist read

iu = lunget()
open (iu, status = 'scratch')
write (iu, *) '&params'
write (iu, *) ' global%' // trim(who) // ' = ' // trim(set_value)
write (iu, *) '/'
rewind (iu)
global = s%global  ! set defaults
read (iu, nml = params, iostat = ios)
close (iu)

if (ios == 0) then
  s%global = global
else
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
endif

end subroutine




