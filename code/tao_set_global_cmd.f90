!+
! Subroutine tao_set_global_cmd (s, who, set_value)
!
! Routine to set global variables
! 
! Input:
!   set_value -- Character(*): Value to set to.
!
!  Output:
!   s         -- tao_super_universe_struct:
!    %global    -- Global variables structure.
!-

subroutine tao_set_global_cmd (s, who, set_value)

use tao_mod
use quick_plot

implicit none

type (tao_super_universe_struct) s

character(*) who, set_value
character(20) :: r_name = 'tao_set_global_cmd'

logical err

!

select case (who)
case ('opt_with_ref') 
  call set_logical (s%global%opt_with_ref, set_value)
case ('opt_with_base') 
  call set_logical (s%global%opt_with_base, set_value)
case default
end select

!----------------------------------------------------------------------------
contains

subroutine set_logical (logic, set_str)

logical logic, err
character(*) set_str

!

logic = eval_logical (set_str, err)
if (err) call out_io (s_error$, r_name, 'BAD LOGICAL VALUE')

end subroutine 

end subroutine




