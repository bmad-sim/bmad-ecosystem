!+
! Subroutine tao_run_cmd (s, which)
!
! Subrutine to minimize the merit function by varying variables until
! the "data" as calculated from the model matches a constraint or measured data.
!
! Input:
!   s     -- Tao_super_universe_struct:
!   which -- Character(*): which optimizer to use. 
!             ' '        -- Same as last time
!             'de'       -- Differential Evolution.
!             'lm'       -- Levenberg - Marquardt (aka lmdif).
!              'custom'  -- Custom routine.
!
! Output:
!   s  -- Tao_super_universe_struct:
!-

subroutine tao_run_cmd (s, which)

use tao_mod
implicit none

type (tao_super_universe_struct) :: s

character(*)  which
character(40) :: r_name = 'tao_run_cmd', my_opti

!

if (.not. any (which /= (/ '      ', 'de    ', 'lm    ', 'custom' /))) then
  call out_io (s_error$, r_name, 'OPTIMIZER NOT RECOGNIZED: ' // which)
  return
endif

if (which /= ' ') s%global%optimizer = which
call out_io (s_blank$, r_name, 'Optimizing with: ' // which)

select case (s%global%optimizer)

case ('de') 
  call tao_de_optimizer (s)

case ('lm') 
  call tao_lm_optimizer (s)

case ('custom')
  call tao_hook_optimizer (s)

end select

end subroutine

