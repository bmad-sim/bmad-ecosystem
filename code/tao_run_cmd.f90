!+
! Subroutine tao_run_cmd (which)
!
! Subrutine to minimize the merit function by varying variables until
! the "data" as calculated from the model matches a constraint or measured data.
!
! Input:
!   which -- Character(*): which optimizer to use. 
!             ' '        -- Same as last time
!             'de'       -- Differential Evolution.
!             'lm'       -- Levenberg - Marquardt (aka lmdif).
!              'custom'  -- Custom routine.
!
! Output:
!-

subroutine tao_run_cmd (which)

use tao_mod
implicit none


character(*)  which
character(40) :: r_name = 'tao_run_cmd', my_opti

!

if (.not. any (which /= (/ '      ', 'de    ', 'lm    ', 'custom' /))) then
  call out_io (s_error$, r_name, 'OPTIMIZER NOT RECOGNIZED: ' // which)
  return
endif

if (which /= ' ') s%global%optimizer = which
call out_io (s_blank$, r_name, 'Optimizing with: ' // which)
call out_io (s_blank$, r_name, &
              "Type ``.'' to stop the optimizer before it's finished.")

select case (s%global%optimizer)

case ('de') 
  call tao_de_optimizer ()

case ('lm') 
  call tao_lm_optimizer ()

case ('custom')
  call tao_hook_optimizer ()

end select

end subroutine

