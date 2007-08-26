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
use tao_var_mod
use tao_lm_optimizer_mod

implicit none

real(rp), allocatable, save :: var_vec(:)
integer n_data, i

character(*)  which
character(40) :: r_name = 'tao_run_cmd', my_opti

logical abort

!

call tao_set_var_useit_opt()
call tao_set_data_useit_opt()

if (.not. any (which /= (/ '      ', 'de    ', 'lm    ', 'lmdif ', 'custom' /))) then
  call out_io (s_error$, r_name, 'OPTIMIZER NOT RECOGNIZED: ' // which)
  return
endif

if (which /= ' ') s%global%optimizer = which
call out_io (s_blank$, r_name, 'Optimizing with: ' // s%global%optimizer)
call out_io (s_blank$, r_name, &
              "Type ``.'' to stop the optimizer before it's finished.")

call tao_get_vars (var_vec)
if (size(var_vec) == 0) then
  call out_io (s_fatal$, r_name, 'No variables to vary!')
  tao_com%optimizer_running = .false.
  return
endif

! See if there are any constraints

n_data = 0
do i = 1, size(s%u)
  n_data = n_data + count(s%u(i)%data(:)%useit_opt)
enddo
if (n_data == 0) then
  call out_io (s_error$, r_name, 'No data constraints defined for the merit function!')
  tao_com%optimizer_running = .false.
  return
endif

! Optimize...

do i = 1, s%global%n_opti_loops

  select case (s%global%optimizer)

  case ('de') 
    call tao_de_optimizer (abort)

  case ('lm') 
    call tao_lm_optimizer (abort)

  case ('lmdif')
    call tao_lmdif_optimizer (abort)

  case ('custom')
    call tao_hook_optimizer (abort)

  end select

  if (abort) exit

enddo

end subroutine

