!+
! Subroutine tao_timer (command)
!
! Routine to stop/start the timer for timing computation speed.
!
!
! Input:
!   command -- Character(*): 'start', 'read', or 'beam'
!-

subroutine tao_timer (command)

use tao_struct

implicit none

real(rp) dtime

character(*) command
character(20), parameter :: r_name = 'tao_timer'

!

select case (command)
case ('start')
  if (s%global%beam_timer_on) call out_io (s_info$, r_name, 'beam timer mode is now off')
  s%global%beam_timer_on = .false.
  call run_timer ('START')

case ('read')
  call run_timer ('READ', dtime)
  if (dtime > 1000) then
    call out_io (s_info$, r_name, 'Time (min): \f11.2\ ', dtime/60)
  else
    call out_io (s_info$, r_name, 'Time (sec): \f11.2\ ', dtime)
  endif

case ('beam')
  s%global%beam_timer_on = .not. s%global%beam_timer_on
  call out_io (s_info$, r_name, 'beam timer mode is now: ' // on_off_logic(s%global%beam_timer_on))

case default
  call out_io (s_error$, r_name, 'UNKNOWN COMMAND: ' // command)

end select

end subroutine
