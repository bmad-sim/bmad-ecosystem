!+
! Subroutine tao_pause_cmd (time)
!
! Routine to pause Tao.
!
! Input:
!   time -- Real(rp): Time to pause in seconds.
!-

subroutine tao_pause_cmd (time)

use tao_struct

implicit none

real(rp) time

integer n_level

character(80) line
character(20) :: r_name = 'tao_pause_cmd'

!

if (time > 0) then
  call milli_sleep (int(1000 * time))
  return
endif

if (time == 0) then
  write (*, '(a)', advance = 'no') 'Press the CR button to continue: '
  read (*, '(a)') line
  return
endif

! 

n_level = s%com%cmd_file_level
if (n_level == 0) then
  call out_io (s_error$, r_name, 'THERE IS NO COMMAND FILE OPEN TO PAUSE')
  return
endif
s%com%cmd_file(n_level)%paused = .true.
call out_io (s_info$, r_name, 'Suspending command file reading...', &
                              'Type "continue" to continue the command file.')

end subroutine
