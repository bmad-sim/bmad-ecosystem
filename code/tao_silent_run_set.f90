!+
! Subroutine tao_silent_run_set (set)
!
! Routine to set on/off silent command file running.
!
! Input:
!   set   -- logical: True is silent running is wanted.
!-

subroutine tao_silent_run_set (set)

use tao_struct
implicit none
logical set

!

if (set) then
  call output_direct (-1, .false., s_blank$, s_success$) ! Do not print 
else
  call output_direct (-1, print_and_capture=s%com%print_to_terminal)
endif

s%global%silent_run = set

end subroutine
