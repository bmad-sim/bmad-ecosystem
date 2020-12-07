!+
! Function tao_user_is_terminating_optimization () result (is_terminating)
!
! Routine to check for keyboard input of a period '.' signaling optimization termination.
!
! Output:
!   is_terminating -- logical: Set True of '.' is detected. False otherwise.
!-

function tao_user_is_terminating_optimization () result (is_terminating)

use tao_interface, dummy => tao_user_is_terminating_optimization
use input_mod

implicit none

logical is_terminating

character(52) :: r_name = 'tao_user_is_terminating_optimization'
character(1) char

!

is_terminating = .false.
if (.not. s%global%optimizer_allow_user_abort) return

do
  call get_tty_char (char, .false., .false.) 
  if (char == '.') then
    call out_io (s_blank$, r_name, 'Optimizer stop signal detected.', 'Stopping now.')
    is_terminating = .true.
    return
  endif
  if (char == achar(0)) return   ! return if there is no more input
enddo

end function tao_user_is_terminating_optimization

