!+
! Subroutine tao_quiet_set (set)
!
! Routine to set silent running.
! Note: s_important$, s_error$ and s_fatal$ levels are never supressed. 
!
! Input:
!   set   -- logical: True is silent running is wanted.
!-

subroutine tao_quiet_set (set)

use tao_struct
implicit none

integer ix, n_level
character(*) set
character(12) name
character(*), parameter :: r_name = 'tao_quiet_set'

!

n_level = s%com%cmd_file_level

if (set == '') then
  name = 'all'

elseif (set == 'cmd-file-end') then
  name = s%com%cmd_file(n_level)%quiet

else
  call match_word (set, [character(12):: 'off', 'all', 'warnings'], ix, .false., .true., name)
  if (ix < 1) then
    call out_io (s_error$, r_name, 'BAD "quiet" COMMAND ARGUMENT: ' // set)
    return
  endif
endif

!

s%com%cmd_file(n_level)%quiet = name
s%global%quiet = name

select case (name)
case ('off')
  call output_direct (-1, .true., s_blank$, s_dwarn$)
case ('warnings')
  call output_direct (-1, .true., s_blank$, s_dwarn$)
  call output_direct (-1, .false., s_warn$, s_dwarn$) ! Do not print 
case ('all')
  call output_direct (-1, .false., s_blank$, s_dwarn$) ! Do not print 
end select

end subroutine
