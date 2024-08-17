!+
! Subroutine tao_quiet_set (set)
!
! Routine to set on/off silent command file running.
!
! Input:
!   set   -- logical: True is silent running is wanted.
!-

subroutine tao_quiet_set (set)

use tao_struct
implicit none

type (out_io_output_direct_struct), save :: out_dir_save
integer ix
character(*) set
character(12) name
character(*), parameter :: r_name = 'tao_quiet_set'

!

if (set == '') then
  name = 'all'

elseif (set == 'cmd-file-end') then
  if (s%global%quiet == 'All') then
    name = 'off'
  else
    name = s%global%quiet
  endif

elseif (set == 'ALL') then
  name = set

else
  call match_word (set, [character(12):: 'off', 'all', 'warnings'], ix, .false., .true., name)
  if (ix < 1) then
    call out_io (s_error$, r_name, 'BAD "quiet" COMMAND ARGUMENT: ' // set)
    return
  endif
endif

!

if (s%global%quiet == 'off') call output_direct(get = out_dir_save)
s%global%quiet = name

select case (name)
case ('off')
  call output_direct (set = out_dir_save)
case ('warnings')
  call output_direct (-1, .false., s_warn$, s_dwarn$) ! Do not print 
case ('all', 'ALL')
  call output_direct (-1, .false., s_blank$, s_dwarn$) ! Do not print 
end select

end subroutine
