module tao_cmd_history_mod

use output_mod
use tao_struct

type cmd_history_struct  ! record the command history
  character(100) cmd     ! the command
  integer :: ix = 0      ! command index (1st command has ix = 1, etc.)
  logical cmd_file       ! Did command come from a command file
end type

type (cmd_history_struct), private, save :: history(100) ! command history
integer, private, save :: ix_history = 0 ! present index to command history array
integer, private, save :: n_history      ! present history index

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_cmd_history_record (cmd)
!
! Subroutine to record a cmd in the command history stack
!-

subroutine tao_cmd_history_record (cmd)

implicit none

character(*) cmd

!

ix_history = ix_history + 1
if (ix_history > size(history)) ix_history = 1
n_history = n_history + 1
history(ix_history)%cmd = cmd
history(ix_history)%ix = n_history

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_cmd_history_print (n_print)
!
! Subroutine to print the command history.
!-

subroutine tao_cmd_history_print (n_print)

implicit none

integer n_print, i
character(120) line
character(16) :: r_name = 'tao_history'

!

if (ix_history == 0) return

if (n_print > 0) then
  i = mod (ix_history - n_print + 1, size(history))
else
  i = mod (ix_history + 1, size(history))
endif

if (i < 1) i = i + size(history)

do
  if (history(i)%ix /= 0) then
    write (line, '(i4, 2a)') history(i)%ix, ': ', trim(history(i)%cmd)
    call out_io (s_blank$, r_name, line)
  endif
  if (i == ix_history) return
  i = i + 1
  if (i > size(history)) i = 1
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_history_cmd (string, err)
!
! Subroutine to print the command history.
!-

subroutine tao_history_cmd (string, err)

implicit none

integer ios, ix1, ix, ix_rec
character(*) string
character(100) cmd_out
character(16) :: r_name = 'tao_history_cmd'
logical err

!

if (string == ' ') then
  call tao_cmd_history_print (50)
  return
endif

!

err = .true.

if (index('-+0123456789', string(1:1)) /= 0) then  ! number
  read (string, *, iostat = ios) ix_rec
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'ERROR READING HISTORY NUMBER')
    return
  endif

  if (ix_rec > 0) then
    if (ix_rec > n_history .or. ix_rec < n_history - (size(history) - 1)) then
      call out_io (s_error$, r_name, 'INVALID INDEX FOR THE HISTORY LIST.')
      return
    endif
    ix = ix_rec + ix_history - n_history
  else
    if (-ix_rec > size(history) - 1 .or. -ix_rec > n_history - 1) then 
      call out_io (s_error$, r_name, 'INVALID INDEX FOR THE HISTORY LIST.')
      return
    endif
    ix = ix_history + ix_rec
  endif

  if (ix < 1) ix = ix + size(history)

!

else

  ix = ix_history
  do

    if (index(history(ix)%cmd, trim(string)) == 1) exit

    ix = ix - 1
    if (ix < 1) ix = ix + size(history)

    if (ix == ix_history .or. history(ix)%ix == 0) then
      call out_io (s_error$, r_name, 'COMMAND NOT FOUND IN THE HISTORY LIST.')
      return
    endif

  enddo

endif

! put the command in the common area so it can be used next.

tao_com%cmd = history(ix)%cmd
tao_com%use_cmd_here = .true.

err = .false.

end subroutine

end module

