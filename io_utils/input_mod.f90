#include "CESR_platform.inc"
                                    
!+
! Modue input_mod
!
! Module for doing single character input from the terminal.
! That is, each keystroke is made available to the program as it occurs.
!-

module input_mod

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_tty_char (this_char, wait, flush)
!
! Subroutine for getting a single character from the terminal.
! Also see: get_a_char
!
! Modules needed:
!   use input_mod
!
! Input:
!   wait      -- Logical: If True then routine will wait until a keystroke
!                  has occured. If False and no keystroke is in the buffer then
!                  achar(0) will be returned as this_char.
!   flush     -- Logical: If True then the keystroke buffer will be cleared 
!                  first before any processing.
!
! Output:
!   this_char -- Character(1): Character returned
!-

subroutine get_tty_char (this_char, wait, flush)

  implicit none

  integer ic

  character this_char
!  character ttychr     ! for VMS

  logical :: wait, flush
  logical :: init_needed = .true.

!

  if (init_needed) then
    call init_tty_char
    init_needed = .false.
  endif

! Unix version

  call get_tty_char_c(ic, wait, flush)
  this_char = achar(ic)

! VMS version

!  call milli_sleep(10)
!
!  do
!    char = ttychr()
!    if (.not. wait) return
!    call milli_sleep(100)
!  enddo

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_a_char (this_char, wait, ignore_this)
!
! Subroutine for getting a single character from the terminal.
! Also see: get_tty_char
!
! Modules needed:
!   use input_mod
!
! Input:
!   wait        -- Logical: If True then routine will wait until a keystroke
!                    has occured. If False and no keystroke is in the buffer 
!                    then achar(0) will be returned as this_char.
!   ignore_this -- Character(*): List of characters to ignore. If a keystroke
!                    matches a character on this list the keystroke is ignored.
!
! Output:
!   this_char -- Character(1): Character returned
!-

subroutine get_a_char (this_char, wait, ignore_this)

  implicit none

  logical :: wait

  character this_char
  character, optional :: ignore_this(:)
  integer ios
!
#ifndef CESR_WINCVF
  do 
    call get_tty_char (this_char, wait, .false.)  ! no flush
    if (.not. wait) return
    if (present(ignore_this)) then
      if (this_char /= achar(0) .and. all(this_char /= ignore_this)) return
    else
      if (this_char /= achar(0)) return  ! finished if not a null char
    endif
  enddo
#else
  do
   read (*, '(a)', iostat = ios) this_char
   if (ios /= 0) cycle
   if (.not. wait) return
   if (present(ignore_this)) then
     if (all(this_char /= ignore_this)) return
   else
     return  ! finished if not a null char
   endif
  enddo
#endif
end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine read_a_line (prompt, line_out, trim_prompt)
!
! Subroutine to read a line of input from the terminal.
! The line is also add to the history buffer so that the up-arrow
! and down-arrow keys can be used to recall past commands.
!
! Modules needed:
!   use input_mod
!
! Input:
!   prompt      -- Character(*): Prompt string to use.
!   trim_prompt -- Logical, optional: If present and True then trim the 
!                   prompt string and add a single blank before printing 
!                   the prompt string. Default is True.
!
! Output:
!   line_out   -- Character(*): Line typed by the user.
!-

subroutine read_a_line (prompt, line_out, trim_prompt)

use cesr_utils

implicit none

character(*) prompt, line_out
logical, optional :: trim_prompt

!

#ifndef CESR_WINCVF

if (logic_option(.true., trim_prompt)) then
  call read_line (trim(prompt) // ' ' // achar(0), line_out)  
else
  call read_line (prompt // achar(0), line_out)
endif

#else

write (*, '(a)', advance = 'NO') trim(prompt)
read (*, '(a)', iostat = ios) line_out

#endif


end subroutine

end module
