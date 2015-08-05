                                    
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
! System Libraries that need to be linked to:
!   readline curses
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

logical :: wait, flush
logical :: init_needed = .true.

! Init

if (init_needed) then
  call init_tty_char
  init_needed = .false.
endif

! Get character

call get_tty_char_c(ic, wait, flush)
this_char = achar(ic)


end subroutine get_tty_char

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
! System Libraries that need to be linked to:
!   readline curses
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

do 
  call get_tty_char (this_char, wait, .false.)  ! no flush
  if (.not. wait) return
  if (present(ignore_this)) then
    if (this_char /= achar(0) .and. all(this_char /= ignore_this)) return
  else
    if (this_char /= achar(0)) return  ! finished if not a null char
  endif
enddo

end subroutine get_a_char

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine read_a_line (prompt, line_out, trim_prompt, prompt_color, prompt_bold)
!
! Subroutine to read a line of input from the terminal.
! The line is also add to the history buffer so that the up-arrow
! and down-arrow keys can be used to recall past commands.
!
! Modules needed:
!   use input_mod
!
! System Libraries that need to be linked to:
!   readline curses
!
! Input:
!   prompt        -- character(*): Prompt string to use.
!   trim_prompt   -- logical, optional: If present and True then trim the 
!                     prompt string and add a single blank before printing 
!                     the prompt string. Default is True.
!   prompt_color  -- character(*), optional: Color of the prompt. Possibilities are: 
!                     'BLACK', 'RED', 'GREEN', 'YELLOW', 'BLUE', 'MAGENTA', 'CYAN', 'GRAY', 'DEFAULT'.
!                     The 'DEFAULT' setting (the default) does not set the prompt color.
!   prompt_bold   -- logical, optional: If present and True then the prompt will be printed in bold.
!
! Output:
!   line_out   -- Character(*): Line typed by the user.
!-

subroutine read_a_line (prompt, line_out, trim_prompt, prompt_color, prompt_bold)

use sim_utils

implicit none

character(*) prompt, line_out
character(*), optional :: prompt_color
character(16) pre, post
character(*), parameter :: r_name = 'read_a_line'

logical, optional :: trim_prompt, prompt_bold

!

pre = ''
post = ''

if (present(prompt_color)) then
  select case (prompt_color)
  case ('BLACK')
    pre = black_color
    post = reset_color
  case ('RED')
    pre = red_color
    post = reset_color
  case ('GREEN')
    pre = green_color
    post = reset_color
  case ('YELLOW')
    pre = yellow_color
    post = reset_color
  case ('BLUE')
    pre = blue_color
    post = reset_color
  case ('MAGENTA')
    pre = magenta_color
    post = reset_color
  case ('CYAN')
    pre = cyan_color
    post = reset_color
  case ('GRAY')
    pre = gray_color
    post = reset_color
  case ('', 'DEFAULT')
  case default
    call out_io (s_warn$, r_name, 'Bad color name: ' // prompt_color)
  end select
endif

if (logic_option(.false., prompt_bold)) then
  pre = trim(pre) // bold_color
  post = reset_color
endif

if (logic_option(.true., trim_prompt)) then
  call read_line (trim(pre) // trim(prompt) // ' ' // trim(post) // achar(0), line_out)  
else
  call read_line (trim(pre) // prompt // trim(post) // achar(0), line_out)
endif

end subroutine read_a_line

end module
