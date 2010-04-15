module tao_command_mod

use output_mod
use tao_mod

type cmd_history_struct  ! record the command history
  character(100) cmd     ! the command
  integer :: ix = 0      ! command index (1st command has ix = 1, etc.)
  logical cmd_file       ! Did command come from a command file
end type

type (cmd_history_struct), private, save :: history(1000) ! command history
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
history(ix_history)%ix = n_history
if (tao_com%cmd_from_cmd_file) then
  history(ix_history)%cmd = '  ! ' // trim(cmd)
else
  history(ix_history)%cmd = trim(cmd)
endif

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

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_cmd_split (cmd_line, n_word, cmd_word, no_extra_words, err, separator)
!
! This routine splits the command line into words.
!
! Input: 
!   cmd_line       -- Character(*): The command line.
!   n_word         -- Integer: Maximum number of words to split command line into.
!   no_extra_words -- Logical: are extra words allowed at the end?
!                        If False then cmd_word(n_word) will contain everything after 
!                        the n_word-1 word.
!   separator      -- Character(*), optional: a list of characters that,
!                        besides a blank space, signify a word boundary. 
!
! Output:
!   cmd_word(n_word) -- Character(*): The individual words.
!   err              -- Logical: error in splitting words
!
! For example: 
!   separator = '-+' 
!   cmd_line = 'model-design'
! Restult:
!   cmd_word(1) = 'model'
!   cmd_word(2) = '-'
!   cmd_word(3) = 'design'
!
! Anything between single or double quotes is treated as a single word.
! Whitespace or a separator inside of "{}", "()", or "[]" is ignored.
!-

subroutine tao_cmd_split (cmd_line, n_word, cmd_word, no_extra_words, err, separator)

integer i, ix, nw, n_word, ix_b1, ix_b2, ix_b3

character(*) cmd_line
character(*), optional :: separator
character(*) cmd_word(:)
character(16) :: r_name = 'tao_cmd_split'
character(200) line
character(1), parameter :: tab = char(9)

logical err
logical no_extra_words

!

err = .false.
line = cmd_line
cmd_word(:) = ''
nw = 0

nw_loop: do 
  call string_trim (line, line, ix)

  if (ix == 0) exit

  ! If extra words allowed, everything left goes into cmd_word(n_word)
  if (nw == n_word - 1 .and. .not. no_extra_words) then
    nw=nw+1; cmd_word(nw) = trim(line)
    return
  endif

  if (line(1:1) == '"') then
    ix = index(line(2:), '"')
    if (ix == 0) ix = len(line)
    nw=nw+1; cmd_word(nw) = line(2:ix)
    line = line(ix+1:)
    cycle
  elseif (line(1:1) == "'") then
    ix = index(line(2:), "'")
    if (ix == 0) ix = len(line)
    nw=nw+1; cmd_word(nw) = line(2:ix)
    line = line(ix+1:)
    cycle
  endif

  ix_b1 = 0; ix_b2 = 0; ix_b3 = 0
  do i = 1, len(line)

    if (line(i:i) == '{') ix_b1 = ix_b1 + 1
    if (line(i:i) == '}') ix_b1 = ix_b1 - 1
    if (line(i:i) == '(') ix_b2 = ix_b2 + 1
    if (line(i:i) == ')') ix_b2 = ix_b2 - 1
    if (line(i:i) == '[') ix_b3 = ix_b3 + 1
    if (line(i:i) == ']') ix_b3 = ix_b3 - 1

    if (ix_b1 /= 0 .or. ix_b2 /= 0 .or. ix_b3 /= 0) cycle

    if (present(separator)) then
      if (index(separator, line(i:i)) /= 0) then
        if (i /= 1) then
          nw=nw+1; cmd_word(nw) = line(1:i-1)
          line = line(i:)
        endif
        nw=nw+1; cmd_word(nw) = line(1:1)
        line = line(2:)
        cycle nw_loop
      endif
    endif

    if (line(i:i) == ' ' .or. line(i:i) == tab) then
      nw=nw+1; cmd_word(nw) = line(1:i-1)
      line = line(i+1:)
      cycle nw_loop
    endif

  enddo

  if (ix_b1 /= 0 .or. ix_b2 /= 0 .or. ix_b3 /= 0) then
    call out_io (s_error$, r_name, 'MISMATCHED "{...}", "(...)", OR "[...]".')
    err = .true.
    return
  endif

  call out_io (s_fatal$, r_name, 'INTERNAL ERROR!')
  call err_exit

enddo nw_loop

!  

if (no_extra_words .and. nw > n_word) then
  call out_io (s_error$, r_name, 'EXTRA STUFF ON COMMAND LINE: ' // line)
  err = .true.
endif

end subroutine tao_cmd_split

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine tao_next_switch (line, switch_list, switch, err, ix_word)
!
! Subroutine return the next switch on the command line.
! A switch is assumed to begin with the '-' character so that
! if the first word on the command line starts with '-' it must
! match something on the switch_list list.
! Switch abbreviations are permitted.
!
! Input:
!   line            -- Character(*): Command line
!   switch_list(*)  -- Character(*): List of valid switches. 
!
! Output:
!   line            -- Character(*): Line with switch word removed. If the first word
!                       does not look like a switch then nothing is removed.
!   switch          -- Character(*): Switch present. This is the full name
!                       even if what was on the command line was an abbreviation.
!                       Return '' if not a switch.
!   err             -- Logical: Set True if switch is not recognized.
!                       An error message will be printed.
!   ix_word         -- Integer: Character length of first word left on line.
!-

subroutine tao_next_switch (line, switch_list, switch, err, ix_word)

character(*) line, switch, switch_list(:)
character(20) :: r_name = 'tao_next_switch'
logical err

integer ix, n, ix_word

!

err = .false.
switch = ''

call string_trim(line, line, ix_word)
if (ix_word == 0) return
if (line(1:1) /= '-') return

call match_word (line(:ix_word), switch_list, n, .true., switch)
if (n < 1) then
  err = .true.
  if (n == 0) then 
    call out_io (s_error$, r_name, 'UNKNOWN SWITCH: ' // line(:ix_word))
  else
    call out_io (s_error$, r_name, 'AMBIGUOUS SWITCH: ' // line(:ix_word))
  endif
  return
endif

call string_trim(line(ix_word+1:), line, ix_word)

end subroutine

end module

