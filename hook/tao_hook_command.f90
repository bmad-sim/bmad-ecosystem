!+
! Subroutine tao_hook_command (cmd_line, found)
!
! Put custom Tao commands here.
!
! This file is already set up so that it is rather simple to add a command. Just
! follow the directions. Keep in mind that you don't have to use the included 
! infrastructure if you feel like doing something really custom -- Just clear
! everything out of here and replace it with whatever you like. All that Tao 
! really cares about is that found is set to either TRUE or FALSE.
!
! Input:
!   cmd_line   -- Character(*): command line
!   found      -- Logical: Set True if the command is handled by this routine.
!                          If false, then the standard commands will be searched
!                          and Tao will spit out an error if the command is not
!                          found
!
!  Output:
!   s          -- tao_super_universe_struct: change whatever you feel like!
!-

subroutine tao_hook_command (command_line, found)

  use tao_mod

  implicit none

  logical found

  integer i, ix, ix_line, ix_cmd, which
  integer int1, int2

  real(rp) value1, value2, this_merit

  character(*) :: command_line
  character(140) :: cmd_line
  character(20) :: r_name = 'tao_hook_command'
  character(40) :: cmd_word(12)
 
  character(16) cmd_name

  !*********
  ! put your list of hook commands in here
  character(16) :: cmd_names(1) = (/ 'hook            '/)
  !*********
  
  logical quit_tao, err
  
!
! found will be set to TRUE if the command is found in here
  found = .false.
  
! blank line => nothing to do

  call string_trim (command_line, cmd_line, ix_line)
  if (ix_line == 0 .or. cmd_line(1:1) == '!') return

! strip the command line of comments

  ix = index(cmd_line, '!')
  if (ix /= 0) cmd_line = cmd_line(:ix-1)        ! strip off comments

! match first word to a command name
! If not found then found = .false.

  call match_word (cmd_line, cmd_names, ix_cmd)
  if (ix_cmd == 0) then
    found = .false.
    return
  elseif (ix_cmd < 0) then
    call out_io (s_error$, r_name, 'AMBIGUOUS HOOK COMMAND')
    found = .true.
    return
  else
    found = .true.
  endif
  cmd_name = cmd_names(ix_cmd)

! Strip off command name from cmd_line and select the appropriate command.

  call string_trim (cmd_line(ix_line+1:), cmd_line, ix_line)
  select case (cmd_name)

!----------------------------------------------------------
! PUT YOUR CUSTOM COMMANDS IN THIS CASE CONSTRUCT

  
!--------------------------------
! HOOK
  
    case ('hook')
      ! split the command line into its separate words
      ! separate words placed in cmd_word(:)
      call cmd_split(10, .true., err); if (err) return

      ! send any output to out_io
      call out_io (s_blank$, r_name, &
       "This is just a dummy command for illustration purposes")
      call out_io (s_blank$, r_name, "I will just echo anything you tell me!")
      call out_io (s_blank$, r_name, "***")
      do i = 1, size(cmd_word)
        if (cmd_word(i) .ne. ' ') &  
          call out_io (s_blank$, r_name, cmd_word(i))
      enddo
      call out_io (s_blank$, r_name, "***")
      

!----------------------------------------------------------
! no case for this command
  
    case default
      found = .false.

  end select

! Do the standard calculations and plotting after command execution.
! See the Tao manual section enititled "Lattice Calculation" in the 
! "Overview" chapter for details on what is performed here.
  call tao_cmd_end_calc

  return

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains

! This routine splits the command into words.
!
! Input: 
!  n_word         -- integer: number of words to split command line into
!  no_extra_words -- logical: are extra words allowed at the end?
!  err            -- logical: error in splitting words
!  separator      -- character(*): a list of characters that, besides a blank space,
!                                  signify a word boundary. 
!
! Output:
!  cmd_word(n_word) -- character(40): the individual words
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
!

subroutine cmd_split (n_word, no_extra_words, err, separator)

  integer i, n, n_word, ix_end_quote
  character(*), optional :: separator
  logical err
  logical no_extra_words

!

  cmd_word(:) = ' '

  if (present(separator)) then
    i = 0
    do 
      if (i == len(cmd_line)) exit
      i = i + 1
      if (index(separator, cmd_line(i:i)) /= 0) then
        cmd_line = cmd_line(:i-1) // ' ' // cmd_line(i:i) // ' ' // &
                                                            cmd_line(i+1:)
        i = i + 3
      endif
    enddo
  endif

  ix_line = 0
  ix_end_quote = 0
  
  do n = 1, n_word
    call string_trim (cmd_line(ix_line+1:), cmd_line, ix_line)
    if (cmd_line(1:1) .eq. '"') then
      ix_end_quote = index(cmd_line(2:), '"')
      if (ix_end_quote .eq. 0) ix_end_quote = len(cmd_line)
      cmd_word(n) = cmd_line(2:ix_end_quote)
    elseif (cmd_line(1:1) .eq. "'") then
      ix_end_quote = index(cmd_line(2:), "'")
      if (ix_end_quote .eq. 0) ix_end_quote = len(cmd_line)
      cmd_word(n) = cmd_line(2:ix_end_quote)
    else
      cmd_word(n) = cmd_line(:ix_line)
    endif 
    if (ix_line == 0) return
  enddo

  cmd_word(n_word) = cmd_line

  if (no_extra_words) then
    call string_trim (cmd_line(ix_line+1:), cmd_line, ix_line)
    if (ix_line /= 0) then
      call out_io (s_error$, r_name, 'EXTRA STUFF ON COMMAND LINE: ' // &
                                                                  cmd_line)
      err = .true.
    endif
  endif

end subroutine cmd_split

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! contains
! This routine converts a string to a real number.

subroutine to_real (str, r_real, err)

character(*) str
real(rp) r_real
integer ios
logical err

!

err = .false.
read (str, *, iostat = ios) r_real

if (ios /= 0) then
  call out_io (s_error$, r_name, 'EXPECTING REAL NUMBER: ' // str)
  err = .true.
  return
endif

end subroutine to_real

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! contains
! This routine converts a string to an integer

subroutine to_int (str, i_int, err)

character(*) str
integer ios, i_int
logical err

!

err = .false.
read (str, *, iostat = ios) i_int

if (ios /= 0) then
  call out_io (s_error$, r_name, 'EXPECTING INTEGER: ' // str)
  err = .true.
  return
endif

end subroutine to_int

end subroutine tao_hook_command
