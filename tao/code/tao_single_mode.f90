!+
! Subroutine tao_single_mode (char)
!
! Subroutine in interpret a command (which is generally a single character) in 
! single character mode.
!
! Input:
!   char -- Character(1): Command. If more characters are needed to complete
!             the command then this routine will get them.
!-

subroutine tao_single_mode (char)

use input_mod
use tao_set_mod, dummy => tao_single_mode
use tao_top10_mod, only:  tao_show_constraints, tao_var_write
use tao_x_scale_mod, only:  tao_x_scale_cmd, tao_x_scale_plot, tao_x_scale_graph
use tao_scale_mod, only:  tao_scale_cmd, tao_scale_graph
use tao_get_user_input_mod, only: tao_get_user_input

implicit none

type (lat_struct), pointer :: lat
type (ele_pointer_struct), allocatable :: eles(:)
type (tao_universe_struct), pointer :: u

integer i, j, ix, ix2, ix_plot, ie, iv, factor, ix_key, ios, ir
integer n1, n2, i_ele, ix_var

real(rp) this_factor, value, this_min, this_max, this_merit
real(rp) m1, m2

real(rp) :: delta_factor(4) = (/ 10, -10, 1, -1 /)

character(10) :: set_char(4) = &
            (/'!@#$%^&*()', 'QWERTYUIOP', '1234567890', 'qwertyuiop' /)

character(1) :: char
character(2) char2
character(80) this_file
character(40) :: this_opt, r_name = 'tao_single_mode'
character(200) line

logical doit, found, err, abort, err_is_fatal



!-------------------------------------------------------
! Look for custom command key binding

do i = 1, size(s%com%key)
  if (char == s%com%key(i)%name) then
    call tao_command (s%com%key(i)%expanded_str, err, err_is_fatal)
    return
  endif
enddo


!-------------------------------------------------------
! is it a variable change?

do i = 1, 4
  if (index(set_char(i), char) == 0) cycle

  ix_key = index(set_char(i), char)
  this_factor = delta_factor(i)
  ix = ix_key + 10*s%com%ix_key_bank
  if (ix > size(s%key)) then
    call out_io (s_error$, r_name, 'KEY NOT BOUND TO VARIABLE.')
    return
  endif

  ix_var = s%key(ix)
  if (ix_var == 0) cycle
  if (ix_var == -1) then
    call out_io (s_error$, r_name, 'KEY NOT BOUND TO A VARIABLE THAT EXISTS.')
    return
  endif

  value = s%var(ix_var)%model_value + this_factor * s%var(ix_var)%key_delta
  call tao_set_var_model_value (s%var(ix_var), value)
  this_merit = tao_merit ()
  return
enddo



!--------------------------------------------------------
! not a variable change so must be something else...

select case (char)

case ('?')
    call tao_help ('keys', '')

! 'g' Go run optimizer.

! 'c' Show constraints

case ('c')

  this_merit = tao_merit ()
  write (*, *)
  write (*, *) 's%global%optimizer:       ', trim(s%global%optimizer)
  write (*, *)
  call tao_show_constraints (0, 'ALL')
  call tao_show_constraints (0, 'TOP10')

case ('g')

  do
    call tao_run_cmd ('', abort)
    if (abort) exit
  enddo

! 'v' Show variables

case ('v')
  write (*, *)
  call tao_var_write ('')

! 'V' Show optimization variables 

case ('V')
  write (*, *)
  call tao_var_write ('', .true.)

! 'Z' Quit single character input mode.

case ('Z')
  s%com%single_mode = .false.
  call out_io (s_blank$, r_name, ' ', 'Entering line mode...')

! 'CR' Do nothing

case (achar(10))
  return

! '<' Scale deltas by factor of 1/2

case ('<')
  s%var(:)%key_delta = s%var(:)%key_delta / 2

! '>' Scale deltas by factor of 2

case ('>')
  s%var(:)%key_delta = s%var(:)%key_delta * 2
        
! "'" Accept command

case ("'")
  s%com%single_mode = .false.
  do
    call tao_get_user_input (line)
    found = .false.
    if (associated(tao_hook_command_ptr)) call tao_hook_command_ptr (line, found)
    if (.not. found) call tao_command (line, err, err_is_fatal)
    if (s%com%saved_cmd_line /= '') cycle ! Use up all commands on line
    ! Keep on going if a command file is open
    ix = s%com%cmd_file_level
    if (ix == 0 .and. .not. s%com%cmd_file(ix)%paused) exit  
  enddo
  s%com%single_mode = .true.

! '<CR>' Just replot.
                               
case(achar(13))  ! Ignore a <CR>

! Error:

case default
  write (*, *) 'What is this you are typing?', iachar(char)
  s%com%single_mode_buffer = ''

!--------------------------------------------------------
! 'Escape' -> Must be an arrow key. Look for the rest of the sequence:

case (achar(27))

  call tao_get_user_input (char2(1:1), wait_flag = .false.)
  call tao_get_user_input (char2(2:2), wait_flag = .false.)

  select case (char2)

  ! '<left_arrow>' Shift key bank down by 10

  case ('[D')
    s%com%ix_key_bank = max(s%com%ix_key_bank-1, 0)

  ! '<right_arrow>': Shift key bank up by 10

  case ('[C')
    s%com%ix_key_bank = min(s%com%ix_key_bank+1, size(s%key-1)/10)

  ! '<up_arrow>': Increase deltas by factor of 10.

  case ('[A')
    s%var(:)%key_delta = s%var(:)%key_delta * 10

  ! '<down_arrow>': Decrease deltas by factor of 10.

  case ('[B')
    s%var(:)%key_delta = s%var(:)%key_delta / 10

  ! 'ESC' Error:

  case default
    write (*, *) 'What is this you are typing?', char2
    s%com%single_mode_buffer = ''

  end select

!--------------------------------------------------------
! '-' Look for rest of sequence

case ('-')

  call tao_get_user_input (char, wait_flag = .true.)

  select case (char)

  ! '-p' Toggle printing.

  case ('p')
    s%global%plot_on = .not. s%global%plot_on
    if (s%global%plot_on) then
      call tao_init_plotting (' ')
      call out_io (s_blank$, r_name, 'Plotting is now: On')
    else
      call out_io (s_blank$, r_name, 'Plotting is now: Off')
    endif

  case default
    write (*, *) 'What is this you are typing?', iachar(char)
    s%com%single_mode_buffer = ''

  end select

!--------------------------------------------------------
! '=' Look for rest of sequence

case ('=')

  call tao_get_user_input (char, wait_flag = .true.)

  select case (char)

  ! '=v' Set variable at given value

  case ('v')
    if (.not. allocated(s%key)) then
      call out_io (s_error$, r_name, 'THE KEY TABLE HAS NOT BEEN SET UP IN THE INIT FILES!')
      return
    endif
    call read_this_input ('Enter Key# and Value: ')
    read (line, *, iostat = ios) ix, value
    if (ios /= 0 .or. .not. is_integer(line)) then
      call out_io (s_error$, r_name, 'I DO NOT UNDERSTAND THIS. NOTHING CHANGED.')
      return
    endif
    if (ix == 0) ix = 10
    ix2 = ix + 10*s%com%ix_key_bank
    if (ix < 0 .or. ix > 10 .or. ix2 > size(s%key)) then
      call out_io (s_error$, r_name, 'KEY NUMBER OUT OF RANGE (0 - 9).')
      return
    endif
    call tao_set_var_model_value (s%var(s%key(ix2)), value)
  
  !----------------------------------
  ! 'Escape' -> Must be an arrow key. Look for the rest of the sequence:

  case (achar(27))

    call tao_get_user_input (char2(1:1), wait_flag = .false.)
    call tao_get_user_input (char2(2:2), wait_flag = .false.)

    select case (char2)

    ! '= <left_arrow>': Paste saved variable value to variable

    case ('[D')
      do i = 1, size(s%key)
        ix_var = s%key(i)
        call tao_set_var_model_value (s%var(ix_var), s%var(ix_var)%key_val0)
      enddo

    ! '= <right_arrow>': Copy variable value to saved

    case ('[C')
      do i = 1, size(s%key)
        ix = s%key(i)
        s%var(ix)%key_val0 = s%var(ix)%model_value
      enddo

    ! '= ESC' Error:

    case default
      write (*, *) 'What is this you are typing?', char2
      s%com%single_mode_buffer = ''

    end select

  case default
    write (*, *) 'What is this you are typing?', iachar(char)
    s%com%single_mode_buffer = ''

  end select

!--------------------------------------------------------
! 'a' Look for rest of sequence

case ('a')

  call tao_get_user_input (char, wait_flag = .true.)

  select case (char)

  !----------------------------------
  ! 'Escape' -> Must be an arrow key. Look for the rest of the sequence:

  case (achar(27))

    call tao_get_user_input (char2(1:1), wait_flag = .false.)
    call tao_get_user_input (char2(2:2), wait_flag = .false.)

    select case (char2)

    ! 'a <left_arrow>': Pan left half screen

    case ('[D')
      call x_pan (-0.5_rp)

    ! 'a <right_arrow>': Pan right half screen

    case ('[C')
      call x_pan (0.5_rp)

    ! 'a <up_arrow>': Pan up half screen

    case ('[A')
      call y_pan (0.5_rp)

    ! 'a <down_arrow>': Pan down half screen

    case ('[B')
      call y_pan (-0.5_rp)

    ! 'a ESC' Error:

    case default
      write (*, *) 'What is this you are typing?', char2
      s%com%single_mode_buffer = ''

    end select

  ! 'a' Error:

  case default
    write (*, *) 'What is this you are typing?', iachar(char)
    s%com%single_mode_buffer = ''

  end select

!--------------------------------------------------------
! 's' Look for rest of sequence

case ('s')

  call tao_get_user_input (char, wait_flag = .true.)

  select case (char)

  !----------------------------------
  ! 'Escape' -> Must be an arrow key. Look for the rest of the sequence:

  case (achar(27))

    call tao_get_user_input (char2(1:1), wait_flag = .false.)
    call tao_get_user_input (char2(2:2), wait_flag = .false.)

    select case (char2)

    ! 's <left_arrow>': X-Scale (scale) plots by a factor of 0.5

    case ('[D')
      call x_scale (0.5_rp)

    ! 's <right_arrow>': X-Scale (scale) plots by a factor of 2/-

    case ('[C')
      call x_scale (2.0_rp)

    ! 's <up_arrow>': Y-Scale (scale) plots by a factor of 2.0

    case ('[A')
      call y_scale (2.0_rp)

    ! 's <down_arrow>': Y-Scale (scale) plots by a factor of 0.5

    case ('[B')
      call y_scale (0.5_rp)

    ! 's ESC' Error:

    case default
      write (*, *) 'What is this you are typing?', char2
      s%com%single_mode_buffer = ''

    end select

  ! 's' Error:

  case default
    write (*, *) 'What is this you are typing?', iachar(char)
    s%com%single_mode_buffer = ''

  end select


!--------------------------------------------------------
! 'z' Look for rest of sequence

case ('z')

  call tao_get_user_input (char, wait_flag = .true.)

  select case (char)

  !----------------------------------
  ! 'Escape' -> Must be an arrow key. Look for the rest of the sequence:

  case (achar(27))

    call tao_get_user_input (char2(1:1), wait_flag = .false.)
    call tao_get_user_input (char2(2:2), wait_flag = .false.)

    select case (char2)

    ! 'z <left_arrow>': Y-Zoom (scale) plots by a factor of 0.5

    case ('[D')
      call x_zoom (0.5_rp)

    ! 'z <right_arrow>': Y-Zoom (scale) plots by a factor of 2.0

    case ('[C')
      call x_zoom (2.0_rp)

    ! 'z <up_arrow>': Y-Zoom (scale) plots by a factor of 2.0

    case ('[A')
      call y_zoom (2.0_rp)

    ! 'z <down_arrow>': Y-Zoom (scale) plots by a factor of 0.5

    case ('[B')
      call y_zoom (0.5_rp)

    ! 'z ESC' Error:

    case default
      write (*, *) 'What is this you are typing?', char2
      s%com%single_mode_buffer = ''

    end select

  ! 'z' Error:

  case default
    write (*, *) 'What is this you are typing?', iachar(char)
    s%com%single_mode_buffer = ''

  end select


!--------------------------------------------------------
! '/' Look for rest of sequence

case ('/')

  call tao_get_user_input (char, wait_flag = .true.)

  select case (char)

  ! '/e' Element

  case ('e')
    call read_this_input (' Element Name or Index: ')
    if (ios /= 0 .or. line == '') return
    call tao_locate_elements (line, 0, eles, err)
    u => tao_pointer_to_universe(0)  
    lat => u%model%lat
    do i = 1, size(eles)
      write (*, *) '!---------------------------------------------------'
      write (*, *) '! Branch Index: ', eles(i)%ele%ix_branch
      write (*, *) '! Element Index:', eles(i)%ele%ix_ele
      call type_ele (eles(i)%ele, .false., 6, .true., radians$, .true.)
    enddo

    if (size(eles) == 0) write (*, *) 'Element not found.'

  ! '/l' Lattice list

  case ('l')
    do i = lbound(s%u, 1), ubound(s%u, 1)
      lat => s%u(i)%model%lat
      write (*, *)
      write (*, *) 'Lat: ', lat%lattice, i
      write (*, *) 'Ix  Name                   S  Beta_a  Beta_b'
      do j = 1, lat%n_ele_track
        write (*, '(i3, 2x, a, f8.2, 2f8.2)') j, lat%ele(j)%name, &
                 lat%ele(j)%s, lat%ele(j)%a%beta, lat%ele(j)%b%beta
      enddo
    enddo

  ! '/b' Switch default branch

  case ('b')

    call read_this_input (' Default Lattice Branch: ')
    call string_trim (line, line, ix)
    call tao_set_default_cmd ('branch', line)

  ! '/p' Just create a postscript file.

  case ('p')
    call tao_write_cmd ('ps')

  ! '/P' Print a hardcopy.

  case ('P')
    call tao_write_cmd ('hard')

  ! '/u' Default universe 

  case ('u')

    call read_this_input (' Default Universe: ')
    call string_trim (line, line, ix)
    call tao_set_default_cmd ('universe', line)

  ! '/v' Output to default file.

  case ('v')
    call tao_var_write (s%global%var_out_file)

  ! '/x' Scale horizontal axis

  case ('x')

    call read_this_input (' Input x-axis min, max: ')
    call string_trim (line, line, ix)
    if (ix == 0) then
      call tao_x_scale_cmd ('*', 0.0_rp, 0.0_rp, err)
    else
      read (line, *, iostat = ios) m1, m2
      if (ios /= 0) then
        write (*, *) 'ERROR READING MIN/MAX.'
        return
      endif
      call tao_x_scale_cmd ('*', m1, m2, err)
    endif

  ! '/y' Scale y-axis

  case ('y')
    call read_this_input (' Input scale min, max: ')
    call string_trim (line, line, ix)
    if (ix == 0) then
      call tao_scale_cmd ('*', 0.0_rp, 0.0_rp)
    else
      read (line, *, iostat = ios) m1, m2
      if (ios /= 0) then
        write (*, *) 'ERROR READING MIN/MAX.'
        return
      endif
      call tao_scale_cmd ('*', m1, m2)
    endif

  ! '/' Error:

  case default
    write (*, *) 'What is this you are typing?', iachar(char)
    s%com%single_mode_buffer = ''

  end select

end select

!-----------------------------------------------------------------------
contains

subroutine x_scale (factor)

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

real(rp) factor
integer i, j

!

do i = 1, size(s%plot_page%region)
  if (.not. s%plot_page%region(i)%visible) cycle
  plot => s%plot_page%region(i)%plot
  if (.not. allocated(plot%graph)) cycle
  if (plot%autoscale_gang_x) then
    call tao_x_scale_plot (plot, plot%graph(1)%x%min * factor, plot%graph(1)%x%max * factor)
  else
    do j = 1, size(plot%graph)
      graph => plot%graph(i)
      call tao_x_scale_graph (graph, graph%x%min * factor, graph%x%max * factor)
    enddo
  endif
enddo

end subroutine

!-----------------------------------------------------------------------
! contains

subroutine y_scale (factor)

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

integer i, j
real(rp) factor

!

do i = 1, size(s%plot_page%region)
  plot => s%plot_page%region(i)%plot
  if (.not. allocated (plot%graph)) cycle
  do j = 1, size(plot%graph)
    graph => plot%graph(j) 
    call tao_scale_graph (graph, graph%y%min * factor, graph%y%max * factor, 'y')
    if (graph%y2%min == graph%y2%max) cycle
    call tao_scale_graph (graph, graph%y2%min * factor, graph%y2%max * factor, 'y2')
  enddo
enddo

end subroutine

!-----------------------------------------------------------------------
! contains

subroutine x_zoom (factor)

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

integer i, j
real(rp) factor, w, c

!

do i = 1, size(s%plot_page%region)
  plot => s%plot_page%region(i)%plot
  if (.not. allocated(plot%graph)) cycle
  do j = 1, size(plot%graph)
    graph => plot%graph(j)
    w = (graph%x%max - graph%x%min) * factor / 2
    c = (graph%x%max + graph%x%min) / 2
    call tao_x_scale_graph (graph, c-w, c+w)
  enddo
enddo

end subroutine

!-----------------------------------------------------------------------
! contains

subroutine y_zoom (factor)

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

integer i, j
real(rp) factor, w, c

!

do i = 1, size(s%plot_page%region)
  plot => s%plot_page%region(i)%plot
  if (.not. allocated (plot%graph)) cycle
  do j = 1, size(plot%graph)
    graph => plot%graph(j) 

    w = (graph%y%max - graph%y%min) * factor / 2
    c = (graph%y%max + graph%y%min) / 2
    call tao_scale_graph (graph, c-w, c+w, 'y')

    w = (graph%y2%max - graph%y2%min) * factor / 2
    c = (graph%y2%max + graph%y2%min) / 2
    if (w == 0) cycle
    call tao_scale_graph (graph, c-w, c+w, 'y2')

  enddo
enddo

end subroutine

!-----------------------------------------------------------------------
! contains

subroutine x_pan (factor)

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

integer i, j
real(rp) factor, w

!

do i = 1, size(s%plot_page%region)
  plot => s%plot_page%region(i)%plot
  if (.not. allocated(plot%graph)) cycle
  do j = 1, size(plot%graph)
    graph => plot%graph(j) 
    w = factor * (graph%x%max - graph%x%min)
    call tao_x_scale_graph (graph, graph%x%min + w, graph%x%max + w)
  enddo
enddo

end subroutine

!-----------------------------------------------------------------------
! contains

subroutine y_pan (factor)

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

integer i, j
real(rp) factor, w

!

do i = 1, size(s%plot_page%region)
  plot => s%plot_page%region(i)%plot
  if (.not. allocated (plot%graph)) cycle
  do j = 1, size(plot%graph)
    graph => plot%graph(j)

    w = factor * (graph%y%max - graph%y%min)
    call tao_scale_graph (graph, graph%y%min + w, graph%y%max + w, 'y')

    w = factor * (graph%y2%max - graph%y2%min)
    if (w == 0) cycle
    call tao_scale_graph (graph, graph%y2%min + w, graph%y2%max + w, 'y2')

  enddo
enddo

end subroutine

!-----------------------------------------------------------------------
! contains

subroutine read_this_input (string)

character(*) string

if (s%global%wait_for_CR_in_single_mode) then
  line = ''
  do i = 1, len(line)
    read '(a)', line(i:i)
    if (line(i:i) == ' ') exit
  enddo

else
  write (*, '(a)', advance = "NO") string
  read (*, '(a)', iostat = ios) line
endif

end subroutine

end subroutine   
