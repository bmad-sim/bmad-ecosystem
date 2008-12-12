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

use tao_struct
use tao_interface
use tao_utils
use quick_plot
use input_mod
use tao_scale_mod
use tao_x_scale_mod
use tao_plot_window_mod
use tao_top10_mod

implicit none

type (lat_struct), pointer :: lat

integer i, j, ix, ix2, ix_plot, ie, iv, factor, ix_key, ios, ir
integer n1, n2, i_ele, ix_var

real(rp) this_factor, value, this_min, this_max, this_merit
real(rp) m1, m2

real(rp) :: delta_factor(4) = (/ 10, -10, 1, -1 /)

character(10) :: set_char(4) = &
            (/'!@#$%^&*()', 'QWERTYUIOP', '1234567890', 'qwertyuiop' /)

character(1) :: char, ignore_space(1) = (/ ' ' /)
character(2) char2
character(80) this_file, line
character(40) :: this_opt, str, r_name = 'tao_single_mode'

logical doit, found, err, abort

! parse input.
! is it a variable change?

do i = 1, 4
  if (index(set_char(i), char) /= 0) then
    ix_key = index(set_char(i), char)
    this_factor = delta_factor(i)
    ix = ix_key + 10*tao_com%ix_key_bank
    if (ix > size(s%key)) then
      call out_io (s_error$, r_name, 'KEY NOT BOUND TO VARIABLE.')
      return
    endif
    ix_var = s%key(ix)
    if (ix_var == 0) cycle
    value = s%var(ix_var)%model_value + this_factor * s%var(ix_var)%key_delta
    call tao_set_var_model_value (s%var(ix_var), value)
    this_merit = tao_merit ()
    return
  endif
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
  tao_com%single_mode = .false.
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
  tao_com%single_mode = .false.
  do
    call tao_get_user_input (line)
    call tao_hook_command (line, found)
    if (.not. found) call tao_command (line, err)
    if (tao_com%multi_commands_here) cycle ! Use up all commands on line
    ! Keep on going if a command file is open
    ix = tao_com%cmd_file_level
    if (ix == 0 .and. .not. tao_com%cmd_file(ix)%paused) exit  
  enddo
  tao_com%single_mode = .true.

! '<CR>' Just replot.
                               
case(achar(13))  ! Ignore a <CR>

! Error:

case default
  write (*, *) 'What is this you are typing?', iachar(char)

!--------------------------------------------------------
! 'Escape' -> Must be an arrow key. Look for the rest of the sequence:

case (achar(27))

  call get_a_char (char2(1:1), .false.)
  call get_a_char (char2(2:2), .false.)

  select case (char2)

  ! '<left_arrow>' Shift key bank down by 10

  case ('[D')
    tao_com%ix_key_bank = max(tao_com%ix_key_bank-1, 0)

  ! '<right_arrow>': Shift key bank up by 10

  case ('[C')
    tao_com%ix_key_bank = min(tao_com%ix_key_bank+1, size(s%key-1)/10)

  ! '<up_arrow>': Increase deltas by factor of 10.

  case ('[A')
    s%var(:)%key_delta = s%var(:)%key_delta * 10

  ! '<down_arrow>': Decrease deltas by factor of 10.

  case ('[B')
    s%var(:)%key_delta = s%var(:)%key_delta / 10

  ! 'ESC' Error:

  case default
    write (*, *) 'What is this you are typing?', char2

  end select

!--------------------------------------------------------
! '-' Look for rest of sequence

case ('-')

  call get_a_char (char, .true., ignore_space)

  select case (char)

  ! '-p' Toggle printing.

  case ('p')
    s%global%plot_on = .not. s%global%plot_on
    if (s%global%plot_on) then
      call tao_init_plotting (' ')
      type *, 'Plotting is now: On'
    else
      type *, 'Plotting is now: Off'
    endif

  case default
    write (*, *) 'What is this you are typing?', iachar(char)

  end select

!--------------------------------------------------------
! '=' Look for rest of sequence

case ('=')

  call get_a_char (char, .true., ignore_space)

  select case (char)

  ! '=v' Set variable at given value

  case ('v')
    if (.not. allocated(s%key)) then
      call out_io (s_error$, r_name, 'THE KEY TABLE HAS NOT BEEN SET UP IN THE INIT FILES!')
      return
    endif
    write (*, '(a)', advance = "NO") 'Enter Key# and Value: '
    read (*, '(a)') line
    read (line, *, iostat = ios) ix, value
    if (ios /= 0 .or. .not. is_integer(line)) then
      call out_io (s_error$, r_name, 'I DO NOT UNDERSTAND THIS. NOTHING CHANGED.')
      return
    endif
    if (ix == 0) ix = 10
    ix2 = ix + 10*tao_com%ix_key_bank
    if (ix < 0 .or. ix > 10 .or. ix2 > size(s%key)) then
      call out_io (s_error$, r_name, 'KEY NUMBER OUT OF RANGE (0 - 9).')
      return
    endif
    call tao_set_var_model_value (s%var(s%key(ix2)), value)
  
  !----------------------------------
  ! 'Escape' -> Must be an arrow key. Look for the rest of the sequence:

  case (achar(27))

    call get_a_char (char2(1:1), .false.)
    call get_a_char (char2(2:2), .false.)

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
      do i = 1, size(s%var)
        s%var(i)%old_value = s%var(i)%model_value
      enddo

    ! '= ESC' Error:

    case default
      write (*, *) 'What is this you are typing?', char2

    end select

  case default
    write (*, *) 'What is this you are typing?', iachar(char)

  end select

!--------------------------------------------------------
! 'a' Look for rest of sequence

case ('a')

  call get_a_char (char, .true., ignore_space)

  select case (char)

  !----------------------------------
  ! 'Escape' -> Must be an arrow key. Look for the rest of the sequence:

  case (achar(27))

    call get_a_char (char2(1:1), .false.)
    call get_a_char (char2(2:2), .false.)

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

    end select

  ! 'a' Error:

  case default
    write (*, *) 'What is this you are typing?', iachar(char)

  end select

!--------------------------------------------------------
! 's' Look for rest of sequence

case ('s')

  call get_a_char (char, .true., ignore_space)

  select case (char)

  !----------------------------------
  ! 'Escape' -> Must be an arrow key. Look for the rest of the sequence:

  case (achar(27))

    call get_a_char (char2(1:1), .false.)
    call get_a_char (char2(2:2), .false.)

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

    end select

  ! 's' Error:

  case default
    write (*, *) 'What is this you are typing?', iachar(char)

  end select


!--------------------------------------------------------
! 'z' Look for rest of sequence

case ('z')

  call get_a_char (char, .true., ignore_space)

  select case (char)

  !----------------------------------
  ! 'Escape' -> Must be an arrow key. Look for the rest of the sequence:

  case (achar(27))

    call get_a_char (char2(1:1), .false.)
    call get_a_char (char2(2:2), .false.)

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

    end select

  ! 'z' Error:

  case default
    write (*, *) 'What is this you are typing?', iachar(char)

  end select


!--------------------------------------------------------
! '/' Look for rest of sequence

case ('/')

  call get_a_char (char, .true., ignore_space)

  select case (char)

! '/e' Element

  case ('e')
    write (*, '(/, a)', advance = "NO") ' Element Name or Index: '
    read (*, '(a)', iostat = ios) str
    if (ios /= 0) return
    call string_trim (str, str, ix) 
    call str_upcase (str, str)
    if (ix == 0) return
    ix = len_trim(str)

    if (str(2:2) == ':') then
      read (str(1:1), *, iostat = ios) n1
      if (n1 > ubound(s%u, 1) .or. n1 < 1) then
        write (*, *) 'ERROR: LAT INDEX NOT VALID.'
        return
      endif
      n2 = n1
      str = str(2:)
    else
      n1 = 1
      n2 = ubound(s%u, 1)
    endif

    ie = -1
    if (index ('01234567890', str(1:1)) /= 0) then
      read (str, *, iostat = ios) ie
      if (ios /= 0) ie = -1
    endif

    found = .false.
    do i = n1, n2
      lat => s%u(i)%model%lat
      do j = 0, lat%n_ele_max
        if (lat%ele(j)%name /= str .and. j /= ie) cycle
        write (*, *) '!---------------------------------------------------'
        write (*, *) '! Lat:', i, ':  ', trim(lat%lattice)
        write (*, *) '! Element Index:', j
        call type_ele (lat%ele(j), .false., 6, .true., radians$, .true., lat)
        found = .true.
      enddo
    enddo

    if (.not. found) write (*, *) 'Element not found.'

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

  ! '/p' Just create a postscript file.

  case ('p')
    call tao_output_cmd ('ps')

  ! '/P' Print a hardcopy.

  case ('P')
    call tao_output_cmd ('hard')

  ! '/u' View universe 

  case ('u')

    write (*, '(a, $)') ' Universe number to view: '
    read '(a)', line
    call string_trim (line, line, ix)
    if (ix == 0) then
      write (*, *) 'ERROR: NO UNIVERSE NUMBER.'
    else
      read (line, *, iostat = ios) i
      if (ios /= 0 .or. i < 0 .or. i > ubound(s%u, 1)) then
        write (*, *) 'ERROR: BAD UNIVERSE NUMBER.'
        return
      endif
      s%global%u_view = i
    endif

  ! '/v' Output to default file.

  case ('v')
    call tao_var_write (s%global%var_out_file)

  ! '/x' Scale horizontal axis

  case ('x')

    write (*, '(a, $)') ' Input x-axis min, max: '
    read '(a)', line
    call string_trim (line, line, ix)
    if (ix == 0) then
      call tao_x_scale_cmd ('all', 0.0_rp, 0.0_rp, err)
    else
      read (line, *, iostat = ios) m1, m2
      if (ios /= 0) then
        write (*, *) 'ERROR READING MIN/MAX.'
        return
      endif
      call tao_x_scale_cmd ('all', m1, m2, err)
    endif

  ! '/y' Scale y-axis

  case ('y')
    write (*, '(a, $)') ' Input scale min, max: '
    read '(a)', line
    call string_trim (line, line, ix)
    if (ix == 0) then
      call tao_scale_cmd ('all', 0.0_rp, 0.0_rp)
    else
      read (line, *, iostat = ios) m1, m2
      if (ios /= 0) then
        write (*, *) 'ERROR READING MIN/MAX.'
        return
      endif
      call tao_scale_cmd ('all', m1, m2)
    endif

  ! '/' Error:

  case default
    write (*, *) 'What is this you are typing?', iachar(char)

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

do i = 1, size(s%plot_region)
  if (.not. s%plot_region(i)%visible) cycle
  plot => s%plot_region(i)%plot
  if (plot%autoscale_gang_x) then
    call tao_x_scale_plot (plot, plot%x%min * factor, plot%x%max * factor)
  else
    if (.not. allocated(plot%graph)) cycle
    do j = 1, size(plot%graph)
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

real(rp) factor

!

do i = 1, size(s%plot_region)
  plot => s%plot_region(i)%plot
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

real(rp) factor, w, c

!

do i = 1, size(s%plot_region)
  plot => s%plot_region(i)%plot
  w = (plot%x%max - plot%x%min) * factor / 2
  c = (plot%x%max + plot%x%min) / 2
  call tao_x_scale_plot (plot, c-w, c+w)
enddo

end subroutine

!-----------------------------------------------------------------------
! contains

subroutine y_zoom (factor)

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

real(rp) factor, w, c

!

do i = 1, size(s%plot_region)
  plot => s%plot_region(i)%plot
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

real(rp) factor, w

!

do i = 1, size(s%plot_region)
  plot => s%plot_region(i)%plot
  w = factor * (plot%x%max - plot%x%min)
  call tao_x_scale_plot (plot, plot%x%min + w, plot%x%max + w)
enddo

end subroutine

!-----------------------------------------------------------------------
! contains

subroutine y_pan (factor)

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

real(rp) factor, w

!

do i = 1, size(s%plot_region)
  plot => s%plot_region(i)%plot
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

end subroutine   
