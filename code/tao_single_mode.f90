subroutine tao_single_mode (char)

  use tao_struct
  use tao_interface
  use tao_utils
  use quick_plot
  use single_char_input_mod
  use tao_scale_mod
  use tao_x_scale_mod
  use tao_plot_window_mod
  use tao_top10_mod

  implicit none

  type (ring_struct), pointer :: ring

  integer i, j, ix, ix_plot, ie, iv, factor, ix_key, ios, ir
  integer n1, n2, i_ele, ix_var

  real(rdef) this_factor, value, this_min, this_max, this_merit
  real(rdef) m1, m2

  real(rdef) :: delta_factor(4) = (/ 10, -10, 1, -1 /)

  character(10) :: set_char(4) = &
            (/'!@#$%^&*()', 'QWERTYUIOP', '1234567890', 'qwertyuiop' /)

  character(1) :: char, ignore_space(1) = (/ ' ' /)
  character(2) char2
  character(80) this_file, line
  character(40) :: this_opt, str, r_name = 'tao_single_mode'

  logical doit, found, err

! parse input.
! is it a variable change?

  do i = 1, 4
    if (index(set_char(i), char) /= 0) then
      ix_key = index(set_char(i), char)
      this_factor = delta_factor(i)
      ix = ix_key + s%global%ix_key_bank
      if (ix > size(s%key)) then
        call out_io (s_error$, r_name, 'KEY NOT BOUND TO VARIABLE.')
        return
      endif
      ix_var = s%key(ix)%ix_var
      if (ix_var == 0) cycle
      value = s%var(ix_var)%model_value + this_factor * s%key(ix)%delta
      call tao_set_var_model_value (s%var(ix_var), value)
      this_merit = tao_merit ()
      return
    endif
  enddo


!----------------------------------
! not a variable change so must be something else...

  select case (char)

  case ('?')
    call tao_help ('keys')

! g: Go run optimizer.

  case ('g')

    s%global%init_opt_wrapper = .true.
    s%global%optimizer_running = .true.
    do
      call tao_run_cmd (' ')
      if (.not. s%global%optimizer_running) exit
    enddo

! s: Show parameters

  case ('s')

    this_merit = tao_merit ()
    write (*, *)
    write (*, *) 's%global%optimizer:       ', trim(s%global%optimizer)
    write (*, *)
    call tao_show_constraints (0, 'ALL')
    call tao_show_constraints (0, 'TOP10')

! v: Set variable at given value

  case ('v')
    write (*, '(a, $)', advance = "NO") ' Enter Key# and Value: '
    read (*, *, iostat = ios) ix, value
    if (ios /= 0) then
      write (*, *) 'ERROR: I DO NOT UNDERSTAND THIS. NOTHING CHANGED.'
      return
    endif
    if (ix == 0) ix = 10
    if (ix < 0 .or. ix > 10) then
      write (*, *) 'ERROR KEY NUMBER OUT OF RANGE (0 - 9).'
      return
    endif
    ix = ix + s%global%ix_key_bank
    ix_var = s%key(ix)%ix_var
    call tao_set_var_model_value (s%var(ix_var), value)

! V: Show variables

  case ('V')
    call tao_var_write (' ')

! z: Quit single character input mode.

  case ('z')
    s%global%single_mode = .false.
    call out_io (s_blank$, r_name, ' ', 'Entering line mode...')

! Z: Quit tao.

  case ('Z')
    doit = .false.
    call tao_query_logical ('y', 'n', 'Quit?', doit)
    if (.not. doit) return
    call tao_destroy_plot_window
    stop

! CR: Do nothing

  case (achar(10))
    return

! <: Scale deltas by factor of 1/2

  case ('<')
    s%key(:)%delta = s%key(:)%delta / 2

! >: Scale deltas by factor of 2

  case ('>')
    s%key(:)%delta = s%key(:)%delta * 2
        
!----------------------------------
! Escape: Look for the rest of the sequence:

  case (achar(27))

    call get_a_char (char2(1:1), .false.)
    call get_a_char (char2(2:2), .false.)

  select case (char2)

! <- (left arrow): Shift key bank down by 10

  case ('[D')
    s%global%ix_key_bank = max(s%global%ix_key_bank-10, 0)

! -> (right arrow): Shift key bank up by 10

  case ('[C')
    s%global%ix_key_bank = min(s%global%ix_key_bank+10, n_key_maxx-10)

! ^ (up arrow): Scale plots by a factor of 0.5

  case ('[A')
    call scale_it_all (0, 0.5_rdef)
                                          
! v (down arrow): Scale plots by a factor of 2

  case ('[B')
    call scale_it_all (0, 2.0_rdef)

! Escape Error:

  case default
    write (*, *) 'What is this you are typing?', char2

  end select

!----------------------------------
! -: Look for rest of sequence

  case ('-')

    call get_a_char (char, .true., ignore_space)

    select case (char)

! -p: Toggle printing.

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

!----------------------------------
! /: Look for rest of sequence

  case ('/')

    call get_a_char (char, .true., ignore_space)

    select case (char)

! /e: Element

    case ('e')
      write (*, '(/, a)', advance = "NO") ' Element Name or Index: '
      read (*, '(a)', iostat = ios) str
      if (ios /= 0) return
      call string_trim (str, str, ix) 
      call str_upcase (str, str)
      if (ix == 0) return
      ix = len_trim(str)

      if (str(ix-1:ix-1) == ':') then
        read (str(ix:), *, iostat = ios) n1
        if (n1 > size(s%u) .or. n1 < 1) then
          write (*, *) 'ERROR: RING INDEX NOT VALID.'
          return
        endif
        n2 = n1
        str = str(:ix-1)
      else
        n1 = 1
        n2 = size(s%u)
      endif

      ie = -1
      if (index ('01234567890', str(1:1)) /= 0) then
        read (str, *, iostat = ios) ie
        if (ios /= 0) ie = -1
      endif

      found = .false.
      do i = n1, n2
        ring => s%u(i)%model%lat
        do j = 0, ring%n_ele_max
          if (ring%ele_(j)%name /= str .and. j /= ie) cycle
          write (*, *) '!---------------------------------------------------'
          write (*, *) '! Ring:', i, ':  ', trim(ring%lattice)
          write (*, *) '! Element Index:', j
          call type_ele (ring%ele_(j), .false., 6, &
                                               .true., radians$, .true., ring)
          found = .true.
        enddo
      enddo

      if (.not. found) write (*, *) 'Element not found.'

! /l: Lattice list

    case ('l')
      do i = 1, size(s%u)
        ring => s%u(i)%model%lat
        write (*, *)
        write (*, *) 'Ring: ', ring%lattice, i
        write (*, *) 'Ix  Name                   S  Beta_x  Beta_y'
        do j = 1, ring%n_ele_ring
          write (*, '(i3, 2x, a, f8.2, 2f8.2)') j, ring%ele_(j)%name, &
                 ring%ele_(j)%s, ring%ele_(j)%x%beta, ring%ele_(j)%y%beta
        enddo
      enddo

! /p: Just create a postscript file.

    case ('p')
      call tao_output_cmd ('ps', ' ')

! /P: Print a hardcopy.

    case ('P')
      call tao_output_cmd ('hard', ' ')
  
! /v: View universe 

    case ('v')

      write (*, '(a, $)') ' Universe number to view: '
      read '(a)', line
      call string_trim (line, line, ix)
      if (ix == 0) then
        write (*, *) 'ERROR: NO UNIVERSE NUMBER.'
      else
        read (line, *, iostat = ios) i
        if (ios /= 0 .or. i < 0 .or. i > size(s%u)) then
          write (*, *) 'ERROR: BAD UNIVERSE NUMBER.'
          return
        endif
        s%global%u_view = i
      endif


! /w: Output to default file.

    case ('w')
      call tao_var_write (s%global%var_out_file)

! /x: Scale horizontal axis

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

!----------------------------------
! /Escape: Look for the rest of the sequence:

    case (achar(27))

      call get_a_char (char2(1:1), .false.)
      call get_a_char (char2(2:2), .false.)

      select case (char2)

! /<- (/left-arrow): Paste saved variable value to variable

      case ('[D')
        do i = 1, size(s%key)
          ix_var = s%key(i)%ix_var
          call tao_set_var_model_value (s%var(ix_var), s%key(i)%val0)
        enddo

! /-> (/right-arrow): Copy variable value to saved

      case ('[C')
        do i = 1, size(s%key)
          ix = s%key(i)%ix_var
          s%key(i)%val0 = s%var(ix)%model_value
        enddo
        do i = 1, size(s%var)
          s%var(i)%old_value = s%var(i)%model_value
        enddo

! /^ (/up-arrow): Increase deltas by factor of 10.

      case ('[A')
        s%key(:)%delta = s%key(:)%delta * 10

! /v (/down-arrow): Increase deltas by factor of 10.

      case ('[B')
        s%key(:)%delta = s%key(:)%delta / 10

! /Escape Error:

      case default
        write (*, *) 'What is this you are typing?', char2

      end select

! / Error:

    case default
      write (*, *) 'What is this you are typing?', iachar(char)

    end select

!----------------------------------

! <CR>: Just replot.
                               
  case(achar(13))  ! Ignore a <CR>

! Error:

  case default
    write (*, *) 'What is this you are typing?', iachar(char)

  end select

!-----------------------------------------------------------------------
contains

subroutine scale_it_all (ix_plot, factor, a_min, a_max)

  type (tao_plot_struct), pointer :: plot
  type (tao_graph_struct), pointer :: graph

  real(rdef) factor
  real(rdef), optional :: a_min, a_max

  integer ix_plot
  character(40) name

! if factor = 0 then use min/max

  if (factor == 0) then
    call tao_scale_cmd ('all', a_min, a_max)
    return
  endif

  do i = 1, size(s%plot_page%region)
    plot => s%plot_page%region(i)%plot
    if (.not. associated (plot%graph)) cycle
    do j = 1, size(plot%graph)
      graph => plot%graph(j) 
      graph%y2%max = graph%y2%max * factor 
      graph%y2%min = graph%y2%min * factor 
      call tao_scale_graph (graph, graph%y%min * factor, graph%y%max * factor)
    enddo
  enddo

end subroutine

end subroutine   
