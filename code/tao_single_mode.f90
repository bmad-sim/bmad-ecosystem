subroutine tao_single_mode (char)

  use tao_mod
  use quick_plot
  use tao_single_mod
  use single_char_input_mod

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

  logical doit, found

! parse input.
! is it a variable change?

  do i = 1, 4
    if (index(set_char(i), char) /= 0) then
      ix_key = index(set_char(i), char)
      this_factor = delta_factor(i)
      ix = ix_key + s%global%ix_key_bank
      ix_var = s%key(ix)%ix_var
      if (ix_var == 0) cycle
      value = s%var(ix_var)%model_value + this_factor * s%key(ix)%delta
      call tao_set_var_model_value (s%var(ix_var), value)
      call tao_merit ()
      return
    endif
  enddo


!----------------------------------
! not a variable change so must be something else...

  select case (char)

  case ('?')
    call fullfilename ('TAO_DIR:doc/tao_single_mode.doc', this_file)
    call type_this_file (this_file)

! z: Quit single character input mode.

  case ('z')
    print *
    s%global%single_mode = .false.

! Z: Quit tao.

  case ('Z')
    doit = .false.
    call tao_query_logical ('y', 'n', 'Quit?', doit)
    if (.not. doit) return
    call qp_close_page
    stop

! v: Set variable at given value

  case ('v')
    print '(a, $)', ' Enter Key# and Value: '
    read (*, *, iostat = ios) ix, value
    if (ios /= 0) then
      print *, 'ERROR: I DO NOT UNDERSTAND THIS. NOTHING CHANGED.'
      return
    endif
    if (ix == 0) ix = 10
    if (ix < 0 .or. ix > 10) then
      print *, 'ERROR KEY NUMBER OUT OF RANGE (0 - 9).'
      return
    endif
    ix = ix + s%global%ix_key_bank
    ix_var = s%key(ix)%ix_var
    call tao_set_var_model_value (s%var(ix_var), value)

! g: Go run optimizer.

  case ('g')

    s%global%init_opt_wrapper = .true.
    s%global%optimizer_running = .true.
    do
      call tao_run_cmd (' ')
      if (.not. s%global%optimizer_running) exit
    enddo

! G: Go run the optimizer.

  case ('G')

    print '(a, $)', ' Optimizer to use: '
    read '(a)', this_opt
    call str_upcase (this_opt, this_opt)
    call string_trim (this_opt, this_opt, ix)

    if (ix == 0) then
      print *, 'No name given. Nothing done.'
      return
    endif

    if (this_opt(1:1) == '@') then
      s%global%optimizer = this_opt(2:)
      this_opt = this_opt(2:)
    endif 

    s%global%optimizer_running = .true.
    s%global%init_opt_wrapper = .true.
    do
      call tao_run_cmd (this_opt)
      if (.not. s%global%optimizer_running) exit
    enddo

! s: Show parameters

  case ('s')

    this_merit = tao_merit ()
    print *
    print *, 's%global%optimizer:       ', trim(s%global%optimizer)
    print *
    call tao_show_constraints (0, 'ALL')
    call tao_show_constraints (0, 'TOP10')

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
    print *, 'What is this your typing?', char2

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
      print *, 'What is this your typing?', iachar(char)

    end select

!----------------------------------
! /: Look for rest of sequence

  case ('/')

    call get_a_char (char, .true., ignore_space)

    select case (char)

! /e: Element

    case ('e')
      print '(/, a, $)', ' Element Name or Index: '
      read (*, '(a)', iostat = ios) str
      if (ios /= 0) return
      call string_trim (str, str, ix) 
      call str_upcase (str, str)
      if (ix == 0) return
      ix = len_trim(str)

      if (str(ix-1:ix-1) == ':') then
        read (str(ix:), *, iostat = ios) n1
        if (n1 > size(s%u) .or. n1 < 1) then
          print *, 'ERROR: RING INDEX NOT VALID.'
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
        ring => s%u(i)%model
        do j = 0, ring%n_ele_max
          if (ring%ele_(j)%name /= str .and. j /= ie) cycle
          print *, '!---------------------------------------------------'
          print *, '! Ring:', i, ':  ', trim(ring%lattice)
          print *, '! Element Index:', j
          call type_ele (ring%ele_(j), .false., 6, &
                                               .true., radians$, .true., ring)
          found = .true.
        enddo
      enddo

      if (.not. found) print *, 'Element not found.'

! /h: Define header title

    case ('h')
      print '(a, $)', ' Header title: '
      read '(a)', s%plot_page%title(1)%string

! /l: Lattice list

    case ('l')
      do i = 1, size(s%u)
        ring => s%u(i)%model
        print *
        print *, 'Ring: ', ring%lattice, i
        ring => s%u(i)%model
        print *, 'Ix  Name                   S  Beta_x  Beta_y'
        do j = 1, ring%n_ele_ring
          print '(i3, 2x, a, f8.2, 2f8.2)', j, ring%ele_(j)%name, &
                 ring%ele_(j)%s, ring%ele_(j)%x%beta, ring%ele_(j)%y%beta
        enddo
      enddo

! /p: Just create a postscript file.

    case ('p')
      call tao_output_cmd ('ps')

! /P: Print a hardcopy.

    case ('P')
      call tao_output_cmd ('hard')
  
! /t: Read from the current tao input file

    case ('t')
      call tao_read_single_input_file (s%global%tao_single_mode_file)
      print *, 'Read in tao file: ', trim(s%global%tao_single_mode_file)

! /T: Read from a tao input file

    case ('T')
      print '(a, $)', ' tao input file name: '
      read '(a)', this_file
      call string_trim (this_file, this_file, ix)
      if (ix == 0) then
        print *, 'No name given. Nothing done.'
        return
      endif
      if (this_file(1:1) == '@') then
        s%global%tao_single_mode_file = this_file(2:)
        this_file = this_file(2:)
      endif
      call tao_read_single_input_file (this_file)
      print *, 'Read in tao file: ', trim(this_file)

! /v: Output variable settings

    case ('v')
!!      call var_print (.true., .true.)

! /w: Output to default file.

    case ('w')
      call tao_var_write ('tao#.manual', .true.)

! /W: Output to a file.

    case ('W')

      do i = 1, size(s%u)

        print '(1x, i1, a, $)', i, ': Output file name: '
        read '(a)', this_file
        call string_trim (this_file, this_file, ix)

        if (ix == 0) then
          print *, 'No name given. Nothing done.'
          return
        endif
        call tao_var_write (this_file, .true.)

      enddo

! /x: Scale horizontal axis

    case ('x')

      print '(a, $)', ' Input x-axis min, max: '
      read '(a)', line
      call string_trim (line, line, ix)
      if (ix == 0) then
        call tao_x_scale_cmd ('all', 0, 1e20_rp)
      else
        read (line, *, iostat = ios) m1, m2
        if (ios /= 0) then
          print *, 'ERROR READING MIN/MAX.'
          return
        endif
        call tao_x_scale_cmd ('all', m1, m2)
      endif

! /number: Scale a single plot

    case ('1', '2', '3', '4', '5', '6', '7', '8', '9')

      read (char, *) ix_plot
      call get_a_char (char, .true., ignore_space)

      select case (char)

! /number:: Get min, max

      case (':')

        print '(a, $)', ' Input plot min, max: '
        read (*, *, iostat = ios) this_min, this_max
        if (ios /= 0) then
          print *, 'ERROR READING MIN/MAX'
          return
        endif
        call scale_it_all (ix_plot, 0.0_rdef, this_min, this_max)

! /number Escape: Look for the rest of the sequence:

      case (achar(27))

        call get_a_char (char2(1:1), .false.)
        call get_a_char (char2(2:2), .false.)

        select case (char2)

! /number ^ (up arrow): scale plot

        case ('[A')
          call scale_it_all (ix_plot, 0.5_rdef)

! /number ^ (up arrow): scale plot

        case ('[B') 
          call scale_it_all (ix_plot, 2.0_rdef)

! /number Escape Error:

        case default
          print *, 'What is this your typing?', char2
  
        end select      ! /number Escape

      end select      ! /number

!----------------------------------
! /Escape: Look for the rest of the sequence:

    case (achar(27))

      call get_a_char (char2(1:1), .false.)
      call get_a_char (char2(2:2), .false.)

      select case (char2)

! /<- (left arrow): Paste saved variable value to variable

      case ('[D')
        do i = 1, size(s%key)
          ix_var = s%key(i)%ix_var
          call tao_set_var_model_value (s%var(ix_var), s%key(i)%val0)
        enddo

! /-> (right arrow): Copy variable value to saved

      case ('[C')
        do i = 1, size(s%key)
          ix = s%key(i)%ix_var
          s%key(i)%val0 = s%var(ix)%model_value
        enddo
        do i = 1, size(s%var)
          s%var(i)%old_value = s%var(i)%model_value
        enddo

! /Escape Error:

      case default
        print *, 'What is this your typing?', char2

      end select

! / Error:

    case default
      print *, 'What is this your typing?', iachar(char)

    end select

!----------------------------------

! <CR>: Just replot.
                               
  case(achar(13))  ! Ignore a <CR>

! Error:

  case default
    print *, 'What is this your typing?', iachar(char)

  end select

!-----------------------------------------------------------------------
contains

subroutine scale_it_all (ix_plot, factor, a_min, a_max)

  type (tao_graph_struct), pointer :: graph

  real(rdef) factor
  real(rdef), optional :: a_min, a_max

  integer ix_plot

! if factor = 0 then use min/max

  if (factor == 0) then
    call tao_scale_cmd ('all', a_min, a_max)
    return
  endif

  do i = 1, size(s%plot_page%plot)
    do j = 1, size(s%plot_page%plot(i)%graph)
      graph => s%plot_page%plot(i)%graph(j) 
      graph%y%max  = graph%y%max * factor 
      graph%y%min  = graph%y%min * factor 
      graph%y2%max = graph%y2%max * factor 
      graph%y2%min = graph%y2%min * factor 
    enddo
  enddo

end subroutine

end subroutine   
