!+
! Subroutine tao_command (command_line, err_flag, err_is_fatal)
!
! Interface to all standard (non hook) tao commands. 
! This routine esentially breaks the command line into words
! and then calls the appropriate routine.
! Commands are case sensitive.
!
! Input:
!   command_line  -- character(*): command line
!
! Output:
!  err_flag       -- logical: Set True on error. False otherwise.
!  err_is_fatal   -- logical: Set True on non-recoverable error. False otherwise
!-

subroutine tao_command (command_line, err_flag, err_is_fatal)

use tao_set_mod, dummy2 => tao_command
use tao_change_mod, only: tao_change_var, tao_change_ele, tao_dmodel_dvar_calc, tao_change_tune, tao_change_z_tune
use tao_command_mod, only: tao_cmd_split, tao_re_execute, tao_next_switch, tao_next_word
use tao_scale_mod, only: tao_scale_cmd
use tao_wave_mod, only: tao_wave_cmd
use tao_x_scale_mod, only: tao_x_scale_cmd
use tao_plot_window_mod, only: tao_destroy_plot_window

! MPI use tao_mpi_mod

implicit none

type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (tao_plot_struct), pointer :: plt

integer i, j, n, iu, ios, n_word, n_eq, stat
integer ix, ix_line, ix_cmd, which
integer int1, int2, uni, wrt, n_level

real(rp) value1, value2, time

character(*) :: command_line
character(len(command_line)) cmd_line
character(*), parameter :: r_name = 'tao_command'
character(1000) :: cmd_word(12)
character(200) list, mask
character(40) gang_str, switch, word, except, branch_str, what
character(16) cmd_name, set_word, axis_name

character(16) :: cmd_names(48) = [character(16):: &
                      'alias', 'call', 'change', 'clear', 'clip', 'continue', 'create', 'cut_ring', 'derivative', &
                      'end_file', 'exit', 'flatten', 'help', 'json', 'ls', 'misalign', 'pause', 'pipe', 'place', &
                      'plot', 'ptc', 'python', 'quit', 're_execute', 'read', 'regression', 'reinitialize', 'reset', &
                      'restore', 'run_optimizer', 'scale', 'set', 'show', 'single_mode', 'spawn', 'taper', &
                      'timer', 'use', 'veto', 'view', 'wave', 'write', 'x_axis', 'x_scale', 'xy_scale', &
                      'debug', 'verbose', 'tree']
character(16) :: cmd_names_old(6) = [&
    'x-scale      ', 'xy-scale     ', 'single-mode  ', 'x-axis       ', 'end-file     ', &
    'output       ']

logical quit_tao, err, err_is_fatal, silent, gang, abort, err_flag, ok
logical include_wall, update, exact, include_this, listing, found

! blank line => nothing to do

err_is_fatal = .false.
err_flag = .false.

call string_trim (command_line, cmd_line, ix_line)
if (ix_line == 0 .or. cmd_line(1:1) == '!') return

! '/' denotes an option so put a space before it so it does not look like part of the command.

ix = index(cmd_line(1:ix_line), '/')
if (ix /= 0) then
  cmd_line = cmd_line(1:ix-1) // ' ' // trim(cmd_line(ix:))
  ix_line = ix - 1
endif

! strip the command line of comments

ix = index(cmd_line, '!')
if (ix /= 0) cmd_line = cmd_line(:ix-1)        ! strip off comments

! match first word to a command name

if (cmd_line(1:5) == 'quiet') then 
  call out_io (s_warn$, r_name, 'The "quiet" command has been replaced by the "set global quiet = <action>" command.')
  return
endif

call match_word (cmd_line, cmd_names, ix_cmd, .true., matched_name = cmd_name)

if (ix_cmd == 0) then  ! Accept old-style names with "-" instead of "_".
  call match_word (cmd_line, cmd_names_old, ix_cmd, .true., matched_name = cmd_name)
  ix = index(cmd_name, '-')
  if (ix /= 0) cmd_name(ix:ix) = '_'
  if (cmd_name == 'output') cmd_name = 'write'
endif

if (ix_cmd == 0) then
  call out_io (s_error$, r_name, 'UNRECOGNIZED COMMAND: ' // cmd_line)
  call tao_abort_command_file()
  return
elseif (ix_cmd < 0) then
  call out_io (s_error$, r_name, 'Ambiguous command (the "help" command will show a list of commands).')
  call tao_abort_command_file()
  return
endif

! Strip off command name from cmd_line 

call string_trim (cmd_line(ix_line+1:), cmd_line, ix_line)

! Something like "set global%rf_on" gets translated to "set global rf_on"

if (cmd_name == 'set') then
  ix = index(cmd_line, '%')
  j = index(cmd_line, ' ')
  if (ix /= 0 .and. ix < j) cmd_line(ix:ix) = ' '
endif

! select the appropriate command.

select case (cmd_name)

!--------------------------------
! ALIAS

case ('alias')

  call tao_cmd_split(cmd_line, 2, cmd_word, .false., err_flag); if (err_flag) return
  call tao_alias_cmd (cmd_word(1), cmd_word(2))
  return

!--------------------------------
! CALL

case ('call')

  call tao_cmd_split(cmd_line, 10, cmd_word, .true., err_flag); if (err_flag) goto 9000
  call tao_call_cmd (cmd_word(1), cmd_word(2:10))
  return

!--------------------------------
! CHANGE

case ('change')

  call tao_cmd_split (cmd_line, 8, cmd_word, .false., err_flag); if (err_flag) goto 9000

  silent = .false.
  update = .false.
  mask = ''
  branch_str = ''
  listing = .false.
  n = size(cmd_word)

  do i = 2, 8
    if (len_trim(cmd_word(i)) < 2) cycle

    if (index('-silent', trim(cmd_word(i))) == 1) then
      silent = .true.
      cmd_word(i:n-1) = cmd_word(i+1:n)

    elseif (index('-update', trim(cmd_word(i))) == 1) then
      update = .true.
      cmd_word(i:n-1) = cmd_word(i+1:n)

    elseif (index('-listing', trim(cmd_word(i))) == 1) then
      listing = .true.
      cmd_word(i:n-1) = cmd_word(i+1:n)

    elseif (index('-branch', trim(cmd_word(i))) == 1) then
      branch_str = cmd_word(i+1)
      cmd_word(i:n-2) = cmd_word(i+2:n)      

    elseif (index('-mask', trim(cmd_word(i))) == 1) then
      mask = cmd_word(i+1)
      cmd_word(i:n-2) = cmd_word(i+2:n)      
    endif
  enddo

  cmd_word(4) = cmd_word(4)//cmd_word(5)//cmd_word(6)//cmd_word(7)//cmd_word(8)

  if (index ('variable', trim(cmd_word(1))) == 1) then
    call tao_change_var (cmd_word(2), cmd_word(3)//cmd_word(4), silent, err_flag)

  elseif (index('element', trim(cmd_word(1))) == 1) then
    call tao_change_ele (cmd_word(2), cmd_word(3), cmd_word(4), update, err_flag)

  elseif (index('tune', trim(cmd_word(1))) == 1) then
    call tao_change_tune (branch_str, mask, listing, cmd_word(2), cmd_word(3)//cmd_word(4), err_flag)

  elseif (index('z_tune', trim(cmd_word(1))) == 1) then
    call tao_change_z_tune (branch_str, cmd_word(2)//cmd_word(3)//cmd_word(4), err_flag)

  elseif (index(trim(cmd_word(1)), 'particle_start') /= 0) then     ! Could be "2@particle_start"
    word = cmd_word(1)
    call tao_change_ele (word, cmd_word(2), cmd_word(3)//cmd_word(4), .false., err_flag)

  else
    call out_io (s_error$, r_name, 'Change who? (should be: "element", "particle_start", or "variable")')
  endif

!--------------------------------
! CLEAR

case ('clear')

  call tao_clear_cmd(cmd_line)

!--------------------------------
! CLIP

case ('clip')

  call tao_cmd_split (cmd_line, 4, cmd_word, .true., err_flag); if (err_flag) return

  gang = .false.
  if (index('-gang', trim(cmd_word(1))) == 1 .and. len_trim(cmd_word(1)) > 1) then
    gang = .true.
    cmd_word(1:3) = cmd_word(2:4)
  endif

  if (cmd_word(2) == ' ') then
    call tao_clip_cmd (gang, cmd_word(1), 0.0_rp, 0.0_rp) 
  else
    call tao_to_real (cmd_word(2), value1, err_flag);  if (err_flag) return
    if (cmd_word(3) /= ' ') then
      call tao_to_real (cmd_word(3), value2, err_flag);  if (err_flag) return
    else
      value2 = value1
      value1 = -value1
    endif
    call tao_clip_cmd (gang, cmd_word(1), value1, value2)
  endif

!--------------------------------
! CONTINUE

case ('continue')

  n_level = s%com%cmd_file_level
  if (s%com%cmd_file(n_level)%paused) then
    s%com%cmd_file(n_level)%paused = .false.
  else
    call out_io (s_error$, r_name, 'NO PAUSED COMMAND FILE HERE.')
  endif

  return

!--------------------------------
! CREATE

case ('create')
  call tao_cmd_split(cmd_line, 3, cmd_word, .false., err_flag)

  if (err_flag) then
     call out_io (s_error$, r_name, 'Error in the create command')
     return
  end if

  select case (trim(downcase(cmd_word(1))))
  case ('data')
    block
      integer, dimension(4) :: id
      integer :: jd,nd,ns
      character(:), allocatable :: ds
      character((len_trim(cmd_word(3))*11)/7+21+len_trim(cmd_word(2))) :: pipe_cmd
      type (tao_d2_data_array_struct), dimension(:), allocatable :: d2_array

      ! Check if the data exists
      call tao_find_data(err, cmd_word(2), d2_array, print_err=.false.)
      if (size(d2_array).ne.0) then
         call out_io (s_error$, r_name, 'data already exists, will not replace')
         return
      end if
      ! Count the number of d1_data
      id(1)=1
      nd = 0
      ds = trim(cmd_word(3))
      ns = len(ds)
      id(2:4) = 0
      do jd=1,ns
         select case (ds(jd:jd))
         case('[')
            id(2) = id(2)+1
         case(':')
            id(3) = id(3)+1
         case(']')
            id(4) = id(4)+1
         end select
      end do
      if (id(2).eq.0.or.id(2).ne.id(3).or.id(3).ne.id(4)) go to 70000
      nd = id(2)
      ! Start constructing the pipe command
      write(pipe_cmd,'(a,i0)') 'data_d2_create '//trim(cmd_word(2))//'^^', nd
      ! Parse the arrays
      id(1)=1
      jd = 1
      do jd=1,nd
         id(2) = index(ds(id(1):),'[') + id(1) - 1
         if (id(2).le.id(1).or.id(2).eq.ns) go to 70000
         pipe_cmd = trim(pipe_cmd)//'^^'//trim(adjustl(ds(id(1):id(2)-1)))
         id(3) = scan(ds(id(2)+1:),':') + id(2)
         if (id(3).le.id(2)+1.or.id(3).eq.ns) go to 70000
         if (.not.is_integer(ds(id(2)+1:id(3)-1))) go to 70000
         pipe_cmd = trim(pipe_cmd)//'^^'//trim(adjustl(ds(id(2)+1:id(3)-1)))
         id(4) = scan(ds(id(3)+1:),']') + id(3)
         if (id(4).le.id(3)+1) go to 70000
         if (.not.is_integer(ds(id(3)+1:id(4)-1))) go to 70000
         pipe_cmd = trim(pipe_cmd)//'^^'//trim(adjustl(ds(id(3)+1:id(4)-1)))
         if (id(4).eq.ns) exit
         if (ds(id(4)+1:id(4)+1).ne.' ') go to 70000
         id(1) = verify(ds(id(4)+1:),' ') + id(4)
      end do
      call tao_pipe_cmd(pipe_cmd)
    end block
    return
    70000 call out_io(s_error$, r_name, 'Correct form is "create data d2_name x[i:j] ..."')

  case default
     call out_io (s_error$, r_name, 'I can only create data')
  end select

  return

!--------------------------------
! CUT_RING

case ('cut_ring')

  u => tao_pointer_to_universe(-1)
  lat => u%model%lat

  what = '-static'

  call tao_cmd_split (cmd_line, 5, cmd_word, .true., err_flag); if (err_flag) goto 9000
  do i = 1, 5
    if (cmd_word(i) == '') exit
    call match_word (cmd_word(i), [character(16):: '-particle_start', '-static', '-zero'], &
                                                                       ix, .true., matched_name=switch)
    select case (switch)
    case ('-particle_start', '-static', '-zero')
      what = switch
    case default
      call out_io (s_error$, r_name, 'Unknown switch: ' // switch, 'Nothing done.')
      return
    end select
  enddo


  if (lat%param%geometry == closed$) then
    lat%param%geometry = open$
  else
    lat%param%geometry = closed$
  endif

  call out_io (s_info$, r_name, 'The lattice geometry is now: ' // geometry_name(lat%param%geometry))

  select case (what)
  case ('-static')
    if (u%model%tao_branch(0)%orbit(0)%state == alive$) then
      u%model%lat%particle_start = u%model%tao_branch(0)%orbit(0)
    endif
  case ('-zero');   u%model%lat%particle_start%vec = 0
  end select
  u%model%tao_branch(0)%orb0 = u%model%lat%particle_start

  u%calc%lattice = .true.
  call tao_lattice_calc (ok)

!--------------------------------
! DEBUG / VERBOSE / TREE

case ('debug')
  s%global%debug_on = (.not. s%global%debug_on)

case ('tree')
  s%global%expression_tree_on = (.not. s%global%expression_tree_on)
  print *, 'Using an expression tree: ', s%global%expression_tree_on

case ('verbose')
  s%global%verbose_on = (.not. s%global%verbose_on)

!--------------------------------
! DERIVATIVE

case ('derivative')

  call tao_dmodel_dvar_calc(.true., err_flag)
  call out_io (s_blank$, r_name, 'Derivative calculated')

  return

!--------------------------------
! END_FILE

case ('end_file')

  n_level = s%com%cmd_file_level
  if (n_level == 0) then
    call out_io (s_error$, r_name, 'END_FILE COMMAND ONLY ALLOWED IN A COMMAND FILE!')
    return
  endif

  call tao_close_command_file()

  if (s%com%cmd_file(n_level-1)%paused) then
    call out_io (s_info$, r_name, 'To continue the paused command file type "continue".')
  endif

  return

!--------------------------------
! EXIT/QUIT

case ('exit', 'quit')

  call string_trim (command_line, cmd_line, ix)
  if (ix < 3) then
    call out_io (s_error$, r_name, &
            'SAFETY FEATURE: YOU NEED TO TYPE AT LEAST THREE CHARACTERS TO QUIT.')
    return
  endif

  if (s%global%plot_on) call tao_destroy_plot_window
  call out_io (s_dinfo$, r_name, "Stopping.")
  !MPI !Finalize MPI if it is on
  !MPI if (s%mpi%on) call tao_mpi_finalize()
  err_is_fatal = .true. ! So Tao will stop.
  return
 
!--------------------------------
! HELP

case ('help')

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err_flag); if (err_flag) return
  call tao_help (cmd_word(1), cmd_word(2))
  return

!--------------------------------
! JSON
! This is experimental. Removal is a possibility if not developed.

case ('json')

  call tao_json_cmd (cmd_line)
  return

!--------------------------------
! LS

case ('ls')
  call system_command ('ls ' // cmd_line, err)
  return

!--------------------------------
! PAUSE

case ('pause')

  time = 0
  call tao_cmd_split (cmd_line, 1, cmd_word, .true., err_flag); if (err_flag) return
  if (cmd_word(1) /= '') then
    read (cmd_word(1), *, iostat = ios) time
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'TIME IS NOT A NUMBER.')
      return
    endif
  endif

  call tao_pause_cmd (time)
  return

!--------------------------------
! PIPE / PYTHON

case ('pipe', 'python')

  call tao_pipe_cmd (cmd_line)
  return

!--------------------------------
! PLACE

case ('place')

  call tao_cmd_split (cmd_line, 3, cmd_word, .true., err_flag); if (err_flag) return

  if (index('-no_buffer', trim(cmd_word(1))) == 1) then
    call tao_place_cmd (cmd_word(2), cmd_word(3), .true.)

  else
    if (cmd_word(3) /= ' ') then
      call out_io (s_error$, r_name, 'BAD PLACE COMMAND: ' // command_line)
      return
    endif
    call tao_place_cmd (cmd_word(1), cmd_word(2))
  endif

!--------------------------------
! PLOT
! NOTE: THIS COMMAND IS DEPRECATED 8/2021.

case ('plot')

  call out_io (s_error$, r_name, 'The "plot" command has been replaced by the "set plot <plot_name> component = ..." command.')
  return

!--------------------------------
! PTC

case ('ptc')

  call tao_cmd_split (cmd_line, 2, cmd_word, .false., err_flag); if (err_flag) goto 9000

  call tao_ptc_cmd (cmd_word(1), cmd_word(2))
  return

!--------------------------------
! RE_EXECUTE

case ('re_execute')

  call tao_re_execute (cmd_line, err_flag)
  return

!--------------------------------
! READ

case ('read')

  silent = .false.
  call tao_cmd_split (cmd_line, 5, cmd_word, .true., err_flag); if (err_flag) goto 9000
  word = ''
  do i = 1, 5
    if (cmd_word(i) == '') exit
    call match_word (cmd_word(i), [character(16):: '-universe', '-silent'], ix, .true., matched_name=switch)
    select case (switch)
    case ('-silent')
      silent = .true.
    case ('-universe')
      word = cmd_word(i+1)
      cmd_word(i:i+1) = cmd_word(i+2:i+3)
      exit
    end select
  enddo

  call tao_read_cmd (cmd_word(1), word, cmd_word(2), silent)

!--------------------------------
! REGRESSION
! This is a private, undocumented command used to produce output for use in regression testing.

case ('regression')
  call tao_regression_test()

!--------------------------------
! RESET
! This is a private, undocumented command used when debugging.

case ('reset')
  call system_command ('reset', err_flag)

!--------------------------------
! RESTORE, USE, VETO

case ('restore', 'use', 'veto')

  call tao_cmd_split(cmd_line, 2, cmd_word, .true., err_flag);  if (err_flag) goto 9000
  
  call match_word (cmd_word(1), [character(8) :: "data", "variable"], which, .true., matched_name = switch)

  select case (switch)
  case ('data')
    call tao_use_data (cmd_name, cmd_word(2))
  case ('variable')
    call tao_use_var (cmd_name, cmd_word(2))
  case default
    call out_io (s_error$, r_name, "Use/veto/restore what? data or variable?")
    return
  end select

!--------------------------------
! REINITIALIZE

case ('reinitialize')

  call tao_cmd_split(cmd_line, 2, cmd_word, .false., err_flag);  if (err_flag) goto 9000

  call match_word (cmd_word(1), ['data', 'tao ', 'beam'], ix, .true., matched_name=word)

  select case (word)

  case ('beam') 
    do i = lbound(s%u, 1), ubound(s%u, 1)
      s%u(i)%model_branch(:)%beam%init_starting_distribution = .true.
      s%u(i)%calc%lattice = .true.
    enddo

  case ('data') 
    s%u(:)%calc%lattice = .true.

  case ('tao') 
    call tao_parse_command_args (err, cmd_word(2));  if (err_flag) goto 9000

    if (s%init%init_file_arg /= '') call out_io (s_info$, r_name, 'Reinitializing with: ' // s%init%init_file_arg)
    call tao_init (err_flag)
    return

  case default
    call out_io (s_error$, r_name, 'Reinit what? Choices are: "beam", "data", or "tao".')
    return
    
  end select

!--------------------------------
! RUN, FLATTEN

case ('run_optimizer', 'flatten')

  call tao_cmd_split (cmd_line, 1, cmd_word, .true., err_flag); if (err_flag) goto 9000
  call tao_run_cmd (cmd_word(1), abort)

!--------------------------------
! SCALE

case ('scale')

  call tao_cmd_split (cmd_line, 7, cmd_word, .true., err_flag); if (err_flag) return

  axis_name = ''
  gang_str = ''
  include_wall = .false.
  exact = .false.

  i = 1
  do
    if (cmd_word(i) == '') exit
    call match_word (cmd_word(i), [character(16):: '-y', '-y2', '-nogang', '-gang', '-include_wall', '-exact'], &
                                                                                     ix, .true., matched_name=switch)
    select case (switch)
    case ('-exact');            exact = .true.
    case ('-y', '-y2');         axis_name = switch(2:)
    case ('-gang', '-nogang');  gang_str = switch(2:)
    case ('-include_wall');     include_wall = .true.
    case default;               i = i + 1;  cycle
    end select

    cmd_word(i:i+6) = cmd_word(i+1:i+7)
  enddo

  if (cmd_word(2) == ' ') then
    call tao_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, axis_name, include_wall, gang_str)
  else
    call tao_to_real (cmd_word(2), value1, err_flag);  if (err_flag) return
    if (cmd_word(3) /= ' ') then
      call tao_to_real (cmd_word(3), value2, err_flag);  if (err_flag) return
    else
      value2 = value1
      value1 = -value1
    endif
    call tao_scale_cmd (cmd_word(1), value1, value2, axis_name, include_wall, gang_str, exact)
  endif

!--------------------------------
! SET

case ('set')
  update = .false.
  set_word = ''
  branch_str = ''
  mask = ''
  listing = .false.
  silent = .false.

  do
    ! "-1" is a universe index and not a switch.
    if (cmd_line(1:1) == '-' .and. cmd_line(1:2) /= '-1') then
      call tao_next_switch (cmd_line, [character(20) :: '-update', '-lord_no_set', '-mask', &
                                    '-branch', '-listing', '-silent'], .true., switch, err_flag)
      if (err_flag) return
      select case (switch)
      case ('-update')
        update = .true.
      case ('-listing')
        listing = .true.
      case ('-lord_no_set')
        call out_io (s_warn$, r_name, 'Note: The "-lord_no_set" no longer exists. This set will be ignored.')
      case ('-branch')
        call tao_next_word(cmd_line, branch_str)
      case ('-mask')
        call tao_next_word(cmd_line, mask)
      case ('-silent')
        silent = .true.
      end select
      cycle
    endif

    if (set_word /= '') exit

    call tao_next_switch (cmd_line, [character(20) :: 'branch', 'data', 'var', 'lattice', &
      'universe', 'curve', 'graph', 'beam_init', 'wave', 'plot', 'bmad_com', 'element', 'opti_de_param', &
      'csr_param', 'floor_plan', 'lat_layout', 'geodesic_lm', 'default', 'key', 'particle_start', &
      'plot_page', 'ran_state', 'symbolic_number', 'beam', 'beam_start', 'dynamic_aperture', &
      'global', 'region', 'calculate', 'space_charge_com', 'ptc_com', 'tune', 'z_tune'], .true., switch, err_flag)
    if (err_flag) return
    set_word = switch
  enddo

  select case (set_word)
  case ('csr_param')
    call out_io (s_warn$, r_name, '"csr_param" structure is now called "space_charge_com"')
    goto 9000
  case ('ran_state'); n_word = 2; n_eq = 1
  case ('beam', 'beam_init', 'bmad_com', 'space_charge_com', 'data', 'global', 'lattice', 'default', &
        'opti_de_param', 'wave', 'floor_plan', 'lat_layout', 'geodesic_lm', 'key', 'symbolic_number', &
        'var', 'beam_start', 'particle_start', 'dynamic_aperture', 'ptc_com'); n_word = 3; n_eq = 2
  case ('universe'); n_word = 4; n_eq = 3
  case ('plot_page'); n_word = 4; n_eq = 2
  case ('branch', 'curve', 'element', 'graph', 'plot', 'region'); n_word = 4; n_eq = 3
  case ('calculate'); n_word = 1; n_eq = 0
  case ('tune'); n_word = 7; n_eq = 0
  case ('z_tune'); n_word = 6; n_eq = 0
  case default
    call out_io (s_error$, r_name, 'SET WHAT? (MUST BE ON OF "branch", "data", "var", ...etc.')
    goto 9000
  end select

  ! Split command line into words. Translate "set ele [1,2]@q[k1]" -> "set ele [1,2]@q k1"

  n = size(cmd_word)
  call tao_cmd_split (cmd_line, n_word, cmd_word, .false., err, '=')

  if (set_word == 'tune' .or. set_word == 'z_tune') then
    j = 1
    do i = 1, n_word
      if (index('-mask', cmd_word(i)) == 1 .and. len_trim(cmd_word(i)) > 1) then
        mask = cmd_word(i+1)
        cmd_word(i:n-2) = cmd_word(i+2:n)
      elseif (index('-branch', cmd_word(i)) == 1 .and. len_trim(cmd_word(i)) > 1) then
        branch_str= cmd_word(i+1)
        cmd_word(i:n-2) = cmd_word(i+2:n)
      elseif (index('-listing', cmd_word(i)) == 1 .and. len_trim(cmd_word(i)) > 1) then
        listing = .true.
        cmd_word(i:n-1) = cmd_word(i+1:n)
      else
        j = j + 1
      endif
    enddo
  endif

  if  (set_word == 'element' .and. index('-update', trim(cmd_word(1))) == 1 .and. len_trim(cmd_word(1)) > 1) then
    update = .true.
    call tao_cmd_split (cmd_line, 5, cmd_word, .false., err, '=')
    cmd_word(1:4) = cmd_word(2:5)
  endif

  ix = str_last_in_set(cmd_word(1), '[')
  if (set_word == 'element' .and. ix /= 0 .and. ix > index(cmd_word(1), '@')) then
    n = len_trim(cmd_word(1)) 
    if (cmd_word(1)(n:n) /= ']') then
      call out_io (s_error$, r_name, 'CANNOT DECODE: ' // cmd_word(1))
      goto 9000
    endif
    cmd_word(3:5) = cmd_word(2:4)
    cmd_word(2) = cmd_word(1)(ix+1:n-1)
    cmd_word(1) = cmd_word(1)(1:ix-1)
  endif

  !

  if (set_word == 'universe' .and. cmd_word(3) /= '=') then  ! Old syntax
    cmd_word(4) = cmd_word(3)
    cmd_word(3) = '='
  endif

  if (n_eq > 0) then
    if (cmd_word(n_eq) /= '=') then
      call out_io (s_error$, r_name, 'SYNTAX PROBLEM. "=" NOT IN CORRECT PLACE.')
      goto 9000
    endif
  endif

  select case (set_word)
  case ('beam')
    call tao_set_beam_cmd (cmd_word(1), unquote(cmd_word(3)), branch_str)
  case ('beam_init')
    call tao_set_beam_init_cmd (cmd_word(1), cmd_word(3), branch_str)
  case ('beam_start', 'particle_start')
    if (set_word == 'beam_start') call out_io (s_warn$, r_name, 'Note: "beam_start" is now named "particle_start".')
    call tao_set_particle_start_cmd (cmd_word(1), cmd_word(3))
  case ('bmad_com')
    call tao_set_bmad_com_cmd (cmd_word(1), cmd_word(3))
  case ('branch')
    call tao_set_branch_cmd (cmd_word(1), cmd_word(2), cmd_word(4)) 
  case ('calculate')
    call tao_set_calculate_cmd (cmd_word(1))
  case ('curve')
    call tao_set_curve_cmd (cmd_word(1), cmd_word(2), cmd_word(4)) 
  case ('data')
    call tao_set_data_cmd (cmd_word(1), cmd_word(3), silent)
  case ('default')
    call tao_set_default_cmd (cmd_word(1), cmd_word(3))
  case ('dynamic_aperture')
    call tao_set_dynamic_aperture_cmd (cmd_word(1), cmd_word(3))
  case ('element')
    call tao_set_elements_cmd (cmd_word(1), cmd_word(2), cmd_word(4), update)
  case ('floor_plan')
    call tao_set_drawing_cmd (s%plot_page%floor_plan, cmd_word(1), cmd_word(3))
  case ('geodesic_lm')
    call tao_set_geodesic_lm_cmd (cmd_word(1), cmd_word(3))
  case ('global')
    call tao_set_global_cmd (cmd_word(1), cmd_word(3))
  case ('graph')
    call tao_set_graph_cmd (cmd_word(1), cmd_word(2), cmd_word(4))
  case ('key')
    call tao_set_key_cmd (cmd_word(1), cmd_word(3))    
  case ('lat_layout')
    call tao_set_drawing_cmd (s%plot_page%lat_layout, cmd_word(1), cmd_word(3))
  case ('lattice')
    call tao_set_lattice_cmd (cmd_word(1), cmd_word(3))
  case ('opti_de_param')
    call tao_set_opti_de_param_cmd (cmd_word(1), cmd_word(3))
  case ('plot ')
    call tao_set_plot_cmd (cmd_word(1), cmd_word(2), cmd_word(4))
  case ('plot_page')
    call tao_set_plot_page_cmd (cmd_word(1), cmd_word(3), cmd_word(4))
  case ('ran_state')
    call tao_set_ran_state_cmd (cmd_word(2))
  case ('ptc_com')
    call tao_set_ptc_com_cmd (cmd_word(1), cmd_word(3))
  case ('region')
    call tao_set_region_cmd (cmd_word(1), cmd_word(2), cmd_word(4))
  case ('space_charge_com')
    call tao_set_space_charge_com_cmd (cmd_word(1), cmd_word(3))
  case ('symbolic_number')
    call tao_set_symbolic_number_cmd(cmd_word(1), cmd_word(3))
  case ('tune')
    if (cmd_word(1) == '=') cmd_word(1:2) = cmd_word(2:3)
    call tao_set_tune_cmd (branch_str, mask, listing, cmd_word(1), cmd_word(2), .false.)
  case ('universe')    
    call tao_set_universe_cmd (cmd_word(1), cmd_word(2), cmd_word(4))
  case ('var')
    call tao_set_var_cmd (cmd_word(1), cmd_word(3))
  case ('wave')
    call tao_set_wave_cmd (cmd_word(1), cmd_word(3), err_flag);  if (err_flag) goto 9000
    call tao_cmd_end_calc
    call tao_show_cmd ('wave')
  case ('z_tune')
    if (cmd_word(1) == '=') cmd_word(1:2) = cmd_word(2:3)
    call tao_set_z_tune_cmd (branch_str, cmd_word(1), .false.)
  end select

!--------------------------------
! SHOW

case ('show')

  call tao_show_cmd (cmd_line)
  return

!--------------------------------
! SINGLE-MODE

case ('single_mode')

  if (cmd_line /= '') then
    call out_io (s_error$, r_name, 'Extra stuff on line: ' // cmd_line)
    return
  endif

  if (s%global%plot_on) then
    found = .false.
    do i = 1, size(s%plot_page%region)
      if (.not. s%plot_page%region(i)%visible) cycle
      plt => s%plot_page%region(i)%plot
      if (.not. allocated(plt%graph)) cycle
      do j = 1, size(plt%graph)
        if (plt%graph(j)%type == 'key_table') found = .true.
      enddo
    enddo
    if (.not. found) call out_io(s_blank$, r_name, '[Note: Use the "place" command in line mode if you want to see the key_table.]')
  endif

  s%com%single_mode = .true.
  call out_io (s_blank$, r_name, 'Entering Single Mode. Waiting for your input...')

  return

!--------------------------------
! SPAWN

case ('spawn')

  call system_command (cmd_line, err_flag)
  if (err_flag) call tao_abort_command_file()
  return

!--------------------------------
! taper

case ('taper')

  except = ''
  word = ''

  call tao_cmd_split (cmd_line, 4, cmd_word, .true., err_flag); if (err_flag) return

  i = 0
  do
    i = i + 1
    if (cmd_word(i) == '') exit
    call match_word (cmd_word(i), [character(20):: '-universe', '-except'], ix, .true., matched_name=switch)

    select case (switch)
    case ('-except')
      i = i + 1
      except = cmd_word(i)
    case ('-universe')
      i = i + 1
      word = cmd_word(i)
    case default
      call out_io (s_error$, r_name, 'UNKNOWN SWITCH: ' // cmd_word(1))
      return
    end select
  enddo

  call tao_taper_cmd(except, word)
  call tao_cmd_end_calc
  return

!--------------------------------
! timer

case ('timer')

  call tao_timer (cmd_line)
  return

!--------------------------------
! view

case ('view')
  call tao_set_default_cmd ('universe', cmd_line)
  call tao_cmd_end_calc
  return 

!--------------------------------
! wave

case ('wave')

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err_flag); if (err_flag) return
  call tao_wave_cmd (cmd_word(1), cmd_word(2), err_flag); if (err_flag) return
  call tao_cmd_end_calc
  call tao_show_cmd ('wave')
  return

!--------------------------------
! write

case ('write')

  call tao_write_cmd (cmd_line)
  return

!--------------------------------
! X_AXIS

case ('x_axis')

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err_flag); if (err_flag) return
  call tao_x_axis_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! X_SCALE

case ('x_scale')

  call tao_cmd_split (cmd_line, 5, cmd_word, .true., err_flag); if (err_flag) return

  gang_str = ''
  include_wall = .false.
  exact = .false.

  i = 1
  do
    if (cmd_word(i) == '') exit
    call match_word (cmd_word(i), [character(16):: '-nogang', '-gang', '-include_wall', '-exact'], &
                                                                               ix, .true., matched_name=switch)

    select case (switch)
    case ('-exact');            exact = .true.
    case ('-gang', '-nogang');  gang_str = switch(2:)
    case ('-include_wall');     include_wall = .true.
    case default;               i = i + 1;  cycle
    end select

    cmd_word(i:i+5) = cmd_word(i+1:i+6)
  enddo

  if (cmd_word(2) == ' ') then
    call tao_x_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, err, include_wall, gang_str)
  else
    call tao_to_real (cmd_word(2), value1, err_flag); if (err_flag) return
    call tao_to_real (cmd_word(3), value2, err_flag); if (err_flag) return
    call tao_x_scale_cmd (cmd_word(1), value1, value2, err, include_wall, gang_str, exact)
  endif

!--------------------------------
! XY_SCALE

case ('xy_scale')

  call tao_cmd_split (cmd_line, 5, cmd_word, .true., err_flag); if (err_flag) return

  include_wall = .false.
  exact = .false.

  i = 1
  do
    if (cmd_word(i) == '') exit
    call match_word (cmd_word(i), [character(16):: '-include_wall', '-exact'], ix, .true., matched_name=switch)

    select case (switch)
    case ('-exact');            exact = .true.
    case ('-include_wall');     include_wall = .true.
    case default;               i = i + 1;  cycle
    end select

    cmd_word(i:i+5) = cmd_word(i+1:i+6)
  enddo


  if (cmd_word(2) == ' ') then
    call tao_x_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, err, include_wall = include_wall)
    call tao_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, include_wall = include_wall) 
  else
    call tao_to_real (cmd_word(2), value1, err_flag);  if (err_flag) return
    if (cmd_word(3) /= ' ') then
      call tao_to_real (cmd_word(3), value2, err_flag);  if (err_flag) return
    else
      value2 = value1
      value1 = -value1
    endif
    call tao_x_scale_cmd (cmd_word(1), value1, value2, err, include_wall = include_wall, exact = exact)
    call tao_scale_cmd (cmd_word(1), value1, value2, include_wall = include_wall, exact = exact)
  endif

!--------------------------------
! DEFAULT

case default

  call out_io (s_error$, r_name, 'INTERNAL COMMAND PARSING ERROR!')
  call err_exit

end select

!------------------------------------------------------------------------
! Do the standard calculations and plotting after command
! Note: wave command bypasses this.

call tao_cmd_end_calc
return

!------------------------------------------------------------------------
! Error:

9000 continue
call tao_abort_command_file()

end subroutine tao_command




