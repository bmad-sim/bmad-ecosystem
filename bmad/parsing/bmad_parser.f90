!+
! Subroutine bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag, parse_lat)
!
! Subroutine to parse a BMAD input file and put the information in lat.
!
! Because of the time it takes to parse a file bmad_parser will save 
! LAT in a "digested" file with the name:
!               <lat_file>.digested_NNN   
! For subsequent calls to the same lat_file, BMAD_PARSER will just read in the
! digested file. bmad_parser will always check to see that the digested file
! is up-to-date and if not the digested file will not be used.
!
! Input:
!   lat_file   -- character(*): Name of the input file.
!   make_mats6 -- logical, optional: Compute the 6x6 transport matrices for the Elements?
!                   Default is True. Do not set False unless you know what you are doing.
!   use_line   -- character(*), optional: If present and not blank, override the use 
!                   statement in the lattice file and use use_line instead.
!
! Output:
!   lat              -- lat_struct: Lat structure. See bmad_struct.f90 for more details.
!     %ele(:)%mat6      -- This is computed assuming an on-axis orbit if make_mats6 = T.
!   digested_read_ok -- logical, optional: Set True if the digested file was
!                        successfully read. False otherwise.
!   err_flag         -- logical, optional: Set true if there is an error, false otherwise.
!                         Note: err_flag does *not* include errors in lat_make_mat6 since
!                         if there is a match element, there is an error raised since
!                         the Twiss parameters have not been set but this is expected. 
!   parse_lat        -- lat_struct, optional: List of elements used to construct the lattice.
!                         Useful if bmad_parser2 will be called. See bmad_parser2 documentation.
!-

subroutine bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag, parse_lat)

use parser_set_attribute_mod, dummy1 => bmad_parser
use wall3d_mod, dummy3 => bmad_parser
use photon_target_mod, dummy4 => bmad_parser
use ptc_interface_mod, only: set_ptc_com_pointers
use random_mod

implicit none

type (lat_struct), target :: lat, in_lat
type (lat_struct), optional :: parse_lat
type (ele_struct) this_ele
type (ele_struct), pointer :: ele, slave, lord, ele2, ele0, param_ele
type (ele_struct), save :: marker_ele
type (seq_struct), target, allocatable :: sequence(:)
type (branch_struct), pointer :: branch0, branch
type (parser_lat_struct), target :: plat
type (parser_ele_struct), pointer :: pele
type (lat_ele_loc_struct) m_slaves(100)
type (ele_pointer_struct), allocatable :: fork_ele_list(:), eles(:)
type (parser_controller_struct), allocatable :: pcon(:)
type (base_line_ele_struct), allocatable :: a_line(:)
type (seq_ele_struct), pointer :: s_ele

real(rp) beta, val

integer, allocatable :: seq_indexx(:)

integer :: ix_param_ele, ix_mad_beam_ele
integer ix_word, i_use, i, j, k, k2, n, ix, ix1, ix2, n_track
integer n_ele_use, digested_version, key, loop_counter, n_ic, n_con
integer iseq_tot, iyy, n_ele_max, n_multi, n0, n_ele, ixc, n_slave, n_loc
integer ib, ie, ib2, ie2, flip, n_branch, n_forks, i_loop, n_branch_max
integer, pointer :: n_max

character(*) lat_file
character(*), optional :: use_line

character(1) delim
character(16), parameter :: r_name = 'bmad_parser'
character(40) word_1, word_2, name, this_name, this_branch_name
character(40), allocatable :: seq_name(:)
character(80) debug_line
character(400) full_lat_file_name, digested_file, call_file
character(280) parse_line_save, line, use_line_str

logical, optional :: make_mats6, digested_read_ok, err_flag
logical delim_found, arg_list_found, wild_here
logical end_of_file, ele_found, match_found, err, finished, exit_on_error
logical multipass, heterogeneous_ele_list
logical auto_bookkeeper_saved, is_photon_fork, created_new_branch

! See if digested file is open and current. If so read in and return.
! Note: The name of the digested file depends upon the real precision.

call cpu_time(bp_com%time0)
call init_bmad()
auto_bookkeeper_saved = bmad_com%auto_bookkeeper
bmad_com%auto_bookkeeper = .false.  

if (present(err_flag)) err_flag = .true.
bp_com = bp_common_struct()
allocate(bp_com%lat_file_names(1))       !! To get around an ifort bug in Versions 18+
bp_com%parser_name = 'bmad_parser'       ! Used for error messages.
debug_line = ''
err = .false.

if (lat_file == '') then
  call parser_error('lattice file name is blank!')
  if (global_com%exit_on_error) then
    call out_io (s_fatal$, r_name, 'BMAD_PARSER FINISHED. EXITING ON ERRORS')
    stop
  endif
  bp_com%parser_name = ''
  return
endif

if (.not. bp_com%always_parse) then
  call form_digested_bmad_file_name (lat_file, digested_file, full_lat_file_name, use_line)
  call read_digested_bmad_file (digested_file, lat, digested_version, err, .true., bp_com%lat_file_names)
  bp_com%num_lat_files = size(bp_com%lat_file_names)
endif

! Must make sure that if use_line is present the digested file has used the 
! correct line

if (present(use_line)) then
  if (use_line /= '') then
    call str_upcase (name, use_line)
    if (name /= lat%use_name) err = .true.
  endif
endif

if (.not. err .and. .not. bp_com%always_parse) then
  if (ptc_private%taylor_order_ptc /= 0 .and. lat%input_taylor_order /= 0 .and. &
                                             lat%input_taylor_order /= ptc_private%taylor_order_ptc) then
     call out_io (s_info$, r_name, 'Taylor_order has changed.', &
           'Taylor_order in digested file: \i4\ ', &
           'Taylor_order now:              \i4\ ', &
           'Will now set to the new Taylor order...', &
           i_array = [lat%input_taylor_order, ptc_private%taylor_order_ptc])
    if (lat%input_taylor_order > ptc_private%taylor_order_ptc) bp_com%write_digested = .false.

  else
    if (present(digested_read_ok)) digested_read_ok = .true.
    call parser_end_stuff (in_lat)
    return
  endif
endif

if (present(digested_read_ok)) digested_read_ok = .false.

! here if not OK bmad_status. So we have to do everything from scratch...
! init variables.

call init_lat (lat, 1)
call init_lat (in_lat, 1000)

call allocate_plat (plat, ubound(in_lat%ele, 1))

do i = 0, ubound(in_lat%ele, 1)
  in_lat%ele(i)%ixx = i   ! Pointer to plat%ele() array
enddo

call out_io (s_info$, r_name, 'Parsing lattice file(s). This might take a minute or so...')
call parser_file_stack('init')
call parser_file_stack('push', lat_file, finished, err)  ! open file on stack
if (err) then
  call parser_end_stuff (in_lat)
  return
endif

iseq_tot = 0                            ! number of sequences encountered
allocate(sequence(500))

call init_bmad_parser_common()
bp_com%extra = extra_parsing_info_struct()
bp_com%input_line_meaningful = .true.

n_max => in_lat%n_ele_max   ! Number of elements encountered
n_max = -1

! Note: The order of def_parameter and def_mad_beam elements is used by parser_set_attribute

n_max = n_max + 1
ele => in_lat%ele(n_max)
call init_ele(ele, beginning_ele$, 0, n_max, in_lat%branch(0))
ele%name = 'BEGINNING'
call set_ele_defaults (ele)   ! Defaults for beginning_ele element
call nametable_add (in_lat%nametable, ele%name, n_max)

n_max = n_max + 1
ele => in_lat%ele(n_max) ! Important: def_parameter must come after def_mad_beam due to overlapping parameters.
call init_ele(ele, def_mad_beam$, 0, n_max, in_lat%branch(0))
ele%name = 'BEAM'                 ! For MAD compatibility.
call nametable_add (in_lat%nametable, ele%name, n_max)
ix_mad_beam_ele = 1

n_max = n_max + 1
ele => in_lat%ele(n_max)  ! Important: def_parameter comes after def_mad_beam.
call init_ele(ele, def_parameter$, 0, n_max, in_lat%branch(0))
ele%name = 'PARAMETER'           ! For parameters 
call nametable_add (in_lat%nametable, ele%name, n_max)
ix_param_ele = 2

n_max = n_max + 1
ele => in_lat%ele(n_max)
call init_ele(ele, def_particle_start$, 0, n_max, in_lat%branch(0))
ele%name = 'PARTICLE_START'           ! For beam starting parameters 
call nametable_add (in_lat%nametable, ele%name, n_max) 

n_max = n_max + 1
ele => in_lat%ele(n_max)
call init_ele(ele, def_ptc_com$, 0, n_max, in_lat%branch(0))
ele%name = 'PTC_COM'           ! Global PTC parameters
call nametable_add (in_lat%nametable, ele%name, n_max)

n_max = n_max + 1
ele => in_lat%ele(n_max)
call init_ele(ele, def_bmad_com$, 0, n_max, in_lat%branch(0))
ele%name = 'BMAD_COM'           ! Global bmad parameters
call nametable_add (in_lat%nametable, ele%name, n_max)

n_max = n_max + 1
ele => in_lat%ele(n_max)
call init_ele(ele, def_space_charge_com$, 0, n_max, in_lat%branch(0))
ele%name = 'SPACE_CHARGE_COM'           ! Space charge parameters
call nametable_add (in_lat%nametable, ele%name, n_max)

lat%n_control_max = 0

call load_parse_line ('init', 1, end_of_file)

!-----------------------------------------------------------
! main parsing loop

loop_counter = 0  ! Used for debugging
parsing_loop: do 

  loop_counter = loop_counter + 1

  ! get a line from the input file and parse out the first word.
  call load_parse_line ('new_command', 1, end_of_file)  ! load an input line
  call get_next_word (word_1, ix_word, '[:](,)= ', delim, delim_found, .true.)
  if (end_of_file) then
    word_1 = 'END_FILE'
    ix_word = 8
  endif

  ! If input line is something like "quadrupole::*[k1] = ..." then shift delim from ":" to "["

  if (delim == ':' .and. bp_com%parse_line(1:1) == ':') then
    ix = index(bp_com%parse_line, '[')
    if (ix /= 0) then
      word_1 = trim(word_1) // ':' // upcase(bp_com%parse_line(:ix-1))
      bp_com%parse_line = bp_com%parse_line(ix+1:)
      delim = '['
      ix_word = len_trim(word_1)
    endif
  endif

  ! Name check. 
  ! If delim = '[' then have attribute redef and things are complicated so do not check.

  if (delim /= '[') then
    if (.not.  verify_valid_name(word_1, ix_word)) cycle
  endif

  !-------------------------------------------
  ! PARSER_DEBUG

  if (word_1(:ix_word) == 'PARSER_DEBUG') then
    debug_line = bp_com%parse_line
    bp_com%parse_line = ''
    call out_io (s_info$, r_name, 'Found in file: "PARSER_DEBUG". Debug is now on')
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! PRINT

  if (word_1(:ix_word) == 'PRINT') then
    call parser_print_line(lat, end_of_file)  ! Put directly in lat.
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! USE_LOCAL_LAT_FILE

  if (word_1(:ix_word) == 'USE_LOCAL_LAT_FILE') then
    bp_com%use_local_lat_file = .true.
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! NO_DIGESTED

  if (word_1(:ix_word) == 'NO_DIGESTED') then
    if (bp_com%write_digested) call out_io (s_info$, r_name, 'Found in file: "NO_DIGESTED". No digested file will be created')
    bp_com%write_digested = .false.
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! WRITE_DIGESTED

  if (word_1(:ix_word) == 'WRITE_DIGESTED') then
    if (.not. bp_com%write_digested) call out_io (s_info$, r_name, 'Found in file: "WRITE_DIGESTED". A digested file will be created')
    bp_com%write_digested = .true.
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! NO_SUPERIMPOSE

  if (word_1(:ix_word) == 'NO_SUPERIMPOSE') then
    bp_com%do_superimpose = .false.
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! DEBUG_MARKER is used to be able to easily set a break within the debugger

  if (word_1(:ix_word) == 'DEBUG_MARKER') then
    word_1 = 'ABC'          ! An executable line to set a break on
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! SLICE_LATTICE, etc.

  select case (word_1(:ix_word))
  case ('SLICE_LATTICE', 'COMBINE_ELEMENTS', 'CALC_REFERENCE_ORBIT', 'START_BRANCH_AT', &
        'COMBINE_CONSECUTIVE_ELEMENTS', 'MERGE_ELEMENTS', 'REMOVE_ELEMENTS')
    call parser_error ('A ' // word_1(:ix_word) // ' COMMAND CAN ONLY BE USED AFTER AN EXPAND_LATTICE COMMAND.')
    cycle parsing_loop
  end select

  !-------------------------------------------
  ! Superimpose statement

  if (word_1(:ix_word) == 'SUPERIMPOSE') then
    call new_element_init('superimpose-command:' // int_str(n_max+1), '', in_lat, plat, err, null_ele$)
    ele => in_lat%ele(n_max)
    call parse_superimpose_command(in_lat, ele, plat%ele(ele%ixx), delim)
    cycle parsing_loop   
  endif

  !-------------------------------------------
  ! USE command...

  if (word_1(:ix_word) == 'USE') then
    lat%use_name = ''
    do
      if (delim /= ',' .and. delim /= ' ') call parser_error ('MISSING COMMA IN "USE" STATEMENT.')
      call get_next_word(word_2, ix_word, ':(=,)', delim, delim_found, .true.)
      if (ix_word == 0) then 
        call parser_error ('CONFUSED "USE" STATEMENT', '')
        cycle parsing_loop
      endif
      if (.not. verify_valid_name(word_2, ix_word)) cycle parsing_loop
      lat%use_name = trim(lat%use_name) // ',' // word_2
      if (.not. delim_found .and. bp_com%parse_line == '') then
        lat%use_name = lat%use_name(2:)   ! Trim initial comma
        cycle parsing_loop
      endif
    enddo
  endif

  !-------------------------------------------
  ! TITLE command

  if (word_1(:ix_word) == 'TITLE') then
    if (delim_found) then
      if (delim /= " " .and. delim /= ",") call parser_error ('BAD DELIMITOR IN "TITLE" COMMAND')
      call bmad_parser_string_attribute_set(this_ele, 'DESCRIP', delim, delim_found)
      lat%title = this_ele%descrip
      deallocate (this_ele%descrip)
    else
      lat%title = bp_com%next_chunk
      bp_com%next_chunk = ''
    endif
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! CALL command

  if (word_1(:ix_word) == 'CALL') then
    call get_called_file(delim, call_file, err)
    if (err) then
      call parser_end_stuff (in_lat)
      return
    endif

    cycle parsing_loop
  endif

  !-------------------------------------------
  ! BEAM command

  if (word_1(:ix_word) == 'BEAM') then
    if (delim /= ',') call parser_error ('"BEAM" NOT FOLLOWED BY COMMA')
    do 
      if (.not. delim_found) exit
      if (delim /= ',') then
        call parser_error ('EXPECTING: "," BUT GOT: ' // delim, 'FOR "BEAM" COMMAND')
        exit
      endif
      call parser_set_attribute (def$, in_lat%ele(ix_mad_beam_ele), delim, delim_found, err)
      if (bp_com%fatal_error_flag) exit parsing_loop
    enddo
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! EXPAND_LATTICE command

  if (word_1(:ix_word) == 'EXPAND_LATTICE') then
    bp_com%detected_expand_lattice_cmd = .true.
    exit parsing_loop
  endif

  !-------------------------------------------
  ! RETURN or END_FILE command

  if (word_1(:ix_word) == 'RETURN' .or.  word_1(:ix_word) == 'END_FILE') then
    call parser_file_stack ('pop', ' ', finished, err)
    if (err) then
      call parser_end_stuff (in_lat)
      return
    endif
    if (finished) exit parsing_loop ! break loop
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! If an element attribute redef

  if (delim == '[') then

    call get_next_word (word_2, ix_word, '][=', delim, delim_found, .true.)
    if (delim /= ']') then
      call parser_error ('OPENING "[" FOUND WITHOUT MATCHING "]"')
      cycle parsing_loop
    endif

    call get_next_word (this_name, ix_word, ':=', delim, delim_found, .true.)
    if (.not. delim_found .or. ix_word /= 0) then
      call parser_error ('MALFORMED ELEMENT ATTRIBUTE REDEFINITION')
      cycle parsing_loop
    endif

    ! If delim is ':' then this is an error since get_next_word treats
    ! a ':=' construction as a '=' 

    if (delim == ':') then
      call parser_error ('MALFORMED ELEMENT ATTRIBUTE REDEF')
      cycle parsing_loop
    endif

    ! Find associated element and evaluate the attribute value...

    if (word_1(1:1) == '+') then 
      word_1 = word_1(2:)
      call parser_error ('"+" sign prefix construct deprecated. Please remove it.', level = s_warn$)
    endif

    ixc = index(word_1, '::')
    wild_here = (index(word_1, '*') /= 0 .or. index(word_1, '%') /= 0) ! Wild card character found
    key = key_name_to_key_index(word_1)
    parse_line_save = trim(word_2) // ' = ' // bp_com%parse_line 

    ! If just a name then we can look this up

    if (word_1 == 'BEAM_START') then
      word_1 = 'PARTICLE_START'   ! For backwards compatibility
      call parser_error ('Beam_start has been renamed particle_start.', &
                         'Note: If you are doing beam tracking, particle_start will now not affect', &
                         'the beam centroid unless beam_init%use_particle_start = True.', level = s_warn$)
    endif

    if (ixc == 0 .and. key == -1 .and. .not. wild_here) then    
      call find_index (word_1, in_lat%nametable, ix)
      if (ix == -1) then
        if (index(name, '##') /= 0) then
          call parser_error ('"ELEMENT##N" CONSTRUCT NOT VALID BEFORE AN "EXPAND_LATTICE" COMMAND: ' // name)
        elseif (is_integer(name)) then
          call parser_error ('USING AN ELEMENT INDEX NOT VALID BEFORE AN "EXPAND_LATTICE" COMMAND: ' // name)
        else
          call parser_error ('ELEMENT NOT FOUND: ' // word_1)
        endif
      else
        ele => in_lat%ele(ix)
        bp_com%parse_line = parse_line_save
        call parser_set_attribute (redef$, ele, delim, delim_found, err, plat%ele(ele%ixx))
        if (bp_com%fatal_error_flag) exit parsing_loop
        if (.not. err .and. delim_found) call parser_error ('BAD DELIMITER: ' // delim)
      endif
      bp_com%parse_line = ''  ! Might be needed in case of error.
      cycle parsing_loop
    endif

    ! Not a simple name so have to loop over all elements and look for a match

    if (key > -1) then   ! If Old style "quadrupole[k1] = ..." syntax
      name = '*'
      wild_here = .true.

    elseif (ixc == 0) then   ! Simple element name: "q01w[k1] = ..."
      key = 0
      name = word_1

    else                    ! "key::name" syntax
      key = key_name_to_key_index(word_1(1:ixc-1), .true.)
      name = word_1(ixc+2:)
      if (name == '') then
        call parser_error('ELEMENT NAME IN "CLASS::NAME" SYNTAX CANNOT BE BLANK', &
                          'NOTE: USE "*" TO MATCH TO ALL NAMES')
        cycle parsing_loop
      endif
      if (key == -1) then
        bp_com%parse_line = ''
        call parser_error ('BAD ELEMENT CLASS NAME: ' // word_1(1:ixc-1))
        cycle parsing_loop
      endif
    endif

    ! When setting an attribute for many elements (EG "1:10[k1] = ..."), only some of the elements may have
    ! the attribute to be set. This is acceptible and should not generate an error.

    ele_found = .false.
    heterogeneous_ele_list = (key == 0 .and. wild_here)

    do i = 0, n_max
      ele => in_lat%ele(i)
      if (key /= 0 .and. ele%key /= key) cycle
      ! No wild card matches permitted for these.
      if (ele%key == beginning_ele$ .or. ele%key == def_mad_beam$ .or. &
          ele%key == def_parameter$ .or. ele%key == def_particle_start$ .or. &
          ele%key == def_bmad_com$ .or. ele%key == def_space_charge_com$) cycle
      if (.not. match_wild(ele%name, trim(name))) cycle
      ! 
      if (heterogeneous_ele_list .and. attribute_index(ele, word_2) < 1) cycle
      bp_com%parse_line = parse_line_save
      ele_found = .true.
      call parser_set_attribute (redef$, ele, delim, delim_found, err, plat%ele(ele%ixx), heterogeneous_ele_list = heterogeneous_ele_list)
      if (bp_com%fatal_error_flag) exit parsing_loop
      if (err .or. delim_found) then
        if (.not. err .and. delim_found) call parser_error ('BAD DELIMITER: ' // delim)
        bp_com%parse_line = '' 
        cycle parsing_loop
      endif
    enddo

    bp_com%parse_line = '' ! Needed if last call to parser_set_attribute did not have a set.

    if (index(name, '##') /= 0) then
      call parser_error ('"ELEMENT##N" CONSTRUCT NOT VALID BEFORE AN "EXPAND_LATTICE" COMMAND: ' // name)
    elseif (index(name, '>>') /= 0) then
      call parser_error ('AN ELEMENT NAME WHICH CONTAINS ">>" NOT VALID BEFORE AN "EXPAND_LATTICE" COMMAND: ' // name, &
                         'THE REASON FOR THIS IS THAT LATTICE BRANCHES ARE NOT FORMED UNTIL LATTICE EXPANSION.')
    elseif (.not. ele_found .and. .not. wild_here) then
      if (index(name, ':') /= 0) then
        call parser_error ('"ELEMENT1:ELEMENT2" CONSTRUCT NOT VALID BEFORE AN "EXPAND_LATTICE" COMMAND: ' // name)
      else
        call parser_error ('ELEMENT NOT FOUND: ' // name)
      endif
    elseif (.not. ele_found .and. attribute_index (key, word_2) == 0) then
      call parser_error ('BAD ATTRIBUTE')
      bp_com%parse_line = '' 
    endif

    cycle parsing_loop
  endif

  !-------------------------------------------
  ! Variable def

  if (word_1 == 'REDEF' .and. delim == ':') then
    call get_next_word (word_1, ix_word, '[:](,)= ', delim, delim_found, .true.)
    call parser_add_constant (word_1, in_lat, .false.)
    cycle parsing_loop
  endif

  if (delim == '=') then
    call parser_add_constant (word_1, in_lat, .true.)
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! If a "(" delimitor then we are looking at a replacement line.

  if (delim == '(') then
    call get_sequence_args (word_1, sequence(iseq_tot+1)%dummy_arg, delim, err)
    ix = size(sequence(iseq_tot+1)%dummy_arg)
    allocate (sequence(iseq_tot+1)%corresponding_actual_arg(ix))
    if (err) cycle parsing_loop
    arg_list_found = .true.
    call get_next_word (word_2, ix_word, '(): =,', delim, delim_found, .true.)
    if (word_2 /= ' ') call parser_error &
                ('":" NOT FOUND AFTER REPLACEMENT LINE ARGUMENT LIST. ' // &
                'FOUND: ' // word_2, 'FOR LINE: ' // word_1)
  else
    arg_list_found = .false.
  endif

  ! must have a ":" delimiter now

  if (delim /= ':') then
    if (delim == '') then
      call parser_error ('UNRECOGNIZED COMMAND: ' // word_1)
    else
      call parser_error ('EXPECTED DELIMITER TO BE ":". BUT IT IS: "' // trim(delim),  '" AFTER WORD: ' // word_1)
    endif
    cycle parsing_loop
  endif

  ! only possibilities left are: element, list, or line
  ! to decide which look at 2nd word

  call get_next_word(word_2, ix_word, ':=,', delim, delim_found, .true.)
  if (ix_word == 0) then
    call parser_error ('NOTHING FOUND AFTER: ' // word_1)
    call parser_end_stuff (in_lat)
    return
  endif

  if (word_2 == 'LINE[MULTIPASS]') then
    word_2 = 'LINE'
    ix_word = 4
    multipass = .true.
  else
    multipass = .false.
  endif

  if (.not. verify_valid_name(word_2, ix_word)) cycle parsing_loop

  ! arg lists are only used with lines

  if (word_2(:ix_word) /= 'LINE' .and. arg_list_found) then
    call parser_error ('ARGUMENTS "XYZ(...):" ARE ONLY USED WITH REPLACEMENT LINES.', &
                                                      'FOR: ' // word_1)
    cycle parsing_loop
  endif

  !-------------------------------------------------------
  ! If line or list

  if (word_2(:ix_word) == 'LINE' .or. word_2(:ix_word) == 'LIST') then
    iseq_tot = iseq_tot + 1
    if (iseq_tot > size(sequence)-1) call reallocate_sequence(sequence, 2*iseq_tot)

    sequence(iseq_tot)%name = word_1
    sequence(iseq_tot)%multipass = multipass

    call new_element_init (word_1, '', in_lat, plat, err)
    ele => in_lat%ele(n_max)

    if (delim /= '=') call parser_error ('EXPECTING: "=" BUT GOT: ' // delim)
    if (word_2(:ix_word) == 'LINE') then
      if (arg_list_found) then
        sequence(iseq_tot)%type = replacement_line$
        ele%key = replacement_line$
      else
        sequence(iseq_tot)%type = line$
        ele%key = def_line$
        call set_ele_defaults (ele)
      endif
    else
      sequence(iseq_tot)%type = list$
      ele%key = list$
    endif
    call parse_line_or_list (sequence, iseq_tot, in_lat, .true.)

    cycle parsing_loop
  endif

  !-------------------------------------------------------
  ! If not line or list then must be an element

  call new_element_init (word_1, word_2, in_lat, plat, err)
  if (err) cycle parsing_loop

  ! Check for valid element key name or if element is part of a element key.
  ! If none of the above then we have an error.

  match_found = .false.  ! found a match?

  call find_index (word_2, in_lat%nametable, i)
  if (i >= 0 .and. i < n_max) then ! i < n_max avoids "abc: abc" construct.
    plat%ele(n_max) = plat%ele(i)
    in_lat%ele(n_max) = in_lat%ele(i)
    in_lat%ele(n_max)%ixx = n_max  ! Restore correct value
    call set_ele_name (in_lat%ele(n_max), word_1)
    match_found = .true.
  endif

  if (.not. match_found .and. word_2 == 'RF') then
    call parser_error ('"RF" ELEMENT TYPE IS AMBIGUOUS. COULD BE "RFCAVITY" OR "RF_BEND".', &
                       'WILL ASSUME THIS IS A "RFCAVITY"', level = s_warn$)
    word_2 = 'RFCAVITY'
  endif

  if (.not. match_found) then
    in_lat%ele(n_max)%key = key_name_to_key_index(word_2, .true.)
    if (key_name_to_key_index(word_1, .false.) > 0) then
      call parser_error ('ELEMENT NAME: ' // word_1, &
                         'IS NOT ALLOWED TO BE THE SAME AS AN ELEMENT CLASS: ' // word_2)
    endif

    if (in_lat%ele(n_max)%key > 0) then
      call set_ele_defaults (in_lat%ele(n_max))
      match_found = .true.
    endif
  endif

  if (.not. match_found) then
    if (index(word_2, '[') /= 0) then
      call parser_error ('"ELEMENT1:ELEMENT2" CONSTRUCT NOT VALID BEFORE AN "EXPAND_LATTICE" COMMAND: ' // name)
    else
      call parser_error ('KEY NAME NOT RECOGNIZED OR AMBIGUOUS: ' // word_2,  &
                         'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
    endif
    cycle parsing_loop
  endif

  ! Element definition...
  ! First parse any overlay/group/girder slave list

  key = in_lat%ele(n_max)%key
  if (key == overlay$ .or. key == group$ .or. key == girder$ .or. key == ramper$) then
    if (delim /= '=') then
      call parser_error ('EXPECTING: "=" BUT GOT: ' // delim,  &
                         'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
      cycle parsing_loop        
    endif

    select case (key)
    case (overlay$); in_lat%ele(n_max)%lord_status = overlay_lord$
    case (group$);   in_lat%ele(n_max)%lord_status = group_lord$
    case (girder$);  in_lat%ele(n_max)%lord_status = girder_lord$
    case (ramper$);  in_lat%ele(n_max)%lord_status = ramper_lord$
    end select

    call get_overlay_group_names(in_lat%ele(n_max), in_lat, plat%ele(n_max), delim, delim_found, .false., err)
    if (err) cycle parsing_loop

    if (key /= girder$ .and. .not. delim_found) then
      call parser_error ('NO CONTROL ATTRIBUTE GIVEN AFTER CLOSING "}"',  &
                    'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
      cycle parsing_loop
    endif
  endif

  if (key == feedback$) in_lat%ele(n_max)%lord_status = control_lord$

  ! Now loop over all attributes...

  do 
    if (.not. delim_found) exit   ! If no delim then we are finished with this element.
    if (delim /= ',') then
      call parser_error ('EXPECTING: "," BUT GOT: ' // delim,  'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
      exit
    endif
    call parser_set_attribute (def$, in_lat%ele(n_max), delim, delim_found, err, plat%ele(n_max))
    if (bp_com%fatal_error_flag) exit parsing_loop
    if (err) cycle parsing_loop
  enddo

enddo parsing_loop       ! main parsing loop

if (bp_com%fatal_error_flag) then
  call parser_end_stuff (in_lat)
  return
endif

!===========================================================================
!---------------------------------------------------------------------------
! We now have read in everything. 

bp_com%input_line_meaningful = .false.
param_ele    => in_lat%ele(ix_param_ele)

! Sort lists and check for duplicates.

allocate (seq_indexx(iseq_tot), seq_name(iseq_tot))
seq_name = sequence(1:iseq_tot)%name
call indexer (seq_name, seq_indexx)

do i = 1, iseq_tot-1
  ix1 = seq_indexx(i)
  ix2 = seq_indexx(i+1)
  if (sequence(ix1)%name == sequence(ix2)%name) call parser_error  ('DUPLICATE LINE NAME ' // sequence(ix1)%name)
enddo

do i = 1, in_lat%nametable%n_max-1
  ix1 = in_lat%nametable%index(i)
  ix2 = in_lat%nametable%index(i+1)
  if (in_lat%ele(ix1)%name == in_lat%ele(ix2)%name) call parser_error ('DUPLICATE ELEMENT NAME ' // in_lat%ele(ix1)%name)
enddo

! Now to expand the lines and lists to find the elements to use.
! First go through the lines and lists and index everything.

do k = 1, iseq_tot
  do i = 1, size(sequence(k)%ele(:))

    s_ele => sequence(k)%ele(i)
    name = s_ele%name

    if (s_ele%ix_arg > 0) then   ! dummy arg
      s_ele%type = element$
      cycle
    endif

    ! Remember: Sequence names also appear in the element list so search the sequence list first.

    call find_index (name, seq_name, seq_indexx, iseq_tot, j)
    if (j == 0) then  ! if not a sequence, it must be an element
      call find_index (name, in_lat%nametable, j)
      if (j < 0) then  ! if not an element, I don't know what it is
        s_ele%ix_ele = -1       ! Struggle on for now...
        s_ele%type = element$
      else
        s_ele%ix_ele = j
        s_ele%type = element$
      endif
    else
      s_ele%ix_ele = j
      s_ele%type = sequence(j)%type
      if (s_ele%type == list$ .and. s_ele%ele_order_reflect) call parser_error ( &
                          'A REFLECTION WITH A LIST IS NOT ALLOWED IN: '  &
                          // sequence(k)%name, 'FOR LIST: ' // s_ele%name, seq = sequence(k))
      if (sequence(k)%type == list$) &
                call parser_error ('A REPLACEMENT LIST: ' // sequence(k)%name, &
                'HAS A NON-ELEMENT MEMBER: ' // s_ele%name)
    endif

  enddo
enddo

!----------------------------------------------------------------------
! Expand all branches...

bp_com%used_line_set_by_calling_routine = .false.

if (present (use_line)) then
  if (use_line /= '') then
    call str_upcase (name, use_line)
    if (lat%use_name /= '' .and. name /= lat%use_name) bp_com%used_line_set_by_calling_routine = .true.
    lat%use_name = name
  endif
endif

if (lat%use_name == blank_name$) then
  if (iseq_tot == 1) then
    call parser_error ('NO "USE" STATEMENT FOUND.', &
                       'However since there is only one line, that will be used.', level = s_warn$)
    lat%use_name = sequence(1)%name
  else
    call parser_error ('NO "USE" STATEMENT FOUND.', 'I DO NOT KNOW WHAT LINE TO USE!')
    call parser_end_stuff (in_lat)
    return
  endif
endif

use_line_str = lat%use_name
n_forks = 0
allocate (fork_ele_list(20))

n_branch_max = 1000
branch_loop: do i_loop = 1, n_branch_max

  ! Expand branches from fork elements before expanding branches from the use command. 

  if (n_forks == 0) then
    if (use_line_str == '') exit
    ix = index(use_line_str, ',')
    if (ix == 0) then
      this_branch_name = use_line_str
      use_line_str = ''
    else
      this_branch_name = use_line_str(1:ix-1)
      use_line_str = use_line_str(ix+1:)
    endif

    call parser_expand_line (1, this_branch_name, sequence, seq_name, seq_indexx, &
                                                    is_true(param_ele%value(no_end_marker$)), n_ele_use, lat, in_lat)
    if (bp_com%fatal_error_flag) then
      call parser_end_stuff (in_lat)
      return
    endif
    is_photon_fork = .false.

  else
    ele => fork_ele_list(1)%ele
    call parser_add_branch (ele, lat, sequence, seq_name, seq_indexx, &
                       is_true(param_ele%value(no_end_marker$)), in_lat, plat, created_new_branch, this_branch_name)
    is_photon_fork = (ele%key == photon_fork$)
    n_forks = n_forks - 1
    fork_ele_list(1:n_forks) = fork_ele_list(2:n_forks+1)
    if (.not. created_new_branch) cycle 
  endif

  n_branch = ubound(lat%branch, 1)
  branch => lat%branch(n_branch)

  call find_index (this_branch_name, in_lat%nametable, ix)
  ele => in_lat%ele(ix) ! line_ele element associated with this branch.

  ele0 => branch%ele(0)
  ele0%value(e_tot$) = -1
  ele0%value(p0c$)   = -1 
  call settable_dep_var_bookkeeping(ele0)

  ! Add energy, species, etc info for all branches except branch(0) which is handled "old style".

  if (n_branch == 0) then
    lat%ele(0)                  = in_lat%ele(0)    ! Beginning element
    lat%ele(0)%orientation      = lat%ele(1)%orientation
    lat%version                 = bmad_inc_version$
    lat%input_file_name         = full_lat_file_name             ! save input file  
    lat%particle_start          = in_lat%particle_start
    lat%a                       = in_lat%a
    lat%b                       = in_lat%b
    lat%z                       = in_lat%z
    lat%photon_type             = in_lat%photon_type
    lat%input_taylor_order      = in_lat%input_taylor_order
    if (allocated(in_lat%custom))   lat%custom   = in_lat%custom
    if (allocated(in_lat%constant)) lat%constant = in_lat%constant

    call mat_make_unit (lat%ele(0)%mat6)
    call clear_lat_1turn_mats (lat)

    lat%param%n_part = param_ele%value(n_part$)

    lat%param%particle = param_ele%ref_species
    if (ele%ref_species /= not_set$) lat%param%particle = ele%ref_species
    if (lat%param%particle == not_set$) lat%param%particle = positron$

    ! The lattice name from a "parameter[lattice] = ..." line is 
    ! stored the param_ele%descrip string

    lat%lattice = in_lat%lattice
    lat%machine = in_lat%machine
    
    ! Set live_branch.

    val = param_ele%value(live_branch$)
    if (val /= real_garbage$) then  ! live_branch has been set.
      lat%param%live_branch = is_true(val)
    endif

    ! Set geometry.

    val = param_ele%value(geometry$)
    if (val /= real_garbage$) then  ! geometry has been set.
      lat%param%geometry = nint(val)
    elseif (lat%param%particle == photon$) then
      lat%param%geometry = open$
    endif

  else  ! n_branch /= 0
    branch%param%particle = ele%ref_species
  endif  

  !----

  if (ele%value(live_branch$) /= real_garbage$) branch%param%live_branch = is_true(ele%value(live_branch$))
  if (ele%value(high_energy_space_charge_on$) /= real_garbage$ .and. is_true(ele%value(high_energy_space_charge_on$))) bmad_com%high_energy_space_charge_on = .true.
  if (ele%value(geometry$) /= real_garbage$) branch%param%geometry = nint(ele%value(geometry$))

  ! Transfer info from line element if parameters have been set.

  branch%param%default_tracking_species = ref_particle$
  val = param_ele%value(default_tracking_species$)
  if (n_branch == 0 .and.  val /= real_garbage$) branch%param%default_tracking_species = nint(val)

  do i = 1, num_ele_attrib$
    if (ele%value(i) == real_garbage$) cycle
    select case (i)
    case (default_tracking_species$);  branch%param%default_tracking_species = nint(ele%value(i))
    case default;                      ele0%value(i) = ele%value(i)
    end select
  enddo

  call set_this_real_val (ele%s,                 ele0%s)
  call set_this_real_val (ele%ref_time,          ele0%ref_time)
  call set_this_real_val (ele%floor%r(1),        ele0%floor%r(1))
  call set_this_real_val (ele%floor%r(2),        ele0%floor%r(2))
  call set_this_real_val (ele%floor%r(3),        ele0%floor%r(3))
  call set_this_real_val (ele%floor%theta,       ele0%floor%theta)
  call set_this_real_val (ele%floor%phi,         ele0%floor%phi)
  call set_this_real_val (ele%floor%psi,         ele0%floor%psi)
  call set_this_twiss_struct (ele%a, ele0%a)
  call set_this_twiss_struct (ele%b, ele0%b)
  call set_this_twiss_struct (ele%z, ele0%z)
  call set_this_xy_disp_struct (ele%x, ele0%x)
  call set_this_xy_disp_struct (ele%y, ele0%y)

  call floor_angles_to_w_mat(ele0%floor%theta, ele0%floor%phi, ele0%floor%psi, ele0%floor%w)

  !

  call create_lat_ele_nametable(lat, lat%nametable)

  if (bp_com%error_flag) then
    call parser_end_stuff (in_lat)
    return
  endif

  ! Go through the IN_LAT elements and put in the superpositions.
  ! If the superposition is a branch element, need to add the branch line.

  call s_calc (lat)              ! calc longitudinal distances
  call control_bookkeeper (lat)

  do i = 1, n_max
    if (in_lat%ele(i)%lord_status /= super_lord$) cycle
    pele => plat%ele(i)
    if (pele%superposition_command_here) then
      call lat_ele_locator (pele%ele_name, in_lat, eles, n_loc, err)
      if (n_loc == 0) then
        call parser_error ('CANNOT FIND ELEMENT FOR SUPERPOSITION: ' // pele%ele_name)
        cycle
      endif
      call parser_add_superimpose (branch, eles(1)%ele, pele, in_lat, plat)
    else
      call parser_add_superimpose (branch, in_lat%ele(i), pele, in_lat, plat)
    endif
    call s_calc (lat)  ! calc longitudinal distances of new branch elements
  enddo

  call adjust_super_slave_names (lat, lat%n_ele_track+1, lat%n_ele_max)

  ! For bookkeeping purposes, null_ele elements with %sub_key = drift$ were created from drifts that were superimposed upon.
  ! Change these to %sub_key = 0 so they will be ignored if parser_add_superimpose is called for other lattice branches.

  do i = lat%n_ele_track+1, lat%n_ele_max
    if (lat%ele(i)%key == null_ele$ .and. lat%ele(i)%sub_key == drift$) lat%ele(i)%sub_key = 0
  enddo

  ! Add branch lines to the list of branches to construct.

  j = 0
  do i = 1, branch%n_ele_max
    if (branch%ele(i)%key /= photon_fork$ .and. branch%ele(i)%key /= fork$) cycle
    j = j + 1
    n_forks = n_forks + 1
    call re_allocate_eles (fork_ele_list, n_forks + 10, .true., .false.)
    fork_ele_list(j+1:n_forks) = fork_ele_list(j:n_forks-1)
    fork_ele_list(j)%ele => branch%ele(i)
  enddo

  if (i_loop == n_branch_max) then
    call parser_error ('1000 BRANCHES GENERATED. LOOKS LIKE AN ENDLESS LOOP')
    call parser_end_stuff (in_lat)
    return
  endif 

enddo branch_loop

! Work on multipass...
! Multipass elements are paired by ele%iyy tag and ele%name must both match.

do ib = 0, ubound(lat%branch, 1)
  do ie = 1, lat%branch(ib)%n_ele_max
    ele => lat%branch(ib)%ele(ie)
    if (ele%iyy == 0) cycle
    if (ele%slave_status == super_slave$) cycle
    if (ele%key == null_ele$) cycle
    n_multi = 0  ! number of elements to slave together
    iyy = ele%iyy
    do ib2 = ib, ubound(lat%branch, 1)
      do ie2 = 1, lat%branch(ib2)%n_ele_max
        if (ib == ib2 .and. ie2 < ie) cycle
        ele2 => lat%branch(ib2)%ele(ie2)
        if (ele2%iyy /= iyy) cycle
        if (ele2%name /= ele%name) cycle
        n_multi = n_multi + 1
        m_slaves(n_multi) = ele_loc (ele2)
        ele2%iyy = 0  ! mark as taken care of
      enddo
    enddo
    call add_this_multipass (lat, m_slaves(1:n_multi))
  enddo
enddo

call drift_multipass_name_correction(lat)

!-------------------------------------
! If a girder element refers to a line then must expand that line.

do i = 1, n_max
  lord => in_lat%ele(i)
  if (lord%key /= girder$) cycle
  pele => plat%ele(i)
  n_slave = size(pele%control)
  j = 0
  do 
    j = j + 1
    if (j > n_slave) exit
    call find_index(pele%control(j)%name, seq_name, seq_indexx, size(seq_name), k, k2)
    if (k == 0) cycle
    call parser_expand_line (-1, pele%control(j)%name, sequence, &
                                         seq_name, seq_indexx, .false., n_ele_use, expanded_line = a_line)
    if (bp_com%fatal_error_flag) then
      call parser_end_stuff (in_lat)
      return
    endif

    ! Put elements from the line expansion into the slave list.

    call move_alloc(pele%control, pcon)
    allocate (pele%control(n_slave+n_ele_use-1))

    pele%control(1:j-1) = pcon(1:j-1)
    pele%control(j:j+n_ele_use-1)%name = a_line(1:n_ele_use)%name
    pele%control(j+n_ele_use:n_slave+n_ele_use-1) = pcon(j:n_slave-1)

    j = j + n_ele_use - 1
  enddo
enddo

! PTC stuff.
! Use arbitrary energy above the rest mass energy since when tracking individual elements the
! true reference energy is used.

if (lat%input_taylor_order /= 0) ptc_private%taylor_order_saved = lat%input_taylor_order
call set_ptc (1000*mass_of(lat%param%particle), lat%param%particle)

! Error check that if a superposition attribute was set that "superposition" was set.

do i = 1, n_max
  pele => plat%ele(i)
  if (in_lat%ele(i)%lord_status == super_lord$ .or. pele%superposition_has_been_set) cycle
  if (pele%ref_name /= blank_name$ .or. pele%offset /= 0 .or. &
      pele%ele_pt /= not_set$ .or. pele%ref_pt /= not_set$) then
    call parser_error ('A SUPERPOSITION PARAMETER HAS BEEN SET BUT "SUPERPOSITION" NOT SPECIFIED FOR: ' // in_lat%ele(i)%name)
  endif
enddo

! Now put in the overlay, girder, and group elements

call parser_add_lords (in_lat, n_max, plat, lat)

! fork element to element bookkeeping

call parser_identify_fork_to_element(lat)

! Remove all null_ele elements

do n = 0, ubound(lat%branch, 1)
  branch => lat%branch(n)
  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    if (ele%key == null_ele$) ele%ix_ele = -1 ! mark for deletion
  enddo
enddo
call remove_eles_from_lat (lat, .false.)

if (bp_com%error_flag) then
  call parser_end_stuff (in_lat)
  return
endif

! Bookkeeping...
! Must do this before calling bmad_parser2 since after an expand_lattice command the lattice 
! file may contain references to dependent element parameters that are computed in lattice_bookkeeper.

call cpu_time(bp_com%time1)

call drift_and_pipe_track_methods_adjustment(lat)
call set_flags_for_changed_attribute(lat)
call parser_init_custom_elements (lat)

call s_calc(lat)
call lattice_bookkeeper (lat, err)
if (err) then
  call parser_end_stuff (in_lat, set_error_flag = .true.)
  return
endif

! Consistancy check

call lat_sanity_check (lat, err)
if (err) then
  call parser_end_stuff (in_lat, set_error_flag = .true.)
  return
endif

! Put in field overlaps

do i = 1, n_max
  lord => in_lat%ele(i)
  pele => plat%ele(i)
  if (.not. allocated(pele%field_overlaps)) cycle
  call lat_ele_locator (lord%name, lat, eles, n)
  if (n < 1) cycle

  do j = 1, size(pele%field_overlaps)
    call create_field_overlap (lat, lord%name, pele%field_overlaps(j), err)
    if (err) then
      call parser_error ('CANNOT FIND ELEMENT: ' // pele%field_overlaps(j), &
                         'WHICH HAS FIELD OVERLAP FROM ELEMENT: ' // lord%name)
    endif
  enddo
enddo

! Do we need to call bmad_parser2?

if (bp_com%detected_expand_lattice_cmd) then
  exit_on_error = global_com%exit_on_error
  global_com%exit_on_error = .false.
  bp_com%bmad_parser_calling = .true.
  bp_com%old_lat => in_lat
  call bmad_parser2 ('FROM: BMAD_PARSER', lat, make_mats6 = .false., parse_lat = in_lat)
  bp_com%bmad_parser_calling = .false.
  global_com%exit_on_error = exit_on_error

else
  ! Apply ramper elements?
  call apply_all_rampers(lat, err)
  if (err) call parser_error ('ERROR APPLYING RAMPERS')
  call lattice_bookkeeper (lat)
endif

!

do n = 0, ubound(lat%branch, 1)
  branch => lat%branch(n)
  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    select case (ele%key)
    case (sample$, diffraction_plate$, photon_init$) 
      call photon_target_setup (ele)
    end select
  enddo
enddo

! Make the transfer matrices.
! Note: The bmad_parser err_flag argument does *not* include errors in 
! lat_make_mat6 since if there is a match element, there is an error raised 
! here since the Twiss parameters have not been set. But this is expected. 

call cpu_time(bp_com%time2)

bmad_com%auto_bookkeeper = auto_bookkeeper_saved  ! potentially saves time with lat_make_mat6
if (logic_option (.true., make_mats6)) call lat_make_mat6(lat, ix_branch = -1) 

call g_integrals_calc(lat)

call cpu_time(bp_com%time3)

!

call create_concatenated_wall3d (lat, err)
if (err) then
  call parser_end_stuff (in_lat)
  return
endif

!-------------------------------------------------------------------------
! Print lattice info if debug is on

if (debug_line /= '') call parser_debug_print_info (lat, debug_line, sequence(1:iseq_tot))

! Write digested file

if (.not. bp_com%error_flag .and. .not. bp_com%always_parse) then
  bp_com%write_digested = bp_com%write_digested .and. digested_version <= bmad_inc_version$
  if (bp_com%write_digested) then
    call write_digested_bmad_file (digested_file, lat, bp_com%num_lat_files, bp_com%lat_file_names, bp_com%extra, err)
    if (.not. err) call out_io (s_info$, r_name, 'Created new digested file')
  endif
endif

call lat_sanity_check (lat, err)
if (err) bp_com%error_flag = .true.

call parser_end_stuff (in_lat)

if (bp_com%extra%undeterministic_ran_function_called) then
  call out_io(s_warn$, r_name, &
                'NOTE: THE RANDOM NUMBER FUNCTION WAS USED IN THE LATTICE FILE SO THIS', &
                '      LATTICE WILL DIFFER FROM OTHER LATTICES GENERATED FROM THE SAME FILE.')
endif

call out_io (s_important$, r_name, 'Lattice parse time(min):\f5.2\ ', r_array = [(bp_com%time1 - bp_com%time0)/60.0_rp])


!---------------------------------------------------------------------
contains

subroutine parser_end_stuff (lat0, set_error_flag)

type (lat_struct) lat0
type (ele_struct), pointer :: ele

logical, optional :: set_error_flag
integer i, j, stat_b(24), stat, ierr
character(400) name

!

if (present(parse_lat) .and. .not. bp_com%error_flag) then
  lat0%ramper_slave_bookkeeping = super_ok$    ! Prevents generation of warnings when ramper_slave_setup is called in next line. 
  parse_lat = lat0
endif

! Calculate the creation hash which can be used by programs to verify that the lattice has not been changed since
! the last time the lattice was read in.

lat%creation_hash = djb_hash(int_str(bmad_inc_version$))
do i = 1, bp_com%num_lat_files
  name = bp_com%lat_file_names(i)
  if (name(1:10) == '!DIGESTED:') cycle  ! Ignore digested file
  stat_b = 0
  ierr = stat(name, stat_b)
  lat%creation_hash = djb_hash(int_str(stat_b(2)) // int_str(stat_b(8)), lat%creation_hash)
enddo

! Restore auto_bookkeeper flag

bmad_com%auto_bookkeeper = auto_bookkeeper_saved

! deallocate pointers

if (allocated (bp_com%lat_file_names))   deallocate (bp_com%lat_file_names)

if (logic_option(.false., set_error_flag)) bp_com%error_flag = .true.
if (bp_com%error_flag) then
  if (global_com%exit_on_error) then
    call out_io (s_fatal$, r_name, 'BMAD_PARSER FINISHED. EXITING ON ERRORS')
    stop
  endif
endif

if (present(err_flag)) err_flag = bp_com%error_flag

! Check wiggler for integer number of periods

do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (ele%key /= wiggler$ .and. ele%key /= undulator$) cycle
  if (ele%field_calc /= planar_model$ .and. ele%field_calc /= helical_model$) cycle
  if (ele%slave_status == super_slave$) cycle
  if (abs(ele%value(n_period$) - nint(ele%value(n_period$))) > 0.01) then
    call out_io (s_warn$, r_name, [&
          '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', &
          '!!!!! WARNING! WIGGLER: ' // ele%name, &
          '!!!!! DOES NOT HAVE AN INTEGER NUMBER OF PERIODS!               ', &
          '!!!!! THIS WILL BE PROBLEMATIC IF YOU ARE USING TAYLOR MAPS!    ', &
          '!!!!! SEE THE BMAD MANUAL FOR MORE DETAILS!                     ', &
          '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' ])
  endif
enddo

call deallocate_lat_pointers(lat0)

bp_com%parser_name = ''

end subroutine parser_end_stuff 

!---------------------------------------------------------------------
! contains

subroutine new_element_init (word1, class_word, lat0, plat, err, ele_key)

type (lat_struct), target :: lat0
type (parser_lat_struct), target :: plat

integer, pointer :: n_max
integer, optional :: ele_key
integer ix
logical err, added, is_marker
character(*) word1, class_word

! "end" is a reserved word except if it is being used to instantiate a marker.

n_max => lat0%n_ele_max

select case (word1)
case ('BEGINNING', 'BEAM', 'PARTICLE_START')
  call parser_error ('NAME OF ELEMENT or LINE or LIST CORRESPONDS TO A RESERVED WORD: ' // word1)
  err = .true.
  return
case ('END')
  if (any(class_word == lat0%ele(1:n_max)%name)) then
    call find_index (class_word, lat0%nametable, ix, add_to_list = .false.)
    is_marker = (lat0%ele(ix)%key == marker$)
  else
    is_marker = (index('MARKER', trim(class_word)) /= 0)
  endif

  if (.not. is_marker) then
    call parser_error ('NAME OF ELEMENT or LINE or LIST CORRESPONDS TO A RESERVED WORD: ' // word1)
    err = .true.
    return
  endif
end select

err = .false.

if (n_max >= ubound(lat0%ele, 1)) then
  call allocate_lat_ele_array (lat0)
  call allocate_plat (plat, ubound(lat0%ele, 1))
endif

n_max = n_max + 1
lat0%ele(n_max)%name = word1
call find_index (word1, lat0%nametable, ix, add_to_list = .true., has_been_added = added)
if (.not. added) then
  call parser_error ('DUPLICATE ELEMENT, LINE, OR LIST NAME: ' // word1)
endif

lat0%ele(n_max)%ixx = n_max  ! Pointer to plat%ele() array
if (present(ele_key)) lat0%ele(n_max)%key = ele_key
plat%ele(n_max)%lat_file = bp_com%current_file%full_name
plat%ele(n_max)%ix_line_in_file = bp_com%current_file%i_line

end subroutine new_element_init

!---------------------------------------------------------------------
! contains

subroutine set_this_twiss_struct (t_in, t_out)
type (twiss_struct) t_in, t_out
call set_this_real_val(t_in%beta,       t_out%beta)
call set_this_real_val(t_in%alpha,      t_out%alpha)
call set_this_real_val(t_in%gamma,      t_out%gamma)
call set_this_real_val(t_in%phi,        t_out%phi)
call set_this_real_val(t_in%eta,        t_out%eta)
call set_this_real_val(t_in%etap,       t_out%etap)
call set_this_real_val(t_in%deta_ds,    t_out%deta_ds)
call set_this_real_val(t_in%sigma,      t_out%sigma)
call set_this_real_val(t_in%sigma_p,    t_out%sigma_p)
call set_this_real_val(t_in%emit,       t_out%emit)
call set_this_real_val(t_in%norm_emit,  t_out%norm_emit)
call set_this_real_val(t_in%dbeta_dpz,  t_out%dbeta_dpz)
call set_this_real_val(t_in%dalpha_dpz, t_out%dalpha_dpz)
call set_this_real_val(t_in%deta_dpz,  t_out%deta_dpz)
call set_this_real_val(t_in%detap_dpz,  t_out%detap_dpz)
end subroutine set_this_twiss_struct

!---------------------------------------------------------------------
! contains

subroutine set_this_xy_disp_struct (t_in, t_out)
type (xy_disp_struct) t_in, t_out
call set_this_real_val(t_in%eta,       t_out%eta)
call set_this_real_val(t_in%etap,      t_out%etap)
call set_this_real_val(t_in%deta_ds,   t_out%deta_ds)
call set_this_real_val(t_in%sigma,     t_out%sigma)
end subroutine set_this_xy_disp_struct

!---------------------------------------------------------------------
! contains

subroutine set_this_real_val (val_in, val_out)
real(rp) val_in, val_out
if (val_in /= real_garbage$) val_out = val_in
end subroutine set_this_real_val


end subroutine bmad_parser
