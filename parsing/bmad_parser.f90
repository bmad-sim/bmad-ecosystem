!+
! Subroutine bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
!
! Subroutine to parse a BMAD input file and put the information in lat.
!
! Because of the time it takes to parse a file bmad_parser will save 
! LAT in a "digested" file with the name:
!               'digested_' // lat_file   
! For subsequent calls to the same lat_file, BMAD_PARSER will just read in the
! digested file. bmad_parser will always check to see that the digested file
! is up-to-date and if not the digested file will not be used.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_file   -- Character(*): Name of the input file.
!   make_mats6 -- Logical, optional: Compute the 6x6 transport matrices for the
!                   Elements and do other bookkeeping like calculate the reference energy? 
!                   Default is True.
!   use_line   -- Character(*), optional: If present and not blank, override the use 
!                   statement in the lattice file and use use_line instead.
!
! Output:
!   lat              -- lat_struct: Lat structure. See bmad_struct.f90 for more details.
!     %ele(:)%mat6      -- This is computed assuming an on-axis orbit 
!     %ele(:)%s         -- This is also computed.
!   digested_read_ok -- Logical, optional: Set True if the digested file was
!                        successfully read. False otherwise.
!   err_flag         -- Logical, optional: Set true if there is an error, false otherwise.
!                         Note: err_flag does *not* include errors in lat_make_mat6 since
!                         if there is a match element, there is an error raised since
!                         the Twiss parameters have not been set but this is expected. 
!
! Defaults:
!   lat%param%particle          = positron$
!   lat%param%lattice_type      = circular_lattice$
!   lat%param%aperture_limit_on = .true.
!-

subroutine bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)

use bmad_parser_mod, except_dummy => bmad_parser
use sim_utils
use multipole_mod
use random_mod

implicit none

type (lat_struct), target :: lat, in_lat, old_lat, lat2
type (ele_struct) this_ele
type (ele_struct), save, pointer :: ele, slave, lord
type (ele_struct), save :: marker_ele
type (seq_struct), save, target :: sequence(1000)
type (branch_struct), pointer :: branch0, branch
type (parser_lat_struct), target :: plat
type (parser_ele_struct), pointer :: pele
type (lat_ele_loc_struct) m_slaves(100)

integer, allocatable :: seq_indexx(:), in_indexx(:)
character(40), allocatable ::  in_name(:), seq_name(:)

integer ix_word, i_use, i, j, k, k2, n, ix, ix1, ix2, n_wall, n_track
integer n_ele_use, digested_version, key, loop_counter, n_ic, n_con
integer  iseq_tot, ix_multipass, n_ele_max, n_multi, n0, n_ele
integer, pointer :: n_max

character(*) lat_file
character(*), optional :: use_line

character(1) delim
character(40) word_2, name
character(16), parameter :: r_name = 'bmad_parser'
character(40) this_name, word_1
character(200) full_lat_file_name, digested_file, call_file
character(280) parse_line_save
character(80) debug_line

logical, optional :: make_mats6, digested_read_ok, err_flag
logical delim_found, arg_list_found, xsif_called, print_err, wild_here
logical end_of_file, found, good_attrib, err, finished, exit_on_error
logical detected_expand_lattice_cmd, multipass, wildcards_permitted, matched
logical auto_bookkeeper_saved

! see if digested file is open and current. If so read in and return.
! Note: The name of the digested file depends upon the real precision.

auto_bookkeeper_saved = bmad_com%auto_bookkeeper
bmad_com%auto_bookkeeper = .true.  

if (present(err_flag)) err_flag = .true.
bp_com%error_flag = .false.              ! set to true on an error
bp_com%parser_name = 'bmad_parser'       ! Used for error messages.
bp_com%write_digested = .true.
bp_com%do_superimpose = .true.
bp_com%input_from_file = .true.
debug_line = ''

call form_digested_bmad_file_name (lat_file, digested_file, full_lat_file_name)
call read_digested_bmad_file (digested_file, lat, digested_version, err_flag = err)

! Must make sure that if use_line is present the digested file has used the 
! correct line

if (present(use_line)) then
  if (use_line /= '') then
    call str_upcase (name, use_line)
    if (name /= lat%use_name) err = .true.
  endif
endif

if (.not. err .and. .not. bp_com%always_parse) then
  call set_taylor_order (lat%input_taylor_order, .false.)
  call set_ptc (1.0e9_rp, lat%param%particle)  ! Energy value used does not matter here
  if (lat%input_taylor_order == bmad_com%taylor_order) then
    if (present(digested_read_ok)) digested_read_ok = .true.
    call parser_end_stuff (.false.)
    return
  else
    if (bmad_status%type_out) then
       call out_io (s_info$, r_name, 'Taylor_order has changed.', &
           'Taylor_order in digested file: \i4\ ', &
           'Taylor_order now:              \i4\ ', &
           i_array = [lat%input_taylor_order, bmad_com%taylor_order ])
    endif
    if (lat%input_taylor_order > bmad_com%taylor_order) bp_com%write_digested = .false.
  endif
endif

if (present(digested_read_ok)) digested_read_ok = .false.

! save lattice for possible Taylor map reuse.

old_lat = lat

! here if not OK bmad_status. So we have to do everything from scratch...
! init variables.

call init_lat (lat, 1)
call init_lat (in_lat, 1000)
allocate (in_indexx(0:1000), in_name(0:1000))

nullify (plat%ele)
call allocate_plat (plat, ubound(in_lat%ele, 1))

do i = 0, ubound(in_lat%ele, 1)
  in_lat%ele(i)%ixx = i   ! Pointer to plat%ele() array
enddo

if (bmad_status%type_out) call out_io (s_info$, r_name, 'Parsing lattice file(s). This might take a few minutes...')
call parser_file_stack('init')
call parser_file_stack('push', lat_file, finished, err)  ! open file on stack
if (err) then
  call parser_end_stuff (.false.)
  return
endif

iseq_tot = 0                            ! number of sequences encountered

call ran_seed_get (state = bp_com%ran%initial_state) ! Get initial random state.
if (bp_com%ran%initial_state%ix == -1) then
  bp_com%ran%deterministic = 0
else
  bp_com%ran%deterministic = 1
endif

bp_com%input_line_meaningful = .true.
bp_com%ran%ran_function_was_called = .false.
bp_com%ran%deterministic_ran_function_was_called = .false.
bp_com%e_tot_set = .false.
bp_com%p0c_set   = .false.

call find_indexx2 (in_lat%ele(0)%name, in_name, in_indexx, 0, -1, ix, add_to_list = .true.)

bp_com%beam_ele => in_lat%ele(1)
bp_com%beam_ele%name = 'BEAM'                 ! For mad compatibility.
bp_com%beam_ele%key = def_beam$               ! "definition of beam"
bp_com%beam_ele%value(particle$) = positron$  ! default
call find_indexx2 (in_lat%ele(1)%name, in_name, in_indexx, 0, 0, ix, add_to_list = .true.)

bp_com%param_ele => in_lat%ele(2)
bp_com%param_ele%name = 'PARAMETER'           ! For parameters 
bp_com%param_ele%key = def_parameter$
bp_com%param_ele%value(lattice_type$) = -1
bp_com%param_ele%value(particle$)     = positron$  ! default
call find_indexx2 (in_lat%ele(2)%name, in_name, in_indexx, 0, 1, ix, add_to_list = .true.)

! The parser actually puts beam_start parameters into in_lat%beam_start so
! the beam_start_ele is actually only needed to act as a place holder.

bp_com%beam_start_ele => in_lat%ele(3)
bp_com%beam_start_ele%name = 'BEAM_START'           ! For beam starting parameters 
bp_com%beam_start_ele%key = def_beam_start$
call find_indexx2 (in_lat%ele(3)%name, in_name, in_indexx, 0, 2, ix, add_to_list = .true.)

bp_com%root_branch_ele => in_lat%ele(4)
bp_com%root_branch_ele%name = 'ROOT_BRANCH'           ! For parameters 
bp_com%root_branch_ele%key = branch$
call find_indexx2 (in_lat%ele(4)%name, in_name, in_indexx, 0, 3, ix, add_to_list = .true.)

n_max => in_lat%n_ele_max
n_max = 4                              ! Number of elements encountered

lat%n_control_max = 0
detected_expand_lattice_cmd = .false.

!-----------------------------------------------------------
! main parsing loop

loop_counter = 0  ! Used for debugging
parsing_loop: do 

  loop_counter = loop_counter + 1

  ! get a line from the input file and parse out the first word

  call load_parse_line ('normal', 1, end_of_file)  ! load an input line
  call get_next_word (word_1, ix_word, '[:](,)= ', delim, delim_found, .true.)
  if (end_of_file) then
    word_1 = 'END_FILE'
    ix_word = 8
  else
    wildcards_permitted = (delim == '[')  ! For 'q*[x_offset] = ...' constructs
    call verify_valid_name(word_1, ix_word, wildcards_permitted)
  endif

  ! PARSER_DEBUG

  if (word_1(:ix_word) == 'PARSER_DEBUG') then
    debug_line = bp_com%parse_line
    if (bmad_status%type_out) call out_io (s_info$, r_name, 'Found in file: "PARSER_DEBUG". Debug is now on')
    cycle parsing_loop
  endif

  ! PRINT

  if (word_1(:ix_word) == 'PRINT') then
    call string_trim (bp_com%input_line2, parse_line_save, ix) ! Will strip off initial "print"
    if (bmad_status%type_out) call out_io (s_dwarn$, r_name, &
                                      'Print Message in Lattice File: ' // parse_line_save(ix+1:))
    cycle parsing_loop
  endif

  ! NO_DIGESTED

  if (word_1(:ix_word) == 'NO_DIGESTED') then
    bp_com%write_digested = .false.
    if (bmad_status%type_out) call out_io (s_info$, r_name, &
                            'Found in file: "NO_DIGESTED". No digested file will be created')
    cycle parsing_loop
  endif

  ! NO_SUPERIMPOSE

  if (word_1(:ix_word) == 'NO_SUPERIMPOSE') then
    bp_com%do_superimpose = .false.
    cycle parsing_loop
  endif

  ! DEBUG_MARKER is used to be able to easily set a break within the debugger

  if (word_1(:ix_word) == 'DEBUG_MARKER') then
    word_1 = 'ABC'          ! An executable line to set a break on
    cycle parsing_loop
  endif

  ! USE command...

  if (word_1(:ix_word) == 'USE') then
    if (delim /= ',') call parser_error ('"USE" NOT FOLLOWED BY COMMA')
    call get_next_word(word_2, ix_word, ':(=,)', delim, delim_found, .true.)
    if (ix_word == 0) then 
      call parser_error ('NO BEAM LINE SPECIFIED WITH "USE"', ' ')
      call parser_end_stuff ()
      return
    endif
    call verify_valid_name(word_2, ix_word)
    lat%use_name = word_2
    cycle parsing_loop
  endif

  ! TITLE command

  if (word_1(:ix_word) == 'TITLE') then
    if (delim_found) then
      if (delim /= " " .and. delim /= ",") call parser_error ('BAD DELIMITOR IN "TITLE" COMMAND')
      call bmad_parser_type_get (this_ele, 'DESCRIP', delim, delim_found)
      lat%title = this_ele%descrip
      deallocate (this_ele%descrip)
    else
      read (bp_com%current_file%f_unit, '(a)') lat%title
      bp_com%current_file%i_line = bp_com%current_file%i_line + 1
    endif
    cycle parsing_loop
  endif

  ! CALL command

  if (word_1(:ix_word) == 'CALL') then
    call get_called_file(delim, call_file, xsif_called, err)
    if (err) then
      call parser_end_stuff ()
      return
    endif

    if (xsif_called) then
      call parser_error ('XSIF_PARSER TEMPORARILY DISABLED. PLEASE SEE DCS.')
      if (bmad_status%exit_on_error) call err_exit
      ! call xsif_parser (call_file, lat, make_mats6, digested_read_ok, use_line) 
      detected_expand_lattice_cmd = .true.
      goto 8000  ! Skip the lattice expansion since xsif_parser does this
    endif

    cycle parsing_loop
  endif

  ! BEAM command

  if (word_1(:ix_word) == 'BEAM') then
    if (delim /= ',') call parser_error ('"BEAM" NOT FOLLOWED BY COMMA')
    do 
      if (.not. delim_found) exit
      if (delim /= ',') then
        call parser_error ('EXPECTING: "," BUT GOT: ' // delim, 'FOR "BEAM" COMMAND')
        exit
      endif
      call parser_set_attribute (def$, bp_com%beam_ele, in_lat, delim, delim_found, err, .true.)
    enddo
    cycle parsing_loop
  endif

  ! EXPAND_LATTICE command

  if (word_1(:ix_word) == 'EXPAND_LATTICE') then
    detected_expand_lattice_cmd = .true.
    exit parsing_loop
  endif

  ! RETURN or END_FILE command

  if (word_1(:ix_word) == 'RETURN' .or.  word_1(:ix_word) == 'END_FILE') then
    call parser_file_stack ('pop', ' ', finished, err)
    if (err) then
      call parser_end_stuff ()
      return
    endif
    if (finished) exit parsing_loop ! break loop
    cycle parsing_loop
  endif

  ! variable definition or element redef...

  ! if an element attribute redef

  if (delim == '[') then

    call get_next_word (word_2, ix_word, ']', delim, delim_found, .true.)
    if (.not. delim_found) then
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

    ! Find associated element and evaluate the attribute value.

    wild_here = .false.
    if (index(word_1, '*') /= 0 .or. index(word_1, '%') /= 0) wild_here = .true.
    matched = .false.

    found = .false.
    good_attrib = .false.

    print_err = .true.
    if (word_1 == '*') print_err = .false.

    if (wild_here .or. any(word_1 == key_name)) then
      do i = 0, n_max

        ele => in_lat%ele(i)

        select case (ele%name)
        case ('BEGINNING', 'BEAM', 'PARAMETER', 'BEAM_START')
          matched = .false.
        case default
          matched = match_wild(ele%name, word_1)
        end select

        if (ele%name /= word_1 .and. key_name(ele%key) /= word_1 .and. .not. matched) cycle

        bp_com%parse_line = trim(word_2) // ' = ' // bp_com%parse_line 
        if (found) then   ! if not first time
          bp_com%parse_line = parse_line_save
        else
          parse_line_save = bp_com%parse_line
        endif
        call parser_set_attribute (redef$, ele, in_lat, delim, delim_found, &
                                                  err, print_err, plat%ele(ele%ixx))
        if (.not. err .and. delim_found) call parser_error ('BAD DELIMITER: ' // delim)
        found = .true.
        if (.not. err) good_attrib = .true.

      enddo

    else  ! Not wild

      call find_indexx2 (word_1, in_name, in_indexx, 0, n_max, ix, ix2)
      do i = ix2, n_max
        if (in_name(in_indexx(i)) /= word_1) exit

        ele => in_lat%ele(in_indexx(i))
        bp_com%parse_line = trim(word_2) // ' = ' // bp_com%parse_line 
        if (found) then   ! if not first time
          bp_com%parse_line = parse_line_save
        else
          parse_line_save = bp_com%parse_line
        endif
        call parser_set_attribute (redef$, ele, in_lat, delim, delim_found, &
                                                    err, print_err, plat%ele(ele%ixx))
        if (.not. err .and. delim_found) call parser_error ('BAD DELIMITER: ' // delim)
        found = .true.
        if (.not. err) good_attrib = .true.
      enddo

    endif

    ! If not found then the parse line must be cleared.
    ! If not found then issue a warning except if a general key redef ("quadrupole[...] = ...").

    if (.not. found) then
      bp_com%parse_line = ''
      if (key_name_to_key_index (word_1, .false.) == -1) call parser_error ('ELEMENT NOT FOUND: ' // word_1)
    endif

    if (found .and. .not. print_err .and. .not. good_attrib) then
      call parser_error ('BAD ATTRIBUTE')
    endif

    cycle parsing_loop

  ! else must be a variable

  elseif (delim == '=') then

    call parser_add_variable (word_1, in_lat)
    cycle parsing_loop

  endif

  ! if a "(" delimitor then we are looking at a replacement line.

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
    call parser_error ('1ST DELIMITER IS NOT ":". IT IS: ' // delim,  &
                                                     'FOR: ' // word_1)
    cycle parsing_loop
  endif

  ! only possibilities left are: element, list, or line
  ! to decide which look at 2nd word

  call get_next_word(word_2, ix_word, ':=,', delim, delim_found, .true.)
  if (ix_word == 0) then
    call parser_error ('NO NAME FOUND AFTER: ' // word_1, ' ')
    call parser_end_stuff 
    return
  endif

  if (word_2 == 'LINE[MULTIPASS]') then
    word_2 = 'LINE'
    ix_word = 4
    multipass = .true.
  else
    multipass = .false.
  endif

  call verify_valid_name(word_2, ix_word)

  ! arg lists are only used with lines

  if (word_2(:ix_word) /= 'LINE' .and. arg_list_found) then
    call parser_error ('ARGUMENTS "XYZ(...):" ARE ONLY USED WITH REPLACEMENT LINES.', &
                                                      'FOR: ' // word_1)
    cycle parsing_loop
  endif

  ! if line or list

  if (word_2(:ix_word) == 'LINE' .or. word_2(:ix_word) == 'LIST') then
    iseq_tot = iseq_tot + 1
    if (iseq_tot > size(sequence)-1) then
      call out_io (s_fatal$, r_name, 'ERROR IN BMAD_PARSER: NEED TO INCREASE LINE ARRAY SIZE!')
      if (bmad_status%exit_on_error) call err_exit
    endif

    sequence(iseq_tot)%name = word_1
    sequence(iseq_tot)%multipass = multipass

    if (delim /= '=') call parser_error ('EXPECTING: "=" BUT GOT: ' // delim)
    if (word_2(:ix_word) == 'LINE') then
      sequence(iseq_tot)%type = line$
      if (arg_list_found) sequence(iseq_tot)%type = replacement_line$
    else
      sequence(iseq_tot)%type = list$
    endif
    call seq_expand1 (sequence, iseq_tot, in_lat, .true.)

  ! if not line or list then must be an element

  else

    if (word_1 == 'BEGINNING' .or. word_1 == 'BEAM' .or. word_1 == 'BEAM_START') then
      call parser_error ('ELEMENT NAME CORRESPONDS TO A RESERVED WORD: ' // word_1)
      cycle parsing_loop
    endif

    n_max = n_max + 1
    if (n_max > ubound(in_lat%ele, 1)) then
      call allocate_lat_ele_array (in_lat)
      call re_allocate2 (in_name, 0, ubound(in_lat%ele, 1))
      call re_allocate2 (in_indexx, 0, ubound(in_lat%ele, 1))
      bp_com%beam_ele => in_lat%ele(1)
      bp_com%param_ele => in_lat%ele(2)
      bp_com%beam_start_ele => in_lat%ele(3)
      call allocate_plat (plat, ubound(in_lat%ele, 1))
    endif

    call init_ele (in_lat%ele(n_max), lat = in_lat)
    in_lat%ele(n_max)%name = word_1
    call find_indexx2 (in_lat%ele(n_max)%name, in_name, in_indexx, 0, n_max-1, ix, add_to_list = .true.)
    in_lat%ele(n_max)%ixx = n_max  ! Pointer to plat%ele() array

    plat%ele(n_max)%lat_file = bp_com%current_file%full_name
    plat%ele(n_max)%ix_line_in_file = bp_com%current_file%i_line

    ! Check for valid element key name or if element is part of a element key.
    ! If none of the above then we have an error.

    found = .false.  ! found a match?

    call find_indexx2 (word_2, in_name, in_indexx, 0, n_max, i)
    if (i >= 0 .and. i < n_max) then ! i < n_max avoids "abc: abc" construct.
      in_lat%ele(n_max) = in_lat%ele(i)
      in_lat%ele(n_max)%ixx = n_max  ! Restore correct value
      in_lat%ele(n_max)%name = word_1
      found = .true.
    endif

    if (.not. found) then
      in_lat%ele(n_max)%key = key_name_to_key_index(word_2, .true.)
      if (in_lat%ele(n_max)%key > 0) then
        call set_ele_defaults (in_lat%ele(n_max))
        found = .true.
      endif
    endif

    if (.not. found) then
      call parser_error ('KEY NAME NOT RECOGNIZED OR AMBIGUOUS: ' // word_2,  &
                    'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
      cycle parsing_loop
    endif

    ! Element definition...
    ! Overlay/group/girder case.

    key = in_lat%ele(n_max)%key
    if (key == overlay$ .or. key == group$ .or. key == girder$) then
      if (delim /= '=') then
        call parser_error ('EXPECTING: "=" BUT GOT: ' // delim,  &
                    'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
        cycle parsing_loop        
      endif

      if (key == overlay$) in_lat%ele(n_max)%lord_status = overlay_lord$
      if (key == group$)   in_lat%ele(n_max)%lord_status = group_lord$
      if (key == girder$)  in_lat%ele(n_max)%lord_status = girder_lord$

      call get_overlay_group_names(in_lat%ele(n_max), in_lat, &
                                                  plat%ele(n_max), delim, delim_found)

      if (key /= girder$ .and. .not. delim_found) then
        call parser_error ('NO CONTROL ATTRIBUTE GIVEN AFTER CLOSING "}"',  &
                      'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
        cycle parsing_loop
      endif

    endif

    ! Not overlay/group/girder case.
    ! Loop over all attributes...

    do 
      if (.not. delim_found) exit   ! If no delim then we are finished with this element.
      if (delim /= ',') then
        call parser_error ('EXPECTING: "," BUT GOT: ' // delim,  &
                      'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
        exit
      endif
      call parser_set_attribute (def$, in_lat%ele(n_max), in_lat, delim, delim_found, &
                      err, .true., plat%ele(n_max))
      if (err) cycle parsing_loop
    enddo

  endif

enddo parsing_loop       ! main parsing loop

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! we now have read in everything. 

bp_com%input_line_meaningful = .false.

! sort elements and lists and check for duplicates
! seq_name(:) and in_name(:) arrays speed up the calls to find_indexx since
! the compiler does not have to repack the memory.

allocate (seq_indexx(iseq_tot), seq_name(iseq_tot))
seq_name = sequence(1:iseq_tot)%name
call indexx (seq_name, seq_indexx)

do i = 1, iseq_tot-1
  ix1 = seq_indexx(i)
  ix2 = seq_indexx(i+1)
  if (sequence(ix1)%name == sequence(ix2)%name) call parser_error  &
                    ('DUPLICATE LINE NAME ' // sequence(ix1)%name)
enddo

do i = 1, n_max-1
  ix1 = in_indexx(i)
  ix2 = in_indexx(i+1)
  if (in_lat%ele(ix1)%name == in_lat%ele(ix2)%name) call parser_error &
                  ('DUPLICATE ELEMENT NAME ' // in_lat%ele(ix1)%name)
enddo

i = 1; j = 1
do
  if (i > iseq_tot) exit
  if (j > n_max) exit
  ix1 = seq_indexx(i)
  ix2 = in_indexx(j)
  if (sequence(ix1)%name == in_lat%ele(ix2)%name) call parser_error  &
        ('LINE AND ELEMENT HAVE THE SAME NAME: ' // sequence(ix1)%name)
  if (sequence(ix1)%name < in_lat%ele(ix2)%name) then
    i = i + 1
  else
    j = j + 1
  endif
enddo

! find line corresponding to the "use" statement and expand the used line.

if (present (use_line)) then
  if (use_line /= '') call str_upcase (lat%use_name, use_line)
endif

if (lat%use_name == blank_name$) then
  call parser_error ('NO "USE" STATEMENT FOUND.', 'I DO NOT KNOW WHAT LINE TO USE!')
  call parser_end_stuff ()
  return
endif

allocate (bp_com%used_line(n_max))

call parser_expand_line (0, lat, lat%use_name, sequence, in_name, in_indexx, &
                                      seq_name, seq_indexx, in_lat, n_ele_use)

if (bp_com%error_flag) then
  call parser_end_stuff ()
  return
endif

lat%n_ele_track = n_ele_use
lat%n_ele_max   = n_ele_use

!---------------------------------------------------------------
! If a girder elements refer to a line then must expand that line.

do i = 1, n_max
  lord => in_lat%ele(i)
  if (lord%lord_status /= girder_lord$) cycle
  pele => plat%ele(i)
  j = 0
  do 
    j = j + 1
    if (j > lord%n_slave) exit
    call find_indexx(pele%name(j), seq_name, seq_indexx, size(seq_name), k, k2)
    if (k == 0) cycle
    call parser_expand_line (0, lat2, pele%name(j), sequence, in_name, in_indexx, &
                                      seq_name, seq_indexx, in_lat, n_ele_use)
    ! Put elements from the line expansion into the slave list.
    ! Remember to ignore drifts.
    lord%n_slave = lord%n_slave - 1   ! Remove beam line name
    pele%name(1:lord%n_slave) = [pele%name(1:j-1), pele%name(j+1:lord%n_slave+1)]
    call re_associate (pele%name, lord%n_slave+n_ele_use)
    do k = 1, n_ele_use
      call find_indexx2 (lat2%ele(k)%name, in_name, in_indexx, 0, n_max, ix, ix2)      
      if (ix /= 0) then
        if (in_lat%ele(ix)%key == drift$) cycle
      endif
      lord%n_slave = lord%n_slave + 1
      pele%name(lord%n_slave) = lat2%ele(k)%name
    enddo
  enddo
enddo

!---------------------------------------------------------------
! We have the line to use in constructing the lat.
! now to put the elements in LAT in the correct order.
! superimpose, overlays, and groups are handled later.
! first load beam parameters.

lat%version = bmad_inc_version$
lat%input_file_name         = full_lat_file_name             ! save input file  

lat%ele(0)     = in_lat%ele(0)    ! Beginning element
lat%beam_start = in_lat%beam_start
lat%a          = in_lat%a
lat%b          = in_lat%b
lat%z          = in_lat%z
lat%absolute_time_tracking      = in_lat%absolute_time_tracking
lat%rf_auto_scale_phase         = in_lat%rf_auto_scale_phase
lat%rf_auto_scale_amp           = in_lat%rf_auto_scale_amp
lat%use_ptc_layout              = in_lat%use_ptc_layout

call mat_make_unit (lat%ele(0)%mat6)
call clear_lat_1turn_mats (lat)

if (bp_com%beam_ele%value(n_part$) /= 0 .and. bp_com%param_ele%value(n_part$) /= 0) &
          call parser_error ('BOTH "PARAMETER[N_PART]" AND "BEAM, N_PART" SET.')
lat%param%n_part = max(bp_com%beam_ele%value(n_part$), bp_com%param_ele%value(n_part$))

ix1 = nint(bp_com%param_ele%value(particle$))
ix2 = nint(bp_com%beam_ele%value(particle$))
if (ix1 /= positron$ .and. ix2 /= positron$) &
        call parser_error ('BOTH "PARAMETER[PARTICLE]" AND "BEAM, PARTICLE" SET.')
lat%param%particle = ix1
if (ix2 /=  positron$) lat%param%particle = ix2

! The lattice name from a "parameter[lattice] = ..." line is 
! stored the bp_com%param_ele%descrip string

if (associated(bp_com%param_ele%descrip)) then
  lat%lattice = bp_com%param_ele%descrip
  deallocate (bp_com%param_ele%descrip)
endif

! Work on multipass before overlays, groups, and superimpose. 
! This is important since the elements in the lattice get
! renamed and if not done first would confuse any overlays, girders, etc.
! Multipass elements are paired by multipass index and multipass line name

do i = 1, lat%n_ele_track
  if (bp_com%used_line(i)%ix_multipass == 0) cycle
  n_multi = 0  ! number of elements to slave together
  ix_multipass = bp_com%used_line(i)%ix_multipass
  do j = i, lat%n_ele_track
    if (bp_com%used_line(j)%ix_multipass /= ix_multipass) cycle
    if (bp_com%used_line(j)%multipass_line /= bp_com%used_line(i)%multipass_line) cycle
    n_multi = n_multi + 1
    m_slaves(n_multi) = ele_to_lat_loc (lat%ele(j))
    bp_com%used_line(j)%ix_multipass = 0  ! mark as taken care of
  enddo
  call add_this_multipass (lat, m_slaves(1:n_multi))
enddo

! A patch element with ref_orb = patch_out$ gets a lat%control element to keep
! track of there the corresponding patch with ref_orb = patch_in$ is.

do i = lat%n_ele_track+1, lat%n_ele_max
  ele => lat%ele(i)
  if (ele%ref_orbit /= patch_out$) cycle
  if (ele%component_name == '') then
    call parser_error ('NO REF_PATCH ELEMENT GIVEN FOR PATCH: ' // ele%name)
    cycle
  endif
  do j = lat%n_ele_track+1, lat%n_ele_max
    if (lat%ele(j)%name == ele%component_name) then
      call reallocate_control(lat, lat%n_control_max+1)

      n_ic = lat%n_ic_max + 1
      lat%n_ic_max = n_ic

      n_con = lat%n_control_max + 1
      lat%n_control_max = n_con

      ele%slave_status = patch_in_slave$
      ele%n_lord = 1
      ele%ic1_lord = n_ic
      ele%ic2_lord = n_ic

      lat%ic(n_ic) = n_con  
      lat%control(n_con)%ix_lord = j
      lat%control(n_con)%ix_slave = i
      lat%control(n_con)%ix_branch = 0
      exit
    endif
  enddo
  if (j == lat%n_ele_max + 1) then
    call parser_error ('CANNOT FIND REF_PATCH FOR PATCH: ' // ele%name, &
                  'CANNOT FIND: ' // ele%component_name)
  endif
enddo

! Set lattice_type.

ix = nint(bp_com%param_ele%value(lattice_type$))
if (ix > 0) then  ! lattice_type has been set.
  lat%param%lattice_type = ix
else              ! else use default
  lat%param%lattice_type = circular_lattice$      ! default 
  if (any(in_lat%ele(:)%key == lcavity$)) then    !   except...
    if (bmad_status%type_out) call out_io (s_warn$, r_name, 'NOTE: THIS LATTICE HAS A LCAVITY.', &
                                  'SETTING THE LATTICE_TYPE TO LINEAR_LATTICE.')
    lat%param%lattice_type = linear_lattice$
  endif
endif

! Set taylor_order

lat%input_taylor_order = nint(bp_com%param_ele%value(taylor_order$))

if (lat%input_taylor_order /= 0) call set_taylor_order (lat%input_taylor_order, .false.)

!-------------------------------------------------------------------------
! energy bookkeeping.

ele => lat%ele(0)

err = .true.
if (bp_com%e_tot_set .and. ele%value(e_tot$) < mass_of(lat%param%particle)) then
  if (bmad_status%type_out) call out_io (s_error$, r_name, 'REFERENCE ENERGY IS SET BELOW MC^2! WILL USE 1000 * MC^2!')
elseif (bp_com%p0c_set .and. ele%value(p0c$) < 0) then
  if (bmad_status%type_out) call out_io (s_error$, r_name, 'REFERENCE MOMENTUM IS NEGATIVE! WILL USE 1000 * MC^2!')
elseif (.not. (bp_com%p0c_set .or. bp_com%e_tot_set)) then
  if (bmad_status%type_out) call out_io (s_warn$, r_name, 'REFERENCE ENERGY IS NOT SET IN LATTICE FILE! WILL USE 1000 * MC^2!')
else
  err = .false.
endif

if (err) then
  ele%value(e_tot$) = 1000 * mass_of(lat%param%particle)
  call convert_total_energy_to (ele%value(e_tot$), lat%param%particle, pc = ele%value(p0c$))
  bp_com%write_digested = .false.
elseif (bp_com%e_tot_set) then
  call convert_total_energy_to (ele%value(e_tot$), lat%param%particle, pc = ele%value(p0c$))
elseif (bp_com%p0c_set) then
  call convert_pc_to (ele%value(p0c$), lat%param%particle, e_tot = ele%value(e_tot$))
endif

ele%value(e_tot_start$) = ele%value(e_tot$)
ele%value(p0c_start$) = ele%value(p0c$)

! Use arbitrary energy above the rest mass energy since when tracking individual elements the
! true reference energy is used.

call set_ptc (1000*mass_of(lat%param%particle), lat%param%particle)

! Add branch lines.
! Branch lines may contain branch elements so this is an iterative process

do i = 1, lat%n_ele_max
  if (lat%ele(i)%key /= photon_branch$ .and. lat%ele(i)%key /= branch$) cycle
  if (lat%ele(i)%slave_status == multipass_slave$) cycle
  call parser_add_branch (lat%ele(i), lat, sequence, in_name, in_indexx, seq_name, seq_indexx, in_lat)
enddo

! Error check that if a superposition attribute was set that "superposition" was set.

do i = 1, n_max
  if (in_lat%ele(i)%lord_status == super_lord$) cycle
  pele => plat%ele(i)
  if (pele%ref_name /= blank_name$ .or. pele%s /= 0 .or. &
      pele%ele_pt /= not_set$ .or. pele%ref_pt /= not_set$) then
    call parser_error ('SUPERPOSITION ATTRIBUTE SET BUT "SUPERPOSITION" NOT SPECIFIED FOR: ' // in_lat%ele(i)%name)
  endif
enddo

! Add end marker. Make sure name is not duplicated. Also do not add marker if one already exists.

if (nint(bp_com%param_ele%value(no_end_marker$)) == 0) then
  call init_ele (marker_ele)
  marker_ele%key = marker$

  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    n_track = branch%n_ele_track

    marker_ele%name = 'END'

    do j = 1, n_track
      if (branch%ele(j)%name == 'END') then
        marker_ele%name = 'END_MARKER'
        exit
      endif
    enddo

    if (branch%ele(n_track)%name /= 'END' .or. branch%ele(n_track)%key /= marker$) then
      if (marker_ele%name /= 'END') then
        call parser_error ('ELEMENT NAMED "END" FOUND IN LATTICE.', & 
                           'ENDING MARKER WILL BE NAMED: ' // trim(marker_ele%name), warn_only = .true.)
      endif
      call insert_element (lat, marker_ele, n_track+1, i)
    endif
  enddo
endif

! Go through the IN_LAT elements and put in the superpositions.
! If the superposition is a branch elemen, need to add the branch line.

call s_calc (lat)              ! calc longitudinal distances
call control_bookkeeper (lat)

do i = 1, n_max

  if (in_lat%ele(i)%lord_status /= super_lord$) cycle
  call add_all_superimpose (lat, in_lat%ele(i), plat%ele(i), in_lat)

  if (in_lat%ele(i)%key /= photon_branch$ .and. in_lat%ele(i)%key /= branch$) cycle
  do j = 0, ubound(lat%branch, 1)
    branch => lat%branch(j)
    n_ele = branch%n_ele_max
    do k = 1, n_ele
      if (branch%ele(k)%name /= in_lat%ele(i)%name) cycle
      call parser_add_branch (branch%ele(k), lat, sequence, in_name, in_indexx, &
                                                              seq_name, seq_indexx, in_lat)
      branch => lat%branch(j)
    enddo
  enddo

  call s_calc (lat)  ! calc longitudinal distances of new branche elements

enddo

! Now put in the overlay_lord, girder, and group elements

call parser_add_lord (in_lat, n_max, plat, lat)

! Consistancy check

call check_lat_controls (lat, err)
if (err) then
  bp_com%error_flag = .true.
  call parser_end_stuff
  return
endif

! Reuse the old taylor series if they exist
! and the old taylor series has the same attributes.

if (logic_option (.true., make_mats6)) then
  call lattice_bookkeeper (lat, err)
  if (err) then
    bp_com%error_flag = .true.
    call parser_end_stuff
    return
  endif
endif

lat%input_taylor_order = bmad_com%taylor_order

! Do we need to expand the lattice and call bmad_parser2?

8000 continue

if (detected_expand_lattice_cmd) then
  exit_on_error = bmad_status%exit_on_error
  bmad_status%exit_on_error = .false.
  bp_com%bmad_parser_calling = .true.
  bp_com%old_lat => in_lat
  call bmad_parser2 ('FROM: BMAD_PARSER', lat, make_mats6 = .false.)
  bp_com%bmad_parser_calling = .false.
  bmad_status%exit_on_error = exit_on_error
endif

! Remove all null_ele elements and init custom stuff.

do n = 0, ubound(lat%branch, 1)
  branch => lat%branch(n)
  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    if (ele%key == null_ele$) ele%key = -1 ! mark for deletion
    if (ele%key == custom$ .or. ele%tracking_method == custom$ .or. &
        ele%mat6_calc_method == custom$ .or. ele%field_calc == custom$ .or. &
        ele%aperture_type == custom$) then
      call init_custom (ele, err)
      if (err) bp_com%error_flag = .true.
    endif
  enddo
enddo
call remove_eles_from_lat (lat, .false.)  

! Make the transfer matrices.
! Note: The bmad_parser err_flag argument does *not* include errors in 
! lat_make_mat6 since if there is a match element, there is an error raised 
! here since the Twiss parameters have not been set. But this is expected. 

call reuse_taylor_elements (old_lat, lat)
if (logic_option (.true., make_mats6)) call lat_make_mat6(lat, -1) 

! Aggragate vacuum chamber wall info for a branch to branch%wall3d structure

! NOTE: This code needs to be modified to take care of continuous aperture walls and 
! wall priorities. OR: Is the branch%wall3d component even needed?

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  cycle

  n_wall = 0
  do j = 0, branch%n_ele_track
    ele => branch%ele(j)
    if (ele%key == capillary$) cycle
    if (.not. associated(ele%wall3d)) cycle
    call wall3d_initializer (ele%wall3d, err_flag)
    n_wall = n_wall + size(ele%wall3d%section)
  enddo
  if (n_wall == 0) cycle

  allocate (branch%wall3d)
  allocate (branch%wall3d%section(n_wall))
  n_wall = 0
  do j = 0, branch%n_ele_track
    ele => branch%ele(j)
    if (ele%key == capillary$) cycle
    if (.not. associated(ele%wall3d)) cycle
    n = size(ele%wall3d%section)
    branch%wall3d%section(n_wall+1:n_wall+n) = ele%wall3d%section
    branch%wall3d%section(n_wall+1:n_wall+n)%s = branch%wall3d%section(n_wall+1:n_wall+n)%s + &
                                                                              ele%s - ele%value(l$)
    n_wall = n_wall + n
  enddo

enddo

! Correct beam_start info

call init_coord (lat%beam_start, lat%beam_start, lat%ele(0), shift_vec6 = .false.)

!-------------------------------------------------------------------------
! write out if debug is on

if (debug_line /= '') call parser_debug_print_info (lat, debug_line)

! write to digested file

if (.not. bp_com%error_flag) then            
  bp_com%write_digested = bp_com%write_digested .and. digested_version <= bmad_inc_version$
  if (bp_com%write_digested) then
    call write_digested_bmad_file (digested_file, lat, bp_com%num_lat_files, &
                                                            bp_com%lat_file_names, bp_com%ran, err)
    if (.not. err .and. bmad_status%type_out) call out_io (s_info$, r_name, 'Created new digested file')
  endif
endif

call parser_end_stuff ()

if (bp_com%ran%ran_function_was_called) then
  if (bmad_status%type_out) call out_io(s_warn$, r_name, &
                'NOTE: THE RANDOM NUMBER FUNCTION WAS USED IN THE LATTICE FILE SO THIS', &
                '      LATTICE WILL DIFFER FROM OTHER LATTICES GENERATED FROM THE SAME FILE.')
endif

!---------------------------------------------------------------------
contains

subroutine parser_end_stuff (do_dealloc)

logical, optional :: do_dealloc
integer i, j

! Restore auto_bookkeeper flag

bmad_com%auto_bookkeeper = auto_bookkeeper_saved

! deallocate pointers

if (logic_option (.true., do_dealloc)) then

  do i = lbound(plat%ele, 1) , ubound(plat%ele, 1)
    if (associated (plat%ele(i)%name)) then
      deallocate(plat%ele(i)%name)
      deallocate(plat%ele(i)%attrib_name)
      deallocate(plat%ele(i)%coef)
    endif
  enddo

  do i = 1, size(sequence(:))
    if (associated (sequence(i)%dummy_arg)) &
              deallocate(sequence(i)%dummy_arg, sequence(i)%corresponding_actual_arg)
    if (associated (sequence(i)%ele)) then
      do j = 1, size(sequence(i)%ele)
        if (associated (sequence(i)%ele(j)%actual_arg)) &
                              deallocate(sequence(i)%ele(j)%actual_arg)
      enddo
      deallocate(sequence(i)%ele)
    endif
  enddo

  if (associated (plat%ele))             deallocate (plat%ele)
  if (allocated (seq_indexx))            deallocate (seq_indexx, seq_name)
  if (allocated (in_indexx))             deallocate (in_indexx, in_name)
  if (allocated (bp_com%used_line))      deallocate (bp_com%used_line)
  if (allocated (bp_com%lat_file_names)) deallocate (bp_com%lat_file_names)

endif

if (bp_com%error_flag) then
  if (bmad_status%exit_on_error) then
    call out_io (s_fatal$, r_name, 'BMAD_PARSER FINISHED. EXITING ON ERRORS')
    stop
  endif
endif

if (present(err_flag)) err_flag = bp_com%error_flag

! Check wiggler for integer number of periods

do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (ele%key /= wiggler$) cycle
  if (ele%sub_key /= periodic_type$) cycle
  if (ele%slave_status == super_slave$) cycle
  if (abs(modulo2(ele%value(n_pole$) / 2, 0.5_rp)) > 0.01) then
    call out_io (s_warn$, r_name, [&
          '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', &
          '!!!!! WARNING! WIGGLER: ' // ele%name, &
          '!!!!! DOES NOT HAVE AN EVEN NUMBER OF POLES!                    ', &
          '!!!!! THIS WILL BE PROBLEMATIC IF YOU ARE USING TAYLOR MAPS!    ', &
          '!!!!! SEE THE BMAD MANUAL FOR MORE DETAILS!                     ', &
          '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' ])
  endif
enddo

call deallocate_lat_pointers(old_lat)
call deallocate_lat_pointers(in_lat)
call deallocate_lat_pointers(lat2)

end subroutine

end subroutine
