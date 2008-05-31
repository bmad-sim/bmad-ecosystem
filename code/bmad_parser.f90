!+
! Subroutine bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line)
!
! Subroutine to parse a BMAD input file and put the information in lat.
!
! Because of the time it takes to parse a file BMAD_PARSER will save 
! LAT in a "digested" file with the name:
!               'digested_' // lat_file   ! for single precision BMAD version
!               'digested8_' // lat_file  ! for double precision BMAD version
! For subsequent calls to the same lat_file, BMAD_PARSER will just read in the
! digested file. BMAD_PARSER will always check to see that the digested file
! is up-to-date and if not the digested file will not be used.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_file   -- Character(*): Name of the input file.
!   make_mats6 -- Logical, optional: Compute the 6x6 transport matrices for the
!                   Elements? Default is True.
!   use_line   -- Character(*), optional: If present then override the use 
!                   statement in the lattice file and use use_line instead.
!
! Output:
!   lat -- lat_struct: Lat structure. See bmad_struct.f90 for more details.
!     %ele(:)%mat6  -- This is computed assuming an on-axis orbit 
!     %ele(:)%s     -- This is also computed.
!   digested_read_ok -- Logical, optional: Set True if the digested file was
!                        successfully read. False otherwise.
!   bmad_status      -- Bmad status common block.
!     %ok              -- Set True if parsing is successful. False otherwise.
!         
! Defaults:
!   lat%param%particle          = positron$
!   lat%param%lattice_type      = circular_lattice$
!   lat%param%aperture_limit_on = .true.
!-

#include "CESR_platform.inc"

subroutine bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line)

  use bmad_parser_mod, except_dummy => bmad_parser
  use cesr_utils
  use multipole_mod
  use random_mod

  implicit none

  type (lat_struct), target :: lat, in_lat
  type (ele_struct) this_ele
  type (seq_struct), save, target :: sequence(1000)
  type (ele_struct), pointer :: beam_ele, param_ele, beam_start_ele
  type (branch_struct), pointer :: branch0, branch
  type (parser_lat_struct) plat
  type (ele_struct), save, pointer :: ele
  type (ele_struct), allocatable, save :: old_ele(:) 
  type (used_seq_struct), allocatable ::  used_line(:)

  integer, allocatable :: seq_indexx(:), in_indexx(:)
  character(40), allocatable ::  in_name(:), seq_name(:)

  integer ix_word, i_use, i, j, k, n, ix, ix1, ix2, ixm(100)
  integer n_ele_use, digested_version, key, loop_counter
  integer  iseq_tot, ix_multipass, n_ele_max, n_multi, n0
  integer, pointer :: n_max

  character(*) lat_file
  character(*), optional :: use_line

  character(1) delim
  character(40) word_2, name
  character(16) :: r_name = 'bmad_parser'
  character(40) this_name, word_1
  character(200) full_lat_file_name, digested_file, call_file
  character(280) parse_line_save
  character(80) debug_line

  real(rp) energy_beam, energy_param, energy_0

  logical, optional :: make_mats6, digested_read_ok
  logical delim_found, arg_list_found, doit, xsif_called
  logical file_end, found, err_flag, finished, exit_on_error
  logical detected_expand_lattice_cmd, multipass, write_digested

! see if digested file is open and current. If so read in and return.
! Note: The name of the digested file depends upon the real precision.

  bp_com%parser_name = 'BMAD_PARSER'  ! Used for error messages.
  write_digested = .true.
  debug_line = ''

  call form_digested_bmad_file_name (lat_file, digested_file, full_lat_file_name)
  call read_digested_bmad_file (digested_file, lat, digested_version)

  ! Must make sure that if use_line is present the digested file has used the 
  ! correct line

  if (present(use_line)) then
    call str_upcase (name, use_line)
    if (name /= lat%name) bmad_status%ok = .false.
  endif

  if (bmad_status%ok) then
    call set_taylor_order (lat%input_taylor_order, .false.)
    call set_ptc (lat%ele(0)%value(e_tot$), lat%param%particle)
    if (lat%input_taylor_order == bmad_com%taylor_order) then
      if (present(digested_read_ok)) digested_read_ok = .true.
      return
    else
      if (bmad_status%type_out) then
         call out_io (s_info$, r_name, 'Taylor_order has changed.', &
             'Taylor_order in digested file: \i4\ ', &
             'Taylor_order now:              \i4\ ', &
             i_array = (/ lat%input_taylor_order, bmad_com%taylor_order /) )
      endif
      if (lat%input_taylor_order > bmad_com%taylor_order) &
                                           write_digested = .false.
    endif
  endif

  if (present(digested_read_ok)) digested_read_ok = .false.

! save all elements that have a taylor series

  call save_taylor_elements (lat, old_ele)

! here if not OK bmad_status. So we have to do everything from scratch...
! init variables.

  nullify (plat%ele)
  call init_lat (in_lat, 1000)
  call init_lat (lat, 1)
  call allocate_plat (in_lat, plat)

  bmad_status%ok = .true.
  if (bmad_status%type_out) &
       call out_io (s_info$, r_name, 'Creating new digested file...')

  bp_com%error_flag = .false.                 ! set to true on an error
  call file_stack('push', lat_file, finished)  ! open file on stack
  if (.not. bmad_status%ok) return
  iseq_tot = 0                            ! number of sequences encountered

  bp_com%input_line_meaningful = .true.
  bp_com%ran_function_was_called = .false.

  call init_ele (in_lat%ele(0))
  in_lat%ele(0)%name = 'BEGINNING'     ! Beginning element
  in_lat%ele(0)%key = init_ele$

  beam_ele => in_lat%ele(1)
  call init_ele (beam_ele)
  beam_ele%name = 'BEAM'                 ! fake beam element
  beam_ele%key = def_beam$               ! "definition of beam"
  beam_ele%value(particle$) = positron$  ! default

  param_ele => in_lat%ele(2)
  call init_ele (param_ele)
  param_ele%name = 'PARAMETER'           ! For parameters 
  param_ele%key = def_parameter$
  param_ele%value(lattice_type$) = -1
  param_ele%value(particle$)     = positron$  ! default

  beam_start_ele => in_lat%ele(3)
  call init_ele (beam_start_ele)
  beam_start_ele%name = 'BEAM_START'           ! For parameters 
  beam_start_ele%key = def_beam_start$

  n_max => in_lat%n_ele_max
  n_max = 3                              ! Number of elements encountered

  lat%n_control_max = 0
  detected_expand_lattice_cmd = .false.

!-----------------------------------------------------------
! main parsing loop

  loop_counter = 0  ! Used for debugging
  parsing_loop: do 

    loop_counter = loop_counter + 1

! get a line from the input file and parse out the first word

    call load_parse_line ('normal', 1, file_end)  ! load an input line
    call get_next_word (word_1, ix_word, '[:](,)= ', delim, delim_found, .true.)
    if (file_end) then
      word_1 = 'END_FILE'
      ix_word = 8
    else
      call verify_valid_name(word_1, ix_word)
    endif

! PARSER_DEBUG

    if (word_1(:ix_word) == 'PARSER_DEBUG') then
      debug_line = bp_com%parse_line
      call out_io (s_info$, r_name, 'FOUND IN FILE: "PARSER_DEBUG". DEBUG IS NOW ON')
      cycle parsing_loop
    endif

! NO_DIGESTED

    if (word_1(:ix_word) == 'NO_DIGESTED') then
      write_digested = .false.
      call out_io (s_info$, r_name, 'FOUND IN FILE: "NO_DIGESTED". NO DIGESTED FILE WILL BE CREATED')
      cycle parsing_loop
    endif

! DEBUG_MARKER is used to be able to easily set a break within the debugger

    if (word_1(:ix_word) == 'DEBUG_MARKER') then
      word_1 = 'ABC'          ! An executable line to set a break on
      cycle parsing_loop
    endif

! USE command...

    if (word_1(:ix_word) == 'USE') then
      if (delim /= ',') call warning ('"USE" NOT FOLLOWED BY COMMA')
      call get_next_word(word_2, ix_word, ':(=,)', delim, delim_found, .true.)
      if (ix_word == 0) call error_exit  &
                                ('NO BEAM LINE SPECIFIED WITH "USE"', ' ')
      call verify_valid_name(word_2, ix_word)
      lat%name = word_2
      cycle parsing_loop
    endif

! TITLE command

    if (word_1(:ix_word) == 'TITLE') then
      if (delim_found) then
        if (delim /= " " .and. delim /= ",") call warning &
                            ('BAD DELIMITOR IN "TITLE" COMMAND')
        call type_get (this_ele, descrip$, delim, delim_found)
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
      call get_called_file(delim, call_file, xsif_called)
      if (.not. bmad_status%ok) return

      if (xsif_called) then
        call warning ('XSIF_PARSER TEMPORARILY DISABLED. PLEASE SEE DCS.')
        call err_exit
        ! call xsif_parser (call_file, lat, make_mats6, digested_read_ok, use_line) 
        detected_expand_lattice_cmd = .true.
        goto 8000  ! Skip the lattice expansion since xsif_parser does this
      endif

      cycle parsing_loop
    endif

! BEAM command

    if (word_1(:ix_word) == 'BEAM') then
      if (delim /= ',') call warning ('"BEAM" NOT FOLLOWED BY COMMA')
      do 
        if (.not. delim_found) exit
        if (delim /= ',') then
          call warning ('EXPECTING: "," BUT GOT: ' // delim, 'FOR "BEAM" COMMAND')
          exit
        endif
        call get_attribute (def$, beam_ele, in_lat, plat, delim, delim_found, err_flag)
      enddo
      cycle parsing_loop
    endif
                   
! LATTICE command

    if (word_1(:ix_word) == 'LATTICE') then
      if (delim /= ':' .and. delim /= '=') then
        call warning ('"LATTICE" NOT FOLLOWED BY ":"')
      else
        if (delim == ':' .and. bp_com%parse_line(1:1) == '=') &
                      bp_com%parse_line = bp_com%parse_line(2:)  ! trim off '='
        call get_next_word (lat%lattice, ix_word, ',', &
                                                   delim, delim_found, .true.)
        print *, 'BMAD_PARSER NOTE:'
        print *, '    DEPRECATED USE OF SYNTAX: "LATTICE = ...".'
        print *, '    USE "PARAMETER[LATTICE] = ..." SYNTAX INSTEAD.'
      endif
      cycle parsing_loop
    endif

! EXPAND_LATTICE command

    if (word_1(:ix_word) == 'EXPAND_LATTICE') then
      detected_expand_lattice_cmd = .true.
      exit parsing_loop
    endif

! RETURN or END_FILE command

    if (word_1(:ix_word) == 'RETURN' .or.  &
                                    word_1(:ix_word) == 'END_FILE') then
      call file_stack ('pop', ' ', finished)
      if (.not. bmad_status%ok) return
      if (finished) exit parsing_loop ! break loop
      cycle parsing_loop
    endif

! variable definition or element redef...

! if an element attribute redef

    found = .false.
    if (delim == '[') then

      call get_next_word (word_2, ix_word, ']', delim, delim_found, .true.)
      if (.not. delim_found) then
        call warning ('OPENING "[" FOUND WITHOUT MATCHING "]"')
        cycle parsing_loop
      endif

      call get_next_word (this_name, ix_word, ':=', delim, delim_found, .true.)
      if (.not. delim_found .or. ix_word /= 0) then
        call warning ('MALFORMED ELEMENT ATTRIBUTE REDEFINITION')
        cycle parsing_loop
      endif

      ! If delim is ':' then this is an error since get_next_word treats
      ! a ':=' construction as a '=' 

      if (delim == ':') then
        call warning ('MALFORMED ELEMENT ATTRIBUTE REDEF')
        cycle parsing_loop
      endif

      ! find associated element and evaluate the attribute value

      do i = 0, n_max

        ele => in_lat%ele(i)

        if (ele%name == word_1 .or. key_name(ele%key) == word_1) then
          bp_com%parse_line = trim(word_2) // ' = ' // bp_com%parse_line 
          if (found) then   ! if not first time
            bp_com%parse_line = parse_line_save
          else
            parse_line_save = bp_com%parse_line
          endif
          call get_attribute (redef$, ele, in_lat, plat, &
                                             delim, delim_found, err_flag)
          if (delim_found) call warning ('BAD DELIMITER: ' // delim)
          found = .true.
        endif

      enddo

      if (.not. found) call warning ('ELEMENT NOT FOUND: ' // word_1)
      cycle parsing_loop

! else must be a variable

    elseif (delim == '=') then

      call parser_add_variable (word_1, in_lat)
      cycle parsing_loop

    endif

! if a "(" delimitor then we are looking at a replacement line.

    if (delim == '(') then
      call get_sequence_args (word_1, sequence(iseq_tot+1)%dummy_arg, &
                                                       delim, err_flag)
      ix = size(sequence(iseq_tot+1)%dummy_arg)
      allocate (sequence(iseq_tot+1)%corresponding_actual_arg(ix))
      if (err_flag) cycle parsing_loop
      arg_list_found = .true.
      call get_next_word (word_2, ix_word, '(): =,', &
                                                 delim, delim_found, .true.)
      if (word_2 /= ' ') call warning &
                  ('":" NOT FOUND AFTER REPLACEMENT LINE ARGUMENT LIST. ' // &
                  'FOUND: ' // word_2, 'FOR LINE: ' // word_1)
    else
      arg_list_found = .false.
    endif

! must have a ":" delimiter now

    if (delim /= ':') then
      call warning ('1ST DELIMITER IS NOT ":". IT IS: ' // delim,  &
                                                       'FOR: ' // word_1)
      cycle parsing_loop
    endif

! only possibilities left are: element, list, or line
! to decide which look at 2nd word

    call get_next_word(word_2, ix_word, ':=,', delim, delim_found, .true.)
    if (ix_word == 0) call error_exit ('NO NAME FOUND AFTER: ' // word_1, ' ')

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
      call warning ('ARGUMENTS "XYZ(...):" ARE ONLY USED WITH REPLACEMENT LINES.', &
                                                        'FOR: ' // word_1)
      cycle parsing_loop
    endif

! if line or list

    if (word_2(:ix_word) == 'LINE' .or. word_2(:ix_word) == 'LIST') then
      iseq_tot = iseq_tot + 1
      if (iseq_tot > size(sequence)-1) then
        print *, 'ERROR IN BMAD_PARSER: NEED TO INCREASE LINE ARRAY!'
        call err_exit
      endif

      sequence(iseq_tot)%name = word_1
      sequence(iseq_tot)%multipass = multipass

      if (delim /= '=') call warning ('EXPECTING: "=" BUT GOT: ' // delim)
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
        call warning ('ELEMENT NAME CORRESPONDS TO A RESERVED WORD: ' // word_1)
        cycle parsing_loop
      endif

      n_max = n_max + 1
      if (n_max > ubound(in_lat%ele, 1)) then
        call allocate_lat_ele_array (in_lat)
        beam_ele => in_lat%ele(1)
        param_ele => in_lat%ele(2)
        beam_start_ele => in_lat%ele(3)
        call allocate_plat (in_lat, plat)
      endif

      call init_ele (in_lat%ele(n_max))
      in_lat%ele(n_max)%name = word_1

      plat%ele(n_max)%lat_file = bp_com%current_file%full_name
      plat%ele(n_max)%ix_line_in_file = bp_com%current_file%i_line

! Check for valid element key name or if element is part of a element key.
! If none of the above then we have an error.

      found = .false.  ! found a match?

      in_lat%ele(n_max)%key = key_name_to_key_index(word_2, .true.)
      if (in_lat%ele(n_max)%key > 0) then
        call preparse_element_init (in_lat%ele(n_max))
        found = .true.
      endif

      if (.not. found) then
        do i = 1, n_max-1
          if (word_2 == in_lat%ele(i)%name) then
            in_lat%ele(n_max) = in_lat%ele(i)
            in_lat%ele(n_max)%name = word_1
            found = .true.
            exit
          endif
        enddo
      endif

      if (.not. found) then
        call warning ('KEY NAME NOT RECOGNIZED: ' // word_2,  &
                       'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
        cycle parsing_loop
      endif

! Element definition...
! First: set defaults.

      call parser_set_ele_defaults (in_lat%ele(n_max))

      key = in_lat%ele(n_max)%key
      if (key == overlay$ .or. key == group$ .or. key == girder$) then
        if (delim /= '=') then
          call warning ('EXPECTING: "=" BUT GOT: ' // delim,  &
                      'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
          cycle parsing_loop        
        endif

        if (key == overlay$) in_lat%ele(n_max)%control_type = overlay_lord$
        if (key == group$)   in_lat%ele(n_max)%control_type = group_lord$
        if (key == girder$)  in_lat%ele(n_max)%control_type = girder_lord$

        call get_overlay_group_names(in_lat%ele(n_max), in_lat, &
                                                    plat, delim, delim_found)

        if (key /= girder$ .and. .not. delim_found) then
          call warning ('NO CONTROL ATTRIBUTE GIVEN AFTER CLOSING "}"',  &
                        'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
          cycle parsing_loop
        endif

      endif

! Second: We need to get the attribute values for the element.

      do 
        if (.not. delim_found) exit          ! if nothing more
        if (delim /= ',') then
          call warning ('EXPECTING: "," BUT GOT: ' // delim,  &
                        'FOR ELEMENT: ' // in_lat%ele(n_max)%name)
          exit
        endif
        call get_attribute (def$, in_lat%ele(n_max), &
                                  in_lat, plat, delim, delim_found, err_flag)
        if (err_flag) cycle parsing_loop
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

  allocate (in_indexx(n_max), in_name(n_max))
  in_name = in_lat%ele(1:n_max)%name
  call indexx (in_name, in_indexx)

  do i = 1, iseq_tot-1
    ix1 = seq_indexx(i)
    ix2 = seq_indexx(i+1)
    if (sequence(ix1)%name == sequence(ix2)%name) call warning  &
                      ('DUPLICATE LINE NAME ' // sequence(ix1)%name)
  enddo

  do i = 1, n_max-1
    ix1 = in_indexx(i)
    ix2 = in_indexx(i+1)
    if (in_lat%ele(ix1)%name == in_lat%ele(ix2)%name) call warning &
                    ('DUPLICATE ELEMENT NAME ' // in_lat%ele(ix1)%name)
  enddo

  i = 1; j = 1
  do
    if (i > iseq_tot) exit
    if (j > n_max) exit
    ix1 = seq_indexx(i)
    ix2 = in_indexx(j)
    if (sequence(ix1)%name == in_lat%ele(ix2)%name) call warning  &
          ('LINE AND ELEMENT HAVE THE SAME NAME: ' // sequence(ix1)%name)
    if (sequence(ix1)%name < in_lat%ele(ix2)%name) then
      i = i + 1
    else
      j = j + 1
    endif
  enddo

! find line corresponding to the "use" statement and expand the used line.

  if (present (use_line)) call str_upcase (lat%name, use_line)
  if (lat%name == blank) call error_exit &
            ('NO "USE" STATEMENT FOUND.', 'I DO NOT KNOW WHAT LINE TO USE!')

  allocate (used_line(n_max))

  call parser_expand_line (0, lat%name, sequence, in_name, in_indexx, &
                  seq_name, seq_indexx, in_lat%ele, lat%ele, used_line, n_ele_use)
  if (bp_com%error_flag) stop

  lat%n_ele_track        = n_ele_use
  lat%n_ele_max          = n_ele_use
  
  call allocate_branch_array (lat%branch, 0, lat)  ! Initial allocation
  lat%branch(0)%name = 'MAIN'

!---------------------------------------------------------------
! we now have the line to use in constructing the lat.
! now to put the elements in LAT in the correct order.
! superimpose, overlays, and groups are handled later.
! first load beam parameters.

  lat%version            = bmad_inc_version$
  lat%input_file_name    = full_lat_file_name             ! save input file  
  lat%n_ic_max           = 0                     
  lat%n_control_max      = 0    
  lat%param%growth_rate  = 0
  lat%param%stable            = .true.
  lat%param%aperture_limit_on = .true.

  lat%ele(0)     = in_lat%ele(0)    ! Beginning element
  lat%beam_start = in_lat%beam_start
  lat%a          = in_lat%a
  lat%b          = in_lat%b
  lat%z          = in_lat%z

  call mat_make_unit (lat%ele(0)%mat6)
  call clear_lat_1turn_mats (lat)

  if (beam_ele%value(n_part$) /= 0 .and. param_ele%value(n_part$) /= 0) &
            call warning ('BOTH "PARAMETER[N_PART]" AND "BEAM, N_PART" SET.')
  lat%param%n_part = max(beam_ele%value(n_part$), param_ele%value(n_part$))

  ix1 = nint(param_ele%value(particle$))
  ix2 = nint(beam_ele%value(particle$))
  if (ix1 /= positron$ .and. ix2 /= positron$) &
          call warning ('BOTH "PARAMETER[PARTICLE]" AND "BEAM, PARTICLE" SET.')
  lat%param%particle = ix1
  if (ix2 /=  positron$) lat%param%particle = ix2

! The lattice name from a "parameter[lattice] = ..." line is 
! stored the param_ele%descrip string

  if (associated(param_ele%descrip)) then
    lat%lattice = param_ele%descrip
    deallocate (param_ele%descrip)
  endif

! Work on multipass before overlays, groups, and usuperimpose. 
! This is important since the elements in the lattice get
! renamed and if not done first would confuse any overlays, girders, etc.
! Multipass elements are paired by multipass index and multipass line name

  do i = 1, lat%n_ele_track
    if (used_line(i)%ix_multipass == 0) cycle
    n_multi = 0  ! number of elements to slave together
    ix_multipass = used_line(i)%ix_multipass
    do j = i, lat%n_ele_track
      if (used_line(j)%ix_multipass == ix_multipass .and. &
          used_line(j)%multipass_line == used_line(i)%multipass_line) then
        n_multi = n_multi + 1
        ixm(n_multi) = j
        used_line(j)%ix_multipass = 0  ! mark as taken care of
      endif
    enddo
    call add_this_multipass (lat, ixm(1:n_multi))
  enddo

! Make sure that taylor order and lattice_type are not being set via
! the old way of doing things.

  do i = 1, bp_com%ivar_tot

    if (bp_com%var_name(i) == 'LATTICE_TYPE') then
      print *, 'BMAD_PARSER NOTE:'
      print *, '   DEPRECATED USE OF SYNTAX: "LATTICE_TYPE = ...".'
      print *, '   USE "PARAMETER[LATTICE_TYPE] = ..." SYNTAX INSTEAD.'
      lat%param%lattice_type = nint(bp_com%var_value(i))
    endif

    if (bp_com%var_name(i) == 'TAYLOR_ORDER') then
      print *, 'BMAD_PARSER NOTE:'
      print *, '   DEPRECATED USE OF SYNTAX: "TAYLOR_ORDER = ...".'
      print *, '   USE "PARAMETER[TAYLOR_ORDER] = ..." SYNTAX INSTEAD.'
      lat%input_taylor_order = nint(bp_com%var_value(i))
    endif

  enddo

! Set lattice_type.

  ix = nint(param_ele%value(lattice_type$))
  if (ix > 0) then  ! lattice_type has been set.
    lat%param%lattice_type = ix
  else              ! else use default
    lat%param%lattice_type = circular_lattice$      ! default 
    if (any(in_lat%ele(:)%key == lcavity$)) then    !   except...
      print *, 'BMAD_PARSER NOTE: THIS LATTICE HAS A LCAVITY.'
      print *, '     SETTING THE LATTICE_TYPE TO LINEAR_LATTICE.'
      lat%param%lattice_type = linear_lattice$
    endif
  endif

! Set taylor_order

  lat%input_taylor_order = nint(param_ele%value(taylor_order$))

  if (lat%input_taylor_order /= 0) &
       call set_taylor_order (lat%input_taylor_order, .false.)

!-------------------------------------------------------------------------
! energy bookkeeping.

  energy_beam  = 1d9 * beam_ele%value(energy_gev$)  
  energy_param = param_ele%value(e_tot$)
  energy_0     = lat%ele(0)%value(e_tot$) 

  if (energy_beam == 0 .and. energy_param == 0 .and. energy_0 == 0) then
    call out_io (s_warn$, r_name, 'e_tot IS 0!')
  elseif (energy_beam /= 0 .and. energy_param == 0 .and. energy_0 == 0) then
    lat%ele(0)%value(e_tot$) = energy_beam
  elseif (energy_beam == 0 .and. energy_param /= 0 .and. energy_0 == 0) then
    lat%ele(0)%value(e_tot$) = energy_param
  elseif (energy_beam == 0 .and. energy_param == 0 .and. energy_0 /= 0) then
    lat%ele(0)%value(e_tot$) = energy_0
  else
    call warning ('BEAM ENERGY SET MULTIPLE TIMES ')
  endif

  call convert_total_energy_to (lat%ele(0)%value(e_tot$), lat%param%particle, &
                                             pc = lat%ele(0)%value(p0c$))

  call set_ptc (lat%ele(0)%value(e_tot$), lat%param%particle)

! Go through the IN_LAT elements and put in the superpositions, groups, etc.
! First put in the superpositions and remove the null_ele elements

  call s_calc (lat)              ! calc longitudinal distances
  do i = 1, n_max
    if (in_lat%ele(i)%control_type /= super_lord$) cycle
    call add_all_superimpose (lat, in_lat%ele(i), plat%ele(i))
  enddo

  call remove_ele_from_lat (lat)  ! remove all null_ele elements.

! Add branch lines

  n0 = 0
  n = 0
  do
    branch0 => lat%branch(n0)
    do i = 1, branch0%n_ele_track
      if (branch0%ele(i)%key /= photon_branch$ .and. branch0%ele(i)%key /= branch$) cycle
      n = n + 1
      call allocate_branch_array (lat%branch, n)
      branch0 => lat%branch(n0)
      branch => lat%branch(n)
      branch%key = branch0%ele(i)%key
      branch%ix_branch = n
      branch%ix_from_branch = 0
      branch%ix_from_ele = i
      call parser_expand_line (n, lat%ele(i)%attribute_name, sequence, in_name, &
        in_indexx, seq_name, seq_indexx, in_lat%ele, branch%ele, used_line, n_ele_use)
      branch%ele(0)%key = init_ele$
      branch%ele(0)%name = 'BEGINNING'
      branch%name        = lat%ele(i)%attribute_name
      if (lat%ele(i)%alias /= '') branch%name = lat%ele(i)%alias
      branch%n_ele_track = n_ele_use
      branch%n_ele_max   = n_ele_use
    enddo
    n0 = n0 + 1
    if (n0 > n) exit
  enddo

! Now put in the overlay_lord, girder, and group elements

  call parser_add_lord (in_lat, n_max, plat, lat)

! global computations
! Reuse the old taylor series if they exist
! and the old taylor series has the same attributes.

  call lattice_bookkeeper (lat)
  lat%input_taylor_order = bmad_com%taylor_order

  call reuse_taylor_elements (lat, old_ele)

! make matrices for entire lat

  doit = .true.
  if (present(make_mats6)) doit = make_mats6
  if (doit) then
    call lat_make_mat6(lat)      ! make 6x6 transport matrices
  endif

! Do we need to expand the lattice and call bmad_parser2?

  8000 continue

  if (detected_expand_lattice_cmd) then
    exit_on_error = bmad_status%exit_on_error
    bmad_status%exit_on_error = .false.
    bp_com%bmad_parser_calling = .true.
    bp_com%old_lat => in_lat
    call bmad_parser2 ('FROM: BMAD_PARSER', lat)
    bp_com%bmad_parser_calling = .false.
    bmad_status%exit_on_error = exit_on_error
  endif

!-------------------------------------------------------------------------
! write out if debug is on

  if (debug_line /= '') call parser_debug_print_info (lat, debug_line)

! deallocate pointers

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

  if (associated (in_lat%ele))     call deallocate_lat_pointers (in_lat)
  if (associated (plat%ele))       deallocate (plat%ele)
  if (allocated (seq_indexx))      deallocate (seq_indexx, seq_name)
  if (allocated (in_indexx))       deallocate (in_indexx, in_name)
  if (allocated (in_lat%control))  deallocate (in_lat%control)
  if (allocated (in_lat%ic))       deallocate (in_lat%ic)
  if (allocated (used_line))       deallocate (used_line)

! error check

  if (bp_com%error_flag) then
    if (bmad_status%exit_on_error) then
       call out_io (s_fatal$, r_name, 'BMAD_PARSER FINISHED. EXITING ON ERRORS')
      stop
    else
      bmad_status%ok = .false.
      return
    endif
  endif
              
  call check_lat_controls (lat, .true.)

! write to digested file

  write_digested = write_digested .and. digested_version <= bmad_inc_version$
  if (write_digested) call write_digested_bmad_file (digested_file, &
                              lat, bp_com%num_lat_files, bp_com%lat_file_names)

end subroutine
