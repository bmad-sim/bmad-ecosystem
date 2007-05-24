!+
! Subroutine bmad_parser2 (lat_file, lat, orbit, make_mats6, 
!                                 digested_file_name, digested_read_ok)
!
! Subroutine parse (read in) a BMAD input file.
! This subrotine assumes that lat already holds an existing lattice.
! To read in a lattice from scratch use BMAD_PARSER or XSIF_PARSER.
!
! With BMAD_PARSER2 you may:
!     a) Modify the attributes of elements.
!     b) Define new overlays and groups.
!     c) Superimpose new elements upon the lattice.
!
! If the digested_file_name argument is present then bmad_parser2 will, If it exists
! and is up-to-date, set lat to the lattice in this digested file instead of modifying 
! lat using lat_file. If the digested file doesn't
! exist or is not up-to-date then a new digested file will be created. If this 
! argument is '*' then the new digested file name is:
!               'digested_' // lat_file   ! for single precision BMAD version
!               'digested8_' // lat_file  ! for double precision BMAD version
!
! Note: If you use the superimpose feature to insert an element into the latttice
!       then the index of a given element already in the lattice may change.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_file    -- Character(*): Input file name
!   lat         -- lat_struct: lattice with existing layout.
!   orbit(0:)  -- Coord_struct, optional: closed orbit for when
!                           bmad_parser2 calls lat_make_mat6
!   make_mats6  -- Logical, optional: Make the 6x6 transport matrices for then
!                   Elements? Default is True.
!   digested_file_name 
!               -- Character(*), optional: Name of a Bmad digested file to use.
!                           
!   digested_read_ok  
!               -- Logical, optinal: Set True if digested file was
!                           successfully read. False otherwise.
!
! Output:
!   lat    -- lat_struct: lattice with modifications.
!-

#include "CESR_platform.inc"

subroutine bmad_parser2 (lat_file, lat, orbit, make_mats6, &
                                           digested_file_name, digested_read_ok)

  use bmad_parser_mod, except => bmad_parser2

  implicit none
    
  type (lat_struct), target :: lat
  type (lat_struct), save :: lat2
  type (ele_struct), pointer :: ele
  type (ele_struct), target, save :: beam_ele, param_ele, beam_start_ele
  type (coord_struct), optional :: orbit(0:)
  type (parser_lat_struct) plat
  type (ele_struct), allocatable, save :: old_ele(:) 

  integer ix_word, i, ix, last_con, ixx, ele_num
  integer key, n_max_old, digested_version
  integer, pointer :: n_max

  character(*) lat_file
  character(*), optional :: digested_file_name
  character(1) delim 
  character(40) word_2, name
  character(16) :: r_name = 'bmad_parser2'
  character(32) word_1
  character(40) this_name
  character(280) parse_line_save, digested_name

  logical, optional :: make_mats6, digested_read_ok
  logical parsing, delim_found, found, doit
  logical file_end, err_flag, finished, write_digested

! init

  bmad_status%ok = .true.
  write_digested = .false.
  bp_com%parser_name = 'BMAD_PARSER2'
  call file_stack('push', lat_file, finished)   ! open file on stack
  if (.not. bmad_status%ok) return

  n_max => lat%n_ele_max
  n_max_old = n_max

  last_con = 0

  call allocate_plat (lat, plat)

  call init_ele(beam_ele)
  beam_ele%name = 'BEAM'              ! fake beam element
  beam_ele%key = def_beam$            ! "definition of beam"
  beam_ele%value(particle$)   = lat%param%particle 
  beam_ele%value(energy_gev$) = 0
  beam_ele%value(n_part$)     = lat%param%n_part

  call init_ele (param_ele)
  param_ele%name = 'PARAMETER'
  param_ele%key = def_parameter$
  param_ele%value(lattice_type$) = lat%param%lattice_type
  param_ele%value(E_TOT$)  = 0
  param_ele%value(taylor_order$) = lat%input_taylor_order
  param_ele%value(n_part$)       = lat%param%n_part

  call init_ele (beam_start_ele)
  beam_start_ele%name = 'BEAM_START'
  beam_start_ele%key = def_beam_start$

! see if a digested bmad file is available

  if (present(digested_file_name)) then
    digested_name = digested_file_name
    if (digested_file_name == '*') &
                  call form_digested_bmad_file_name (lat_file, digested_name)
    call read_digested_bmad_file (digested_name, lat2, digested_version)
    if (bmad_status%ok) then
      lat = lat2
      call deallocate_lat_pointers (lat2)
      return
    endif
    call save_taylor_elements (lat2, old_ele)
    call deallocate_lat_pointers (lat2)
    write_digested = .true.
    if (bmad_status%type_out) &
             call out_io (s_info$, r_name, 'Creating new digested file...')
  endif

  bmad_status%ok = .true.

!-----------------------------------------------------------
! main parsing loop

  bp_com%input_line_meaningful = .true.

  parsing_loop: do

! get a line from the input file and parse out the first word

    call load_parse_line ('normal', 1, file_end)  ! load an input line
    call get_next_word (word_1, ix_word, '[:](,)=', delim, delim_found, .true.)
    if (file_end) then
      word_1 = 'END_FILE'
      ix_word = 8
    else
      call verify_valid_name(word_1, ix_word)
    endif

! PARSER_DEBUG

    if (word_1(:ix_word) == 'PARSER_DEBUG') then
      bp_com%parser_debug = .true.
      bp_com%debug_line = bp_com%parse_line
      print *, 'FOUND IN FILE: "PARSER_DEBUG". DEBUG IS NOW ON'
      cycle parsing_loop
    endif

! NO_DIGESTED

    if (word_1(:ix_word) == 'NO_DIGESTED') then
      write_digested = .false.
      print *, 'FOUND IN FILE: "NO_DIGESTED". NO DIGESTED FILE WILL BE CREATED'
      cycle parsing_loop
    endif

! CALL command

    if (word_1(:ix_word) == 'CALL') then
      call get_called_file(delim)
      if (.not. bmad_status%ok) return
      cycle parsing_loop

    endif

! BEAM command

    if (word_1(:ix_word) == 'BEAM') then
      if (delim /= ',')  &
              call warning ('"BEAM" NOT FOLLOWED BY COMMA', ' ')

      parsing = .true.
      do while (parsing)
        if (.not. delim_found) then
          parsing = .false.
        elseif (delim /= ',') then
          call warning ('EXPECTING: "," BUT GOT: ' // delim,  &
                                             'FOR "BEAM" COMMAND')
          parsing = .false.
        else
          call get_attribute (def$, beam_ele, lat, plat, &
                                          delim, delim_found, err_flag)
          if (err_flag) cycle parsing_loop
        endif
      enddo

      cycle parsing_loop

    endif

! LATTICE command

    if (word_1(:ix_word) == 'LATTICE') then
      if ((delim /= ':' .or. bp_com%parse_line(1:1) /= '=') &
                                       .and. (delim /= '=')) then
        call warning ('"LATTICE" NOT FOLLOWED BY ":="', ' ')
      else
        if (delim == ':') bp_com%parse_line = bp_com%parse_line(2:)  ! trim off '='
        call get_next_word (lat%lattice, ix_word, ',', &
                                            delim, delim_found, .true.)
      endif
      cycle parsing_loop
    endif

! RETURN or END_FILE command

    if (word_1(:ix_word) == 'RETURN' .or.  &
                                    word_1(:ix_word) == 'END_FILE') then
      call file_stack ('pop', ' ', finished)
      if (.not. bmad_status%ok .and. bmad_status%exit_on_error) call err_exit
      if (.not. bmad_status%ok) return
      if (finished) then
        exit parsing_loop
      else
        cycle parsing_loop
      endif
    endif

! Variable definition or element redef...

! if an element attribute redef.

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
        if (lat%ele(i)%name == word_1 .or. &
                          key_name(lat%ele(i)%key) == word_1) then
          bp_com%parse_line = trim(word_2) // ' = ' // bp_com%parse_line 
          if (found) then   ! if not first time
            bp_com%parse_line = parse_line_save
          else
            parse_line_save = bp_com%parse_line
          endif
          call get_attribute (redef$, lat%ele(i), lat, plat, &
                                               delim, delim_found, err_flag)
          if (delim_found) call warning ('BAD DELIMITER: ' // delim, ' ')
          found = .true.
        endif
      enddo

      ! If bmad_parser2 has been called from bmad_parser then check if the
      ! element was just not used in the lattice. If so then just ignore it.
      if (.not. found) then
        if (bp_com%bmad_parser_calling) then
          do i = 0, bp_com%old_lat%n_ele_max
            if (bp_com%old_lat%ele(i)%name == word_1) then
              bp_com%parse_line = ' '  ! discard rest of statement
              cycle parsing_loop       ! goto next statement
            endif
          enddo
        endif
        call warning ('ELEMENT NOT FOUND: ' // word_1)
      endif

      cycle parsing_loop

! else must be a variable

    elseif (delim == '=') then

      call parser_add_variable (word_1, lat)
      cycle parsing_loop

    endif

! bad delimiter

    if (delim /= ':') then
      call warning ('1ST DELIMITER IS NOT ":". IT IS: ' // delim,  &
                                                       'FOR: ' // word_1)
      cycle parsing_loop
    endif

! only possibilities left are: element, list, or line
! to decide which look at 2nd word

    call get_next_word(word_2, ix_word, ':=,', delim, delim_found, .true.)
    if (ix_word == 0) then
      call error_exit ('NO NAME FOUND AFTER: ' // word_1, ' ')
    endif

    call verify_valid_name(word_2, ix_word)

! if line or list then this is an error for bmad_parser2

    if (word_2(:ix_word) == 'LINE' .or. word_2(:ix_word) == 'LIST') then
      call warning ('LINES OR LISTS NOT PERMITTED: ' // word_1, ' ')

! if not line or list then must be an element

    else

      n_max = n_max + 1
      if (n_max > ubound(lat%ele, 1)) then
        call allocate_lat_ele(lat)
        call allocate_plat (lat, plat)
      endif

      lat%ele(n_max)%name = word_1
      last_con = last_con + 1     ! next free slot
      lat%ele(n_max)%ixx = last_con

      do i = 1, n_max-1
        if (lat%ele(n_max)%name == lat%ele(i)%name) then
          call warning ('DUPLICATE ELEMENT NAME ' // lat%ele(n_max)%name, ' ')
          exit
        endif
      enddo

! Check for valid element key name or if element is part of a element key.
! If none of the above then we have an error.

      found = .false.  ! found a match?

      lat%ele(n_max)%key = key_name_to_key_index(word_2, .true.)
      if (lat%ele(n_max)%key > 0) then
        call preparse_element_init (lat%ele(n_max))
        found = .true.
      endif

      if (.not. found) then
        do i = 1, n_max-1
          if (word_2 == lat%ele(i)%name) then
            lat%ele(n_max) = lat%ele(i)
            lat%ele(n_max)%name = word_1
            found = .true.
            exit
          endif
        enddo
      endif

      if (.not. found) then
        call warning ('KEY NAME NOT RECOGNIZED: ' // word_2,  &
                       'FOR ELEMENT: ' // lat%ele(n_max)%name)
        lat%ele(n_max)%key = 1       ! dummy value
      endif

! now get the attribute values.
! For control elements lat%ele()%IXX temporarily points to
! the plat structure where storage for the control lists is
                   
      call parser_set_ele_defaults (lat%ele(n_max))

      key = lat%ele(n_max)%key
      if (key == overlay$ .or. key == group$ .or. key == i_beam$) then
        if (delim /= '=') then
          call warning ('EXPECTING: "=" BUT GOT: ' // delim,  &
                      'FOR ELEMENT: ' // lat%ele(n_max)%name)
        else
          if (key == overlay$) lat%ele(n_max)%control_type = overlay_lord$
          if (key == group$)   lat%ele(n_max)%control_type = group_lord$
          if (key == i_beam$)  lat%ele(n_max)%control_type = i_beam_lord$
          call get_overlay_group_names(lat%ele(n_max), lat,  &
                                              plat, delim, delim_found)
        endif
        if (key /= i_beam$ .and. .not. delim_found) then
          call warning ('NO CONTROL ATTRIBUTE GIVEN AFTER CLOSING "}"',  &
                        'FOR ELEMENT: ' // lat%ele(n_max)%name)
          n_max = n_max - 1
          cycle parsing_loop
        endif
      endif

      parsing = .true.
      do while (parsing)
        if (.not. delim_found) then          ! if nothing more
          parsing = .false.           ! break loop
        elseif (delim /= ',') then
          call warning ('EXPECTING: "," BUT GOT: ' // delim,  &
                        'FOR ELEMENT: ' // lat%ele(n_max)%name)
          n_max = n_max - 1
          cycle parsing_loop
        else
          call get_attribute (def$, lat%ele(n_max), &
                                  lat, plat, delim, delim_found, err_flag)
          if (err_flag) then
            n_max = n_max - 1
            cycle parsing_loop
          endif
        endif
      enddo

! Element must be a group, overlay, or superimpose element

      if (key /= overlay$ .and. key /= group$ .and. &
              lat%ele(n_max)%control_type /= super_lord$) then
        call warning ('ELEMENT MUST BE AN OVERLAY, SUPERIMPOSE, ' //  &
                                             'OR GROUP: ' // word_1, ' ')
        n_max = n_max - 1
        cycle parsing_loop
      endif

    endif

  enddo parsing_loop

!---------------------------------------------------------------
! Now we have read everything in

  bp_com%input_line_meaningful = .false.

  lat%param%particle    = nint(beam_ele%value(particle$))
  lat%param%lattice_type = nint(param_ele%value(lattice_type$))
  lat%input_taylor_order = nint(param_ele%value(taylor_order$))

  if (beam_ele%value(energy_gev$) /= 0) then
    lat%ele(0)%value(E_TOT$) = beam_ele%value(energy_gev$) * 1e9
  elseif (param_ele%value(E_TOT$) /= 0) then
    lat%ele(0)%value(E_TOT$) = param_ele%value(E_TOT$)
  endif

  if (lat%param%n_part /= param_ele%value(n_part$)) then
    lat%param%n_part = param_ele%value(n_part$)
  else
    lat%param%n_part = beam_ele%value(n_part$)
  endif

! Transfer the new elements to a safe_place

  ele_num = n_max - n_max_old
  allocate (lat2%ele(1:ele_num))
  lat2%ele(1:ele_num) = lat%ele(n_max_old+1:n_max)
  n_max = n_max_old

! Do bookkeeping for settable dependent variables.

  do i = 1, ele_num
    ele => lat2%ele(i)
    call settable_dep_var_bookkeeping (ele)
  enddo

! Put in the new elements...
! First put in superimpose elements

  do i = 1, ele_num
    ele => lat2%ele(i)
    if (ele%control_type /= super_lord$) cycle

    select case (ele%key)
    case (wiggler$)
      if (ele%sub_key == periodic_type$) then
        if (ele%value(l_pole$) == 0 .and. ele%value(n_pole$) /= 0) then
          ele%value(l_pole$) = ele%value(l$) / ele%value(n_pole$) 
        endif
      endif
    end select

    ixx = ele%ixx
    call add_all_superimpose (lat, ele, plat%ele(ixx))
  enddo

  call remove_all_null_ele_elements (lat)

! Go through and create the overlay, i_beam, and group lord elements.

  call parser_add_lord (lat2, ele_num, plat, lat)

! Reuse the old taylor series if they exist
! and the old taylor series has the same attributes.

  call reuse_taylor_elements (lat, old_ele)

! make matrices for entire lat

  call lattice_bookkeeper (lat)
  call compute_reference_energy (lat, .true.)
  doit = .true.
  if (present(make_mats6)) doit = make_mats6
  if (doit) call lat_make_mat6(lat, -1, orbit)  ! make transport matrices
  call s_calc (lat)                       ! calc loginitudinal distances
  call lat_geometry (lat)                ! lat layout

!-----------------------------------------------------------------------------
! error check

  if (bp_com%error_flag .and. bmad_status%exit_on_error) then
    print *, 'BMAD_PARSER2 FINISHED. EXITING ON ERRORS'
    stop
  endif

  call check_lat_controls (lat, .true.)

  do i = lbound(plat%ele, 1) , ubound(plat%ele, 1)
    if (associated (plat%ele(i)%name)) then
      deallocate(plat%ele(i)%name)
      deallocate(plat%ele(i)%attrib_name)
      deallocate(plat%ele(i)%coef)
    endif
  enddo

  if (associated (plat%ele))        deallocate (plat%ele)

  do i = 1, size(lat2%ele)
    call deallocate_ele_pointers (lat2%ele(i))
  enddo

! write to digested file

  if (write_digested .and. .not. bp_com%parser_debug .and. &
      digested_version <= bmad_inc_version$) call write_digested_bmad_file  &
             (digested_name, lat, bp_com%num_lat_files, bp_com%lat_file_names)

end subroutine
