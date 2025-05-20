!+
! Subroutine bmad_parser2 (lat_file, lat, orbit, make_mats6, err_flag, parse_lat)
!
! Subroutine parse (read in) a BMAD input file.
! This subrotine assumes that lat already holds an existing lattice.
! To read in a lattice from scratch use bmad_parser.
!
! With bmad_parser2 you may:
!     a) Modify the attributes of elements.
!     b) Define new overlays and groups.
!     c) Superimpose new elements upon the lattice.
!
! Note: Unlike bmad_parser, no digested file will be created.
!
! Note: If you use the superimpose feature to insert an element into the latttice
!       then the index of a given element already in the lattice may change.
!
! Input:
!   lat_file    -- Character(*): Input file name.
!   lat         -- lat_struct: lattice with existing layout.
!   orbit(0:)   -- Coord_struct, optional: closed orbit for when
!                           bmad_parser2 calls lat_make_mat6
!   make_mats6  -- Logical, optional: Make the 6x6 transport matrices for then
!                   Elements? Default is True.
!   parse_lat   -- lat_struct, optional: Used by bmad_parser to pass to bmad_parser2 a 
!                    list of elements that were defined in the lattice file but not used.
!                    This is useful in preventing errors being generated if group/overlay
!                    elements definded by lat_file refer to unused slaves in parse_lat.
!
! Output:
!   lat    -- lat_struct: lattice with modifications.
!-

subroutine bmad_parser2 (lat_file, lat, orbit, make_mats6, err_flag, parse_lat)

use parser_set_attribute_mod, except_dummy => bmad_parser2
use twiss_and_track_mod, except_dummy2 => bmad_parser2

implicit none
  
type (lat_struct), target :: lat
type (lat_struct), optional :: parse_lat
type (lat_struct) :: lat2
type (ele_struct), pointer :: ele, param_ele, lord, slave, slave2
type (ele_pointer_struct), allocatable :: eles(:)
type (parser_ele_struct), pointer :: pele
type (coord_struct), optional :: orbit(0:)
type (parser_lat_struct), target :: plat
type (branch_struct), pointer :: branch
type (coord_array_struct), allocatable :: orb_array(:)
real(rp) v1, v2

integer ix_word, i, j, n, ix, ix1, ix2, n_plat_ele, ixx, ix_word_1
integer key, n_max_old, n_loc, n_def_ele, is, is2, ib, ie, why_not_free, status
integer, pointer :: n_max

character(*) lat_file
character(1) delim 
character(16) :: r_name = 'bmad_parser2'
character(40) word_1, slice_start, slice_end, temp_ele_name
character(40) word_2, name, this_name, old_parser_name
character(80) debug_line
character(280) parse_line_save, string, extra_ele_names

logical, optional :: make_mats6, err_flag
logical parsing, found, delim_found, err, key_here, move
logical end_of_file, finished, good_attrib, wildcards_permitted, integer_permitted
logical multiple_eles_here, heterogeneous_ele_list

! Init...

if (present(err_flag)) err_flag = .true.
bp_com%write_digested2 = (.not. bp_com%always_parse)
old_parser_name = bp_com%parser_name
bp_com%parser_name = 'bmad_parser2'     ! Note: Code in parser_set_attribute depends upon this string being what it is.
bp_com%input_from_file = .true.
bp_com%fatal_error_flag = .false.       ! Set True on fatal (must abort now) error 
bp_com%calc_reference_orbit = .false.

debug_line = ''

! If lat_file = 'FROM: BMAD_PARSER' then bmad_parser2 has been called by 
! bmad_parser (after an expand_lattice command). 
! In this case we just read from the current open file.

if (lat_file /= 'FROM: BMAD_PARSER') then
  bp_com%do_superimpose = .true.
  call parser_file_stack('init')
  call parser_file_stack('push', lat_file, finished, err)   ! open file on stack
  if (err) return
endif

call init_bmad_parser_common(lat)

! Note: The order of def_parameter and def_mad_beam elements is used by parser_set_attribute
! due to overlapping parameters

n_max => lat%n_ele_max
n_def_ele = 7 + size(lat%branch)
call allocate_plat (plat, n_def_ele)
if (ubound(lat%ele, 1) < n_max + n_def_ele) call allocate_lat_ele_array(lat, n_max+n_def_ele+100, do_ramper_slave_setup = .true.)

ele => lat%ele(n_max+1)  ! Important: def_parameter comes after def_mad_beam.
call init_ele(ele, def_mad_beam$, 0, n_max+1, lat%branch(0))
ele%name = 'BEAM'              ! For MAD compatibility.
ele%value(n_part$)  = real_garbage$
ele%ref_species     = not_set$
ele%ixx = 1                    ! Pointer to plat%ele() array
extra_ele_names = ele%name 
call nametable_add(lat%nametable, ele%name, n_max+1)

ele => lat%ele(n_max+2)  ! Important: def_parameter comes after def_mad_beam
call init_ele (ele, def_parameter$, 0, n_max+2, lat%branch(0))
ele%name = 'PARAMETER'
ele%value(n_part$)    = real_garbage$
ele%ref_species       = not_set$
ele%ixx = 2                    ! Pointer to plat%ele() array
extra_ele_names = trim(extra_ele_names) // ', ' // ele%name
call nametable_add(lat%nametable, ele%name, n_max+2)

ele => lat%ele(n_max+3)
call init_ele (ele, def_particle_start$, 0, n_max+3, lat%branch(0))
ele%name = 'PARTICLE_START'
ele%ixx  = 3                    ! Pointer to plat%ele() array
extra_ele_names = trim(extra_ele_names) // ', ' // ele%name
call nametable_add(lat%nametable, ele%name, n_max+3)

ele => lat%ele(n_max+4)
call init_ele (ele, def_ptc_com$, 0, n_max+4, lat%branch(0))
ele%name = 'PTC_COM'
ele%ixx  = 4                    ! Pointer to plat%ele() array
extra_ele_names = trim(extra_ele_names) // ', ' // ele%name
call nametable_add(lat%nametable, ele%name, n_max+4)

ele => lat%ele(n_max+5)
call init_ele (ele, def_bmad_com$, 0, n_max+5, lat%branch(0))
ele%name = 'BMAD_COM'
ele%ixx  = 5                    ! Pointer to plat%ele() array
extra_ele_names = trim(extra_ele_names) // ', ' // ele%name
call nametable_add(lat%nametable, ele%name, n_max+5)

ele => lat%ele(n_max+6)
call init_ele (ele, def_space_charge_com$, 0, n_max+6, lat%branch(0))
ele%name = 'SPACE_CHARGE_COM'
ele%ixx  = 6                    ! Pointer to plat%ele() array
extra_ele_names = trim(extra_ele_names) // ', ' // ele%name
call nametable_add(lat%nametable, ele%name, n_max+6)

ele => lat%ele(n_max+7)
call init_ele (ele, null_ele$, 0, n_max+7, lat%branch(0))
ele%name = '<TEMP_ELE>'  ! Used for new element parameter storage
extra_ele_names = trim(extra_ele_names) // ', ' // ele%name
call nametable_add(lat%nametable, ele%name, n_max+7)
temp_ele_name = ele%name

do i = 0, ubound(lat%branch, 1)
  ele => lat%ele(n_max+8+i)
  call init_ele(ele, def_line$, 0, n_max+8+i, lat%branch(0))
  ele%name = lat%branch(i)%name
  ele%value(ix_branch$) = i
  ele%ixx = 8 + i
  extra_ele_names = trim(extra_ele_names) // ', ' // ele%name
  call nametable_add(lat%nametable, ele%name, n_max+8+i)
enddo

n_plat_ele = n_def_ele
n_max_old = n_max
n_max = n_max + n_def_ele

if (.not. bp_com%bmad_parser_calling) call load_parse_line ('init', 1, end_of_file)
call init_lat (lat2, 10)


!-----------------------------------------------------------
! main parsing loop

bp_com%input_line_meaningful = .true.

parsing_loop: do

  ! get a line from the input file and parse out the first word

  call load_parse_line ('new_command', 1, end_of_file)  ! load an input line
  call get_next_word (word_1, ix_word, '[:](,)= ', delim, delim_found, .true.)
  if (end_of_file) then
    word_1 = 'END_FILE'
    ix_word = 8
  endif

  ! If input line is something like "quadrupole::*[k1] = ..." or "q1:q2[x_limit] = ..." then shift delim from ":" to "[".
  ! Careful: Do not want to shift delim for something like "q:quad, l = d[l]"

  if (delim == ':' .or. delim == ',') then
    ix = index(bp_com%parse_line, '[')
    ix2 = index(bp_com%parse_line, '=')
    if (ix /= 0 .and. ix2 /= 0 .and. ix < ix2) then
      word_1 = trim(word_1) // ':' // bp_com%parse_line(:ix-1)
      bp_com%parse_line = bp_com%parse_line(ix+1:)
      delim = '['
      ix_word = len_trim(word_1)
    endif
  endif

  ! Name check. 
  ! If delim = '[' then have attribute redef and things are complicated so do not check.

  if (delim /= '[') then
    if (.not. verify_valid_name(word_1, ix_word)) cycle parsing_loop
  endif

  !-------------------------------------------
  ! CALC_REFERENCE_ORBIT

  if (word_1(:ix_word) == 'CALC_REFERENCE_ORBIT') then
    bp_com%calc_reference_orbit = .true.
    if (bp_com%parse_line /= '') then
      call parser_error ('EXTRA STUFF ON CALC_REFERENCE_ORBIT LINE: ' // bp_com%parse_line)
    endif
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! COMBINE_CONSECUTIVE_ELEMENTS

  if (word_1(:ix_word) == 'COMBINE_CONSECUTIVE_ELEMENTS') then
    call combine_consecutive_elements (lat, err)
    if (err) call parser_error ('ERROR COMBINING ELEMENTS.')
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! MERGE_ELEMENTS

  if (word_1(:ix_word) == 'MERGE_ELEMENTS' .or. word_1(:ix_word) == 'COMBINE_ELEMENTS') then
    string = trim(bp_com%parse_line) // ', ' // trim(extra_ele_names)
    call lat_ele_locator(string, lat, eles, n_loc, err)
    if (err) then
      call parser_error ('ERROR IN ELEMENT LIST FOR MERGE_ELEMENTS: ' // bp_com%parse_line)
      bp_com%parse_line = ''
      cycle parsing_loop
    endif

    bp_com%parse_line = ''

    do ib = 0, ubound(lat%branch, 1)
      lat%branch(ib)%ele%select = .false.
    enddo

    do i = 1, n_loc
      eles(i)%ele%select = .true.
    enddo

    if (bp_com%calc_reference_orbit) then
      call twiss_and_track(lat, orb_array, status)
      if (status /= ok$) then
        call parser_error ('ERROR CALCULATING REFERENCE ORBIT')
        cycle parsing_loop
      endif
      call make_hybrid_lat (lat, lat2, .true., orb_array)
    else
      call make_hybrid_lat (lat, lat2, .true.)
    endif

    lat = lat2
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! REMOVE_ELEMENTS

  if (word_1(:ix_word) == 'REMOVE_ELEMENTS') then
    call lat_ele_locator(bp_com%parse_line, lat, eles, n_loc, err)
    if (err) then
      call parser_error ('ERROR IN ELEMENT LIST FOR REMOVE_ELEMENTS: ' // bp_com%parse_line)
      bp_com%parse_line = ''
      cycle parsing_loop
    endif

    ! Need to call lattice_bookkeeper in the case where an overlay is to be removed since slave parameter
    ! values needs to be set first.
    call lattice_bookkeeper(lat)

    do i = 1, n_loc
      eles(i)%ele%ix_ele = -1
    enddo

    call remove_eles_from_lat(lat)
    bp_com%parse_line = ''

    cycle parsing_loop
  endif

  !-------------------------------------------
  ! PARSER_DEBUG

  if (word_1(:ix_word) == 'PARSER_DEBUG') then
    debug_line = bp_com%parse_line
    call out_io (s_info$, r_name, 'FOUND IN FILE: "PARSER_DEBUG". DEBUG IS NOW ON')
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! PRINT

  if (word_1(:ix_word) == 'PRINT') then
    call parser_print_line(lat, end_of_file)
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! NO_DIGESTED

  if (word_1(:ix_word) == 'NO_DIGESTED') then
    if (bp_com%write_digested2) call out_io (s_info$, r_name, 'Found in file: "NO_DIGESTED". No digested file will be created')
    bp_com%write_digested  = .false.
    bp_com%write_digested2 = .false.
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! SLICE_LATTICE

  if (word_1(:ix_word) == 'SLICE_LATTICE') then
    string = trim(bp_com%parse_line) // ', ' // trim(extra_ele_names)
    call slice_lattice (lat, string, err)
    if (err) call parser_error ('ERROR SLICING LATTICE USING: ' // bp_com%parse_line)
    bp_com%parse_line = ''
    cycle parsing_loop    
  endif

  !-------------------------------------------
  ! START_BRANCH_AT

  if (word_1(:ix_word) == 'START_BRANCH_AT') then
    move = .false.
    if (delim == ',') then
      call get_next_word (word_1, ix_word, ' ', delim, delim_found, .true.)
      if (word_1 /= 'MOVE_END_MARKER') then
        call parser_error('EXPECTED "MOVE_END_MARKER" AFTER COMMA IN START_BRANCH_AT STATEMENT.')
        cycle parsing_loop
      endif
      move = .true.
    endif

    if (delim /= ' ') then
      call parser_error('MALFORMED START_BRANCH_AT STATEMENT.')
      cycle parsing_loop
    endif

    call start_branch_at (lat, trim(bp_com%parse_line), move, err)
    if (err) call parser_error ('ERROR STARTING BRANCH AT: ' // bp_com%parse_line)
    bp_com%parse_line = ''
    cycle parsing_loop    
  endif

  !-------------------------------------------
  ! Superimpose statement

  if (word_1(:ix_word) == 'SUPERIMPOSE') then
    n_plat_ele = n_plat_ele + 1     ! next free slot
    if (n_plat_ele > ubound(plat%ele, 1)) call allocate_plat (plat, 2*size(plat%ele))
    pele => plat%ele(n_plat_ele)
    call parse_superimpose_command(lat, ele, pele, delim)
    call lat_ele_locator (pele%ele_name, lat, eles, n_loc, err)
    if (n_loc == 0 .or. err) then
      call parser_error ('CANNOT FIND ELEMENT FOR SUPERPOSITION: ' // pele%ele_name)
      cycle parsing_loop    
    endif
    call parser2_add_superimpose (lat, eles(1)%ele, pele, parse_lat)
    cycle parsing_loop    
  endif

  !-------------------------------------------
  ! WRITE_DIGESTED

  if (word_1(:ix_word) == 'WRITE_DIGESTED') then
    if (.not. bp_com%write_digested2) call out_io (s_info$, r_name, 'Found in file: "WRITE_DIGESTED". A digested file will be created')
    bp_com%write_digested = .true.
    bp_com%write_digested2 = .true.
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! CALL command

  if (word_1(:ix_word) == 'CALL') then
    call get_called_file(delim, string, err)
    if (err) return
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! BEAM command

  if (word_1(:ix_word) == 'BEAM') then
    if (delim /= ',')  call parser_error ('"BEAM" NOT FOLLOWED BY COMMA')

    parsing = .true.
    do while (parsing)
      if (.not. delim_found) then
        parsing = .false.
      elseif (delim /= ',') then
        call parser_error ('EXPECTING: "," BUT GOT: ' // delim, 'FOR "BEAM" COMMAND')
        parsing = .false.
      else
        call lat_ele_locator ('BEAM', lat, eles, n_loc, err)
        if (n_loc /= 1 .or. err) call err_exit
        call parser_set_attribute (def$, eles(1)%ele, delim, delim_found, err, check_free = .true.)
        if (err) cycle parsing_loop
      endif
    enddo

    cycle parsing_loop
  endif

  !-------------------------------------------
  ! EXPAND_LATTICE command: This does nothing since the lattice has already been expanded.

  if (word_1(:ix_word) == 'EXPAND_LATTICE') then
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! RETURN or END_FILE command

  if (word_1(:ix_word) == 'RETURN' .or.  word_1(:ix_word) == 'END_FILE') then
    call parser_file_stack ('pop', ' ', finished, err)
    if (err) return
    if (finished) then
      exit parsing_loop
    else
      cycle parsing_loop
    endif
  endif

  !---------------------------------------
  ! if an element attribute redef.

  if (delim == '[') then

    ! Old style "SBEND" is equiv to "SBEND::*"
    ix = index(word_1, '::')
    if (ix == 0 .and. key_name_to_key_index(word_1, .false.) > -1) word_1 = trim(word_1) // '::*'

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

    ! Find associated element and evaluate the attribute value

    if (word_1(1:1) == '+') then
      word_1 = word_1(2:)
      call parser_error ('"+" sign prefix construct deprecated. Please remove it.', level = s_warn$)
    endif

    parse_line_save = trim(word_2) // ' = ' // bp_com%parse_line 

    ! 

    if (word_1 == 'BEAM_START') word_1 = 'PARTICLE_START'   ! For backwards compatibility

    if (any(word_1 == key_name)) then   ! If Old style "quadrupole[k1] = ..." syntax
      name = word_1 // '::*'
      multiple_eles_here = .true.
    else
      name = word_1
    endif

    multiple_eles_here = (index(word_1, '*') /= 0 .or. index(word_1, '%') /= 0 .or. &
                 index(word_1, ',') /= 0 .or. index(word_1, ':') /= 0) ! Wild card character or range found
    key = 0
    ix = index(name, '::')
    if (ix /= 0) key = key_name_to_key_index (word_1(1:ix-1))
    ! OK if a homogeneous list is labeled heterogeneous. But not vice versa.
    heterogeneous_ele_list = (multiple_eles_here .and. index(word_1, ',') == 0 .and. key == 0)    

    ! When setting an attribute for all elements then suppress error printing.

    call lat_ele_locator (name, lat, eles, n_loc, err)
    if (err) then
      call parser_error ('MALFORMED ELEMENT(S) NAME: ' // quote(name))
      cycle parsing_loop
    endif

    found = .false.

    do i = 1, n_loc
      ele => eles(i)%ele

      ! No wild card matches permitted for these
      if (ele%name == 'BEGINNING' .or. ele%name == 'BEAM' .or. ele%name == 'SPACE_CHARGE_COM' .or. &
          ele%name == 'PARAMETER' .or. ele%name == 'PARTICLE_START' .or. ele%name == 'BMAD_COM') then
        if (word_1 /= ele%name) cycle
      endif

      ! When setting a hetero set of eles, cycle if not a valid attrib for this ele
      if (heterogeneous_ele_list .and. attribute_index(ele, word_2) == 0) cycle  

      if (multiple_eles_here .and. (ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$)) then
        if (.not. attribute_free(ele, word_2, .false., why_not_free = why_not_free)) then
          if (why_not_free == super_slave$ .or. why_not_free == multipass_slave$) then
            ele => pointer_to_lord(ele, 1)
            if (ele%slave_status == multipass_slave$) ele => pointer_to_lord(ele, 1)
          endif
        endif
      endif

      bp_com%parse_line = parse_line_save
      found = .true.
      call parser_set_attribute (redef$, ele, delim, delim_found, err, check_free = .true., heterogeneous_ele_list = heterogeneous_ele_list)
      if (.not. err .and. delim_found) then
        call parser_error ('BAD DELIMITER: ' // delim)
        exit
      endif

      call set_flags_for_changed_attribute (ele)
    enddo

    bp_com%parse_line = '' ! Needed if last call to parser_set_attribute did not have a set.

    ! If bmad_parser2 has been called from bmad_parser then check if the
    ! element was just not used in the lattice. If so then just ignore it.

    if (.not. found .and. .not. multiple_eles_here) then
      if (bp_com%bmad_parser_calling) then
        do i = 0, bp_com%old_lat%n_ele_max
          if (bp_com%old_lat%ele(i)%name == word_1) cycle parsing_loop
        enddo
      endif

      call parser_error ('ELEMENT NOT FOUND: ' // word_1)
    endif

    cycle parsing_loop

  endif

  !---------------------------------------
  ! Variable definition 

  if (word_1 == 'REDEF' .and. delim == ':') then
    call get_next_word (word_1, ix_word, '[:](,)= ', delim, delim_found, .true.)
    call parser_add_constant (word_1, lat, .false.)
    cycle parsing_loop
  endif

  if (delim == '=') then
    call parser_add_constant (word_1, lat, .true.)
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! bad delimiter

  if (delim /= ':') then
    if (delim == '') then
      call parser_error ('UNRECOGNIZED COMMAND: ' // word_1)
    else
      call parser_error ('EXPECTED DELIMITER TO BE ":". BUT IT IS: "' // trim(delim),  '" AFTER WORD: ' // word_1)
    endif
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! Only possibilities left are: element, list, or line
  ! to decide which look at 2nd word

  call get_next_word(word_2, ix_word, ':=,', delim, delim_found, .true.)
  if (ix_word == 0) then
    call parser_error ('NOTHING FOUND AFTER: ' // word_1)
    if (global_com%exit_on_error) call err_exit
  endif

  if (.not. verify_valid_name(word_2, ix_word)) cycle parsing_loop

  !-------------------------------------------------------
  ! if line or list then this is an error for bmad_parser2

  if (word_2(:ix_word) == 'LINE' .or. word_2(:ix_word) == 'LIST') then
    call parser_error ('DEFINING LINES OR LISTS NOT PERMITTED AFTER LATTICE EXPLANSION: ' // word_1, &
                       'I DO NOT KNOW TO DO WITH THIS SINCE THE LATTICE ALREADY EXISTS!')
    cycle parsing_loop
  endif

  !-------------------------------------------------------
  ! if none of the above then must be an element

  call lat_ele_locator (word_1, lat, eles, n_loc, err)
  if (n_loc > 0) then
    call parser_error ('DUPLICATE ELEMENT NAME ' // ele%name)
    exit
  endif

  ! Element must be stored tempararily in lat to keep set_element_attribute happy.

  call lat_ele_locator (temp_ele_name, lat, eles, n_loc, err)
  ele => eles(1)%ele
  call init_ele(ele, null_ele$, 0, ele%ix_ele, lat%branch(0))
  ele%name = word_1

  n_plat_ele = n_plat_ele + 1     ! next free slot
  ele%ixx = n_plat_ele
  if (n_plat_ele > ubound(plat%ele, 1)) call allocate_plat (plat, 2*size(plat%ele))
  pele => plat%ele(n_plat_ele)

  pele%lat_file = bp_com%current_file%full_name
  pele%ix_line_in_file = bp_com%current_file%i_line

  ! Check for valid element key name or if element is part of a element key.
  ! If none of the above then we have an error.

  found = .false.  ! found a match?

  call lat_ele_locator (word_2, lat, eles, n_loc, err)
  if (n_loc > 0) then
    ixx = ele%ixx  ! save
    ele = eles(1)%ele
    ele%ixx = ixx   ! Restore correct value
    call set_ele_name (ele, word_1)
    found = .true.
    exit
  endif

  if (.not. found) then
    ele%key = key_name_to_key_index(word_2, .true.)
    if (key_name_to_key_index(word_1, .false.) > 0) then
      call parser_error ('ELEMENT NAME: ' // word_1, &
                         'IS NOT ALLOWED TO BE THE SAME AS AN ELEMENT CLASS: ' // word_2)
    endif

    if (ele%key > 0) then
      call set_ele_defaults (ele)
      found = .true.
    endif
  endif

  if (.not. found) then
    call parser_error ('KEY NAME NOT RECOGNIZED OR AMBIGUOUS: ' // word_2,  &
                       'FOR ELEMENT: ' // ele%name)
    ele%key = 1       ! dummy value
  endif

  ! now get the attribute values.
  ! For control elements lat%ele()%ixx temporarily points to
  ! the plat structure where storage for the control lists is
               
  key = ele%key
  if (key == overlay$ .or. key == group$ .or. key == girder$ .or. key == ramper$) then
    if (delim /= '=') then
      call parser_error ('EXPECTING: "=" BUT GOT: ' // delim,  &
                  'FOR ELEMENT: ' // ele%name)

    else
      select case (key)
      case (overlay$); ele%lord_status = overlay_lord$
      case (group$);   ele%lord_status = group_lord$
      case (girder$);  ele%lord_status = girder_lord$
      case (ramper$);  ele%lord_status = ramper_lord$
      end select

      call get_overlay_group_names(ele, lat,  pele, delim, delim_found, .false., err)
      if (err) cycle parsing_loop
    endif
    if (key /= girder$ .and. .not. delim_found) then
      call parser_error ('NO CONTROL ATTRIBUTE GIVEN AFTER CLOSING "}"',  &
                    'FOR ELEMENT: ' // ele%name)
      cycle parsing_loop
    endif
  endif

  parsing = .true.
  do while (parsing)
    if (.not. delim_found) then          ! if nothing more
      parsing = .false.           ! break loop
    elseif (delim /= ',') then
      call parser_error ('EXPECTING: "," BUT GOT: ' // delim,  &
                    'FOR ELEMENT: ' // ele%name)
      cycle parsing_loop
    else
      call parser_set_attribute (def$, ele, delim, delim_found, err, pele)
      call set_flags_for_changed_attribute (ele)
      if (err) then
        cycle parsing_loop
      endif
    endif
  enddo

  ! Element must be a group, overlay, or superimpose element

  call settable_dep_var_bookkeeping (ele)

  select case (ele%key)
  case (wiggler$, undulator$)
    if (ele%field_calc == planar_model$ .or. ele%field_calc == helical_model$) then
      if (ele%value(l_period$) == 0 .and. ele%value(n_period$) /= 0) then
        ele%value(l_period$) = ele%value(l$) / ele%value(n_period$) 
      endif
    endif
  end select

  call drift_multipass_name_correction(lat) ! In case superimposing upon multipass elements.

  lat2%ele(1) = ele
  call set_ele_name (ele, temp_ele_name)
  ele => lat2%ele(1) ! Needed since call to add superimpose will shift elements in lat.

  if (ele%lord_status == super_lord$) then
    ixx = ele%ixx
    call parser2_add_superimpose (lat, lat2%ele(1), plat%ele(ixx), parse_lat)

  elseif (key == overlay$ .or. key == group$ .or. key == girder$ .or. key == ramper$) then
    call parser_add_lords (lat2, 1, plat, lat, parse_lat)

  else
    call parser_error ('ERROR FOR ELEMENT: ' // word_1, &
                       'ELEMENTS DEFINED AFTER LATTICE EXPANSION MUST BE SUPERIMPOSED OR BE A CONTROLLER', &
                       'TYPE ELEMENT (OVERLAY, GROUP, GIRDER, OR RAMPER ELEMENT).')
    cycle parsing_loop
  endif

  pele => plat%ele(ele%ixx)
  if (.not. allocated(pele%field_overlaps)) cycle
  call lat_ele_locator (ele%name, lat, eles, n)
  if (n < 1) cycle

  do j = 1, size(pele%field_overlaps)
    call create_field_overlap (lat, ele%name, pele%field_overlaps(j), err)
    if (err) then
      call parser_error ('CANNOT FIND ELEMENT: ' // pele%field_overlaps(j), &
                         'WHICH HAS FIELD OVERLAP FROM ELEMENT: ' // ele%name)
    endif
  enddo

enddo parsing_loop

call adjust_super_slave_names (lat, lat%n_ele_track+1, lat%n_ele_max)


!---------------------------------------------------------------
! Now we have read everything in

bp_com%input_line_meaningful = .false.

lat%input_taylor_order = bmad_com%taylor_order
call set_ptc()   ! Will set Taylor_order

! Must deal with fact that some parameters (geometry, etc) can be set in different places.

call lat_ele_locator ('BEAM', lat, eles, n_loc, err)
if (n_loc /= 1 .or. err) call err_exit
call lat_ele_locator ('PARAMETER', lat, eles, n_loc, err)
if (n_loc /= 1 .or. err) call err_exit
param_ele => eles(1)%ele

if (param_ele%value(geometry$) /= real_garbage$)    lat%param%geometry = nint(param_ele%value(geometry$))
if (param_ele%value(live_branch$) /= real_garbage$) lat%param%live_branch = is_true(param_ele%value(live_branch$))
if (param_ele%value(default_tracking_species$) /= real_garbage$) lat%param%default_tracking_species = nint(param_ele%value(default_tracking_species$))
if (param_ele%value(high_energy_space_charge_on$) /= real_garbage$ .and. is_true(param_ele%value(high_energy_space_charge_on$))) &
                                                                                                         bmad_com%high_energy_space_charge_on = .true.

if (param_ele%ref_species /= not_set$)     lat%param%particle = param_ele%ref_species

if (associated(param_ele%descrip)) then
  lat%lattice = param_ele%descrip
  deallocate (param_ele%descrip)
endif

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  ele => branch%ele(0)
  if (ele%value(e_tot$) < 0) then
    call convert_pc_to (ele%value(p0c$), branch%param%particle, e_tot = ele%value(e_tot$))
  elseif (ele%value(p0c$) < 0) then
    call convert_total_energy_to (ele%value(e_tot$), branch%param%particle, pc = ele%value(p0c$))
  endif
enddo

if (param_ele%value(n_part$) /= real_garbage$)     lat%param%n_part = param_ele%value(n_part$)

! Remove BEAM and other elements that were temporarily put in as well as null elements.

call lat_ele_locator(extra_ele_names, lat, eles, n_loc, err)
do i = 1, n_loc
  eles(i)%ele%ix_ele = -1
enddo

do i = 1, lat%n_ele_max
  if (lat%ele(i)%key == null_ele$) lat%ele(i)%ix_ele = -1 ! mark for deletion
enddo

call remove_eles_from_lat(lat)
call ramper_slave_setup(lat)

! Some bookkeeping

call drift_and_pipe_track_methods_adjustment(lat)

do ib = 0, ubound(lat%branch, 1)
  ele => lat%branch(ib)%ele(0)
  call floor_angles_to_w_mat(ele%floor%theta, ele%floor%phi, ele%floor%psi, ele%floor%w)
enddo

call apply_all_rampers(lat, err)
if (err) call parser_error ('ERROR APPLYING RAMPERS')

call set_flags_for_changed_attribute(lat)
call lattice_bookkeeper (lat)

! make matrices for entire lat

if (logic_option (.true., make_mats6)) call lat_make_mat6(lat, -1, orbit) 

!-----------------------------------------------------------------------------
! error check

if (debug_line /= '') call parser_debug_print_info (lat, debug_line)

call lat_sanity_check (lat, err)
if (err) bp_com%error_flag = .true.

if (bp_com%error_flag .and. global_com%exit_on_error) then
  call out_io (s_info$, r_name, 'FINISHED. EXITING ON ERRORS')
  stop
endif

call deallocate_lat_pointers (lat2)

if (present(err_flag)) err_flag = bp_com%error_flag

bp_com%parser_name = old_parser_name

! And cleanup

if (lat_file /= 'FROM: BMAD_PARSER') then
  if (allocated (bp_com%lat_file_names)) deallocate (bp_com%lat_file_names)
endif

end subroutine
