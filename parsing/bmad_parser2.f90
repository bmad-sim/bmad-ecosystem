!+
! Subroutine bmad_parser2 (lat_file, lat, orbit, make_mats6, err_flag)
!
! Subroutine parse (read in) a BMAD input file.
! This subrotine assumes that lat already holds an existing lattice.
! To read in a lattice from scratch use bmad_parser or xsif_parser.
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
! Modules needed:
!   use bmad
!
! Input:
!   lat_file    -- Character(*): Input file name.
!   lat         -- lat_struct: lattice with existing layout.
!   orbit(0:)   -- Coord_struct, optional: closed orbit for when
!                           bmad_parser2 calls lat_make_mat6
!   make_mats6  -- Logical, optional: Make the 6x6 transport matrices for then
!                   Elements? Default is True.
!
! Output:
!   lat    -- lat_struct: lattice with modifications.
!-

subroutine bmad_parser2 (lat_file, lat, orbit, make_mats6, err_flag)

use bmad_parser_mod, except_dummy => bmad_parser2

implicit none
  
type (lat_struct), target :: lat
type (lat_struct) :: lat2
type (ele_struct), pointer :: ele, mad_beam_ele, param_ele, lord, slave, slave2
type (ele_pointer_struct), allocatable :: eles(:)
type (parser_ele_struct), pointer :: pele
type (coord_struct), optional :: orbit(0:)
type (parser_lat_struct), target :: plat
type (branch_struct), pointer :: branch

real(rp) v1, v2

integer ix_word, i, j, n, ix, ix1, ix2, n_plat_ele, ixx, ele_num, ix_word_1
integer key, n_max_old, n_loc, n_def_ele, is, is2, ib, ie
integer, pointer :: n_max, n_ptr
integer, allocatable :: lat_indexx(:)

character(*) lat_file
character(1) delim 
character(16) :: r_name = 'bmad_parser2'
character(32) word_1
character(40) word_2, name, this_name
character(40), allocatable :: lat_name(:)
character(80) debug_line
character(280) parse_line_save, call_file

logical, optional :: make_mats6, err_flag
logical parsing, found, delim_found, xsif_called, err, key_here
logical end_of_file, finished, good_attrib, wildcards_permitted, integer_permitted
logical wild_here, wild_and_key0

! Init...

if (present(err_flag)) err_flag = .true.
bp_com%write_digested2 = (.not. bp_com%always_parse)
bp_com%parser_name = 'bmad_parser2'
bp_com%input_from_file = .true.
bp_com%fatal_error_flag = .false.       ! Set True on fatal (must abort now) error 
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

! Mark multipass elements in case there is a superposition

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  branch%ele%iyy = 0
enddo

do ie = lat%n_ele_track+1, lat%n_ele_max
  ele => lat%ele(ie)
  if (ele%lord_status /= multipass_lord$) cycle
  do is = 1, ele%n_slave
    slave => pointer_to_slave(ele, is)
    slave%iyy = ie
    if (slave%lord_status /= super_lord$) cycle
    do is2 = 1, slave%n_slave
      slave2 => pointer_to_slave(slave, is2)
      slave2%iyy = ie
    enddo
  enddo
enddo

!

n_max => lat%n_ele_max
n_def_ele = 5 + size(lat%branch)
call allocate_plat (plat, n_def_ele)
if (ubound(lat%ele, 1) < n_max + n_def_ele) call allocate_lat_ele_array(lat, n_max+n_def_ele+100)

ele => lat%ele(n_max+1)
call init_ele(ele, def_mad_beam$, 0, n_max+1, lat%branch(0))
ele%name = 'BEAM'              ! For MAD compatibility.
ele%value(n_part$)     = lat%param%n_part
ele%value(particle$)   = lat%param%particle
ele%ixx = 1                    ! Pointer to plat%ele() array

ele => lat%ele(n_max+2)
call init_ele (ele, def_parameter$, 0, n_max+2, lat%branch(0))
ele%name = 'PARAMETER'
ele%value(geometry$)     = lat%param%geometry
ele%value(live_branch$)  = logic_to_int(lat%param%live_branch)
ele%value(n_part$)       = lat%param%n_part
ele%value(particle$)     = lat%param%particle
ele%ixx = 2                    ! Pointer to plat%ele() array

ele => lat%ele(n_max+3)
call init_ele (ele, def_beam_start$, 0, n_max+3, lat%branch(0))
ele%name = 'BEAM_START'
ele%ixx  = 3                    ! Pointer to plat%ele() array

ele => lat%ele(n_max+4)
call init_ele (ele, def_bmad_com$, 0, n_max+4, lat%branch(0))
ele%name = 'BMAD_COM'
ele%ixx  = 4                    ! Pointer to plat%ele() array

ele => lat%ele(n_max+5)
call init_ele (ele, beginning_ele$, 0, n_max+5, lat%branch(0))
ele%name = 'BEGINNING'
ele%ixx  = 5                    ! Pointer to plat%ele() array

do i = 0, ubound(lat%branch, 1)
  ele => lat%ele(n_max+6+i)
  call init_ele(ele, line_ele$, 0, n_max+6+i, lat%branch(0))
  ele%name = lat%branch(i)%name
  ele%ixx = 6 + i
enddo

n_plat_ele = n_def_ele
n_max_old = n_max
n_max = n_max + n_def_ele

if (.not. bp_com%bmad_parser_calling) call load_parse_line ('init', 1, end_of_file)

!-----------------------------------------------------------
! main parsing loop

bp_com%input_line_meaningful = .true.

parsing_loop: do

  ! get a line from the input file and parse out the first word

  call load_parse_line ('new_command', 1, end_of_file)  ! load an input line
  call get_next_word (word_1, ix_word, '[:](,)=', delim, delim_found, .true.)
  if (end_of_file) then
    word_1 = 'END_FILE'
    ix_word = 8
  endif

  ! If input line is something like "quadrupole::*[k1] = ..." then shift delim from ":" to "["

  if (delim == ':' .and. bp_com%parse_line(1:1) == ':') then
    ix = index(bp_com%parse_line, '[')
    if (ix /= 0) then
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
  ! PARSER_DEBUG

  if (word_1(:ix_word) == 'PARSER_DEBUG') then
    debug_line = bp_com%parse_line
    call out_io (s_info$, r_name, 'FOUND IN FILE: "PARSER_DEBUG". DEBUG IS NOW ON')
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! PRINT

  if (word_1(:ix_word) == 'PRINT') then
    call string_trim (bp_com%input_line2, parse_line_save, ix) ! so can strip off initial "print"
    n_ptr => bp_com%num_lat_files
    if (size(bp_com%lat_file_names) < n_ptr + 1) call re_allocate(bp_com%lat_file_names, n_ptr+100)
    n_ptr = n_ptr + 1
    bp_com%lat_file_names(n_ptr) = '!PRINT:' // trim(parse_line_save(ix+2:)) ! To save in digested
    call out_io (s_dwarn$, r_name, 'Print Message in Lattice File: ' // parse_line_save(ix+2:))
    ! This prevents bmad_parser from thinking print string is a command.
    call load_parse_line ('init', 1, end_of_file)
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! NO_DIGESTED

  if (word_1(:ix_word) == 'NO_DIGESTED') then
    bp_com%write_digested  = .false.
    bp_com%write_digested2 = .false.
    call out_io (s_info$, r_name, 'FOUND IN FILE: "NO_DIGESTED". NO DIGESTED FILE WILL BE CREATED')
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! CALL command

  if (word_1(:ix_word) == 'CALL') then
    call get_called_file(delim, call_file, xsif_called, err)
    if (err) return
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! BEAM command

  if (word_1(:ix_word) == 'BEAM') then
    if (delim /= ',')  call parser_error ('"BEAM" NOT FOLLOWED BY COMMA', ' ')

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
        call parser_set_attribute (def$, eles(1)%ele, lat, delim, delim_found, err, check_free = .true.)
        if (err) cycle parsing_loop
      endif
    enddo

    cycle parsing_loop
  endif

  !-------------------------------------------
  ! LATTICE command

  if (word_1(:ix_word) == 'LATTICE') then
    if ((delim /= ':' .or. bp_com%parse_line(1:1) /= '=') .and. (delim /= '=')) then
      call parser_error ('"LATTICE" NOT FOLLOWED BY ":="', ' ')
    else
      if (delim == ':') bp_com%parse_line = bp_com%parse_line(2:)  ! trim off '='
      call get_next_word (lat%lattice, ix_word, ',', delim, delim_found, .true.)
    endif
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

    ! Find associated element and evaluate the attribute value

    if (word_1(1:1) == '+') then
      word_1 = word_1(2:)
      call parser_error ('"+" sign prefix construct deprecated. Please remove it.', level = s_warn$)
    endif

    parse_line_save = trim(word_2) // ' = ' // bp_com%parse_line 

    ! 

    if (any(word_1 == key_name)) then   ! If Old style "quadrupole[k1] = ..." syntax
      name = word_1 // '::*'
      wild_here = .true.
    else
      name = word_1
    endif

    wild_here = (index(word_1, '*') /= 0 .or. index(word_1, '%') /= 0) ! Wild card character found
    key = 0
    ix = index(name, '::')
    if (ix /= 0) key = key_name_to_key_index (word_1(1:ix-1))
    wild_and_key0 = (wild_here .and. key == 0)

    ! When setting an attribute for all elements then suppress error printing.

    call lat_ele_locator (name, lat, eles, n_loc, err)
    if (err) then
      bp_com%error_flag = .true.
      bp_com%parse_line = ''
      cycle parsing_loop
    endif

    found = .false.

    do i = 1, n_loc
      ele => eles(i)%ele

      ! No wild card matches permitted for these
      if (ele%name == 'BEGINNING' .or. ele%name == 'BEAM' .or. &
          ele%name == 'PARAMETER' .or. ele%name == 'BEAM_START' .or. ele%name == 'BMAD_COM') then
        if (word_1 /= ele%name) cycle
      endif

      if (wild_and_key0 .and. attribute_index(ele, word_2) == 0) cycle

      bp_com%parse_line = parse_line_save
      found = .true.
      call parser_set_attribute (redef$, ele, lat, delim, delim_found, err, check_free = .true., wild_and_key0 = wild_and_key0)
      if (.not. err .and. delim_found) then
        call parser_error ('BAD DELIMITER: ' // delim)
        exit
      endif

      call set_flags_for_changed_attribute (ele)
    enddo

    bp_com%parse_line = '' ! Needed if last call to parser_set_attribute did not have a set.

    ! If bmad_parser2 has been called from bmad_parser then check if the
    ! element was just not used in the lattice. If so then just ignore it.

    if (.not. found .and. .not. wild_here) then
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

  if (delim == '=') then
    call parser_add_variable (word_1, lat)
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! bad delimiter

  if (delim /= ':') then
    call parser_error ('1ST DELIMITER IS NOT ":". IT IS: ' // delim,  'FOR: ' // word_1)
    cycle parsing_loop
  endif

  !-------------------------------------------
  ! Only possibilities left are: element, list, or line
  ! to decide which look at 2nd word

  call get_next_word(word_2, ix_word, ':=,', delim, delim_found, .true.)
  if (ix_word == 0) then
    call parser_error ('NO NAME FOUND AFTER: ' // word_1, ' ')
    if (global_com%exit_on_error) call err_exit
  endif

  if (.not. verify_valid_name(word_2, ix_word)) cycle parsing_loop

  !-------------------------------------------------------
  ! if line or list then this is an error for bmad_parser2

  if (word_2(:ix_word) == 'LINE' .or. word_2(:ix_word) == 'LIST') then
    call parser_error ('LINES OR LISTS NOT PERMITTED: ' // word_1, ' ')
    cycle parsing_loop
  endif

  !-------------------------------------------------------
  ! if none of the above then must be an element

  n_max = n_max + 1
  if (n_max > ubound(lat%ele, 1)) call allocate_lat_ele_array(lat)
  ele => lat%ele(n_max)

  ele%name = word_1

  n_plat_ele = n_plat_ele + 1     ! next free slot
  ele%ixx = n_plat_ele
  if (n_plat_ele > ubound(plat%ele, 1)) call allocate_plat (plat, 2*size(plat%ele))
  pele => plat%ele(n_plat_ele)

  pele%lat_file = bp_com%current_file%full_name
  pele%ix_line_in_file = bp_com%current_file%i_line

  do i = 1, n_max-1
    if (ele%name == lat%ele(i)%name) then
      call parser_error ('DUPLICATE ELEMENT NAME ' // ele%name, ' ')
      exit
    endif
  enddo

  ! Check for valid element key name or if element is part of a element key.
  ! If none of the above then we have an error.

  found = .false.  ! found a match?

  do i = 1, n_max-1
    if (word_2 == lat%ele(i)%name) then
      ixx = ele%ixx  ! save
      ele = lat%ele(i)
      ele%ixx = ixx   ! Restore correct value
      ele%name = word_1
      found = .true.
      exit
    endif
  enddo

  if (.not. found) then
    ele%key = key_name_to_key_index(word_2, .true.)
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
  if (key == overlay$ .or. key == group$ .or. key == girder$) then
    if (delim /= '=') then
      call parser_error ('EXPECTING: "=" BUT GOT: ' // delim,  &
                  'FOR ELEMENT: ' // ele%name)
    else
      if (key == overlay$) ele%lord_status = overlay_lord$
      if (key == group$)   ele%lord_status = group_lord$
      if (key == girder$)  ele%lord_status = girder_lord$
      call get_overlay_group_names(ele, lat,  pele, delim, delim_found, .false.)
    endif
    if (key /= girder$ .and. .not. delim_found) then
      call parser_error ('NO CONTROL ATTRIBUTE GIVEN AFTER CLOSING "}"',  &
                    'FOR ELEMENT: ' // ele%name)
      n_max = n_max - 1
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
      n_max = n_max - 1
      cycle parsing_loop
    else
      call parser_set_attribute (def$, ele, lat, delim, delim_found, err, pele)
      call set_flags_for_changed_attribute (ele)
      if (err) then
        n_max = n_max - 1
        cycle parsing_loop
      endif
    endif
  enddo

  ! Element must be a group, overlay, or superimpose element

  if (key /= overlay$ .and. key /= group$ .and. ele%lord_status /= super_lord$) then
    call parser_error ('ELEMENT MUST BE AN OVERLAY, SUPERIMPOSE, ' //  &
                                         'OR GROUP: ' // word_1, ' ')
    n_max = n_max - 1
    cycle parsing_loop
  endif

enddo parsing_loop

!---------------------------------------------------------------
! Now we have read everything in

bp_com%input_line_meaningful = .false.

call lat_ele_locator ('BEAM', lat, eles, n_loc, err)
if (n_loc /= 1 .or. err) call err_exit
mad_beam_ele => eles(1)%ele
call lat_ele_locator ('PARAMETER', lat, eles, n_loc, err)
if (n_loc /= 1 .or. err) call err_exit
param_ele => eles(1)%ele

lat%param%geometry = nint(param_ele%value(geometry$))
lat%param%live_branch = is_true(param_ele%value(live_branch$))
lat%input_taylor_order = bmad_com%taylor_order

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


v1 = param_ele%value(n_part$)
v2 = mad_beam_ele%value(n_part$)
if (lat%param%n_part /= v1 .and. lat%param%n_part /= v2) then
  call parser_error ('BOTH "PARAMETER[N_PART]" AND "BEAM, N_PART" SET.')
else if (v1 /= lat%param%n_part) then
  lat%param%n_part = v1
else
  lat%param%n_part = v2
endif

ix1 = nint(param_ele%value(particle$))
ix2 = nint(mad_beam_ele%value(particle$))
if (ix1 /= lat%param%particle .and. ix2 /= lat%param%particle) &
        call parser_error ('BOTH "PARAMETER[PARTICLE]" AND "BEAM, PARTICLE" SET.')
lat%param%particle = ix1
if (ix2 /=  lat%param%particle) lat%param%particle = ix2

! Spin

call lat_ele_locator ('BEAM_START', lat, eles, n_loc, err)
if (n_loc /= 1 .or. err) call err_exit
call parser_set_spin (eles(1)%ele, lat%beam_start)

! Transfer the new elements to a safe_place and reset lat%n_max

ele_num = n_max - n_max_old - n_def_ele
allocate (lat2%ele(1:ele_num))
lat2%ele(1:ele_num) = lat%ele(n_max_old+n_def_ele+1:n_max)
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
  if (ele%lord_status /= super_lord$) cycle

  select case (ele%key)
  case (wiggler$, undulator$)
    if (ele%sub_key == periodic_type$) then
      if (ele%value(l_pole$) == 0 .and. ele%value(n_pole$) /= 0) then
        ele%value(l_pole$) = ele%value(l$) / ele%value(n_pole$) 
      endif
    endif
  end select

  ixx = ele%ixx

  do j = 0, ubound(lat%branch, 1)
    call parser_add_superimpose (lat%branch(j), ele, plat%ele(ixx))
    call parser_check_superimpose_valid_ref (ele, lat, plat%ele(ixx))
  enddo
enddo

! Go through and create the overlay, girder, and group lord elements.

call parser_add_lord (lat2, ele_num, plat, lat)

! put in field overlap

do i = 1, ele_num
  lord => lat2%ele(i)
  pele => plat%ele(lord%ixx)
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


! make matrices for entire lat

do i = 1, lat%n_ele_max
  if (lat%ele(i)%key == null_ele$) lat%ele(i)%key = -1 ! mark for deletion
enddo
call remove_eles_from_lat (lat)  ! remove all null_ele elements.

call set_flags_for_changed_attribute(lat)
call lattice_bookkeeper (lat)

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

if (allocated(lat_name))        deallocate (lat_name, lat_indexx)

call deallocate_lat_pointers (lat2)

if (present(err_flag)) err_flag = bp_com%error_flag

end subroutine
