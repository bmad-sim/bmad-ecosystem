!+
! Subroutine bmad_parser2 (in_file, ring, orbit_, make_mats6)
!
! Subroutine parse (read in) a BMAD input file.
! This subrotine assumes that ring already holds an existing lattice.
! To read in a lattice from scratch use BMAD_PARSER.
!
! With BMAD_PARSER2 you may:
!     a) Modify the attributes of elements.
!     b) Define new overlays and groups.
!     c) Superimpose new elements upon the ring.
!
! Note: If you use the superimpose feature to insert an element into the ring
!       then the index of a given element already in the ring may change.
!
! Modules needed:
!   use bmad
!
! Input:
!   in_file     -- Character(*): Input file name
!   ring        -- Ring_struct: Ring with existing layout
!   orbit_(0:)  -- Coord_struct, optional: closed orbit for when
!                           bmad_parser2 calls ring_make_mat6
!   make_mats6  -- Logical, optional: Make the 6x6 transport matrices for then
!                   Elements? Default is True.
!
! Output:
!   ring    -- Ring_struct: Ring with modifications
!-

#include "CESR_platform.inc"

subroutine bmad_parser2 (in_file, ring, orbit_, make_mats6)

  use bmad_parser_mod, except => bmad_parser2

  implicit none
    
  type (ring_struct), target :: ring, r_temp
  type (ele_struct), pointer :: ele
  type (coord_struct), optional :: orbit_(0:)
  type (parser_ring_struct) pring

  integer ix_word, i, ix, last_con, ixx, ele_num
  integer key, n_max_old
  integer, pointer :: n_max

  character(*) in_file
  character(1) delim 
  character(16) word_2, name
  character(32) word_1
  character(40) this_name
  character(280) parse_line_save

  logical, optional :: make_mats6
  logical parsing, delim_found, found, doit
  logical file_end, err_flag, finished

! init

  bmad_status%ok = .true.
  bp_com%parser_name = 'BMAD_PARSER2'
  call file_stack('push', in_file, finished)   ! open file on stack
  if (.not. bmad_status%ok) return

  n_max => ring%n_ele_max
  n_max_old = n_max

  last_con = 0

  call allocate_pring (ring, pring)

  beam_ele%name = 'BEAM'              ! fake beam element
  beam_ele%key = def_beam$            ! "definition of beam"
  beam_ele%value(particle$)   = ring%param%particle 
  beam_ele%value(energy_gev$) = 0
  beam_ele%value(n_part$)     = ring%param%n_part

  param_ele%name = 'PARAMETER'
  param_ele%key = def_parameter$
  param_ele%value(lattice_type$) = ring%param%lattice_type
  param_ele%value(beam_energy$)  = 0
  param_ele%value(taylor_order$) = ring%input_taylor_order
  param_ele%value(n_part$)       = ring%param%n_part

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
      bp_com%write_digested = .false.
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
          call get_attribute (def$, beam_ele, ring, pring, &
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
        call get_next_word (ring%lattice, ix_word, ',', &
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

      do i = 0, n_max
        if (ring%ele_(i)%name == word_1 .or. &
                          key_name(ring%ele_(i)%key) == word_1) then
          bp_com%parse_line = trim(word_2) // ' = ' // bp_com%parse_line 
          if (found) then   ! if not first time
            bp_com%parse_line = parse_line_save
          else
            parse_line_save = bp_com%parse_line
          endif
          call get_attribute (redef$, ring%ele_(i), ring, pring, &
                                               delim, delim_found, err_flag)
          if (delim_found) call warning ('BAD DELIMITER: ' // delim, ' ')
          found = .true.
        endif
      enddo

      ! If bmad_parser2 has been called from bmad_parser then check if the
      ! element was just not used in the lattice. If so then just ignore it.
      if (.not. found) then
        if (bp_com%bmad_parser_calling) then
          do i = 0, bp_com%old_ring%n_ele_max
            if (bp_com%old_ring%ele_(i)%name == word_1) then
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

      call parser_add_variable (word_1, ring)
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
      if (n_max > ubound(ring%ele_, 1)) then
        call allocate_ring_ele_(ring)
        call allocate_pring (ring, pring)
      endif

      ring%ele_(n_max)%name = word_1
      last_con = last_con + 1     ! next free slot
      ring%ele_(n_max)%ixx = last_con

      do i = 1, n_max-1
        if (ring%ele_(n_max)%name == ring%ele_(i)%name) then
          call warning ('DUPLICATE ELEMENT NAME ' // ring%ele_(n_max)%name, ' ')
          exit
        endif
      enddo

! check for valid element key name

      found = .false.

      do i = 1, n_key
        if (word_2(:ix_word) == key_name(i)(:ix_word)) then
          ring%ele_(n_max)%key = i
          call preparse_element_init (ring%ele_(n_max))
          found = .true.
          exit
        endif
      enddo

! check if element part of a element class

      if (.not. found) then
        do i = 1, n_max-1
          if (word_2 == ring%ele_(i)%name) then
            ring%ele_(n_max) = ring%ele_(i)
            ring%ele_(n_max)%name = word_1
            found = .true.
            exit
          endif
        enddo
      endif

! if none of the above then we have an error

      if (.not. found) then
        call warning ('KEY NAME NOT RECOGNIZED: ' // word_2,  &
                       'FOR ELEMENT: ' // ring%ele_(n_max)%name)
        ring%ele_(n_max)%key = 1       ! dummy value
      endif

! now get the attribute values.
! For control elements RING%ELE_()%IXX temporarily points to
! the pring structure where storage for the control lists is
                   
      call parser_set_ele_defaults (ring%ele_(n_max))

      key = ring%ele_(n_max)%key
      if (key == overlay$ .or. key == group$ .or. key == i_beam$) then
        if (delim /= '=') then
          call warning ('EXPECTING: "=" BUT GOT: ' // delim,  &
                      'FOR ELEMENT: ' // ring%ele_(n_max)%name)
        else
          if (key == overlay$) ring%ele_(n_max)%control_type = overlay_lord$
          if (key == group$)   ring%ele_(n_max)%control_type = group_lord$
          if (key == i_beam$)  ring%ele_(n_max)%control_type = i_beam_lord$
          call get_overlay_group_names(ring%ele_(n_max), ring,  &
                                              pring, delim, delim_found)
        endif
        if (key /= i_beam$ .and. .not. delim_found) then
          call warning ('NO CONTROL ATTRIBUTE GIVEN AFTER CLOSING "}"',  &
                        'FOR ELEMENT: ' // ring%ele_(n_max)%name)
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
                        'FOR ELEMENT: ' // ring%ele_(n_max)%name)
          n_max = n_max - 1
          cycle parsing_loop
        else
          call get_attribute (def$, ring%ele_(n_max), &
                                  ring, pring, delim, delim_found, err_flag)
          if (err_flag) then
            n_max = n_max - 1
            cycle parsing_loop
          endif
        endif
      enddo

! Element must be a group, overlay, or superimpose element

      if (key /= overlay$ .and. key /= group$ .and. &
              ring%ele_(n_max)%control_type /= super_lord$) then
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

  ring%param%particle    = nint(beam_ele%value(particle$))
  ring%param%lattice_type = nint(param_ele%value(lattice_type$))
  ring%input_taylor_order = nint(param_ele%value(taylor_order$))

  if (beam_ele%value(energy_gev$) /= 0) then
    ring%ele_(0)%value(beam_energy$) = beam_ele%value(energy_gev$) * 1e9
  elseif (param_ele%value(beam_energy$) /= 0) then
    ring%ele_(0)%value(beam_energy$) = param_ele%value(beam_energy$)
  endif

  if (ring%param%n_part /= param_ele%value(n_part$)) then
    ring%param%n_part = param_ele%value(n_part$)
  else
    ring%param%n_part = beam_ele%value(n_part$)
  endif

! Transfer the new elements to a safe_place

  ele_num = n_max - n_max_old
  allocate (r_temp%ele_(1:ele_num))
  r_temp%ele_(1:ele_num) = ring%ele_(n_max_old+1:n_max)
  n_max = n_max_old

! Do bookkeeping for settable dependent variables.

  do i = 1, ele_num
    ele => r_temp%ele_(i)
    call settable_dep_var_bookkeeping (ele)
  enddo

! Put in the new elements...
! First put in superimpose elements

  do i = 1, ele_num
    ele => r_temp%ele_(i)
    if (ele%control_type /= super_lord$) cycle

    select case (ele%key)
    case (wiggler$)
      if (ele%sub_key == periodic_type$) then
        if (ele%value(kz$) == 0 .and. ele%value(l$) /= 0) then
          ele%value(kz$) = pi * ele%value(n_pole$) / ele%value(l$)
        endif
      endif
    end select

    ixx = ele%ixx
    call add_all_superimpose (ring, ele, pring%ele(ixx))
  enddo

! Go through and create the overlay, i_beam, and group lord elements.

  call parser_add_lord (r_temp, ele_num, pring, ring)

! make matrices for entire ring

  call lattice_bookkeeper (ring)
  call compute_reference_energy (ring, .true.)
  doit = .true.
  if (present(make_mats6)) doit = make_mats6
  if (doit) call ring_make_mat6(ring, -1, orbit_)  ! make transport matrices
  call s_calc (ring)                       ! calc loginitudinal distances
  call ring_geometry (ring)                ! ring layout

!-----------------------------------------------------------------------------
! error check

  if (bp_com%error_flag .and. bmad_status%exit_on_error) then
    print *, 'BMAD_PARSER2 FINISHED. EXITING ON ERRORS'
    stop
  endif

  call check_ring_controls (ring, .true.)

  do i = lbound(pring%ele, 1) , ubound(pring%ele, 1)
    if (associated (pring%ele(i)%name_)) then
      deallocate(pring%ele(i)%name_)
      deallocate(pring%ele(i)%attrib_name_)
      deallocate(pring%ele(i)%coef_)
    endif
  enddo

  if (associated (pring%ele))        deallocate (pring%ele)

  do i = 1, size(r_temp%ele_)
    call deallocate_ele_pointers (r_temp%ele_(i))
  enddo

end subroutine
