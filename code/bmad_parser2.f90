!+
! Subroutine bmad_parser2 (in_file, ring, orbit_)
!
! Subroutine to parse (read in) a BMAD input file to modify an existing lattice.
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
!   in_file -- Character*(*): Input file name
!   ring    -- Ring_struct: Ring with existing layout
!   orbit_(0:n_ele_maxx) -- [Optional] Coord_struct: closed orbit for when
!                           bmad_parser2 calls ring_make_mat6
!
! Output:
!   ring    -- Ring_struct: Ring with modifications
!-

!$Id$
!$Log$
!Revision 1.5  2002/06/13 14:54:22  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:32:11  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:36  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:48  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine bmad_parser2 (in_file, ring, orbit_)

  use bmad_parser_mod

  implicit none
    
  type (ring_struct), target :: ring, r_temp
  type (ele_struct)  beam_ele, ele
  type (coord_struct), optional :: orbit_(0:n_ele_maxx)
  type (parser_ring_struct) pring

  integer ix_word, ick, i, j, k, ix, ixe, ix_lord
  integer jmax, i_key, last_con, ic, ixx, ele_num
  integer key, ix_super, ivar, n_max_old
  integer, pointer :: n_max

  character*(*) in_file
  character*16 word_2, name, a_name
  character*16 name1, name2
  character delim*1, word_1*32, call_file*200

  logical parsing, delim_found, found, matched_delim
  logical file_end, match_found, err_flag, finished

! init

  pcom%parser_name = 'BMAD_PARSER2'
  pcom%n_files = 0
  pcom%error_flag = .false.                  ! set to true on an error
  call file_stack('init', in_file, finished)    ! open file on stack
  call file_stack('push', in_file, finished)    ! open file on stack
  if (.not. bmad_status%ok) return
  call load_parse_line ('init', 0, file_end) ! initialize subroutine
  call init_bmad_parser_common

  n_max => ring%n_ele_max
  n_max_old = n_max

  last_con = 0

  pring%ele(:)%ref_name = blank
  pring%ele(:)%ref_pt  = center$
  pring%ele(:)%ele_pt  = center$
  pring%ele(:)%s       = 0
  pring%ele(:)%common_lord = .false.
  ring%ele_(:)%ixx = 0

  beam_ele%name = 'BEAM'              ! fake beam element
  beam_ele%key = def_beam$            ! "definition of beam"
  beam_ele%value(particle$) = ring%param%particle  ! default
  beam_ele%value(energy$)   = ring%param%energy
  beam_ele%value(n_part$)   = ring%param%n_part

!-----------------------------------------------------------
! main parsing loop

  parsing_loop: do

! get a line from the input file and parse out the first word

    call load_parse_line ('normal', 1, file_end)  ! load an input line
    call get_next_word (word_1, ix_word, ':(,)=', delim, delim_found, .true.)
    if (file_end) then
      word_1 = 'END_FILE'
      ix_word = 8
    else
      call verify_valid_name(word_1, ix_word)
    endif

! CALL command

    if (word_1(:ix_word) == 'CALL') then

      if (delim /= ',')  &
              call warning ('"CALL" NOT FOLLOWED BY COMMA', ' ')
      call get_next_word(call_file, ix_word, ':=,', delim, delim_found, .true.)
      if (ix_word == 0) then
        call warning ('NOTHING AFTER "CALL"', ' ')
      elseif (index('FILENAME', call_file(:ix_word)) /= 1) then
        call warning ('INVALID "CALL" COMMAND', ' ')
      elseif (delim /= '=') then
        call warning ('NO "=" AFTER "FILENAME"', ' ')
      else
        call get_next_word(call_file, ix_word, '=,', &
                                       delim, delim_found, .false.)
        if (ix_word == 0) then
          call warning ('NO FILE NAME SPECIFIED', ' ')
        else
          call file_stack ('push', call_file, finished)
          if (.not. bmad_status%ok) return
        endif
      endif
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
      if ((delim /= ':' .or. pcom%parse_line(1:1) /= '=') &
                                       .and. (delim /= '=')) then
        call warning ('"LATTICE" NOT FOLLOWED BY ":="', ' ')
      else
        if (delim == ':') pcom%parse_line = pcom%parse_line(2:)  ! trim off '='
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

! variable definition or element redef
! Note: "var := num" is old-style variable definition syntax.

    matched_delim = .false.
    if (delim == ':' .and. pcom%parse_line(1:1) == '=') then  ! old style
      matched_delim = .true.
      pcom%parse_line = pcom%parse_line(2:)      ! trim off "="
      ix = index(word_1, '[')
    elseif (delim == '=') then
      matched_delim = .true.
      ix = index(word_1, '[')
    endif

! if an element attribute redef...

    found = .false.
    if (matched_delim .and. ix /= 0) then
      name = word_1(:ix-1)
      do i = 1, n_max
        if (ring%ele_(i)%name == name) then
          name = word_1(ix+1:)    ! name of attribute
          ix = index(name, ']')
          name = name(:ix-1)
          pcom%parse_line = name // ' = ' // pcom%parse_line 
          call get_attribute (redef$, ring%ele_(i), ring, pring, &
                                               delim, delim_found, err_flag)
          if (delim_found) call warning ('BAD DELIMITER: ' // delim, ' ')
          found = .true.
        endif
      enddo

      if (.not. found) call warning ('ELEMENT NOT FOUND: ' // name, ' ')

      cycle parsing_loop

! else must be a variable

    elseif (matched_delim) then

      found = .false.
      do i = 1, pcom%ivar_tot-1
        if (word_1 == var_(i)%name) then
          ivar = i
          found = .true.
        endif
      enddo

      if (.not. found) then
        pcom%ivar_tot = pcom%ivar_tot + 1
        ivar = pcom%ivar_tot
        if (pcom%ivar_tot > ivar_maxx) then
          type *, 'ERROR IN BMAD_PARSER2: NEED TO INCREASE IVAR_MAXX!'
          call err_exit
        endif
      endif

      var_(ivar)%name = word_1
      call evaluate_value (var_(ivar)%name, var_(ivar)%value, &
                                    ring, delim, delim_found, err_flag)
      if (delim /= ' ' .and. .not. err_flag) call warning  &
            ('EXTRA CHARACTERS ON RHS: ' // pcom%parse_line,  &
             'FOR VARIABLE: ' // var_(ivar)%name)
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
      if (n_max > n_ele_maxx) then
        type *, 'ERROR IN BMAD_PARSER2: NEED TO INCREASE ELEMENT ARRAY!'
        call err_exit
      endif

      call init_ele (ring%ele_(n_max))    ! init element
      ring%ele_(n_max)%name = word_1
      last_con = last_con + 1     ! next free slot
      ring%ele_(n_max)%ixx = last_con

      do i = 1, n_max-1
        if (ring%ele_(n_max)%name == ring%ele_(i)%name)  &
                         call warning ('DUPLICATE ELEMENT NAME ' //  &
                         ring%ele_(n_max)%name, ' ')
      enddo

! check if element part of a element class

      do i = 1, n_max-1
        if (word_2 == ring%ele_(i)%name) then
          i_key = ring%ele_(i)%key
          ring%ele_(n_max)%key =  i_key
          ring%ele_(n_max)%type = ring%ele_(i)%type
          ring%ele_(n_max)%value = ring%ele_(i)%value
          exit
        endif
      enddo

! check for valid element key name
! if none of the above then we have an error

      if (i == n_max) then
        do i = 1, n_key-1
          if (word_2(:ix_word) == key_name(i)(:ix_word)) then
            ring%ele_(n_max)%key = i
            exit
          endif
        enddo
        if (i == n_key) then
          call warning ('KEY NAME NOT RECOGNIZED: ' // word_2,  &
                       'FOR ELEMENT: ' // ring%ele_(n_max)%name)
          ring%ele_(n_max)%key = 1       ! dummy value
        endif
      endif

! now get the attribute values.
! For control elements RING.ELE_().IXX temporarily points to
! the pring structure where storage for the control lists is
                   
      key = ring%ele_(n_max)%key
      if (key == overlay$ .or. key == group$) then
        if (delim /= '=') then
          call warning ('EXPECTING: "=" BUT GOT: ' // delim,  &
                      'FOR ELEMENT: ' // ring%ele_(n_max)%name)
        else
          ring%ele_(n_max)%control_type = key
          call get_overlay_group_names(ring%ele_(n_max), ring,  &
                                              pring, delim, delim_found)
        endif
        if (.not. delim_found) then
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

  ring%param%particle = nint(beam_ele%value(particle$))
  ring%param%energy   = beam_ele%value(energy$)
  ring%param%n_part    = beam_ele%value(n_part$)

! Transfer the new elements to a safe_place

  ele_num = n_max - n_max_old
  r_temp%ele_(1:ele_num) = ring%ele_(n_max_old+1:n_max)
  n_max = n_max_old

! Put in the new elements...
! First put in superimpose elements

  do i = 1, ele_num
    if (r_temp%ele_(i)%control_type == super_lord$) then
      ixx = r_temp%ele_(i)%ixx
      call add_all_superimpose (ring, r_temp%ele_(i), pring%ele(ixx))
    endif
  enddo

! Find controlled elements for new overlay and group lords

  do ixe = 1, ele_num

    if (r_temp%ele_(ixe)%control_type == super_lord$) cycle

    ic = r_temp%ele_(ixe)%ixx
    jmax = r_temp%ele_(ixe)%n_slave

    do j = 1, jmax

      match_found = .false.

      do k = 1, ring%n_ele_max
                 
        name1 = ring%ele_(k)%name
        name2 = null_name
        ick = index(name1, '\')              !'
        if (ick /= 0) name2 = name1(:ick-1)

        if (pring%ele(ic)%name_(j) == name1 .or.  &
                pring%ele(ic)%name_(j) == name2 .and.  &
                (ring%ele_(k)%control_type /= super_slave$)) then

          if (match_found) &
              call warning ('MULTIPLE CONTROL MATCHES FOUND FOR: '  &
              //  r_temp%ele_(ixe)%name, &
              'MATCHING STRING: ' // pring%ele(ic)%name_(j))

          match_found = .true.
          pring%ele(ic)%cs_(j)%ix_slave = k

          a_name = pring%ele(ic)%attrib_name_(j)
          ix = attribute_index(ring%ele_(k), a_name)
          if (a_name /= blank) then
            pring%ele(ic)%cs_(j)%ix_attrib = ix
            if (ix < 1) call warning ('BAD ATTRIBUTE NAME: ' // a_name,  &
                     'FOR OVERLAY/GROUP ELEMENT: ' // r_temp%ele_(ixe)%name)
          else
            pring%ele(ic)%cs_(j)%ix_attrib = r_temp%ele_(ixe)%ix_value
          endif

        endif ! pring%ele(ic).name(j) == name1 or name2

      enddo   ! k = 1, ring.n_ele_ring+ring.n_ele_control

      if (.not. match_found)  &
             call warning ('NO SLAVE ELEMENT FOUND FOR: ' //  &
             r_temp%ele_(ixe)%name, 'FAILED TO FIND: ' // pring%ele(ic)%name_(j))

    enddo   ! j = 1, jmax

  enddo

! Go through the overlay elements.
! Transfer control information from PRING%ELE array to RING.CONTROL_ array.

  do i = 1, ele_num
    ele = r_temp%ele_(i)
    if (ele%control_type == overlay_lord$) then
      call new_control (ring, ix_lord)
      ring%ele_(ix_lord) = ele
      call create_overlay (ring, ix_lord, ele%ix_value, ele%n_slave, &
                                                     pring%ele(ele%ixx)%cs_)
    endif
  enddo

  if (ring%n_control_array > n_control_maxx) then
    type *, 'ERROR IN BMAD_PARSER2: NOT ENOUGH CONTROL ELEMENTS !!!'
    type *, '      YOU NEED TO INCREASE N_CONTROL_MAXX IN BMAD_STRUCT !!!'
    call err_exit
  endif

! make up groups
! Note that CREATE_GROUP needs RING.PARAM.TOTAL_LENGTH

  do i = 1, ele_num
    ele = r_temp%ele_(i)
    if (ele%control_type == group_lord$ .and. ic > 0) then
      call new_control (ring, ix_lord)
      ring%ele_(ix_lord) = ele
      call create_group (ring, ix_lord, ele%n_slave, pring%ele(ele%ixx)%cs_)
    endif
  enddo

! make matrices for entire ring

  call ring_make_mat6 (ring, -1, orbit_)
  call s_calc (ring)                       ! calc loginitudinal distances
  call ring_geometry (ring)                ! ring layout

!-------------------------------------------------------------------------
! write out if debug is on

  if (pcom%parser_debug) then

    type *
    type *, '----------------------------------------'
    type *, 'Number of Elements in the Regular Ring:', ring%n_ele_ring
    do i = 1, ring%n_ele_ring
      type *, '-------------'
      type *, 'Ele #', i
      call type_ele (ring%ele_(i), .false., 0, .false., 0, .true., ring)
    enddo

    type *
    type *, '----------------------------------------'
    type *, 'Control elements: ', ring%n_ele_max - ring%n_ele_ring
    do i = ring%n_ele_ring+1, ring%n_ele_max
      type *, '-------------'
      type *, 'Ele #', i
      call type_ele (ring%ele_(i), .false., 0, .false., 0, .true., ring)
    enddo


    type *
    type *, '----------------------------------------'
    type *, 'Ring Used: ', ring%name
    type *, 'Number of ring elements:', ring%n_ele_ring
    type *, 'List:                               Key      Length         S'
    do i = 1, ring%n_ele_ring
      type '(3x, i3, 2a, 3x, a, 2f10.2)', i, ') ', ring%ele_(i)%name,  &
        key_name(ring%ele_(i)%key), ring%ele_(i)%value(l$), ring%ele_(i)%s
    enddo

  endif

!-----------------------------------------------------------------------------
! error check

  if (pcom%error_flag .and. bmad_status%exit_on_error) then
    type *, 'BMAD_PARSER2 FINISHED. EXITING ON ERRORS'
    call exit
  endif

  call check_ring_controls (ring, .true.)

end subroutine
