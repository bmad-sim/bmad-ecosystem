!+
! Subroutine bmad_parser (in_file, ring)
!
! Subroutine to parse a BMAD input file.
!
! Because of the time it takes to parse a file BMAD_PARSER will save 
! RING in a "digested" file with the name:
!                 'DIGESTED_' // IN_FILE
! For subsequent calls to the same IN_FILE, BMAD_PARSER will just read in the
! digested file. BMAD_PARSER will always check to see that the digested file
! is up-to-date.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   in_file -- Character: Name of the input file.
!
! Output:
!   ring -- Ring_struct: Ring structure. See bmad_struct.f90 for more details.
!     %ele_(:)%mat6 -- This is computed assuming an on-axis orbit 
!     %ele_(:)%s    -- This is also computed.
!         
! Defaults:
!   ring%n_ele_symm         = 0
!   ring%param%particle     = positron$
!   ring%param%symmetry     = no_symmetry$
!   ring%param%z_decoupled  = .false.
!   ring%param%aperture_limit_on = .true.
!   ring%param%lattice_type = circular_lattice$
!
! For more info on the parser see the BMAD documentation.
! DCS 10/6/97
!-

!$Id$
!$Log$
!Revision 1.6  2002/01/08 21:44:36  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.5  2001/10/12 20:53:34  rwh24
!DCS changes and two files added
!
!Revision 1.4  2001/10/05 18:23:57  rwh24
!Bug Fixes
!
!Revision 1.3  2001/10/02 18:49:11  rwh24
!More compatibility updates; also added many explicit variable declarations.
!
!Revision 1.2  2001/09/27 18:31:48  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine bmad_parser (in_file, ring)

  use local_bmad_struct
  use local_bmad_interface
  use cesr_utils
  use nr, only: indexx
  
  implicit none

  type (ring_struct), target :: ring, in_ring
  type (seq_struct), target :: seq_(n_ele_maxx)
  type (seq_struct), pointer :: seq, seq2
  type (seq_ele_struct), pointer :: s_ele
  type (parser_ring_struct) pring
  type (seq_stack_struct) stack(20)
  type (ele_struct), pointer :: ele

  integer ix_word, i_use, i, j, k, n, ix, ixr, ixs, i_ring
  integer i_lev, i_key, ic, ix_lord
  integer iseq, ix_ring(n_ele_maxx), n_ele_ring
  integer ix_super, digested_version, key
  integer ix1, ix2, iv, n_arg
  integer ivar, ixx, j_lord, n_slave
  integer seq_indexx(n_ele_maxx), in_indexx(n_ele_maxx), r_indexx(n_ele_maxx)
  integer, pointer :: n_max

  character*(*) in_file
  character*16 word_2, name, a_name, dummy_name(20)
  character*16 name_(n_ele_maxx)
  character delim*1, word_1*40
  character*200 path, basename, full_name, digested_file, call_file
  
  real angle

  logical parsing, delim_found, matched_delim, arg_list_found
  logical file_end, found, err_flag, finished
  logical, save :: init_needed = .true.

! see if digested file is open and current.
! if so read in and return

  pcom%parser_name = 'BMAD_PARSER'

!  ix = index(in_file, ']')
!  if (ix /= 0) then
!    digested_file = in_file(:ix) // 'digested_' // in_file(ix+1:)
!  else
!    digested_file = 'digested_' // in_file
!  endif

  ix = SplitFileName(in_file, path, basename)
  digested_file = in_file(:ix) // 'digested_' // in_file(ix+1:)

  call read_digested_bmad_file (digested_file, ring, digested_version)
  if (bmad_status%ok) return

! init

  if (init_needed) then
    do i = lbound(pring%ele, 1) , ubound(pring%ele, 1)
      nullify (pring%ele(i)%name_)
    enddo
    do i = 1, size(seq_(:))
      nullify (seq_(i)%arg, seq_(i)%ele)
    enddo
    init_needed = .false.
  endif

! here if not OK bmad_status. So we have to do everything from scratch...
! init variables.

  bmad_status%ok = .true.
  if (bmad_status%type_out) print *, 'BMAD_PARSER: CREATING NEW DIGESTED FILE...'

  pcom%n_files = 0
  pcom%error_flag = .false.                 ! set to true on an error
  call file_stack('init', in_file, finished)   ! init stack
  call file_stack('push', in_file, finished)   ! open file on stack
  if (.not. bmad_status%ok) return
  call load_parse_line ('init', 0, file_end) ! initialize subroutine
  pcom%iseq_tot = 0                     ! number of sequences encountered
  pcom%ivar_tot = 0                     ! number of variables encountered
  call init_bmad_parser_common
  ring%name = ' '

  n_max => in_ring%n_ele_max
  n_max = 0                         ! number of elements encountered

  pring%ele(:)%ref_name = blank
  pring%ele(:)%ref_pt  = center$
  pring%ele(:)%ele_pt  = center$
  pring%ele(:)%s       = 0
  pring%ele(:)%common_lord = .false.

  in_ring%ele_(0)%name = 'BEAM'                 ! fake beam element
  in_ring%ele_(0)%key = def_beam$               ! "definition of beam"
  in_ring%ele_(0)%value(particle$) = positron$  ! default

! %ixx is used as a pointer from the in_ring%ele_ array to the pring%ele array

  forall (i = 0:n_ele_maxx) in_ring%ele_(i)%ixx = i

  ring%n_control_array = 0

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

! USE command...

    if (word_1(:ix_word) == 'USE') then
      if (delim /= ',') call warning ('"USE" NOT FOLLOWED BY COMMA')
      call get_next_word(word_2, ix_word, ':(=,)', delim, delim_found, .true.)
      if (ix_word == 0) call error_exit  &
                                ('NO BEAM LINE SPECIFIED WITH "USE"', ' ')
      call verify_valid_name(word_2, ix_word)
      ring%name = word_2
      cycle parsing_loop
    endif

! TITLE command

    if (word_1(:ix_word) == 'TITLE') then
      if (delim_found) call warning ('EXTRA STUFF FOUND AFTER "TITLE"')
      read (pcom%f_unit, '(a)') ring%title
      pcom%i_line = pcom%i_line + 1
      cycle parsing_loop
    endif

! CALL command

    if (word_1(:ix_word) == 'CALL') then
      if (delim /= ',')  call warning ('"CALL" NOT FOLLOWED BY COMMA')
      call get_next_word(call_file, ix_word, ':=,', delim, delim_found, .true.)
      if (ix_word == 0) then
        call warning ('NOTHING AFTER "CALL"')
      elseif (index('FILENAME', call_file(:ix_word)) /= 1) then
        call warning ('INVALID "CALL" COMMAND')
      elseif (delim /= '=') then
        call warning ('NO "=" AFTER "FILENAME"')
      else
        call get_next_word(call_file, ix_word, '=,', &
                                       delim, delim_found, .false.)
        if (ix_word == 0) then
          call warning ('NO FILE NAME SPECIFIED')
        else
          call fullfilename(call_file, call_file)
          call file_stack ('push', call_file, finished)
          if (.not. bmad_status%ok) return
        endif
      endif     
      cycle parsing_loop
    endif

! BEAM command

    if (word_1(:ix_word) == 'BEAM') then
      if (delim /= ',') call warning ('"BEAM" NOT FOLLOWED BY COMMA')
      parsing = .true.
      do while (parsing)
        if (.not. delim_found) then
          parsing = .false.
        elseif (delim /= ',') then
          call warning ('EXPECTING: "," BUT GOT: ' // delim,  &
                                                     'FOR "BEAM" COMMAND')
          parsing = .false.
        else
          call get_attribute (def$, in_ring%ele_(0), &
                                 in_ring, pring, delim, delim_found, err_flag)
        endif
      enddo
      cycle parsing_loop
    endif
                   
! LATTICE command

    if (word_1(:ix_word) == 'LATTICE') then
      if ((delim /= ':' .or. pcom%parse_line(1:1) /= '=') &
                                       .and. (delim /= '=')) then
        call warning ('"LATTICE" NOT FOLLOWED BY ":="')
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
      if (finished) exit ! break loop
      cycle parsing_loop
    endif

! variable definition or element redef.
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

! if an element attribute redef

    if (matched_delim .and. ix /= 0) then
      name = word_1(:ix-1)  
      do i = 1, n_max
        if (in_ring%ele_(i)%name == name) then
          name = word_1(ix+1:)    ! name of attribute
          ix = index(name, ']')
          name = name(:ix-1)
          pcom%parse_line = name // ' = ' // pcom%parse_line 
          call get_attribute (redef$, in_ring%ele_(i), in_ring, pring, &
                                          delim, delim_found, err_flag)
          if (delim_found) call warning ('BAD DELIMITER: ' // delim)
          cycle parsing_loop
        endif
      enddo

      call warning ('ELEMENT NOT FOUND: ' // name)
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
          print *, 'ERROR IN BMAD_PARSER: NEED TO INCREASE IVAR_MAXX!'
          call err_exit
        endif
      endif

      var_(ivar)%name = word_1

      call evaluate_value (var_(ivar), in_ring, delim, delim_found, err_flag)
      if (delim /= ' ' .and. .not. err_flag) call warning  &
                    ('EXTRA CHARACTERS ON RHS: ' // pcom%parse_line,  &
                     'FOR VARIABLE: ' // var_(ivar)%name)
      cycle parsing_loop

    endif

! if a "(" delimitor then we are looking at a replacement line.

    if (delim == '(') then
      n_arg = 0
      do
        call get_next_word (word_2, ix_word, '(): =,', &
                                                 delim, delim_found, .true.)
        if (ix_word == 0 .or. delim == '( :=') then
          call warning ('BAD ARGUMENT LIST FOR: ', word_1)
          cycle parsing_loop
        endif
        n_arg = n_arg + 1
        dummy_name(n_arg) = word_2
        if (delim == ')') then
          call get_next_word (word_2, ix_word, '(): =,', &
                                                 delim, delim_found, .true.)
          if (word_2 /= ' ') call warning &
                  ('":" NOT FOUND AFTER REPLACEMENT LINE ARGUMENT LIST. ' // &
                  'FOUND: ' // word_2, 'FOR LINE: ' // word_1)
          exit
        endif
      enddo
      allocate (seq_(pcom%iseq_tot+1)%arg(n_arg))
      seq_(pcom%iseq_tot+1)%arg(:)%dummy_name = dummy_name(1:n_arg)
      arg_list_found = .true.
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

    call verify_valid_name(word_2, ix_word)

! arg lists are only used with lines

    if (word_2(:ix_word) /= 'LINE' .and. arg_list_found) then
      call warning ('ARGUMENT LISTS "(...)" ARE ONLY USED WITH LINES', &
                                                        'FOR: ' // word_1)
      cycle parsing_loop
    endif

! if line or list

    if (word_2(:ix_word) == 'LINE' .or. word_2(:ix_word) == 'LIST') then
      pcom%iseq_tot = pcom%iseq_tot + 1
      if (pcom%iseq_tot > n_ele_maxx-1) then
        print *, 'ERROR IN BMAD_PARSER: NEED TO INCREASE LINE ARRAY!'
        call err_exit
      endif
      seq_(pcom%iseq_tot)%name = word_1

      if (delim /= '=') call warning ('EXPECTING: "=" BUT GOT: ' // delim)
      if (word_2(:ix_word) == 'LINE') then
        seq_(pcom%iseq_tot)%type = line$
        if (arg_list_found) seq_(pcom%iseq_tot)%type = replacement_line$
      else
        seq_(pcom%iseq_tot)%type = list$
      endif
      call seq_expand1 (seq_, pcom%iseq_tot, .true.)

! if not line or list then must be an element

    else

      n_max = n_max + 1
      if (n_max > n_ele_maxx) then
        print *, 'ERROR IN BMAD_PARSER: NEED TO INCREASE ELEMENT ARRAY!'
        call err_exit
      endif

      call init_ele (in_ring%ele_(n_max))
      in_ring%ele_(n_max)%name = word_1

! check for valid element key name

      do i = 1, n_key-1
        if (word_2(:ix_word) == key_name(i)(:ix_word)) then
          in_ring%ele_(n_max)%key = i
          exit
        endif
      enddo

! check if element part of a element class
! if none of the above then we have an error

      if (i == n_key) then
        do i = 1, n_max-1
          if (word_2 == in_ring%ele_(i)%name) then
            i_key = in_ring%ele_(i)%key
            in_ring%ele_(n_max)%key = i_key
            in_ring%ele_(n_max)%type = in_ring%ele_(i)%type
            in_ring%ele_(n_max)%value = in_ring%ele_(i)%value
            exit
          endif
        enddo

        if (i == n_max) then
          call warning ('KEY NAME NOT RECOGNIZED: ' // word_2,  &
                       'FOR ELEMENT: ' // in_ring%ele_(n_max)%name)
          cycle parsing_loop
        endif

      endif

! now get the attribute values.
! For control elements RING.ELE_()%IXX temporarily points to
! the PRING%ELE() array where storage for the control lists is

      key = in_ring%ele_(n_max)%key
      if (key == overlay$ .or. key == group$ .or. key == coil$) then
        if (delim /= '=') then
          call warning ('EXPECTING: "=" BUT GOT: ' // delim,  &
                      'FOR ELEMENT: ' // in_ring%ele_(n_max)%name)
          cycle parsing_loop        
        endif

        if (key == coil$) then
          in_ring%ele_(n_max)%control_type = container_slave$
        else
          in_ring%ele_(n_max)%control_type = key
        endif
        call get_overlay_group_names(in_ring%ele_(n_max), in_ring, &
                                                    pring, delim, delim_found)

        if (key /= coil$ .and. .not. delim_found) then
          call warning ('NO CONTROL ATTRIBUTE GIVEN AFTER CLOSING "}"',  &
                        'FOR ELEMENT: ' // in_ring%ele_(n_max)%name)
          cycle parsing_loop
        endif

      endif

      parsing = .true.
      do while (parsing)
        if (.not. delim_found) then          ! if nothing more
          parsing = .false.
        elseif (delim /= ',') then
          call warning ('EXPECTING: "," BUT GOT: ' // delim,  &
                        'FOR ELEMENT: ' // in_ring%ele_(n_max)%name)
          cycle parsing_loop
        else
          call get_attribute (def$, in_ring%ele_(n_max), &
                                  in_ring, pring, delim, delim_found, err_flag)
          if (err_flag) cycle parsing_loop
        endif
      enddo

    endif

  enddo parsing_loop       ! main parsing loop

!---------------------------------------------------------------
! we now have read in everything. 

! sort elements and lists and check for duplicates

  call indexx (seq_(1:pcom%iseq_tot)%name, seq_indexx(1:pcom%iseq_tot))
  call indexx (in_ring%ele_(1:n_max)%name, in_indexx(1:n_max))

  do i = 1, pcom%iseq_tot-1
    ix1 = seq_indexx(i)
    ix2 = seq_indexx(i+1)
    if (seq_(ix1)%name == seq_(ix2)%name) call warning  &
                      ('DUPLICATE LINE NAME ' // seq_(ix1)%name)
  enddo

  do i = 1, n_max-1
    ix1 = in_indexx(i)
    ix2 = in_indexx(i+1)
    if (in_ring%ele_(ix1)%name == in_ring%ele_(ix2)%name) call warning &
                    ('DUPLICATE ELEMENT NAME ' // in_ring%ele_(ix1)%name)
  enddo

  i = 1; j = 1
  do
    if (i > pcom%iseq_tot) exit
    if (j > n_max) exit
    ix1 = seq_indexx(i)
    ix2 = in_indexx(j)
    if (seq_(ix1)%name == in_ring%ele_(ix2)%name) call warning  &
          ('LINE AND ELEMENT HAVE THE SAME NAME: ' // seq_(i)%name)
    if (seq_(ix1)%name < in_ring%ele_(ix2)%name) then
      i = i + 1
    else
      j = j + 1
    endif
  enddo

! find line corresponding to the "use" statement.

  if (ring%name == blank) call err_exit &
            ('NO "USE" COMMAND FOUND.', 'I DO NOT KNOW WHAT LINE TO USE!')

  call find_indexx (ring%name, seq_(:)%name, seq_indexx, pcom%iseq_tot, i_use)

  if (i_use == 0) call err_exit &
            ('CANNOT FIND DEFINITION FOR "USE" LINE: ' // ring%name, ' ')

  if (seq_(i_use)%type /= line$) call error_exit  &
                              ('NAME AFTER "USE" IS NOT A LINE!', ' ')

! Now to expand the lines and lists to find the elements to use
! first go through the lines and lists and index everything

  do k = 1, pcom%iseq_tot
    do i = 1, size(seq_(k)%ele(:))

      s_ele => seq_(k)%ele(i)

      ix = index(s_ele%name, '\')
      if (ix /= 0) then
        name = s_ele%name(:ix-1)
      else
        name = s_ele%name
      endif
  
      if (s_ele%ix_arg > 0) cycle  ! dummy arg

      call find_indexx (name, in_ring%ele_(1:)%name, in_indexx, n_max, j)

      if (j == 0) then  ! if not an element it must be a sequence
        call find_indexx (name, seq_(:)%name, seq_indexx, pcom%iseq_tot, j)
        if (j == 0) then  ! if not a sequence then I don't know what it is
          call warning  &
              ('NOT A DEFINED ELEMENT, LINE, OR LIST: ' // s_ele%name)
          s_ele%ix_array = 1
          s_ele%type = element$
        else
          s_ele%ix_array = j
          s_ele%type = seq_(j)%type
          if (seq_(k)%type == list$) call warning ('LIST: ' // seq_(k)%name,  &
                         'CONTAINS A NON-ELEMENT: ' // s_ele%name)
        endif
        if (s_ele%type == list$ .and. s_ele%reflect) call warning ( &
                          'A REFLECTION WITH A LIST IS NOT ALLOWED IN: '  &
                          // seq_(k)%name, 'FOR LIST: ' // s_ele%name)

      else    ! if an element...
        s_ele%ix_array = j
        s_ele%type = element$
      endif

    enddo
  enddo

! to expand the ring we use a stack for nested sublines
! IX_RING is the expanded array of elements in the ring.
! init stack

  i_lev = 1                          ! level on the stack
  stack(1)%ix_seq  = i_use           ! which sequence to use for the ring
  stack(1)%ix_ele  =  1              ! we start at the beginning
  stack(1)%reflect = +1              ! and move forward

  seq => seq_(i_use)

  n_ele_ring = 0
             
! expand line

  parsing = .true.
  line_expansion: do while (parsing)

    s_ele => seq%ele(stack(i_lev)%ix_ele)

! if an element

    if (s_ele%type == element$) then
      call pushit (ix_ring, n_ele_ring, s_ele%ix_array)
      name_(n_ele_ring) = s_ele%name
      call increment_pointer (stack(i_lev)%ix_ele, stack(i_lev)%reflect)

! if a list

    elseif (s_ele%type == list$) then
      seq2 => seq_(s_ele%ix_array)
      j = seq2%ix
      call pushit (ix_ring, n_ele_ring, seq2%ele(j)%ix_array)
      name_(n_ele_ring) = s_ele%name
      seq2%ix = seq2%ix + 1
      if (seq2%ix > size(seq2%ele(:))) seq2%ix = 1
      call increment_pointer (stack(i_lev)%ix_ele, stack(i_lev)%reflect)

! if a line:
!     a) move pointer on current level past line element
!     b) go to the next higher level
!     c) initialize pointers for the higher level to use the line

    elseif (s_ele%type == line$) then
      call increment_pointer (stack(i_lev)%ix_ele, stack(i_lev)%reflect)
      i_lev = i_lev + 1
      seq => seq_(s_ele%ix_array)
      stack(i_lev)%ix_seq = s_ele%ix_array
      stack(i_lev)%reflect = stack(i_lev-1)%reflect
      if (s_ele%reflect) stack(i_lev)%reflect = -stack(i_lev)%reflect
      if (stack(i_lev)%reflect == 1) then
        stack(i_lev)%ix_ele = 1
      else
        stack(i_lev)%ix_ele = size(seq%ele(:))
      endif

    else
      call warning ('INTERNAL SEQUENCE ERROR!')

    endif

! if we have got to the end of the current line then pop the stack back to
! the next lower level.
! Also check if we have gotten to level 0 which says that we are done

    do
      if (stack(i_lev)%ix_ele < 1 .or. &
                        stack(i_lev)%ix_ele > size(seq%ele(:))) then
        i_lev = i_lev - 1
        if (i_lev == 0) exit line_expansion
        seq => seq_(stack(i_lev)%ix_seq)
      else
        exit
      endif
    enddo

  enddo line_expansion

!---------------------------------------------------------------
! we now have expanded the lines and lists.
! now to put the elements in RING in the correct order.
! superimpose, overlays, and groups are handled later.
! first load beam parameters.

  ring%version            = bmad_inc_version$
  ring%param%particle     = nint(in_ring%ele_(0)%value(particle$))
  ring%param%energy       = in_ring%ele_(0)%value(energy$)
  ring%param%n_part       = in_ring%ele_(0)%value(n_part$)
  ring%param%symmetry     = no_symmetry$
  ring%param%lattice_type = circular_lattice$
  ring%param%z_decoupled  = .false.                ! default
  ring%n_ele_ring         = n_ele_ring
  ring%n_ele_use          = n_ele_ring
  ring%n_ele_max          = n_ele_ring
  ring%param%aperture_limit_on = .true.
  ring%n_ele_symm         = 0                     ! no symmetry point
  inquire (file = in_file, name = full_name)      !! full input file_name
  ix = index(full_name, ';')
  if (ix /= 0) full_name = full_name(:ix-1)
  ring%input_file_name    = full_name             ! save input file
  call set_symmetry (ring%param%symmetry, ring)   ! set ring.n_ele_use
  ring%ele_(0)%name = 'BEGINNING'                 ! Just a name

  do i = 1, pcom%ivar_tot
    if (var_(i)%name == 'SYMMETRY') ring%param%symmetry = var_(i)%value
    if (var_(i)%name == 'LATTICE_TYPE')  &
                              ring%param%lattice_type = var_(i)%value
  enddo

! Use name as given in sequence lists since elements can
! have different names from the defining elements. 
! Eg: B01\H2 gets its definition from B01.
! Also: convert rbends to sbends
             
  if (n_ele_ring > n_ele_maxx) then
    print *, 'ERROR IN BMAD_PARSER: NUMBER OF ELEMENTS EXCEEDS ELEMENT ARRAY'
    print *, '       SIZE THE RING STRUCTURE:', n_ele_ring, n_ele_maxx
    call err_exit
  endif

  do i = 1, n_ele_ring

    ele => ring%ele_(i)

    ele = in_ring%ele_(ix_ring(i))
    ele%name = name_(i)

    if (ele%key == sbend$ .or. ele%key == rbend$) then

      if (ele%value(rho$) /= 0 .and. ele%value(rho_design$) == 0) &
                                ele%value(rho_design$) = ele%value(rho$)
      if (ele%value(rho_design$) /= 0 .and. ele%value(rho$) == 0) &
                                ele%value(rho$) = ele%value(rho_design$)

      if (ele%key == rbend$) then
        ele%value(l_chord$) = ele%value(l$)
        angle = ele%value(angle$) 
        if (ele%value(angle$) /= 0) then
          ele%value(l$) = ele%value(l_chord$) * angle / (2 * sin(angle/2))
        elseif (ele%value(rho_design$) /= 0) then
          ele%value(l$) = angle * ele%value(rho_design$)
        endif
        ele%value(e1$) = ele%value(e1$) + angle / 2
        ele%value(e2$) = ele%value(e2$) + angle / 2
        ele%key = sbend$
      endif

      if (ele%value(rho$) == 0 .and. ele%value(angle$) /= 0) then
        ele%value(rho$)        = ele%value(l$) / ele%value(angle$) 
        ele%value(rho_design$) = ele%value(l$) / ele%value(angle$) 
      elseif (ele%value(angle$) == 0 .and. ele%value(rho_design$) /= 0) then
        ele%value(angle$) = ele%value(l$) / ele%value(rho$)
        ele%value(rho_design$) = ele%value(l$) / ele%value(angle$) 
      elseif (ele%value(rho_design$) /= 0 .and. ele%value(angle$) /= 0)  then
        call warning ('BOTH RHO AND ANGLE SPECIFIED FOR BEND: ' // ele%name)
      endif

      if (ele%value(hgapx$) == 0) ele%value(hgapx$) = ele%value(hgap$)
      if (ele%value(fintx$) == 0) ele%value(fintx$) = ele%value(fint$)

    endif


    if (ele%control_type == container_slave$) then
      ele%n_lord = ele%n_slave
      ele%n_slave = 0
      ele%ic1_lord = ring%n_ic_array + 1
      ele%ic2_lord = ring%n_ic_array + ele%n_lord
      ring%n_ic_array = ele%ic2_lord
      do k = 1, ele%n_lord
        name = pring%ele(ix_ring(i))%name_(k)
        call find_indexx (name, in_ring%ele_(1:)%name, in_indexx, n_max, j)
        if (j == 0) exit
        ele%value(l$) = ele%value(l$) + in_ring%ele_(j)%value(l$)
        in_ring%ele_(j)%n_slave = in_ring%ele_(j)%n_slave + 1
      enddo
    endif

  enddo

!-------------------------------------------------------------------------
! Go through the IN_RING elements and put in the superpositions, groups, etc.

  call s_calc (ring)              ! calc longitudinal distances

! superpositions first

  do i = 1, n_max
    if (in_ring%ele_(i)%control_type /= super_lord$) cycle
    call add_all_superimpose (ring, in_ring%ele_(i), pring%ele(i))
  enddo

! Now put in the overlay_lord elements

  do i = 1, n_max
    if (in_ring%ele_(i)%control_type == overlay_lord$) then
      ixr = ring%n_ele_max
      call indexx (ring%ele_(1:ixr)%name, r_indexx(1:ixr))

      ring%n_ele_max = ring%n_ele_max + 1
      ixs = ring%n_ele_max
      ring%ele_(ixs) = in_ring%ele_(i)
      n_slave = in_ring%ele_(i)%n_slave

      do j = 1, n_slave
        call find_indexx (pring%ele(i)%name_(j), ring%ele_(1:)%name, &
                                                          r_indexx, ixr, k)
        if (k == 0) then
          call warning ('CANNOT FIND OVERLAY_SLAVE FOR: ' // &
                  ring%ele_(ixs)%name, 'CANNOT FIND: ' // pring%ele(i)%name_(j))
          cycle
        endif
        pring%ele(i)%cs_(j)%ix_slave = k
        a_name = pring%ele(i)%attrib_name_(j)
        if (a_name == blank) then
          pring%ele(i)%cs_(j)%ix_attrib = ring%ele_(ixs)%ix_value
        else
          ix = attribute_index(ring%ele_(k), a_name)
          pring%ele(i)%cs_(j)%ix_attrib = ix
          if (ix < 1) call warning('BAD ATTRIBUTE NAME: ' // a_name,  &
                       'IN OVERLAY ELEMENT: ' // ring%ele_(ixs)%name)
        endif
      enddo

      iv = ring%ele_(ixs)%ix_value
      call create_overlay (ring, ixs, iv, n_slave, pring%ele(i)%cs_)
      cycle

    endif  

! put in groups
! Note that CREATE_GROUP needs RING%PARAM%TOTAL_LENGTH

    if (in_ring%ele_(i)%control_type == group_lord$) then
      ixr = ring%n_ele_max
      call indexx (ring%ele_(1:ixr)%name, r_indexx(1:ixr))

      ring%n_ele_max = ring%n_ele_max + 1
      ixs = ring%n_ele_max
      ring%ele_(ixs) = in_ring%ele_(i)
      n_slave = in_ring%ele_(i)%n_slave              

      do j = 1, n_slave
        call find_indexx (pring%ele(i)%name_(j), ring%ele_(1:)%name, &
                                                        r_indexx, ixr, k)
        if (k == 0) then
          call warning ('CANNOT FIND GROUP_SLAVE FOR: ' // &
                ring%ele_(ixs)%name, 'CANNOT FIND: ' // pring%ele(i)%name_(j))
          cycle
        endif
        pring%ele(i)%cs_(j)%ix_slave = k
        a_name = pring%ele(i)%attrib_name_(j)
        if (a_name == blank) then
          pring%ele(i)%cs_(j)%ix_attrib = ring%ele_(ixs)%ix_value
        else
          ix = attribute_index(ring%ele_(k), a_name)
          pring%ele(i)%cs_(j)%ix_attrib = ix
          if (ix < 1) call warning('BAD ATTRIBUTE NAME: ' // a_name,  &
                       'IN OVERLAY ELEMENT: ' // ring%ele_(ixs)%name)
        endif
      enddo

      call create_group (ring, ixs, n_slave, pring%ele(i)%cs_)

    endif

  enddo

! put the components for the containers in ring.

  do i = 1, ring%n_ele_ring
    if (ring%ele_(i)%control_type == container_slave$) then
      ixx = ring%ele_(i)%ixx
      n = ring%n_ele_max
      call indexx (ring%ele_(1:n)%name, r_indexx(1:n))
      do k = 1, ring%ele_(i)%n_lord
        name = pring%ele(ixx)%name_(k)
        call find_indexx (name, ring%ele_(1:n)%name, r_indexx, n, j_lord)
        if (j_lord == 0) then
          call find_indexx (name, in_ring%ele_(:)%name, in_indexx, n_max, j)
          if (j == 0) then
            call warning ('UNABLE TO FIND COMPONENT FOR: ' // &
                            ring%ele_(i)%name, 'CANNOT FIND: ' // name)
            cycle
          endif
          ring%n_ele_max = ring%n_ele_max + 1
          j_lord = ring%n_ele_max 
          ring%ele_(j_lord) = in_ring%ele_(j)
          ring%ele_(j_lord)%control_type = component_lord$
          ring%ele_(j_lord)%ix1_slave = ring%n_control_array + 1
          ring%ele_(j_lord)%ix2_slave = ring%n_control_array + &
                                                    ring%ele_(j_lord)%n_slave
          ring%n_control_array = ring%ele_(j_lord)%ix2_slave
          ix2 = ring%ele_(j_lord)%ixx
          pring%ele(ix2)%ix_count = 0
          n = ring%n_ele_max
          call indexx (ring%ele_(1:n)%name, r_indexx(1:n))
        endif

        ix2 = ring%ele_(j_lord)%ixx
        pring%ele(ix2)%ix_count = pring%ele(ix2)%ix_count + 1
        ix = ring%ele_(j_lord)%ix1_slave + pring%ele(ix2)%ix_count - 1
        ring%control_(ix)%ix_slave = i
        ring%control_(ix)%ix_lord = j_lord

        ic = ring%ele_(i)%ic1_lord + k - 1
        ring%ic_(ic) = ix

      enddo
    endif
  enddo

! make matrices for entire ring

  call s_calc (ring)                       ! calc longitudinal distances
  call ring_geometry (ring)                ! ring layout
  call ring_make_mat6(ring, -1)            ! make 6x6 transport matrices

!-------------------------------------------------------------------------
! write out if debug is on

  if (pcom%parser_debug) then
    
    if (index(pcom%debug_line, 'VAR') /= 0) then
      print *
      print *, '----------------------------------------'
      print *, 'Number of Variables:', pcom%ivar_tot - pcom%ivar_init
      do i = pcom%ivar_init+1, pcom%ivar_tot
        print *
        print *, 'Var #', i-pcom%ivar_tot
        print *, 'Name: ', var_(i)%name
        print *, 'Value:', var_(i)%value
      enddo
    endif

    if (index(pcom%debug_line, 'SEQ') /= 0) then
      print *
      print *, '----------------------------------------'
      print *, 'Number of Lines/Lists defined:', pcom.iseq_tot
      do i = 1, pcom.iseq_tot
        print *
        print *, 'Sequence #', i
        print *, 'Name: ', seq_(i)%name
        print *, 'Type:', seq_(i)%type
        print *, 'Number of elements:', size(seq_(i)%ele)
        print *, 'List:'
        do j = 1, size(seq_(i)%ele(:))
          print '(4x, a, l3, i3)', seq_(i)%ele(j)%name, &
              seq_(i)%ele(j)%reflect, seq_(i)%ele(j)%ix_arg
        enddo
      enddo
    endif

    if (index(pcom%debug_line, 'SLAVE') /= 0) then
      print *
      print *, '----------------------------------------'
      print *, 'Number of Elements in Regular Ring:', ring%n_ele_ring
      do i = 1, ring%n_ele_ring
        print *, '-------------'
        print *, 'Ele #', i
        call type_ele (ring%ele_(i), .false., 0, .false., .true., ring)
      enddo
    endif

    if (index(pcom%debug_line, 'LORD') /= 0) then
      print *
      print *, '----------------------------------------'
      print *, 'LORD elements: ', ring%n_ele_max - ring%n_ele_ring
      do i = ring%n_ele_ring+1, ring%n_ele_max
        print *, '-------------'
        print *, 'Ele #', i
        call type_ele (ring%ele_(i), .false., 0, .false., .true., ring)
      enddo
    endif

    if (index(pcom%debug_line, 'RING') /= 0) then  
      print *
      print *, '----------------------------------------'
      print *, 'Ring Used: ', ring%name
      print *, 'Number of ring elements:', ring%n_ele_ring
      print *, 'List:                               Key      Length         S'
      do i = 1, ring%n_ele_ring
        type '(3x, i3, 2a, 3x, a, 2f10.2)', i, ') ', ring%ele_(i)%name,  &
          key_name(ring%ele_(i)%key), ring%ele_(i)%value(l$), ring%ele_(i)%s
      enddo
    endif

  endif

!-----------------------------------------------------------------------------
! deallocate pointers

  do i = lbound(pring%ele, 1) , ubound(pring%ele, 1)
    if (associated (pring%ele(i)%name_)) then
      deallocate(pring%ele(i)%name_)
      deallocate(pring%ele(i)%attrib_name_)
      deallocate(pring%ele(i)%cs_)
    endif
  enddo

  do i = 1, size(seq_(:))
    if (associated (seq_(i)%arg)) deallocate(seq_(i)%arg)
    if (associated (seq_(i)%ele)) deallocate(seq_(i)%ele)
  enddo



! error check

  if (pcom%error_flag) then
    if (bmad_status%exit_on_error) then
      print *, 'BMAD_PARSER FINISHED. EXITING ON ERRORS'
      call exit
    else
      bmad_status%ok = .false.
      return
    endif
  endif
              
  call check_ring_controls (ring, .true.)

! write to digested file

  if (.not. pcom%no_digested .and. .not. pcom%parser_debug .and. &
      digested_version <= bmad_inc_version$) call write_digested_bmad_file  &
               (digested_file, ring, pcom%n_files, pcom%file_name_)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_next_word (word, ix_word,
!                        delim_list, delim, delim_found, upper_case_word)
!
! Subroutine to get the next word from the input stream.
!
! Input:
!   word       -- Character*(*): Word returned
!   delim_list -- Character*(*): List of valid delimiters
!   upper_case_word -- Logical: if True then convert word to upper case.
!
! Output
!   ix_word     -- Integer: length of WORD
!   delim       -- Character*1: Actual delimiter found
!   delim_found -- Logical: Set true if a delimiter found. A delimiter
!                    may not be found if the end of the line is reached first.
!-


subroutine get_next_word (word, ix_word, delim_list, &
                    delim, delim_found, upper_case_word)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  integer ix_a, ix_word

  character*(*) word, delim_list, delim
                           
  logical delim_found, file_end, upper_case_word

! check for continuation character and if found then load more characters
! into the parse line.
! after that get the first word in PCOM.PARSE_LINE

  do
    ix_a = index(pcom%parse_line, '&')
    if (ix_a == 0 .or. ix_a > 70) exit
    call load_parse_line('continue', ix_a, file_end)
  enddo

  call word_read (pcom%parse_line, delim_list,  &
                         word, ix_word, delim, delim_found, pcom%parse_line)

  if (upper_case_word) call str_upcase (word, word)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine file_stack (how, file_name, finished)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  integer lunget, i_line_(10), f_unit_(10), i_level

  character*(*) how, file_name
  character*200 stack_file_name_(10)

  logical finished

  save i_level

!

  finished = .false.

  if (how == 'init') then
    i_level = 0
    return

  elseif (how == 'push') then
    i_level = i_level + 1
    if (i_level > 10) then
      print *, 'ERROR: CALL NESTING GREATER THAN 10 LEVELS'
      call err_exit
    endif
    stack_file_name_(i_level) = file_name
    pcom%current_file_name = file_name
    f_unit_(i_level) = lunget()
    pcom%f_unit = f_unit_(i_level)
    if (i_level /= 1) i_line_(i_level-1) = pcom%i_line
    pcom%i_line = 0
    open (unit = pcom%f_unit, file = pcom%current_file_name,  &
                                     status = 'old', readonly, err = 9000)
    pcom%n_files = pcom%n_files + 1
    inquire (file = file_name, name = pcom%file_name_(pcom%n_files))
  elseif (how == 'pop') then
    close (unit = pcom%f_unit)
    i_level = i_level - 1
    if (i_level < 0) then
      call error_exit ('BAD "RETURN"', ' ')
    elseif (i_level > 0) then
      pcom%current_file_name = stack_file_name_(i_level)
      pcom%f_unit = f_unit_(i_level)
      pcom%i_line = i_line_(i_level)
    else    ! i_level == 0
      finished = .true.
    endif
  else
    print *, 'BMAD_PARSER: INTERNAL ERROR IN FILE_STACK SUBROUTINE'
    call err_exit
  endif

  return

9000  continue
  if (bmad_status%type_out .or. bmad_status%exit_on_error) print *,  &
      'ERROR IN ', trim(pcom%parser_name), ': UNABLE TO OPEN FILE: ', &
                                       trim(pcom%current_file_name)
  if (bmad_status%exit_on_error) call err_exit
  bmad_status%ok = .false.

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine load_parse_line (how, ix_cmd, file_end)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  integer ix_cmd, ix

  character*(*) how
  character*80 line, pending_line

  logical cmd_pending, file_end


!

  file_end = .false.

! init

  if (how == 'init') then
    pcom%parser_debug = .false.
    pcom%no_digested = .false.
    read (pcom%f_unit, '(a)') line
    if (index(line, '!PARSER_DEBUG')) then
      pcom%parser_debug = .true.
      pcom%debug_line = line
      print *, 'FOUND IN FILE: "!PARSER_DEBUG". DEBUG IS NOW ON'
    elseif (index(line, '!NO_DIGESTED')) then
      pcom%no_digested = .true.
      print *, 'FOUND IN FILE: "!NO_DIGESTED". NO DIGESTED FILE WILL BE CREATED'
    endif
    rewind (unit = pcom%f_unit)
    return
  endif

  if (how /= 'normal' .and. how /= 'continue') call error_exit  &
                                     ('INTERNAL ERROR #4: CALL HELP', ' ')

!

1000    continue
    if (cmd_pending) then
      line = pending_line
      cmd_pending = .false.
    else
      read (pcom%f_unit, '(a)', end = 9000) line
      pcom%i_line = pcom%i_line + 1
    endif

! strip off comments

    ix = index(line, '!')
    if (ix == 1) then
      line = ' '
    elseif (ix > 1) then
      line = line(:ix-1)
    endif

! semi-colon delimiter means that we need to split the line
! and save the 2nd piece for the next time around.

    ix = index(line, ';')
    if (ix == 1) then
      cmd_pending = .true.
      pending_line = line(ix+1:)
      line = ' '
    elseif (ix > 1) then
      cmd_pending = .true.
      pending_line = line(ix+1:)
      line = line(:ix-1)
    else
      cmd_pending = .false.
    endif

! if the command line is blank then go back for more input

  call string_trim (line, line, ix)
  if (ix == 0 .and. how == 'normal') goto 1000

  pcom%parse_line(ix_cmd:) = line

  return

9000  continue
  file_end = .true.
  pcom%parse_line = ' '

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine evaluate_value (var, ring, final_delim, final_delim_found, err_flag)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  type (parser_var_struct)  var
  type (ring_struct)  ring

  integer i_lev, i_op, i

  integer plus$ / 1 /, minus$ / 2 /, times$ / 3 /, divide$ / 4 /
  integer l_parens$ / 5 /, power$ / 7 /, chs$ / 8 /, none$ / 9 /
  integer numeric$ / 100 /

  integer level(9) / 1, 1, 2, 2, 0, 0, 4, 3, -1 /
  character*1 op_name(9) / '+', '-', '*', '/', '(', ')', '^', '-', ' ' /

  integer stk_type(200), op_(200), ix_word, i_delim, i2, ix0

  real stk_value(200), value

  character line*70, delim*1, word*40, final_delim*1, word0*40

  logical parsing, delim_found, final_delim_found, op_pending, split
  logical err_flag

! init

  err_flag = .false.
  i_lev = 0
  i_op = 0

! get line to parse

  call get_next_word (line, ix_word, ',}', final_delim, final_delim_found, .true.)
  if (ix_word == 0) call warning  &
                         ('NO VALUE FOUND FOR: ' // var%name)

! parsing loop to build up the stack

  parsing = .true.
  do while(parsing)

! get a word

    call word_read (line, '+-*/()^',  &
                         word, ix_word, delim, delim_found, line)

    if (delim == '*' .and. line(1:1) == '*') then
      call warning ('EXPONENTIATION SYMBOL IS "^" AS OPPOSED TO "**"!',  &
                    'for: ' // var%name)
      err_flag = .true.
      return
    endif

! just make sure we are not chopping a number in two, e.g. "3.5e-7" should not
! get split at the "-" even though "-" is a delimiter

    split = .true.         ! assume initially that we have a split number
    if (ix_word == 0) then
      split = .false.
    elseif (word(ix_word:ix_word) /= 'E') then
      split = .false.
    endif
    if (delim(1:1) /= '-' .and. delim(1:1) /= '+') split = .false.
    do i = 1, ix_word-1
      if (index('.0123456789', word(i:i)) == 0) split = .false.
    enddo

! If still SPLIT = .TRUE. then we need to unsplit

    if (split) then
      word0 = word(:ix_word) // delim
      ix0 = ix_word + 1
      call word_read (line, '+-*/()^',  &
                         word, ix_word, delim, delim_found, line)
      word = word0(:ix0) // word
      ix_word = ix_word + ix0
    endif

! Now see what we got...

! for a "(" delim

    if (delim == '(') then
      if (ix_word /= 0) then
        call warning ('UNEXPECTED CHARACTERS ON RHS BEFORE "(": ' // word,  &
                                                  'FOR: ' // var%name)
        err_flag = .true.
        return
      endif
      call pushit (op_, i_op, l_parens$)
      op_pending = .false.

! for a unary "-"

    elseif (delim == '-' .and. ix_word == 0) then
      call pushit (op_, i_op, chs$)
      op_pending = .false.

! for a ")" delim

    elseif (delim == ')') then
      if (ix_word == 0) call error_exit  &
            ('CONSTANT OR VARIABLE MISSING FOR: ' // var%name, ' ')
      call word_to_value (word, ring, value)
      call pushit (stk_type, i_lev, numeric$)
      stk_value(i_lev) = value

2500  do i = i_op, 1, -1     ! release pending ops
        if (op_(i) == l_parens$) goto 3000  ! break do loop
        call pushit (stk_type, i_lev, op_(i))
      enddo
      call warning ('UNMATCHED ")" ON RHS', 'FOR: ' // var%name)
      err_flag = .true.
      return

3000  continue
      i_op = i - 1
      call word_read (line, '+-*/()^',  &
                         word, ix_word, delim, delim_found, line)
      if (ix_word /= 0) then
        call warning ('UNEXPECTED CHARACTERS ON RHS AFTER ")"',  &
                                                  'FOR: ' // var%name)
        err_flag = .true.
        return
      endif

      if (delim == ')') then
        goto 2500     ! if more ')' need to release more
      elseif (delim == '(') then
        call warning  &
            ('")(" CONSTRUCT DOES NOT MAKE SENSE FOR: ' // var%name)
        err_flag = .true.
        return
      else
        op_pending = .true.
      endif

! For binary "+-/*^" delims

    else
      if (ix_word == 0) then
        call warning ('CONSTANT OR VARIABLE MISSING FOR: ' // var%name)
        err_flag = .true.
        return
      endif
      call word_to_value (word, ring, value)
      call pushit (stk_type, i_lev, numeric$)
      stk_value(i_lev) = value
      op_pending = .true.
    endif

! OP_PENDING = .true. means that we have an operation that is waiting to
! be identified

    if (op_pending) then
      do i = 1, 7
        if (delim == op_name(i)) then      ! op identified
          i_delim = i                        ! op id number
          goto 4000
        endif
      enddo

      if (delim_found) then   ! how could this be?
        call error_exit ('INTERNAL ERROR #01: GET HELP', ' ')
      else                    ! must be that we are at the end of the line
        i_delim = none$
      endif

! now see if there are operations on the OP_ stack that need to be transferred
! to the STK_ stack

4000  continue

      do i = i_op, 1, -1
        if (level(op_(i)) >= level(i_delim)) then
          call pushit (stk_type, i_lev, op_(i))
        else
          goto 5000
        endif
      enddo

! put the pending operation on the OP_ stack

5000  continue
      i_op = i
      if (i_delim == none$) then
        goto 6000
      else
        call pushit (op_, i_op, i_delim)
      endif

    endif

  enddo

! now go through the stack and perform the operations

6000  continue
  if (i_op /= 0) then
    call warning ('UNMATCHED "(" ON RHS', 'FOR: ' // var%name)
    err_flag = .true.
    return
  endif

  i2 = 0
  do i = 1, i_lev
    if (stk_type(i) == numeric$) then
      i2 = i2 + 1
      stk_value(i2) = stk_value(i)
    elseif (stk_type(i) == chs$) then
      stk_value(i2) = -stk_value(i2)
    elseif (stk_type(i) == plus$) then
      stk_value(i2-1) = stk_value(i2-1) + stk_value(i2)
      i2 = i2 - 1
    elseif (stk_type(i) == minus$) then
      stk_value(i2-1) = stk_value(i2-1) - stk_value(i2)
      i2 = i2 - 1
    elseif (stk_type(i) == times$) then
      stk_value(i2-1) = stk_value(i2-1) * stk_value(i2)
      i2 = i2 - 1
    elseif (stk_type(i) == divide$) then
      if (stk_value(i2) == 0) call error_exit  &
              ('DIVIDE BY 0 ON RHS', 'FOR: ' // var%name)
      stk_value(i2-1) = stk_value(i2-1) / stk_value(i2)
      i2 = i2 - 1
    elseif (stk_type(i) == power$) then
      stk_value(i2-1) = stk_value(i2-1)**stk_value(i2)
      i2 = i2 - 1
    else
      call error_exit ('INTERNAL ERROR #02: GET HELP', ' ')
    endif
  enddo


  if (i2 /= 1) call error_exit ('INTERNAL ERROR #03: GET HELP', ' ')

  var%value = stk_value(1)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine increment_pointer (ix, reflect)

  implicit none

  integer ix, reflect

!

  if (reflect > 0) then
    ix = ix + 1
  else
    ix = ix - 1
  endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine pushit (stack, i_lev, value)

  implicit none

  integer stack(:), i_lev, value

!

  i_lev = i_lev + 1
  stack(i_lev) = value

end subroutine
                       
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine word_to_value (word, ring, value)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  type (ring_struct)  ring

  integer i, i_ele, ix1, ix2
  real value
  character*(*) word
  character*16 name

! see if this is numeric

  if (index('-+.0123456789', word(1:1)) /= 0) then
    read (word, *, err = 9000) value
    return
  endif

!

  call verify_valid_name(word, len_trim(word))

! If word has a "[...]" then it is a element attribute

  ix1 = index(word, '[')
  if (ix1 /= 0) then   

    name = word(:ix1-1)    ! name of attribute
    i_ele = 0
    do i = 0, ring%n_ele_max
      if (ring%ele_(i)%name == name) then
        i_ele = i
        goto 1000
      endif
    enddo

    call warning ('ELEMENT NOT DEFINED: ' // name)
    value = 0
    return

1000    continue

    ix2 = index(word, ']')
    name = word(ix1+1:ix2-1)

    if (name == 'S') then
      if (pcom%parser_name == 'BMAD_PARSER2') then
        value = ring%ele_(i_ele)%s
      else
        call warning ('"S" ATTRIBUTE CAN ONLY BE USED WITH BMAD_PARSER2')
      endif
    else
      i = attribute_index(ring%ele_(i_ele), name)
      if (i < 1) call warning('BAD ATTRIBUTE NAME: ' // word)
      value = ring%ele_(i_ele)%value(i)
    endif

    return
  endif

! None of the above? must be a variable

  do i = 1, pcom%ivar_tot
    if (word == var_(i)%name) then
      value = var_(i)%value
      return
    endif
  enddo

  call warning ('VARIABLE USED BUT NOT YET DEFINED: ' // word)
  return

9000  call warning ('BAD VARIABLE: ' // word)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine get_attribute (how, ele, ring, pring, &
                                             delim, delim_found, err_flag)
                            
  use local_bmad_struct
  use local_bmad_interface

  implicit none

  type (parser_var_struct)  var
  type (ring_struct)  ring
  type (parser_ring_struct) pring
  type (ele_struct)  ele, ele0

  integer i, ic, ix_word, how

  character*16 word, tilt_word / 'TILT' /
  character delim*1
  character*16 super_names(11) / 'SUPERIMPOSE', 'OFFSET', 'REFERENCE',  &
                'ELE_BEGINNING', 'ELE_CENTER', 'ELE_END',  &
                'REF_BEGINNING', 'REF_CENTER', 'REF_END', &
                'COMMON_LORD', ' ' /

  logical delim_found, err_flag

! Get next WORD.
! If an overlay or group element then word is just an attribute to control
! [except for a "GROUP[COMMAND] = 0.343" redef construct]

  err_flag = .false.
  call get_next_word (word, ix_word, ':, =', delim, delim_found, .true.)

  if (ele%key == overlay$) then
    i = attribute_index(ele, word)       ! general attribute search
    if (i < 1) then
      call warning ('BAD OVERLAY ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
      err_flag = .true.
      return
    endif
    if (i == type$ .or. i == alias$) then
      call type_get (ele, i, delim, delim_found)
    else
      if (how == def$) then
        ele%ix_value = i
        ele%attribute_name = word
      endif
      if (delim == '=') then  ! value
        var%name = ele%name(:8) // ':' // word ! in case there is an error
        call evaluate_value (var, ring, delim, delim_found, err_flag)
        if (err_flag) return
        ele%value(i) = var%value
      else
        ele%value(i) = 0
      endif
    endif
    return

  elseif (ele%key == group$) then

    if (how == def$) then
      ele0%key = overlay$
      i = attribute_index(ele0, word)       ! general attribute search
    else   ! how == redef$
      i = attribute_index(ele, word)
    endif

    if (i < 1) then
      call warning ('BAD GROUP ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
      err_flag = .true.
      return
    endif

    if (i == type$ .or. i == alias$) then
      call type_get (ele, i, delim, delim_found)
    else
      if (how == def$) ele%ix_value = i
      if (delim == '=') then  ! value
        var%name = ele%name(:8) // ':' // word ! in case there is an error
        call evaluate_value (var, ring, delim, delim_found, err_flag)
        if (err_flag) return
        if (how == def$) then
          ele%value(command$) = var%value
        else
          ele%value(i) = var%value
        endif
      elseif (how == redef$) then
        call warning ('NO VALUE GIVEN FOR ATTRIBUTE FOR: ' // ele%name)
        err_flag = .true.
      endif
    endif

    return

  endif

! if not an overlay then see if it is an ordinary attribute.
! if not an ordinary attribute then might be a superimpose switch

  i = attribute_index(ele, word)

  if (i < 1) then          ! if not an ordinary attribute...
    if (ix_word == 0) then  ! no word
      call warning  &
            ('"," NOT FOLLOWED BY ATTRIBUTE NAME FOR: ' // ele%name)
      err_flag = .true.
      return
    else
      if (word(:ix_word) == 'REF') word = 'REFERENCE' ! allowed abbrev
      call match_word (word, super_names, i)
      if (i < 1) then
        call warning  &
            ('BAD ATTRIBUTE NAME: ' // word, 'FOR ELEMENT: ' // ele%name)

        err_flag = .true.
        return
      else    ! valid superimpose switch

        if (ele%ixx == 0) then
          call warning ('ELEMENT HAS NO ASSOCIATED INFO: ' // ele%name) 
          return
        endif
        ic = ele%ixx
        if (super_names(i) == 'SUPERIMPOSE') then
          ele%control_type = super_lord$
        elseif (super_names(i) == 'REF_BEGINNING') then
          pring%ele(ic)%ref_pt = begin$
        elseif (super_names(i) == 'REF_CENTER') then
          pring%ele(ic)%ref_pt = center$
        elseif (super_names(i) == 'REF_END') then
          pring%ele(ic)%ref_pt = end$
        elseif (super_names(i) == 'ELE_BEGINNING') then
          pring%ele(ic)%ele_pt = begin$
        elseif (super_names(i) == 'ELE_CENTER') then
          pring%ele(ic)%ele_pt = center$
        elseif (super_names(i) == 'ELE_END') then
          pring%ele(ic)%ele_pt = end$
        elseif (super_names(i) == 'COMMON_LORD') then
          pring%ele(ic)%common_lord = .true.
        elseif (super_names(i) == 'REFERENCE') then
          call get_next_word(pring%ele(ic)%ref_name, ix_word,  &
                                             ':=,', delim, delim_found, .true.)
        elseif (super_names(i) == 'OFFSET') then
          var%name = ele%name(:8) // ':' // word   ! in case of error
          call evaluate_value (var, ring, delim, delim_found, err_flag)
          if (err_flag) return
          pring%ele(ic)%s = var%value
        else
          print *, 'ERROR IN BMAD_PARSER: INTERNAL ERROR. PLEASE GET HELP!'
          call err_exit
        endif
      endif
      return
    endif
  endif

! check that next delim is a "=". If not check for a possible default value
! otherwise it is an error

1000  continue

  if (delim /= '=')  then
    if (word == tilt_word) then
      if (ele%key == quadrupole$) then
        ele%value(tilt$) = pi / 4
      elseif (ele%key == sextupole$) then
        ele%value(tilt$) = pi / 6
      elseif (ele%key == octupole$) then
        ele%value(tilt$) = pi / 8
      else
        call warning ('SORRY I''M NOT PROGRAMMED TO USE A "TILT" DEFAULT',  &
                'FOR A: ' // key_name(ele%key))
        err_flag = .true.
        return
      endif
      return
    else
      call warning ('EXPECTING "=" AFTER ATTRIBUTE: ' // word,  &
                         'FOR ELEMENT: ' // ele%name)
      err_flag = .true.
      return
    endif
  endif

! get the value of the attribute.
! The TYPE attribute is special because its "value" is a character string

  if (i == type$ .or. i == alias$) then         ! if TYPE attribute
    call type_get (ele, i, delim, delim_found)
  else    ! normal attribute
    var%name = ele%name(:8) // ':' // word   ! in case there is an error
    call evaluate_value (var, ring, delim, delim_found, err_flag)
    if (err_flag) return
    if (i == mat6_calc_method$) then
      ele%mat6_calc_method = nint(var%value)
    elseif (i == tracking_method$) then   
      ele%tracking_method = nint(var%value)
    elseif (i == num_steps$) then    
      ele%num_steps = nint(var%value)
    else
      ele%value(i) = var%value
    endif
  endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine type_get (ele, ix_type, delim, delim_found)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  type (ele_struct)  ele

  integer ix, ix_word, ix_type
  character word*16, delim*1, type_name*16
  logical delim_found

!

  call string_trim(pcom%parse_line, pcom%parse_line, ix)
  ix = index(pcom%parse_line(2:), '"')

  if (pcom%parse_line(1:1) /= '"' .or. ix == 0) then
    call warning ('MISSING QUOTE MARK (") FOR TYPE = "attribute"',  &
                          'FOR ELEMENT: ' // ele%name)
    if (ix /= 0) then
      pcom%parse_line = pcom%parse_line(ix+2:)
    elseif (pcom%parse_line(1:1) == '"') then
      pcom%parse_line = pcom%parse_line(2:)
    endif
  endif

  if (ix == 1) then
    type_name = ' '
  else
    type_name = pcom%parse_line(2:ix)
    call str_upcase (type_name, type_name)
  endif

  if (ix_type == type$) then
    ele%type = type_name
  elseif (ix_type == alias$) then
    ele%alias = type_name
  endif

  pcom%parse_line = pcom%parse_line(ix+2:)
  call get_next_word (word, ix_word, ',=', delim, delim_found, .true.)
  if (ix_word /= 0) call warning (  &
                'EXTRA CHARACTERS FOUND AFTER TYPE ATTRIBUTE: ' // word,  &
                'FOR ELEMENT: ' // ele%name)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine GET_OVERLAY_GROUP_NAMES (ELE, RING, PRING, DELIM, DELIM_FOUND)
!-
        
subroutine get_overlay_group_names (ele, ring, pring, delim, delim_found)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  type (ele_struct)  ele
  type (parser_ring_struct) pring
  type (ring_struct)  ring
  type (parser_var_struct)  var
  type (control_struct) cs_(100)

  integer ic, ix_word, ixs, j, k
                             
  character delim*1, word_in*40, word*40
  character*16 name_(100), attrib_name_(100)

  logical delim_found, parsing, err_flag, file_end
                      
!

  call get_next_word (word_in, ix_word, '{,}', delim, delim_found, .true.)
  if (delim /= '{' .or. ix_word /= 0) call warning  &
          ('BAD ' // control_name(ele%control_type) // 'SPEC: ' // word_in,  &
          'FOR ELEMENT: ' // ele%name)

!

  parsing = .true.
  do while (parsing)

    call get_next_word (word_in, ix_word, '{,}/', delim, delim_found, .true.)
    ele%n_slave = ele%n_slave + 1
    ixs = ele%n_slave
    word = word_in

    j = index(word, '[')
    if (j > 1) then
      k = index(word, ']')
      if (k <= j+1) then
        call warning ('BAD ATTRIBUTE SPEC: ' // word_in, 'FOR: ' // ele%name)
        word = word(:k-1) // word(j+1:)
      else
        attrib_name_(ixs) = word(j+1:k-1)
        word = word(:j-1) // word(k+1:)
      endif
    else
      attrib_name_(ixs) = blank
    endif

    name_(ixs) = word

    if (delim == '/') then
      var%name = ele%name
      call evaluate_value (var, ring, delim, delim_found, err_flag)
      if (err_flag) then
        call warning ('BAD COEFFICIENT: ' // word_in,  &
                                          'FOR ELEMENT: ' // ele%name)
        call load_parse_line ('normal', 1, file_end)         ! next line
        return
      endif
      cs_(ixs)%coef = var%value
    else
      cs_(ixs)%coef = 1.0
    endif

    if (delim == '}') then
      parsing = .false.
      call get_next_word (word, ix_word, ',=:', delim, delim_found, .true.)
    elseif (delim /= ',') then
      call warning ('BAD ' // control_name(ele%control_type) //  &
              'SPEC: ' // word_in, 'FOR: ' // ele%name)
      parsing = .false.
    endif
                          
  enddo

!

  ic = ele%ixx
  allocate (pring%ele(ic)%cs_(ixs), pring%ele(ic)%name_(ixs), &
                                       pring%ele(ic)%attrib_name_(ixs))
  pring%ele(ic)%cs_ = cs_(1:ixs)
  pring%ele(ic)%name_ = name_(1:ixs)
  pring%ele(ic)%attrib_name_ = attrib_name_(1:ixs)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine verify_valid_name (name, ix_name)

  use local_bmad_interface

  implicit none

  integer i, ix_name, len, ix1, ix2

  character*(*) name
  character*27 letters / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\' /
  character*40   valid_chars / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\0123456789_[]' /
  character*1 tab

  logical OK

  parameter (tab = char(9))

! check for blank spaces

  do i = 1, min(ix_name, len(name))
    if (name(i:i) == ' ' .or. name(i:i) == tab) call warning  &
                          ('NO DELIMITER BETWEEN NAMES: ' // name)
  enddo

! check for name too long

  if (ix_name > len(name)) then
     call warning ('NAME TOO LONG: ' // name)
     ix_name = len(name)      ! chop name
  endif

! check for name too short

  if (ix_name == 0) call warning ('BLANK NAME')

! check for invalid characters in name

  OK = .true.
  if (index(letters, name(1:1)) == 0) OK = .false.
  do i = 2, ix_name
    if (index(valid_chars, name(i:i)) == 0) OK = .false.
  enddo

  if (.not. OK) call warning ('INVALID NAME: UNRECOGNIZED CHARACTERS IN: '  &
                                   // name)

! check for non matched "[" "]" pairs

  ix1 = index(name, '[')
  ix2 = index(name, ']')
  if (ix1 /= 0 .or. ix2 /= 0) then
    if (ix1 == 0) call warning ('UNMATCHED BRACKET: ' // name)
    if (ix2 <= ix1+1) call warning  &
                    ('INVALID: REVERSED BRACKETS: ' // name)
    if (index(name(ix1+1:), '[') /= 0 .or. index(name(ix2+1:), ']') /=  &
                   0) call warning ('INVALID: BAD BRACKETS: ' // name)
    if (ix2 /= len(name)) then
      if (name(ix2+1:ix2+1) /= ' ') call warning  &
                    ('INVALID: NOTHING IN BRACKETS: ' // name)
    endif
  endif

! check for more than 16 characters

  if (ix1 == 0 .and. ix_name > 16)  &
            call warning ('NAME HAS > 16 CHARACTERS: ' // name)

  if (ix1 > 17 .or. ix2 - ix1 > 17)  &
            call warning ('NAME HAS > 16 CHARACTERS: ' // name)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine error_exit (what1, what2)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  character*(*) what1, what2

!

  print *, 'FATAL ERROR IN BMAD_PARSER: ', what1
  if (what2 /= ' ') type '(22x, a)', what2
  print *, '      IN FILE: ', pcom%current_file_name(:60)
  print *, '      AT OR BEFORE LINE:', pcom%i_line
  call err_exit

end subroutine


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine warning (what1, what2)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  integer ix
  character*(*) what1
  character*(*), optional :: what2

! PCOM.ERROR_FLAG is a common logical used so program will stop at end of parsing

  if (bmad_status%type_out) then
    print *, 'ERROR IN ', trim(pcom%parser_name), ': ', what1
    if (present(what2)) type '(22x, a)', what2
    print *, '      IN FILE: ', pcom%current_file_name(:60)
    print *, '      AT OR BEFORE LINE:', pcom%i_line
  endif

  pcom%error_flag = .true.
  bmad_status%ok = .false.

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine init_bmad_parser_common

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  var_(1)%name = 'PI'
  var_(1)%value = pi

  var_(2)%name = 'TWOPI'
  var_(2)%value = twopi

  var_(3)%name = 'DEGRAD'
  var_(3)%value= 180 / pi

  var_(4)%name = 'RADDEG'
  var_(4)%value = pi / 180

  var_(5)%name = 'E'
  var_(5)%value = 2.718281828459

  var_(6)%name = 'E_MASS'
  var_(6)%value = e_mass

  var_(7)%name = 'C_LIGHT'
  var_(7)%value = c_light

  var_(8)%name = 'POSITRON'
  var_(8)%value = positron$

  var_(9)%name = 'ELECTRON'
  var_(9)%value = electron$

  var_(10)%name = 'MOBIUS'
  var_(10)%value = mobius_symmetry$

  var_(11)%name = 'E_CHARGE'
  var_(11)%value = e_charge

  var_(12)%name = 'EMASS'      ! old style
  var_(12)%value = e_mass

  var_(13)%name = 'CLIGHT'     ! old style
  var_(13)%value = c_light

  var_(14)%name = 'LINAC_LATTICE'
  var_(14)%value = linac_lattice$

  var_(15)%name = 'LINEAR_LATTICE'
  var_(15)%value = linear_lattice$

  var_(16)%name = 'CIRCULAR_LATTICE'
  var_(16)%value = circular_lattice$

  var_(17)%name = 'BMAD_STANDARD'
  var_(17)%value = bmad_standard$

  var_(18)%name = 'DA_MAP'
  var_(18)%value = DA_map$

  var_(19)%name = 'CUSTOM_CALC'
  var_(19)%value = custom_calc$

  var_(20)%name = 'RUNGE_KUTTA'
  var_(20)%value = runge_kutta$

  var_(21)%name = 'PROTON'
  var_(21)%value = proton$

  var_(22)%name = 'ANTIPROTON'
  var_(22)%value = antiproton$

  pcom%ivar_tot = 20
  pcom%ivar_init = 20

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine compute_super_lord_s (ele, ring, pring)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  type (ring_struct)   ring
  type (ele_struct)    ele
  type (parser_ring_struct) pring

  integer ic, k
  logical found

! compute position of super_lord element
! simple case

  ic = ele%ixx

  if (ring%n_ele_max < ring%n_ele_ring) then
    call warning ('N_ELE_MAX LESS THAN N_ELE_RING!')
    call err_exit
  endif

  if (pring%ele(ic)%ref_name == blank) then
    call compute2_super_lord_s (ring, 0, ele, pring%ele(ic))
    return
  endif

! search

  found = .false.

  do k = 1, ring%n_ele_max
    if (pring%ele(ic)%ref_name == ring%ele_(k)%name) then
      if (found) then
        call warning ('MULTIPLE NAME MATCHES FOR' //  &
                  ' REFERENCE OF SUPERIMPOSE: ' // ele%name)
        return
      else
        call compute2_super_lord_s (ring, k, ele, pring%ele(ic))
        found = .true.
      endif
    endif
  enddo

  if (.not. found) call warning ('UNKNOWN REFERENCE ELEMENT: ' //  &
        pring%ele(ic)%ref_name, 'FOR SUPERIMPOSE ELEMENT: ' // ele%name)

end subroutine
                                    
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine add_all_superimpose (ring, ele_in, pele)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele, ele_in, ele2, this_ele
  type (parser_ele_struct) pele

  integer ix, i, j, it, nic, nn, i_sup, i_ele, ic
  integer n_inserted, ix_lord, n_con

  character*16 matched_name(200)

  logical match_wild, have_inserted, create_new

! init

  ele = ele_in   ! in case ele changes
  ele2 = ele
  n_inserted = 0
  ring%ele_%iyy = 0    ! to keep track of where we have inserted

! If no refrence point then superposition is simple

  if (pele%ref_name == blank) then
    call compute2_super_lord_s (ring, 0, ele2, pele)
    call add_superimpose (ring, ele2, i_sup)
    return
  endif

! insert ele in the ring
! do not insert twice at the same spot

  do 

    have_inserted = .false.

    ele_loop: do i_ele = 1, ring%n_ele_max

      if (ring%ele_(i_ele)%control_type /= free$)  &
                                   call control_bookkeeper (ring, i_ele)
      this_ele = ring%ele_(i_ele)
      ic = this_ele%control_type
       
      if (ic == group_lord$ .or. ic == super_slave$) cycle
      if (this_ele%iyy == 1) cycle

      if (match_wild(this_ele%name, pele%ref_name)) then

        do i = 1, n_inserted
          if (this_ele%name == matched_name(i)) cycle ele_loop
        enddo
       
        ring%ele_(i_ele)%iyy = 1
        call compute2_super_lord_s (ring, i_ele, ele2, pele)
        call string_trim(ele%name, ele%name, ix)
        ele2%name = ele%name(:ix)            
        call add_superimpose (ring, ele2, i_sup)

        do i = i_ele, ring%n_ele_ring
          ring%ele_(i)%s = ring%ele_(i-1)%s + ring%ele_(i)%value(l$)
        enddo

        n_inserted = n_inserted + 1
        matched_name(n_inserted) = ele2%name
        have_inserted = .true.   

      endif

    enddo ele_loop

    if (.not. have_inserted) exit

  enddo

! error check

  if (n_inserted == 0) call warning ('NO MATCH FOR REFERENCE ELEMENT: ' //  &
            pele%ref_name, 'FOR SUPERPOSITION OF: ' // ele%name)


! if there is to be no common lord then we are done

  if (.not. pele%common_lord) return

! here for common_lord, not scalled multipoles

  if (ele%key /= multipole$ .and. ele%key /= ab_multipole$) then
    type *, 'ERROR IN INSERT_MUTIPLE: ELEMENT ', ring%ele_(i)%name
    type *, '      IS USED WITH THE "COMMON_LORD" ATTRIBUTE BUT'
    type *, '      THIS ELEMENT IS NOT A MULIPOLE OR AB_MULTIPOLE'
    call err_exit
  endif

  ring%n_ele_max = ring%n_ele_max + 1
  nn = ring%n_ele_max 

  n_con = ring%n_control_array 
  ring%n_control_array = n_con + n_inserted

  ring%ele_(nn) = ele
  ring%ele_(nn)%control_type = super_lord$
  ring%ele_(nn)%n_slave = n_inserted
  ring%ele_(nn)%ix1_slave = n_con + 1
  ring%ele_(nn)%ix2_slave = n_con + n_inserted

  do i = n_con + 1, n_con + n_inserted
    ring%control_(i)%ix_lord = nn
    ring%control_(i)%ix_attrib = 0
  enddo

  j = 0
  do i = 1, ring%n_ele_max-1
    if (any (matched_name(1:n_inserted) == ring%ele_(i)%name)) then
      it = ring%ele_(i)%control_type
      if (it /= free$) then
        type *, 'ERROR IN INSERT_MUTIPLE: SLAVE ', ring%ele_(i)%name
        type *, '      OF LORD ', ele%name
        type *, '      IS NOT A "FREE" ELEMENT BUT IS: ', control_name(it)
        call err_exit
      endif
      j = j + 1
      ring%ele_(i)%control_type = super_slave$
      nic = ring%n_ic_array + 1
      ring%ele_(i)%n_lord = 1
      ring%ele_(i)%ic1_lord = nic
      ring%ele_(i)%ic2_lord = nic
      ring%ic_(nic) = n_con + j
      ring%control_(n_con+j)%ix_slave = i
      ring%n_ic_array = nic
    endif
  enddo

  if (j /= n_inserted) then
    type *, 'ERROR IN INSERT_MUTIPLE: SLAVE NUMBER MISMATCH', j, n_inserted
    call err_exit
  endif


end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine compute2_super_lord_s (ring, i_ref, ele, pele)

  use local_bmad_struct
  use local_bmad_interface

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele
  type (parser_ele_struct) pele

  integer i_ref, i, ix

  real s_ref_begin, s_ref_end

!

  ele%s = pele%s

  if (pele%ele_pt == begin$) then
    ele%s = ele%s + ele%value(l$)
  elseif (pele%ele_pt == center$) then
    ele%s = ele%s + ele%value(l$) / 2
  elseif (pele%ele_pt /= end$) then
    print *, 'ERROR IN COMPUTE2_SUPER_LORD_S: CONTROL #1 INTERNAL ERROR!'
    call err_exit
  endif

!

  if (ring%ele_(i_ref)%control_type == overlay_lord$) then
    s_ref_begin = 1e10
    s_ref_end = 0
    do i = ring%ele_(i_ref)%ix1_slave, ring%ele_(i_ref)%ix2_slave
      ix = ring%control_(i)%ix_slave
      s_ref_begin = min(s_ref_begin,  &
                         ring%ele_(ix)%s - ring%ele_(ix)%value(l$))
      s_ref_end = max(s_ref_end, ring%ele_(ix)%s)
    enddo
  elseif (ring%ele_(i_ref)%control_type == group_lord$) then
    call warning ('SUPERPOSING: ' // ele%name, 'UPON GROUP' // pele%ref_name)
    return
  else
    s_ref_begin = ring%ele_(i_ref)%s - ring%ele_(i_ref)%value(l$)
    s_ref_end = ring%ele_(i_ref)%s
  endif

  if (pele%ref_pt == begin$) then
    ele%s = ele%s + s_ref_begin
  elseif (pele%ref_pt == center$) then
    ele%s = ele%s + (s_ref_begin + s_ref_end) / 2
  elseif (pele%ref_pt == end$) then
    ele%s = ele%s + s_ref_end
  else
    print *, 'ERROR IN COMPUTE_SUPER_LORD_S: CONTROL #2 INTERNAL ERROR!'
    call err_exit
  endif

!

  if (ele%s > ring%param%total_length) then
    ele%s = ele%s - ring%param%total_length
  elseif (ele%s < 0) then
    ele%s = ele%s + ring%param%total_length
  endif

end subroutine
        
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine find_indexx (name, names, an_indexx, n_max, ix_match)


  implicit none

  integer ix1, ix2, ix3, n_max, ix_match
  integer an_indexx(:)

  character*16 name, names(:)
  character*16 this_name

! simple case

  if (n_max == 0) then
    ix_match = 0
    return
  endif

!

  ix1 = 1
  ix3 = n_max

  do

    ix2 = (ix1 + ix3) / 2 
    this_name = names(an_indexx(ix2))

    if (this_name == name) then
      ix_match = an_indexx(ix2)
      return
    elseif (this_name < name) then
      ix1 = ix2 + 1
    else
      ix3 = ix2 - 1
    endif
                       
    if (ix1 > ix3) then
      ix_match = 0
      return
    endif

  enddo

end subroutine
