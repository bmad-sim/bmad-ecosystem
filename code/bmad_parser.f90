!+
! Subroutine bmad_parser (in_file, ring, make_mats6, digested_read_ok)
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
!   use bmad
!
! Input:
!   in_file    -- Character: Name of the input file.
!   make_mats6 -- Logical, optional: Make the 6x6 transport matrices for then
!                   Elements? Default is True.
!
! Output:
!   ring -- Ring_struct: Ring structure. See bmad_struct.f90 for more details.
!     %ele_(:)%mat6  -- This is computed assuming an on-axis orbit 
!     %ele_(:)%s     -- This is also computed.
!   digested_read_ok -- Logical, optional: Set True if the digested file was
!                        successfully read. False otherwise.
!         
! Defaults:
!   ring%n_ele_symm              = 0
!   ring%param%particle          = positron$
!   ring%param%symmetry          = no_symmetry$
!   ring%param%lattice_type      = circular_lattice$
!   ring%param%aperture_limit_on = .true.
!   ring%param%damping_on        = .false.
!   ring%param%radiation_on      = .false.
!
! For more info on the parser see the BMAD documentation.
! DCS 10/6/97
!-

#include "CESR_platform.inc"

subroutine bmad_parser (in_file, ring, make_mats6, digested_read_ok)

  use bmad_parser_mod
  use cesr_utils
  
  implicit none

  type (ring_struct), target :: ring, in_ring
  type (seq_struct), target :: this_seq(200)
  type (seq_struct), pointer :: seq, seq2
  type (seq_ele_struct), pointer :: s_ele
  type (parser_ring_struct) pring
  type (seq_stack_struct) stack(20)
  type (ele_struct), save, pointer :: ele, old_ele(:) => null()
  type (control_struct), pointer :: cs_(:) => null()

  integer ix_word, i_use, i, j, k, n, ix, ixr, ixs, i_ring, it
  integer i_lev, i_key, ic, ix_lord
  integer iseq, ix_ring(n_ele_maxx), n_ele_ring
  integer ix_super, digested_version, key, ct
  integer ix1, ix2, iv, n_arg
  integer ivar, ixx, j_lord, n_slave
  integer seq_indexx(n_ele_maxx), in_indexx(n_ele_maxx)
  integer, pointer :: n_max

  character*(*) in_file
  character(16) word_2, name, a_name, dummy_name(20)
  character(16) name_(n_ele_maxx)
  character(1) delim*1
  character(200) path, basename, full_name, digested_file, call_file
  character(40) this_name, word_1
  character(280) parse_line_save

  real(rdef) angle, old_energy

  logical, optional :: make_mats6, digested_read_ok
  logical parsing, delim_found, matched_delim, arg_list_found, doit
  logical file_end, found, err_flag, finished, save_taylor
  logical vmask(n_attrib_maxx)
  logical, save :: init_needed = .true.

! see if digested file is open and current. If so read in and return.
! Note: The name of the digested file depends upon the real precision.

  bp_com%parser_name = 'BMAD_PARSER'  ! Used for error messages.
  bp_com%write_digested = .true.


  inquire (file = in_file, name = full_name)      ! full input file_name
  ix = index(full_name, ';')
  if (ix /= 0) full_name = full_name(:ix-1)
  ring%input_file_name = full_name      ! needed by read_digested_bmad_file

  ix = SplitFileName(in_file, path, basename)
  if (rdef == 8) then
    digested_file = in_file(:ix) // 'digested8_' // in_file(ix+1:)
  else
    digested_file = in_file(:ix) // 'digested_' // in_file(ix+1:)
  endif

  call read_digested_bmad_file (digested_file, ring, digested_version)

  if (bmad_status%ok) then
    call set_taylor_order (ring%input_taylor_order, .false.)
    call set_ptc (ring%param)
    if (ring%input_taylor_order == bmad_com%taylor_order) then
      if (present(digested_read_ok)) digested_read_ok = .true.
      return
    else
      if (bmad_status%type_out) &
                  print *, 'BMAD_PARSER: Taylor_order has changed.'
      if (ring%input_taylor_order > bmad_com%taylor_order) &
                                           bp_com%write_digested = .false.
    endif
  endif

  if (present(digested_read_ok)) digested_read_ok = .false.

! save all elements that have a taylor series

  old_energy = ring%param%beam_energy

  ix = 0
  do i = 1, ring%n_ele_max
    if (associated(ring%ele_(i)%taylor(1)%term)) ix = ix +1
  enddo

  if (ix /= 0) then
    if (associated(old_ele)) deallocate(old_ele)
    allocate(old_ele(ix))
    ix = 0
    do i = 1, ring%n_ele_max
      if (associated(ring%ele_(i)%taylor(1)%term)) then
        ix = ix +1
        old_ele(ix) = ring%ele_(i)
      endif
    enddo
  endif  

! init

  if (init_needed) then
    do i = lbound(pring%ele, 1) , ubound(pring%ele, 1)
      nullify (pring%ele(i)%name_)
    enddo
    do i = 1, size(this_seq(:))
      nullify (this_seq(i)%arg, this_seq(i)%ele)
    enddo
    init_needed = .false.
  endif

! here if not OK bmad_status. So we have to do everything from scratch...
! init variables.

  bmad_status%ok = .true.
  if (bmad_status%type_out) &
                        print *, 'BMAD_PARSER: Creating new digested file...'

  bp_com%n_files = 0
  bp_com%error_flag = .false.                 ! set to true on an error
  call file_stack('init', in_file, finished)   ! init stack
  call file_stack('push', in_file, finished)   ! open file on stack
  if (.not. bmad_status%ok) return
  call load_parse_line ('init', 0, file_end) ! initialize subroutine
  bp_com%iseq_tot = 0                     ! number of sequences encountered
  bp_com%ivar_tot = 0                     ! number of variables encountered
  call init_bmad_parser_common
  ring%name = ' '
  ring%lattice = ' '

  pring%ele(:)%ref_name = blank
  pring%ele(:)%ref_pt  = center$
  pring%ele(:)%ele_pt  = center$
  pring%ele(:)%s       = 0
  pring%ele(:)%common_lord = .false.

  call init_ele (in_ring%ele_(0))
  in_ring%ele_(0)%name = 'BEGINNING'     ! Beginning element
  in_ring%ele_(0)%key = init_ele$

  n_max => in_ring%n_ele_max
  n_max = 0                              ! Number of elements encountered

  call init_ele (beam_ele)
  beam_ele%name = 'BEAM'                 ! fake beam element
  beam_ele%key = def_beam$               ! "definition of beam"
  beam_ele%value(particle$) = positron$  ! default

  call init_ele (param_ele)
  param_ele%name = 'PARAMETER'    ! For parameters 
  param_ele%key = def_parameter$
  param_ele%value(lattice_type$) = circular_lattice$  ! Default
  param_ele%value(symmetry$)     = no_symmetry$       ! Default

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
      read (bp_com%f_unit, '(a)') ring%title
      bp_com%i_line = bp_com%i_line + 1
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
          call get_attribute (def$, beam_ele, &
                                 in_ring, pring, delim, delim_found, err_flag)
        endif
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
        call get_next_word (ring%lattice, ix_word, ',', &
                                                   delim, delim_found, .true.)
      endif
      cycle parsing_loop
    endif

! RETURN or END_FILE command

    if (word_1(:ix_word) == 'RETURN' .or.  &
                                    word_1(:ix_word) == 'END_FILE') then
      call file_stack ('pop', ' ', finished)
      if (.not. bmad_status%ok) return
      if (finished) exit ! break loop
      cycle parsing_loop
    endif

! variable definition or element redef.
! Note: "var := num" is old-style variable definition syntax.

    matched_delim = .false.
    if (delim == ':' .and. bp_com%parse_line(1:1) == '=') then  ! old style
      matched_delim = .true.
      bp_com%parse_line = bp_com%parse_line(2:)      ! trim off "="
    elseif (delim == '=') then
      matched_delim = .true.
    endif

! if an element attribute redef

    found = .false.
    ix = index(word_1, '[')

    if (matched_delim .and. ix /= 0) then
      name = word_1(:ix-1)  
      do i = 0, n_max+1

        if (i == n_max+1) then
          ele => param_ele
        else
          ele => in_ring%ele_(i)
        endif

        if (ele%name == name) then
          ix = index(word_1, '[')
          this_name = word_1(ix+1:)    ! name of attribute
          ix = index(this_name, ']')
          this_name = this_name(:ix-1)
          bp_com%parse_line = trim(this_name) // ' = ' // bp_com%parse_line 
          if (found) then   ! if not first time
            bp_com%parse_line = parse_line_save
          else
            parse_line_save = bp_com%parse_line
          endif
          call get_attribute (redef$, ele, in_ring, pring, &
                                             delim, delim_found, err_flag)
          if (delim_found) call warning ('BAD DELIMITER: ' // delim)
          found = .true.
        endif

      enddo

      if (.not. found) call warning ('ELEMENT NOT FOUND: ' // name)
      cycle parsing_loop

! else must be a variable

    elseif (matched_delim) then

      found = .false.
      do i = 1, bp_com%ivar_tot-1
        if (word_1 == var_(i)%name) then
          ivar = i
          found = .true.
        endif
      enddo

      if (.not. found) then
        bp_com%ivar_tot = bp_com%ivar_tot + 1
        ivar = bp_com%ivar_tot
        if (bp_com%ivar_tot > ivar_maxx) then
          print *, 'ERROR IN BMAD_PARSER: NEED TO INCREASE IVAR_MAXX!'
          call err_exit
        endif
      endif

      var_(ivar)%name = word_1

      call evaluate_value (var_(ivar)%name, var_(ivar)%value, &
                                       in_ring, delim, delim_found, err_flag)
      if (delim /= ' ' .and. .not. err_flag) call warning  &
                    ('EXTRA CHARACTERS ON RHS: ' // bp_com%parse_line,  &
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
      allocate (this_seq(bp_com%iseq_tot+1)%arg(n_arg))
      this_seq(bp_com%iseq_tot+1)%arg(:)%dummy_name = dummy_name(1:n_arg)
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
      call warning ('ARGUMENTS "XYZ(...):" ARE ONLY USED WITH REPLACEMENT LINES.', &
                                                        'FOR: ' // word_1)
      cycle parsing_loop
    endif

! if line or list

    if (word_2(:ix_word) == 'LINE' .or. word_2(:ix_word) == 'LIST') then
      bp_com%iseq_tot = bp_com%iseq_tot + 1
      if (bp_com%iseq_tot > size(this_seq)-1) then
        print *, 'ERROR IN BMAD_PARSER: NEED TO INCREASE LINE ARRAY!'
        call err_exit
      endif
      this_seq(bp_com%iseq_tot)%name = word_1

      if (delim /= '=') call warning ('EXPECTING: "=" BUT GOT: ' // delim)
      if (word_2(:ix_word) == 'LINE') then
        this_seq(bp_com%iseq_tot)%type = line$
        if (arg_list_found) this_seq(bp_com%iseq_tot)%type = replacement_line$
      else
        this_seq(bp_com%iseq_tot)%type = list$
      endif
      call seq_expand1 (this_seq, bp_com%iseq_tot, .true., in_ring)

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

      found = .false.  ! found a match?

      do i = 1, n_key
        if (word_2(:ix_word) == key_name(i)(:ix_word)) then
          in_ring%ele_(n_max)%key = i
          call preparse_element_init (in_ring%ele_(n_max))
          found = .true.
          exit
        endif
      enddo

! check if element is part of a element class
! if none of the above then we have an error

      if (.not. found) then
        do i = 1, n_max-1
          if (word_2 == in_ring%ele_(i)%name) then
            in_ring%ele_(n_max) = in_ring%ele_(i)
            in_ring%ele_(n_max)%name = word_1
            found = .true.
            exit
          endif
        enddo
      endif

      if (.not. found) then
        call warning ('KEY NAME NOT RECOGNIZED: ' // word_2,  &
                       'FOR ELEMENT: ' // in_ring%ele_(n_max)%name)
        cycle parsing_loop
      endif

! Element definition.
! we need to get the attribute values for the element.
! For control elements RING.ELE_()%IXX temporarily points to
! the PRING%ELE() array where storage for the control lists is

      key = in_ring%ele_(n_max)%key

      if (key == overlay$ .or. key == group$) then
        if (delim /= '=') then
          call warning ('EXPECTING: "=" BUT GOT: ' // delim,  &
                      'FOR ELEMENT: ' // in_ring%ele_(n_max)%name)
          cycle parsing_loop        
        endif

        in_ring%ele_(n_max)%control_type = key
        call get_overlay_group_names(in_ring%ele_(n_max), in_ring, &
                                                    pring, delim, delim_found)

        if (.not. delim_found) then
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

!---------------------------------------------------------------------------
! we now have read in everything. 

! sort elements and lists and check for duplicates

  call indexx (this_seq(1:bp_com%iseq_tot)%name, seq_indexx(1:bp_com%iseq_tot))
  call indexx (in_ring%ele_(1:n_max)%name, in_indexx(1:n_max))

  do i = 1, bp_com%iseq_tot-1
    ix1 = seq_indexx(i)
    ix2 = seq_indexx(i+1)
    if (this_seq(ix1)%name == this_seq(ix2)%name) call warning  &
                      ('DUPLICATE LINE NAME ' // this_seq(ix1)%name)
  enddo

  do i = 1, n_max-1
    ix1 = in_indexx(i)
    ix2 = in_indexx(i+1)
    if (in_ring%ele_(ix1)%name == in_ring%ele_(ix2)%name) call warning &
                    ('DUPLICATE ELEMENT NAME ' // in_ring%ele_(ix1)%name)
  enddo

  i = 1; j = 1
  do
    if (i > bp_com%iseq_tot) exit
    if (j > n_max) exit
    ix1 = seq_indexx(i)
    ix2 = in_indexx(j)
    if (this_seq(ix1)%name == in_ring%ele_(ix2)%name) call warning  &
          ('LINE AND ELEMENT HAVE THE SAME NAME: ' // this_seq(i)%name)
    if (this_seq(ix1)%name < in_ring%ele_(ix2)%name) then
      i = i + 1
    else
      j = j + 1
    endif
  enddo

! find line corresponding to the "use" statement.

  if (ring%name == blank) call error_exit &
            ('NO "USE" STATEMENT FOUND.', 'I DO NOT KNOW WHAT LINE TO USE!')

  call find_indexx (ring%name, this_seq(:)%name, &
                                    seq_indexx, bp_com%iseq_tot, i_use)
  if (i_use == 0) call error_exit &
      ('CANNOT FIND DEFINITION OF LINE IN "USE" STATEMENT: ' // ring%name, ' ')

  if (this_seq(i_use)%type /= line$) call error_exit  &
                      ('NAME IN "USE" STATEMENT IS NOT A LINE!', ' ')

! Now to expand the lines and lists to find the elements to use.
! First go through the lines and lists and index everything.

  do k = 1, bp_com%iseq_tot
    do i = 1, size(this_seq(k)%ele(:))

      s_ele => this_seq(k)%ele(i)

      ix = index(s_ele%name, '\')   ! '
      if (ix /= 0) then
        name = s_ele%name(:ix-1)
      else
        name = s_ele%name
      endif
  
      if (s_ele%ix_arg > 0) cycle  ! dummy arg

      call find_indexx (name, in_ring%ele_(1:)%name, in_indexx, n_max, j)
      if (j == 0) then  ! if not an element it must be a sequence
        call find_indexx (name, this_seq(:)%name, seq_indexx, bp_com%iseq_tot, j)
        if (j == 0) then  ! if not a sequence then I don't know what it is
          call warning ('NOT A DEFINED ELEMENT, LINE, OR LIST: ' // &
                 s_ele%name, 'IN THE LINE/LIST: ' // this_seq(k)%name, &
                 seq = this_seq(k))
          s_ele%ix_array = 1
          s_ele%type = element$
        else
          s_ele%ix_array = j
          s_ele%type = this_seq(j)%type
          if (this_seq(k)%type == list$) call warning ('LIST: ' // & 
                  this_seq(k)%name, 'CONTAINS A NON-ELEMENT: ' // s_ele%name, &
                  seq = this_seq(k))
        endif
        if (s_ele%type == list$ .and. s_ele%reflect) call warning ( &
                          'A REFLECTION WITH A LIST IS NOT ALLOWED IN: '  &
                          // this_seq(k)%name, 'FOR LIST: ' // s_ele%name, &
                          seq = this_seq(k))

      else    ! if an element...
        s_ele%ix_array = j
        s_ele%type = element$
      endif

    enddo
  enddo

! to expand the ring we use a stack for nested sublines.
! IX_RING is the expanded array of elements in the ring.
! init stack

  i_lev = 1                          ! level on the stack
  stack(1)%ix_seq  = i_use           ! which sequence to use for the ring
  stack(1)%ix_ele  =  1              ! we start at the beginning
  stack(1)%reflect = +1              ! and move forward

  seq => this_seq(i_use)

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
      seq2 => this_seq(s_ele%ix_array)
      j = seq2%ix
      call pushit (ix_ring, n_ele_ring, seq2%ele(j)%ix_array)
      name_(n_ele_ring) = seq2%ele(j)%name
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
      seq => this_seq(s_ele%ix_array)
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
        seq => this_seq(stack(i_lev)%ix_seq)
      else
        exit
      endif
    enddo

  enddo line_expansion

!---------------------------------------------------------------
! we now have the line to use in constructing the ring.
! now to put the elements in RING in the correct order.
! superimpose, overlays, and groups are handled later.
! first load beam parameters.

  ring%version            = bmad_inc_version$
  ring%input_file_name    = full_name             ! save input file
  ring%param%particle     = nint(beam_ele%value(particle$))
  ring%param%n_part       = beam_ele%value(n_part$)
  ring%n_ele_ring         = n_ele_ring
  ring%n_ele_use          = n_ele_ring
  ring%n_ele_max          = n_ele_ring
  ring%param%aperture_limit_on  = .true.
  ring%param%damping_on         = .false.
  ring%param%fluctuations_on    = .false.
  ring%n_ele_symm         = 0                     ! no symmetry point
  ring%n_ic_array         = 0                     
  ring%n_control_array    = 0    
  call set_symmetry (ring%param%symmetry, ring)   ! set ring.n_ele_use

  ring%ele_(0) = in_ring%ele_(0)    ! Beginning element
  call init_ele (ring%ele_init)

! New way of doing things

  ring%param%lattice_type = nint(param_ele%value(lattice_type$))
  ring%param%symmetry = nint(param_ele%value(symmetry$))

  if (nint(param_ele%value(taylor_order$)) /= 0) &
            ring%input_taylor_order = nint(param_ele%value(taylor_order$))

! old way of doing things

  do i = 1, bp_com%ivar_tot
    if (var_(i)%name == 'SYMMETRY') ring%param%symmetry = nint(var_(i)%value)
    if (var_(i)%name == 'LATTICE_TYPE')  &
                              ring%param%lattice_type = nint(var_(i)%value)
    if (var_(i)%name == 'TAYLOR_ORDER') &
                              ring%input_taylor_order = nint(var_(i)%value)
  enddo

!

  ring%param%beam_energy  = param_ele%value(beam_energy$)
  ring%param%energy       = beam_ele%value(energy$)

  if (ring%param%beam_energy /= 0 .and. ring%param%energy /= 0) call warning &
         ('BOTH "BEAM, ENERGY" AND "PARAMETER[BEAM_ENERGY]" CONSTRUCTS USED!')

  if (ring%param%beam_energy /= 0) then
    ring%param%energy  = 1d-9 * ring%param%beam_energy
  else
    ring%param%beam_energy  = 1d9 * ring%param%energy
  endif

!

  if (ring%input_taylor_order /= 0) &
                  call set_taylor_order (ring%input_taylor_order, .false.)

  if (n_ele_ring > n_ele_maxx) then
    print *, 'ERROR IN BMAD_PARSER: NUMBER OF ELEMENTS EXCEEDS ELEMENT ARRAY'
    print *, '       SIZE THE RING STRUCTURE:', n_ele_ring, n_ele_maxx
    call err_exit
  endif

  if (any(in_ring%ele_(:)%key == lcavity$) .and. &
                          ring%param%lattice_type /= linac_lattice$) then
    print *, 'ERROR IN BMAD_PARSER: THERE IS A LCAVITY BUT THE'
    print *, '      LATTICE_TYPE IS NOT SET TO LINAC_LATTICE!'
    call err_exit
  endif

! transfer the ele information from the in_ring to ring

  do i = 1, n_ele_ring

! Use the name as given in sequence lists since elements can have different 
! names from the defining elements. Eg: B01\H2 gets its definition from B01.

    ele => ring%ele_(i)
    ele = in_ring%ele_(ix_ring(i)) 
    ele%name = name_(i)

! Convert rbends to sbends and evaluate G if needed.
! Needed is the length and either: angle, G, or rho.

    if (ele%key == sbend$ .or. ele%key == rbend$) then

      angle = ele%value(angle$) 

      if (ele%value(g$) /= 0 .and. ele%value(rho$) /= 0) &
            call warning ('BOTH G AND RHO SPECIFIED FOR BEND: ' // ele%name)
      if (ele%value(rho$) /= 0) ele%value(g$) = 1 / ele%value(rho$)

      if (ele%value(g$) /= 0 .and. angle /= 0) call warning &
            ('BOTH ANGLE AND G OR RHO SPECIFIED FOR BEND: ' // ele%name)

      if (ele%key == rbend$) then
        ele%value(l_chord$) = ele%value(l$)
        
        if (angle /= 0) then
          ele%value(l$) = ele%value(l_chord$) * angle / (2 * sin(angle/2))
        elseif (ele%value(g$) /= 0 .and. ele%value(l_chord$) == 0) then
          angle = ele%value(g$) * ele%value(l_chord$) / 2
          ele%value(l$) = ele%value(l_chord$) * asin(angle)/ angle
        endif
        ele%value(e1$) = ele%value(e1$) + angle / 2
        ele%value(e2$) = ele%value(e2$) + angle / 2
        ele%key = sbend$
      endif

      if (ele%value(angle$) /= 0) ele%value(g$) = &
                                        ele%value(angle$) / ele%value(l$) 

      if (ele%value(hgapx$) == 0) ele%value(hgapx$) = ele%value(hgap$)
      if (ele%value(fintx$) == 0) ele%value(fintx$) = ele%value(fint$)

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

! Now put in the overlay_lord and group elements

  do i = 1, n_max
    ct = in_ring%ele_(i)%control_type
    if (ct /= overlay_lord$ .and. ct /= group_lord$) cycle
    ring%n_ele_max = ring%n_ele_max + 1
    ixs = ring%n_ele_max
    ring%ele_(ixs) = in_ring%ele_(i)
    call find_slaves_for_parser (ring, pring%ele(i)%name_, &
                   pring%ele(i)%attrib_name_, pring%ele(i)%coef_, cs_)
    ele => ring%ele_(ixs)
    if (ct == overlay_lord$) then
      call create_overlay (ring, ixs, ele%ix_value, ele%n_slave, cs_)
    else
      call create_group (ring, ixs, ele%n_slave, cs_)
    endif
  enddo

! make matrices for entire ring

  do i = ring%n_ele_ring+1, ring%n_ele_max
    call control_bookkeeper (ring, i)      ! need this for ring_geometry
  enddo

  call s_calc (ring)                       ! calc longitudinal distances
  call ring_geometry (ring)                ! ring layout
  call set_ptc (ring%param)

! Reuse the old taylor series if they exist
! and the old taylor series has the same attributes.

  if (associated(old_ele)) then

    do i = 1, ring%n_ele_max

      if (old_energy /= ring%param%beam_energy) exit
      
      ele => ring%ele_(i)
      call attribute_bookkeeper (ele, ring%param) ! so value arrays are equal

      vmask = .true.
      if (ele%key == wiggler$) vmask((/k1$, rho$, b_max$/)) = .false.
      do j = 1, size(old_ele)
        if (old_ele(j)%key /= ele%key) cycle
        if (old_ele(j)%name /= ele%name) cycle
        if (any(old_ele(j)%value /= ele%value .and. vmask)) cycle
        if (old_ele(j)%num_steps /= ele%num_steps) cycle
        if (old_ele(j)%integration_order /= ele%integration_order) cycle
        if (associated(old_ele(j)%wig_term) .and. &
                                               associated(ele%wig_term)) then
          if (size(old_ele(j)%wig_term) /= size(ele%wig_term)) cycle
          do it = 1, size(old_ele(j)%wig_term)
            if (old_ele(j)%wig_term(it)%coef /= ele%wig_term(it)%coef) cycle
            if (old_ele(j)%wig_term(it)%kx /= ele%wig_term(it)%kx) cycle
            if (old_ele(j)%wig_term(it)%ky /= ele%wig_term(it)%ky) cycle
            if (old_ele(j)%wig_term(it)%kz /= ele%wig_term(it)%kz) cycle
            if (old_ele(j)%wig_term(it)%phi_z /= ele%wig_term(it)%phi_z) cycle
          enddo
        elseif (associated(old_ele(j)%wig_term) .xor. &
                                            associated(ele%wig_term)) then
          cycle
        endif
        if (any(old_ele(j)%taylor(:)%ref /= 0)) cycle
        if (bmad_com%taylor_order > old_ele(j)%taylor_order) cycle
        exit
      enddo

      if (j == size(old_ele) + 1) cycle

      print *, 'BMAD_PARSER: Reusing Taylor for: ', old_ele(j)%name

      do it = 1, 6
        ix = 0
        do k = 1, size(old_ele(j)%taylor(it)%term)
         if (sum(old_ele(j)%taylor(it)%term(k)%exp(:)) <= &
                                      bmad_com%taylor_order) ix = ix + 1
        enddo
        allocate (ele%taylor(it)%term(ix))
        ix = 0
        do k = 1, size(old_ele(j)%taylor(it)%term)
         if (sum(old_ele(j)%taylor(it)%term(k)%exp(:)) <= &
                                      bmad_com%taylor_order) then
            ix = ix + 1
            ele%taylor(it)%term(ix) = old_ele(j)%taylor(it)%term(k)
          endif      
        enddo
      enddo
      ele%taylor_order = bmad_com%taylor_order
      ele%taylor(:)%ref = 0
    enddo

    deallocate(old_ele)

  endif

!

  if (ring%param%beam_energy == 0 .and. ring%ele_(0)%value(energy$) == 0) &
                    print *, 'WARNING FROM BMAD_PARSER: BEAM_ENERGY IS 0!'

  doit = .true.
  if (present(make_mats6)) doit = make_mats6
  call compute_element_energy (ring)
  if (doit) call ring_make_mat6(ring, -1)      ! make 6x6 transport matrices

!-------------------------------------------------------------------------
! write out if debug is on

  if (bp_com%parser_debug) then
    
    if (index(bp_com%debug_line, 'VAR') /= 0) then
      print *
      print *, '----------------------------------------'
      print *, 'Number of Variables:', bp_com%ivar_tot - bp_com%ivar_init
      do i = bp_com%ivar_init+1, bp_com%ivar_tot
        print *
        print *, 'Var #', i-bp_com%ivar_tot
        print *, 'Name: ', var_(i)%name
        print *, 'Value:', var_(i)%value
      enddo
    endif

    if (index(bp_com%debug_line, 'SEQ') /= 0) then
      print *
      print *, '----------------------------------------'
      print *, 'Number of Lines/Lists defined:', bp_com%iseq_tot
      do i = 1, bp_com%iseq_tot
        print *
        print *, 'Sequence #', i
        print *, 'Name: ', this_seq(i)%name
        print *, 'Type:', this_seq(i)%type
        print *, 'Number of elements:', size(this_seq(i)%ele)
        print *, 'List:'
        do j = 1, size(this_seq(i)%ele(:))
          print '(4x, a, l3, i3)', this_seq(i)%ele(j)%name, &
              this_seq(i)%ele(j)%reflect, this_seq(i)%ele(j)%ix_arg
        enddo
      enddo
    endif

    if (index(bp_com%debug_line, 'SLAVE') /= 0) then
      print *
      print *, '----------------------------------------'
      print *, 'Number of Elements in Regular Ring:', ring%n_ele_ring
      do i = 1, ring%n_ele_ring
        print *, '-------------'
        print *, 'Ele #', i
        call type_ele (ring%ele_(i), .false., 0, .false., 0, .true., ring)
      enddo
    endif

    if (index(bp_com%debug_line, 'LORD') /= 0) then
      print *
      print *, '----------------------------------------'
      print *, 'LORD elements: ', ring%n_ele_max - ring%n_ele_ring
      do i = ring%n_ele_ring+1, ring%n_ele_max
        print *, '-------------'
        print *, 'Ele #', i
        call type_ele (ring%ele_(i), .false., 0, .false., 0, .true., ring)
      enddo
    endif

    if (index(bp_com%debug_line, 'RING') /= 0) then  
      print *
      print *, '----------------------------------------'
      print *, 'Ring Used: ', ring%name
      print *, 'Number of ring elements:', ring%n_ele_ring
      print *, 'List:                               Key      Length         S'
      do i = 1, ring%n_ele_ring
        print '(3x, i3, 2a, 3x, a, 2f10.2)', i, ') ', ring%ele_(i)%name,  &
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
      deallocate(pring%ele(i)%coef_)
    endif
  enddo

  do i = 1, size(this_seq(:))
    if (associated (this_seq(i)%arg)) deallocate(this_seq(i)%arg)
    if (associated (this_seq(i)%ele)) deallocate(this_seq(i)%ele)
  enddo

! error check

  if (bp_com%error_flag) then
    if (bmad_status%exit_on_error) then
      print *, 'BMAD_PARSER FINISHED. EXITING ON ERRORS'
      stop
    else
      bmad_status%ok = .false.
      return
    endif
  endif
              
  call check_ring_controls (ring, .true.)

! write to digested file

  if (bp_com%write_digested .and. .not. bp_com%parser_debug .and. &
      digested_version <= bmad_inc_version$) call write_digested_bmad_file  &
               (digested_file, ring, bp_com%n_files, bp_com%file_name_)

end subroutine
