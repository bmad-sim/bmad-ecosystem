!+
! Subroutine bmad_parser (in_file, ring, make_mats6, digested_read_ok)
!
! Subroutine to parse a BMAD input file and put the information in ring.
!
! Because of the time it takes to parse a file BMAD_PARSER will save 
! RING in a "digested" file with the name:
!               'digested_' // IN_FILE   ! for single precision BMAD version
!               'digested8_' // IN_FILE  ! for double precision BMAD version
! For subsequent calls to the same IN_FILE, BMAD_PARSER will just read in the
! digested file. BMAD_PARSER will always check to see that the digested file
! is up-to-date and if not the digested file will not be used.
!
! Modules needed:
!   use bmad
!
! Input:
!   in_file    -- Character: Name of the input file.
!   make_mats6 -- Logical, optional: Compute the 6x6 transport matrices for the
!                   Elements? Default is True.
!
! Output:
!   ring -- Ring_struct: Ring structure. See bmad_struct.f90 for more details.
!     %ele_(:)%mat6  -- This is computed assuming an on-axis orbit 
!     %ele_(:)%s     -- This is also computed.
!   digested_read_ok -- Logical, optional: Set True if the digested file was
!                        successfully read. False otherwise.
!   bmad_status      -- Bmad status common block.
!     %ok              -- Set True if parsing is successful. False otherwise.
!         
! Defaults:
!   ring%param%particle          = positron$
!   ring%param%lattice_type      = circular_lattice$
!   ring%param%aperture_limit_on = .true.
!-

#include "CESR_platform.inc"

subroutine bmad_parser (in_file, ring, make_mats6, digested_read_ok)

  use bmad_parser_mod
  use cesr_utils
  
  implicit none

  type (ring_struct), target :: ring, in_ring
  type (ele_struct) this_ele
  type (seq_struct), save, target :: sequence_(1000)
  type (seq_struct), pointer :: seq, seq2
  type (seq_ele_struct), target :: this_s_ele
  type (seq_ele_struct), pointer :: s_ele
  type (parser_ring_struct) pring
  type (seq_stack_struct) stack(20)
  type (ele_struct), save, pointer :: ele, old_ele(:) => null()

  integer, allocatable :: ix_ring(:)
  integer, allocatable :: seq_indexx(:), in_indexx(:)
  character(16), allocatable ::  name_(:), in_name(:), seq_name(:)

  integer ix_word, i_use, i, j, k, n, ix, ixr, ixs, i_ring, it
  integer i_lev, i_key, ic, ix_lord
  integer iseq, n_ele_ring, temp_maxx
  integer ix_super, digested_version, key, ct
  integer ix1, ix2, iv, n_arg, iseq_tot
  integer ivar, ixx, j_lord, n_slave
  integer, pointer :: n_max

  character*(*) in_file
  character(16) word_2, name, a_name
  character(1) delim*1
  character(200) path, basename, full_name, digested_file
  character(40) this_name, word_1
  character(280) parse_line_save

  real(rp) angle, old_energy

  logical, optional :: make_mats6, digested_read_ok
  logical parsing, delim_found, matched_delim, arg_list_found, doit
  logical file_end, found, err_flag, finished, save_taylor
  logical, save :: init_needed = .true.

! see if digested file is open and current. If so read in and return.
! Note: The name of the digested file depends upon the real precision.

  bp_com%parser_name = 'BMAD_PARSER'  ! Used for error messages.
  bp_com%write_digested = .true.
  bp_com%input_line_meaningful = .true.

  inquire (file = in_file, name = full_name)      ! full input file_name
  ix = index(full_name, ';')
  if (ix /= 0) full_name = full_name(:ix-1)
  ring%input_file_name = full_name      ! needed by read_digested_bmad_file

  ix = SplitFileName(in_file, path, basename)
  if (rp == 8) then
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
      if (bmad_status%type_out) then
        print *, 'BMAD_PARSER: Taylor_order has changed.'
        print *, '    Taylor_order in digested file:', ring%input_taylor_order
        print *, '    Taylor_order now:             ', bmad_com%taylor_order
      endif
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

! here if not OK bmad_status. So we have to do everything from scratch...
! init variables.

  nullify (pring%ele)
  call init_ring (in_ring, 1000)
  call allocate_pring (in_ring, pring)

  bmad_status%ok = .true.
  if (bmad_status%type_out) &
                        print *, 'BMAD_PARSER: Creating new digested file...'

  bp_com%n_files = 0
  bp_com%error_flag = .false.                 ! set to true on an error
  call file_stack('init', in_file, finished)   ! init stack
  call file_stack('push', in_file, finished)   ! open file on stack
  if (.not. bmad_status%ok) return
  bp_com%parser_debug = .false.
  iseq_tot = 0                            ! number of sequences encountered
  bp_com%ivar_tot = 0                     ! number of variables encountered
  call init_bmad_parser_common

  ring%name = ' '
  ring%lattice = ' '

  call init_ele (in_ring%ele_(0))
  in_ring%ele_(0)%name = 'BEGINNING'     ! Beginning element
  in_ring%ele_(0)%key = init_ele$

  call mat_make_unit (in_ring%ele_(0)%mat6)
  in_ring%param%t1_mat6 = 0
  in_ring%param%t1_mat4 = 0

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

  ring%n_control_max = 0

!-----------------------------------------------------------
! main parsing loop

  parsing_loop: do 

! get a line from the input file and parse out the first word

    call load_parse_line ('normal', 1, file_end)  ! load an input line
    call get_next_word (word_1, ix_word, ':(,)= ', delim, delim_found, .true.)
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
      call str_upcase (bp_com%debug_line, bp_com%debug_line)
      print *, 'FOUND IN FILE: "PARSER_DEBUG". DEBUG IS NOW ON'
      cycle parsing_loop
    endif

! NO_DIGESTED

    if (word_1(:ix_word) == 'NO_DIGESTED') then
      bp_com%write_digested = .false.
      print *, 'FOUND IN FILE: "NO_DIGESTED". NO DIGESTED FILE WILL BE CREATED'
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
      ring%name = word_2
      cycle parsing_loop
    endif

! TITLE command

    if (word_1(:ix_word) == 'TITLE') then
      if (delim_found) then
        if (delim /= " " .and. delim /= ",") call warning &
                            ('BAD DELIMITOR IN "TITLE" COMMAND')
        call type_get (this_ele, descrip$, delim, delim_found)
        ring%title = this_ele%descrip
        deallocate (this_ele%descrip)
      else
        read (bp_com%current_file%f_unit, '(a)') ring%title
        bp_com%current_file%i_line = bp_com%current_file%i_line + 1
      endif
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

      call parser_add_variable (word_1, in_ring)
      cycle parsing_loop

    endif

! if a "(" delimitor then we are looking at a replacement line.

    if (delim == '(') then
      call get_sequence_args (word_1, sequence_(iseq_tot+1)%arg, &
                                                       delim, err_flag)
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
      if (iseq_tot > size(sequence_)-1) then
        print *, 'ERROR IN BMAD_PARSER: NEED TO INCREASE LINE ARRAY!'
        call err_exit
      endif
      sequence_(iseq_tot)%name = word_1

      if (delim /= '=') call warning ('EXPECTING: "=" BUT GOT: ' // delim)
      if (word_2(:ix_word) == 'LINE') then
        sequence_(iseq_tot)%type = line$
        if (arg_list_found) sequence_(iseq_tot)%type = replacement_line$
      else
        sequence_(iseq_tot)%type = list$
      endif
      call seq_expand1 (sequence_, iseq_tot, in_ring, .true.)

! if not line or list then must be an element

    else

      n_max = n_max + 1
      if (n_max > in_ring%n_ele_maxx) then
        call allocate_ring_ele_ (in_ring)
        call allocate_pring (in_ring, pring)
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
! For control elements RING%ELE_()%IXX temporarily points to
! the PRING%ELE() array where storage for the control lists is

      key = in_ring%ele_(n_max)%key

      if (key == custom$) then
        in_ring%ele_(n_max)%mat6_calc_method = custom$
        in_ring%ele_(n_max)%tracking_method  = custom$
        in_ring%ele_(n_max)%field_calc       = custom$
      endif

      if (key == taylor$) then   ! start with unit matrix
        call add_taylor_term (in_ring%ele_(n_max), 1, 1.0_rp, (/ 1, 0, 0, 0, 0, 0 /)) 
        call add_taylor_term (in_ring%ele_(n_max), 2, 1.0_rp, (/ 0, 1, 0, 0, 0, 0 /)) 
        call add_taylor_term (in_ring%ele_(n_max), 3, 1.0_rp, (/ 0, 0, 1, 0, 0, 0 /)) 
        call add_taylor_term (in_ring%ele_(n_max), 4, 1.0_rp, (/ 0, 0, 0, 1, 0, 0 /)) 
        call add_taylor_term (in_ring%ele_(n_max), 5, 1.0_rp, (/ 0, 0, 0, 0, 1, 0 /)) 
        call add_taylor_term (in_ring%ele_(n_max), 6, 1.0_rp, (/ 0, 0, 0, 0, 0, 1 /)) 
      endif

      if (key == overlay$ .or. key == group$ .or. key == i_beam$) then
        if (delim /= '=') then
          call warning ('EXPECTING: "=" BUT GOT: ' // delim,  &
                      'FOR ELEMENT: ' // in_ring%ele_(n_max)%name)
          cycle parsing_loop        
        endif

        if (key == overlay$) in_ring%ele_(n_max)%control_type = overlay_lord$
        if (key == group$)   in_ring%ele_(n_max)%control_type = group_lord$
        if (key == i_beam$)  in_ring%ele_(n_max)%control_type = i_beam_lord$

        call get_overlay_group_names(in_ring%ele_(n_max), in_ring, &
                                                    pring, delim, delim_found)

        if (key /= i_beam$ .and. .not. delim_found) then
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

  bp_com%input_line_meaningful = .false.

! sort elements and lists and check for duplicates
! seq_name(:) and in_name(:) arrays speed up the calls to find_indexx since
! the compiler does not have to repack the memory.

  allocate (seq_indexx(iseq_tot), seq_name(iseq_tot))
  seq_name = sequence_(1:iseq_tot)%name
  call indexx (seq_name, seq_indexx(1:iseq_tot))

  allocate (in_indexx(n_max), in_name(n_max))
  in_name = in_ring%ele_(1:n_max)%name
  call indexx (in_name, in_indexx(1:n_max))

  do i = 1, iseq_tot-1
    ix1 = seq_indexx(i)
    ix2 = seq_indexx(i+1)
    if (sequence_(ix1)%name == sequence_(ix2)%name) call warning  &
                      ('DUPLICATE LINE NAME ' // sequence_(ix1)%name)
  enddo

  do i = 1, n_max-1
    ix1 = in_indexx(i)
    ix2 = in_indexx(i+1)
    if (in_ring%ele_(ix1)%name == in_ring%ele_(ix2)%name) call warning &
                    ('DUPLICATE ELEMENT NAME ' // in_ring%ele_(ix1)%name)
  enddo

  i = 1; j = 1
  do
    if (i > iseq_tot) exit
    if (j > n_max) exit
    ix1 = seq_indexx(i)
    ix2 = in_indexx(j)
    if (sequence_(ix1)%name == in_ring%ele_(ix2)%name) call warning  &
          ('LINE AND ELEMENT HAVE THE SAME NAME: ' // sequence_(ix1)%name)
    if (sequence_(ix1)%name < in_ring%ele_(ix2)%name) then
      i = i + 1
    else
      j = j + 1
    endif
  enddo

! find line corresponding to the "use" statement.

  if (ring%name == blank) call error_exit &
            ('NO "USE" STATEMENT FOUND.', 'I DO NOT KNOW WHAT LINE TO USE!')

  call find_indexx (ring%name, seq_name, seq_indexx, iseq_tot, i_use)
  if (i_use == 0) call error_exit &
      ('CANNOT FIND DEFINITION OF LINE IN "USE" STATEMENT: ' // ring%name, ' ')

  if (sequence_(i_use)%type /= line$) call error_exit  &
                      ('NAME IN "USE" STATEMENT IS NOT A LINE!', ' ')

! Now to expand the lines and lists to find the elements to use.
! First go through the lines and lists and index everything.

  do k = 1, iseq_tot
    do i = 1, size(sequence_(k)%ele(:))

      s_ele => sequence_(k)%ele(i)

      ix = index(s_ele%name, '\')   ! '
      if (ix /= 0) then
        name = s_ele%name(:ix-1)
      else
        name = s_ele%name
      endif
  
      if (s_ele%ix_arg > 0) then   ! dummy arg
        s_ele%type = element$
        cycle
      endif

      call find_indexx (name, in_name, in_indexx, n_max, j)
      if (j == 0) then  ! if not an element it must be a sequence
        call find_indexx (name, seq_name, seq_indexx, iseq_tot, j)
        if (j == 0) then  ! if not a sequence then I don't know what it is
          s_ele%ix_array = -1
          s_ele%type = element$
        else
          s_ele%ix_array = j
          s_ele%type = sequence_(j)%type
        endif
        if (s_ele%type == list$ .and. s_ele%reflect) call warning ( &
                          'A REFLECTION WITH A LIST IS NOT ALLOWED IN: '  &
                          // sequence_(k)%name, 'FOR LIST: ' // s_ele%name, &
                          seq = sequence_(k))
        if (sequence_(k)%type == list$) &
                call warning ('A REPLACEMENT LIST: ' // sequence_(k)%name, &
                'HAS A NON-ELEMENT MEMBER: ' // s_ele%name)
 
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
  seq => sequence_(i_use)

  stack(1)%ix_seq  = i_use           ! which sequence to use for the ring
  stack(1)%ix_ele  =  1              ! we start at the beginning
  stack(1)%reflect = +1              ! and move forward
  stack(1)%rep_count = seq%ele(1)%rep_count

  n_ele_ring = 0
           
  allocate (name_(in_ring%n_ele_maxx))
  allocate (ix_ring(in_ring%n_ele_maxx))
  ix_ring = -1
  sequence_(:)%ix = 1  ! Init. Used for replacement list index

! expand line

  parsing = .true.
  line_expansion: do while (parsing)

    s_ele => seq%ele(stack(i_lev)%ix_ele)

    ix = s_ele%ix_arg
    if (ix /= 0) then  ! it is a dummy argument.
      name = seq%arg(ix)%actual_name
      s_ele => this_s_ele
      s_ele%name = name
      call find_indexx (name, in_name, in_indexx, n_max, j)
      if (j == 0) then  ! if not an element it must be a sequence
        call find_indexx (name, seq_name, seq_indexx, iseq_tot, j)
        if (j == 0) then  ! if not a sequence then I don't know what it is
          call warning ('CANNOT FIND DEFINITION FOR: ' // name, &
                          'IN LINE: ' // seq%name, seq)
          call err_exit
        endif
        s_ele%ix_array = j
        s_ele%type = sequence_(j)%type
      else
        s_ele%ix_array = j 
        s_ele%type = element$
      endif
      
    endif

    stack(i_lev)%rep_count = stack(i_lev)%rep_count - 1
    if (stack(i_lev)%rep_count == 0) then
      call increment_pointer (stack(i_lev)%ix_ele, stack(i_lev)%reflect)
      ix = stack(i_lev)%ix_ele
      if (ix > 0 .and. ix <= size(seq%ele)) &
                            stack(i_lev)%rep_count = seq%ele(ix)%rep_count
    endif

! if an element

    select case (s_ele%type)
    case (element$) 
      if (n_ele_ring+1 > size(ix_ring)) then
        n = 1.5*n_ele_ring
        call reallocate_integer (ix_ring, n)
        call reallocate_string (name_, len(name_(1)), n)
      endif
      call pushit (ix_ring, n_ele_ring, s_ele%ix_array)
      name_(n_ele_ring) = s_ele%name

      if (s_ele%ix_array < 1) call warning('NOT A DEFINED ELEMENT: ' // &
                          s_ele%name, 'IN THE LINE/LIST: ' // seq%name, seq)

! if a list

    case (list$) 
      seq2 => sequence_(s_ele%ix_array)
      j = seq2%ix
      if (n_ele_ring+1 > size(ix_ring)) then
        n = 1.5*n_ele_ring
        call reallocate_integer (ix_ring, n)
        call reallocate_string (name_, len(name_(1)), n)
      endif
      call pushit (ix_ring, n_ele_ring, seq2%ele(j)%ix_array)
      name_(n_ele_ring) = seq2%ele(j)%name
      seq2%ix = seq2%ix + 1
      if (seq2%ix > size(seq2%ele(:))) seq2%ix = 1

! if a line:
!     a) move pointer on current level past line element
!     b) go to the next higher level
!     c) initialize pointers for the higher level to use the line

    case (line$, replacement_line$)
      i_lev = i_lev + 1
      if (s_ele%type == replacement_line$) then
        seq2 => sequence_(s_ele%ix_array)
        if (size(seq2%arg) /= size(s_ele%arg)) then
          call warning ('WRONG NUMBER OF ARGUMENTS FORREPLACEMENT LINE: ' // &
              s_ele%name, 'WHEN USED IN LINE: ' // seq%name, seq)
          call err_exit
        endif
        arg_loop: do i = 1, size(seq2%arg)
          seq2%arg(i)%actual_name = s_ele%arg(i)%actual_name
          if (associated(seq%arg)) then
            do j = 1, size(seq%arg)
              if (seq2%arg(i)%actual_name == seq%arg(j)%dummy_name) then
                seq2%arg(i)%actual_name = seq%arg(j)%actual_name
                seq2%arg(i)%ix_array = seq%arg(j)%ix_array
                cycle arg_loop
              endif
            enddo
          endif
          name = seq2%arg(i)%actual_name
        enddo arg_loop
      endif
      seq => sequence_(s_ele%ix_array)
      stack(i_lev)%ix_seq = s_ele%ix_array
      stack(i_lev)%reflect = stack(i_lev-1)%reflect
      if (s_ele%reflect) stack(i_lev)%reflect = -stack(i_lev)%reflect
      if (stack(i_lev)%reflect == 1) then
        ix = 1
      else
        ix = size(seq%ele(:))
      endif
      stack(i_lev)%ix_ele = ix
      stack(i_lev)%rep_count = seq%ele(ix)%rep_count

    case default
      call warning ('INTERNAL SEQUENCE ERROR!')

    end select

! if we have got to the end of the current line then pop the stack back to
! the next lower level.
! Also check if we have gotten to level 0 which says that we are done

    do
      if (stack(i_lev)%ix_ele < 1 .or. &
                        stack(i_lev)%ix_ele > size(seq%ele(:))) then
        i_lev = i_lev - 1
        if (i_lev == 0) exit line_expansion
        seq => sequence_(stack(i_lev)%ix_seq)
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

  call init_ring (ring, n_ele_ring+100)

  ring%version            = bmad_inc_version$
  ring%input_file_name    = full_name             ! save input file
  ring%param%particle     = nint(beam_ele%value(particle$))
  ring%param%n_part       = beam_ele%value(n_part$)
  ring%param%charge       = 0
  ring%n_ele_ring         = n_ele_ring
  ring%n_ele_use          = n_ele_ring
  ring%n_ele_max          = n_ele_ring
  ring%param%aperture_limit_on  = .true.
  ring%n_ic_max           = 0                     
  ring%n_control_max      = 0    
  ring%n_ele_use          = ring%n_ele_ring

  ring%ele_(0) = in_ring%ele_(0)    ! Beginning element

! New way of doing things

  ring%param%lattice_type = nint(param_ele%value(lattice_type$))

  if (nint(param_ele%value(taylor_order$)) /= 0) &
            ring%input_taylor_order = nint(param_ele%value(taylor_order$))

! old way of doing things

  do i = 1, bp_com%ivar_tot
    if (bp_com%var_name(i) == 'LATTICE_TYPE')  &
                          ring%param%lattice_type = nint(bp_com%var_value(i))
    if (bp_com%var_name(i) == 'TAYLOR_ORDER') &
                          ring%input_taylor_order = nint(bp_com%var_value(i))
  enddo

! Beam energy

  if (param_ele%value(beam_energy$) == 0) then
    ring%param%beam_energy = 1d9 * beam_ele%value(energy_gev$)  
  elseif (beam_ele%value(energy_gev$) == 0) then
    ring%param%beam_energy = param_ele%value(beam_energy$)
  else
    call warning &
         ('BOTH "BEAM, ENERGY" AND "PARAMETER[BEAM_ENERGY]" CONSTRUCTS USED!')
  endif

! Set taylor order and lattice_type

  if (ring%input_taylor_order /= 0) &
       call set_taylor_order (ring%input_taylor_order, .false.)

  if (any(in_ring%ele_(:)%key == lcavity$) .and. &
                          ring%param%lattice_type /= linac_lattice$) then
    print *, 'Note in BMAD_PARSER: This lattice has a LCAVITY.'
    print *, '     Setting the LATTICE_TYPE to LINAC_LATTICE.'
    ring%param%lattice_type = linac_lattice$
  endif

! transfer the ele information from the in_ring to ring

  do i = 1, n_ele_ring

! Use the name as given in sequence lists since elements can have different 
! names from the defining elements. Eg: B01\H2 gets its definition from B01.

    ele => ring%ele_(i)
    ele = in_ring%ele_(ix_ring(i)) 
    ele%name = name_(i)

    select case (ele%key)

! Convert rbends to sbends and evaluate G if needed.
! Needed is the length and either: angle, G, or rho.

    case (sbend$, rbend$) 

      ele%sub_key = ele%key  ! save input format.
      angle = ele%value(angle$) 

      if (ele%value(b_field$) /= 0 .and. ele%value(g$) /= 0) call warning &
          ('BOTH G AND B_FIELD SET FOR A BEND: ' // ele%name)

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

! Accept Use of Delta_E for lcavities

    case (lcavity$) 

      if (ele%value(delta_e$) /= 0) then
        if (ele%value(gradient$) /= 0) call warning &
                ('BOTH DELTA_E AND gradient NON-ZERO FOR A LCAVITY:', ele%name)
        ele%value(gradient$) = ele%value(delta_e$) / ele%value(l$)
      endif

! check for inconsistancies

    case (solenoid$)
      if (ele%value(b_field$) /= 0 .and. ele%value(ks$) /= 0) call warning &
          ('BOTH KS AND B_FIELD SET FOR A SOLENOID: ' // ele%name)

    case (quadrupole$)
      if (ele%value(b_gradient$) /= 0 .and. ele%value(k1$) /= 0) call warning &
          ('BOTH K1 AND B_GRADIENT SET FOR A QUAD: ' // ele%name)

    case (sextupole$)
      if (ele%value(b_gradient$) /= 0 .and. ele%value(k2$) /= 0) call warning &
          ('BOTH K2 AND B_GRADIENT SET FOR A SEXTUPOLE: ' // ele%name)

    case (octupole$)
      if (ele%value(b_gradient$) /= 0 .and. ele%value(k3$) /= 0) call warning &
          ('BOTH K3 AND B_GRADIENT SET FOR A OCTUPOLE: ' // ele%name)

    end select

! aperture

    if (ele%value(aperture$) /= 0) then
      ele%value(x_limit$) = ele%value(aperture$)
      ele%value(y_limit$) = ele%value(aperture$)
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

! Now put in the overlay_lord, i_beam, and group elements

  call parser_add_lord (in_ring, 1, n_max, pring, ring)

! make matrices for entire ring

  do i = ring%n_ele_ring+1, ring%n_ele_max
    call control_bookkeeper (ring, i)      ! need this for ring_geometry
  enddo

  call s_calc (ring)                       ! calc longitudinal distances
  call ring_geometry (ring)                ! ring layout
  call set_ptc (ring%param)
  ring%input_taylor_order = bmad_com%taylor_order

! Reuse the old taylor series if they exist
! and the old taylor series has the same attributes.

  if (associated(old_ele)) then

    do i = 1, ring%n_ele_max

      if (old_energy /= ring%param%beam_energy) exit
      
      ele => ring%ele_(i)
      call attribute_bookkeeper (ele, ring%param) ! so value arrays are equal

      do j = 1, size(old_ele)
        if (any(old_ele(j)%taylor(:)%ref /= 0)) cycle
        if (bmad_com%taylor_order > old_ele(j)%taylor_order) cycle
        if (.not. equivalent_eles (old_ele(j), ele)) cycle
        exit
      enddo

      if (j == size(old_ele) + 1) cycle

      print *, 'BMAD_PARSER: Reusing Taylor for: ', old_ele(j)%name
      call transfer_ele_taylor (old_ele(j), ele, bmad_com%taylor_order)
    enddo

    deallocate(old_ele)

  endif

!

  if (ring%param%beam_energy == 0 .and. &
                          ring%ele_(0)%value(beam_energy$) == 0) then
    print *, 'WARNING FROM BMAD_PARSER: BEAM_ENERGY IS 0!'
  elseif (ring%param%beam_energy /= 0 .and. &
                          ring%ele_(0)%value(beam_energy$) /= 0) then
    print *, &
        'WARNING FROM BMAD_PARSER: BEAM_ENERGY AND BEGINNING[ENERGY] BOTH SET!'
  endif

  doit = .true.
  if (present(make_mats6)) doit = make_mats6
  if (doit) then
    call ring_make_mat6(ring, -1)      ! make 6x6 transport matrices
  else
    call compute_element_energy (ring)
    do i = 1, ring%n_ele_max
      call attribute_bookkeeper (ring%ele_(i), ring%param)
    enddo
  endif

!-------------------------------------------------------------------------
! write out if debug is on

  if (bp_com%parser_debug) then
    
    if (index(bp_com%debug_line, 'VAR') /= 0) then
      print *
      print *, '----------------------------------------'
      print *, 'Number of Defined Variables:', &
                                      bp_com%ivar_tot - bp_com%ivar_init
      do i = bp_com%ivar_init+1, bp_com%ivar_tot
        print *
        print *, 'Var #', i
        print *, 'Name: ', bp_com%var_name(i)
        print *, 'Value:', bp_com%var_value(i)
      enddo
    endif

    if (index(bp_com%debug_line, 'SEQ') /= 0) then
      print *
      print *, '----------------------------------------'
      print *, 'Number of Lines/Lists defined:', iseq_tot
      do i = 1, iseq_tot
        print *
        print *, 'Sequence #', i
        print *, 'Name: ', sequence_(i)%name
        print *, 'Type:', sequence_(i)%type
        print *, 'Number of elements:', size(sequence_(i)%ele)
        print *, 'List:'
        do j = 1, size(sequence_(i)%ele(:))
          print '(4x, a, l3, 2i3)', sequence_(i)%ele(j)%name, &
              sequence_(i)%ele(j)%reflect, sequence_(i)%ele(j)%rep_count, &
              sequence_(i)%ele(j)%ix_arg
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

  do i = 1, size(sequence_(:))
    if (associated (sequence_(i)%arg)) deallocate(sequence_(i)%arg)
    if (associated (sequence_(i)%ele)) then
      do j = 1, size(sequence_(i)%ele)
        if (associated (sequence_(i)%ele(j)%arg)) &
                                        deallocate(sequence_(i)%ele(j)%arg)
      enddo
      deallocate(sequence_(i)%ele)
    endif
  enddo

  if (associated (in_ring%ele_))     call deallocate_ring_pointers (in_ring)
  if (associated (pring%ele))        deallocate (pring%ele)
  if (allocated (ix_ring))           deallocate (ix_ring)
  if (allocated (seq_indexx))        deallocate (seq_indexx, seq_name)
  if (allocated (in_indexx))         deallocate (in_indexx, in_name)
  if (allocated (name_))             deallocate (name_)
  if (associated (in_ring%control_)) deallocate (in_ring%control_)
  if (associated (in_ring%ic_))      deallocate (in_ring%ic_)

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
