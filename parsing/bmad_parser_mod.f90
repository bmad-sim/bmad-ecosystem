!+
! Module bmad_parser_mod
!
! This module is a collection of helper routines used by bmad_parser and bmad_parser2.
! The routines in this module are specifically taylored for bmad_parser and
! bmad_parser2 and cannot, in general, be used otherwise.
!-

module bmad_parser_mod

use ptc_interface_mod
use bookkeeper_mod
use wake_mod
use attribute_mod
use add_superimpose_mod

private parse_grid, parse_map

! A "sequence" is a line or a list.
! The information about a sequence is stored in a seq_struct.

! A seq_struct has an array of seq_ele_struct structures.
! Each seq_ele_struct represents an individual element in a sequence and, 
! since sequences can be nested, can itself be a line or a list.

type seq_ele_struct
  character(40) name             ! name of element, subline, or sublist
  character(40), pointer :: actual_arg(:) => null()
  character(40) :: tag = ''      ! tag name.
  integer type                   ! LINE$, REPLACEMENT_LINE$, LIST$, ELEMENT$
  integer ix_ele                 ! if an element: pointer to ELE array
                                 ! if a list: pointer to SEQ array
  integer ix_arg                 ! index in arg list (for replacement lines)
  integer rep_count              ! how many copies of an element
  logical reflect                ! reflection sequence
end type

type seq_struct
  character(40) name                 ! name of sequence
  type (seq_ele_struct), pointer :: ele(:) => null()
  character(40), pointer :: dummy_arg(:) => null()
  character(40), pointer :: corresponding_actual_arg(:) => null()
  integer type                       ! LINE$, REPLACEMENT_LINE$ or LIST$
  integer ix                         ! current index of element in %ELE
  integer indexx                     ! alphabetical order sorted index
  character(200) file_name     ! file where sequence is defined
  integer ix_line              ! line number in filewhere sequence is defined
  logical multipass
end type

type used_seq_struct
  character(40) :: name = ''           ! name of sequence
  character(40) :: multipass_line = '' ! name of root multipass line
  character(40) :: tag = ''             ! tag name.
  integer :: ix_multipass = 0           ! index used to sort elements
end type    

! A LIFO stack structure is used in the final evaluation of the line that is
! used to form a lattice

type seq_stack_struct
  integer ix_seq                ! index to seq(:) array
  integer ix_ele                ! index to seq%ele(:) array
  integer rep_count             ! repetition count
  integer direction             ! +1 => forwad, -1 => back reflection.
  character(40) :: tag = ''
  logical multipass
end type

! A LIFO stack structure is used to hold the list of input lattice files
! that are currently open.

type stack_file_struct
  character(200) logical_name
  character(200) :: full_name = ''
  character(200) :: dir = './'
  character(200) parse_line_saved
  integer i_line
  integer f_unit
  logical inline_call_active
end type

!-----------------------------------------------------------

integer, parameter, private :: plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
integer, parameter, private :: l_parens$ = 5, r_parens$ = 6, power$ = 7
integer, parameter, private :: unary_minus$ = 8, unary_plus$ = 9, no_delim$ = 10
integer, parameter, private :: sin$ = 11, cos$ = 12, tan$ = 13
integer, parameter, private :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17, sqrt$ = 18
integer, parameter, private :: log$ = 19, exp$ = 20, ran$ = 21, ran_gauss$ = 22
integer, parameter, private :: numeric$ = 100

integer, parameter, private :: eval_level(22) = [1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
                            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9]

type eval_stack_struct
  integer type
  real(rp) value
end type

! structure for holding the control names and pointers for superimpose and overlay elements

type parser_ele_struct
  character(40) ref_name
  character(40), pointer :: name(:) => null()
  character(40), pointer :: attrib_name(:) => null()
  character(200) lat_file    ! File where element was defined.
  integer ix_line_in_file    ! Line in file where element was defined.
  real(rp), pointer :: coef(:) => null()
  real(rp) s
  integer ix_count
  integer ele_pt, ref_pt
  integer indexx
  logical create_em_field_slave
end type

type parser_lat_struct
  type (parser_ele_struct), pointer :: ele(:) => null()
end type

!

integer, parameter :: line$ = 1, list$ = 2, element$ = 3
integer, parameter :: replacement_line$ = 4
integer, parameter :: begin$ = -1, center$ = 0, end$ = 1
integer, parameter :: def$ = 1, redef$ = 2

!------------------------------------------------
! common stuff

integer, parameter :: n_parse_line = 280

type bp_var_struct
  character(40) name      ! variable name
  real(rp) value          ! variable value
  integer :: indexx = 0   ! variable sort index
end type

type bp_common_struct
  type (ele_struct), pointer :: param_ele, beam_ele, beam_start_ele, root_branch_ele
  type (stack_file_struct), pointer :: current_file
  type (stack_file_struct), pointer :: calling_file
  type (lat_struct), pointer :: old_lat
  type (used_seq_struct), allocatable ::  used_line(:)
  type (bp_var_struct), allocatable :: var(:)   ! variable name
  type (ran_parsing_struct) ran
  integer num_lat_files               ! Number of files opened
  integer ivar_tot, ivar_init
  character(200), allocatable :: lat_file_names(:) ! List of all files used to create lat
  character(n_parse_line) parse_line
  character(n_parse_line) input_line1          ! For debug messages
  character(n_parse_line) input_line2          ! For debug messages
  character(40) parser_name
  character(200) :: dirs(2) 
  logical :: bmad_parser_calling = .false.     ! used for expand_lattice
  logical error_flag     
  logical input_line_meaningful
  logical do_superimpose
  logical write_digested      ! For bmad_parser
  logical write_digested2     ! For bmad_parser2
  logical input_from_file     ! Input is from a lattice file?
  logical inline_call_active
  logical e_tot_set, p0c_set
  logical :: always_parse = .false. ! For debugging to force parsing
end type

!

type (bp_common_struct), save :: bp_com

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_set_attribute (how, ele, lat, delim, delim_found, 
!                                            err_flag, print_err, pele, check_free)
!
! Subroutine used by bmad_parser and bmad_parser2 to get the value of
! an attribute from the input file and set the appropriate value in an element.
!
! This subroutine is not intended for general use.
!
! Input:
!   how -- Integer: Either def$ if the element is being construct from scratch or
!             redef$ if the element has already been formed and this is part of a
!             "ele_name[attrib_name] = value" construct.
!   lat -- lat_struct: Lattice. Needed if the attribute value is an expression
!             that uses values of other elements.
!   print_err    -- Logical: If False then do not print error messages.
!   check_free   -- Logical, optional: If present and True then an error will be generated
!                     if the attribute is not free to vary. Used by bmad_parser2.
!
! Output
!   ele          -- ele_struct: Element whos attribute this is.
!   delim        -- Character(1): Delimiter found where the parsing of the input line stops.
!   delim_found  -- Logical: Delimiter found? False if end of input command.
!   err_flag     -- Logical: Set True if there is a problem parsing the input.
!   pele         -- parser_ele_struct, optional: Structure to hold additional 
!                     information that cannot be stored in the ele argument.
!-

subroutine parser_set_attribute (how, ele, lat, delim, delim_found, &
                                             err_flag, print_err, pele, check_free)

use random_mod
       
implicit none

type (lat_struct), target :: lat
type (parser_ele_struct), optional :: pele
type (ele_struct), target ::  ele
type (ele_struct), target, save ::  ele0
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: bele
type (wig_term_struct), pointer :: wig_term(:)
type (real_pointer_struct), allocatable, save :: r_ptrs(:)
type (wall3d_section_struct), pointer :: section
type (wall3d_vertex_struct), pointer :: v_ptr
type (em_field_mode_struct), pointer :: em_modes(:)
type (em_field_mode_struct), pointer :: em_mode

real(rp) kx, ky, kz, tol, value, coef
real(rp), pointer :: r_ptr

integer i, j, ix_word, how, ix_word1, ix_word2, ios, ix, i_out, ix_coef
integer expn(6), ix_attrib, i_section, ix_v, ix_sec, i_mode, i_term, ib, ie, im

character(40) :: word, str_ix, attrib_word, word2
character(1) delim, delim1, delim2
character(80) str, err_str, line

logical delim_found, err_flag, logic, print_err, set_done, end_of_file, do_evaluate
logical, optional :: check_free

! Get next WORD.
! If an overlay or group element then word is just an attribute to control
! [except for a "GROUP[COMMAND] = 0.343" redef construct]

err_flag = .true.  ! assume the worst
call get_next_word (word, ix_word, ':, =()', delim, delim_found, call_check = .true.)

! taylor

if (ele%key == taylor$ .and. word(1:1) == '{') then

  word = word(2:)             ! strip off '{'
  read (word, *, iostat = ios) i_out
  if (delim /= ':' .or. ix_word == 0 .or. ios /= 0) then
    call parser_error ('BAD "OUT" IN TERM FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  call evaluate_value (str, coef, lat, delim, delim_found, err_flag)
  if (err_flag) return

  call get_next_word (line, ix_word, '},', delim, delim_found)
  read (line, *, iostat = ios) expn
  if (delim /= '}' .or. ix_word == 0 .or. ios /= 0) then
    call parser_error ('BAD "EXPONENT" IN TERM FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  call add_this_taylor_term (ele, i_out, coef, expn)
  call get_next_word (word, ix_word, '},', delim, delim_found)

  if (ix_word /= 0 .or. (delim_found .and. delim /= ',')) then
    call parser_error ('BAD TERM ENDING FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  return
endif

! overlay

if (ele%key == overlay$) then
  i = attribute_index(ele, word)       ! general attribute search

  if (i == type$ .or. i == alias$) then
    call bmad_parser_type_get (ele, word, delim, delim_found)

  else

    if (i < 1) then
      if (print_err) call parser_error ('BAD OVERLAY ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
      return
    endif

    if (how == def$) then
      ele%ix_value = i
      ele%component_name = word
    endif

    if (ele%ix_value /= i) then
      if (print_err) call parser_error ('BAD OVERLAY ATTRIBUTE SET FOR: ' // ele%name, &
            'YOU ARE TRYING TO SET: ' // word, &
            'BUT YOU SHOULD BE SETTING: ' // ele%component_name)
      return
    endif

    value = 0
    if (delim == '=') then  ! value
      call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
      if (err_flag) return
    endif
    call pointer_to_indexed_attribute (ele, i, .true., r_ptr, err_flag, .true.)
    r_ptr = value

    if (attrib_free_problem(word)) return

  endif

  err_flag = .false.
  return
endif

! group...

if (ele%key == group$) then

  if (how == def$) then
    ele0%key = overlay$
    i = attribute_index(ele0, word)       ! general attribute search
  else   ! how == redef$
    i = attribute_index(ele, word)
  endif

  if (i < 1) then
    if (print_err) call parser_error ('BAD GROUP ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
    return
  endif

  if (i == type$ .or. i == alias$) then
    call bmad_parser_type_get (ele, word, delim, delim_found)
  else
    if (how == def$) then
      ele%ix_value = i
      ele%component_name = word
    endif
    if (delim == '=') then  ! value
      call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
      if (err_flag) return
      if (how == def$) then
        ele%value(command$) = value
      else
        ele%value(i) = value
      endif
    elseif (how == redef$) then
      call parser_error ('NO VALUE GIVEN FOR ATTRIBUTE FOR: ' // ele%name)
    endif
    if (attrib_free_problem(word)) return
  endif

  err_flag = .false.
  return

endif

! beginning element or beam_start element

if (ele%key == init_ele$ .or. ele%key == def_beam_start$) then
  call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag) 
  if (err_flag) return
  call pointers_to_attribute (lat, ele%name, word, .false., r_ptrs, err_flag, .false.)
  if (err_flag .or. size(r_ptrs) /= 1) then
    if (print_err) call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  r_ptrs(1)%r = value
  if (word == 'E_TOT') then
    bp_com%e_tot_set = .true.
    bp_com%p0c_set   = .false.
  elseif (word == 'P0C') then
    bp_com%e_tot_set = .false.
    bp_com%p0c_set   = .true.
  endif

  return
endif

! Long-range wake

if (word == 'LR' .and. delim == '(') then

  call get_next_word (word, ix_word, '=', delim, delim_found)
  if (.not. delim_found) then
    call parser_error ('NO "=" SIGN FOUND', 'FOR ELEMENT: ' // ele%name)
    return
  endif
  call pointer_to_attribute (ele, 'LR(' // word, .false., r_ptr, err_flag, .false.)
  if (err_flag) then
    call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif
  call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
  r_ptr = value
  return
endif

! "WALL." redef

if (word(1:5) == 'WALL.') then
  if (delim /= '(') then
    call parser_error ('MALFORMED WALL COMPONENT REDEF IN ELEMENT: ' // ele%name)
    return
  endif

  ix_sec = evaluate_array_index (err_flag, ')', word2, '(=', delim)
  if (err_flag .or. .not. associated(ele%wall3d%section) .or. ix_sec < 0 .or. ix_sec > size(ele%wall3d%section)) then
    call parser_error('BAD ' // trim(word) // ' INDEX', 'FOR ELEMENT: ' // ele%name)
    return
  endif
  section => ele%wall3d%section(ix_sec)

  select case (word(6:))
  case ('SECTION')

    if (word2 == '.S' .and. delim == '=') then
      r_ptr => section%s

    elseif (word2 == '.DR_DS' .and. delim == '=') then
      r_ptr => section%dr_ds

    elseif (word2 == '.V' .and. delim == '(') then
      ix_v = evaluate_array_index (err_flag, ')', word, '=', delim)
      if (err_flag .or. ix_v < 0 .or. ix_v > size(section%v)) then
        call parser_error('BAD VERTEX INDEX',  'FOR ELEMENT: ' // ele%name)
        return
      endif
      v_ptr => section%v(ix_v)

      select case (word)
      case ('.X');        r_ptr => v_ptr%x
      case ('.Y');        r_ptr => v_ptr%y
      case ('.RADIUS_X'); r_ptr => v_ptr%radius_x
      case ('.RADIUS_Y'); r_ptr => v_ptr%radius_y
      case ('.TILT');     r_ptr => v_ptr%tilt
      case default
        call parser_error('BAD WALL SECTION VERTEX COMPONENT: ' // word, 'FOR ELEMENT: ' // ele%name)
        return
      end select
    else
      call parser_error('BAD WALL SECTION COMPONENT: ' // word2, 'FOR ELEMENT: ' // ele%name)
      return
    endif

  case default
    call parser_error ('BAD WALL COMPONENT REDEF: ' // word, 'IN ELEMENT: ' // ele%name)
  end select

  call evaluate_value (trim(ele%name), r_ptr, lat, delim, delim_found, err_flag)

  return
endif

! if not an overlay then see if it is an ordinary attribute.
! if not an ordinary attribute then might be a superimpose switch

if (word(:ix_word) == 'REF') word = 'REFERENCE' ! allowed abbrev

ix_attrib = attribute_index(ele, word)
attrib_word = word

if (attrib_free_problem(attrib_word)) return

if (ix_word == 0) then  ! no word
  call parser_error  ('"," NOT FOLLOWED BY ATTRIBUTE NAME FOR: ' // ele%name)
  return
endif

if (ix_attrib < 1) then
  if (ele%key == drift$ .and. (word == 'HKICK' .or. word == 'VKICK' .or. &
        word == 'BL_HKICK' .or. word == 'BL_VKICK')) then
    if (print_err) call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name, &
                      'ONE SOLUTION IS TO MAKE THIS DRIFT A "PIPE" ELEMENT.')
  else
    if (print_err) call parser_error ('BAD ATTRIBUTE NAME: ' // word, 'FOR ELEMENT: ' // ele%name)
  endif
  return
endif

! wall cross-section definition

if (attrib_word == 'WALL') then
  if (associated (ele%wall3d%section)) then
    call parser_error ('MULTIPLE WALL DEFINITIONS FOR ELEMENT: ' // ele%name)
    return
  endif

  ! Expect "= {"
  call get_next_word (word, ix_word, '{,()', delim2, delim_found, call_check = .true.)
  if (delim /= '=' .or. delim2 /= '{' .or. word /= '') then
    call parser_error ('NO "= {" FOUND AFTER "WALL"', 'FOR ELEMENT: ' // ele%name)
    return
  endif

  call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

  ! Loop over all sections or ele_anchor_pt

  i_section = 0
  section_loop: do    

    ! Possible "}" is end of wall 
    if (delim /= '}' .and. word == '') exit

    ! "ele_anchor_pt ="

    if (word == 'ELE_ANCHOR_PT') then
      if (delim /= '=') then
        call parser_error ('NO "=" FOUND AFTER WALL ELE_ANCHOR_PT FOR ELEMENT: ' // ele%name)
        return
      endif
      call get_next_word (word2, ix_word, ',}', delim, delim_found)
      call match_word(word2, anchor_pt_name, ele%wall3d%ele_anchor_pt, can_abbreviate = .false.)
      if (ele%wall3d%ele_anchor_pt < 1) then
        call parser_error ('BAD WALL3D ELE_ANCHOR_PT: ' // word2, 'FOR ELEMENT: ' // ele%name)
        return
      endif
      ! delim is parsed below so just put it back on the parse line.
      bp_com%parse_line = delim // bp_com%parse_line

    ! Must be section
    ! Expect "section = {" 

    else
      call get_next_word (word2, ix_word, '{},()', delim2, delim_found)

      if (word /= 'SECTION' .or. delim /= '=' .or. word2 /= '' .or. delim2 /= '{') then
        call parser_error ('NO "SECTION = {" SIGN FOUND IN WALL STRUCTURE', 'FOR ELEMENT: ' // ele%name)
        return
      endif

      ! Read in section

      i_section = i_section + 1
      call re_associate (ele%wall3d%section, i_section)
      section => ele%wall3d%section(i_section)

      ! Expect "S ="
      call get_next_word (word, ix_word, '{},()=', delim, delim_found)

      if (word /= 'S' .or. delim /= '=') then
        call parser_error ('EXPECTED "S =" AT START OF WALL SECTION BUT GOT: ' // trim(word) // delim, &
                             'FOR: ' // ele%name)
        return
      endif

      call evaluate_value (trim(ele%name), section%s, lat, delim, delim_found, err_flag, ',')
      if (err_flag) return
      if (ele%key == capillary$) ele%value(l$) = section%s

      ! Parse "V() = ..." constructs.

      ix_v = 0

      do
        ! Expect "V(" or "dr_ds ="
        call get_next_word (word, ix_word, '{},()=', delim, delim_found)

        if (word == 'DR_DS') then
          if (delim /= '=') then
            call parser_error ('NO "=" AFTER "DR_DS" IN WALL SECTION FOR:' // ele%name)
            return
          endif
          call evaluate_value (trim(ele%name), section%dr_ds, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return

        !

        elseif (word == 'V' .and. delim == '(') then

          ix_v = ix_v + 1
          section%n_vertex_input = ix_v
          call re_allocate (section%v, ix_v)

          call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
          read (word, *, iostat = ios) j 
          if (ios /= 0 .or. ix_v /= j) then
            call parser_error ('BAD OR OUT OF ORDER WALL SECTION VERTEX INDEX NUMBER FOR: ' // ele%name)
            return
          endif

          call get_next_word (word, ix_word, '{},()', delim2, delim_found)
          if (delim /= ')' .or. word /= '=' .or. delim2 /= '{') then        
            call parser_error ('MALFORMED ORDER WALL SECTION VERTEX FOR: ' // ele%name)
            return
          endif

          call evaluate_value (trim(ele%name), section%v(ix_v)%x, lat, delim, delim_found, err_flag, ',')
          if (err_flag) return

          call evaluate_value (trim(ele%name), section%v(ix_v)%y, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return

          if (delim == ',') then
            call evaluate_value (trim(ele%name), section%v(ix_v)%radius_x, lat, delim, delim_found, err_flag, ',}')
            if (err_flag) return
          endif

          if (delim == ',') then
            call evaluate_value (trim(ele%name), section%v(ix_v)%radius_y, lat, delim, delim_found, err_flag, ',}')
            if (err_flag) return
          endif

          if (delim == ',') then
            call evaluate_value (trim(ele%name), section%v(ix_v)%tilt, lat, delim, delim_found, err_flag, '}')
            if (err_flag) return
          endif

          call get_next_word (word, ix_word, '{},()=', delim, delim_found)
          if (word /= '' .or. (delim /= '}' .and. delim /= ',')) then
            call parser_error ('BAD SYNTAX IN WALL SECTION DEFINITION FOR ELEMENT: ' // ele%name)
            return
          endif

        else
          call parser_error ('EXPECTED "V(" BUT GOT: ' // trim(word) // delim, &
                               'IN WALL SECTION DEFINITION FOR ELEMENT: ' // ele%name)
          return
        endif

        if (delim == '}') exit

      enddo

    endif

    call get_next_word (word, ix_word, '{},()=', delim, delim_found)

    ! Possible "}" is end of wall structure
    if (delim == '}' .and. word == '') exit

    ! Must be ","
    if (word /= '' .or. delim /= ',') then
      call parser_error ('BAD SYNTAX IN WALL DEFINITION FOR ELEMENT: ' // ele%name)
    endif


    ! Expect "section" or "ele_anchor_pt"

    call get_next_word (word, ix_word, '{},()=', delim, delim_found)

    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER: ' // word, 'IN WALL STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif

    if (word /= 'SECTION' .and. word /= 'ELE_ANCHOR_PT') then
      call parser_error('DO NOT UNDERSTAND: ' // word, 'IN WALL STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif

  enddo section_loop

  ! Check for next thing on line and return

  call get_next_word (word, ix_word, '{},()', delim, delim_found)
  if (word /= '') call parser_error('EXTRA CHARACTERS AT END OF WALL SPECIFICATION IN ELEMENT: ' // ele%name) 
  return

endif

! rf field

if (attrib_word == 'FIELD') then

  ! Expect "= {"
  call get_next_word (word, ix_word, '{,()', delim2, delim_found, call_check = .true.)
  if (delim /= '=' .or. delim2 /= '{' .or. word /= '') then
    call parser_error ('NO "= {" FOUND AFTER "FIELD"', 'FOR ELEMENT: ' // ele%name)
    return
  endif

  ! Loop over all modes

  do

    if (associated(ele%em_field)) then
      i_mode = size(ele%em_field%mode) + 1
      em_modes => ele%em_field%mode
      nullify(ele%em_field)
      call init_em_field (ele%em_field, i_mode)
      ele%em_field%mode(1:i_mode-1) = em_modes 
      deallocate(em_modes)
    else
      call init_em_field (ele%em_field, 1)
      i_mode = 1
    endif

    em_mode => ele%em_field%mode(i_mode)
    if (ele%key == lcavity$ .or. ele%key == rfcavity$) em_mode%harmonic = 1 ! Default

    ! Expect "MODE = {"

    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    call get_next_word (word2, ix_word, '{},()', delim2, delim_found)
    if (word /= 'MODE' .or. delim /= '=' .or. word2 /= '' .or. delim2 /= '{') then
      call parser_error ('NO "MODE = {" SIGN FOUND IN FIELD STRUCTURE', 'FOR ELEMENT: ' // ele%name)
      return
    endif

    ! Read in mode...

    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    do

      ! Expect "<component> = "

      if (delim /= '=') then
        call parser_error ('NO "=" SIGN FOUND IN MODE DEFINITION', 'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      endif

      do_evaluate = .true.
      word2 = word

      select case (word2)

      case ('F_DAMP');         r_ptr => em_mode%f_damp
      case ('DPHI0_REF');      r_ptr => em_mode%dphi0_ref
      case ('STORED_ENERGY');  r_ptr => em_mode%stored_energy
      case ('PHI0_AZIMUTH');   r_ptr => em_mode%phi0_azimuth
      case ('FIELD_SCALE');    r_ptr => em_mode%field_scale

      case ('GRID') 
        call parse_grid(em_mode%grid, ele, lat, delim, delim_found, err_flag, print_err)
        if (err_flag) return
        do_evaluate = .false.

      case ('MAP') 
        call parse_map(em_mode%map, ele, lat, delim, delim_found, err_flag, print_err)
        if (err_flag) return
        do_evaluate = .false.

      case ('M', 'HARMONIC')
        call get_next_word (word, ix_word, ',}', delim, delim_found)
        if (.not. is_integer(word) .or. (delim /= ',' .and. delim /= '}')) then
          call parser_error ('BAD "M = <INTEGER>" CONSTRUCT', &
                               'FOUND IN MODE DEFINITION IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
          return
        endif

        if (word2 == 'M') then;  read (word, *) em_mode%m
        else;                    read (word, *) em_mode%harmonic
        endif

        do_evaluate = .false.

      case ('MASTER_SCALE')
        call get_next_word (word, ix_word, ',}', delim, delim_found)
        ix = attribute_index(ele, word)
        if (ix < 1 .or. ix > n_attrib_maxx) then
          call parser_error ('BAD NAME FOR "MASTER_SCALE = <NAME>" CONSTRUCT', &
                               'FOUND IN MODE DEFINITION IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
          return
        endif
        em_mode%master_scale = ix
        do_evaluate = .false.

      case default
        call parser_error ('UNKNOWN MODE COMPONENT: ' // word, &
                             'FOUND IN MODE DEFINITION IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      end select

      if (do_evaluate) call evaluate_value (trim(ele%name), r_ptr, lat, delim, delim_found, err_flag, ',}')

      ! Possible "}" is end of mode
      if (delim == '}') exit

      call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    enddo

    ! Check if map data has already been read in for another element.
    ! If so, save space by pointing to the existing map.

    branch_loop: do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      do ie = 1, branch%n_ele_max
        bele => branch%ele(ie)
        if (.not. associated(bele%em_field)) cycle    
        if (bele%ix_ele == ele%ix_ele .and. bele%ix_branch == ele%ix_branch) cycle
        do im = 1, size(bele%em_field%mode)
          if (.not. associated(bele%em_field%mode(im)%map)) cycle
          if (bele%em_field%mode(im)%map%file /= em_mode%map%file) cycle
          deallocate(em_mode%map)
          em_mode%map => bele%em_field%mode(im)%map
          em_mode%map%n_link = em_mode%map%n_link + 1        
          exit branch_loop
        end do
      enddo 
    enddo branch_loop

    ! Expect "," or "}"
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    if (word /= '' .or. (delim /= '}' .and. delim /= ',')) then
      call parser_error ('NO "," OR "}" FOUND AFTER MODE DEFINITION', &
                           'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif

    if (delim == '}') exit

  enddo

  call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
  err_flag = .false.
  return

endif

! wiggler term attribute

if (ix_attrib == term$ .and. ele%key == wiggler$) then

  err_flag = .true. ! assume the worst

  if (delim /= '(') then   ! ) then
    call parser_error ('"TERM" FOR A WIGGLER NOT FOLLOWED BY A "(" FOR: ' // ele%name)  ! )
    return
  endif

  call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.) ! (
  if (delim /= ')') then
    call parser_error ('CANNOT FIND CLOSING ")" for a "TERM(i)" FOR A WIGGLER"', 'FOR: ' // ele%name)
    return
  endif

  read (word, *, iostat = ios) ix
  if (ix < 1 .or. ios /= 0) then
    call parser_error ('BAD TERM NUMBER FOR A WIGGLER FOR: ' // ele%name)
    return
  endif

  write (str_ix, '(a, i3, a)') 'TERM(', ix, ')'

  if (.not. associated(ele%wig)) then
    allocate(ele%wig)
    allocate(ele%wig%term(ix))
  elseif (ix > size(ele%wig%term)) then
    allocate (wig_term(size(ele%wig%term)))
    wig_term = ele%wig%term
    deallocate (ele%wig%term)
    allocate (ele%wig%term(ix))
    ele%wig%term(1:size(wig_term)) = wig_term
    deallocate (wig_term)
  endif

! 1) chop "=", 2) chop to "{", 3) chop to "}", 4) chop to "," if it exists

  call get_next_word (word, ix_word1, ':,={}', delim1, delim_found, .true.) 
  call get_next_word (word, ix_word2, ':,={}', delim2, delim_found, .true., call_check = .true.)  

  if (delim1 /= '=' .or. delim2 /= '{' .or. ix_word1 /= 0 .or. ix_word2 /= 0) then
    call parser_error ('CONFUSED SYNTAX FOR TERM IN WIGGLER: ' // ele%name, str_ix)
    return
  endif

  err_str = trim(ele%name) // ' ' // str_ix

  call evaluate_value (err_str, ele%wig%term(ix)%coef, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return
 
  call evaluate_value (err_str, ele%wig%term(ix)%kx, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return

  call evaluate_value (err_str, ele%wig%term(ix)%ky, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return

  call evaluate_value (err_str, ele%wig%term(ix)%kz, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return

  call evaluate_value (err_str, ele%wig%term(ix)%phi_z, lat, delim, delim_found, err_flag, '}')
  if (err_flag) return

  kx = ele%wig%term(ix)%kx
  ky = ele%wig%term(ix)%ky
  kz = ele%wig%term(ix)%kz
  tol = 1e-5 * (kx**2 + ky**2 + kz**2)

  if (abs(ky**2 - kx**2 - kz**2) < tol) then
    ele%wig%term(ix)%type = hyper_y$
  elseif (abs(ky**2 + kx**2 - kz**2) < tol) then
    ele%wig%term(ix)%type = hyper_xy$
  elseif (abs(ky**2 - kx**2 + kz**2) < tol) then
    ele%wig%term(ix)%type = hyper_x$
  else
    call parser_error ('WIGGLER TERM DOES NOT HAVE CONSISTANT Kx, Ky, and Kz', &
                  'FOR WIGGLER: ' // ele%name // '  ' // str_ix)
    err_flag = .true.
    return
  endif

  call get_next_word (word, ix_word,  ':,=()', delim,  delim_found, .true.)  
  if (ix_word /= 0) then
    call parser_error ('BAD SYNTAX FOR WIGGLER: ' // ele%name, str_ix)
    err_flag = .true.
    return
  endif

  ele%sub_key = map_type$
  return

endif

! Check that next delim is a "=". 
! If not, it might be a flag attribute or an attribute that has a default value.

if (delim /= '=')  then
  err_flag = .false.

  if (ele%key == multipole$ .and. ix_attrib >= t0$) then
    if (.not. associated(ele%a_pole)) call multipole_init (ele)
    ele%b_pole(ix_attrib-t0$) = pi / (2*(ix_attrib-t0$) + 2)
    return
  endif

  if (ix_attrib == tilt$) then
    select case (ele%key)
    case (sbend$, rbend$, mirror$)
      ele%value(tilt$) = pi / 2
      return
    case (quadrupole$, sol_quad$) 
      ele%value(tilt$) = pi / 4
      return
    case (sextupole$) 
      ele%value(tilt$) = pi / 6
      return
    case (octupole$) 
      ele%value(tilt$) = pi / 8
      return
    case default
      if (attribute_name(ele, tilt$) == 'TILT') then
        call parser_error ('SORRY I''M NOT PROGRAMMED TO USE A "TILT" DEFAULT' // &
                      'FOR A: ' // key_name(ele%key), 'FOR: ' // ele%name)
        err_flag = .true.
        return
      endif
    end select
  endif

  if (ele%key == sbend$ .or. ele%key == rbend$) then
    select case (ix_attrib)
    case (fint$)
      ele%value(fint$) = 0.5
      return
    case (fintx$)
      ele%value(fintx$) = 0.5
      return
    end select
  endif

  select case (ix_attrib)

  case (superimpose$)
    ele%lord_status = super_lord$

  case (ref_beginning$)
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ref_pt = begin$

  case (create_em_field_slave$)
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%create_em_field_slave = .true.

  case (ref_center$)
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ref_pt = center$

  case (ref_end$)
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ref_pt = end$

  case (ele_beginning$)
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ele_pt = begin$

  case (ele_center$)
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ele_pt = center$

  case (ele_end$)
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ele_pt = end$

  case default
    call parser_error ('EXPECTING "=" AFTER ATTRIBUTE: ' // word,  'FOR ELEMENT: ' // ele%name)
    err_flag = .true.
  end select

  return
endif

! get the value of the attribute.
! The TYPE, ALIAS, and DESCRIP attributes are special because their "values"
! are character strings

select case (attrib_word)

case ('REFERENCE')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  call get_next_word(pele%ref_name, ix_word,  ':=,', delim, delim_found, .true.)

case ('OFFSET')
  call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
  if (err_flag) return
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%s = value

case('TYPE', 'ALIAS', 'DESCRIP', 'SR_WAKE_FILE', 'LR_WAKE_FILE', 'LATTICE', 'TO', &
     'REF_PATCH', 'CRYSTAL_TYPE')
  call bmad_parser_type_get (ele, attrib_word, delim, delim_found)

case ('SYMPLECTIFY') 
  if (how == def$ .and. (delim == ',' .or. .not. delim_found)) then
    ele%symplectify = .true.
  else
    call get_logical (attrib_word, ele%symplectify, err_flag)
  endif
  
case ('IS_ON')
  call get_logical (attrib_word, ele%is_on, err_flag)

case ('APERTURE_LIMIT_ON') 
  call get_logical (attrib_word, lat%param%aperture_limit_on, err_flag)

case ('ABSOLUTE_TIME_TRACKING')
  call get_logical (attrib_word, lat%absolute_time_tracking, err_flag)

case ('USE_PTC_LAYOUT')
  call get_logical (attrib_word, lat%use_ptc_layout, err_flag)

case ('RF_AUTO_SCALE_PHASE')
  call get_logical (attrib_word, lat%rf_auto_scale_phase, err_flag)

case ('RF_AUTO_SCALE_AMP')
  call get_logical (attrib_word, lat%rf_auto_scale_amp, err_flag)

case ('CSR_CALC_ON')
  call get_logical (attrib_word, ele%csr_calc_on, err_flag)

case ('MAP_WITH_OFFSETS')
  call get_logical (attrib_word, ele%map_with_offsets, err_flag)

case ('OFFSET_MOVES_APERTURE')
  call get_logical (attrib_word, ele%offset_moves_aperture, err_flag)

case ('FIELD_MASTER')
  call get_logical (attrib_word, ele%field_master, err_flag)

case ('SCALE_MULTIPOLES')
  call get_logical (attrib_word, ele%scale_multipoles, err_flag)

case ('DIFFRACTION_TYPE')
  call get_switch (attrib_word, diffraction_type_name(1:), ele%sub_key, err_flag)

case ('FIELD_CALC')
  call get_switch (attrib_word, field_calc_name(1:), ele%field_calc, err_flag)

case ('ROOT_BRANCH_NAME')
  call get_next_word(bp_com%root_branch_ele%name, ix_word,  ':=,', delim, delim_found, .true.)

case default   ! normal attribute

  if (attribute_type(attrib_word) == is_logical$) then
    call get_logical_real (attrib_word, ele%value(ix_attrib), err_flag)


  else
    call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
    if (err_flag) return

    if (ix_attrib >= a0$ .and. ix_attrib <= b20$) then  ! multipole attribute
        if (.not. associated(ele%a_pole)) call multipole_init (ele)
        if (ix_attrib >= b0$) then
          ele%b_pole(ix_attrib-b0$) = value
        else
          ele%a_pole(ix_attrib-a0$) = value
        endif
    elseif (ix_attrib == mat6_calc_method$) then
      ele%mat6_calc_method = nint(value)
    elseif (ix_attrib == tracking_method$) then
      ele%tracking_method = nint(value)
    elseif (ix_attrib == spin_tracking_method$) then
      ele%spin_tracking_method = nint(value)
    elseif (ix_attrib == ref_orbit$) then
      ele%ref_orbit = nint(value)
    elseif (ix_attrib == aperture_at$) then
      ele%aperture_at = nint(value)
    elseif (ix_attrib == aperture_type$) then
      ele%aperture_type = nint(value)
    elseif (ix_attrib == ran_seed$) then
      call ran_seed_put (nint(value))  ! init random number generator
      if (nint(value) == 0) then  ! Using system clock -> Not determinisitc.
        bp_com%ran%deterministic = 0
      else
        bp_com%ran%deterministic = 2
      endif
    elseif (ix_attrib == aperture$) then
      ele%value(x1_limit$) = value
      ele%value(x2_limit$) = value
      ele%value(y1_limit$) = value
      ele%value(y2_limit$) = value
    elseif (ix_attrib == x_limit$) then
      ele%value(x1_limit$) = value
      ele%value(x2_limit$) = value
    elseif (ix_attrib == y_limit$) then
      ele%value(y1_limit$) = value
      ele%value(y2_limit$) = value
    else
      ele%value(ix_attrib) = value

      ix = len_trim(attrib_word)
      if (ix > 9 .and. index(attrib_word, '_GRADIENT') == ix-8) ele%field_master = .true.
      if (ix > 6 .and. index(attrib_word, '_FIELD') == ix-5) ele%field_master = .true.
      if (ix > 10 .and. index(attrib_word, '_FIELD_ERR') == ix-9) ele%field_master = .true.
      if (attrib_word(1:3) == 'BL_') ele%field_master = .true.

      select case (ix_attrib)

      case (num_steps$)
        ele%value(ds_step$) = abs(ele%value(l$) * nint(ele%value(num_steps$)))

      case (e_tot$)
        if (ele%key == def_beam$ .or. ele%key == def_parameter$) then
          lat%ele(0)%value(e_tot$) = value
          if (ele%key == def_beam$) lat%ele(0)%value(e_tot$) = 1d9 * value
          bp_com%e_tot_set = .true.
          bp_com%p0c_set   = .false.
        endif

      case (p0c$)
        if (ele%key == def_beam$ .or. ele%key == def_parameter$) then
          lat%ele(0)%value(p0c$) = value
          if (ele%key == def_beam$) lat%ele(0)%value(p0c$) = 1d9 * value
          bp_com%e_tot_set = .false.
          bp_com%p0c_set   = .true.
        endif

      case (lr_freq_spread$)
        call randomize_lr_wake_frequencies (ele, set_done)
        if (set_done) call bp_set_ran_status

      end select

    endif

  endif

end select

err_flag = .false.

!--------------------------------------------------------
contains

function attrib_free_problem (attrib_name) result (is_problem)

character(*) attrib_name
logical is_problem, is_free

!

is_problem = .false.

if (logic_option(.false., check_free)) then
  is_free = attribute_free (ele, attrib_name, lat, print_err)
  if (.not. is_free) then
    if (print_err) call parser_error ('ATTRIBUTE NOT FREE TO BE SET: ' // attrib_name, 'FOR: ' // ele%name)
    err_flag = .true.
    is_problem = .true.
  endif
endif

end function attrib_free_problem

!--------------------------------------------------------
! contains

subroutine get_logical (attrib_name, this_logic, err)

character(*) attrib_name
logical this_logic, err

!

call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.)
this_logic = evaluate_logical (word, ios)
if (ios /= 0 .or. ix_word == 0) then
  call parser_error ('BAD "' // trim(attrib_name) // '" SWITCH FOR: ' // ele%name, 'I DO NOT UNDERSTAND: ' // word)
  err = .true.
else
  err = .false.
endif

end subroutine get_logical

!--------------------------------------------------------
! contains

subroutine get_switch (name, name_list, this_switch, err)

character(*) name, name_list(:)
integer this_switch
logical err

!

call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.)
call match_word (word, name_list, this_switch, can_abbreviate = .false.)
if (this_switch < 1) then
  call parser_error ('BAD "' // trim(name) // '" SWITCH FOR: ' // ele%name, 'I DO NOT UNDERSTAND: ' // word)
  err = .true.
else
  err = .false.
endif

end subroutine get_switch

!--------------------------------------------------------
! contains

subroutine get_logical_real (name, logic_real, err)

character(*) name
real(rp) logic_real
logical this_logical, err

!

call get_logical (name, this_logical, err)
if (err) return

if (this_logical) then
  logic_real = 1
else
  logic_real = 0
endif

err = .false.

end subroutine get_logical_real

end subroutine parser_set_attribute 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine get_called_file (delim, call_file, xsif_called, err)

implicit none

character(1) delim
character(*) call_file

integer ix_word, ix, n
logical delim_found, finished, xsif_called, err

!

err = .true.

if (delim /= ',')  call parser_error ('"CALL" NOT FOLLOWED BY COMMA', stop_here = .true.)
call get_next_word(call_file, ix_word, ':=,', delim, delim_found, .true.)

if (ix_word == 0) then
  call parser_error ('NOTHING AFTER "CALL"', stop_here = .true.)
  return
elseif (index('FILENAME', call_file(:ix_word)) /= 1) then
  call parser_error ('INVALID "CALL" COMMAND', stop_here = .true.)
  return
elseif (delim /= '=') then
  call parser_error ('NO "=" AFTER "FILENAME"', stop_here = .true.)
  return
endif

call get_next_word(call_file, ix_word, ',', delim, delim_found, .false.)
if (ix_word == 0) then
  call parser_error ('NO FILE NAME SPECIFIED', stop_here = .true.)
  return
endif

if (call_file(1:1) == '"') then
  call_file = call_file(2:)
  ix = index(call_file, '"')
  if (ix == 0 .or. ix /= len_trim(call_file)) then
    call parser_error ('MISSING DOUBLE QUOTE MARK (") FOR CALL STATEMENT', stop_here = .true.)
    return
  endif
  call_file(ix:ix) = ' '
endif

if (call_file(1:1) == "'") then
  call_file = call_file(2:)
  ix = index(call_file, "'")
  if (ix == 0 .or. ix /= len_trim(call_file)) then
    call parser_error ("MISSING SINGLE QUOTE MARK (') FOR CALL STATEMENT", stop_here = .true.)
    return
  endif
  call_file(ix:ix) = ' '
endif

if (call_file(1:6) == 'xsif::') then
  call_file = call_file(7:)
  n = size(bp_com%lat_file_names)
  if (n < bp_com%num_lat_files + 1) &
              call re_allocate (bp_com%lat_file_names, n + 100)
  bp_com%num_lat_files = bp_com%num_lat_files + 1 
  inquire (file = call_file, name = bp_com%lat_file_names(bp_com%num_lat_files))
  xsif_called = .true.
  err = .false.
  return
endif

xsif_called = .false.
call parser_file_stack ('push', call_file, finished, err) ! err gets set here

end subroutine get_called_file

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine add_this_taylor_term (ele, i_out, coef, expn) 
!
! Subroutine used by bmad_parser and bmad_parser2 to parse the input file.
! This subroutine is not intended for general use.
!-

subroutine add_this_taylor_term (ele, i_out, coef, expn)

implicit none

type (ele_struct) ele
real(rp) coef
integer i_out, expn(6)

!

if (.not. associated(ele%taylor(1)%term)) then
  allocate (ele%taylor(1)%term(0), ele%taylor(2)%term(0))
  allocate (ele%taylor(3)%term(0), ele%taylor(4)%term(0))
  allocate (ele%taylor(5)%term(0), ele%taylor(6)%term(0))
endif

if (i_out < 1 .or. i_out > 6) then
  call parser_error ('"OUT" VALUE IN TAYLOR TERM NOT IN RANGE (1 - 6)', &
                'FOR TAYLOR ELEMENT: ' // ele%name)
  return
endif

call add_taylor_term (ele%taylor(i_out), coef, expn, .true.)

end subroutine add_this_taylor_term

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_next_word (word, ix_word, delim_list, delim, delim_found, upper_case_word)
!
! Subroutine to get the next word from the input stream.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   word            -- Character(*): Word returned
!   delim_list      -- Character(*): List of valid delimiters
!   upper_case_word -- Logical, optional: if True then convert word to 
!                       upper case. Default is True.
!   call_check      -- Logical, optional: If present and True then check for 'call::<filename>' construct.
!
! Output
!   ix_word     -- Integer: length of word argument
!   delim       -- Character1: Actual delimiter found
!   delim_found -- Logical: Set true if a delimiter found. A delimiter
!                    may not be found if the end of the line is reached first.
!-


subroutine get_next_word (word, ix_word, delim_list, delim, delim_found, upper_case_word, call_check)

implicit none

integer ix_a, ix_word

character(*) word, delim_list, delim

integer n, ix

logical delim_found, end_of_file
logical, optional :: upper_case_word, call_check

character(100) line
character(6) str

! Possible inline call...

if (logic_option(.false., call_check)) then
  call string_trim(bp_com%parse_line, bp_com%parse_line, ix)
  call str_upcase (str, bp_com%parse_line(1:6))
  if (str == 'CALL::') then
    bp_com%parse_line = bp_com%parse_line(7:)
    call word_read (bp_com%parse_line, ',} ',  line, ix_word, delim, delim_found, bp_com%parse_line)
    bp_com%parse_line = delim // bp_com%parse_line  ! put delim back on parse line.
    call parser_file_stack ('push_inline', line)
  endif
endif

! check for continuation character and, if found, then load more characters
! into the parse line from the lattice file. 
! If the input is not from a file then skip this.

if (bp_com%input_from_file) then 
  do
    n = len_trim(bp_com%parse_line)
    if (n == 0 .or. n > 60) exit

    select case (bp_com%parse_line(n:n))
    case (',', '+', '-', '*', '/', '(', '{', '[', '=')
      call load_parse_line('continue', n+2, end_of_file)
    case ('&')
      call load_parse_line('continue', n, end_of_file)

    case default
      ! If in an inline called file then make sure the rest of the file is blank and
      ! return to the calling file
      if (bp_com%inline_call_active) then
        call load_parse_line('continue', n+2, end_of_file)
        if (bp_com%parse_line(n+1:) /= '') then
          call string_trim (bp_com%parse_line(n+1:), line, ix)
          call str_upcase (line(1:10), line(1:10))
          if (line /= 'END_FILE') THEN
            call parser_error ('EXTRA STUFF IN FILE')
          endif
        endif
        bp_com%parse_line(n+1:) = ''
        call parser_file_stack ('pop')
      else
        exit
      endif

    end select

  enddo
endif

! Get the first word in bp_com%parse_line

call word_read (bp_com%parse_line, delim_list,  word, ix_word, delim, delim_found, bp_com%parse_line)

if (len(word) < ix_word) then
  call parser_error ('BAD WORD: ' // bp_com%parse_line)
  ix_word = len(word)
endif

if (present(upper_case_word)) then
  if (upper_case_word) call str_upcase (word, word)
else
  call str_upcase (word, word)
endif

! Note: "var := num" is old-style variable definition syntax.
! If delim is ":" and next char is "=" then use "=" as the delim

if (delim == ':' .and. index(delim_list, '=') /= 0 .and. bp_com%parse_line(1:1) == '=') then
  delim = '='
  bp_com%parse_line = bp_com%parse_line(2:)
endif

end subroutine get_next_word

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_file_stack (how, file_name_in, finished, err)
!
! Subroutine to keep track of the files that are opened for reading.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_file_stack (how, file_name_in, finished, err)

implicit none

integer, parameter :: f_maxx = 20
type (stack_file_struct), save, target :: file(0:f_maxx)

integer, save :: i_level
integer i, ix, ios, n

character(*) how
character(*), optional :: file_name_in
character(200) file_name, basename, file_name2

logical, optional :: finished, err
logical found_it, is_relative, valid, err_flag

! "Init" means init

if (present(err)) err = .true.

if (how == 'init') then
  i_level = 0
  call fullfilename('./', file_name)
  bp_com%dirs(2) = file_name
  file(:)%dir = file_name
  if (present(err)) err = .false.
  if (.not. allocated(bp_com%lat_file_names)) allocate(bp_com%lat_file_names(100))
  bp_com%inline_call_active = .false.
  return
endif

! "push" means open a file and put its name on the stack.

if (present(finished)) finished = .false.

select case (how)
case ('push', 'push_inline')

  i_level = i_level + 1    ! number of files currently open
  if (i_level > f_maxx) then
    call parser_error ('CALL NESTING GREATER THAN 20 LEVELS')
    if (bmad_status%exit_on_error) call err_exit
  endif

  if (how == 'push_inline') then
    file(i_level)%parse_line_saved = bp_com%parse_line
    file(i_level)%inline_call_active = .true.    
    bp_com%parse_line = '&'
    bp_com%inline_call_active = .true.
  endif

  bp_com%current_file => file(i_level)
  bp_com%calling_file => file(i_level-1)

  if (i_level == 1) then   ! if we are just starting out then init some vars.
    bp_com%num_lat_files = 0           ! total number of files opened
    bp_com%error_flag = .false.  ! set to true on an error
    bp_com%current_file%full_name = ' '
    bp_com%input_line_meaningful = .false.
    call init_bmad_parser_common
  endif

  call fullfilename (file_name_in, file_name2, valid)
  if (.not. valid) then
    call parser_error ('MALFORMED FILE NAME: ' // file_name_in, stop_here = .true.)
    if (bmad_status%exit_on_error) call err_exit
    do i = 1, i_level-1
      close (file(i_level)%f_unit)
    enddo
    return
  endif

  ix = splitfilename (file_name2, file(i_level)%dir, basename, is_relative)
  if (is_relative) then
    call append_subdirectory (trim(file(i_level-1)%dir), file(i_level)%dir, file(i_level)%dir, err_flag)
    if (err_flag) call parser_error ('BAD DIRECTORY SYNTAX FOR: ' // file_name, stop_here = .true.)
  endif
  bp_com%dirs(1) = file(i_level-1)%dir
  call find_file (file_name2, found_it, file_name, bp_com%dirs)
  file(i_level)%logical_name = file_name2
  file(i_level)%full_name = file_name
  file(i_level)%f_unit = lunget()

  open (file(i_level)%f_unit, file = file_name, status = 'OLD', action = 'READ', iostat = ios)
  if (ios /= 0 .or. .not. found_it) then
    bp_com%current_file => file(i_level-1)  ! For warning
    if (file_name2 == file_name)  then
      call parser_error ('UNABLE TO OPEN FILE: ' // file_name, stop_here = .true.)
    else
      call parser_error ('UNABLE TO OPEN FILE: ' // file_name, &
                    'THIS FROM THE LOGICAL FILE NAME: ' // file_name2, &
                    stop_here = .true.)
    endif
    do i = 1, i_level-1
      close (file(i_level)%f_unit)
    enddo
    return
  endif

  bp_com%current_file%i_line = 0

  n = size(bp_com%lat_file_names)
  if (n < bp_com%num_lat_files + 1) call re_allocate (bp_com%lat_file_names, n + 100)
  bp_com%num_lat_files = bp_com%num_lat_files + 1 
  inquire (file = file_name, name = bp_com%lat_file_names(bp_com%num_lat_files))

! "pop" means close the current file and pop its name off the stack

case ('pop')
  close (unit = bp_com%current_file%f_unit)
  i_level = i_level - 1
  if (i_level < 0) then
    call parser_error ('BAD "RETURN"')
    return
  elseif (i_level > 0) then
    bp_com%current_file => file(i_level)
    bp_com%calling_file => file(i_level-1)
  else    ! i_level == 0
    if (present(finished)) finished = .true.
  endif

  if (bp_com%inline_call_active) then
    bp_com%parse_line = trim(bp_com%parse_line) // ' ' // file(i_level+1)%parse_line_saved
    bp_com%inline_call_active = file(i_level+1)%inline_call_active
  endif

  bp_com%inline_call_active = .false.

! Programming error

case default
  call parser_error ('INTERNAL ERROR IN PARSER_FILE_STACK SUBROUTINE!')
  if (bmad_status%exit_on_error) call err_exit
end select

if (present(err)) err = .false.

end subroutine parser_file_stack

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine load_parse_line (load_type, ix_start, end_of_file) 
!
! Subroutine to load characters from the input file.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   load_type -- Character(*): 'continue' or 'normal'
!   ix_start  -- Integer: index in bp_com%parse_line string where to append stuff.
!
! Output:
!   end_of_file       -- Logical: 
!   bp_com%parse_line -- String to append to.
!-

subroutine load_parse_line (load_type, ix_start, end_of_file)

implicit none

integer ix_start, ix

character(*) load_type
character(n_parse_line+20), save :: line, saved_line
character(1), parameter :: tab = achar(9)

logical :: have_saved_line = .false., end_of_file

!

end_of_file = .false.

do

  ! Read a line or use saved_line if it exists

  if (have_saved_line) then
    line = saved_line
    have_saved_line = .false.
  else
    read (bp_com%current_file%f_unit, '(a)', end = 9000) line
    bp_com%current_file%i_line = bp_com%current_file%i_line + 1
    if (line(n_parse_line-ix_start-20:) /= ' ') &
      call parser_error ('INPUT LINE HAS TOO MANY CHARACTERS:', line)
  endif

  ! %input_line1 and %input_line2 are for error messages if needed.
  ! Only the input string being parsed is saved in these lines.
  ! 'normal' load_type means we are loading a new input string so start from scratch.
  ! 'continue' load_type means keep the existing input string.

  if (load_type == 'continue') then
    bp_com%input_line1 = bp_com%input_line2
    bp_com%input_line2 = line
  elseif (load_type == 'normal') then
    bp_com%input_line1 = ' '
    bp_com%input_line2 = line
  else
    call parser_error ('INTERNAL ERROR #4: CALL HELP')
    if (bmad_status%exit_on_error) call err_exit
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
    have_saved_line = .true.
    saved_line = line(ix+1:)
    line = ' '
  elseif (ix > 1) then
    have_saved_line = .true.
    saved_line = line(ix+1:)
    line = line(:ix-1)
  else
    have_saved_line = .false.
  endif

  ! if the command line is blank then go back for more input

  call string_trim (line, line, ix)
  if (ix /= 0 .or. have_saved_line) exit
enddo

! now simply append the line to %parse_line starting at ix_start

call str_substitute (line)
bp_com%parse_line(ix_start:) = line

return

!

9000  continue
end_of_file = .true.

if (bp_com%parse_line /= '' .and. .not. bp_com%inline_call_active) then
  call parser_error ('FILE ENDED BEFORE PARSING FINISHED', stop_here = .true.)
endif

end subroutine load_parse_line

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function evaluate_array_index (err_flag, delim_list1, word2, delim_list2, delim2) result (this_index)
!
! Function of evaluate the index of an array. Typically the text peing parsed looks like:
!      "5) = ..."  or
!      "6).COMP = ..."
!
! Input:
!   delim_list1 -- Character(1): Delimitor after the integer. Normally ')'.
!   delim_list2 -- Character(*): Delimitor list to mark the end of word2. Normally '='.
!
! Output:
!   err_flag    -- Logical: Set True if there is an error. False otherwise.
!   word2       -- Character(*): Word found after delim1. Normally this should be blank.
!   delim2      -- Character(1): Actual delimitor found after word2.
!   this_index  -- Integer: Integer value
!-

function evaluate_array_index (err_flag, delim_list1, word2, delim_list2, delim2) result (this_index)

implicit none

character(*) delim_list1, word2, delim_list2, delim2
integer this_index, ix_word
character(1) delim
character(20) word
logical err_flag, delim_found


! Init

err_flag = .true.
this_index = -1

! Get integer

call get_next_word (word, ix_word, delim_list1, delim, delim_found)
if (.not. delim_found) return

if (.not. is_integer(word)) return
read (word, *) this_index

! Get word after integer

call get_next_word (word2, ix_word, delim_list2, delim2, delim_found)
if (delim_found) err_flag = .false.

end function evaluate_array_index

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function evaluate_logical (word, iostat) result (this_logic)
!
! Function of convert a string into a logical value.
! Accepted possibilities are:
!   .TRUE.  .FALSE. 
!    TRUE    FALSE
!    T       F
!
! Input:
!   word   -- Character(*): Input string.
!
! Output:
!   this_logic -- Logical: Result.
!   iostat -- Integer: Status: Returns 0 if conversion successful. 
!-

function evaluate_logical (word, iostat) result (this_logic)

implicit none

character(*), intent(in) :: word
character(len(word)+8) :: wd
logical this_logic
integer, intent(out) :: iostat
integer i

!

iostat = -1

call str_upcase(wd, word)
do i = 1, len(word)
  if (wd(i:i) /= ' ') then
    if (wd(i:i+6) == '.TRUE. ' .or. wd(i:i+4) == 'TRUE ' .or. &
                                          wd(i:i+1) == 'T ') then
      this_logic = .true.
      iostat = 0
    elseif (wd(i:i+7) == '.FALSE. ' .or. wd(i:i+5) == 'FALSE ' .or. &
                                          wd(i:i+1) == 'F ') then
      this_logic = .false.
      iostat = 0
    endif
    return
  endif
enddo

end function evaluate_logical

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine evaluate_value (err_str, value, lat, delim, delim_found, err_flag, end_delims)
!
! This routine creates an "evaluation stack" structure which can be used 
! to evaluate an arithmethic expression.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine evaluate_value (err_str, value, lat, delim, delim_found, err_flag, end_delims)

use random_mod

implicit none

type (lat_struct)  lat
type (eval_stack_struct) stk(200)

integer i_lev, i_op, i

integer op(200), ix_word, i_delim, i2, ix_word2

real(rp) value

character(*) err_str
character(*), optional :: end_delims
character(1) delim
character(80) word, word2

logical delim_found, split, ran_function_pending, first_get_next_word_call
logical err_flag

! The general idea is to rewrite the expression on a stack in reverse polish.
! Reverse polish means that the operand goes last so that 2 * 3 is writen 
! on the stack as: [2, 3, *]

! The stack is called: stk
! Since operations move towards the end of the stack we need a separate
! stack called op which keeps track of what operations have not yet
! been put on stk.

! init

err_flag = .true.
i_lev = 0
i_op = 0
ran_function_pending = .false.
first_get_next_word_call = .true.

! parsing loop to build up the stack.

parsing_loop: do

! get a word

  if (first_get_next_word_call) then
    call get_next_word (word, ix_word, '+-*/()^,:} ', delim, delim_found, call_check = .true.)
    first_get_next_word_call = .false.
  else
    call get_next_word (word, ix_word, '+-*/()^,:} ', delim, delim_found)
  endif

  if (delim == '*' .and. word(1:1) == '*') then
    call parser_error ('EXPONENTIATION SYMBOL IS "^" AS OPPOSED TO "**"!',  &
                  'for: ' // err_str)
    return
  endif

  if (ran_function_pending .and. (ix_word /= 0 .or. delim /= ')')) then
    call parser_error ('RAN AND RAN_GAUSS DO NOT TAKE AN ARGUMENT', 'FOR: ' // err_str)
    return
  endif

!--------------------------
! Preliminary: If we have split up something that should have not been split
! then put it back together again...

! just make sure we are not chopping a number in two, e.g. "3.5e-7" should not
! get split at the "-" even though "-" is a delimiter

  split = .true.         ! assume initially that we have a split number
  if (ix_word == 0) then
    split = .false.
  elseif (word(ix_word:ix_word) /= 'E' .and. word(ix_word:ix_word) /= 'D') then
    split = .false.
  endif
  if (delim(1:1) /= '-' .and. delim(1:1) /= '+') split = .false.
  do i = 1, ix_word-1
    if (index('.0123456789', word(i:i)) == 0) split = .false.
  enddo

! If still SPLIT = .TRUE. then we need to unsplit

  if (split) then
    word = word(:ix_word) // delim
    call get_next_word (word2, ix_word2, '+-*/()^,:}', delim, delim_found)
    word = word(:ix_word+1) // word2
    ix_word = ix_word + ix_word2
  endif

! Something like "lcav[lr(2).freq]" will get split on the "("

  if (delim == '(' .and. index(word, '[LR') /= 0) then
    call get_next_word (word2, ix_word2, '+-*/(^,:}', delim, delim_found)
    word = word(:ix_word) // '(' // word2
    ix_word = ix_word + ix_word2 + 1
  endif

!---------------------------
! Now see what we got...

! For a "(" delim we must have a function

  if (delim == '(') then

    ran_function_pending = .false.
    if (ix_word /= 0) then
      select case (word)
      case ('SIN') 
        call pushit (op, i_op, sin$)
      case ('COS') 
        call pushit (op, i_op, cos$)
      case ('TAN') 
        call pushit (op, i_op, tan$)
      case ('ASIN') 
        call pushit (op, i_op, asin$)
      case ('ACOS') 
        call pushit (op, i_op, acos$)
      case ('ATAN') 
        call pushit (op, i_op, atan$)
      case ('ABS') 
        call pushit (op, i_op, abs$)
      case ('SQRT') 
        call pushit (op, i_op, sqrt$)
      case ('LOG') 
        call pushit (op, i_op, log$)
      case ('EXP') 
        call pushit (op, i_op, exp$)
      case ('RAN') 
        call pushit (op, i_op, ran$)
        ran_function_pending = .true.
        call bp_set_ran_status
      case ('RAN_GAUSS') 
        call pushit (op, i_op, ran_gauss$)
        ran_function_pending = .true.
        call bp_set_ran_status
      case default
        call parser_error ('UNEXPECTED CHARACTERS ON RHS BEFORE "(": ' // word,  &
                                                'FOR: ' // err_str)
        return
      end select
    endif

    call pushit (op, i_op, l_parens$)
    cycle parsing_loop

! for a unary "-"

  elseif (delim == '-' .and. ix_word == 0) then
    call pushit (op, i_op, unary_minus$)
    cycle parsing_loop

! for a unary "+"

  elseif (delim == '+' .and. ix_word == 0) then
    call pushit (op, i_op, unary_plus$)
    cycle parsing_loop

! for a ")" delim

  elseif (delim == ')') then
    if (ix_word == 0) then
      if (.not. ran_function_pending) then
        call parser_error  ('CONSTANT OR VARIABLE MISSING BEFORE ")"', 'FOR: ' // err_str)
        return
      endif
      ran_function_pending = .false.
    else
      call word_to_value (word, lat, value)
      call pushit (stk%type, i_lev, numeric$)
      stk(i_lev)%value = value
    endif

    do
      do i = i_op, 1, -1     ! release pending ops
        if (op(i) == l_parens$) exit          ! break do loop
        call pushit (stk%type, i_lev, op(i))
      enddo

      if (i == 0) then
        call parser_error ('UNMATCHED ")" ON RHS', 'FOR: ' // err_str)
        return
      endif

      i_op = i - 1

      call get_next_word (word, ix_word, '+-*/()^,:}', delim, delim_found)
      if (ix_word /= 0) then
        call parser_error ('UNEXPECTED CHARACTERS ON RHS AFTER ")"',  &
                                                  'FOR: ' // err_str)
        return
      endif

      if (delim /= ')') exit  ! if no more ')' then no need to release more
    enddo


    if (delim == '(') then
      call parser_error ('")(" CONSTRUCT DOES NOT MAKE SENSE FOR: ' // err_str)
      return
    endif

! For binary "+-/*^" delims

  else
    if (ix_word == 0) then
      call parser_error ('CONSTANT OR VARIABLE MISSING IN EVALUATING: ' // err_str)
      return
    endif
    call word_to_value (word, lat, value)
    call pushit (stk%type, i_lev, numeric$)
    stk(i_lev)%value = value
  endif

! If we are here then we have an operation that is waiting to be identified

  if (.not. delim_found) delim = ':'

  select case (delim)
  case ('+')
    i_delim = plus$
  case ('-')
    i_delim = minus$
  case ('*')
    i_delim = times$
  case ('/')
    i_delim = divide$
  case (')')
    i_delim = r_parens$
  case ('^')
    i_delim = power$
  case (',', '}', ':')
    i_delim = no_delim$
  case default
    call parser_error ('MALFORMED EXPRESSION')
    bp_com%parse_line = ' '
    return
  end select

! now see if there are operations on the OP stack that need to be transferred
! to the STK stack

  do i = i_op, 1, -1
    if (eval_level(op(i)) >= eval_level(i_delim)) then
      if (op(i) == l_parens$) then
        call parser_error ('UNMATCHED "(" IN EVALUATING: ' // err_str)
        return
      endif
      call pushit (stk%type, i_lev, op(i))
    else
      exit
    endif
  enddo

! put the pending operation on the OP stack

  i_op = i
  if (i_delim == no_delim$) then
    exit parsing_loop
  else
    call pushit (op, i_op, i_delim)
  endif

enddo parsing_loop

!------------------------------------------------------------------
! now go through the stack and perform the operations

if (i_op /= 0) then
  call parser_error ('UNMATCHED "(" IN EVALUATING: ' // err_str)
  return
endif

if (i_lev == 0) call parser_error ('NO VALUE FOUND FOR: ' // err_str)

i2 = 0
do i = 1, i_lev
  if (stk(i)%type == numeric$) then
    i2 = i2 + 1
    stk(i2)%value = stk(i)%value
  elseif (stk(i)%type == unary_minus$) then
    stk(i2)%value = -stk(i2)%value
  elseif (stk(i)%type == unary_plus$) then
    stk(i2)%value = stk(i2)%value
  elseif (stk(i)%type == plus$) then
    stk(i2-1)%value = stk(i2-1)%value + stk(i2)%value
    i2 = i2 - 1
  elseif (stk(i)%type == minus$) then
    stk(i2-1)%value = stk(i2-1)%value - stk(i2)%value
    i2 = i2 - 1
  elseif (stk(i)%type == times$) then
    stk(i2-1)%value = stk(i2-1)%value * stk(i2)%value
    i2 = i2 - 1
  elseif (stk(i)%type == divide$) then
    if (stk(i2)%value == 0) then
      call parser_error ('DIVIDE BY 0 ON RHS', 'FOR: ' // err_str)
      return
    endif
    stk(i2-1)%value= stk(i2-1)%value / stk(i2)%value
    i2 = i2 - 1
  elseif (stk(i)%type == power$) then
    stk(i2-1)%value = stk(i2-1)%value**stk(i2)%value
    i2 = i2 - 1
  elseif (stk(i)%type == sin$) then
    stk(i2)%value = sin(stk(i2)%value)
  elseif (stk(i)%type == cos$) then
    stk(i2)%value = cos(stk(i2)%value)
  elseif (stk(i)%type == tan$) then
    stk(i2)%value = tan(stk(i2)%value)
  elseif (stk(i)%type == asin$) then
    stk(i2)%value = asin(stk(i2)%value)
  elseif (stk(i)%type == acos$) then
    stk(i2)%value = acos(stk(i2)%value)
  elseif (stk(i)%type == atan$) then
    stk(i2)%value = atan(stk(i2)%value)
  elseif (stk(i)%type == abs$) then
    stk(i2)%value = abs(stk(i2)%value)
  elseif (stk(i)%type == sqrt$) then
    stk(i2)%value = sqrt(stk(i2)%value)
  elseif (stk(i)%type == log$) then
    stk(i2)%value = log(stk(i2)%value)
  elseif (stk(i)%type == exp$) then
    stk(i2)%value = exp(stk(i2)%value)
  elseif (stk(i)%type == ran$) then
    i2 = i2 + 1
    call ran_uniform(stk(i2)%value)
  elseif (stk(i)%type == ran_gauss$) then
    i2 = i2 + 1
    call ran_gauss(stk(i2)%value)
  else
    call parser_error ('INTERNAL ERROR #02: GET HELP')
    if (bmad_status%exit_on_error) call err_exit
  endif
enddo

if (i2 /= 1) then
  call parser_error ('INTERNAL ERROR #03: GET HELP')
  if (bmad_status%exit_on_error) call err_exit
endif

value = stk(1)%value
err_flag = .false.

! Check that final delim matches.

if (present(end_delims)) then
  if (.not. delim_found .or. index(end_delims, delim) == 0) then
    call parser_error ('BAD DELIMITOR AFTER VALUE FOR: ' // err_str)
    err_flag = .true.
  endif
endif


end subroutine evaluate_value

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine pushit (stack, i_lev, value)

implicit none

integer stack(:), i_lev, value

!

i_lev = i_lev + 1

if (i_lev > size(stack)) then
  call parser_error ('STACK OVERFLOW.', 'EXPERT HELP IS NEEDED!')
  if (bmad_status%exit_on_error) call err_exit
endif

stack(i_lev) = value

end subroutine pushit
                     
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine word_to_value (word, lat, value)

implicit none

type (lat_struct), target ::  lat
type (ele_struct), pointer :: ele
type (real_pointer_struct), allocatable, save :: ptr(:)
type (ele_pointer_struct), allocatable, save :: eles(:)

integer i, ix1, ix2, ix_word, ios, ix, n_loc
real(rp) value
real(rp), pointer :: v(:)
character(*) word
character(40) attrib_name, ele_name
logical err_flag

! see if this is numeric

if (index('-+.0123456789', word(1:1)) /= 0) then
  read (word, *, iostat = ios) value
  if (ios /= 0) call parser_error ('BAD VARIABLE: ' // word)
  return
endif

! If not numeric...

ix_word = len_trim(word)
call verify_valid_name (word, ix_word)

! If word does not have a "[...]" then it must be a variable

ix1 = index(word, '[')
if (ix1 == 0) then   
  call find_indexx (word, bp_com%var%name, bp_com%var%indexx, bp_com%ivar_tot, i)
  if (i == 0) then
    call parser_error ('VARIABLE USED BUT NOT YET DEFINED: ' // word)
    value = 0
  else
    value = bp_com%var(i)%value
  endif
  return
endif

! Here if word does have a "[...]" then is a element attribute

ele_name = word(:ix1-1)    ! name of attribute

ix2 = index(word, ']')
attrib_name = word(ix1+1:ix2-1)

if (attrib_name == 'S' .and. bp_com%parser_name /= 'bmad_parser2') then
  call parser_error ('"S" ATTRIBUTE CAN ONLY BE USED WITH BMAD_PARSER2')
endif

! Apertures, etc.

select case (attrib_name)
case ('APERTURE', 'X_LIMIT', 'Y_LIMIT')
  call lat_ele_locator (ele_name, lat, eles, n_loc, err_flag)
  if (.not. err_flag .and. size(eles) > 0) then
    v => eles(1)%ele%value
    if (attrib_name == 'APERTURE') then
      if (v(x1_limit$) /= v(x2_limit$) .or. v(x1_limit$) /= v(y1_limit$) .or. &
          v(x1_limit$) /= v(y2_limit$)) then
        err_flag = .true.
      else
        value = v(x1_limit$)
      endif
    elseif (attrib_name == 'X_LIMIT') then
      if (v(x1_limit$) /= v(x2_limit$)) then
        err_flag = .true.
      else
        value = v(x1_limit$)
      endif
    elseif (attrib_name == 'Y1_LIMIT') then
      if (v(y1_limit$) /= v(y2_limit$)) then
        err_flag = .true.
      else
        value = v(y1_limit$)
      endif
    endif
  endif

  ! If size(eles) > 1 then there must be more than one element of the same name.

  if (size(eles) > 1) call parser_error (&
            'MULTIPLE ELEMENTS OF THE SAME NAME REFERENCED IN ATTRIBUTE: ' // word, warn_only = .true.)

! Everything else

case default
  call pointers_to_attribute (lat, ele_name, attrib_name, .false., ptr, err_flag, .false.)
  if (err_flag .or. size(ptr) == 0) then
    call parser_error('BAD ATTRIBUTE: ' // word)
  else
    value = ptr(1)%r
  endif

  ! If size(ptr) > 1 then there must be more than one element of the same name.

  if (size(ptr) > 1) call parser_error (&
            'MULTIPLE ELEMENTS OF THE SAME NAME REFERENCED IN ATTRIBUTE: ' // word, warn_only = .true.)

end select

end subroutine word_to_value

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_variable (word, lat)

implicit none

type (lat_struct) lat
character(*) word
character(1) delim
integer i, ivar, ixm, ixm2
logical delim_found, err_flag

!

call find_indexx (word, bp_com%var%name, bp_com%var%indexx, bp_com%ivar_tot, i)
if (i /= 0) then
  call parser_error ('VARIABLES ARE NOT ALLOWED TO BE REDEFINED: ' // word)
  call evaluate_value (word, bp_com%var(i)%value, lat, delim, delim_found, err_flag)
  return
endif

bp_com%ivar_tot = bp_com%ivar_tot + 1
if (bp_com%ivar_tot > size(bp_com%var%name)) call reallocate_bp_com_var()
ivar = bp_com%ivar_tot
bp_com%var(ivar)%name = word
call evaluate_value (bp_com%var(ivar)%name, bp_com%var(ivar)%value, &
                                     lat, delim, delim_found, err_flag)
if (delim_found .and. .not. err_flag) call parser_error  &
                  ('EXTRA CHARACTERS ON RHS: ' // bp_com%parse_line,  &
                   'FOR VARIABLE: ' // bp_com%var(ivar)%name)

! Reindex.

call find_indexx (word, bp_com%var%name, bp_com%var%indexx, ivar-1, ixm, ixm2)

bp_com%var(ixm2+1:ivar)%indexx = bp_com%var(ixm2:ivar-1)%indexx
bp_com%var(ixm2)%indexx = ivar

end subroutine parser_add_variable

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine bmad_parser_type_get (ele, attrib_name, delim, delim_found)

implicit none

type (ele_struct)  ele

integer ix, ix_word

character(*) attrib_name
character(40)  word
character(1)   delim, str_end
character(200) type_name

logical delim_found

!

call string_trim(bp_com%parse_line, bp_com%parse_line, ix)

str_end = bp_com%parse_line(1:1)
if (str_end == '"' .or. str_end == "'") then
  bp_com%parse_line = bp_com%parse_line(2:)
  ix = index(bp_com%parse_line, str_end)
  if (ix == 0) then
    call parser_error ('MISSING ENDING QUOTE MARK FOR TYPE = "attribute"',  &
                        'FOR ELEMENT: ' // ele%name)
    type_name = ' '
  else
    type_name = bp_com%parse_line(1:ix-1)
    bp_com%parse_line = bp_com%parse_line(ix+1:)
    call get_next_word (word, ix_word, ',=', delim, delim_found, .true.)
    if (ix_word /= 0) call parser_error ( &
              'EXTRA CHARACTERS FOUND AFTER TYPE ATTRIBUTE: ' // word,  &
              'FOR ELEMENT: ' // ele%name)
  endif
else
  call get_next_word (type_name, ix_word, ',= ', delim, delim_found, .false.)
endif

select case (attrib_name)
case ('TYPE')
  ele%type = type_name
case ('ALIAS')
  ele%alias = type_name
case ('DESCRIP', 'LATTICE')
  if (.not. associated(ele%descrip)) allocate (ele%descrip) 
  ele%descrip = type_name
case ('SR_WAKE_FILE') 
  call read_sr_wake (ele, type_name)
case ('LR_WAKE_FILE') 
  call read_lr_wake (ele, type_name)
case ('TO', 'REF_PATCH')
  ele%component_name = type_name
  call upcase_string (ele%component_name)
case ('CRYSTAL_TYPE')
  ele%component_name = type_name
case default
  call parser_error ('INTERNAL ERROR IN BMAD_PARSER_TYPE_GET: I NEED HELP!')
  if (bmad_status%exit_on_error) call err_exit
end select

end subroutine bmad_parser_type_get

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine read_lr_wake (ele, lr_file_name)
!
! Subroutine to read in a long-range wake field from an external file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   ele          -- Ele_struct: Element containing wake structure.
!   lr_file_name -- Character(*):  Name of long-range wake field file.
!
! Output:
!   ele -- Ele_struct: Element with wake information.
!     %rf_wake%lr(:)       -- Long-range wake potential.
!-
      
subroutine read_lr_wake (ele, lr_file_name)

implicit none

type lr_wake_input_struct
  real(rp) freq         ! Actual Frequency in Hz.
  real(rp) R_over_Q     ! Strength in V/C/m^2.
  real(rp) Q            ! Quality factor.
  integer m             ! Order (1 = dipole, 2 = quad, etc.)
  character(16) angle   ! polarization angle (radians/2pi).
  real(rp) b_sin, b_cos, a_sin, a_cos, t_ref
end type

type (ele_struct) ele
type (lr_wake_input_struct) lr(500)
integer n_row, iu, i, j, ios
character(*) lr_file_name
character(200) full_file_name
logical set_done
namelist / long_range_modes / lr

! Init

if (.not. associated(ele%rf_wake)) allocate (ele%rf_wake)
if (.not. associated(ele%rf_wake%sr_table))       allocate (ele%rf_wake%sr_table(0))
if (.not. associated(ele%rf_wake%sr_mode_long))  allocate (ele%rf_wake%sr_mode_long(0))
if (.not. associated(ele%rf_wake%sr_mode_trans)) allocate (ele%rf_wake%sr_mode_trans(0))
if (associated(ele%rf_wake%lr)) deallocate (ele%rf_wake%lr)

! get data

call find_this_file (iu, lr_file_name, full_file_name)
if (iu < 0) return

ele%rf_wake%lr_file = lr_file_name

lr%freq = -1
lr%angle = ''
lr%b_sin = 0
lr%b_cos = 0
lr%a_sin = 0
lr%a_cos = 0
lr%t_ref = 0

read (iu, nml = long_range_modes, iostat = ios)
close (iu)
if (ios /= 0) then
  call parser_error ('CANNOT READ LONG_RANGE_MODES NAMELIST FOR ELEMENT: ' // ele%name, & 
                'FROM FILE: '// full_file_name)
  return
endif

n_row = count(lr%freq /= -1)
allocate (ele%rf_wake%lr(n_row))
j = 0
do i = 1, size(lr)
  if (lr(i)%freq == -1) cycle

  j = j + 1
  ele%rf_wake%lr(j)%freq_in   = lr(i)%freq
  ele%rf_wake%lr(j)%freq      = lr(i)%freq
  ele%rf_wake%lr(j)%r_over_q  = lr(i)%r_over_q
  ele%rf_wake%lr(j)%q         = lr(i)%q
  ele%rf_wake%lr(j)%m         = lr(i)%m
  ele%rf_wake%lr(j)%b_sin     = lr(i)%b_sin
  ele%rf_wake%lr(j)%b_cos     = lr(i)%b_cos
  ele%rf_wake%lr(j)%a_sin     = lr(i)%a_sin
  ele%rf_wake%lr(j)%a_cos     = lr(i)%a_cos
  ele%rf_wake%lr(j)%t_ref     = lr(i)%t_ref

  call downcase_string(lr(i)%angle)
  if (lr(i)%angle == '') then
    call parser_error ('LONG_RANGE_MODE ANGLE IS MISSING. MUST BE NUMBER OR "UNPOLARIZED"', & 
                  'FOR ELEMENT: ' // ele%name, &
                  'IN FILE: ' // full_file_name)
    cycle
  endif

  if (index('unpolarized', trim(lr(j)%angle)) == 1) then
    ele%rf_wake%lr(j)%polarized = .false.
    ele%rf_wake%lr(j)%angle     = 0
  else
    ele%rf_wake%lr(j)%polarized = .true.
    read (lr(j)%angle, *, iostat = ios) ele%rf_wake%lr(j)%angle
    if (ios /= 0) then
      call parser_error ('BAD LONG_RANGE_MODE ANGLE.', &
                    'FOR ELEMENT: ' // ele%name, &
                    'IN FILE: ' // full_file_name)
      cycle
    endif
  endif
enddo

call randomize_lr_wake_frequencies (ele, set_done)
if (set_done) call bp_set_ran_status

end subroutine read_lr_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine read_sr_wake (ele, sr_file_name)
!
! Subroutine to read in a short-range wake field from an external file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   ele          -- Ele_struct: Element containing wake structure.
!   sr_file_name -- Character(*):  Name of short-range wake field file.
!
! Output:
!   ele -- Ele_struct: Element with wake information.
!     %rf_wake%sr_table(:)       -- Short-range wake potential.
!     %rf_wake%sr_mode_long(:)  -- Short-range wake potential.
!     %rf_wake%sr_mode_trans(:) -- Short-range wake potential.
!-

subroutine read_sr_wake (ele, sr_file_name)

implicit none

type (ele_struct) ele
type (rf_wake_sr_mode_struct) longitudinal(100), transverse(100)

real(rp) dz, z_max
real(rp), allocatable :: col1(:), col2(:), col3(:)
integer n_row, n, j, iu, ios, ix, i

character(*) sr_file_name
character(80) line
character(200) full_file_name

logical found_it

namelist / short_range_modes / z_max, longitudinal, transverse

! init

if (.not. associated(ele%rf_wake)) allocate (ele%rf_wake)
if (.not. associated(ele%rf_wake%lr)) allocate (ele%rf_wake%lr(0))
if (associated(ele%rf_wake%sr_table))       deallocate (ele%rf_wake%sr_table)
if (associated(ele%rf_wake%sr_mode_long))  deallocate (ele%rf_wake%sr_mode_long)
if (associated(ele%rf_wake%sr_mode_trans)) deallocate (ele%rf_wake%sr_mode_trans)

allocate (ele%rf_wake%sr_table(0), ele%rf_wake%sr_mode_long(0), ele%rf_wake%sr_mode_trans(0))

! get sr_table data

iu = 0
ele%rf_wake%sr_file = sr_file_name
call find_this_file (iu, sr_file_name, full_file_name)
if (iu < 0) return


! count number of lines in the file

i = 0
do
  read (iu, '(a)', iostat = ios) line
  if (ios < 0) then   ! end-of-file
    allocate (col1(i), col2(i), col3(i))
    exit
  endif
  i = i + 1
enddo
rewind (iu)

!

i = 0

do
  read (iu, '(a)', iostat = ios) line
  if (ios < 0) then   ! end-of-file
    close (iu)
    iu = 0
    exit
  endif
  if (ios > 0) then
    call parser_error ('ERROR READING WAKE FILE: ' // full_file_name)
    return
  endif
  call string_trim (line, line, ix)
  if (line(1:1) == '!') cycle  ! skip comments.
  if (ix == 0) cycle          ! skip blank lines.
  call str_upcase (line, line)
  if (line(1:) == 'END_SECTION') exit
  i = i + 1
  n_row = i
  read (line, *, iostat = ios) col1(i), col2(i), col3(i)

  if (ios /= 0) then
    call parser_error ('ERROR PARSING WAKE FILE: ' // full_file_name, &
                                        'CANNOT READ LINE: ' // line)
    return
  endif

enddo

allocate (ele%rf_wake%sr_table(0:n_row-1))
ele%rf_wake%sr_table%z     = col1(1:n_row)
ele%rf_wake%sr_table%long  = col2(1:n_row)
ele%rf_wake%sr_table%trans = col3(1:n_row)

deallocate (col1, col2, col3)

! err check

if (n_row > 1) then
  if (ele%rf_wake%sr_table(0)%z /= 0) then
    call parser_error ('WAKEFIELDS DO NOT START AT Z = 0!', &
                                  'IN FILE: ' // ele%rf_wake%sr_file)
    return
  endif

  n = n_row - 1
  dz = ele%rf_wake%sr_table(n)%z / n

  do j = 1, n
    if (abs(ele%rf_wake%sr_table(j)%z - dz * j) > 1e-4 * abs(dz)) then
      write (line, '(a, i5)') &
               'WAKEFIELD POINTS DO NOT HAVE UNIFORM DZ FOR POINT:', j
      call parser_error (line, 'IN FILE: ' // ele%rf_wake%sr_file)
      return
    endif
  enddo               

  ! if dz > 0 means that an old-style file is being used.

  if (dz > 0) call parser_error ( &
          'SHORT-RANGE WAKEFIELD FILE TABLES NOW MUST HAVE Z < 0! ' // full_file_name, &
          'REMEMBER THAT Wt NEEDS TO BE NEGATIVE ALSO!')

  if (ele%rf_wake%sr_table(1)%trans > 0) call parser_error ( &
           'POSITIVE Wt IN WAKEFIELD FILE INDICATES SIGN ERROR! ' // full_file_name)

endif

if (iu == 0) return  ! end of file reached

! Get sr_mode_long data

longitudinal(:)%phi = real_garbage$
transverse(:)%phi = real_garbage$
z_max = real_garbage$

read (iu, nml = short_range_modes, iostat = ios)
close (1)
if (ios /= 0) then
  call parser_error ('CANNOT READ SHORT_RANGE_MODES NAMELIST FROM FILE: ' & 
                    // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

n = count(longitudinal%phi /= real_garbage$)
allocate (ele%rf_wake%sr_mode_long(n))
ele%rf_wake%sr_mode_long = longitudinal(1:n)
if (any(longitudinal(1:n)%phi == real_garbage$)) call parser_error ( &
    'JUMBLED INDEX FOR LONGITUDINAL SHORT_RANGE_MODES FROM FILE: ' &
    // full_file_name, 'FOR ELEMENT: ' // ele%name)

n = count(transverse%phi /= real_garbage$)
allocate (ele%rf_wake%sr_mode_trans(n))
ele%rf_wake%sr_mode_trans = transverse(1:n)
if (any(transverse(1:n)%phi == real_garbage$)) call parser_error ( &
    'JUMBLED INDEX FOR TRANSVERSE SHORT_RANGE_MODES FROM FILE: ' &
    // full_file_name, 'FOR ELEMENT: ' // ele%name)


ele%rf_wake%z_sr_mode_max = z_max
if (z_max == real_garbage$) call parser_error ( &
    'Z_MAX NOT SET FOR SHORT_RANGE_MODES FROM FILE: ' &
    // full_file_name, 'FOR ELEMENT: ' // ele%name)

end subroutine read_sr_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine find_this_file (iu, file, full_file_name)
!
! Subroutine to open a file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   file -- Character(*): Name of wake field file
!
! Output:
!   iu             -- Integer: Open file unit. 
!                       If negative then there was an error. 
!   full_file_name -- Character(*): Full name with directory spec.
!-

subroutine find_this_file (iu, file, full_file_name)

implicit none

character(*) file, full_file_name
integer iu, ios, n
logical found_it

! open file

iu = lunget()
bp_com%dirs(2) = bp_com%calling_file%dir
call find_file (file, found_it, full_file_name, bp_com%dirs)
open (iu, file = full_file_name, status = 'OLD', action = 'READ', iostat = ios)
if (ios /= 0) then
call parser_error ('CANNOT OPEN WAKE FILE: ' // file)
iu = -1
return
endif

! If we have not read in this file before then add this to the list of files
! that are used to create the lattice.

n = size(bp_com%lat_file_names)
if (n < bp_com%num_lat_files + 1) call re_allocate (bp_com%lat_file_names, n + 100)
n = bp_com%num_lat_files
inquire (file = full_file_name, name = bp_com%lat_file_names(n+1))
if (all(bp_com%lat_file_names(n+1) /= bp_com%lat_file_names(1:n))) &
                                              bp_com%num_lat_files = n + 1
end subroutine find_this_file

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_overlay_group_names (ele, lat, pele, delim, delim_found)
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-
      
subroutine get_overlay_group_names (ele, lat, pele, delim, delim_found)

implicit none

type (ele_struct)  ele
type (parser_ele_struct) pele
type (lat_struct)  lat

real(rp) coef(200)
real(rp) value

integer ix_word, ixs, j, k
                           
character(1) delim
character(40) word_in, word
character(40) name(200), attrib_name(200)

logical delim_found, err_flag, end_of_file
                    
!

call get_next_word (word_in, ix_word, '{,}', delim, delim_found, .true.)
if (delim /= '{' .or. ix_word /= 0) call parser_error  &
        ('BAD ' // control_name(ele%lord_status) // 'SPEC: ' // word_in,  &
        'FOR ELEMENT: ' // ele%name)

! loop over all names in "{...}" list

do 

  call get_next_word (word_in, ix_word, '{,}/:', delim, delim_found, .true.)

  ! If "{}" with no slaves... 
  if (delim == '}' .and. ix_word == 0 .and. ele%n_slave == 0) then
    call get_next_word (word, ix_word, ',=:', delim, delim_found, .true.)
    exit
  endif

  ele%n_slave = ele%n_slave + 1
  ixs = ele%n_slave
  word = word_in

  j = index(word, '[')
  if (j > 1) then
    k = index(word, ']')
    if (k <= j+1) then
      call parser_error ('BAD ATTRIBUTE SPEC: ' // word_in, 'FOR: ' // ele%name)
      word = word(:k-1) // word(j+1:)
    else
      attrib_name(ixs) = word(j+1:k-1)
      word = word(:j-1) // word(k+1:)
    endif
  else
    attrib_name(ixs) = blank_name$
  endif

  name(ixs) = word

  if (delim == '/' .or. delim == ':') then
    call evaluate_value (trim(ele%name), value, lat, delim, delim_found, err_flag)
    if (err_flag) then
      call parser_error ('BAD COEFFICIENT: ' // word_in,  &
                                        'FOR ELEMENT: ' // ele%name)
      call load_parse_line ('normal', 1, end_of_file)         ! next line
      return
    endif
    coef(ixs) = value
  else
    coef(ixs) = 1.0
  endif

  if (delim == '}') then
    call get_next_word (word, ix_word, ',=:', delim, delim_found, .true.)
    exit
  elseif (delim /= ',') then
    call parser_error ('BAD ' // control_name(ele%lord_status) //  &
            'SPEC: ' // word_in, 'FOR: ' // ele%name)
    exit
  endif
                        
enddo

!

ixs = ele%n_slave

! if (ixs == 0) call parser_error ( &
!        'NO SLAVE ELEMENTS ASSOCIATED WITH GROUP/OVERLAY ELEMENT: ' // ele%name)

allocate (pele%coef(ixs), pele%name(ixs), pele%attrib_name(ixs))
pele%coef = coef(1:ixs)
pele%name = name(1:ixs)
pele%attrib_name = attrib_name(1:ixs)

end subroutine get_overlay_group_names

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-


subroutine verify_valid_name (name, ix_name, wildcards_permitted, integer_permitted)

implicit none

integer i, ix_name, ix1, ix2

character(*) name
character(27), parameter :: letters = '\ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
character(44), parameter :: valid_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\0123456789_[]().#'
character(1), parameter :: tab = achar(9)

logical, optional :: wildcards_permitted, integer_permitted
logical OK, wild_permit

! init

wild_permit = logic_option (.false., wildcards_permitted)

! check for blank spaces

do i = 1, min(ix_name, len(name))
  if (name(i:i) == ' ' .or. name(i:i) == tab) call parser_error  &
                        ('NO DELIMITER BETWEEN NAMES: ' // name)
enddo

! check for name too long

if (ix_name > len(name)) then
   call parser_error ('NAME TOO LONG: ' // name)
   ix_name = len(name)      ! chop name
endif

! check for name too short

if (ix_name == 0) call parser_error ('BLANK NAME')

! Check if integer

if (logic_option(.false., integer_permitted)) then
  if (is_integer(name)) return
endif

! check for invalid characters in name

OK = .true.
if (index(letters, name(1:1)) == 0 .and. &
          .not. (wild_permit .and. (name(1:1) == '*' .or. name(1:1) == '%'))) OK = .false.

do i = 1, ix_name
  if (index(valid_chars, name(i:i)) == 0 .and. &
         .not. (wild_permit .and. (name(i:i) == '*' .or. name(i:i) == '%'))) OK = .false.
enddo

if (.not. OK) call parser_error ('INVALID NAME: UNRECOGNIZED CHARACTERS IN: ' // name)

! Check for non-matched "(" ")" pairs

ix1 = index(name, '(')
ix2 = index(name, ')')
if (ix1 /= 0 .or. ix2 /= 0) then
  if (ix1 == 0) call parser_error ('UNMATCHED PARENTHESIS: ' // name)
  if (ix2 <= ix1+1) call parser_error  ('INVALID: REVERSED PARENTHESES: ' // name)
  if (index(name(ix1+1:), '(') /= 0 .or. index(name(ix2+1:), ')') /=  &
                 0) call parser_error ('INVALID: BAD PARENTHESES: ' // name)
endif

! Check for non matched "[" "]" pairs

ix1 = index(name, '[')
ix2 = index(name, ']')
if (ix1 /= 0 .or. ix2 /= 0) then
  if (ix1 == 0) call parser_error ('UNMATCHED BRACKET: ' // name)
  if (ix2 <= ix1+1) call parser_error  ('INVALID: REVERSED BRACKETS: ' // name)
  if (index(name(ix1+1:), '[') /= 0 .or. index(name(ix2+1:), ']') /=  &
                 0) call parser_error ('INVALID: BAD BRACKETS: ' // name)
  if (ix2 /= len(name)) then
    if (name(ix2+1:ix2+1) /= ' ') call parser_error  &
                  ('INVALID: SOMETHING AFTER CLOSING "]" BRACKET: ' // name)
  endif
endif

! check for more than 40 characters

if ((ix1 == 0 .and. ix_name > 40) .or. (ix1 > 41 .or. ix2 - ix1 > 41)) &
                          call parser_error ('NAME HAS > 40 CHARACTERS: ' // name)

end subroutine verify_valid_name

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_error (what1, what2, what3, seq, pele, stop_here, warn_only)

implicit none

type (seq_struct), optional :: seq
type (parser_ele_struct), optional :: pele

character(*) what1
character(*), optional :: what2, what3
character(160) lines(12)
character(16), parameter :: r_name = 'parser_error'
integer nl
logical, optional :: stop_here, warn_only

! bp_com%error_flag is a common logical used so program will stop at end of parsing

if (bmad_status%type_out) then

  nl = 0

  if (logic_option(.false., warn_only)) then
    nl=nl+1; lines(nl) = 'WARNING IN ' // trim(bp_com%parser_name) // ': ' // trim(what1)
  else
    nl=nl+1; lines(nl) = 'ERROR IN ' // trim(bp_com%parser_name) // ': ' // trim(what1)
  endif

  if (present(what2)) then
    nl=nl+1; lines(nl) = '     ' // trim(what2)
  endif

  if (present(what3)) then
    nl=nl+1; lines(nl) = '     ' // trim(what3)
  endif

  if (present(seq)) then
    nl=nl+1; lines(nl) = '      IN FILE: ' // trim(seq%file_name)
    nl=nl+1; write (lines(nl), '(a, i0)') '      AT LINE: ', seq%ix_line
  elseif (bp_com%current_file%full_name /= ' ') then
    if (bp_com%input_line_meaningful) then
      nl=nl+1; lines(nl) = '      IN FILE: ' // trim(bp_com%current_file%full_name)
      nl=nl+1; write (lines(nl), '(a, i0)') '      AT OR BEFORE LINE: ', bp_com%current_file%i_line
    else
      nl=nl+1; lines(nl) = '      ROOT FILE: ' // trim(bp_com%current_file%full_name)
    endif
  endif

  if (bp_com%input_line_meaningful) then
    if (len_trim(bp_com%input_line1) /= 0) then
      nl=nl+1; lines(nl) = '     ' // trim(bp_com%input_line1)
    endif
    if (len_trim(bp_com%input_line2) /= 0) then
      nl=nl+1; lines(nl) = '     ' // trim(bp_com%input_line2)
    endif
  endif

  if (present(pele)) then
    nl=nl+1; lines(nl) = '      ELEMENT DEFINED IN FILE: ' // trim(pele%lat_file)
    nl=nl+1; write (lines(nl), '(a, i0)') '      AT LINE: ', pele%ix_line_in_file
  endif

  nl=nl+1; lines(nl) = ''

  call out_io (s_error$, r_name, lines(1:nl))

endif

! Warnings do not result in bp_com%error_flag being set

if (.not. logic_option(.false., warn_only)) then
  bp_com%error_flag = .true.
  if (logic_option(.false., stop_here) .and. bmad_status%exit_on_error) stop
endif

end subroutine parser_error

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine init_bmad_parser_common

implicit none

integer nn, nt, i

! 

if (allocated(bp_com%var)) deallocate (bp_com%var)

nn = 17  ! number of "constant" variables
bp_com%ivar_init = nn + ubound(calc_method_name, 1) + &
           ubound(ref_orbit_name, 1) + ubound(element_end_name, 1) + ubound(aperture_type_name, 1) + &
           size(particle_name) + ubound(polarization_name, 1)
bp_com%ivar_tot = bp_com%ivar_init

nt = bp_com%ivar_tot
allocate (bp_com%var(nt))

bp_com%var( 1) = bp_var_struct('PI', pi, 0)
bp_com%var( 2) = bp_var_struct('TWOPI', twopi, 0)
bp_com%var( 3) = bp_var_struct('DEGRAD', 180 / pi, 0)
bp_com%var( 4) = bp_var_struct('RADDEG', pi / 180, 0)
bp_com%var( 5) = bp_var_struct('E_LOG', 2.718281828459_rp, 0)
bp_com%var( 6) = bp_var_struct('E_MASS', e_mass, 0)
bp_com%var( 7) = bp_var_struct('C_LIGHT', c_light, 0)
bp_com%var( 8) = bp_var_struct('M_ELECTRON', m_electron, 0)
bp_com%var( 9) = bp_var_struct('M_PROTON', m_proton, 0)
bp_com%var(10) = bp_var_struct('R_P', r_p, 0)
bp_com%var(11) = bp_var_struct('E_CHARGE', e_charge, 0)
bp_com%var(12) = bp_var_struct('EMASS', e_mass, 0)
bp_com%var(13) = bp_var_struct('CLIGHT', c_light, 0)
bp_com%var(14) = bp_var_struct('LINEAR_LATTICE', real(linear_lattice$, rp), 0)
bp_com%var(15) = bp_var_struct('CIRCULAR_LATTICE', real(circular_lattice$, rp), 0)
bp_com%var(16) = bp_var_struct('R_E', r_e, 0)
bp_com%var(17) = bp_var_struct('DEGREES', 180 / pi, 0)

do i = lbound(particle_name, 1), ubound(particle_name, 1)
  nn = nn + 1
  call str_upcase (bp_com%var(nn)%name, particle_name(i))
  bp_com%var(nn)%value = i
enddo

do i = 1, ubound(calc_method_name, 1)
  nn = nn + 1
  call str_upcase (bp_com%var(nn)%name, calc_method_name(i))
  bp_com%var(nn)%value = i
enddo

do i = 1, ubound(element_end_name, 1)
  nn = nn + 1
  call str_upcase (bp_com%var(nn)%name, element_end_name(i))
  bp_com%var(nn)%value = i
enddo

do i = 1, ubound(ref_orbit_name, 1)
  nn = nn + 1
  call str_upcase (bp_com%var(nn)%name, ref_orbit_name(i))
  bp_com%var(nn)%value = i
enddo

do i = 1, ubound(aperture_type_name, 1)
  nn = nn + 1
  call str_upcase (bp_com%var(nn)%name, aperture_type_name(i))
  bp_com%var(nn)%value = i
enddo

call indexx (bp_com%var(1:nt)%name, bp_com%var(1:nt)%indexx)

end subroutine init_bmad_parser_common

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine add_this_multipass (lat, m_slaves, lord_in)

use multipass_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: slave, lord, slave2, lord2, ele
type (ele_struct), optional :: lord_in
type (branch_struct), pointer :: branch
type (lat_ele_loc_struct) m_slaves(:)

integer i, j, k, n, i1, ix, ixc, ixic, ix_lord, ixb, ixb2, ix_n
integer n_multipass, ic, ix_l1, ix_l0, ix_pass, n_links

character(40) base_name

! Count slaves.
! If i > lat%n_ele_track we are looking at cloning a super_lord which should
! not happen.

n_multipass = size(m_slaves)

! setup multipass_lord

call new_control (lat, ix_lord)
lord => lat%ele(ix_lord)

if (present(lord_in)) then
  lord = lord_in   ! Use lord_in as template
else
  lord = pointer_to_ele(lat, m_slaves(1))  ! Set attributes equal to first slave.
endif

lord%lord_status = multipass_lord$
lord%n_slave = n_multipass
lord%ix1_slave = 0
lord%ix2_slave = -1
call add_lattice_control_structs (lat, lord)
if (lord%key == sbend$ .and. lord%ref_orbit == 0) lord%ref_orbit = single_ref$

! Setup bookkeeping between lord and slaves

do i = 1, n_multipass
  slave => pointer_to_ele (lat, m_slaves(i))
  ixc = i + lord%ix1_slave - 1
  lat%control(ixc)%ix_lord = ix_lord
  lat%control(ixc)%ix_slave = slave%ix_ele
  lat%control(ixc)%ix_branch = slave%ix_branch
  if (slave%n_lord /= 0) then
    call parser_error ('INTERNAL ERROR: CONFUSED MULTIPASS SETUP.', &
                  'PLEASE GET EXPERT HELP!')
    if (bmad_status%exit_on_error) call err_exit
  endif
  slave%n_lord = 1
  write (slave%name, '(2a, i1)') trim(slave%name), '\', i   ! '
  call add_lattice_control_structs (lat, slave)
  slave%slave_status = multipass_slave$
  ixic = slave%ic1_lord
  lat%ic(ixic) = ixc
  ! If slave is a super_lord then create the appropriate super_slave names
  i1 = 0
  do j = 1, slave%n_slave
    slave2 => pointer_to_slave(slave, j)
    if (slave2%n_lord == 1) then
      i1 = i1 + 1
      write (slave2%name, '(2a, i0, a, i0)') trim(lord%name), '#', i1, '\', i      ! '
    else
      slave2%name = ''
      do k = 1, slave2%n_lord
        lord2 => pointer_to_lord(slave2, k)
        lord2 => pointer_to_lord(lord2, 1)
        slave2%name = trim(slave2%name) // trim(lord2%name) // '\'     ! '
      enddo
      write (slave2%name, '(a, i0)') trim(slave2%name), i  
    endif
  enddo
enddo

! If this is a drift multipass whose multipass_slave elements are the result
! of splitting a drift with superposition then make sure that all split drift elements 
! of the lattice with the same base name have a name of the form "<base_name>#<n>" where
! <n> is an index from 1 for the first split drift.

ixb = index(lord%name, '#') - 1
if (lord%key == drift$ .and. ixb > 0) then
  ix_n = 0
  base_name = lord%name(1:ixb) 

  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    do j = 1, branch%n_ele_track
      ele => branch%ele(j)
      if (ele%key /= drift$) cycle
      ixb2 = index(ele%name, '#') - 1
      if (ixb2 /= ixb) cycle
      if (base_name(1:ixb) /= ele%name(1:ixb)) cycle
      ! super_slave drifts are temporary constructs that need to be ignored.
      ! This routine will be called later to correct the name of such elements.
      if (ele%slave_status == super_slave$) cycle 
      if (ele%slave_status == multipass_slave$) then
        call multipass_chain (ele, lat, ix_pass, n_links)
        if (ix_pass /= 1) cycle  ! Only do renaming once
        lord2 => pointer_to_lord(ele, 1)
        ix_n = ix_n + 1
        write (lord2%name, '(2a, i0)') base_name(1:ixb), '#', ix_n
        do k = 1, lord2%n_slave
          slave => pointer_to_slave(lord2, k)
          write (slave%name, '(2a, i0, a, i0)') base_name(1:ixb), '#', ix_n, '\', k      !'
        enddo
      else
        ix_n = ix_n + 1
        write (ele%name, '(2a, i0)') base_name(1:ixb), '#', ix_n
      endif
    enddo
  enddo

endif

!

call control_bookkeeper (lat, lord)

end subroutine add_this_multipass

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine reallocate_bp_com_var()

implicit none

type (bp_var_struct) :: var_temp(size(bp_com%var))

integer n

!

var_temp = bp_com%var

deallocate (bp_com%var)

n = bp_com%ivar_tot+200
allocate (bp_com%var(n))

n = size(var_temp)
bp_com%var(1:n) = var_temp


end subroutine reallocate_bp_com_var

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine add_all_superimpose (lat, super_ele_in, pele, in_lat)

use multipass_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct) super_ele_in
type (ele_struct), save :: super_ele_saved, super_ele
type (ele_struct), pointer :: ref_ele, ele, slave, lord, super_ele_out
type (parser_ele_struct) pele
type (multipass_all_info_struct) m_info
type (lat_struct), optional :: in_lat
type (lat_ele_loc_struct), allocatable :: m_slaves(:)
type (branch_struct), pointer :: branch

integer ix, i, j, k, it, nic, nn, i_ele, ib
integer n_inserted, n_con, i_br, ix_branch

character(40) matched_name(200), num, name
character(40), allocatable :: multi_name(:)
character(80) line

logical have_inserted, found, err_flag

! init

if (.not. bp_com%do_superimpose) return

call settable_dep_var_bookkeeping (super_ele_in)

call init_ele(super_ele_saved)
call init_ele(super_ele)

super_ele_saved = super_ele_in      ! in case super_ele_in changes
super_ele = super_ele_saved        ! 
n_inserted = 0

! If no refrence point then superposition is simple

if (pele%ref_name == blank_name$) then
  call compute_super_lord_s (lat, lat%ele(0), super_ele, pele)
  call add_superimpose (lat, super_ele, 0, err_flag, save_null_drift = .true., &
                                 create_em_field_slave = pele%create_em_field_slave)
  if (err_flag) bp_com%error_flag = .true.
  return
endif

! Insert ele in the lat.
! Do not insert twice at the same spot.

do i_br = 0, ubound(lat%branch, 1)
  branch => lat%branch(i_br)
  branch%ele%old_is_on = .false.    ! to keep track of where we have inserted
  branch%ele%ix_pointer = -1
enddo


do 

  have_inserted = .false.

  do i_br = 0, ubound(lat%branch, 1)
    branch => lat%branch(i_br)

    ele_loop: do i_ele = 1, branch%n_ele_max

      ref_ele => branch%ele(i_ele)
       
      if (ref_ele%lord_status == group_lord$ .or. ref_ele%slave_status == super_slave$) cycle
      if (ref_ele%lord_status == girder_lord$) cycle
      if (ref_ele%old_is_on) cycle
      if (.not. match_wild(ref_ele%name, pele%ref_name)) cycle

      do i = 1, n_inserted
        if (ref_ele%name == matched_name(i)) cycle ele_loop
      enddo
     
      ref_ele%old_is_on = .true.

      ! If superimposing on a multipass_lord then the superposition
      ! must be done at all multipass locations.

      if (ref_ele%lord_status == multipass_lord$) then
        allocate (m_slaves(ref_ele%n_slave), multi_name(ref_ele%n_slave))
        do i = 1, ref_ele%n_slave
          slave => pointer_to_slave(ref_ele, i)
          slave%ix_pointer = i  ! tag ref element
        enddo
        ix_branch = slave%ix_branch  ! Branch of slaves
        branch => lat%branch(ix_branch)
        call string_trim(super_ele_saved%name, super_ele_saved%name, ix)
        super_ele%name = super_ele_saved%name(:ix)

        ! Put in the superposition at the multipass locations.
        ! Since elements get shuffled around, tag the superimposed elements 
        !     with "temp_name!" to identify them later.

        i = 0  ! Element index
        j = 0  ! Number of superpositions done.
        do while (i < branch%n_ele_max)
          i = i + 1
          if (branch%ele(i)%ix_pointer /= j+1) cycle
          j = j + 1
          call compute_super_lord_s (lat, branch%ele(i), super_ele, pele)
          ! Don't need to save drifts since a multipass_lord drift already exists.
          call add_superimpose (lat, super_ele, ix_branch, err_flag, super_ele_out, &
                        save_null_drift = .false., create_em_field_slave = pele%create_em_field_slave)
          if (err_flag) bp_com%error_flag = .true.
          super_ele_out%name = 'temp_name!'
        enddo

        ! Mark any super_lord drifts (along with any lords that control the super_lord) for 
        ! future deletion at the end of lattice parsing.

        call multipass_all_info (lat, m_info) ! Save multipass info for later.

        do i = lat%n_ele_track+1, lat%n_ele_max 
          ele => lat%ele(i)
          if (ele%key /= drift$) cycle
          if (ele%lord_status /= super_lord$) cycle 
          ele%key = null_ele$ ! mark for deletion

          do j = 1, ele%n_lord
            lord => pointer_to_lord(ele, j)
            lord%key = null_ele$  ! Mark lord for deletion
          enddo

          ! Need to remove super_lord/super_slave links otherwise the code below gets confused
          ! when it tries to connect the former super_slave drifts.
          do while (ele%n_slave /= 0)
            call remove_lord_slave_link (ele, pointer_to_slave(ele, 1))
          enddo

        enddo

        ! Add a multipass_lord to control the created super_lords.
        ! If the super_lords have a single super_slave and the super_slave
        ! has only a single super_lord, the super_lords
        ! can be eliminated and the created multipass_lord can control the
        ! super_slaves directly.

        j = 0

        do i = 1, branch%n_ele_track
          if (branch%ele(i)%name == 'temp_name!') then
            branch%ele(i)%name = super_ele_saved%name
            j = j + 1
            m_slaves(j) = ele_to_lat_loc (branch%ele(i))
          endif
        enddo

        do i = lat%n_ele_track+1, lat%n_ele_max
          if (lat%ele(i)%name == 'temp_name!') then
            lat%ele(i)%name = super_ele_saved%name
            j = j + 1
            m_slaves(j) = ele_to_lat_loc (lat%ele(i))
          endif
        enddo

        ele => pointer_to_ele (lat, m_slaves(1))
        if (ele%lord_status == super_lord$ .and. ele%n_slave == 1) then
          slave => pointer_to_slave(ele, 1)
          if (slave%n_lord == 1) then
            do i = 1, size(m_slaves)
              ele => pointer_to_ele (lat, m_slaves(i))
              ele%key = -1 ! Mark for deletion
              ele => pointer_to_slave(ele, 1)
              ele%name = super_ele_saved%name
              m_slaves(i) = ele_to_lat_loc (ele)
            enddo
            call remove_eles_from_lat (lat, .false.)
          endif
        endif

        call add_this_multipass (lat, m_slaves, super_ele_saved) 

        ! Reconnect drifts that were part of the multipass region.

        do i = 1, size(m_info%top)
          do j = 1, size(m_info%top(i)%slave, 2)
            slave => m_info%top(i)%slave(1, j)%ele
            if (slave%key /= drift$) cycle
            if (slave%slave_status == multipass_slave$) cycle
            do k = 1, size(m_info%top(i)%slave(:, j))
              ele => m_info%top(i)%slave(k, j)%ele
              m_slaves(k) = ele_to_lat_loc (ele)
              ib = index(ele%name, '\') ! '
              if (ib /= 0) ele%name = ele%name(1:ib-1) // ele%name(ib+2:)
            enddo
            call add_this_multipass (lat, m_slaves)
          enddo
        enddo

        call deallocate_multipass_all_info_struct (m_info)

      ! Else not superimposing on a multipass_lord ...

      else
        call compute_super_lord_s (lat, branch%ele(i_ele), super_ele, pele)
        call string_trim(super_ele_saved%name, super_ele_saved%name, ix)
        super_ele%name = super_ele_saved%name(:ix)            
        call add_superimpose (lat, super_ele, i_br, err_flag, super_ele_out, &
                    save_null_drift = .true., create_em_field_slave = pele%create_em_field_slave)
        if (err_flag) bp_com%error_flag = .true.
        call control_bookkeeper (lat, super_ele_out)
      endif

      call s_calc (lat)

      n_inserted = n_inserted + 1
      matched_name(n_inserted) = super_ele%name
      have_inserted = .true.   

    enddo ele_loop

  enddo

  if (.not. have_inserted) exit

enddo

! Error check. If the reference element has been defined but not used in the lattice
! then this is not an error

if (n_inserted == 0) then
  do i = 1, in_lat%n_ele_max
    found = match_wild(in_lat%ele(i)%name, pele%ref_name)
    if (found) exit
  enddo
  if (.not. found) call parser_error ('NO MATCH FOR REFERENCE ELEMENT: ' //  &
          pele%ref_name, 'FOR SUPERPOSITION OF: ' // super_ele_saved%name, pele = pele)
endif

end subroutine add_all_superimpose

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine compute_super_lord_s (lat, ref_ele, super_ele, pele)

implicit none

type (lat_struct), target :: lat
type (ele_struct) ref_ele, super_ele
type (ele_struct), pointer :: slave
type (parser_ele_struct) pele
type (branch_struct), pointer :: branch

integer i, ix, ct

real(rp) s_ref_begin, s_ref_end

! Find the reference point on the element being superimposed.

super_ele%s = pele%s

if (pele%ele_pt == begin$) then
  super_ele%s = super_ele%s + super_ele%value(l$)
elseif (pele%ele_pt == center$ .or. pele%ele_pt == not_set$) then
  super_ele%s = super_ele%s + super_ele%value(l$) / 2
elseif (pele%ele_pt /= end$) then
  call parser_error ('ERROR IN COMPUTE_SUPER_LORD_S: CONTROL #1 INTERNAL ERROR!')
  if (bmad_status%exit_on_error) call err_exit
endif

! Find the refernce point in the lattice.

ct = ref_ele%lord_status
if (ct == overlay_lord$ .or. ct == girder_lord$) then
  s_ref_begin = 1e10
  s_ref_end = 0
  do i = 1, ref_ele%n_slave
    slave => pointer_to_slave(ref_ele, i)
    s_ref_begin = min(s_ref_begin, slave%s - slave%value(l$))
    s_ref_end = max(s_ref_end, slave%s)
  enddo
elseif (ct == group_lord$) then
  call parser_error ('SUPERPOSING: ' // super_ele%name, 'UPON GROUP' // pele%ref_name)
  return
else
  s_ref_begin = ref_ele%s - ref_ele%value(l$)
  s_ref_end = ref_ele%s
endif

! Now compute the s position at the end of the element and put it in ele%s.

if (pele%ref_pt == begin$) then
  super_ele%s = super_ele%s + s_ref_begin
elseif (pele%ref_pt == center$ .or. pele%ref_pt == not_set$) then
  super_ele%s = super_ele%s + (s_ref_begin + s_ref_end) / 2
elseif (pele%ref_pt == end$) then
  super_ele%s = super_ele%s + s_ref_end
else
  call parser_error ('ERROR IN COMPUTE_SUPER_LORD_S: CONTROL #2 INTERNAL ERROR!')
  if (bmad_status%exit_on_error) call err_exit
endif

! For circular lattices a superimpose can wrap around the beginning or 
! the end of the lattice.

branch => lat%branch(ref_ele%ix_branch)
if (branch%param%lattice_type == circular_lattice$) then
  if (super_ele%s > branch%ele(branch%n_ele_track)%s) then
    super_ele%s = super_ele%s - branch%param%total_length
  elseif (super_ele%s < 0) then
    super_ele%s = super_ele%s + branch%param%total_length
  endif
endif

end subroutine compute_super_lord_s

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_sequence_args (seq_name, arg_list, delim, err_flag)
!
! Subroutine to get the argument list for a replacement_line or a list.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine get_sequence_args (seq_name, arg_list, delim, err_flag)

implicit none

integer n_arg, ix_word

character(*), pointer :: arg_list(:)
character(*) seq_name
character(1) delim
character(40) name(20), word

logical delim_found, err_flag

!

n_arg = 0
err_flag = .true.

do
  call get_next_word (word, ix_word, '(): =,', delim, delim_found, .true.)
  if (ix_word == 0 .or. delim == '( :=') then
    call parser_error ('BAD ARGUMENT LIST FOR: ', seq_name)
    return
  endif
  n_arg = n_arg + 1
  name(n_arg) = word
  if (delim == ')') exit
enddo

err_flag = .false.
allocate (arg_list(n_arg))
arg_list = name(1:n_arg)

end subroutine get_sequence_args

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine seq_expand1 (sequence, iseq_tot, lat, top_level)
!
! Subroutine to expand a sequence.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

recursive subroutine seq_expand1 (sequence, iseq_tot, lat, top_level)

implicit none

type (seq_struct), target :: sequence(:)
type (seq_struct), pointer :: seq, sub_seq
type (seq_ele_struct), pointer :: s_ele(:)
type (seq_ele_struct), pointer :: s_ele2(:)
type (seq_ele_struct), pointer :: this_ele
type (lat_struct) lat

integer ix_ele, iseq_tot, ix_word, ix, i, n
integer, save :: ix_internal = 0

real(rp) rcount

character(40) word
character(1) delim, c_delim
character(40) str, name
character(n_parse_line) parse_line_saved

logical delim_found, replacement_line_here, c_delim_found
logical err_flag, top_level

! init

allocate (s_ele(ubound(lat%ele, 1)))
s_ele%type = 0
s_ele%ix_ele = 0
s_ele%ix_arg = 0

! save info on what file we are parsing for error messages.

seq => sequence(iseq_tot)
seq%file_name = bp_com%current_file%full_name 
seq%ix_line = bp_com%current_file%i_line

! first thing should be a "("

call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)

if (delim /= '(') call parser_error  &
      ('EXPECTING "(", GOT: ' // delim, 'FOR LINE: ' // seq%name)
if (ix_word /= 0)  call parser_error  &
      ('EXTRANEOUS STUFF BEFORE "(", GOT: ' // word,  &
      'FOR LINE: ' // seq%name)

! now parse list proper

ix_ele = 1
this_ele => s_ele(ix_ele)

do 

  call get_next_word (word, ix_word, ':=(,)[]', delim, delim_found, .true.)

  ix = index(word, '*')          ! E.g. word = '-3*LINE'
  if (ix /= 0) then
    ! Evaluate the rep count. The ":" is to prevent confusion if the rep count
    ! has an ending like "," that would signal line continuation.
    parse_line_saved = bp_com%parse_line
    bp_com%parse_line = word(:ix-1) // ":" 
    call evaluate_value (trim(seq%name) // ' Repetition Count', rcount, &
                            lat, c_delim, c_delim_found, err_flag)
    this_ele%rep_count = nint(rcount)
    if (err_flag) return
    if (bp_com%parse_line /= '' .or. c_delim /= ":") then
      call parser_error ('MALFORMED REPETION COUNT FOUND IN SEQUENCE: ' // seq%name)
      return
    endif
    bp_com%parse_line = parse_line_saved
    this_ele%name = word(ix+1:)
    if (this_ele%rep_count < 0) then
      this_ele%reflect = .true.
    else
      this_ele%reflect = .false.
    endif
    this_ele%rep_count = abs(this_ele%rep_count)
    ix_word = ix_word - ix
  elseif (word(1:1) == '-') then
    this_ele%reflect = .true.
    this_ele%rep_count = 1
    this_ele%name = word(2:)
    ix_word = ix_word - 1
  else
    this_ele%reflect = .false.
    this_ele%rep_count = 1
    this_ele%name = word
  endif

  name = this_ele%name
  if (name /= ' ') call verify_valid_name (name, ix_word)

  ! Check for line tag

  if (delim == '[') then
    call get_next_word (this_ele%tag, ix_word, '[]:=(,)', delim, delim_found, .true.)
    if (delim /= ']') then
      call parser_error ('NO MATCHING "]" FOUND FOR OPENING "[" IN SEQUENCE: ' // seq%name)
      return
    endif
    call get_next_word (word, ix_word, '[]:=(,)', delim, delim_found, .true.)
    if (ix_word > 0) then
      call parser_error ('ILLEGAL CHARACTERS AFTER CLOSING "]" FOUND IN SEQUENCE: ' // seq%name)
      return
    endif   
  endif

  ! Check for a subline or replacement line.
  ! If there is one then save as an internal sequence.

  replacement_line_here = .false.

  if (delim == '(') then ! subline or replacement line

    ! if a subline...
    if (name == ' ') then  
      ix_internal = ix_internal + 1
      write (str, '(a, i3.3)') '#Internal', ix_internal   ! unique name 
      this_ele%name = str
      iseq_tot = iseq_tot + 1
      sub_seq => sequence(iseq_tot) 
      sub_seq%name = str
      sub_seq%type = seq%type
      sub_seq%multipass = seq%multipass
      if (sub_seq%type == replacement_line$) then
        ix = size (seq%dummy_arg)
        allocate (sub_seq%dummy_arg(ix), &
              sub_seq%corresponding_actual_arg(ix), this_ele%actual_arg(ix))
        sub_seq%dummy_arg = seq%dummy_arg
        this_ele%actual_arg = seq%dummy_arg
      endif
      bp_com%parse_line = '(' // bp_com%parse_line
      call seq_expand1 (sequence, iseq_tot, lat, .false.)

    ! else this is a replacement line
    else    
      replacement_line_here = .true.
      call get_sequence_args (name, this_ele%actual_arg, delim, err_flag)
      if (err_flag) return
    endif

    call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)
    if (word /= ' ') call parser_error &
              ('NO COMMA AFTER SUBLINE OR REPLACEMENT LINE. FOUND: ' // &
               word, 'IN THE SEQUENCE: ' // seq%name)
  endif

  if (this_ele%name == ' ') call parser_error &
            ('SUB-ELEMENT NAME IS BLANK FOR LINE/LIST: ' // seq%name)

  ! if a replacement line then look for element in argument list

  this_ele%ix_arg = 0
  if (seq%type == replacement_line$) then
    do i = 1, size(seq%dummy_arg)
      if (seq%dummy_arg(i) == this_ele%name) then
        this_ele%ix_arg = i
        exit
      endif
    enddo
  endif

! 

  n = size(s_ele)
  ix_ele = ix_ele + 1

  if (ix_ele > n) then
    allocate (s_ele2(n))      
    s_ele2 = s_ele(1:n)
    deallocate (s_ele)
    allocate (s_ele(n+1000))
    s_ele(1:n) = s_ele2
    deallocate(s_ele2)
  endif

  this_ele => s_ele(ix_ele)

  if (delim == ')') exit

  if (delim /= ',') then
    call parser_error ('EXPECTING "," GOT: ' // delim, 'FOR LINE: ' // seq%name)
    exit
  endif
         
enddo

! make sure there is nothing else if at top level

if (top_level) then
  call get_next_word(word, ix_word, ':=() ', delim, delim_found, .true.)
  if (delim_found .or. ix_word /= 0) call parser_error  &
        ('EXTRA CHARACTERS AFTER CLOSING ")"',  'FOR LINE: ' // seq%name)
endif

! transfer

ix_ele = ix_ele - 1
allocate (seq%ele(ix_ele))

do i = 1, ix_ele
  seq%ele(i) = s_ele(i)
enddo

deallocate (s_ele)

end subroutine seq_expand1

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine allocate_plat (plat, n_ele_max) 
!
! Subroutine to allocate allocatable array sizes.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

Subroutine allocate_plat (plat, n_ele_max)

implicit none

type (parser_lat_struct) plat
type (parser_ele_struct), pointer :: temp_pele(:)

integer i, n_now, n_ele_max

! assume all the arrays have the same size

if (associated(plat%ele)) then
  n_now = ubound(plat%ele, 1)
  allocate (temp_pele(0:n_now))
  temp_pele = plat%ele
  deallocate (plat%ele)
  allocate (plat%ele(0:n_ele_max))
  plat%ele(0:n_now) = temp_pele
  deallocate (temp_pele)

else
  allocate (plat%ele(0:n_ele_max))
  n_now = -1
endif

! %ixx is used as a pointer from the in_lat%ele array to the plat%ele array

do i = n_now+1, ubound(plat%ele, 1)
  nullify (plat%ele(i)%name)

  plat%ele(i)%ref_name = blank_name$
  plat%ele(i)%ref_pt  = not_set$
  plat%ele(i)%ele_pt  = not_set$
  plat%ele(i)%s       = 0
  plat%ele(i)%create_em_field_slave = .false.
enddo

end subroutine allocate_plat

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_add_lord (in_lat, n2, plat, lat)
!
! Subroutine to add overlay, group, and girder lords to the lattice.
! For overlays and groups: If multiple elements have the same name then 
! use all of them.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_lord (in_lat, n2, plat, lat)

implicit none

type (lat_struct), target :: in_lat, lat
type (ele_struct), pointer :: lord, lord2, slave, ele
type (parser_lat_struct), target :: plat
type (parser_ele_struct), pointer :: pele
type (control_struct), pointer, save :: cs(:) => null()
type (branch_struct), pointer :: branch

integer i, ic, n, n2, k, k2, ix, j, ib, ie, n_list, ix2, ns, ixs, ii, ix_end
integer ix_lord, k_slave, ix_ele_now, ix_girder_end, ix_super_lord_end
integer, allocatable :: r_indexx(:), ix_ele(:), ix_branch(:)

character(40), allocatable :: name_list(:)
character(40) name, input_slave_name, attrib_name, missing_slave_name

logical err, slave_not_in_lat, found_ele_match, created_girder_lord

! Setup...
! in_lat has the lords that are to be added to lat.
! We add an extra 1000 places to the arrays to give us some overhead.

n = n2 + 1000
do i = 0, ubound(lat%branch, 1)
    n = n + lat%branch(i)%n_ele_max
enddo

allocate (name_list(n), ix_ele(n), ix_branch(n), r_indexx(n))
allocate (cs(1000))

n_list = 0
do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  n = branch%n_ele_max
  ix2 = n_list + n
  name_list(n_list+1:ix2) = branch%ele(1:n)%name
  ix_ele(n_list+1:ix2)    = branch%ele(1:n)%ix_ele
  ix_branch(n_list+1:ix2) = branch%ele(1:n)%ix_branch
  n_list = ix2
enddo

call indexx (name_list(1:n_list), r_indexx(1:n_list)) ! get sorted list

! loop over elements

main_loop: do n = 1, n2

  lord => in_lat%ele(n)  ! next lord to add
  pele => plat%ele(lord%ixx)
  
  !-----------------------------------------------------
  ! overlay and groups

  select case (lord%lord_status)
  case (overlay_lord$, group_lord$)
 
    call new_control (lat, ix_lord)  ! get index in lat where lord goes
    lat%ele(ix_lord) = lord

    ! Find where the slave elements are. 
    ! If a slave element is not in lat but is in in_lat then the slave has 
    ! not yet been used in the lattice list. 
    ! In this case do not add the lord to the lattice.

    j = 0 ! number of slaves found
    slave_not_in_lat = .false.  ! Is there a slave that is not in the lattice?

    do i = 1, lord%n_slave

      name = pele%name(i)
      call find_indexx (name, name_list, r_indexx, n_list, k, k2)

      if (k == 0) then
        slave_not_in_lat = .true.
        missing_slave_name = name
      endif

      if ((k == 0 .and. j > 0) .or. (k > 0 .and. slave_not_in_lat) .or. &
          (k == 0 .and. all(in_lat%ele(1:n2)%name /= name))) then
        call parser_error ('CANNOT FIND SLAVE FOR: ' // lord%name, &
                      'CANNOT FIND: '// missing_slave_name, pele = pele)
        lat%n_ele_max = lat%n_ele_max - 1 ! Undo new_control call
        cycle main_loop
      endif

      if (k == 0) cycle

      ! There might be more than 1 element with %name = name. 
      ! Loop over all elements whose name matches name.
      ! Put the info into the cs structure.

      do 
        j = j + 1
        k = r_indexx(k2)
        cs(j)%coef = pele%coef(i)
        cs(j)%ix_slave = ix_ele(k)
        cs(j)%ix_branch = ix_branch(k)
        cs(j)%ix_lord = -1             ! dummy value
        attrib_name = pele%attrib_name(i)
        if (attrib_name == blank_name$) attrib_name = lord%component_name
        slave => pointer_to_ele (lat, ix_ele(k), ix_branch(k))
        ix = attribute_index(slave, attrib_name)
        ! If attribute not found it may be a special attribute like accordian_edge$.
        ! A special attribute will have ix > n_attrib_maxx
        if (ix < 1 .and. lord%lord_status == group_lord$) then
          ix = attribute_index(lord, attrib_name)
          if (ix <= n_attrib_maxx) ix = 0  ! Mark as not valid
        endif
        cs(j)%ix_attrib = ix
        if (ix < 1) then
          call parser_error ('IN OVERLAY OR GROUP ELEMENT: ' // lord%name, &
                        'ATTRIBUTE: ' // attrib_name, &
                        'IS NOT A VALID ATTRIBUTE OF: ' // slave%name, &
                        pele = pele)
          cycle main_loop
        endif
        k2 = k2 + 1
        if (k2 > n_list) exit
        k = r_indexx(k2)
        slave => pointer_to_ele (lat, ix_ele(k), ix_branch(k))
        ! exit loop if no more matches
        if (slave%name /= name) exit 
      enddo

    enddo

    lord%n_slave = j

    ! If the lord has no slaves then discard it

    if (j == 0) then
      lat%n_ele_max = lat%n_ele_max - 1 ! Undo new_control call
      cycle main_loop
    endif

    ! put the element name in the list r_indexx list

    call find_indexx (lord%name, name_list, r_indexx, n_list, k, add_to_list = .true.)
    n_list = n_list + 1
    ix_ele(k) = ix_lord
    ix_branch(k) = 0

    do ii = 1, n_list-1
      if (name_list(r_indexx(ii)) > name_list(r_indexx(ii+1))) then
        call parser_error ('PARER_ADD_LORD INTERNAL ERROR!')
        if (bmad_status%exit_on_error) call err_exit
      endif
    enddo

    ! create the lord

    ns = lord%n_slave

    select case (lord%lord_status)
    case (overlay_lord$)
      call create_overlay (lat, ix_lord, lord%component_name, cs(1:ns), err)
    case (group_lord$)
      call create_group (lat, ix_lord, cs(1:ns), err)
    end select
    if (err) call parser_error ('ELEMENT OR GROUP: ' // lord%name, &
                           'IS TRYING TO CONTROL AN ATTRIBUTE THAT IS NOT FREE TO VARY!', &
                           pele = pele)

  !-----------------------------------------------------
  ! girder
  ! Create an girder element for each lattice section that matches the slave list names.
  ! If no lattice element names match any names on the girder slave list then assume 
  ! this girder is for a different lattice and ignore the girder. If some do match and 
  ! some don't then flag this as an error.

  case (girder_lord$)

    ! Loop over all elements in the lattice.

    created_girder_lord = .false.

    branch_loop: do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      ix_girder_end = -1

      ele_loop: do ie = 1, branch%n_ele_track

        if (ie <= ix_girder_end) cycle

        ! Loop over girder slave list and see if this section matches.

        ix_ele_now = ie
        ix_super_lord_end = -1   ! Not in a super_lord
        ixs = 1       ! Index of girder slave element we are looking for

        slave_loop: do            ! loop over all girder slaves

          if (ix_ele_now > branch%n_ele_track) then
            if (branch%param%lattice_type == linear_lattice$) cycle ele_loop
            ix_ele_now = ix_ele_now - branch%n_ele_track
          endif

          if (ixs > lord%n_slave) exit
          ele => pointer_to_ele (lat, ix_ele_now, ib)
          input_slave_name = pele%name(ixs)

          if (ele%key == drift$) then
            ix_ele_now = ix_ele_now + 1
            cycle
          endif

          ! If a super_slave then go to the lords to match the name.
          ! The logic here is complicated by the fact that elements can, for example,
          ! be completely contained within another element.

          if (ele%slave_status == super_slave$) then
            do ic = 1, ele%n_lord
              lord2 => pointer_to_lord(ele, ic)
              slave => pointer_to_slave(lord2, lord2%n_slave)
              ix_end = slave%ix_ele
              if (lord2%slave_status == multipass_slave$) lord2 => pointer_to_lord(lord2, 1)
              if (match_wild(lord2%name, input_slave_name)) then
                cs(ixs)%ix_slave  = lord2%ix_ele
                cs(ixs)%ix_branch = lord2%ix_branch
                ixs = ixs + 1  ! Next girder slave
                ix_super_lord_end = ix_far_index(ix_ele_now, ix_super_lord_end, ix_end)
                cycle slave_loop
              endif
            enddo
            ! Here if no match. 
            ! If not previously in a super lord then there is no overall match to the slave list.
            if (ix_super_lord_end == -1) cycle ele_loop
            ix_ele_now = ix_ele_now + 1
            cycle
          endif

          ix_super_lord_end = -1  ! Not in a super_lord

          ! If a multipass_slave

          if (ele%slave_status == multipass_slave$) then
            lord2 => pointer_to_lord(ele, 1)
            if (match_wild(lord2%name, input_slave_name)) then
              cs(ixs)%ix_slave  = lord2%ix_ele
              cs(ixs)%ix_branch = lord2%ix_branch
              ixs = ixs + 1  ! Next girder slave
              ix_ele_now = ix_ele_now + 1
              cycle
            endif
            ! No match to the slave list 
            cycle ele_loop
          endif

          ! Regular element

          if (match_wild(ele%name, input_slave_name)) then
            cs(ixs)%ix_slave  = ele%ix_ele
            cs(ixs)%ix_branch = ele%ix_branch
            ixs = ixs + 1  ! Next girder slave
            ix_ele_now = ix_ele_now + 1
            cycle
          endif

          
          ! No match to the slave list. If a marker then ignore

          if (ele%key == marker$) then
            ix_ele_now = ix_ele_now + 1
            cycle
          endif

          ! Match failed. Start again 

          cycle ele_loop

        enddo slave_loop

        ! create the girder element

        call new_control (lat, ix_lord)
        call create_girder (lat, ix_lord, cs(1:lord%n_slave), lord)
        created_girder_lord = .true.

        ix_girder_end = cs(lord%n_slave)%ix_slave
        if (ix_girder_end < ie) exit ele_loop  ! And on to next branch

      enddo ele_loop

    enddo branch_loop

    if (.not. created_girder_lord) then
      call parser_error ('FAILED TO FIND REGION IN LATTICE FOR CRATING GIRDER: ' // lord%name, warn_only = .true.)
    endif

  end select

enddo main_loop

! cleanup

deallocate (r_indexx)
deallocate (name_list)
deallocate (cs)

!-------------------------------------------------------------------------
contains

! Function to return the index that is farthest (reached last) when moving
! from ix_now in a positive direction. Tricky part is if there is wrap around.
! An index of -1 means that the corresponding point does not exist.

function ix_far_index (ix_now, ix1, ix2) result (ix_far)

implicit none

integer ix_now, ix1, ix2, ix_far

!

if (ix1 < 0) then      ! Point 1 does not exist so point 2 wins by default
  ix_far = ix2
elseif (ix2 < 0) then  ! Point 2 does not exist so point 1 wins by default
  ix_far = ix1
elseif (ix1 < ix_now .and. ix2 < ix_now) then  ! both wrapped case
  ix_far = max(ix1, ix2)
elseif (ix1 > ix_now .and. ix2 > ix_now) then ! No wrap ("normal") case
  ix_far = max(ix1, ix2)
else                      ! One is wrapped but not the other case
  ix_far = min(ix1, ix2)
endif

end function ix_far_index

end subroutine parser_add_lord

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_set_ele_defaults (ele)
!
! Subroutine initialize an element given its key (key).
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_set_ele_defaults (ele)

type (ele_struct) ele

select case (ele%key)

case (bend_sol_quad$) 
  ele%mat6_calc_method = symp_lie_bmad$
  ele%tracking_method  = symp_lie_bmad$

case (branch$, photon_branch$)
  ele%value(direction$) = 1
  ele%value(particle$) = real_garbage$
  ele%value(lattice_type$) = linear_lattice$

case (crystal$)
  ele%value(follow_diffracted_beam$) = 1  ! True
  ele%value(ref_polarization$) = sigma_polarization$ 

case (custom$)  
  ele%mat6_calc_method = custom$
  ele%tracking_method  = custom$
  ele%field_calc       = custom$

case (ecollimator$)
  ele%offset_moves_aperture = .true.
  ele%aperture_type = elliptical$

case (lcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_scale$) = 1
  ele%value(n_cell$) = 1

case (multilayer_mirror$)
  ele%value(ref_polarization$) = sigma_polarization$  

case (rbend$, sbend$)
  ele%value(fintx$) = real_garbage$
  ele%value(hgapx$) = real_garbage$

case (rcollimator$)
  ele%offset_moves_aperture = .true.

case (rfcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_scale$) = 1
  ele%value(n_cell$) = 1

case (taylor$)   ! start with unit matrix
  ele%tracking_method = taylor$  
  ele%mat6_calc_method = taylor$ 
  call taylor_make_unit (ele%taylor)

case (wiggler$) 
  ele%sub_key = periodic_type$   
  ele%value(polarity$) = 1.0     

case (e_gun$)
  ele%tracking_method = time_runge_kutta$
  ele%mat6_calc_method = tracking$

end select

end subroutine parser_set_ele_defaults

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine settable_dep_var_bookkeeping (ele)
!
! Subroutine initialize dependent variables in an element.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine settable_dep_var_bookkeeping (ele)

use random_mod

implicit none

type (ele_struct) ele

real(rp) angle, rr
integer n
logical kick_set, length_set, set_done

!

kick_set = (ele%value(hkick$) /= 0) .or. (ele%value(vkick$) /= 0)

select case (ele%key)

! Convert rbends to sbends and evaluate G if needed.
! Needed is the length and either: angle, G, or rho.

case (sbend$, rbend$) 

  ele%sub_key = ele%key  ! save input format.
  angle = ele%value(angle$) 

  ! Only one of b_field, g, or rho may be set.

  if (ele%value(b_field$) /= 0 .and. ele%key == rbend$) call parser_error &
          ("B_FIELD NOT SETTABLE FOR AN RBEND (USE AN SBEND INSTEAD): " // ele%name)

  if (ele%value(b_field$) /= 0 .and. ele%value(g$) /= 0) call parser_error &
          ('BOTH G AND B_FIELD SET FOR A BEND: ' // ele%name)

  if (ele%value(b_field$) /= 0 .and. ele%value(rho$) /= 0) call parser_error &
          ('BOTH RHO AND B_FIELD SET FOR A BEND: ' // ele%name)

  if (ele%value(g$) /= 0 .and. ele%value(rho$) /= 0) &
            call parser_error ('BOTH G AND RHO SPECIFIED FOR BEND: ' // ele%name)

  ! if rho is set then this gives g

  if (ele%value(rho$) /= 0) ele%value(g$) = 1 / ele%value(rho$)

  ! If g and angle are set then this determines l

  length_set = .false.
  if (ele%value(g$) /= 0 .and. angle /= 0) then
    if (ele%value(l$) /= 0) call parser_error ('ANGLE, G/RHO, AND L SPECIFIED FOR BEND: ' // ele%name)
    ele%value(l$) = angle / ele%value(g$)
    length_set = .true.
  endif

  ! Convert an rbend

  if (ele%key == rbend$) then

    if (.not. length_set) then
      ele%value(l_chord$) = ele%value(l$)
      if (ele%value(l_chord$) == 0) then
        angle = 0
        ele%value(l$) = 0
      elseif (angle /= 0) then
        ele%value(l$) = ele%value(l_chord$) * angle / (2 * sin(angle/2))
      elseif (ele%value(g$) /= 0) then
        angle = 2 * asin(ele%value(l_chord$) * ele%value(g$) / 2)
        ele%value(l$) = ele%value(l_chord$) * angle / (2 * sin(angle/2))
      endif
    endif
    
    ele%value(e1$) = ele%value(e1$) + angle / 2
    ele%value(e2$) = ele%value(e2$) + angle / 2
    ele%key = sbend$

  endif

  ! 

  if (ele%value(angle$) /= 0 .and. ele%value(l$) == 0) then
    call parser_error ('THE BENDING ANGLE IS NONZERO IN A ZERO LENGTH BEND! ' // ele%name)
  elseif (ele%value(angle$) /= 0) then
    ele%value(g$) = ele%value(angle$) / ele%value(l$) 
  endif

  ! If fintx or hgapx are real_garbage then they have not been set.
  ! If so, set their valuse to fint and hgap.

  if (ele%value(hgapx$) == real_garbage$) ele%value(hgapx$) = ele%value(hgap$)
  if (ele%value(fintx$) == real_garbage$) ele%value(fintx$) = ele%value(fint$)

! Accept Use of Voltage for lcavities and vary the mode frequencies.

case (lcavity$) 

  if (ele%value(voltage$) /= 0) then
    if (ele%value(gradient$) /= 0) call parser_error &
                ('BOTH VOLTAGE AND GRADIENT NON-ZERO FOR A LCAVITY:', ele%name)
    ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)
  endif

! for a periodic_type wiggler n_pole is a dependent attribute

case (wiggler$)
  if (ele%sub_key == periodic_type$) then

    if (ele%value(l_pole$) == 0 .and. ele%value(n_pole$) /= 0) then
      ele%value(l_pole$) = ele%value(l$) / ele%value(n_pole$) 
    endif

  endif

! check for inconsistancies

case (solenoid$)
  if (ele%field_master .and. (ele%value(ks$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (KS, HKICK, ETC.) AND FIELD SET FOR A SOLENOID.')

case (sol_quad$)
  if (ele%field_master .and. (ele%value(ks$) /= 0 .or. &
                            ele%value(k1$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K1, HKICK, ETC.) AND FIELD SET FOR A SOL_QUAD.')

case (quadrupole$)
  if (ele%field_master .and. (ele%value(k1$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K1, HKICK, ETC.) AND FIELD SET FOR A QUAD.')

case (sextupole$)
  if (ele%field_master .and. (ele%value(k2$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K2, HKICK, ETC.) AND FIELD SET FOR A SEXTUPOLE.')

case (octupole$)
  if (ele%field_master .and. (ele%value(k3$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K3, HKICK, ETC.) AND FIELD SET FOR A OCTUPOLE.')

case (hkicker$, vkicker$)
  if (ele%field_master .and. (ele%value(kick$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH AND BL_KICK SET FOR A H/VKICKER.')

case (e_gun$)
  if (ele%value(gradient$) == 0 .and. ele%value(l$) /= 0) ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)

end select

! set ds_step if not already set.

if (attribute_index(ele, 'DS_STEP') > 0) then  ! If this is an attribute for this element...
  if (ele%value(num_steps$) > 0) then
    ele%value(ds_step$) = abs(ele%value(l$) / ele%value(num_steps$))
  elseif (ele%value(ds_step$) == 0) then
    if (ele%key == wiggler$ .and. ele%value(l_pole$) /= 0) then
      ele%value(ds_step$) = ele%value(l_pole$) / 10
    else
      ele%value(ds_step$) = bmad_com%default_ds_step
    endif
  endif
endif

! multipass lord needs to have field_master = T or must define a reference energy.

if (ele%lord_status == multipass_lord$ .and. .not. ele%field_master) then
  select case (ele%key)
  case (quadrupole$, sextupole$, octupole$, solenoid$, sol_quad$, sbend$, &
        hkicker$, vkicker$, kicker$, elseparator$, bend_sol_quad$)
    n = 0
    if (nint(ele%value(n_ref_pass$)) /= 0) n = n + 1 
    if (ele%value(p0c$) /= 0) n = n + 1
    if (ele%value(e_tot$) /= 0) n = n + 1
    if (n == 0) call parser_error ( &
          'FOR MULTIPASS LORD: ' // ele%name, &
          'N_REF_PASS, E_TOT, AND P0C ARE ALL ZERO AND FIELD_MASTER = FALSE!')
    if (n > 1) call parser_error ( &
          'FOR MULTIPASS LORD: ' // ele%name, &
          'MORE THAN ONE OF N_REF_PASS, E_TOT, AND P0C ARE SET NON-ZERO!')
  end select
endif


end subroutine settable_dep_var_bookkeeping 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine form_digested_bmad_file_name (lat_file, digested_file, full_lat_file)
!
! Subroutine to form the standard name of the Bmad digested file. 
! The standard digested file name has 'digested_' prepended to the file name.
!
! Modules needed:
!   use bmad_parser_mod
!
! Input:
!   lat_file -- Character(200): Input lattice file name.
!
! Output:
!   digested_file -- Character(200): Name of the digested file.
!   full_lat_file  -- Character(200), optional: Input lattice file name with full directory
!
!-

subroutine form_digested_bmad_file_name (lat_file, digested_file, full_lat_file)

character(*) lat_file, digested_file
character(*), optional :: full_lat_file
character(200) full_name, path, basename

integer ix

!

call fullfilename (lat_file, full_name)
inquire (file = full_name, name = full_name)  ! full input file_name
ix = index(full_name, ';')
if (ix /= 0) full_name = full_name(:ix-1)

if (present (full_lat_file)) full_lat_file = full_name

ix = SplitFileName(full_name, path, basename)
digested_file = trim(path) // 'digested_' // basename

! This only affects VMS programs.
! What we want to do is change the directory for lattice files in
! $CESR_ONLINE/machine_data/lattice...
! However 'CESR_ONLINE' is a logical that will get translated by the inquire function so
! we only check that 'lattice' is the top directory.

ix = max (index_nocase(digested_file, '[lattice.'), &
          index_nocase(digested_file, '[000000.lattice.'))
if (ix /= 0) then
  ix = index_nocase(digested_file, 'lattice.')
  digested_file = 'CESR_ONLINE:[machine_data.lattice.' // digested_file(ix+8:)
endif

end subroutine form_digested_bmad_file_name

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine reuse_taylor_elements (old_lat, lat)
!
! Subroutine to reuse saved taylor maps in a lattice.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine reuse_taylor_elements (old_lat, lat)

implicit none

type (lat_struct), target :: old_lat, lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: old_ele(:)

integer i, j, k, m, ib, n_old
integer, allocatable :: order(:)

! old_lat may be uninitalized (For example if the digested file did not exist.)
! In this case there is nothing to be done

if (.not. allocated (old_lat%branch)) return

! find the taylor order for each taylor map

if (.not. allocated(old_lat%branch)) return

n_old = 0
do i = 0, ubound(old_lat%branch, 1)
  branch => old_lat%branch(i)
  do j = 1, branch%n_ele_max
    if (associated(branch%ele(j)%taylor(1)%term)) n_old = n_old + 1
  enddo
enddo

allocate (old_ele(n_old), order(n_old))

n_old = 0
do i = 0, ubound(old_lat%branch, 1)
  branch => old_lat%branch(i)
  do j = 1, branch%n_ele_max
    ele => branch%ele(j)
    if (.not. associated(ele%taylor(1)%term)) cycle
    n_old = n_old + 1
    old_ele(n_old)%ele => ele
    order(n_old) = 0
    do k = 1, 6
      do m = 1, size(ele%taylor(k)%term)
        order(n_old) = max(order(n_old), maxval(ele%taylor(k)%term(m)%expn(:)))
      enddo
    enddo
  enddo
enddo

! Reuse the old taylor series for an element if: 
!   1) The old map exists and 
!   2) The new element has need of the map and
!   3) The old element and new element have the same attributes.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do i = 1, branch%n_ele_max

    ele => branch%ele(i)
    if (ele%tracking_method /= taylor$ .and. ele%tracking_method /= symp_map$ .and. &
        ele%mat6_calc_method /= taylor$ .and. ele%mat6_calc_method /= symp_map$) cycle

    call attribute_bookkeeper (ele, branch%param) ! for equivalent_taylor_attributes test

    do j = 1, n_old
      if (bmad_com%taylor_order > order(j)) cycle
      if (.not. equivalent_taylor_attributes (old_ele(j)%ele, ele)) cycle
      call transfer_ele_taylor (old_ele(j)%ele, ele, bmad_com%taylor_order)
      exit
    enddo

  enddo

enddo

deallocate(old_ele)

end subroutine reuse_taylor_elements

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_add_branch (branch_ele, lat, sequence, in_name, in_indexx, &
!                                                        seq_name, seq_indexx, in_lat)
!
! Subroutine to do line expansion.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

recursive subroutine parser_add_branch (branch_ele, lat, sequence, in_name, in_indexx, &
                                                        seq_name, seq_indexx, in_lat)

implicit none

type (lat_struct), target :: lat, in_lat
type (ele_struct) branch_ele
type (ele_struct), pointer :: b_ele, ele2
type (seq_struct), target :: sequence(:)
type (branch_struct), pointer :: branch

integer, allocatable :: seq_indexx(:), in_indexx(:)
integer j, nb, n_ele_use, n

character(*), allocatable ::  in_name(:), seq_name(:)

!

nb = ubound(lat%branch, 1) + 1
call allocate_branch_array (lat, nb)
branch_ele%value(ix_branch_to$) = nb
branch => lat%branch(nb)

! If branch_ele%value(particle$) has not been set then use default.
! Default for branch$ is just the same particle as the current branch.

if (branch_ele%value(particle$) == real_garbage$) then 
  if (branch_ele%key == branch$) then
    branch_ele%value(particle$) = lat%branch(branch_ele%ix_branch)%param%particle
  else
    branch_ele%value(particle$) = photon$
  endif
endif

!

branch%param%particle = nint(branch_ele%value(particle$))
branch%param%lattice_type = nint(branch_ele%value(lattice_type$))
branch%ix_branch      = nb
branch%ix_from_branch = branch_ele%ix_branch
branch%ix_from_ele    = branch_ele%ix_ele
branch%name           = branch_ele%name
call parser_expand_line (nb, lat, branch_ele%component_name, sequence, in_name, &
                              in_indexx, seq_name, seq_indexx, in_lat, n_ele_use)
branch%n_ele_track = n_ele_use
branch%n_ele_max   = n_ele_use

do j = 1, n_ele_use
  ele2 => lat%branch(nb)%ele(j)
  ! Make sure we don't have an endless loop
  n = nb
  do
    if (n == 0) exit
    b_ele => pointer_to_ele(lat, lat%branch(n)%ix_from_ele, lat%branch(n)%ix_from_branch)
    if (b_ele%name == ele2%name) then
      call parser_error ('ENDLESS BRANCHING LOOP DETECTED. BRANCHING ON BRANCHING ELEMENT: ' // ele2%name)
      if (bmad_status%exit_on_error) call err_exit
    endif
    n = b_ele%ix_branch
  enddo
  ! Now add the branch
  if (ele2%key == photon_branch$ .or. ele2%key == branch$) then
    call parser_add_branch (ele2, lat, sequence, in_name, in_indexx, &
                                                    seq_name, seq_indexx, in_lat)
  endif
enddo

end subroutine parser_add_branch

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_expand_line (ix_branch, lat, use_name, sequence, in_name, &
!                               in_indexx, seq_name, seq_indexx, in_lat, n_ele_use)
!
! Subroutine to do line expansion.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_expand_line (ix_branch, lat, use_name, sequence, in_name, &
                               in_indexx, seq_name, seq_indexx, in_lat, n_ele_use)

implicit none

type (lat_struct), target :: lat, in_lat
type (ele_struct), pointer :: ele_line(:)
type (seq_struct), target :: sequence(:)
type (seq_ele_struct), pointer :: s_ele, this_seq_ele
type (seq_stack_struct) stack(40)
type (seq_struct), pointer :: seq, seq2
type (used_seq_struct), allocatable ::  used2(:)
type (seq_ele_struct), target :: dummy_seq_ele

integer, allocatable :: ix_lat(:)
integer, allocatable :: seq_indexx(:), in_indexx(:)
integer iseq_tot, i_lev, i_use, n0_multi, n_ele_use, n_max
integer i, j, k, n, ix, ix_multipass, ix_branch

character(*), allocatable ::  in_name(:), seq_name(:)
character(*) use_name
character(40) name, multipass_line

! find line corresponding to the "use" statement.

iseq_tot = size(seq_indexx)
n_max = in_lat%n_ele_max

call find_indexx (use_name, seq_name, seq_indexx, iseq_tot, i_use)
if (i_use == 0) then
  call parser_error ('CANNOT FIND DEFINITION OF LINE IN "USE" STATEMENT: ' // use_name, ' ')
  return
endif

if (sequence(i_use)%type /= line$) then
  call parser_error ('NAME IN "USE" STATEMENT IS NOT A LINE!', ' ')
  return
endif

! Now to expand the lines and lists to find the elements to use.
! First go through the lines and lists and index everything.

do k = 1, iseq_tot
  do i = 1, size(sequence(k)%ele(:))

    s_ele => sequence(k)%ele(i)
    name = s_ele%name

    if (s_ele%ix_arg > 0) then   ! dummy arg
      s_ele%type = element$
      cycle
    endif

    call find_indexx2 (name, in_name, in_indexx, 0, n_max, j)
    if (j < 0) then  ! if not an element it must be a sequence
      call find_indexx (name, seq_name, seq_indexx, iseq_tot, j)
      if (j == 0) then  ! if not a sequence then I don't know what it is
        s_ele%ix_ele = -1
        s_ele%type = element$
      else
        s_ele%ix_ele = j
        s_ele%type = sequence(j)%type
      endif
      if (s_ele%type == list$ .and. s_ele%reflect) call parser_error ( &
                          'A REFLECTION WITH A LIST IS NOT ALLOWED IN: '  &
                          // sequence(k)%name, 'FOR LIST: ' // s_ele%name, &
                          seq = sequence(k))
      if (sequence(k)%type == list$) &
                call parser_error ('A REPLACEMENT LIST: ' // sequence(k)%name, &
                'HAS A NON-ELEMENT MEMBER: ' // s_ele%name)
 
    else    ! if an element...
      s_ele%ix_ele = j
      s_ele%type = element$
    endif

  enddo
enddo

! to expand the "used" line we use a stack for nested sublines.
! IX_LAT is the expanded array of elements in the lat.
! init stack

i_lev = 1                          ! level on the stack
seq => sequence(i_use)

stack(1)%ix_seq    = i_use           ! which sequence to use for the lat
stack(1)%ix_ele    =  1              ! we start at the beginning
stack(1)%direction = +1              ! and move forward
stack(1)%rep_count = seq%ele(1)%rep_count
stack(1)%multipass = seq%multipass
stack(1)%tag = ''

n_ele_use = 0
         
allocate (ix_lat(n_max))
ix_lat = -1
sequence(:)%ix = 1  ! Init. Used for replacement list index

if (stack(1)%multipass) then
  call parser_error ('"USE"D LINE FOR LATTICE EXPANSION IS MARKED MULTIPASS!')
  if (bmad_status%exit_on_error) call err_exit
endif

!-------------------------------------------------------------------------
! Expand "used" line...

line_expansion: do

  ! If rep_count is zero then we are finished with the current element.

  if (stack(i_lev)%rep_count == 0) then      ! goto next element in the sequence
    ! goto the next element by changing %ix_ele index by +/- 1 
    stack(i_lev)%ix_ele = stack(i_lev)%ix_ele + stack(i_lev)%direction 
    ix = stack(i_lev)%ix_ele

    ! Check if off the end of the current line...
    if (ix > 0 .and. ix <= size(seq%ele)) then  ! Nope. Still have more element to process
      stack(i_lev)%rep_count = seq%ele(ix)%rep_count  ! set the rep_count for the next ele.
      if (stack(i_lev)%rep_count == 0) cycle          ! For "0*sub_line" construct.

    ! If we have got to the end of the current line then pop the stack back to
    ! the next lower level.
    else  
      i_lev = i_lev - 1
      if (i_lev == 0) exit line_expansion    ! level 0 -> we are done.
      seq => sequence(stack(i_lev)%ix_seq)

      ! If we have stepped out of a multipass line which has been transversed in reverse
      !   then we need to do some bookkeeping to keep the elements straight.
      if (.not. stack(i_lev)%multipass .and. stack(i_lev+1)%multipass) then
        if (stack(i_lev+1)%direction == -1) then
          bp_com%used_line(n0_multi:n_ele_use)%ix_multipass = &
                        bp_com%used_line(n_ele_use:n0_multi:-1)%ix_multipass
        endif
      endif
      cycle
    endif

  endif

  ! 

  s_ele => seq%ele(stack(i_lev)%ix_ele)  ! next element, line, or list
  stack(i_lev)%rep_count = stack(i_lev)%rep_count - 1

  ! if s_ele is a dummy arg then get corresponding actual arg.

  ix = s_ele%ix_arg
  if (ix /= 0) then  ! it is a dummy argument.
    name = seq%corresponding_actual_arg(ix)
    s_ele => dummy_seq_ele
    s_ele%name = name
    call find_indexx2 (name, in_name, in_indexx, 0, n_max, j)
    if (j < 0) then  ! if not an element it must be a sequence
      call find_indexx (name, seq_name, seq_indexx, iseq_tot, j)
      if (j == 0) then  ! if not a sequence then I don't know what it is
        call parser_error ('CANNOT FIND DEFINITION FOR: ' // name, &
                          'IN LINE: ' // seq%name, seq = seq)
        if (bmad_status%exit_on_error) call err_exit
        return
      endif
      s_ele%ix_ele = j
      s_ele%type = sequence(j)%type
    else
      s_ele%ix_ele = j 
      s_ele%type = element$
    endif
    
  endif

  ! if an element

  select case (s_ele%type)

  case (element$, list$) 

    if (s_ele%type == list$) then
      seq2 => sequence(s_ele%ix_ele)
      j = seq2%ix
      this_seq_ele => seq2%ele(j)
      seq2%ix = seq2%ix + 1
      if (seq2%ix > size(seq2%ele(:))) seq2%ix = 1
    else
      if (s_ele%tag /= '') then
        call parser_error ('ELEMENTS IN A LINE OR LIST ARE NOT ALLOWED TO HAVE A TAG.', &
                      'FOUND ILLEGAL TAG FOR ELEMENT: ' // s_ele%name, &
                      'IN THE LINE/LIST: ' // seq%name, seq)
      endif
      this_seq_ele => s_ele
    endif

    if (this_seq_ele%ix_ele < 1) call parser_error('NOT A DEFINED ELEMENT: ' // &
                          s_ele%name, 'IN THE LINE/LIST: ' // seq%name, seq = seq)


    if (n_ele_use+1 > size(ix_lat)) then
      n = 1.5*n_ele_use
      call re_allocate (ix_lat, n)
      ix = size(bp_com%used_line) 
      allocate (used2(ix))
      used2(1:ix) = bp_com%used_line(1:ix)
      deallocate (bp_com%used_line)
      allocate (bp_com%used_line(1:n))
      bp_com%used_line(1:ix) = used2(1:ix)
      deallocate (used2)
    endif

    call pushit (ix_lat, n_ele_use, this_seq_ele%ix_ele)

    bp_com%used_line(n_ele_use)%name = this_seq_ele%name

    if (stack(i_lev)%tag /= '' .and. s_ele%tag /= '') then
      bp_com%used_line(n_ele_use)%tag =  trim(stack(i_lev)%tag) // '.' // s_ele%tag
    elseif (s_ele%tag /= '') then
      bp_com%used_line(n_ele_use)%tag = s_ele%tag
    else
      bp_com%used_line(n_ele_use)%tag =  stack(i_lev)%tag
    endif

    if (stack(i_lev)%multipass) then
      ix_multipass = ix_multipass + 1
      bp_com%used_line(n_ele_use)%ix_multipass = ix_multipass
      bp_com%used_line(n_ele_use)%multipass_line = multipass_line
    else
      bp_com%used_line(n_ele_use)%ix_multipass = 0
    endif


  ! if a line:
  !     a) move pointer on current level past line element
  !     b) go to the next higher level
  !     c) initialize pointers for the higher level to use the line

  case (line$, replacement_line$)
    i_lev = i_lev + 1
    if (i_lev > size(stack)) then
      call parser_error ('NESTED LINES EXCEED STACK DEPTH! SUSPECT INFINITE LOOP!')
      if (bmad_status%exit_on_error) call err_exit
    endif
    if (s_ele%type == replacement_line$) then
      seq2 => sequence(s_ele%ix_ele)
      if (size(seq2%dummy_arg) /= size(s_ele%actual_arg)) then
        call parser_error ('WRONG NUMBER OF ARGUMENTS FORREPLACEMENT LINE: ' // &
            s_ele%name, 'WHEN USED IN LINE: ' // seq%name, seq = seq)
        return
      endif
      arg_loop: do i = 1, size(seq2%dummy_arg)
        seq2%corresponding_actual_arg(i) = s_ele%actual_arg(i)
        if (associated(seq%dummy_arg)) then
          do j = 1, size(seq%dummy_arg)
            if (seq2%corresponding_actual_arg(i) == seq%dummy_arg(j)) then
              seq2%corresponding_actual_arg(i) = seq%corresponding_actual_arg(j)
              cycle arg_loop
            endif
          enddo
        endif
        name = seq2%corresponding_actual_arg(i)
      enddo arg_loop
    endif

    seq => sequence(s_ele%ix_ele)
    stack(i_lev)%ix_seq = s_ele%ix_ele
    stack(i_lev)%direction = stack(i_lev-1)%direction
    stack(i_lev)%multipass = (stack(i_lev-1)%multipass .or. seq%multipass)
    if (stack(i_lev-1)%tag /= '' .and. s_ele%tag /= '') then
       stack(i_lev)%tag = trim(stack(i_lev-1)%tag) // '.' // s_ele%tag
    elseif (stack(i_lev-1)%tag /= '') then
       stack(i_lev)%tag = trim(stack(i_lev-1)%tag)
    else
       stack(i_lev)%tag = s_ele%tag
    endif
    if (s_ele%reflect) stack(i_lev)%direction = -stack(i_lev)%direction

    if (stack(i_lev)%direction == 1) then
      ix = 1
    else
      ix = size(seq%ele(:))
    endif

    stack(i_lev)%ix_ele = ix
    stack(i_lev)%rep_count = seq%ele(ix)%rep_count

    if (stack(i_lev)%multipass .and. .not. stack(i_lev-1)%multipass) then
      ix_multipass = 1
      n0_multi = n_ele_use + 1
      multipass_line = sequence(stack(i_lev)%ix_seq)%name
    endif

  case default
    call parser_error ('INTERNAL SEQUENCE ERROR!')

  end select

enddo line_expansion

! Stop here if there has been an error

if (bp_com%error_flag) return

! Transfer the ele information from the in_lat to lat and
! do the bookkeeping for settable dependent variables.

if (ix_branch == 0) then  ! Main branch
  call allocate_lat_ele_array(lat, n_ele_use)
  ele_line => lat%ele
else                    ! branch line
  call allocate_lat_ele_array(lat, n_ele_use, ix_branch)
  ele_line => lat%branch(ix_branch)%ele
endif

ele_line(0)%ix_branch = ix_branch

do i = 1, n_ele_use
  ele_line(i) = in_lat%ele(ix_lat(i)) 
  ele_line(i)%name = bp_com%used_line(i)%name
  if (bp_com%used_line(i)%tag /= '') ele_line(i)%name = &
                trim(bp_com%used_line(i)%tag) // '.' // ele_line(i)%name
  call settable_dep_var_bookkeeping (ele_line(i))
enddo

! Cleanup

deallocate (ix_lat)

end subroutine parser_expand_line

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine bp_set_ran_status
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine bp_set_ran_status

if (bp_com%ran%deterministic == 0) then
  bp_com%ran%ran_function_was_called = .true.
elseif (bp_com%ran%deterministic == 1) then
  bp_com%ran%deterministic_ran_function_was_called = .true.
endif

end subroutine bp_set_ran_status
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_debug_print_info (debug_line)
!
! Subroutine to remove all null_ele elements.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_debug_print_info (lat, debug_line)

type (lat_struct) lat
character(*) debug_line
integer i, ix

!

call str_upcase (debug_line, debug_line)

if (index(debug_line, 'VAR') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'Number of Defined Variables:', bp_com%ivar_tot - bp_com%ivar_init
  do i = bp_com%ivar_init+1, bp_com%ivar_tot
    print *
    print *, 'Var #', i
    print *, 'Name: ', bp_com%var(i)%name
    print *, 'Value:', bp_com%var(i)%value
  enddo
endif

if (index(debug_line, 'SLAVE') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'Number of Elements in Tracking Lattice:', lat%n_ele_track
  do i = 1, lat%n_ele_track
    print *, '-------------'
    print *, 'Ele #', i
    call type_ele (lat%ele(i), .false., 0, .false., 0, .true., lat, .true., .false., .true., .true.)
  enddo
endif

if (index(debug_line, 'LORD') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'LORD elements: ', lat%n_ele_max - lat%n_ele_track
  do i = lat%n_ele_track+1, lat%n_ele_max
    print *, '-------------'
    print *, 'Ele #', i
    call type_ele (lat%ele(i), .false., 0, .false., 0, .true., lat, .true., .false., .true., .true.)
  enddo
endif

if (index(debug_line, 'LATTICE') /= 0) then  
  print *
  print *, '----------------------------------------'
  print *, 'Lattice Used: ', lat%use_name
  print *, 'Number of lattice elements:', lat%n_ele_track
  print *, 'List:                                 Key                 Length         S'
  do i = 1, lat%n_ele_track
    print '(i4, 2a, 3x, a, 2f10.2)', i, ') ', lat%ele(i)%name(1:30),  &
      key_name(lat%ele(i)%key), lat%ele(i)%value(l$), lat%ele(i)%s
  enddo
  print *, '---- Lord Elements ----'
  do i = lat%n_ele_track+1, lat%n_ele_max
    print '(2x, i4, 2a, 3x, a, 2f10.2)', i, ') ', lat%ele(i)%name(1:30),  &
           key_name(lat%ele(i)%key), lat%ele(i)%value(l$), lat%ele(i)%s
  enddo
endif

ix = index(debug_line, 'ELE')
if (ix /= 0) then
  print *
  print *, '----------------------------------------'
  call string_trim (debug_line(ix+3:), debug_line, ix)
  do
    if (ix == 0) exit
    read (debug_line, *) i
    print *
    print *, '----------------------------------------'
    print *, 'Element #', i
    call type_ele (lat%ele(i), .false., 0, .true., 0, .true., lat, .true., .true., .true., .true.)
    call string_trim (debug_line(ix+1:), debug_line, ix)
  enddo
endif

if (index(debug_line, 'BEAM_START') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'beam_start:'
  print '(3x, 6es13.4)', lat%beam_start%vec      
endif

end subroutine parser_debug_print_info

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! parse_map (grid, ele, lat, delim, delim_found, err_flag, print_err)
!
! Subroutine to parse a "map = {}" construct
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is private to bmad_parser_mod.
! This must read in:
! {type = ,
!    dr = , 
!    r0 = , 
!    pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) ) 
!    .
!    .
!    . ) },
!-

subroutine parse_map (map, ele, lat, delim, delim_found, err_flag, print_err)

implicit none

type (em_field_map_struct), pointer :: map
type (ele_struct), target :: ele
type (lat_struct) lat

real(rp), allocatable :: array(:)

complex(rp), pointer :: c_ptr(:)

integer ix_word, i_term

character(1) delim, delim2
character(40) word, word2, name

logical err_flag, print_err, delim_found

!

err_flag = .true.
if (.not. associated(map)) allocate (map)
allocate (array(1024))

! Expect {

call get_next_word (word, ix_word, '{', delim, delim_found, call_check = .true.)
if ((word /= '') .or. (delim /= '{')) then
  call parser_error ('NO { SIGN FOUND IN MAP DEFINITION',  &
                    'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
  return
endif

!

do

  ! Read attriubute
  call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

  select case (word)

  case ('ELE_ANCHOR_PT')
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER MAP ' // word,  &
                         'IN MAP STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif
    call get_next_word (word2, ix_word, ',}', delim, delim_found)

    ! Evaluate string into integer.

    if (word == 'ELE_ANCHOR_PT') then
      call match_word(word2, anchor_pt_name, map%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
    endif
  
    if (name == '') then
      call parser_error ('UNKNKOWN MAP ' // trim(word) // ': ' // word2, &
                         'FOUND IN MODE MAP DEFINITION FOR ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('E_COEF_RE', 'E_COEF_IM', 'B_COEF_RE', 'B_COEF_IM')

    ! Expect "("
    call get_next_word (word, ix_word, ',({', delim, delim_found)
    if (word /= '' .or. delim /= '(') then
      call parser_error ('NO "(" FOUND AFTER "' // trim(word2) // ' =" ', &
                           'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif

    ! Read list of values.
    call re_allocate(array, 1024, .false.)
    do i_term = 1, 100000
      call get_next_word (word, ix_word, '{},()', delim, delim_found)
      if ((delim /= ',' .and. delim /= ')') .or. .not. is_real(word)) then
        call parser_error ('ERROR PARSING ARRAY FOR: ' // word2, &
                             'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      endif
      if (i_term > size(array)) call re_allocate(array, 2*size(array))
      read (word, *) array(i_term)
      if (delim == ')') exit
    enddo

    if (allocated(map%term)) then
      if (size(map%term) /= i_term) then
        call parser_error ('ARRAY SIZE MISMATCH FOR: ' // word2, &
                           'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      endif
    else
      allocate(map%term(i_term))
    endif

    select case (word2)
    case ('E_COEF_RE', 'E_COEF_IM'); c_ptr => map%term%e_coef 
    case ('B_COEF_RE', 'B_COEF_IM'); c_ptr => map%term%b_coef
    end select

    if (word2(8:9) == 'RE') then
      if (any(real(c_ptr) /= 0)) then
        call parser_error ('DUPLICATE ARRAY FOR: ' // word2, &
                           'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      endif
      c_ptr = c_ptr + array(1:i_term)

    else
      if (any(aimag(c_ptr) /= 0)) then
        call parser_error ('DUPLICATE ARRAY FOR: ' // word2, &
                           'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      endif
      c_ptr = c_ptr + i_imaginary * array(1:i_term)
    endif

    ! Expect "," or "}"
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    if (word /= '' .or. (delim /= ',' .and. delim /= '}')) then
      call parser_error ('BAD ' // trim(word2) // ' = (...) CONSTRUCT', &
                           'FOUND IN MODE DEFINITION IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif

  case ('DZ');            
    if (.not. associated(map)) allocate (map)
    call evaluate_value (trim(ele%name), map%dz, lat, delim, delim_found, err_flag, ',}')

  end select

  ! Possible "}" is end of mode
  if (delim == '}') exit

  call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

enddo

! Get final separator after grid construct.
 
call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

if (ix_word /= 0) then
    call parser_error ('UNKNOWN TEXT FOUND AFTER MODE GRID DEFINITION: ' // word, &
                         'FOR ELEMENT: ' // ele%name)
    return 
endif

! Deallocate

deallocate(array)
err_flag = .false.

end subroutine parse_map

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! parse_grid (grid, ele, lat, delim, delim_found, err_flag, print_err)
!
! Subroutine to parse a "grid = {}" construct
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is private to bmad_parser_mod.
! This must read in:
! {type = ,
!    dr = , 
!    r0 = , 
!    pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) ) 
!    .
!    .
!    . ) },
!-

subroutine parse_grid(grid, ele, lat, delim, delim_found, err_flag, print_err)


type grid_pt_struct
  integer :: ix1 , ix2, ix3
  complex(rp) :: field(6) = 0
end type


type (em_field_grid_struct), pointer :: grid
type (ele_struct) :: ele
type (ele_struct), pointer :: bele
type (lat_struct),  target :: lat
type (branch_struct), pointer :: branch
type(grid_pt_struct), allocatable :: array(:), array2(:)

character(1) delim, delim2
character(40) :: word, word2, name

real(rp), allocatable :: dr(:), r0(:)

integer ix_word, ix_word2
integer pt_counter, n, i, ib, ie, im, ix1, ix2, ix3, max_ix1, max_ix2, max_ix3
integer grid_dim,  num_dr, num_r0

logical delim_found, delim_found2, err_flag, print_err, err_flag2

!

if (.not. associated(grid)) allocate (grid)

! Expect {

call get_next_word (word, ix_word, '{', delim, delim_found, call_check = .true.)
if ((word /= '') .or. (delim /= '{')) then
  call parser_error ('NO { SIGN FOUND IN GRID DEFINITION',  &
                    'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
  return
endif

! Set %file to be the last called file with 

write(grid%file, '(2a, i0)') trim(bp_com%current_file%full_name),  ':', bp_com%current_file%i_line

! Init

allocate(array(1024))
pt_counter = 0
err_flag = .true.

do    

  ! Read attriubute
  call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

  select case (word)

  case ('TYPE', 'ELE_ANCHOR_PT')
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER GRID ' // word,  &
                         'IN GRID STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif
    call get_next_word (word2, ix_word, ',}', delim, delim_found)
    ! Check to see if this is a valid type by checking against em_grid_type_name(:)

    if (word == 'TYPE') then
      call match_word(word2, em_grid_type_name, grid%type, can_abbreviate = .false., matched_name = name)
    else
      call match_word(word2, anchor_pt_name, grid%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
    endif
  
    if (name == '') then
      call parser_error ('UNKNKOWN GRID ' // trim(word) // ': ' // word2, &
                         'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('R0')       
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER GRID R0',  &
                      'IN GRID STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif
    ! expect ( 1. ) or (1. , 2.) or (1., 2., 3.)
    call parse_real_list(err_flag2, r0, num_r0, num_expected = 3)
    if (err_flag2) return
    grid%r0 = r0
    !Expect , or }
    call get_next_word (word, ix_word, ',}', delim, delim_found)     
    if (word /= '') then
      call parser_error ('BAD INPUT AFTER R0 DEFINITION: ' // word , &
                                 'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if

    case ('DR')      
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER GRID DR',  &
                      'IN GRID STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif      
    ! expect ( 1.) or (1. , 2.) or (1., 2., 3.)
    call parse_real_list(err_flag2, dr, num_dr, num_expected = 3)
    if (err_flag2) return
    grid%dr = dr
    call get_next_word (word, ix_word, ',}', delim, delim_found)     
    if (word /= '') then
      call parser_error ('BAD INPUT AFTER DR DEFINITION: ' // word , &
                                 'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if

  case ('PT')
    !Increment 
    pt_counter = pt_counter +1
      !Reallocate temporary structure if needed
    if (pt_counter > size(array)) then
      n = size(array)
      allocate( array2(n))
      array2 = array
      deallocate(array)
      allocate(array(2*n))
      array(1:n) = array2
      deallocate(array2)
    end if

    !Expect "("
    if (delim /= '(') then
      call parser_error ('BAD GRID PT CONSTRUCT', &
                           'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if
     
    !Get indices
    call parse_pt_indices(array(pt_counter)%ix1, delim, err_flag2)
    if (err_flag2) return
    if ( delim == ',') call parse_pt_indices(array(pt_counter)%ix2, delim, err_flag2)
    if (delim == ',') call parse_pt_indices(array(pt_counter)%ix3, delim, err_flag2) 

    !Expect last delim was ")"
    if (delim /= ')') then
      call parser_error ('BAD GRID PT CONSTRUCT, NO CLOSING ")" FOR INDICES', &
      'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if
      
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    call get_next_word (word2, ix_word2, '{}=,()', delim2, delim_found2)
    if ((word /= '') .or. (word2 /= '') &
         .or. (delim /= '=') .or. (delim2 /= '(')) then
      call parser_error ('BAD GRID PT CONSTRUCT, NO  = "(" ', &
                 'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if
    !Get as many field components as listed
    do i = 1, 6
      call parse_complex_component(array(pt_counter)%field(i), delim, err_flag2)
      if (err_flag2) return
      if (delim == ')') exit
      if (delim /= ',') then
        call parser_error ('BAD GRID PT CONSTRUCT, NO "," BETWEEN FIELD COMPONENTS', &
            'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
        return
      end if    
    end do
    !Expect , or }
    call get_next_word (word, ix_word, ',}', delim, delim_found) 
    if (word /= '') then
      call parser_error ('BAD INPUT AFTER PT DEFINITION: ' // word , &
                           'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if
  
  case default
    call parser_error ('UNKNOWN GRID ATTRIBUTE: ' // word, &
                         'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
    return 
        
  end select 

  if (delim == '}') exit   

enddo

! Get final separator after grid construct.
 
call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

if (ix_word /= 0) then
    call parser_error ('UNKNOWN TEXT FOUND AFTER MODE GRID DEFINITION: ' // word, &
                         'FOR ELEMENT: ' // ele%name)
    return 
endif

!Clear pts

if (allocated(grid%pt)) deallocate(grid%pt)

!Allocate grid for different dimensions

grid_dim = em_grid_dimension(grid%type)
select case(grid_dim)
  case (1)
    max_ix1 = maxval(array(1:pt_counter)%ix1)
    ix2 = 1
    ix3 = 1
    allocate(grid%pt(0:max_ix1, 1:1, 1:1))
  case (2)
    max_ix1 = maxval(array(1:pt_counter)%ix1)
    max_ix2 = maxval(array(1:pt_counter)%ix2)
    ix3 = 1
    allocate(grid%pt(0:max_ix1, 0:max_ix2, 1:1))
  case (3)
    max_ix1 = maxval(array(1:pt_counter)%ix1)
    max_ix2 = maxval(array(1:pt_counter)%ix2)
    max_ix3 = maxval(array(1:pt_counter)%ix3)
    allocate(grid%pt(0:max_ix1, 0:max_ix2, 0:max_ix3))
  case default
    call parser_error ('BAD GRID DIMENSION', &
               'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
    return
end select

!Assign grid values
do i = 1, pt_counter
  ix1 = array(i)%ix1
  if (grid_dim >1)   ix2 = array(i)%ix2
  if (grid_dim == 3) ix3 = array(i)%ix3
  grid%pt(ix1, ix2, ix3)%E(1:3) = array(i)%field(1:3)
  grid%pt(ix1, ix2, ix3)%B(1:3) = array(i)%field(4:6)
end do

! Clear temporary array

deallocate(array)

! Check if grid data has already been read in for another element.
! If so, save space by pointing to the existing grid.

branch_loop: do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 1, branch%n_ele_max
    bele => branch%ele(ie)
    if (.not. associated(bele%em_field)) cycle    
    if (bele%ix_ele == ele%ix_ele .and. bele%ix_branch == ele%ix_branch) cycle
    do im = 1, size(bele%em_field%mode)
      if (.not. associated(bele%em_field%mode(im)%grid)) cycle
      if (bele%em_field%mode(im)%grid%file /= grid%file) cycle
      deallocate(grid)
      grid => bele%em_field%mode(im)%grid
      grid%n_link = grid%n_link + 1        
      exit branch_loop
    end do
  enddo 
enddo branch_loop

err_flag = .false.

!------------------------------------------------------------

contains

subroutine parse_pt_indices (ix, delim, err_flag2)

character(1) delim
character(40) :: word
integer ix, ix_word
logical delim_found, err_flag2

!

err_flag2 = .true.

!Expect integer followed by "," or ")"

call get_next_word (word, ix_word, ',)', delim, delim_found)
if (.not. is_integer(word)) then
  call parser_error ('BAD GRID PT INDEX, NOT AN INTEGER', &
                         'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
  return
end if       
read (word, *) ix
err_flag2 = .false.

end subroutine parse_pt_indices

!-----------------------------------------------------------
! contains 
! subroutine parse_complex_component(complex_component, delim, err_flag2)
! looks for (x, y) or x followed by , or ) 
! returns complex field_component and next delim, which should be , or )
subroutine parse_complex_component(complex_component, delim, err_flag2)

character(1) delim, delim2
character(40) :: word, word2
integer  ix_word, ix_word2
logical delim_found, delim_found2, err_flag2
real(rp) x, y
complex(rp) complex_component

!

err_flag2 = .true.

! Expect "(" for complex, "," for real in the middle of the list, and ")" at the end of the list

call get_next_word (word, ix_word, ',()', delim, delim_found)
if (is_real(word)) then
  read (word, *) x
  complex_component = cmplx(x, 0.0_rp)
else if (delim == '(') then
  call get_next_word (word, ix_word, ',', delim, delim_found)
  call get_next_word (word2, ix_word2, ')', delim2, delim_found2)
  !Expect: real, real ) 
  if ((.not. is_real(word)) .or. (.not. is_real(word2)) &
    .or. (delim /= ',') .or. (delim2 /= ')') &
    .or. (.not. delim_found) .or. (.not. delim_found2)) then
    call parser_error ('BAD COMPLEX COMPONENT CONSTRUCT', &
                         'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
    return
   end if
   !
   read (word, *) x
      read (word2, *) y
      complex_component = cmplx(x,y)
      !Look for "," or end of list ")"
   call get_next_word (word, ix_word, ',)', delim, delim_found)
end if

err_flag2 = .false.

end subroutine parse_complex_component


end subroutine parse_grid


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! parse_integer_list (integer_array, num_found, 
!           num_expected = 1, open_delim = '(', 
!           separator = ',', close_delim = ')', 
!           do_resize = .false.)
!
! subroutine to parse a list of integers of the form
!    open_delim integer_1 separator integer_2 . . . close_delim
!
!       example:   (1, 2, 4, 8) 
!
! Note: nearly identical to subroutine parse_real_list
!
! Input:
!  integer_array -- Integer, allocatable: the array to be read in 
!
!   Optional: 
!   num_expected = 1     -- integer : number of expected arguments
!              Used to initialize integer_array
!   open_delim   = '('   -- character(1) : opening delimeter
!   separator    = ','   -- character(1) : separating character
!   close_delim  = ')'   -- character(1) : closing delimeter
!   do_resize    = .false.  -- logical : resize integer_array to num_found
!
! Output:
!   integer_array(1:num_found) --integer(rp) : Array of values
!   num_found    -- integer : number of elements


subroutine parse_integer_list(integer_array, num_found, &
          num_expected, open_delim, separator, close_delim, &
          do_resize)

!Arguments
integer, allocatable :: integer_array(:)
integer :: num_found
integer, optional :: num_expected
character(1), optional :: open_delim, close_delim, separator
logical, optional :: do_resize

!Local
integer num_expect
character(1) delim, op_delim, cl_delim, sep
character(40) :: word
integer  ix_word
logical delim_found, resize

character(20), parameter :: r_name =  'parse_integer_list'

!Optional arguments
num_expect = 1
op_delim = '('
cl_delim = ')'
sep      = ','
resize = .false. 
if (present(num_expected)) num_expect = num_expected
if (present(open_delim)) op_delim = open_delim
if (present(close_delim)) cl_delim = close_delim
if (present(separator)) sep = separator
if (present(do_resize)) resize = do_resize


!Expect op_delim
call get_next_word (word, ix_word, op_delim, delim, delim_found)
if ((word /= '') .or. (delim /= op_delim)) then
  call parser_error ('BAD OPENING DELIMITER', &
                         'IN: ' // r_name)
    return
end if

!Initial allocation
call re_allocate(integer_array, num_expected, .false.)

!counter
num_found = 0

!Get integers
do 
  call get_next_word (word, ix_word, sep // cl_delim, delim, delim_found)
  if (.not. is_integer(word) ) then 
    call parser_error ('BAD ELEMENT: ' // word, &
                         'IN : ' // r_name)
    return
   end if    
  !integer is found
  num_found = num_found + 1
  !reallocate if needed  
  if (size(integer_array) < num_found) call re_allocate (integer_array, 2*num_found, .false.)
  
  !Read value
   read (word, *)  integer_array(num_found) 
  
  !Exit if cl_delim is found
  if (delim == cl_delim) exit
  
  !Check separator
  if (delim /= sep) then
    call parser_error ('BAD SEPARATOR', &
                         'IN: ' // r_name)
    return  
  end if
  
end do

!Resize if asked
if (resize) call re_allocate(integer_array, num_found, .true.)

end subroutine parse_integer_list

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! parse_real_list (err_flag, real_array, num_found, 
!           num_expected = 1, open_delim = '(', 
!           separator = ',', close_delim = ')', 
!           do_resize = .false., default_value = 0.0_rp)
!
! subroutine to parse a list of reals of the form
!    open_delim real_1 separator real_2 . . . close_delim
!
!       example:   (1.2, 2.3, 4.4, 8.5) 
!
! Note: nearly identical to subroutine parse_integer_list
!
! Input:
!  real_array -- real(rp), allocatable: the array to be read in 
!
!   Optional: 
!   num_expected = 1     -- integer : number of expected arguments
!              Used to initialize integer_array
!   open_delim   = '('   -- character(1) : opening delimeter
!   separator    = ','   -- character(1) : separating character
!   close_delim  = ')'   -- character(1) : closing delimeter
!   do_resize    = .false.  -- logical : resize integer_array to num_found
!   default_value = 0.0_rp  -- real(rp) : inital assignment of real_array elements
!                  useful when do_resize = .false. 
!
! Output:
!   real_array(1:num_found) -- real(rp) : Array of values
!   num_found    -- integer : number of elements


subroutine parse_real_list (err_flag, real_array, num_found, &
          num_expected, open_delim, separator, close_delim, &
          do_resize, default_value)

! Arguments

real(rp), allocatable :: real_array(:)
integer :: num_found
integer, optional :: num_expected
character(1), optional :: open_delim, close_delim, separator
logical err_flag
logical, optional :: do_resize
real(rp), optional :: default_value

! Local
integer num_expect
character(1) delim, op_delim, cl_delim, sep
character(40) :: word
integer  ix_word
logical delim_found, resize
real(rp) :: default_val

character(20), parameter :: r_name =  'parse_real_list'

! Optional arguments

err_flag = .true.
num_expect = 1
op_delim = '('
cl_delim = ')'
sep      = ','
resize = .false. 
default_val = 0.0_rp
if (present(num_expected)) num_expect = num_expected
if (present(open_delim)) op_delim = open_delim
if (present(close_delim)) cl_delim = close_delim
if (present(separator)) sep = separator
if (present(do_resize)) resize = do_resize
if (present(default_value)) default_val = default_value


!Expect op_delim
call get_next_word (word, ix_word, op_delim, delim, delim_found)
if ((word /= '') .or. (delim /= op_delim)) then
  call parser_error ('BAD OPENING DELIMITER', 'IN: ' // r_name)
  return
end if

!Initial allocation
call re_allocate(real_array, num_expected, .false.)
real_array = default_val

!counter
num_found = 0

!Get integers
do 
  call get_next_word (word, ix_word, sep // cl_delim, delim, delim_found)
  if (.not. is_real(word) ) then 
    call parser_error ('BAD ELEMENT: ' // word, 'IN : ' // r_name)
    return
   end if    
  !real is found
  num_found = num_found + 1
  !reallocate if needed  
  if (size(real_array) < num_found) call re_allocate (real_array, 2*num_found, .false.)
  
  !Read value
   read (word, *)  real_array(num_found) 
  
  !Exit if cl_delim is found
  if (delim == cl_delim) exit
  
  !Check separator
  if (delim /= sep) then
    call parser_error ('BAD SEPARATOR', 'IN: ' // r_name)
    return  
  end if

end do

!Resize if asked
if (resize) call re_allocate(real_array, num_found, .true.)

err_flag = .false.

end subroutine parse_real_list

end module
