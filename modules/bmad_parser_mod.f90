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
use cross_section_mod

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
  integer i_line
  integer f_unit
end type

!-----------------------------------------------------------
! structure for holding the control names and pointers for
! superimpose and overlay elements

integer, private :: plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
integer, private :: l_parens$ = 5, r_parens$ = 6, power$ = 7
integer, private :: unary_minus$ = 8, unary_plus$ = 9, no_delim$ = 10
integer, private :: sin$ = 11, cos$ = 12, tan$ = 13
integer, private :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17, sqrt$ = 18
integer, private :: log$ = 19, exp$ = 20, ran$ = 21, ran_gauss$ = 22
integer, private :: numeric$ = 100

integer, private :: eval_level(22) = (/ 1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
                            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 /)


type eval_stack_struct
  integer type
  real(rp) value
end type

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
  logical common_lord
end type

type parser_lat_struct
  type (parser_ele_struct), pointer :: ele(:) => null()
end type

!

integer, parameter :: line$ = 1, list$ = 2, element$ = 3
integer, parameter :: replacement_line$ = 4

integer begin$, center$, end$
parameter (begin$ = -1)
parameter (center$ = 0)
parameter (end$ = 1)

integer def$, redef$
parameter (def$ = 1)
parameter (redef$ = 2)

!------------------------------------------------
! common stuff

integer, parameter :: n_parse_line = 280

type bp_var_struct
  character(40) name   ! variable name
  real(rp) value       ! variable value
  integer indexx       ! variable sort index
end type

type bp_common_struct
  type (ele_struct), pointer :: param_ele, beam_ele, beam_start_ele
  logical e_tot_set, p0c_set
  type (stack_file_struct), pointer :: current_file
  type (stack_file_struct), pointer :: calling_file
  type (lat_struct), pointer :: old_lat
  type (used_seq_struct), allocatable ::  used_line(:)
  type (bp_var_struct), allocatable :: var(:)   ! variable name
  integer num_lat_files               ! Number of files opened
  integer ivar_tot, ivar_init
  character(200), allocatable :: lat_file_names(:) ! List of all files used to create lat
  character(n_parse_line) parse_line
  character(n_parse_line) input_line1          ! For debug messages
  character(n_parse_line) input_line2          ! For debug messages
  character(40) parser_name
  character(200) :: dirs(2) 
  logical :: bmad_parser_calling = .false.     ! used for expand_lattice
  logical error_flag     ! Needed since bmad_status%ok gets set by many routines.
  logical input_line_meaningful
  logical ran_function_was_called
  logical do_superimpose
  logical write_digested      ! For bmad_parser
  logical write_digested2     ! For bmad_parser2
  logical input_from_file     ! Input is from a lattice file?
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

type (lat_struct)  lat
type (parser_ele_struct), optional :: pele
type (ele_struct), target ::  ele
type (ele_struct), target, save ::  ele0
type (wig_term_struct), pointer :: wig_term(:)
type (real_pointer_struct), allocatable, save :: r_ptrs(:)
type (cross_section_struct), pointer :: cross

real(rp) kx, ky, kz, tol, value, coef
real(rp), pointer :: r_ptr

integer i, j, ix_word, how, ix_word1, ix_word2, ios, ix, i_out
integer expn(6), ix_attrib, i_cross, ix_v

character(40) :: word, str_ix, attrib_word
character(1) delim, delim1, delim2
character(80) str, err_str, line

logical delim_found, err_flag, logic, print_err, set_done, end_of_file
logical, optional :: check_free

! Get next WORD.
! If an overlay or group element then word is just an attribute to control
! [except for a "GROUP[COMMAND] = 0.343" redef construct]

err_flag = .true.  ! assume the worst
call get_next_word (word, ix_word, ':, =()', delim, delim_found, .true.)

! taylor

if (ele%key == taylor$ .and. word(1:1) == '{') then

  word = word(2:)             ! strip off '{'
  read (word, *, iostat = ios) i_out
  if (delim /= ':' .or. ix_word == 0 .or. ios /= 0) then
    call parser_warning ('BAD "OUT" IN TERM FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  call evaluate_value (str, coef, lat, delim, delim_found, err_flag)
  if (err_flag) return

  call get_next_word (line, ix_word, '},', delim, delim_found, .true.)
  read (line, *, iostat = ios) expn
  if (delim /= '}' .or. ix_word == 0 .or. ios /= 0) then
    call parser_warning ('BAD "EXPONENT" IN TERM FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  call add_this_taylor_term (ele, i_out, coef, expn)
  call get_next_word (word, ix_word, '},', delim, delim_found, .true.)

  if (ix_word /= 0 .or. (delim_found .and. delim /= ',')) then
    call parser_warning ('BAD TERM ENDING FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  return
endif

if (ele%key == overlay$) then
  i = attribute_index(ele, word)       ! general attribute search
  if (word == 'B_GRADIENT') call b_grad_warning(ele)

  if (i == type$ .or. i == alias$) then
    call bmad_parser_type_get (ele, word, delim, delim_found)

  else

    if (i < 1) then
      if (print_err) call parser_warning ('BAD OVERLAY ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
      return
    endif

    if (how == def$) then
      ele%ix_value = i
      ele%component_name = word
    endif

    if (ele%ix_value /= i) then
      if (print_err) call parser_warning ('BAD OVERLAY ATTRIBUTE SET FOR: ' // ele%name, &
            'YOU ARE TRYING TO SET: ' // word, &
            'BUT YOU SHOULD BE SETTING: ' // ele%component_name)
      return
    endif

    value = 0
    if (delim == '=') then  ! value
      call evaluate_value (trim(ele%name) // ' ' // word, value, &
                                    lat, delim, delim_found, err_flag)
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
    if (word == 'B_GRADIENT') call b_grad_warning(ele0)
  else   ! how == redef$
    i = attribute_index(ele, word)
    if (word == 'B_GRADIENT') call b_grad_warning(ele)
  endif

  if (i < 1) then
    if (print_err) call parser_warning ('BAD GROUP ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
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
      call parser_warning ('NO VALUE GIVEN FOR ATTRIBUTE FOR: ' // ele%name)
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
    if (print_err) call parser_warning ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  r_ptrs(1)%r = value
  call changed_attribute_bookkeeper (lat, ele, r_ptrs(1)%r)
  return
endif

! Long-range wake

if (word == 'LR' .and. delim == '(') then

  call get_next_word (word, ix_word, '=', delim, delim_found, .true.)
  if (.not. delim_found) then
    call parser_warning ('NO "=" SIGN FOUND', 'FOR ELEMENT: ' // ele%name)
    return
  endif
  call pointer_to_attribute (ele, 'LR(' // word, .false., r_ptr, err_flag, .false.)
  if (err_flag) then
    call parser_warning ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif
  call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
  r_ptr = value
  return
endif

! if not an overlay then see if it is an ordinary attribute.
! if not an ordinary attribute then might be a superimpose switch

if (word(:ix_word) == 'REF') word = 'REFERENCE' ! allowed abbrev
if (word == 'B_GRADIENT') call b_grad_warning (ele)

ix_attrib = attribute_index(ele, word)
attrib_word = word

if (attrib_free_problem(attrib_word)) return

if (ix_word == 0) then  ! no word
  call parser_warning  ('"," NOT FOLLOWED BY ATTRIBUTE NAME FOR: ' // ele%name)
  return
endif
if (ix_attrib < 1) then
  if (print_err) call parser_warning  ('BAD ATTRIBUTE NAME: ' // word, 'FOR ELEMENT: ' // ele%name)
  return
endif

! Capillary cross-section definition

if (ix_attrib == cross$ .and. ele%key == capillary$) then
  if (associated (ele%cross_section)) then
    call re_associate (ele%cross_section, size(ele%cross_section) + 1)
  else
    allocate (ele%cross_section(1))
  endif

  i_cross = size(ele%cross_section)
  cross => ele%cross_section(i_cross)
  ix_v = 0

  if (delim /= '=') then
    call parser_warning ('NO "=" SIGN FOUND AFTER "CROSS"', 'FOR ELEMENT: ' // ele%name)
    return
  endif
  ! Expect "{"
  call get_next_word (word, ix_word, '{,()', delim, delim_found)
  if (delim /= '{') then
    call parser_warning ('NO "{" SIGN FOUND AFTER "CROSS ="', 'FOR ELEMENT: ' // ele%name)
    return
  endif

  do
    ! Expect "s", "v"
    call get_next_word (word, ix_word, '{},()=', delim, delim_found)

    if (word == 's' .and. delim == '=') then
      call evaluate_value (trim(ele%name), cross%s, lat, delim, delim_found, err_flag, ',}')
      if (err_flag) return
      if (delim == '}') exit

    elseif (word == 'v' .and. delim == '(') then
      call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
      ix_v = ix_v + 1
      read (j, '(i)', iostat = ios) word
      if (ios /= 0 .or. ix_v /= j) then
        call parser_warning ('BAD OR OUT OF ORDER CROSS-SECTION VERTEX INDEX NUMBER FOR: ' // ele%name)
        return
      endif

      call get_next_word (word, ix_word, '{},()', delim2, delim_found)
      if (delim /= ')' .or. word /= '=' .or. delim2 /= '{') then        
        call parser_warning ('MALFORMED ORDER CROSS-SECTION VERTEX FOR: ' // ele%name)
        return
      endif

      call evaluate_value (trim(ele%name), cross%v(ix_v)%x, lat, delim, delim_found, err_flag, ',')
      if (err_flag) return

      call evaluate_value (trim(ele%name), cross%v(ix_v)%y, lat, delim, delim_found, err_flag, ',}')
      if (err_flag) return

      if (delim == ',') then
        call evaluate_value (trim(ele%name), cross%v(ix_v)%radius_x, lat, delim, delim_found, err_flag, ',}')
        if (err_flag) return
      endif

      if (delim == ',') then
        call evaluate_value (trim(ele%name), cross%v(ix_v)%radius_y, lat, delim, delim_found, err_flag, ',}')
        if (err_flag) return
      endif

      if (delim == ',') then
        call evaluate_value (trim(ele%name), cross%v(ix_v)%tilt, lat, delim, delim_found, err_flag, '}')
        if (err_flag) return
      endif

      call get_next_word (word, ix_word, '{},()=', delim, delim_found)
      if (word /= '' .or. (delim /= '}' .and. delim /= ',')) then
        call parser_warning ('BAD SYNTAX IN CROSS DEFINITION FOR ELEMENT: ' // ele%name)
        return
      endif
      if (delim == '}') exit

    else
      call parser_warning ('BAD SYNTAX IN CROSS DEFINITION FOR ELEMENT: ' // ele%name)
      return
    endif

  enddo

  call get_next_word (word, ix_word, '{},()=', delim, delim_found)
  if (word /= '' .or. (delim /= ' ' .and. delim /= ',')) then
    call parser_warning ('BAD SYNTAX IN CROSS DEFINITION FOR ELEMENT: ' // ele%name)
    return
  endif

endif

! Capillary s_spline definition

if (ix_attrib == s_spline$ .and. ele%key == capillary$) then

  if (delim /= '=') then
    call parser_warning ('NO "=" SIGN FOUND AFTER "CROSS"', 'FOR ELEMENT: ' // ele%name)
    return
  endif
  ! Expect "{"
  call get_next_word (word, ix_word, '{,()', delim, delim_found)
  if (delim /= '{') then
    call parser_warning ('NO "{" SIGN FOUND AFTER "CROSS ="', 'FOR ELEMENT: ' // ele%name)
    return
  endif

  i_cross = size(ele%cross_section)
  cross => ele%cross_section(i_cross)

  call evaluate_value (trim(ele%name), cross%s_spline(1), lat, delim, delim_found, err_flag, ',}')
  if (err_flag) return

  if (delim == ',') then
    call evaluate_value (trim(ele%name), cross%s_spline(2), lat, delim, delim_found, err_flag, '}')
    if (err_flag) return
  endif

  return
endif

! wiggler term attribute

if (ix_attrib == term$ .and. ele%key == wiggler$) then

  err_flag = .true. ! assume the worst

  if (delim /= '(') then   ! ) then
    call parser_warning ('"TERM" FOR A WIGGLER NOT FOLLOWED BY A "(" FOR: ' // ele%name)  ! )
    return
  endif

  call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.) ! (
  if (delim /= ')') then
    call parser_warning ('CANNOT FIND CLOSING ")" for a "TERM(i)" FOR A WIGGLER"', 'FOR: ' // ele%name)
    return
  endif

  read (word, *, iostat = ios) ix
  if (ix < 1 .or. ios /= 0) then
    call parser_warning ('BAD TERM NUMBER FOR A WIGGLER FOR: ' // ele%name)
    return
  endif

  write (str_ix, '(a, i3, a)') 'TERM(', ix, ')'

  if (.not. associated(ele%wig_term)) then
    allocate(ele%wig_term(ix))
  elseif (ix > size(ele%wig_term)) then
    allocate (wig_term(size(ele%wig_term)))
    wig_term = ele%wig_term
    deallocate (ele%wig_term)
    allocate (ele%wig_term(ix))
    ele%wig_term(1:size(wig_term)) = wig_term
    deallocate (wig_term)
  endif

! 1) chop "=", 2) chop to "{", 3) chop to "}", 4) chop to "," if it exists

  call get_next_word (word, ix_word1, ':,={}', delim1, delim_found, .true.)  
  call get_next_word (word, ix_word2, ':,={}', delim2, delim_found, .true.)  

  if (delim1 /= '=' .or. delim2 /= '{' .or. ix_word1 /= 0 .or. ix_word2 /= 0) then
    call parser_warning ('CONFUSED SYNTAX FOR TERM IN WIGGLER: ' // ele%name, str_ix)
    return
  endif

  err_str = trim(ele%name) // ' ' // str_ix

  call evaluate_value (err_str, ele%wig_term(ix)%coef, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return
 
  call evaluate_value (err_str, ele%wig_term(ix)%kx, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return

  call evaluate_value (err_str, ele%wig_term(ix)%ky, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return

  call evaluate_value (err_str, ele%wig_term(ix)%kz, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return

  call evaluate_value (err_str, ele%wig_term(ix)%phi_z, lat, delim, delim_found, err_flag, '}')
  if (err_flag) return

  kx = ele%wig_term(ix)%kx
  ky = ele%wig_term(ix)%ky
  kz = ele%wig_term(ix)%kz
  tol = 1e-5 * (kx**2 + ky**2 + kz**2)

  if (abs(ky**2 - kx**2 - kz**2) < tol) then
    ele%wig_term(ix)%type = hyper_y$
  elseif (abs(ky**2 + kx**2 - kz**2) < tol) then
    ele%wig_term(ix)%type = hyper_xy$
  elseif (abs(ky**2 - kx**2 + kz**2) < tol) then
    ele%wig_term(ix)%type = hyper_x$
  else
    call parser_warning ('WIGGLER TERM DOES NOT HAVE CONSISTANT Kx, Ky, and Kz', &
                  'FOR WIGGLER: ' // ele%name // '  ' // str_ix)
    err_flag = .true.
    return
  endif

  call get_next_word (word, ix_word,  ':,=()', delim,  delim_found, .true.)  
  if (ix_word /= 0) then
    call parser_warning ('BAD SYNTAX FOR WIGGLER: ' // ele%name, str_ix)
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
        call parser_warning ('SORRY I''M NOT PROGRAMMED TO USE A "TILT" DEFAULT' // &
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
    if (.not. present(pele)) call parser_warning ('INTERNAL ERROR...')
    pele%ref_pt = begin$

  case (ref_center$)
    if (.not. present(pele)) call parser_warning ('INTERNAL ERROR...')
    pele%ref_pt = center$

  case (ref_end$)
    if (.not. present(pele)) call parser_warning ('INTERNAL ERROR...')
    pele%ref_pt = end$

  case (ele_beginning$)
    if (.not. present(pele)) call parser_warning ('INTERNAL ERROR...')
    pele%ele_pt = begin$

  case (ele_center$)
    if (.not. present(pele)) call parser_warning ('INTERNAL ERROR...')
    pele%ele_pt = center$

  case (ele_end$)
    if (.not. present(pele)) call parser_warning ('INTERNAL ERROR...')
    pele%ele_pt = end$

  case (common_lord$)
    if (.not. present(pele)) call parser_warning ('INTERNAL ERROR...')
    pele%common_lord = .true.

  case default
    call parser_warning ('EXPECTING "=" AFTER ATTRIBUTE: ' // word,  'FOR ELEMENT: ' // ele%name)
    err_flag = .true.
  end select

  return
endif

! get the value of the attribute.
! The TYPE, ALIAS, and DESCRIP attributes are special because their "values"
! are character strings

select case (attrib_word)

case ('REFERENCE')
  if (.not. present(pele)) call parser_warning ('INTERNAL ERROR...')
  call get_next_word(pele%ref_name, ix_word,  ':=,', delim, delim_found, .true.)

case ('OFFSET')
  call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
  if (err_flag) return
  if (.not. present(pele)) call parser_warning ('INTERNAL ERROR...')
  pele%s = value

case('TYPE', 'ALIAS', 'DESCRIP', 'SR_WAKE_FILE', 'LR_WAKE_FILE', 'LATTICE', 'TO', &
     'REF_PATCH', 'CRYSTAL_TYPE')
  call bmad_parser_type_get (ele, attrib_word, delim, delim_found)

case ('SYMPLECTIFY') 
  if (how == def$ .and. (delim == ',' .or. .not. delim_found)) then
    ele%symplectify = .true.
  else
    call get_logical ('SYMPLECTIFY', ele%symplectify)
    if (ios /= 0 .or. ix_word == 0) return
  endif
  
case ('IS_ON')
  call get_logical ('IS_ON', ele%is_on)
  if (ios /= 0 .or. ix_word == 0) return

case ('MATCH_END')
  call get_logical ('MATCH_END', logic)
  if (logic) then; ele%value(match_end$) = 1
  else;            ele%value(match_end$) = 0
  endif
  if (ios /= 0 .or. ix_word == 0) return

case ('PATCH_END')
  call get_logical ('PATCH_END', logic)
  if (logic) then; ele%value(patch_end$) = 1
  else;            ele%value(patch_end$) = 0
  endif
  if (ios /= 0 .or. ix_word == 0) return

case ('TRANSLATE_AFTER')
  call get_logical ('TRANSLATE_AFTER', logic)
  if (logic) then; ele%value(translate_after$) = 1
  else;            ele%value(translate_after$) = 0
  endif
  if (ios /= 0 .or. ix_word == 0) return

case ('MATCH_END_ORBIT')
  call get_logical ('MATCH_END_ORBIT', logic)
  if (logic) then; ele%value(match_end_orbit$) = 1
  else;            ele%value(match_end_orbit$) = 0
  endif
  if (ios /= 0 .or. ix_word == 0) return

case ('APERTURE_LIMIT_ON') 
  call get_logical ('APERTURE_LIMIT_ON', lat%param%aperture_limit_on)
  if (ios /= 0 .or. ix_word == 0) return

case ('CSR_CALC_ON')
  call get_logical ('CSR_CALC_ON', ele%csr_calc_on)
  if (ios /= 0 .or. ix_word == 0) return

case ('MAP_WITH_OFFSETS')
  call get_logical ('MAP_WITH_OFFSETS', ele%map_with_offsets)
  if (ios /= 0 .or. ix_word == 0) return

case ('OFFSET_MOVES_APERTURE')
  call get_logical ('OFFSET_MOVES_APERTURE', ele%offset_moves_aperture)
  if (ios /= 0 .or. ix_word == 0) return

case ('FIELD_MASTER')
  call get_logical ('FIELD_MASTER', ele%field_master)
  if (ios /= 0 .or. ix_word == 0) return

case ('SCALE_MULTIPOLES')
  call get_logical ('SCALE_MULTIPOLES', ele%scale_multipoles)
  if (ios /= 0 .or. ix_word == 0) return

case ('DIFFRACTION_TYPE')
  call match_word (attrib_word, diffraction_type_name(1:), ele%sub_key)
  if (ele%sub_key < 1) call parser_warning ('BAD DIFFRACTION_TYPE ATTRIBUTE.')

case default   ! normal attribute

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
  elseif (ix_attrib == field_calc$) then
    ele%field_calc = nint(value)
  elseif (ix_attrib == tracking_method$) then
    ele%tracking_method = nint(value)
  elseif (ix_attrib == ref_orbit$) then
    ele%ref_orbit = nint(value)
  elseif (ix_attrib == aperture_at$) then
    ele%aperture_at = nint(value)
  elseif (ix_attrib == aperture_type$) then
    ele%aperture_type = nint(value)
  elseif (ix_attrib == ran_seed$) then
    call ran_seed_put (nint(value))  ! init random number generator
    bp_com%ran_function_was_called = .true.
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
    if (ele%key == custom$ .and. ix_attrib == l$) ele%value(l_original$) = value

    select case (ix_attrib)

    case (num_steps$)
      ele%value(ds_step$) = abs(ele%value(l$) * nint(ele%value(num_steps$)))

    case (e_tot$)
      select case (ele%key)
      case (def_beam$)
        lat%ele(0)%value(e_tot$) = 1d9 * value
        bp_com%param_ele%value(e_tot$) = 1d9 * value
      case (def_parameter$, init_ele$)
        lat%ele(0)%value(e_tot$) = value
        bp_com%param_ele%value(e_tot$) = value
        bp_com%beam_ele%value(e_tot$) = 1d-9 * value
      end select
      bp_com%e_tot_set = .true.
      bp_com%p0c_set   = .false.

    case (lr_freq_spread$)
      call randomize_lr_wake_frequencies (ele, set_done)
      if (set_done) bp_com%ran_function_was_called = .true.

    case (p0c$)
      select case (ele%key)
      case (def_beam$)
        lat%ele(0)%value(p0c$) = 1d9 * value
        bp_com%param_ele%value(p0c$) = 1d9 * value
      case (def_parameter$, init_ele$)
        lat%ele(0)%value(p0c$) = value
        bp_com%param_ele%value(p0c$) = value
        bp_com%beam_ele%value(p0c$) = 1d-9 * value
      end select
      bp_com%e_tot_set = .false.
      bp_com%p0c_set   = .true.
    end select

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
    if (print_err) call parser_warning ('ATTRIBUTE NOT FREE TO BE SET: ' // attrib_name, 'FOR: ' // ele%name)
    err_flag = .true.
    is_problem = .true.
  endif
endif

end function attrib_free_problem

!--------------------------------------------------------
! contains

subroutine get_logical (name, this_logic)

character(*) name
logical this_logic

!

call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.)
this_logic = evaluate_logical (word, ios)
if (ios /= 0 .or. ix_word == 0) then
  call parser_warning ('BAD "' // trim(name) // '" SWITCH FOR: ' // ele%name)
endif

end subroutine get_logical

!--------------------------------------------------------
! contains

subroutine b_grad_warning (ele)

type (ele_struct) ele
integer kk

!

do kk = 1, 20
  print *, trim(bp_com%parser_name), ' WARNING!'
enddo

if (ele%key == group$ .or. ele%key == overlay$) then
  print *, '      B_GRADIENT NEEDS TO BE REPLACED BY B1_GRADIENT, B2_GRADIENT, OR B3_GRAIENT.'
else
  print *, '      B_GRADIENT NEEDS TO BE REPLACED BY: ', trim(attribute_name(ele, b_gradient$))
endif

if (bp_com%current_file%full_name /= ' ') then
  if (bp_com%input_line_meaningful) then
    print *, '      IN FILE: ', trim(bp_com%current_file%full_name)
    print *, '      AT OR BEFORE LINE:', bp_com%current_file%i_line
  else
    print *, '      ROOT FILE: ', trim(bp_com%current_file%full_name)
  endif
endif

do kk = 1, 10
  print *, 'NO DIGESTED FILE WILL BE MADE BECAUSE OF THIS!'
enddo
bp_com%write_digested = .false.
bp_com%write_digested2 = .false.

end subroutine

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

if (delim /= ',')  call parser_warning ('"CALL" NOT FOLLOWED BY COMMA', stop_here = .true.)
call get_next_word(call_file, ix_word, ':=,', delim, delim_found, .true.)

if (ix_word == 0) then
  call parser_warning ('NOTHING AFTER "CALL"', stop_here = .true.)
  return
elseif (index('FILENAME', call_file(:ix_word)) /= 1) then
  call parser_warning ('INVALID "CALL" COMMAND', stop_here = .true.)
  return
elseif (delim /= '=') then
  call parser_warning ('NO "=" AFTER "FILENAME"', stop_here = .true.)
  return
endif

call get_next_word(call_file, ix_word, ',', delim, delim_found, .false.)
if (ix_word == 0) then
  call parser_warning ('NO FILE NAME SPECIFIED', stop_here = .true.)
  return
endif

if (call_file(1:1) == '"') then
  call_file = call_file(2:)
  ix = index(call_file, '"')
  if (ix == 0 .or. ix /= len_trim(call_file)) then
    call parser_warning ('MISSING DOUBLE QUOTE MARK (") FOR CALL STATEMENT', stop_here = .true.)
    return
  endif
  call_file(ix:ix) = ' '
endif

if (call_file(1:1) == "'") then
  call_file = call_file(2:)
  ix = index(call_file, "'")
  if (ix == 0 .or. ix /= len_trim(call_file)) then
    call parser_warning ("MISSING SINGLE QUOTE MARK (') FOR CALL STATEMENT", stop_here = .true.)
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
call file_stack ('push', call_file, finished, err) ! err gets set here

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
  call parser_warning ('"OUT" VALUE IN TAYLOR TERM NOT IN RANGE (1 - 6)', &
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
!
! Output
!   ix_word     -- Integer: length of word argument
!   delim       -- Character1: Actual delimiter found
!   delim_found -- Logical: Set true if a delimiter found. A delimiter
!                    may not be found if the end of the line is reached first.
!-


subroutine get_next_word (word, ix_word, delim_list, delim, delim_found, upper_case_word)

implicit none

integer ix_a, ix_word

character(*) word, delim_list, delim

integer n

logical delim_found, end_of_file
logical, optional :: upper_case_word

! check for continuation character and, if found, then load more characters
! into the parse line from the lattice file. 
! If the input is not from a file then skip this.

if (bp_com%input_from_file) then 
  n = len_trim(bp_com%parse_line)
  if (n > 0 .and. n < len(bp_com%parse_line)/2) then
    select case (bp_com%parse_line(n:n))
    case (',', '+', '-', '*', '/', '(', '{', '[', '=')
      call load_parse_line('continue', n+2, end_of_file)
    case ('&')
      call load_parse_line('continue', n, end_of_file)
    end select
  endif
endif

! Get the first word in bp_com%parse_line

call word_read (bp_com%parse_line, delim_list,  &
                       word, ix_word, delim, delim_found, bp_com%parse_line)

if (len(word) < ix_word) then
  call parser_warning ('BAD WORD: ' // bp_com%parse_line)
  ix_word = len(word)
endif

if (present(upper_case_word)) then
  if (upper_case_word) call str_upcase (word, word)
else
  call str_upcase (word, word)
endif

! Note: "var := num" is old-style variable definition syntax.
! If delim is ":" and next char is "=" then use "=" as the delim

if (delim == ':' .and. index(delim_list, '=') /= 0 .and. &
                                        bp_com%parse_line(1:1) == '=') then
  delim = '='
  bp_com%parse_line = bp_com%parse_line(2:)
endif

end subroutine get_next_word

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine file_stack (how, file_name_in, finished, err)
!
! Subroutine to keep track of the files that are opened for reading.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine file_stack (how, file_name_in, finished, err)

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
  return
endif

! "push" means open a file and put its name on the stack.

finished = .false.

if (how == 'push') then

  i_level = i_level + 1    ! number of files currently open
  if (i_level > f_maxx) then
    print *, 'ERROR: CALL NESTING GREATER THAN 20 LEVELS'
    call err_exit
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
    call parser_warning ('MALFORMED FILE NAME: ' // file_name_in, stop_here = .true.)
    if (bmad_status%exit_on_error) call err_exit
    do i = 1, i_level-1
      close (file(i_level)%f_unit)
    enddo
    return
  endif

  ix = splitfilename (file_name2, file(i_level)%dir, basename, is_relative)
  if (is_relative) then
    call append_subdirectory (trim(file(i_level-1)%dir), &
                                           file(i_level)%dir, file(i_level)%dir, err_flag)
    if (err_flag) call parser_warning ('BAD DIRECTORY SYNTAX FOR: ' // file_name, stop_here = .true.)
  endif
  bp_com%dirs(1) = file(i_level-1)%dir
  call find_file (file_name2, found_it, file_name, bp_com%dirs)
  file(i_level)%logical_name = file_name2
  file(i_level)%full_name = file_name
  file(i_level)%f_unit = lunget()

  open (file(i_level)%f_unit, file = file_name,  &
                               status = 'OLD', action = 'READ', iostat = ios)
  if (ios /= 0 .or. .not. found_it) then
    bp_com%current_file => file(i_level-1)  ! For warning
    if (file_name2 == file_name)  then
      call parser_warning ('UNABLE TO OPEN FILE: ' // file_name, stop_here = .true.)
    else
      call parser_warning ('UNABLE TO OPEN FILE: ' // file_name, &
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

elseif (how == 'pop') then
  close (unit = bp_com%current_file%f_unit)
  i_level = i_level - 1
  if (i_level < 0) then
    call parser_warning ('BAD "RETURN"')
    return
  elseif (i_level > 0) then
    bp_com%current_file => file(i_level)
    bp_com%calling_file => file(i_level-1)
  else    ! i_level == 0
    finished = .true.
  endif
else
  print *, 'BMAD_PARSER: INTERNAL ERROR IN FILE_STACK SUBROUTINE'
  call err_exit
endif

if (present(err)) err = .false.

end subroutine file_stack

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine load_parse_line (load_type, ix_start, end_of_file) 
!
! Subroutine to load characters from the input file.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine load_parse_line (load_type, ix_start, end_of_file)

implicit none

integer ix_start, ix

character(*) load_type
character(n_parse_line+20), save :: line, saved_line

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
      call parser_warning ('INPUT LINE HAS TOO MANY CHARACTERS:', line)
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
    call error_exit ('INTERNAL ERROR #4: CALL HELP')    
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

bp_com%parse_line(ix_start:) = line

return

!

9000  continue
end_of_file = .true.
bp_com%parse_line = ' '

end subroutine load_parse_line

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
!   iostat -- Integer: Status: Returns 0 if conversion successful. 
!
! Output:
!   this_logic -- Logical: Result.
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

logical delim_found, split, ran_function_pending
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

! parsing loop to build up the stack.

parsing_loop: do

! get a word

  call get_next_word (word, ix_word, '+-*/()^,:} ', delim, delim_found)

  if (delim == '*' .and. word(1:1) == '*') then
    call parser_warning ('EXPONENTIATION SYMBOL IS "^" AS OPPOSED TO "**"!',  &
                  'for: ' // err_str)
    return
  endif

  if (ran_function_pending .and. (ix_word /= 0 .or. delim /= ')')) then
    call parser_warning ('RAN AND RAN_GAUSS DO NOT TAKE AN ARGUMENT', 'FOR: ' // err_str)
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
        bp_com%ran_function_was_called = .true.
      case ('RAN_GAUSS') 
        call pushit (op, i_op, ran_gauss$)
        ran_function_pending = .true.
        bp_com%ran_function_was_called = .true.
      case default
        call parser_warning ('UNEXPECTED CHARACTERS ON RHS BEFORE "(": ' // word,  &
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
        call parser_warning  ('CONSTANT OR VARIABLE MISSING BEFORE ")"', 'FOR: ' // err_str)
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
        call parser_warning ('UNMATCHED ")" ON RHS', 'FOR: ' // err_str)
        return
      endif

      i_op = i - 1

      call get_next_word (word, ix_word, '+-*/()^,:}', delim, delim_found)
      if (ix_word /= 0) then
        call parser_warning ('UNEXPECTED CHARACTERS ON RHS AFTER ")"',  &
                                                  'FOR: ' // err_str)
        return
      endif

      if (delim /= ')') exit  ! if no more ')' then no need to release more
    enddo


    if (delim == '(') then
      call parser_warning ('")(" CONSTRUCT DOES NOT MAKE SENSE FOR: ' // err_str)
      return
    endif

! For binary "+-/*^" delims

  else
    if (ix_word == 0) then
      call parser_warning ('CONSTANT OR VARIABLE MISSING IN EVALUATING: ' // err_str)
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
    call parser_warning ('MALFORMED EXPRESSION')
    bp_com%parse_line = ' '
    return
  end select

! now see if there are operations on the OP stack that need to be transferred
! to the STK stack

  do i = i_op, 1, -1
    if (eval_level(op(i)) >= eval_level(i_delim)) then
      if (op(i) == l_parens$) then
        call parser_warning ('UNMATCHED "(" IN EVALUATING: ' // err_str)
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
  call parser_warning ('UNMATCHED "(" IN EVALUATING: ' // err_str)
  return
endif

if (i_lev == 0) call parser_warning ('NO VALUE FOUND FOR: ' // err_str)

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
      call parser_warning ('DIVIDE BY 0 ON RHS', 'FOR: ' // err_str)
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
    call error_exit ('INTERNAL ERROR #02: GET HELP')
  endif
enddo

if (i2 /= 1) call error_exit ('INTERNAL ERROR #03: GET HELP')

value = stk(1)%value
err_flag = .false.

! Check that final delim matches.

if (present(end_delims)) then
  if (.not. delim_found .or. index(end_delims, delim) == 0) then
    call parser_warning ('BAD DELIMITOR AFTER VALUE FOR: ' // err_str)
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
  print *, 'ERROR IN ', trim(bp_com%parser_name), ': STACK OVERFLOW.'
  print *, '      EXPERT HELP IS NEEDED!'
  call err_exit
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

integer i, ix1, ix2, ix_word, ios, ix
real(rp) value
character(*) word
character(40) attrib_name, ele_name
logical err_flag

! see if this is numeric

if (index('-+.0123456789', word(1:1)) /= 0) then
  read (word, *, iostat = ios) value
  if (ios /= 0) call parser_warning ('BAD VARIABLE: ' // word)
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
    call parser_warning ('VARIABLE USED BUT NOT YET DEFINED: ' // word)
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

if (attrib_name == 'S' .and. bp_com%parser_name /= 'BMAD_PARSER2') then
  call parser_warning ('"S" ATTRIBUTE CAN ONLY BE USED WITH BMAD_PARSER2')
endif

call pointers_to_attribute (lat, ele_name, attrib_name, .false., ptr, err_flag, .false.)
if (err_flag .or. size(ptr) == 0) then
  call parser_warning('BAD ATTRIBUTE: ' // word)
else
  value = ptr(1)%r
endif

! If size(ptr) > 1 then there must be more than one element of the same name.

if (size(ptr) > 1) call parser_warning (&
            'MULTIPLE ELEMENTS OF THE SAME NAME REFERENCED IN ATTRIBUTE: ' // word)

return

end subroutine

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
  call parser_warning ('VARIABLES ARE NOT ALLOWED TO BE REDEFINED: ' // word)
  call evaluate_value (word, bp_com%var(i)%value, lat, delim, delim_found, err_flag)
  return
endif

bp_com%ivar_tot = bp_com%ivar_tot + 1
if (bp_com%ivar_tot > size(bp_com%var%name)) call reallocate_bp_com_var()
ivar = bp_com%ivar_tot
bp_com%var(ivar)%name = word
call evaluate_value (bp_com%var(ivar)%name, bp_com%var(ivar)%value, &
                                     lat, delim, delim_found, err_flag)
if (delim_found .and. .not. err_flag) call parser_warning  &
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
    call parser_warning ('MISSING ENDING QUOTE MARK FOR TYPE = "attribute"',  &
                        'FOR ELEMENT: ' // ele%name)
    type_name = ' '
  else
    type_name = bp_com%parse_line(1:ix-1)
    bp_com%parse_line = bp_com%parse_line(ix+1:)
    call get_next_word (word, ix_word, ',=', delim, delim_found, .true.)
    if (ix_word /= 0) call parser_warning (  &
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
  call upcase_string (ele%component_name)
case default
  print *, 'INTERNAL ERROR IN BMAD_PARSER_TYPE_GET: I NEED HELP!'
  call err_exit
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
!     %wake%lr(:)       -- Long-range wake potential.
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

if (.not. associated(ele%wake)) allocate (ele%wake)
if (.not. associated(ele%wake%sr_table))       allocate (ele%wake%sr_table(0))
if (.not. associated(ele%wake%sr_mode_long))  allocate (ele%wake%sr_mode_long(0))
if (.not. associated(ele%wake%sr_mode_trans)) allocate (ele%wake%sr_mode_trans(0))
if (associated(ele%wake%lr)) deallocate (ele%wake%lr)

! get data

call find_this_file (iu, lr_file_name, full_file_name)
if (iu < 0) return

ele%wake%lr_file = lr_file_name

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
  call parser_warning ('CANNOT READ LONG_RANGE_MODES NAMELIST FOR ELEMENT: ' // ele%name, & 
                'FROM FILE: '// full_file_name)
  return
endif

n_row = count(lr%freq /= -1)
allocate (ele%wake%lr(n_row))
j = 0
do i = 1, size(lr)
  if (lr(i)%freq == -1) cycle

  j = j + 1
  ele%wake%lr(j)%freq_in   = lr(i)%freq
  ele%wake%lr(j)%freq      = lr(i)%freq
  ele%wake%lr(j)%r_over_q  = lr(i)%r_over_q
  ele%wake%lr(j)%q         = lr(i)%q
  ele%wake%lr(j)%m         = lr(i)%m
  ele%wake%lr(j)%b_sin     = lr(i)%b_sin
  ele%wake%lr(j)%b_cos     = lr(i)%b_cos
  ele%wake%lr(j)%a_sin     = lr(i)%a_sin
  ele%wake%lr(j)%a_cos     = lr(i)%a_cos
  ele%wake%lr(j)%t_ref     = lr(i)%t_ref

  call downcase_string(lr(i)%angle)
  if (lr(i)%angle == '') then
    call parser_warning ('LONG_RANGE_MODE ANGLE IS MISSING. MUST BE NUMBER OR "UNPOLARIZED"', & 
                  'FOR ELEMENT: ' // ele%name, &
                  'IN FILE: ' // full_file_name)
    cycle
  endif

  if (index('unpolarized', trim(lr(j)%angle)) == 1) then
    ele%wake%lr(j)%polarized = .false.
    ele%wake%lr(j)%angle     = 0
  else
    ele%wake%lr(j)%polarized = .true.
    read (lr(j)%angle, *, iostat = ios) ele%wake%lr(j)%angle
    if (ios /= 0) then
      call parser_warning ('BAD LONG_RANGE_MODE ANGLE.', &
                    'FOR ELEMENT: ' // ele%name, &
                    'IN FILE: ' // full_file_name)
      cycle
    endif
  endif
enddo

call randomize_lr_wake_frequencies (ele, set_done)
if (set_done) bp_com%ran_function_was_called = .true.
  
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
!     %wake%sr_table(:)       -- Short-range wake potential.
!     %wake%sr_mode_long(:)  -- Short-range wake potential.
!     %wake%sr_mode_trans(:) -- Short-range wake potential.
!-
      
subroutine read_sr_wake (ele, sr_file_name)

implicit none

type (ele_struct) ele
type (sr_mode_wake_struct) longitudinal(100), transverse(100)

real(rp) dz, z_max
real(rp), allocatable :: col1(:), col2(:), col3(:)
integer n_row, n, j, iu, ios, ix, i

character(*) sr_file_name
character(80) line
character(200) full_file_name

logical found_it

namelist / short_range_modes / z_max, longitudinal, transverse

! init

if (.not. associated(ele%wake)) allocate (ele%wake)
if (.not. associated(ele%wake%lr)) allocate (ele%wake%lr(0))
if (associated(ele%wake%sr_table))       deallocate (ele%wake%sr_table)
if (associated(ele%wake%sr_mode_long))  deallocate (ele%wake%sr_mode_long)
if (associated(ele%wake%sr_mode_trans)) deallocate (ele%wake%sr_mode_trans)

allocate (ele%wake%sr_table(0), ele%wake%sr_mode_long(0), ele%wake%sr_mode_trans(0))

! get sr_table data

iu = 0
ele%wake%sr_file = sr_file_name
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
    call parser_warning ('ERROR READING WAKE FILE: ' // full_file_name)
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
    call parser_warning ('ERROR PARSING WAKE FILE: ' // full_file_name, &
                                        'CANNOT READ LINE: ' // line)
    return
  endif

enddo

allocate (ele%wake%sr_table(0:n_row-1))
ele%wake%sr_table%z     = col1(1:n_row)
ele%wake%sr_table%long  = col2(1:n_row)
ele%wake%sr_table%trans = col3(1:n_row)

deallocate (col1, col2, col3)

! err check

if (n_row > 1) then
  if (ele%wake%sr_table(0)%z /= 0) then
    call parser_warning ('WAKEFIELDS DO NOT START AT Z = 0!', &
                                  'IN FILE: ' // ele%wake%sr_file)
    return
  endif

  n = n_row - 1
  dz = ele%wake%sr_table(n)%z / n

  do j = 1, n
    if (abs(ele%wake%sr_table(j)%z - dz * j) > 1e-4 * abs(dz)) then
      write (line, '(a, i5)') &
               'WAKEFIELD POINTS DO NOT HAVE UNIFORM DZ FOR POINT:', j
      call parser_warning (line, 'IN FILE: ' // ele%wake%sr_file)
      return
    endif
  enddo               

  ! if dz > 0 means that an old-style file is being used.

  if (dz > 0) call parser_warning ( &
          'SHORT-RANGE WAKEFIELD FILE TABLES NOW MUST HAVE Z < 0! ' // full_file_name, &
          'REMEMBER THAT Wt NEEDS TO BE NEGATIVE ALSO!')

  if (ele%wake%sr_table(1)%trans > 0) call parser_warning ( &
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
  call parser_warning ('CANNOT READ SHORT_RANGE_MODES NAMELIST FROM FILE: ' & 
                    // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

n = count(longitudinal%phi /= real_garbage$)
allocate (ele%wake%sr_mode_long(n))
ele%wake%sr_mode_long = longitudinal(1:n)
if (any(longitudinal(1:n)%phi == real_garbage$)) call parser_warning ( &
    'JUMBLED INDEX FOR LONGITUDINAL SHORT_RANGE_MODES FROM FILE: ' &
    // full_file_name, 'FOR ELEMENT: ' // ele%name)

n = count(transverse%phi /= real_garbage$)
allocate (ele%wake%sr_mode_trans(n))
ele%wake%sr_mode_trans = transverse(1:n)
if (any(transverse(1:n)%phi == real_garbage$)) call parser_warning ( &
    'JUMBLED INDEX FOR TRANSVERSE SHORT_RANGE_MODES FROM FILE: ' &
    // full_file_name, 'FOR ELEMENT: ' // ele%name)


ele%wake%z_sr_mode_max = z_max
if (z_max == real_garbage$) call parser_warning ( &
    'Z_MAX NOT SET FOR SHORT_RANGE_MODES FROM FILE: ' &
    // full_file_name, 'FOR ELEMENT: ' // ele%name)

end subroutine

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
call parser_warning ('CANNOT OPEN WAKE FILE: ' // file)
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
end subroutine

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
if (delim /= '{' .or. ix_word /= 0) call parser_warning  &
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
      call parser_warning ('BAD ATTRIBUTE SPEC: ' // word_in, 'FOR: ' // ele%name)
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
    call evaluate_value (trim(ele%name), value, &
                                lat, delim, delim_found, err_flag)
    if (err_flag) then
      call parser_warning ('BAD COEFFICIENT: ' // word_in,  &
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
    call parser_warning ('BAD ' // control_name(ele%lord_status) //  &
            'SPEC: ' // word_in, 'FOR: ' // ele%name)
    exit
  endif
                        
enddo

!

ixs = ele%n_slave

! if (ixs == 0) call parser_warning ( &
!        'NO SLAVE ELEMENTS ASSOCIATED WITH GROUP/OVERLAY ELEMENT: ' // ele%name)

allocate (pele%coef(ixs), pele%name(ixs), pele%attrib_name(ixs))
pele%coef = coef(1:ixs)
pele%name = name(1:ixs)
pele%attrib_name = attrib_name(1:ixs)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-


subroutine verify_valid_name (name, ix_name, wildcards_permitted)

implicit none

integer i, ix_name, ix1, ix2

character(*) name
character(27), parameter :: letters = '\ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
character(44), parameter :: valid_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\0123456789_[]().#'
character(1), parameter :: tab = achar(9)

logical, optional :: wildcards_permitted
logical OK, wild_permit

! init

wild_permit = logic_option (.false., wildcards_permitted)

! check for blank spaces

do i = 1, min(ix_name, len(name))
  if (name(i:i) == ' ' .or. name(i:i) == tab) call parser_warning  &
                        ('NO DELIMITER BETWEEN NAMES: ' // name)
enddo

! check for name too long

if (ix_name > len(name)) then
   call parser_warning ('NAME TOO LONG: ' // name)
   ix_name = len(name)      ! chop name
endif

! check for name too short

if (ix_name == 0) call parser_warning ('BLANK NAME')

! check for invalid characters in name

OK = .true.
if (index(letters, name(1:1)) == 0 .and. &
          .not. (wild_permit .and. (name(1:1) == '*' .or. name(1:1) == '%'))) OK = .false.

do i = 1, ix_name
  if (index(valid_chars, name(i:i)) == 0 .and. &
         .not. (wild_permit .and. (name(i:i) == '*' .or. name(i:i) == '%'))) OK = .false.
enddo

if (.not. OK) call parser_warning ('INVALID NAME: UNRECOGNIZED CHARACTERS IN: ' // name)

! Check for non-matched "(" ")" pairs

ix1 = index(name, '(')
ix2 = index(name, ')')
if (ix1 /= 0 .or. ix2 /= 0) then
  if (ix1 == 0) call parser_warning ('UNMATCHED PARENTHESIS: ' // name)
  if (ix2 <= ix1+1) call parser_warning  ('INVALID: REVERSED PARENTHESES: ' // name)
  if (index(name(ix1+1:), '(') /= 0 .or. index(name(ix2+1:), ')') /=  &
                 0) call parser_warning ('INVALID: BAD PARENTHESES: ' // name)
endif

! Check for non matched "[" "]" pairs

ix1 = index(name, '[')
ix2 = index(name, ']')
if (ix1 /= 0 .or. ix2 /= 0) then
  if (ix1 == 0) call parser_warning ('UNMATCHED BRACKET: ' // name)
  if (ix2 <= ix1+1) call parser_warning  ('INVALID: REVERSED BRACKETS: ' // name)
  if (index(name(ix1+1:), '[') /= 0 .or. index(name(ix2+1:), ']') /=  &
                 0) call parser_warning ('INVALID: BAD BRACKETS: ' // name)
  if (ix2 /= len(name)) then
    if (name(ix2+1:ix2+1) /= ' ') call parser_warning  &
                  ('INVALID: SOMETHING AFTER CLOSING "]" BRACKET: ' // name)
  endif
endif

! check for more than 40 characters

if ((ix1 == 0 .and. ix_name > 40) .or. (ix1 > 41 .or. ix2 - ix1 > 41)) &
                          call parser_warning ('NAME HAS > 40 CHARACTERS: ' // name)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine error_exit (what1, what2)

implicit none

character(*) what1
character(*), optional :: what2

!

print *, 'FATAL ERROR IN ', trim(bp_com%parser_name) 
print '(5x, a)', trim(what1)
if (present(what2)) then
  if (what2 /= ' ') print '(5x, a)', trim(what2)
endif

if (bp_com%input_line_meaningful) then
  print *, '      IN FILE: ', trim(bp_com%current_file%full_name)
  print *, '      AT OR BEFORE LINE:', bp_com%current_file%i_line
endif

call err_exit

end subroutine


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_warning (what1, what2, what3, seq, pele, stop_here)

implicit none

type (seq_struct), optional :: seq
type (parser_ele_struct), optional :: pele

character(*) what1
character(*), optional :: what2, what3

logical, optional :: stop_here

! BP_COM%ERROR_FLAG is a common logical used so program will stop at end of parsing

if (bmad_status%type_out) then

  print *, 'ERROR IN ', trim(bp_com%parser_name), ': ', trim(what1)

  if (present(what2)) print '(22x, a)', trim(what2)
  if (present(what3)) print '(22x, a)', trim(what3)

  if (present(seq)) then
    print *, '      IN FILE: ', trim(seq%file_name)
    print *, '      AT LINE:', seq%ix_line
  elseif (bp_com%current_file%full_name /= ' ') then
    if (bp_com%input_line_meaningful) then
      print *, '      IN FILE: ', trim(bp_com%current_file%full_name)
      print *, '      AT OR BEFORE LINE:', bp_com%current_file%i_line
    else
      print *, '      ROOT FILE: ', trim(bp_com%current_file%full_name)
    endif
  endif

  if (bp_com%input_line_meaningful) then
     if (len_trim(bp_com%input_line1) /= 0) print '(5x, a)', trim(bp_com%input_line1)
     if (len_trim(bp_com%input_line2) /= 0) print '(5x, a)', trim(bp_com%input_line2)
  endif

  if (present(pele)) then
    print *, '      ELEMENT DEFINED IN FILE: ', trim(pele%lat_file)
    print *, '      AT LINE: ', pele%ix_line_in_file
  endif

  print *

endif

bp_com%error_flag = .true.
bmad_status%ok = .false.

if (logic_option(.false., stop_here) .and. bmad_status%exit_on_error) stop

end subroutine parser_warning

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
           ubound(ref_orbit_name, 1) + ubound(element_end_name, 1) + ubound(shape_name, 1) + &
           size(particle_name)
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

do i = 1, ubound(shape_name, 1)
  nn = nn + 1
  call str_upcase (bp_com%var(nn)%name, shape_name(i))
  bp_com%var(nn)%value = i
enddo

call indexx (bp_com%var(1:nt)%name, bp_com%var(1:nt)%indexx)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine add_this_multipass (lat, m_slaves, lord_in)

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: slave, lord, slave2, lord2
type (ele_struct), optional :: lord_in
type (lat_ele_loc_struct) m_slaves(:)

integer i, j, k, n, i1, ix, ixc, ixic, ix_lord
integer n_multipass, ic, ix_l1, ix_l0

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
    call parser_warning ('INTERNAL ERROR: CONFUSED MULTIPASS SETUP.', &
                  'PLEASE GET EXPERT HELP!')
    call err_exit
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
    slave2 => pointer_to_slave(lat, slave, j)
    if (slave2%n_lord == 1) then
      i1 = i1 + 1
      write (slave2%name, '(2a, i0, a, i0)') trim(lord%name), '#', i1, '\', i      ! '
    else
      slave2%name = ''
      do k = 1, slave2%n_lord
        lord2 => pointer_to_lord(lat, slave2, k)
        lord2 => pointer_to_lord(lat, lord2, 1)
        slave2%name = trim(slave2%name) // trim(lord2%name) // '\'     ! '
      enddo
      write (slave2%name, '(a, i0)') trim(slave2%name), i  
    endif
  enddo
enddo

call control_bookkeeper (lat, lord)

end subroutine

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


end subroutine

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

logical have_inserted, found

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
  call add_superimpose (lat, super_ele, 0)
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
          slave => pointer_to_slave (lat, ref_ele, i)
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
          call add_superimpose (lat, super_ele, ix_branch, super_ele_out)
          super_ele_out%name = 'temp_name!'
        enddo

        ! Remove any multipass_lord drifts that no longer do anything.
        ! We can recognize these elements since they control super_slaves.
        ! Also mark any of their lords for deletion.

        do i = lat%n_ele_track+1, lat%n_ele_max 
          ele => lat%ele(i)
          if (ele%key /= drift$) cycle
          slave => pointer_to_slave (lat, ele, 1)
          if (slave%slave_status /= super_slave$) cycle
          ele%key = -1 ! mark for deletion
          do j = 1, ele%n_lord
            lord => pointer_to_lord (lat, ele, j)
            lord%key = -1  ! Mark lord for deletion
          enddo
        enddo
        call remove_eles_from_lat (lat, .false.) ! and delete

        ! Add a multipass_lord to control the created super_lords.
        ! If the super_lords have a single super_slave and the super_slave
        ! has only a single super_lord, the super_lords
        ! can be eliminated and the created multipass_lord can control the
        ! super_slaves directly

        j = 0

        do i = 1, branch%n_ele_track
          if (lat%ele(i)%name == 'temp_name!') then
            lat%ele(i)%name = super_ele_saved%name
            j = j + 1
            m_slaves(j) = ele_to_lat_loc (lat%ele(i))
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
          slave => pointer_to_slave (lat, ele, 1)
          if (slave%n_lord == 1) then
            do i = 1, size(m_slaves)
              ele => pointer_to_ele (lat, m_slaves(i))
              ele%key = -1 ! Mark for deletion
              ele => pointer_to_slave (lat, ele, 1)
              ele%name = super_ele_saved%name
              m_slaves(i) = ele_to_lat_loc (ele)
            enddo
            call remove_eles_from_lat (lat, .false.)
          endif
        endif

        call add_this_multipass (lat, m_slaves, super_ele_saved) 

        ! Reconnect drifts that were part of the multipass region.

        call multipass_all_info (lat, m_info) ! Save multipass info for later.

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
        call add_superimpose (lat, super_ele, i_br, super_ele_out)
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
  if (.not. found) call parser_warning ('NO MATCH FOR REFERENCE ELEMENT: ' //  &
          pele%ref_name, 'FOR SUPERPOSITION OF: ' // super_ele_saved%name, pele = pele)
endif

! if there is to be no common lord then we are done

if (.not. pele%common_lord) return

! here for common_lord, not scalled multipoles

if (super_ele_saved%key /= multipole$ .and. super_ele_saved%key /= ab_multipole$) then
  call parser_warning ( &
          'ELEMENT ' // super_ele_saved%name, &
          'IS USED WITH THE "COMMON_LORD" ATTRIBUTE BUT', &
          'THIS ELEMENT IS NOT A MULTIPOLE OR AB_MULTIPOLE', pele = pele)
  return
endif

lat%n_ele_max = lat%n_ele_max + 1
if (lat%n_ele_max > ubound(lat%ele, 1)) call allocate_lat_ele_array(lat)

nn = lat%n_ele_max 

n_con = lat%n_control_max 
lat%n_control_max = n_con + n_inserted
if (lat%n_control_max > size(lat%control)) call reallocate_control(lat, lat%n_control_max+100)


lat%ele(nn) = super_ele_saved
lat%ele(nn)%lord_status = super_lord$
lat%ele(nn)%n_slave = n_inserted
lat%ele(nn)%ix1_slave = n_con + 1
lat%ele(nn)%ix2_slave = n_con + n_inserted

do i = n_con + 1, n_con + n_inserted
  lat%control(i)%ix_lord = nn
  lat%control(i)%ix_attrib = 0
enddo

j = 0
do i_br = 0, ubound(lat%branch, 1)
  branch => lat%branch(i_br)

  do i = 1, branch%n_ele_max-1
    if (any (matched_name(1:n_inserted) == branch%ele(i)%name)) then
      it = branch%ele(i)%slave_status
      if (it /= free$) then
        call parser_warning ('SLAVE: ' // branch%ele(i)%name, &
                      'OF LORD: ' // super_ele_saved%name, &
                      'IS NOT A "FREE" ELEMENT BUT IS: ' // control_name(it), pele = pele)
        return
      endif
      j = j + 1
      branch%ele(i)%slave_status = super_slave$
      nic = lat%n_ic_max + 1
      branch%ele(i)%n_lord = 1
      branch%ele(i)%ic1_lord = nic
      branch%ele(i)%ic2_lord = nic
      lat%ic(nic) = n_con + j
      lat%control(n_con+j)%ix_slave  = i
      lat%control(n_con+j)%ix_branch = i_br
      lat%n_ic_max = nic
    endif
  enddo

enddo

if (j /= n_inserted) then
  call parser_warning ('INTERNAL ERROR! SLAVE NUMBER MISMATCH!')
  call err_exit
endif


end subroutine

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
elseif (pele%ele_pt == center$) then
  super_ele%s = super_ele%s + super_ele%value(l$) / 2
elseif (pele%ele_pt /= end$) then
  print *, 'ERROR IN COMPUTE_SUPER_LORD_S: CONTROL #1 INTERNAL ERROR!'
  call err_exit
endif

! Find the refernce point in the lattice.

ct = ref_ele%lord_status
if (ct == overlay_lord$ .or. ct == girder_lord$) then
  s_ref_begin = 1e10
  s_ref_end = 0
  do i = 1, ref_ele%n_slave
    slave => pointer_to_slave (lat, ref_ele, i)
    s_ref_begin = min(s_ref_begin, slave%s - slave%value(l$))
    s_ref_end = max(s_ref_end, slave%s)
  enddo
elseif (ct == group_lord$) then
  call parser_warning ('SUPERPOSING: ' // super_ele%name, 'UPON GROUP' // pele%ref_name)
  return
else
  s_ref_begin = ref_ele%s - ref_ele%value(l$)
  s_ref_end = ref_ele%s
endif

! Now compute the s position at the end of the element and put it in ele%s.

if (pele%ref_pt == begin$) then
  super_ele%s = super_ele%s + s_ref_begin
elseif (pele%ref_pt == center$) then
  super_ele%s = super_ele%s + (s_ref_begin + s_ref_end) / 2
elseif (pele%ref_pt == end$) then
  super_ele%s = super_ele%s + s_ref_end
else
  print *, 'ERROR IN COMPUTE_SUPER_LORD_S: CONTROL #2 INTERNAL ERROR!'
  call err_exit
endif

! For circular lattices a superimpose can wrap around the beginning or 
! the end of the lattice.

branch => lat%branch(ref_ele%ix_branch)
if (branch%param%lattice_type == circular_lattice$) then
  if (super_ele%s > branch%ele(branch%n_ele_track)%s) then
    super_ele%s = super_ele%s - lat%param%total_length
  elseif (super_ele%s < 0) then
    super_ele%s = super_ele%s + branch%param%total_length
  endif
endif

end subroutine

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
    call parser_warning ('BAD ARGUMENT LIST FOR: ', seq_name)
    return
  endif
  n_arg = n_arg + 1
  name(n_arg) = word
  if (delim == ')') exit
enddo

err_flag = .false.
allocate (arg_list(n_arg))
arg_list = name(1:n_arg)

end subroutine

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

if (delim /= '(') call parser_warning  &
      ('EXPECTING "(", GOT: ' // delim, 'FOR LINE: ' // seq%name)
if (ix_word /= 0)  call parser_warning  &
      ('EXTRANEOUS STUFF BEFORE "(", GOT: ' // word,  &
      'FOR LINE: ' // seq%name)

! now parse list proper

ix_ele = 1
this_ele => s_ele(ix_ele)

do 

  call get_next_word (word, ix_word, ':=(,)[]', delim, delim_found, .true.)

  ix = index(word, '*')          ! E.g. word = '-3*LINE'
  if (ix /= 0) then
    bp_com%parse_line = word(:ix-1) // "," // bp_com%parse_line
    call evaluate_value (trim(seq%name) // ' Repetition Count', rcount, &
                            lat, c_delim, c_delim_found, err_flag)
    this_ele%rep_count = nint(rcount)
    if (err_flag) return
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
      call parser_warning ('NO MATCHING "]" FOUND FOR OPENING "[" IN SEQUENCE: ' // seq%name)
      return
    endif
    call get_next_word (word, ix_word, '[]:=(,)', delim, delim_found, .true.)
    if (ix_word > 0) then
      call parser_warning ('ILLEGAL CHARACTERS AFTER CLOSING "]" FOUND IN SEQUENCE: ' // seq%name)
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
    if (word /= ' ') call parser_warning &
              ('NO COMMA AFTER SUBLINE OR REPLACEMENT LINE. FOUND: ' // &
               word, 'IN THE SEQUENCE: ' // seq%name)
  endif

  if (this_ele%name == ' ') call parser_warning &
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
    call parser_warning ('EXPECTING "," GOT: ' // delim, 'FOR LINE: ' // seq%name)
    exit
  endif
         
enddo

! make sure there is nothing else if at top level

if (top_level) then
  call get_next_word(word, ix_word, ':=() ', delim, delim_found, .true.)
  if (delim_found .or. ix_word /= 0) call parser_warning  &
        ('EXTRA CHARACTERS AFTER CLOSING ")"',  'FOR LINE: ' // seq%name)
endif

! transfer

ix_ele = ix_ele - 1
allocate (seq%ele(ix_ele))

do i = 1, ix_ele
  seq%ele(i) = s_ele(i)
enddo

deallocate (s_ele)

end subroutine

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
  plat%ele(i)%ref_pt  = center$
  plat%ele(i)%ele_pt  = center$
  plat%ele(i)%s       = 0
  plat%ele(i)%common_lord = .false.
enddo

end subroutine

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
type (ele_struct), pointer :: lord, lord2, slave
type (parser_lat_struct) plat
type (control_struct), pointer, save :: cs(:) => null()
type (branch_struct), pointer :: branch

integer ixx, i, ic, n, n2, k, k2, ix, j, ie, ix1, ix2, ns, ixs
integer ix_lord, k_slave, ix_ele_original
integer, allocatable :: r_indexx(:), ix_ele(:), ix_branch(:)

character(40), allocatable :: name_list(:)
character(40) name, name1, slave_name, attrib_name, missing_slave_name

logical err, slave_not_in_lat

! Setup...
! in_lat has the lords that are to be added to lat.
! We add an extra 1000 places to the arrays to give us some overhead.

n = n2 + 1000
do i = 0, ubound(lat%branch, 1)
    n = n + lat%branch(i)%n_ele_max
enddo

allocate (name_list(n), ix_ele(n), ix_branch(n), r_indexx(n))
allocate (cs(1000))

ix1 = 0
do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  n = branch%n_ele_max
  ix2 = ix1 + n
  name_list(ix1+1:ix2) = branch%ele(1:n)%name
  ix_ele(ix1+1:ix2)    = branch%ele(1:n)%ix_ele
  ix_branch(ix1+1:ix2) = branch%ele(1:n)%ix_branch
  ix1 = ix2
enddo

call indexx (name_list(1:ix1), r_indexx(1:ix1)) ! get sorted list

! loop over elements

main_loop: do n = 1, n2

  lord => in_lat%ele(n)  ! next lord to add

  !-----------------------------------------------------
  ! overlay and groups

  select case (lord%lord_status)
  case (overlay_lord$, group_lord$)
 
    call new_control (lat, ix_lord)  ! get index in lat where lord goes
    lat%ele(ix_lord) = lord
    ixx = lord%ixx

    ! Find where the slave elements are. 
    ! If a slave element is not in lat but is in in_lat then the slave has 
    ! not yet been used in the lattice list. 
    ! In this case do not add the lord to the lattice.

    j = 0 ! number of slaves found
    slave_not_in_lat = .false.  ! Is there a slave that is not in the lattice?

    do i = 1, lord%n_slave

      name = plat%ele(ixx)%name(i)
      call find_indexx (name, name_list, r_indexx, ix1, k, k2)

      if (k == 0) then
        slave_not_in_lat = .true.
        missing_slave_name = name
      endif

      if ((k == 0 .and. j > 0) .or. (k > 0 .and. slave_not_in_lat) .or. &
          (k == 0 .and. all(in_lat%ele(1:n2)%name /= name))) then
        call parser_warning ('CANNOT FIND SLAVE FOR: ' // lord%name, &
                      'CANNOT FIND: '// missing_slave_name, pele = plat%ele(ixx))
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
        cs(j)%coef = plat%ele(ixx)%coef(i)
        cs(j)%ix_slave = ix_ele(k)
        cs(j)%ix_branch = ix_branch(k)
        cs(j)%ix_lord = -1             ! dummy value
        attrib_name = plat%ele(ixx)%attrib_name(i)
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
          call parser_warning ('IN OVERLAY OR GROUP ELEMENT: ' // lord%name, &
                        'ATTRIBUTE: ' // attrib_name, &
                        'IS NOT A VALID ATTRIBUTE OF: ' // slave%name, &
                        pele = plat%ele(ixx))
          cycle main_loop
        endif
        k2 = k2 + 1
        if (k2 > ix1) exit
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

    call find_indexx (lord%name, lat%ele(1:ix1)%name, r_indexx(1:ix1), ix1, k, k2)
    ix1 = ix1 + 1
    r_indexx(k2+1:ix1) = r_indexx(k2:ix1-1)
    r_indexx(k2) = ix1
    name_list(ix1) = lord%name
    ix_ele(ix1) = ix_lord
    ix_branch(ix1) = 0

    ! create the lord

    ns = lord%n_slave

    select case (lord%lord_status)
    case (overlay_lord$)
      call create_overlay (lat, ix_lord, lord%component_name, cs(1:ns), err)
    case (group_lord$)
      call create_group (lat, ix_lord, cs(1:ns), err)
    end select
    if (err) call parser_warning ('ELEMENT OR GROUP: ' // lord%name, &
                           'IS TRYING TO CONTROL AN ATTRIBUTE THAT IS NOT FREE TO VARY!', &
                           pele = plat%ele(ixx))

  !-----------------------------------------------------
  ! girder
  ! Create an girder element for each element whose name matches the
  ! first name in the slave list.

  case (girder_lord$)

    ixx = lord%ixx
    name1 = plat%ele(ixx)%name(1)

    call find_indexx (name1, name_list, r_indexx, ix1, k_slave, k2)
    if (k_slave == 0) then
      call parser_warning ('CANNOT FIND START ELEMENT FOR GIRDER: ' // lord%name, &
                    'CANNOT FIND: '// name, pele = plat%ele(ixx))
      cycle
    endif

    slave => pointer_to_ele (lat, ix_ele(k_slave), ix_branch(k_slave))

    if (slave%lord_status == super_lord$) then 
      slave => pointer_to_slave (lat, slave, 1)
    endif

    branch => lat%branch(slave%ix_branch)
    ix_ele_original = slave%ix_ele

    ! Loop over all matches to the first name.

    do 

      ixs = 0       ! Index of slave element we are looking for

      slave_loop: do            ! loop over all slaves
        ixs = ixs + 1
        if (ixs > lord%n_slave) exit
        slave_name = plat%ele(ixx)%name(ixs)

        do  ! loop over all lattice elements
          if (slave%slave_status == super_slave$) then
            do ic = 1, slave%n_lord
              lord2 => pointer_to_lord (lat, slave, ic)
              if (match_wild(lord2%name, slave_name)) then
                cs(ixs)%ix_slave  = lord2%ix_ele
                cs(ixs)%ix_branch = lord2%ix_branch
                cycle slave_loop
              endif
            enddo
          else
            if (match_wild(slave%name, slave_name)) then
              cs(ixs)%ix_slave  = slave%ix_ele
              cs(ixs)%ix_branch = slave%ix_branch
              cycle slave_loop
            endif
          endif
          if (slave%ix_ele == branch%n_ele_track) then
            slave => branch%ele(1)
          else
            slave => branch%ele(slave%ix_ele + 1)
          endif
          if (slave%ix_ele == ix_ele_original) then
            call parser_warning ('CANNOT FIND END ELEMENT FOR GIRDER: ' // lord%name, &
                          'CANNOT FIND: ' // slave_name, pele = plat%ele(ixx))
            cycle main_loop
          endif
        enddo 
      enddo slave_loop

      ! create the girder element

      call new_control (lat, ix_lord)
      call create_girder (lat, ix_lord, cs(1:lord%n_slave), lord)

      k2 = k2 + 1
      k_slave = r_indexx(k2)
      if (slave%name /= name1) exit
      if (k_slave > lat%n_ele_track) then ! must be a super_lord.
        ix = slave%ix1_slave
        k_slave = lat%control(ix)%ix_slave
      endif

    enddo 

  end select

enddo main_loop

! cleanup

deallocate (r_indexx)
deallocate (name_list)
deallocate (cs)

end subroutine

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

case (wiggler$) 
  ele%sub_key = periodic_type$   ! default
  ele%value(polarity$) = 1.0     ! default

case (custom$)  ! set defaults
  ele%mat6_calc_method = custom$
  ele%tracking_method  = custom$
  ele%field_calc       = custom$

case (bend_sol_quad$) ! set defaults
  ele%mat6_calc_method = symp_lie_bmad$
  ele%tracking_method  = symp_lie_bmad$

case (taylor$)   ! start with unit matrix
  ele%tracking_method = taylor$  ! default
  ele%mat6_calc_method = taylor$ ! default
  ele%taylor_order = 999         ! make large.
  call taylor_make_unit (ele%taylor)

case (rbend$, sbend$)
  ele%value(fintx$) = real_garbage$
  ele%value(hgapx$) = real_garbage$

case (branch$, photon_branch$)
  ele%value(direction$) = 1

case (rcollimator$)
  ele%offset_moves_aperture = .true.

case (ecollimator$)
  ele%offset_moves_aperture = .true.
  ele%aperture_type = elliptical$

case (lcavity$)
  ele%value(coupler_at$) = exit_end$

end select

end subroutine

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

  if (ele%value(b_field$) /= 0 .and. ele%key == rbend$) call parser_warning &
          ("B_FIELD NOT SETTABLE FOR AN RBEND (USE AN SBEND INSTEAD): " // ele%name)

  if (ele%value(b_field$) /= 0 .and. ele%value(g$) /= 0) call parser_warning &
          ('BOTH G AND B_FIELD SET FOR A BEND: ' // ele%name)

  if (ele%value(b_field$) /= 0 .and. ele%value(rho$) /= 0) call parser_warning &
          ('BOTH RHO AND B_FIELD SET FOR A BEND: ' // ele%name)

  if (ele%value(g$) /= 0 .and. ele%value(rho$) /= 0) &
            call parser_warning ('BOTH G AND RHO SPECIFIED FOR BEND: ' // ele%name)

  ! if rho is set then this gives g

  if (ele%value(rho$) /= 0) ele%value(g$) = 1 / ele%value(rho$)

  ! If g and angle are set then this determines l

  length_set = .false.
  if (ele%value(g$) /= 0 .and. angle /= 0) then
    if (ele%value(l$) /= 0) call parser_warning ('ANGLE, G/RHO, AND L SPECIFIED FOR BEND: ' // ele%name)
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
    call parser_warning ('THE BENDING ANGLE IS NONZERO IN A ZERO LENGTH BEND! ' // ele%name)
  elseif (ele%value(angle$) /= 0) then
    ele%value(g$) = ele%value(angle$) / ele%value(l$) 
  endif

  ! If fintx or hgapx are real_garbage then they have not been set.
  ! If so, set their valuse to fint and hgap.

  if (ele%value(hgapx$) == real_garbage$) ele%value(hgapx$) = ele%value(hgap$)
  if (ele%value(fintx$) == real_garbage$) ele%value(fintx$) = ele%value(fint$)

! Accept Use of Delta_E for lcavities and vary the mode frequencies.

case (lcavity$) 

  if (ele%value(delta_e$) /= 0) then
    if (ele%value(gradient$) /= 0) call parser_warning &
                ('BOTH DELTA_E AND GRADIENT NON-ZERO FOR A LCAVITY:', ele%name)
    ele%value(gradient$) = ele%value(delta_e$) / ele%value(l$)
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
  if (ele%field_master .and. (ele%value(ks$) /= 0 .or. kick_set)) call parser_warning &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (KS, HKICK, ETC.) AND FIELD SET FOR A SOLENOID.')

case (sol_quad$)
  if (ele%field_master .and. (ele%value(ks$) /= 0 .or. &
                            ele%value(k1$) /= 0 .or. kick_set)) call parser_warning &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K1, HKICK, ETC.) AND FIELD SET FOR A SOL_QUAD.')

case (quadrupole$)
  if (ele%field_master .and. (ele%value(k1$) /= 0 .or. kick_set)) call parser_warning &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K1, HKICK, ETC.) AND FIELD SET FOR A QUAD.')

case (sextupole$)
  if (ele%field_master .and. (ele%value(k2$) /= 0 .or. kick_set)) call parser_warning &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K2, HKICK, ETC.) AND FIELD SET FOR A SEXTUPOLE.')

case (octupole$)
  if (ele%field_master .and. (ele%value(k3$) /= 0 .or. kick_set)) call parser_warning &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K3, HKICK, ETC.) AND FIELD SET FOR A OCTUPOLE.')

case (hkicker$, vkicker$)
  if (ele%field_master .and. (ele%value(kick$) /= 0 .or. kick_set)) call parser_warning &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH AND BL_KICK SET FOR A H/VKICKER.')

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
    if (n == 0) call parser_warning ( &
          'FOR MULTIPASS LORD: ' // ele%name, &
          'N_REF_PASS, E_TOT, AND P0C ARE ALL ZERO AND FIELD_MASTER = FALSE!')
    if (n > 1) call parser_warning ( &
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
! The standard digested file name has 'digested_' (single precision Bmad) or 
! 'digested8_' (standard double precision Bmad) prepended to the file name.
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
if (rp == 8) then
  digested_file = trim(path) // 'digested8_' // basename
else
  digested_file = trim(path) // 'digested_' // basename
endif

! This only affects VMS programs.
! What we want to do is change the directory for lattice files in CESR_MNT:[lattice...].
! However 'CESR_MNT' is a logical that will get translated by the inquire function so
! we only check that 'lattice' is the top directory.

ix = max (index_nocase(digested_file, '[lattice.'), &
          index_nocase(digested_file, '[000000.lattice.'))
if (ix /= 0) then
  ix = index_nocase(digested_file, 'lattice.')
  digested_file = 'U:[cesr.lattice.' // digested_file(ix+8:)
endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine save_taylor_elements (lat, ele_array)
!
! Subroutine to save the taylor maps in a lattice.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine save_taylor_elements (lat, ele_array)

implicit none

type (lat_struct) lat
type (ele_struct), allocatable :: ele_array(:) 

integer i, ix

!

if (.not. associated(lat%ele)) return

ix = 0
do i = 1, lat%n_ele_max
  if (associated(lat%ele(i)%taylor(1)%term)) ix = ix + 1
enddo

if (ix /= 0) then
  if (allocated(ele_array)) deallocate(ele_array)
  allocate(ele_array(ix))
  ix = 0
  do i = 1, lat%n_ele_max
    if (associated(lat%ele(i)%taylor(1)%term)) then
      ix = ix + 1
      ele_array(ix) = lat%ele(i)
    endif
  enddo
endif  

end subroutine save_taylor_elements

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine reuse_taylor_elements (lat, ele_array)
!
! Subroutine to reuse saved taylor maps in a lattice.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine reuse_taylor_elements (lat, ele_array)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_struct), allocatable :: ele_array(:) 

integer i, j

! Reuse the old taylor series if they exist
! and the old taylor series has the same attributes.

if (.not. allocated(ele_array)) return

do i = 1, lat%n_ele_max

  ele => lat%ele(i)
  call attribute_bookkeeper (ele, lat%param) ! for equivalent_taylor_attributes test

  do j = 1, size(ele_array)
    if (any(ele_array(j)%taylor(:)%ref /= 0)) cycle
    if (bmad_com%taylor_order > ele_array(j)%taylor_order) cycle
    if (.not. equivalent_taylor_attributes (ele_array(j), ele)) cycle
    exit
  enddo

  if (j == size(ele_array) + 1) cycle
  call out_io (s_info$, bp_com%parser_name, 'Reusing Taylor for: ' // ele_array(j)%name)
  call transfer_ele_taylor (ele_array(j), ele, bmad_com%taylor_order)
enddo

deallocate(ele_array)

end subroutine

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
type (ele_struct), pointer :: ele2
type (seq_struct), target :: sequence(:)
type (branch_struct), pointer :: branch

integer, allocatable :: seq_indexx(:), in_indexx(:)
integer j, nb, n_ele_use

character(*), allocatable ::  in_name(:), seq_name(:)

!

nb = ubound(lat%branch, 1) + 1
call allocate_branch_array (lat, nb)
branch_ele%value(ix_branch_to$) = nb
branch => lat%branch(nb)
if (branch_ele%key == branch$) then
  branch%param%particle = lat%branch(branch_ele%ix_branch)%param%particle
else
  branch%param%particle = photon$
endif
branch%param%lattice_type = linear_lattice$
branch%key            = branch_ele%key
branch%ix_branch      = nb
branch%ix_from_branch = branch_ele%ix_branch
branch%ix_from_ele    = branch_ele%ix_ele
branch%name           = branch_ele%name
!! if (branch_ele%alias /= '') branch%name = branch_ele%alias
call parser_expand_line (nb, lat, branch_ele%component_name, sequence, in_name, &
                              in_indexx, seq_name, seq_indexx, in_lat, n_ele_use)
branch%ele(0)%key     = init_ele$
branch%ele(0)%name    = 'BEGINNING'
branch%n_ele_track = n_ele_use
branch%n_ele_max   = n_ele_use

do j = 1, n_ele_use
  ele2 => lat%branch(nb)%ele(j)
  if (ele2%key == photon_branch$ .or. ele2%key == branch$) then
    call parser_add_branch (ele2, lat, sequence, in_name, in_indexx, &
                                                    seq_name, seq_indexx, in_lat)
  endif
enddo

end subroutine

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
n_max = size(in_name)

call find_indexx (use_name, seq_name, seq_indexx, iseq_tot, i_use)
if (i_use == 0) then
  call parser_warning ('CANNOT FIND DEFINITION OF LINE IN "USE" STATEMENT: ' // use_name, ' ')
  return
endif

if (sequence(i_use)%type /= line$) then
  call parser_warning ('NAME IN "USE" STATEMENT IS NOT A LINE!', ' ')
  return
endif

! Now to expand the lines and lists to find the elements to use.
! First go through the lines and lists and index everything.

do k = 1, iseq_tot
  do i = 1, size(sequence(k)%ele(:))

    s_ele => sequence(k)%ele(i)
    name = s_ele%name

!      ix = index(name, '\')   ! ' 
!      if (ix /= 0) name = name(:ix-1) ! strip off everything after \

    if (s_ele%ix_arg > 0) then   ! dummy arg
      s_ele%type = element$
      cycle
    endif

    call find_indexx (name, in_name, in_indexx, n_max, j)
    if (j == 0) then  ! if not an element it must be a sequence
      call find_indexx (name, seq_name, seq_indexx, iseq_tot, j)
      if (j == 0) then  ! if not a sequence then I don't know what it is
        s_ele%ix_ele = -1
        s_ele%type = element$
      else
        s_ele%ix_ele = j
        s_ele%type = sequence(j)%type
      endif
      if (s_ele%type == list$ .and. s_ele%reflect) call parser_warning ( &
                          'A REFLECTION WITH A LIST IS NOT ALLOWED IN: '  &
                          // sequence(k)%name, 'FOR LIST: ' // s_ele%name, &
                          seq = sequence(k))
      if (sequence(k)%type == list$) &
                call parser_warning ('A REPLACEMENT LIST: ' // sequence(k)%name, &
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
  call parser_warning ('"USE"D LINE FOR LATTICE EXPANSION IS MARKED MULTIPASS!')
  call err_exit
endif

! Expand "used" line...

line_expansion: do

  ! if rep_count is zero then change %ix_ele index by +/- 1 and reset the rep_count.
  ! if we have got to the end of the current line then pop the stack back to
  ! the next lower level.
  ! Also check if we have gotten to level 0 which says that we are done.
  ! If we have stepped out of a multipass line which has been transversed in reverse
  !   then we need to do some bookkeeping to keep the elements straight.

  if (stack(i_lev)%rep_count == 0) then      ! goto next element in the sequence
    stack(i_lev)%ix_ele = stack(i_lev)%ix_ele + stack(i_lev)%direction 
    ix = stack(i_lev)%ix_ele

    if (ix > 0 .and. ix <= size(seq%ele)) then
      stack(i_lev)%rep_count = seq%ele(ix)%rep_count
    else
      i_lev = i_lev - 1
      if (i_lev == 0) exit line_expansion
      seq => sequence(stack(i_lev)%ix_seq)
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
    call find_indexx (name, in_name, in_indexx, n_max, j)
    if (j == 0) then  ! if not an element it must be a sequence
      call find_indexx (name, seq_name, seq_indexx, iseq_tot, j)
      if (j == 0) then  ! if not a sequence then I don't know what it is
        call parser_warning ('CANNOT FIND DEFINITION FOR: ' // name, &
                          'IN LINE: ' // seq%name, seq = seq)
        call err_exit
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
        call parser_warning ('ELEMENTS IN A LINE OR LIST ARE NOT ALLOWED TO HAVE A TAG.', &
                      'FOUND ILLEGAL TAG FOR ELEMENT: ' // s_ele%name, &
                      'IN THE LINE/LIST: ' // seq%name, seq)
      endif
      this_seq_ele => s_ele
    endif

    if (this_seq_ele%ix_ele < 1) call parser_warning('NOT A DEFINED ELEMENT: ' // &
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
      call parser_warning ('NESTED LINES EXCEED STACK DEPTH!')
      call err_exit
    endif
    if (s_ele%type == replacement_line$) then
      seq2 => sequence(s_ele%ix_ele)
      if (size(seq2%dummy_arg) /= size(s_ele%actual_arg)) then
        call parser_warning ('WRONG NUMBER OF ARGUMENTS FORREPLACEMENT LINE: ' // &
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
    call parser_warning ('INTERNAL SEQUENCE ERROR!')

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
    call type_ele (lat%ele(i), .false., 0, .false., 0, .true., lat)
  enddo
endif

if (index(debug_line, 'LORD') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'LORD elements: ', lat%n_ele_max - lat%n_ele_track
  do i = lat%n_ele_track+1, lat%n_ele_max
    print *, '-------------'
    print *, 'Ele #', i
    call type_ele (lat%ele(i), .false., 0, .false., 0, .true., lat)
  enddo
endif

if (index(debug_line, 'LATTICE') /= 0) then  
  print *
  print *, '----------------------------------------'
  print *, 'Lattice Used: ', lat%name
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
    call type_ele (lat%ele(i), .false., 0, .true., 0, .true., lat)
    call string_trim (debug_line(ix+1:), debug_line, ix)
  enddo
endif

if (index(debug_line, 'BEAM_START') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'beam_start:'
  print '(3x, 6es13.4)', lat%beam_start%vec      
endif

end subroutine

end module


