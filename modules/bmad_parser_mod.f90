!+
! Module bmad_parser_mod
!
! This module is a collection of helper routines used by bmad_parser and bmad_parser2.
! The routines in this module are specifically taylored for bmad_parser and
! bmad_parser2 and cannot, in general, be used otherwise.
!-

#include "CESR_platform.inc"

module bmad_parser_mod

  use ptc_interface_mod
  use bookkeeper_mod

! A "sequence" is a line or a list.
! The information about a sequence is stored in a seq_struct.

! A seq_struct has an array of seq_ele_struct structures.
! Each seq_ele_struct represents an individual element in a sequence and, 
! since sequences can be nested, can itself be a line or a list.

  type seq_ele_struct
    character(16) name             ! name of element, subline, or sublist
    character(16), pointer :: actual_arg(:) => null()
    integer type                   ! LINE$, REPLACEMENT_LINE$, LIST$, ELEMENT$
    integer ix_array               ! if an element: pointer to ELE_ array
                                   ! if a list: pointer to SEQ_ array
    integer ix_arg                 ! index in arg list (for replacement lines)
    integer rep_count              ! how many copies of an element
    logical reflect                ! reflection of subline?
  end type

  type seq_struct
    character(16) name                 ! name of sequence
    type (seq_ele_struct), pointer :: ele(:) => null()
    character(16), pointer :: dummy_arg(:) => null()
    character(16), pointer :: corresponding_actual_arg(:) => null()
    integer type                       ! LINE$, REPLACEMENT_LINE$ or LIST$
    integer ix                         ! current index of element in %ELE
    integer indexx                     ! alphabetical order sorted index
    character(200) file_name     ! file where sequence is defined
    integer ix_line              ! line number in filewhere sequence is defined
  end type


! A LIFO stack structure is used in the final evaluation of the line that is
! used to form a lattice

  type seq_stack_struct
    integer ix_seq                ! index to seq_(:) array
    integer ix_ele                ! index to seq%ele(:) array
    integer rep_count             ! repetition count
    integer reflect               ! reflection sequence?
  end type

! A LIFO stack structure is used to hold the list of input lattice files
! that are currently open.

  type stack_file_struct
    character(200) logical_name
    character(200) full_name
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
    character(16) ref_name
    character(16), pointer :: name_(:) => null()
    character(16), pointer :: attrib_name_(:) => null()
    real(rp), pointer :: coef_(:) => null()
    integer ix_count
    integer ele_pt, ref_pt
    logical common_lord
    real(rp) s
    integer indexx
  end type

  type parser_ring_struct
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

  type bp_com_struct
    type (stack_file_struct) current_file
    type (stack_file_struct) calling_file
    character(16), pointer :: var_name(:) => null()    ! variable name
    real(rp), pointer :: var_value(:) => null()        ! variable value
    integer, pointer :: var_indexx(:) => null()        ! variable sort index
    integer n_files
    character(200) file_name_(50)        ! List of files all opened.
    character(280) parse_line
    character(140) input_line1          ! For debug messages
    character(140) input_line2          ! For debug messages
    character(16) parser_name
    character(72) debug_line
    logical parser_debug, write_digested, error_flag
    integer ivar_tot, ivar_init
    logical input_line_meaningful
    character(200) :: dirs(3) = (/ &
                      './           ', './           ', '$BMAD_LAYOUT:' /)
  end type

!

  type (bp_com_struct), save :: bp_com
  type (ele_struct), target, save :: beam_ele, param_ele

  character(16) :: blank = ' '

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_attribute (how, ele, ring, pring,
!                                             delim, delim_found, err_flag)
!
! Subroutine used by bmad_parser and bmad_parser2 to get the value of
! an attribute from the input file.
!
! This subroutine is not intended for general use.
!-


subroutine get_attribute (how, ele, ring, pring, &
                                             delim, delim_found, err_flag)

  use random_mod
         
  implicit none

  type (ring_struct)  ring
  type (parser_ring_struct) pring
  type (ele_struct), target ::  ele, ele0
  type (wig_term_struct), pointer :: wig_term(:)

  real(rp) kx, ky, kz, tol, value, coef
  real(rp), pointer :: r_ptr

  integer i, ic, ix_word, how, ix_word1, ix_word2, ix_word3, ios, ix, i_out
  integer expn(6), ix_attrib

  character(16) :: word, tilt_word = 'TILT', str_ix
  character delim*1, delim1*1, delim2*1, str*80, err_str*40, line*80
  character(16) :: super_names(11) = (/ &
                'SUPERIMPOSE  ', 'OFFSET       ', 'REFERENCE    ',  &
                'ELE_BEGINNING', 'ELE_CENTER   ', 'ELE_END      ',  &
                'REF_BEGINNING', 'REF_CENTER   ', 'REF_END      ', &
                'COMMON_LORD  ', '             ' /)

  logical delim_found, err_flag

! taylor

  if (ele%key == taylor$) then

    call get_next_word (str, ix_word, '}', delim, delim_found, .true.)
    str = trim(str) // delim
    bp_com%parse_line = str // bp_com%parse_line

    call get_next_word (word, ix_word, ',{}', delim, delim_found, .true.)
    if (delim /= '{' .or. ix_word /= 0) then
      call warning ('BAD TERM FOR TAYLOR ELEMENT: ' // ele%name, &
                                              'CANNOT PARSE: ' // str)
      return
    endif

    call get_next_word (word, ix_word, ',:}', delim, delim_found, .true.)
    read (word, *, iostat = ios) i_out
    if (delim /= ':' .or. ix_word == 0 .or. ios /= 0) then
      call warning ('BAD "OUT" IN TERM FOR TAYLOR ELEMENT: ' // ele%name, &
                                                      'CANNOT PARSE: ' // str)
      return
    endif

    call evaluate_value (str, coef, ring, delim, delim_found, err_flag)
    if (err_flag) return

    call get_next_word (line, ix_word, '},', delim, delim_found, .true.)
    read (line, *, iostat = ios) expn
    if (delim /= '}' .or. ix_word == 0 .or. ios /= 0) then
      call warning ('BAD "EXPONENT" IN TERM FOR TAYLOR ELEMENT: ' // &
                                            ele%name, 'CANNOT PARSE: ' // str)
      return
    endif

    call add_taylor_term (ele, i_out, coef, expn)
    call get_next_word (word, ix_word, '},', delim, delim_found, .true.)

    if (ix_word /= 0 .or. (delim_found .and. delim /= ',')) then
      call warning ('BAD TERM ENDING FOR TAYLOR ELEMENT: ' // ele%name, &
                                                      'CANNOT PARSE: ' // str)
      return
    endif

    return
  endif

! Get next WORD.
! If an overlay or group element then word is just an attribute to control
! [except for a "GROUP[COMMAND] = 0.343" redef construct]

  err_flag = .true.  ! assume the worst
  call get_next_word (word, ix_word, ':, =()', delim, delim_found, .true.)

  if (ele%key == overlay$) then
    i = attribute_index(ele, word)       ! general attribute search
    if (i < 1) then
      call warning ('BAD OVERLAY ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
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
        call evaluate_value (trim(ele%name) // ' ' // word, value, &
                                      ring, delim, delim_found, err_flag)
        if (err_flag) return
        ele%value(i) = value
      else
        ele%value(i) = 0
      endif
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
      call warning ('BAD GROUP ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
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
        call evaluate_value (trim(ele%name) // ' ' // word, value, &
                                ring, delim, delim_found, err_flag)
        if (err_flag) return
        if (how == def$) then
          ele%value(command$) = value
        else
          ele%value(i) = value
        endif
      elseif (how == redef$) then
        call warning ('NO VALUE GIVEN FOR ATTRIBUTE FOR: ' // ele%name)
      endif
    endif

    err_flag = .false.
    return

  endif

! beginning element

  if (ele%key == init_ele$) then
    call evaluate_value (trim(ele%name) // ' ' // word, value, &
                                      ring, delim, delim_found, err_flag) 
    if (err_flag) return
    call pointer_to_attribute (ele, word, .false., r_ptr, ix_attrib, err_flag)
    if (err_flag) then
      bp_com%error_flag = .true.
      return
    endif

    r_ptr = value
    
    if (ele%x%beta /= 0) ele%x%gamma = (1 + ele%x%alpha**2) / ele%x%beta
    if (ele%y%beta /= 0) ele%y%gamma = (1 + ele%y%alpha**2) / ele%y%beta
    ele%gamma_c = sqrt(1 - ele%c_mat(1,1)*ele%c_mat(2,2) + &
                                              ele%c_mat(1,2)*ele%c_mat(2,1))
    return
  endif

! if not an overlay then see if it is an ordinary attribute.
! if not an ordinary attribute then might be a superimpose switch

  i = attribute_index(ele, word)

  if (i < 1) then          ! if not an ordinary attribute...
    if (ix_word == 0) then  ! no word
      call warning  &
            ('"," NOT FOLLOWED BY ATTRIBUTE NAME FOR: ' // ele%name)
      return
    else
      if (word(:ix_word) == 'REF') word = 'REFERENCE' ! allowed abbrev
      call match_word (word, super_names, i)
      if (i < 1) then
        call warning  &
            ('BAD ATTRIBUTE NAME: ' // word, 'FOR ELEMENT: ' // ele%name)

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
          call evaluate_value (trim(ele%name) // ' ' // word, value, &
                                          ring, delim, delim_found, err_flag)
          if (err_flag) return
          pring%ele(ic)%s = value
        else
          print *, 'ERROR IN BMAD_PARSER: INTERNAL ERROR. PLEASE GET HELP!'
          call err_exit
        endif
      endif

      err_flag = .false.
      return

    endif
  endif

! wiggler term attribute

  if (i == term$) then

    err_flag = .true. ! assume the worst

    if (delim /= '(') then   ! ) then
      call warning ('"TERM" FOR A WIGGLER NOT FOLLOWED BY A "(" FOR: ' // &
                                                              ele%name)  ! )
      return
    endif

    call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.) ! (
    if (delim /= ')') then
      call warning ('CANNOT FIND CLOSING ")" for a "TERM(i)" FOR A WIGGLER"', &
                    'FOR: ' // ele%name)
      return
    endif

    read (word, *, iostat = ios) ix
    if (ix < 1 .or. ios /= 0) then
      call warning ('BAD TERM NUMBER FOR A WIGGLER FOR: ' // ele%name)
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

    if (delim1 /= '=' .or. delim2 /= '{' .or. &
                                    ix_word1 /= 0 .or. ix_word2 /= 0) then
      call warning ('CONFUSED SYNTAX FOR TERM IN WIGGLER: ' // ele%name, &
                                                                      str_ix)
      return
    endif

    err_str = trim(ele%name) // ' ' // str_ix

    call evaluate_value (err_str, ele%wig_term(ix)%coef, &
                              ring, delim, delim_found, err_flag)
    if (err_flag) return
 
    call evaluate_value (err_str, ele%wig_term(ix)%kx, &
                              ring, delim, delim_found, err_flag)
    if (err_flag) return

    call evaluate_value (err_str, ele%wig_term(ix)%ky, &
                              ring, delim, delim_found, err_flag)
    if (err_flag) return

    call evaluate_value (err_str, ele%wig_term(ix)%kz, &
                              ring, delim, delim_found, err_flag)
    if (err_flag) return

    call evaluate_value (err_str, ele%wig_term(ix)%phi_z, &
                              ring, delim, delim_found, err_flag)
    if (err_flag) return


    if (delim /= '}') then
      call warning ('ENDING "}" NOT FOUND FOR WIGGLER: ' // ele%name, str_ix)
      err_flag = .true.
      return
    endif

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
      call warning ('WIGGLER TERM DOES NOT HAVE CONSISTANT Kx, Ky, and Kz', &
                    'FOR WIGGLER: ' // ele%name // '  ' // str_ix)
      err_flag = .true.
      return
    endif

    call get_next_word (word, ix_word,  ':,=()', delim,  delim_found, .true.)  
    if (ix_word /= 0) then
      call warning ('BAD SYNTAX FOR WIGGLER: ' // ele%name, str_ix)
      err_flag = .true.
      return
    endif

    ele%sub_key = map_type$
    return

  endif

! check that next delim is a "=". If not check for a possible default value
! otherwise it is an error

  if (delim /= '=')  then
    if (word == tilt_word) then
      select case (ele%key)
      case (sbend$, rbend$)
        ele%value(tilt$) = pi / 2
      case (quadrupole$, sol_quad$) 
        ele%value(tilt$) = pi / 4
      case (sextupole$) 
        ele%value(tilt$) = pi / 6
      case (octupole$) 
        ele%value(tilt$) = pi / 8
      case default
        call warning ('SORRY I''M NOT PROGRAMMED TO USE A "TILT" DEFAULT' // &
                'FOR A: ' // key_name(ele%key), 'FOR: ' // ele%name)
      end select
    elseif (word == 'FINT') then
      ele%value(fint$) = 0.5
    elseif (word == 'FINTX') then
      ele%value(fintx$) = 0.5
    elseif (ele%key == multipole$) then
      if (i >= t0$) then
        ele%b(i-t0$) = pi / (2*(i-t0$) + 2)
      else
        call warning ('EXPECTING "=" AFTER MULTIPOLE ATTRIBUTE: ' // word,  &
                         'FOR ELEMENT: ' // ele%name)
        err_flag = .true.
      endif
    else
      call warning ('EXPECTING "=" AFTER ATTRIBUTE: ' // word,  &
                         'FOR ELEMENT: ' // ele%name)
      err_flag = .true.
    endif
    return
  endif

! get the value of the attribute.
! The TYPE, ALIAS, and DESCRIP attributes are special because their "values"
! are character strings

  select case (i)

  case(type$, alias$, descrip$, sr_wake_file$, lr_wake_file$)
    call type_get (ele, i, delim, delim_found)

  case (symplectify$) 
    if (how == def$ .and. (delim == ',' .or. .not. delim_found)) then
      ele%symplectify = .true.
    else
      call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.)
      ele%symplectify = evaluate_logical (word, ios)
      if (ios /= 0 .or. ix_word == 0) then
        call warning ('BAD "SYMPLECTIFY" SWITCH FOR: ' // ele%name)
        return
      endif
    endif
    
  case (is_on$)
    call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.)
    ele%is_on = evaluate_logical (word, ios)
    if (ios /= 0 .or. ix_word == 0) then
      call warning ('BAD "IS_ON" SWITCH FOR: ' // ele%name)
      return
    endif

  case default   ! normal attribute

    call evaluate_value (trim(ele%name) // ' ' // word, value, &
                                      ring, delim, delim_found, err_flag)
    if (err_flag) return

    if (i >= a0$ .and. i <= b20$) then  ! multipole attribute
        if (.not. associated(ele%a)) call multipole_init (ele)
        if (i >= b0$) then
          ele%b(i-b0$) = value
        else
          ele%a(i-a0$) = value
        endif
    elseif (i == mat6_calc_method$) then
      ele%mat6_calc_method = nint(value)
    elseif (i == field_calc$) then
      ele%field_calc = nint(value)
    elseif (i == tracking_method$) then
      ele%tracking_method = nint(value)
    elseif (i == num_steps$) then
      ele%num_steps = nint(value)
    elseif (i == integration_ord$) then
      ele%integration_ord = nint(value)
    elseif (i == ptc_kind$) then
      ele%ptc_kind = nint(value)
    elseif (i == aperture_at$) then
      ele%aperture_at = nint(value)
    elseif (i == ran_seed$) then
      call ran_seed (nint(value))  ! init random number generator
    else
      ele%value(i) = value
      if (i == b_field$) ele%field_master = .true.
    endif

  end select

  err_flag = .false.

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine get_called_file (delim)

  implicit none

  character(1) delim
  character(200) call_file

  integer ix_word, ix
  logical delim_found, finished

!

  if (delim /= ',')  call warning ('"CALL" NOT FOLLOWED BY COMMA', ' ')
  call get_next_word(call_file, ix_word, ':=,', delim, delim_found, .true.)
  if (ix_word == 0) then
    call warning ('NOTHING AFTER "CALL"', ' ')
  elseif (index('FILENAME', call_file(:ix_word)) /= 1) then
    call warning ('INVALID "CALL" COMMAND', ' ')
  elseif (delim /= '=') then
    call warning ('NO "=" AFTER "FILENAME"', ' ')
  else
    call get_next_word(call_file, ix_word, ',', &
                                       delim, delim_found, .false.)
    if (ix_word == 0) then
      call warning ('NO FILE NAME SPECIFIED', ' ')
    else
      if (call_file(1:1) == '"') then
        call_file = call_file(2:)
        ix = index(call_file, '"')
        if (ix == 0 .or. ix /= len_trim(call_file)) then
          call warning ('MISSING DOUBLE QUOTE MARK (") FOR CALL STATEMENT')
          return
        endif
        call_file(ix:ix) = ' '
      endif
      call file_stack ('push', call_file, finished)
      if (.not. bmad_status%ok) return
    endif
  endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine add_taylor_term (ele, i_out, coef, expn) 
!
! Subroutine used by bmad_parser and bmad_parser2 to parse the input file.
! This subroutine is not intended for general use.
!-

subroutine add_taylor_term (ele, i_out, coef, expn)

  implicit none

  type (ele_struct) ele
  type (taylor_term_struct), allocatable :: term(:)

  real(rp) coef

  integer i_out, expn(6), n, i


!

  if (.not. associated(ele%taylor(1)%term)) then
    allocate (ele%taylor(1)%term(0), ele%taylor(2)%term(0))
    allocate (ele%taylor(3)%term(0), ele%taylor(4)%term(0))
    allocate (ele%taylor(5)%term(0), ele%taylor(6)%term(0))
  endif

  if (i_out < 1 .or. i_out > 6) then
    call warning ('"OUT" VALUE IN TAYLOR TERM NOT IN RANGE (1 - 6)', &
                  'FOR TAYLOR ELEMENT: ' // ele%name)
    return
  endif

! find if we already have a taylor term like this.
! if so just substitute in for the new coef

  n = size(ele%taylor(i_out)%term)

  do i = 1, n
    if (all(ele%taylor(i_out)%term(i)%exp == expn)) then
      ele%taylor(i_out)%term(i)%coef = coef
      return
    endif
  enddo

! new term

  allocate (term(n+1))
  term(1:n) = ele%taylor(i_out)%term
  term(n+1)%coef = coef
  term(n+1)%exp = expn

  deallocate (ele%taylor(i_out)%term)
  allocate (ele%taylor(i_out)%term(n+1))
  ele%taylor(i_out)%term = term

  deallocate(term)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_next_word (word, ix_word,
!                        delim_list, delim, delim_found, upper_case_word)
!
! Subroutine to get the next word from the input stream.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   word            -- Character*(*): Word returned
!   delim_list      -- Character*(*): List of valid delimiters
!   upper_case_word -- Logical, optional: if True then convert word to 
!                       upper case. Default is True.
!
! Output
!   ix_word     -- Integer: length of WORD
!   delim       -- Character*1: Actual delimiter found
!   delim_found -- Logical: Set true if a delimiter found. A delimiter
!                    may not be found if the end of the line is reached first.
!-


subroutine get_next_word (word, ix_word, delim_list, &
                                    delim, delim_found, upper_case_word)

  implicit none

  integer ix_a, ix_word

  character*(*) word, delim_list, delim
                           
  logical delim_found, file_end, to_upper
  logical, optional :: upper_case_word

! check for continuation character and if found then load more characters
! into the parse line.
! after that get the first word in BP_COM%PARSE_LINE

  do
    ix_a = index(bp_com%parse_line, '&')
    if (ix_a == 0 .or. ix_a > 140) exit
    call load_parse_line('continue', ix_a, file_end)
  enddo

  call word_read (bp_com%parse_line, delim_list,  &
                         word, ix_word, delim, delim_found, bp_com%parse_line)

  if (present(upper_case_word)) then
    if (upper_case_word) call str_upcase (word, word)
  else
    call str_upcase (word, word)
  endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine file_stack (how, file_name_in, finished)
!
! Subroutine to keep track of the files that are opened for reading.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine file_stack (how, file_name_in, finished)

  implicit none

  integer, parameter :: f_maxx = 10
  type (stack_file_struct), save :: file(0:f_maxx)

  integer ix, i_level, ios

  character(*) how, file_name_in
  character(200) file_name, basename
  logical finished, found_it

  save i_level

!

  finished = .false.

  if (how == 'init') then
    i_level = 0
    return

  elseif (how == 'push') then
    i_level = i_level + 1
    if (i_level > f_maxx) then
      print *, 'ERROR: CALL NESTING GREATER THAN 10 LEVELS'
      call err_exit
    endif
    ix = splitfilename (file_name_in, file(i_level)%dir, basename)
    bp_com%dirs(2) = file(i_level-1)%dir
    call find_file (file_name_in, found_it, file_name, bp_com%dirs)
    file(i_level)%logical_name = file_name_in
    file(i_level)%full_name = file_name
    file(i_level)%f_unit = lunget()
    bp_com%current_file = file(i_level)
    bp_com%calling_file = file(i_level-1)
    if (i_level /= 1) file(i_level-1)%i_line = bp_com%current_file%i_line
    bp_com%current_file%i_line = 0

    open (bp_com%current_file%f_unit, file = file_name,  &
                                 status = 'OLD', action = 'READ', iostat = ios)
    if (ios /= 0 .or. .not. found_it) then
      print *, 'ERROR IN ', trim(bp_com%parser_name)
      print *, '      UNABLE TO OPEN FILE: ', trim(file_name)
      if (file_name_in /= file_name)  print *, &
              '       THIS FROM THE LOGICAL FILE NAME: ', trim(file_name_in)
      if (bmad_status%exit_on_error) call err_exit
      bmad_status%ok = .false.
      return
    endif

    bp_com%n_files = bp_com%n_files + 1
    inquire (file = file_name, name = bp_com%file_name_(bp_com%n_files))

  elseif (how == 'pop') then
    close (unit = bp_com%current_file%f_unit)
    i_level = i_level - 1
    if (i_level < 0) then
      call error_exit ('BAD "RETURN"', ' ')
    elseif (i_level > 0) then
      bp_com%current_file = file(i_level)
    else    ! i_level == 0
      finished = .true.
    endif
  else
    print *, 'BMAD_PARSER: INTERNAL ERROR IN FILE_STACK SUBROUTINE'
    call err_exit
  endif


end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine load_parse_line (how, ix_cmd, file_end) 
!
! Subroutine to load characters from the input file.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine load_parse_line (how, ix_cmd, file_end)

  implicit none

  integer ix_cmd, ix, ios

  character*(*) how
  character*140 line, pending_line

  logical :: cmd_pending = .false., file_end

!

  file_end = .false.

1000    continue
    if (cmd_pending) then
      line = pending_line
      cmd_pending = .false.
    else
      read (bp_com%current_file%f_unit, '(a)', end = 9000) line
      bp_com%current_file%i_line = bp_com%current_file%i_line + 1
    endif

    if (how == 'continue') then
      bp_com%input_line1 = bp_com%input_line2
      bp_com%input_line2 = line
    elseif (how == 'normal') then
      bp_com%input_line1 = ' '
      bp_com%input_line2 = line
    else
      call error_exit ('INTERNAL ERROR #4: CALL HELP', ' ')    
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
  if (ix == 0 .and. .not. cmd_pending) goto 1000

  bp_com%parse_line(ix_cmd:) = line

  return

!

  9000  continue
  file_end = .true.
  bp_com%parse_line = ' '

end subroutine

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
  
end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine evaluate_value (err_str, value, ring, delim, delim_found, err_flag)
!
! This routine creates an "evaluation stack" structure which can be used 
! to evaluate an arithmethic expression.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine evaluate_value (err_str, value, ring, delim, delim_found, err_flag)

  use random_mod

  implicit none

  type (ring_struct)  ring
  type (eval_stack_struct) stk(200)

  integer i_lev, i_op, i

  integer op_(200), ix_word, i_delim, i2, ix0

  real(rp) value

  character(*) err_str
  character(1) delim
  character(40) word, word0

  logical delim_found, split, ran_function_pending
  logical err_flag, op_found

! The general idea is to rewrite the expression on a stack in reverse polish.
! Reverse polish means that the operand goes last so that 2 * 3 is writen 
! on the stack as: [2, 3, *]

! The stack is called: stk
! Since operations move towards the end of the stack we need a separate
! stack called op_ which keeps track of what operations have not yet
! been put on stk.

! init

  err_flag = .false.
  i_lev = 0
  i_op = 0
  ran_function_pending = .false.

! parsing loop to build up the stack.

  parsing_loop: do

! get a word

    call get_next_word (word, ix_word, '+-*/()^,:}', delim, delim_found)

    if (delim == '*' .and. word(1:1) == '*') then
      call warning ('EXPONENTIATION SYMBOL IS "^" AS OPPOSED TO "**"!',  &
                    'for: ' // err_str)
      err_flag = .true.
      return
    endif

    if (ran_function_pending .and. (ix_word /= 0 .or. delim /= ')')) &
          call error_exit ('RAN AND RAN_GAUSS DO NOT TAKE AN ARGUMENT', 'FOR: ' // err_str)


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
      call get_next_word (word, ix_word, '+-*/()^,:}', delim, delim_found)
      word = word0(:ix0) // word
      ix_word = ix_word + ix0
    endif

! Now see what we got...

! for a "(" delim

    if (delim == '(') then

      ran_function_pending = .false.
      if (ix_word /= 0) then
        if (word == 'SIN') then
          call pushit (op_, i_op, sin$)
        elseif (word == 'COS') then
          call pushit (op_, i_op, cos$)
        elseif (word == 'TAN') then
          call pushit (op_, i_op, tan$)
        elseif (word == 'ASIN') then
          call pushit (op_, i_op, asin$)
        elseif (word == 'ACOS') then
          call pushit (op_, i_op, acos$)
        elseif (word == 'ATAN') then
          call pushit (op_, i_op, atan$)
        elseif (word == 'ABS') then
          call pushit (op_, i_op, abs$)
        elseif (word == 'SQRT') then
          call pushit (op_, i_op, sqrt$)
        elseif (word == 'LOG') then
          call pushit (op_, i_op, log$)
        elseif (word == 'EXP') then
          call pushit (op_, i_op, exp$)
        elseif (word == 'RAN') then
          call pushit (op_, i_op, ran$)
          ran_function_pending = .true.
        elseif (word == 'RAN_GAUSS') then
          call pushit (op_, i_op, ran_gauss$)
          ran_function_pending = .true.
        else
          call warning ('UNEXPECTED CHARACTERS ON RHS BEFORE "(": ' // word,  &
                                                  'FOR: ' // err_str)
          err_flag = .true.
          return
        endif
      endif

      call pushit (op_, i_op, l_parens$)
      cycle parsing_loop

! for a unary "-"

    elseif (delim == '-' .and. ix_word == 0) then
      call pushit (op_, i_op, unary_minus$)
      cycle parsing_loop

! for a unary "+"

    elseif (delim == '+' .and. ix_word == 0) then
      call pushit (op_, i_op, unary_plus$)
      cycle parsing_loop

! for a ")" delim

    elseif (delim == ')') then
      if (ix_word == 0) then
        if (.not. ran_function_pending) call error_exit  &
              ('CONSTANT OR VARIABLE MISSING BEFORE ")"', 'FOR: ' // err_str)
        ran_function_pending = .false.
      else
        call word_to_value (word, ring, value)
        call pushit (stk%type, i_lev, numeric$)
        stk(i_lev)%value = value
      endif

      do
        do i = i_op, 1, -1     ! release pending ops
          if (op_(i) == l_parens$) exit          ! break do loop
          call pushit (stk%type, i_lev, op_(i))
        enddo

        if (i == 0) then
          call warning ('UNMATCHED ")" ON RHS', 'FOR: ' // err_str)
          err_flag = .true.
          return
        endif

        i_op = i - 1

        call get_next_word (word, ix_word, '+-*/()^,:}', delim, delim_found)
        if (ix_word /= 0) then
          call warning ('UNEXPECTED CHARACTERS ON RHS AFTER ")"',  &
                                                    'FOR: ' // err_str)
          err_flag = .true.
          return
        endif

        if (delim /= ')') exit  ! if no more ')' then no need to release more
      enddo


      if (delim == '(') then
        call warning  &
                    ('")(" CONSTRUCT DOES NOT MAKE SENSE FOR: ' // err_str)
        err_flag = .true.
        return
      endif

! For binary "+-/*^" delims

    else
      if (ix_word == 0) then
        call warning ('CONSTANT OR VARIABLE MISSING FOR: ' // err_str)
        err_flag = .true.
        return
      endif
      call word_to_value (word, ring, value)
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
        call error_exit ('INTERNAL ERROR #01: GET HELP', ' ')
    end select

! now see if there are operations on the OP_ stack that need to be transferred
! to the STK_ stack

    do i = i_op, 1, -1
      if (eval_level(op_(i)) >= eval_level(i_delim)) then
        call pushit (stk%type, i_lev, op_(i))
      else
        exit
      endif
    enddo

! put the pending operation on the OP_ stack

    i_op = i
    if (i_delim == no_delim$) then
      exit parsing_loop
    else
      call pushit (op_, i_op, i_delim)
    endif

  enddo parsing_loop

!------------------------------------------------------------------
! now go through the stack and perform the operations

  if (i_op /= 0) then
    call warning ('UNMATCHED "(" ON RHS', 'FOR: ' // err_str)
    err_flag = .true.
    return
  endif

  if (i_lev == 0) call warning ('NO VALUE FOUND FOR: ' // err_str)

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
      if (stk(i2)%value == 0) call error_exit  &
              ('DIVIDE BY 0 ON RHS', 'FOR: ' // err_str)
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
      call random_number(stk(i2)%value)
    elseif (stk(i)%type == ran_gauss$) then
      i2 = i2 + 1
      call ran_gauss(stk(i2)%value)
    else
      call error_exit ('INTERNAL ERROR #02: GET HELP', ' ')
    endif
  enddo


  if (i2 /= 1) call error_exit ('INTERNAL ERROR #03: GET HELP', ' ')

  value = stk(1)%value

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

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

end subroutine
                       
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine word_to_value (word, ring, value)

  implicit none

  type (ring_struct), target ::  ring
  type (ele_struct), pointer :: ele

  integer i, ix1, ix2, ix_word, ios, ix
  real(rp) value
  real(rp), pointer :: ptr
  character*(*) word
  character(16) name
  logical err_flag

! see if this is numeric

  if (index('-+.0123456789', word(1:1)) /= 0) then
    read (word, *, iostat = ios) value
    if (ios /= 0) call warning ('BAD VARIABLE: ' // word)
    return
  endif

! If not numeric...

  ix_word = len_trim(word)
  call verify_valid_name (word, ix_word)

! If word does not have a "[...]" then it must be a variable

  ix1 = index(word, '[')
  if (ix1 == 0) then   
    call find_indexx (word, bp_com%var_name, &
                                    bp_com%var_indexx, bp_com%ivar_tot, i)
    if (i == 0) then
      call warning ('VARIABLE USED BUT NOT YET DEFINED: ' // word)
    else
      value = bp_com%var_value(i)
    endif
    return
  endif

! Here if word does have a "[...]" then is a element attribute

  name = word(:ix1-1)    ! name of attribute

  if (name == beam_ele%name) then
    ele => beam_ele

  else
    do i = 0, ring%n_ele_max
      if (ring%ele_(i)%name == name) then
        ele => ring%ele_(i)
        exit
      endif
    enddo

    if (i == ring%n_ele_max + 1) then
      call warning ('ELEMENT NOT DEFINED: ' // name)
      value = 0
      return
    endif

  endif

  ix2 = index(word, ']')
  name = word(ix1+1:ix2-1)

  if (name == 'S' .and. bp_com%parser_name /= 'BMAD_PARSER2') then
    call warning ('"S" ATTRIBUTE CAN ONLY BE USED WITH BMAD_PARSER2')
  endif

  call pointer_to_attribute (ele, name, .false., ptr, ix, err_flag, .false.)
  if (err_flag) then
    call warning('BAD ATTRIBUTE NAME: ' // word)
  else
    value = ptr
  endif

  return

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_variable (word, ring)

  implicit none

  type (ring_struct) ring
  character(*) word
  character(1) delim
  integer i, ivar, n, ixm, ixm2
  logical delim_found, err_flag

!

  call find_indexx (word, bp_com%var_name, &
                                    bp_com%var_indexx, bp_com%ivar_tot, i)
  if (i /= 0) then
    call warning ('VARIABLES ARE NOT ALLOWED TO BE REDEFINED: ' // word)
    call evaluate_value (word, bp_com%var_value(i), ring, &
                                              delim, delim_found, err_flag)
    return
  endif

  bp_com%ivar_tot = bp_com%ivar_tot + 1
  if (bp_com%ivar_tot > size(bp_com%var_name)) call reallocate_bp_com_var()
  ivar = bp_com%ivar_tot
  bp_com%var_name(ivar) = word
  call evaluate_value (bp_com%var_name(ivar), bp_com%var_value(ivar), &
                                       ring, delim, delim_found, err_flag)
  if (delim_found .and. .not. err_flag) call warning  &
                    ('EXTRA CHARACTERS ON RHS: ' // bp_com%parse_line,  &
                     'FOR VARIABLE: ' // bp_com%var_name(ivar))

! Reindex.

  call find_indexx (word, bp_com%var_name, bp_com%var_indexx, ivar-1, ixm, ixm2)

  bp_com%var_indexx(ixm2+1:ivar) = bp_com%var_indexx(ixm2:ivar-1)
  bp_com%var_indexx(ixm2) = ivar

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine type_get (ele, ix_type, delim, delim_found)

  implicit none

  type (ele_struct)  ele

  integer ix, ix_word, ix_type
  character word*16, delim*1, type_name*200
  logical delim_found

!

  call string_trim(bp_com%parse_line, bp_com%parse_line, ix)

  if (bp_com%parse_line(1:1) == '"') then
    bp_com%parse_line = bp_com%parse_line(2:)
    ix = index(bp_com%parse_line, '"')
    if (ix == 0) then
      call warning ('MISSING DOUBLE QUOTE MARK (") FOR TYPE = "attribute"',  &
                          'FOR ELEMENT: ' // ele%name)
      type_name = ' '
    else
      type_name = bp_com%parse_line(1:ix-1)
      bp_com%parse_line = bp_com%parse_line(ix+1:)
      call get_next_word (word, ix_word, ',=', delim, delim_found, .true.)
      if (ix_word /= 0) call warning (  &
                'EXTRA CHARACTERS FOUND AFTER TYPE ATTRIBUTE: ' // word,  &
                'FOR ELEMENT: ' // ele%name)
    endif
  else
    call get_next_word (type_name, ix_word, ',= ', delim, delim_found, .true.)
  endif

  select case (ix_type)
  case (type$)
    ele%type = type_name
  case (alias$)
    ele%alias = type_name
  case (descrip$)
    if (.not. associated(ele%descrip)) allocate (ele%descrip) 
    ele%descrip = type_name
  case (sr_wake_file$) 
    if (.not. associated(ele%wake%sr_file)) allocate (ele%wake%sr_file)
    ele%wake%sr_file = type_name
    call read_sr_wake (ele)
  case (lr_wake_file$) 
    if (.not. associated(ele%wake%lr_file)) allocate (ele%wake%lr_file)
    ele%wake%lr_file = type_name
    print *, 'ERROR: LR_WAKES NOT YET IMPLEMENTED!'
    call err_exit
    call read_lr_wake (ele)
  case default
    print *, 'INTERNAL ERROR IN TYPE_GET: I NEED HELP!'
    call err_exit
  end select

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine read_sr_wake (ele)
!
! Subroutine to read in a wake field from an external file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   ele -- Ele_struct: Element
!     %wake%sr_file -- Name of wake field file
!
! Output:
!   ele -- Ele_struct: Element with wake information.
!     %wake%sr(:) -- Short-range wake potential.
!-
        
subroutine read_sr_wake (ele)

  implicit none

  type (ele_struct) ele
  type (sr_wake_struct), allocatable :: sr(:), sr2(:)
  real(rp) dz
  integer i, j, ix, iu, ios

  character(200) file_name, line

  logical found_it

! open file

  iu = lunget()
  bp_com%dirs(2) = bp_com%calling_file%dir
  call find_file (ele%wake%sr_file, found_it, file_name, bp_com%dirs)
  open (iu, file = file_name, status = 'OLD', action = 'READ', iostat = ios)
  if (ios /= 0) then
    call warning ('CANNOT OPEN WAKE FILE: ' // ele%wake%sr_file, &
                            'FOR LCAVITY: ' // ele%name)
    return
  endif

! read

  allocate (sr(0:500))

  i = -1

  do
    read (iu, '(a)', iostat = ios) line
    if (ios < 0) exit  ! end-of-file
    if (ios > 0) then
      call warning ('ERROR READING WAKE FILE: ' // ele%wake%sr_file, &
                            'FOR LCAVITY: ' // ele%name)
      return
    endif
    call string_trim (line, line, ix)
    if (line(1:1) == '!') cycle  ! skip comments.
    if (ix == 0) cycle          ! skip blank lines.

    if (i == ubound(sr, 1)) then
      allocate (sr2(0:i))
      sr2 = sr
      deallocate (sr)
      allocate (sr(0:i+500))
      sr(0:i) = sr2
      deallocate(sr2)
    endif

    i = i + 1
    read (line, *, iostat = ios) sr(i)%z, sr(i)%long, sr(i)%trans

    if (ios /= 0) then
      call warning ('ERROR PARSING WAKE FILE: ' // ele%wake%sr_file, &
                           'CANNOT READ LINE: ' // line)
      return
    endif

  enddo

  close (iu)
  if (associated(ele%wake%sr)) deallocate (ele%wake%sr)
  allocate (ele%wake%sr(0:i))
  ele%wake%sr = sr(0:i)

! err check

  if (ele%wake%sr(0)%z /= 0) then
    call warning ('WAKEFIELDS DO NOT START AT Z = 0!', &
                                    'IN FILE: ' // ele%wake%sr_file)
    return
  endif

  dz = ele%wake%sr(i)%z / i
  do j = 1, i
    if (abs(ele%wake%sr(j)%z - dz * j) > 1e-4 * dz) then
      write (line, '(a, i5)') &
                      'WAKEFIELD POINTS DO NOT HAVE UNIFORM DZ FOR POINT:', j
      call warning (line, 'IN FILE: ' // ele%wake%sr_file)
      return
    endif
  enddo               

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine read_lr_wake (ele)
!
! Subroutine to read in a wake field from an external file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   ele -- Ele_struct: Element
!     %wake%lr_file -- Name of wake field file
!
! Output:
!   ele -- Ele_struct: Element with wake information.
!     %wake%lr(:) -- Short-range wake potential.
!-
        
subroutine read_lr_wake (ele)

  implicit none

  type (ele_struct) ele

!

  print *, 'ERROR: READ_LR_WAKE NOT YET IMPLEMENTED!'
  call err_exit

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine GET_OVERLAY_GROUP_NAMES (ELE, RING, PRING, DELIM, DELIM_FOUND)
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-
        
subroutine get_overlay_group_names (ele, ring, pring, delim, delim_found)

  implicit none

  type (ele_struct)  ele
  type (parser_ring_struct) pring
  type (ring_struct)  ring

  real(rp) coef_(200)
  real(rp) value
  
  integer ic, ix_word, ixs, j, k
                             
  character delim*1, word_in*40, word*40
  character(16) name_(200), attrib_name_(200)

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
      call evaluate_value (trim(ele%name), value, &
                                  ring, delim, delim_found, err_flag)
      if (err_flag) then
        call warning ('BAD COEFFICIENT: ' // word_in,  &
                                          'FOR ELEMENT: ' // ele%name)
        call load_parse_line ('normal', 1, file_end)         ! next line
        return
      endif
      coef_(ixs) = value
    else
      coef_(ixs) = 1.0
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
  allocate (pring%ele(ic)%coef_(ixs), pring%ele(ic)%name_(ixs), &
                                       pring%ele(ic)%attrib_name_(ixs))
  pring%ele(ic)%coef_ = coef_(1:ixs)
  pring%ele(ic)%name_ = name_(1:ixs)
  pring%ele(ic)%attrib_name_ = attrib_name_(1:ixs)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-


subroutine verify_valid_name (name, ix_name)

  implicit none

  integer i, ix_name, len, ix, ix1, ix2

  character*(*) name
  character*27 letters / '\ABCDEFGHIJKLMNOPQRSTUVWXYZ' / 
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
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine preparse_element_init (ele)

  implicit none

  type (ele_struct) ele

!

  if (ele%key == wiggler$) then
    ele%sub_key = periodic_type$   ! default
    ele%value(polarity$) = 1.0     ! default
  endif

  if (ele%key == taylor$) then
    ele%tracking_method = taylor$  ! default
    ele%mat6_calc_method = taylor$ ! default
    call add_taylor_term (ele, 1, 1.0_rp, (/ 1, 0, 0, 0, 0, 0 /))
    call add_taylor_term (ele, 2, 1.0_rp, (/ 0, 1, 0, 0, 0, 0 /))
    call add_taylor_term (ele, 3, 1.0_rp, (/ 0, 0, 1, 0, 0, 0 /))
    call add_taylor_term (ele, 4, 1.0_rp, (/ 0, 0, 0, 1, 0, 0 /))
    call add_taylor_term (ele, 5, 1.0_rp, (/ 0, 0, 0, 0, 1, 0 /))
    call add_taylor_term (ele, 6, 1.0_rp, (/ 0, 0, 0, 0, 0, 1 /))
  endif


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

  character*(*) what1
  character*(*), optional :: what2

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

subroutine warning (what1, what2, seq)

  implicit none

  integer ix
  character*(*) what1
  character*(*), optional :: what2
  type (seq_struct), optional :: seq

! BP_COM%ERROR_FLAG is a common logical used so program will stop at end of parsing

  if (bmad_status%type_out) then

    print *, 'ERROR IN ', trim(bp_com%parser_name), ': ', trim(what1)

    if (present(what2)) print '(22x, a)', trim(what2)

    if (present(seq)) then
      print *, '      IN FILE: ', trim(seq%file_name)
      print *, '      AT LINE:', seq%ix_line
    else
      print *, '      IN FILE: ', trim(bp_com%current_file%full_name)
      print *, '      AT OR BEFORE LINE:', bp_com%current_file%i_line
    endif

    if (bp_com%input_line_meaningful) then
       if (len_trim(bp_com%input_line1) /= 0) print '(5x, a)', trim(bp_com%input_line1)
       if (len_trim(bp_com%input_line2) /= 0) print '(5x, a)', trim(bp_com%input_line2)
    endif

    print *

  endif

  bp_com%error_flag = .true.

end subroutine

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

  nn = 20  ! number of "constant" variables
  bp_com%ivar_init = nn + ubound(calc_method_name, 1) + &
                                             ubound(aperture_at_name, 1)
  bp_com%ivar_tot = bp_com%ivar_init

  nt = bp_com%ivar_tot
  allocate (bp_com%var_name(nt), bp_com%var_value(nt), bp_com%var_indexx(nt))

  bp_com%var_name(1)  = 'PI'
  bp_com%var_value(1) = pi

  bp_com%var_name(2)  = 'TWOPI'
  bp_com%var_value(2) = twopi

  bp_com%var_name(3)  = 'DEGRAD'
  bp_com%var_value(3) = 180 / pi

  bp_com%var_name(4)  = 'RADDEG'
  bp_com%var_value(4) = pi / 180

  bp_com%var_name(5)  = 'E'
  bp_com%var_value(5) = 2.718281828459

  bp_com%var_name(6)  = 'E_MASS'
  bp_com%var_value(6) = e_mass

  bp_com%var_name(7)  = 'C_LIGHT'
  bp_com%var_value(7) = c_light

  bp_com%var_name(8)  = 'POSITRON'
  bp_com%var_value(8) = positron$

  bp_com%var_name(9)  = 'ELECTRON'
  bp_com%var_value(9) = electron$

  bp_com%var_name(10)  = 'R_P'
  bp_com%var_value(10) = r_p

  bp_com%var_name(11)  = 'E_CHARGE'
  bp_com%var_value(11) = e_charge

  bp_com%var_name(12)  = 'EMASS'      ! old style
  bp_com%var_value(12) = e_mass

  bp_com%var_name(13)  = 'CLIGHT'     ! old style
  bp_com%var_value(13) = c_light

  bp_com%var_name(14)  = 'LINEAR_LATTICE'
  bp_com%var_value(14) = linear_lattice$

  bp_com%var_name(15)  = 'CIRCULAR_LATTICE'
  bp_com%var_value(15) = circular_lattice$

  bp_com%var_name(16)  = 'R_E'
  bp_com%var_value(16) = r_e

  bp_com%var_name(17)  = 'PROTON'
  bp_com%var_value(17) = proton$

  bp_com%var_name(18)  = 'ANTIPROTON'
  bp_com%var_value(18) = antiproton$

  bp_com%var_name(19)  = 'M_ELECTRON'
  bp_com%var_value(19) = m_electron

  bp_com%var_name(20)  = 'M_PROTON'
  bp_com%var_value(20) = m_proton

  do i = 1, ubound(calc_method_name, 1)
    nn = nn + 1
    call str_upcase (bp_com%var_name(nn), calc_method_name(i))
    bp_com%var_value(nn) = i
  enddo

  do i = 1, ubound(aperture_at_name, 1)
    nn = nn + 1
    call str_upcase (bp_com%var_name(nn), aperture_at_name(i))
    bp_com%var_value(nn) = i
  enddo

  call indexx (bp_com%var_name(1:nt), bp_com%var_indexx(1:nt))

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

  character(16) :: var_name_temp(size(bp_com%var_name))
  real(rp) :: var_value_temp(size(bp_com%var_value))
  integer :: var_indexx_temp(size(bp_com%var_indexx))

  integer n

!

  var_name_temp = bp_com%var_name
  var_value_temp = bp_com%var_value
  var_indexx_temp = bp_com%var_indexx

  deallocate (bp_com%var_name, bp_com%var_value, bp_com%var_indexx)

  n = bp_com%ivar_tot+200
  allocate (bp_com%var_name(n), bp_com%var_value(n), bp_com%var_indexx(n))

  n = size(var_indexx_temp)
  bp_com%var_name(1:n) = var_name_temp
  bp_com%var_value(1:n) = var_value_temp
  bp_com%var_indexx(1:n) = var_indexx_temp


end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine compute_super_lord_s (ele, ring, pring)

  implicit none

  type (ring_struct)   ring
  type (ele_struct)    ele
  type (parser_ring_struct) pring

  integer ic, k
  logical found

! compute position of super_lord element
! simple case

  ic = ele%ixx

  if (ring%n_ele_max < ring%n_ele_use) then
    call warning ('N_ELE_MAX LESS THAN n_ele_use!')
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
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine add_all_superimpose (ring, ele_in, pele)

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele, ele_in, ele2, this_ele
  type (parser_ele_struct) pele

  integer ix, i, j, it, nic, nn, i_sup, i_ele, ic
  integer n_inserted, ix_lord, n_con

  character(16) matched_name(200)

  logical match_wild, have_inserted, create_new

! init

  ele = ele_in   ! in case ele changes
  ele2 = ele
  n_inserted = 0
  ring%ele_%internal_logic = .false.    ! to keep track of where we have inserted

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
      if (ic == i_beam_lord$) cycle
      if (this_ele%internal_logic) cycle

      if (match_wild(this_ele%name, pele%ref_name)) then

        do i = 1, n_inserted
          if (this_ele%name == matched_name(i)) cycle ele_loop
        enddo
       
        ring%ele_(i_ele)%internal_logic = .true.
        call compute2_super_lord_s (ring, i_ele, ele2, pele)
        call string_trim(ele%name, ele%name, ix)
        ele2%name = ele%name(:ix)            
        call add_superimpose (ring, ele2, i_sup)

        do i = i_ele, ring%n_ele_use
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
    print *, 'ERROR IN INSERT_MUTIPLE: ELEMENT ', ring%ele_(i)%name
    print *, '      IS USED WITH THE "COMMON_LORD" ATTRIBUTE BUT'
    print *, '      THIS ELEMENT IS NOT A MULIPOLE OR AB_MULTIPOLE'
    call err_exit
  endif

  ring%n_ele_max = ring%n_ele_max + 1
  if (ring%n_ele_max > ubound(ring%ele_, 1)) call allocate_ring_ele_(ring)

  nn = ring%n_ele_max 

  n_con = ring%n_control_max 
  ring%n_control_max = n_con + n_inserted
  if (ring%n_control_max > size(ring%control_)) &
   ring%control_ => reallocate_control_(ring%control_, ring%n_control_max+1000)


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
        print *, 'ERROR IN INSERT_MUTIPLE: SLAVE ', ring%ele_(i)%name
        print *, '      OF LORD ', ele%name
        print *, '      IS NOT A "FREE" ELEMENT BUT IS: ', control_name(it)
        call err_exit
      endif
      j = j + 1
      ring%ele_(i)%control_type = super_slave$
      nic = ring%n_ic_max + 1
      if (nic > size(ring%ic_)) ring%ic_ => reallocate (ring%ic_, nic+500)
      ring%ele_(i)%n_lord = 1
      ring%ele_(i)%ic1_lord = nic
      ring%ele_(i)%ic2_lord = nic
      ring%ic_(nic) = n_con + j
      ring%control_(n_con+j)%ix_slave = i
      ring%n_ic_max = nic
    endif
  enddo

  if (j /= n_inserted) then
    print *, 'ERROR IN INSERT_MUTIPLE: SLAVE NUMBER MISMATCH', j, n_inserted
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

subroutine compute2_super_lord_s (ring, i_ref, ele, pele)

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele
  type (parser_ele_struct) pele

  integer i_ref, i, ix, ct

  real(rp) s_ref_begin, s_ref_end

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

  ct = ring%ele_(i_ref)%control_type
  if (ct == overlay_lord$ .or. ct == i_beam_lord$) then
    s_ref_begin = 1e10
    s_ref_end = 0
    do i = ring%ele_(i_ref)%ix1_slave, ring%ele_(i_ref)%ix2_slave
      ix = ring%control_(i)%ix_slave
      s_ref_begin = min(s_ref_begin,  &
                         ring%ele_(ix)%s - ring%ele_(ix)%value(l$))
      s_ref_end = max(s_ref_end, ring%ele_(ix)%s)
    enddo
  elseif (ct == group_lord$) then
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
!+
! Subroutine find_indexx (name, names, an_indexx, n_max, ix_match, ix2_match)
!
! Subroutine to find a matching name in a list of names
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   name         -- Character(16): Name to match to.
!   names(:)     -- Character(16): Array of names.
!   an_indexx(:) -- Integer: Sorted index for names(:) array.
!                     names(an_indexx(i)) is in alphabetical order.
!   n_max        -- Integer: Use only names(1:n_max) part of array.
!
! Output:
!   ix_match  -- Integer: If a match is found then:
!                             names(ix_match) = name
!                  If no match is found then ix_match = 0.
!   ix2_match -- Integer, optional: 
!                  If a match is found then
!                              an_indexx(ix2_match) = ix_match
!                              names(an_indexx(ix2_match-1)) /= name
!                  If no match is found then 
!                    for j = an_indexx(ix2_match):
!                              names(j) > name
!                    and if ix2_match > 1 then for j = an_indexx(ix2_match-1):
!                              names(j) < name
!-

subroutine find_indexx (name, names, an_indexx, n_max, ix_match, ix2_match)

  implicit none

  integer ix1, ix2, ix3, n_max, ix_match
  integer, optional :: ix2_match
  integer an_indexx(:)

  character(16) name, names(:)
  character(16) this_name

! simple case

  if (n_max == 0) then
    ix_match = 0
    if (present(ix2_match)) ix2_match = 0
    return
  endif

!

  ix1 = 1
  ix3 = n_max

  do

    ix2 = (ix1 + ix3) / 2 
    this_name = names(an_indexx(ix2))

    if (this_name == name) then
      do ! if there are duplicate names in the list choose the first one
        if (ix2 == 1) exit
        if (names(an_indexx(ix2-1)) /= this_name) exit
        ix2 = ix2 - 1
      enddo
      ix_match = an_indexx(ix2)
      if (present(ix2_match)) ix2_match = ix2
      return
    elseif (this_name < name) then
      ix1 = ix2 + 1
    else
      ix3 = ix2 - 1
    endif
                       
    if (ix1 > ix3) then
      ix_match = 0
      if (present(ix2_match)) then
        if (this_name < name) then
          ix2_match = ix2 + 1
        else
          ix2_match = ix2
        endif
      endif
      return
    endif

  enddo

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
  character(16) name(20), word

  logical delim_found, err_flag

!

  n_arg = 0
  err_flag = .true.

  do
    call get_next_word (word, ix_word, '(): =,', delim, delim_found, .true.)
    if (ix_word == 0 .or. delim == '( :=') then
      call warning ('BAD ARGUMENT LIST FOR: ', seq_name)
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
! Subroutine seq_expand1 (sequence_, iseq_tot, ring, top_level)
!
! Subroutine to expand a sequence.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

recursive subroutine seq_expand1 (sequence_, iseq_tot, ring, top_level)

  implicit none

  type (seq_struct), target :: sequence_(:)
  type (seq_struct), pointer :: seq, sub_seq
  type (seq_ele_struct), allocatable, target :: s_ele(:)
  type (seq_ele_struct), allocatable :: s_ele2(:)
  type (seq_ele_struct), pointer :: this_ele
  type (ring_struct) ring

  integer ix_ele, iseq_tot, ix_word, ix, i, j, n, ios, i_rl
  integer, save :: ix_internal = 0

  real(rp) rcount

  character*20 word
  character delim*1, str*16, name*16, c_delim*1
              
  logical delim_found, replacement_line_here, c_delim_found
  logical err_flag, top_level

! init

  allocate (s_ele(ring%n_ele_maxx))
  s_ele%type = 0
  s_ele%ix_array = 0
  s_ele%ix_arg = 0

! save info on what file we are parsing for error messages.

  seq => sequence_(iseq_tot)
  seq%file_name = bp_com%current_file%full_name 
  seq%ix_line = bp_com%current_file%i_line

! first thing should be a "("

  call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)

  if (delim /= '(') call warning  &
        ('EXPECTING "(", GOT: ' // delim, 'FOR LINE: ' // seq%name)
  if (ix_word /= 0)  call warning  &
        ('EXTRANEOUS STUFF BEFORE "(", GOT: ' // word,  &
        'FOR LINE: ' // seq%name)

! now parse list proper

  ix_ele = 1
  this_ele => s_ele(ix_ele)

  do 

    call get_next_word (word, ix_word, ':=(,)', delim, delim_found, .true.)

    ix = index(word, '*')          ! E.g. word = '-3*LINE'
    if (ix /= 0) then
      bp_com%parse_line = word(:ix-1) // "," // bp_com%parse_line
      call evaluate_value (trim(seq%name) // ' Repetition Count', rcount, &
                              ring, c_delim, c_delim_found, err_flag)
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

! Check for a subline or replacement line.
! If there is one then save as an internal sequence.

    name = this_ele%name
    if (name /= ' ') call verify_valid_name (name, ix_word)

    replacement_line_here = .false.

    if (delim == '(') then ! subline or replacement line
      if (name == ' ') then
        ix_internal = ix_internal + 1
        write (str, '(a, i3.3)') '#Internal', ix_internal   ! unique name 
        this_ele%name = str
        iseq_tot = iseq_tot + 1
        sub_seq => sequence_(iseq_tot) 
        sub_seq%name = str
        sub_seq%type = seq%type
        if (sub_seq%type == replacement_line$) then
          ix = size (seq%dummy_arg)
          allocate (sub_seq%dummy_arg(ix), &
                sub_seq%corresponding_actual_arg(ix), this_ele%actual_arg(ix))
          sub_seq%dummy_arg = seq%dummy_arg
          this_ele%actual_arg = seq%dummy_arg
        endif
        bp_com%parse_line = '(' // bp_com%parse_line
        call seq_expand1 (sequence_, iseq_tot, ring, .false.)
      else
        replacement_line_here = .true.
        call get_sequence_args (name, this_ele%actual_arg, delim, err_flag)
        if (err_flag) return
      endif
      call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)
      if (word /= ' ') call warning &
                ('NO COMMA AFTER SUBLINE OR REPLACEMENT LINE. FOUND: ' // &
                 word, 'IN THE SEQUENCE: ' // seq%name)
    endif

    if (this_ele%name == ' ') call warning &
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
    this_ele => s_ele(ix_ele)

    if (ix_ele > n) then
      allocate (s_ele2(n))      
      s_ele2 = s_ele(1:n)
      deallocate (s_ele)
      allocate (s_ele(n+1000))
      s_ele(1:n) = s_ele2
      deallocate(s_ele2)
    endif

    if (delim == ')') exit

    if (delim /= ',') then
      call warning ('EXPECTING "," GOT: ' // delim, 'FOR LINE: ' // seq%name)
      exit
    endif
           
  enddo

! make sure there is nothing else if at top level

  if (top_level) then
    call get_next_word(word, ix_word, ':=() ', delim, delim_found, .true.)
    if (delim_found .or. ix_word /= 0) call warning  &
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
! Subroutine allocate_pring (ring, pring) 
!
! Subroutine to allocate allocatable array sizes.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

Subroutine allocate_pring (ring, pring)

  implicit none

  type (ring_struct) ring
  type (parser_ring_struct) pring
  type (parser_ele_struct), allocatable :: temp_pele(:)

  integer i, n_now

! assume all the arrays have the same size

  if (associated(pring%ele)) then
    n_now = ubound(pring%ele, 1)
    allocate (temp_pele(0:n_now))
    temp_pele = pring%ele
    deallocate (pring%ele)
    allocate (pring%ele(0:ring%n_ele_maxx))
    pring%ele(0:n_now) = temp_pele
    deallocate (temp_pele)

  else
    allocate (pring%ele(0:ring%n_ele_maxx))
    n_now = -1
  endif

! %ixx is used as a pointer from the in_ring%ele_ array to the pring%ele array

  do i = n_now+1, ubound(pring%ele, 1)
    nullify (pring%ele(i)%name_)
    pring%ele(i)%ref_name = blank
    pring%ele(i)%ref_pt  = center$
    pring%ele(i)%ele_pt  = center$
    pring%ele(i)%s       = 0
    pring%ele(i)%common_lord = .false.
    ring%ele_(i)%ixx = i
  enddo

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_add_lord (in_ring, n2, pring, ring)
!
! Subroutine to add overlay, group, and i_beam lords to the lattice.
! For overlays and groups: If multiple elements have the same name then 
! use all of them.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_lord (in_ring, n2, pring, ring)

  implicit none

  type (ring_struct), target :: in_ring, ring
  type (ele_struct), pointer :: lord
  type (parser_ring_struct) pring
  type (control_struct), pointer, save :: cs_(:) => null()

  integer ixx, i, ic, n, n2, k, k2, ix, j, ie, ij, ix1, ns, ixs
  integer ix_lord, ix_slave(1000), k_slave, k_slave_original
  integer, allocatable, save :: r_indexx(:)

  character(16), allocatable :: name_(:)
  character(16) name, name1, slave_name, attrib_name

  logical delete

! setup
! in_ring has the lords that are to be added to ring.
! we add an extra 1000 places to the arrays to give us some overhead.

  n = ring%n_ele_max + n2 + 1000

  allocate (r_indexx(n))
  allocate (name_(n))
  allocate (cs_(1000))

  ix1 = ring%n_ele_max
  name_(1:ix1) = ring%ele_(1:ix1)%name
  call indexx (name_(1:ix1), r_indexx(1:ix1)) ! get sorted list

! loop over elements

  main_loop: do n = 1, n2

    lord => in_ring%ele_(n)  ! next lord to add

!-----------------------------------------------------
! overlay and groups

    select case (lord%control_type)
    case (overlay_lord$, group_lord$)
 
      call new_control (ring, ix_lord)  ! get index in ring where lord goes
      ring%ele_(ix_lord) = lord
      ixx = lord%ixx

! Find where the slave elements are. 
! If a slave element is not in ring but is in in_ring then the slave has 
! not been used in the lattice list. In this case do not add the lord to 
! the lattice.

      j = 0 ! number of slaves found

      do i = 1, lord%n_slave

        name = pring%ele(ixx)%name_(i)
        call find_indexx (name, name_, r_indexx, ix1, k, k2)
        if (k == 0) then  ! not in ring
          if (all(in_ring%ele_(1:n2)%name /= name)) then ! Not in in_ring.
            call warning ('CANNOT FIND SLAVE FOR: ' // lord%name, &
                          'CANNOT FIND: '// name)
          endif
          cycle
        endif

! There might be more than 1 element with %name = name. 
! Loop over all elements whose name matches name.
! Put the info into the cs_ structure.

        do 
          j = j + 1
          k = r_indexx(k2)
          cs_(j)%coef = pring%ele(ixx)%coef_(i)
          cs_(j)%ix_slave = k
          cs_(j)%ix_lord = -1             ! dummy value
          attrib_name = pring%ele(ixx)%attrib_name_(i)
          if (attrib_name == blank) then
            cs_(j)%ix_attrib = lord%ix_value
          else
            ix = attribute_index(ring%ele_(k), attrib_name)
            cs_(j)%ix_attrib = ix
            if (ix < 1) then
              call warning ('BAD ATTRIBUTE NAME: ' // attrib_name, &
                            'IN ELEMENT: ' // lord%name)
            endif
          endif
          k2 = k2 + 1
          if (k2 > ix1-1) exit
          k = r_indexx(k2)
          if (ring%ele_(k)%name /= name) exit ! exit loop if no more matches
        enddo

      enddo

      lord%n_slave = j

! put the element name in the list r_indexx list

      call find_indexx (lord%name, ring%ele_(1:ix1)%name, &
                                             r_indexx(1:ix1), ix1-1, k, k2)
      ix1 = ix1 + 1
      r_indexx(k2+1:ix1) = r_indexx(k2:ix1-1)
      r_indexx(k2) = ix1
      name_(ix1) = lord%name

! create the lord

      ns = lord%n_slave

      select case (lord%control_type)
      case (overlay_lord$)
        call create_overlay (ring, ix_lord, lord%attribute_name, cs_(1:ns))
      case (group_lord$)
        call create_group (ring, ix_lord, cs_(1:ns))
      end select

!-----------------------------------------------------
! i_beam
! Create an i_beam element for each element whose name matches the
! first name in the slave list.

    case (i_beam_lord$) 

      ixx = lord%ixx
      name1 = pring%ele(ixx)%name_(1)

      call find_indexx (name1, name_, r_indexx, ix1, k_slave, k2)
      if (k_slave == 0) then
        call warning ('CANNOT FIND START ELEMENT FOR I_BEAM: ' // lord%name, &
                      'CANNOT FIND: '// name)
        cycle
      endif

      if (k_slave > ring%n_ele_use) then ! must be a super_lord.
        ix = ring%ele_(k)%ix1_slave
        k_slave = ring%control_(ix)%ix_slave
      endif

      k_slave_original = k_slave

! Loop over all matches to the first name.

      do 

        ixs = 0       ! Index of slave element we are looking for

        slave_loop: do            ! loop over all slaves
          ixs = ixs + 1
          if (ixs > lord%n_slave) exit
          slave_name = pring%ele(ixx)%name_(ixs)

          do  ! loop over all lattice elements
            if (ring%ele_(k_slave)%control_type == super_slave$) then
              do ic = ring%ele_(k_slave)%ic1_lord, ring%ele_(k_slave)%ic2_lord
                ie = ring%control_(ic)%ix_lord
                if (match_wild(ring%ele_(ie)%name, slave_name)) then
                  ix_slave(ixs) = ie
                  cycle slave_loop
                endif
              enddo
            else
              if (match_wild(ring%ele_(k_slave)%name, slave_name)) then
                ix_slave(ixs) = k_slave
                cycle slave_loop
              endif
            endif
            k_slave = k_slave + 1  
            if (k_slave == ring%n_ele_use + 1) k_slave = 1
            if (k_slave == k_slave_original) then
              call warning ('CANNOT FIND END ELEMENT FOR I_BEAM: ' // &
                                     lord%name, 'CANNOT FIND: ' // slave_name)
              cycle main_loop
            endif
          enddo 
        enddo slave_loop

! create the i_beam element

        call new_control (ring, ix_lord)
        call create_i_beam (ring, ix_lord, ix_slave(1:lord%n_slave), lord)

        k2 = k2 + 1
        k_slave = r_indexx(k2)
        if (ring%ele_(k_slave)%name /= name1) exit
        if (k_slave > ring%n_ele_use) then ! must be a super_lord.
          ix = ring%ele_(k_slave)%ix1_slave
          k_slave = ring%control_(ix)%ix_slave
        endif

      enddo 

    end select

  enddo main_loop

! cleanup

  deallocate (r_indexx)
  deallocate (name_)
  deallocate (cs_)

end subroutine

end module
