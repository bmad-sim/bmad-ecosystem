#include "CESR_platform.inc"

module bmad_parser_mod

! This is for bmad_parser and bmad_parser2

  use ptc_interface_mod
  use bookkeeper_mod

! structure for a decleared variable

  type parser_var_struct
    character(16) name                    ! variable name
    real(rdef) value                     ! variable value
  end type                      

! structure for holding the contents of lines and lists (sequences)

  type seq_ele_struct
    character(16) name              ! name of element, subline, or sublist
    integer type                   ! LINE$, LIST$, ELEMENT$
    integer ix_array               ! if an element: pointer to ELE_ array
                                   ! if a list: pointer to SEQ_ array
    integer ix_arg                 ! for replacement lines
    logical reflect                ! reflection of subline?
  end type

  type replacement_arg_struct
    character(16) dummy_name
    character(16) actual_name
  end type

! Head structure for lines and lists (sequences)

  type seq_struct
    character(16) name                  ! name of sequence
    integer type                       ! LINE$, REPLACEMENT_LINE$ or LIST$
    integer ix                         ! current index of element in %ELE
    integer indexx                     ! alphabetical order sorted index
    type (seq_ele_struct), pointer :: ele(:)
    type (replacement_arg_struct), pointer :: arg(:)
    character(200) file_name           ! file where sequence is defined
    integer ix_line                    ! line where sequence is defined
  end type

! stack structure

  type seq_stack_struct
    integer ix_seq                ! index to seq_(:) array
    integer ix_ele                ! index to seq%ele(:) array
    integer reflect               ! reflection sequence?
  end type

!-----------------------------------------------------------
! structure for holding the control names and pointers for
! superimpose and overlay elements

  type parser_ele_struct
    character(16) ref_name
    character(16), pointer :: name_(:)
    character(16), pointer :: attrib_name_(:)
    real(rdef), pointer :: coef_(:)
    integer ix_count
    integer ele_pt, ref_pt
    logical common_lord
    real(rdef) s
    integer indexx
  end type

  type parser_ring_struct
    type (parser_ele_struct) ele(0:n_ele_maxx)
  end type

! component_struct

  type component_struct
    character(16) name
    integer ix_ele
  end type

  integer ivar_maxx
  parameter (ivar_maxx =  500)

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

!

  type bp_com_struct
    integer i_line, f_unit, n_files
    character*200 current_file_name, file_name_(20)
    character*280 parse_line
    character(16) parser_name
    character*72 debug_line
    logical parser_debug, no_digested, error_flag
    integer iseq_tot, ivar_tot, ivar_init
  end type

! common stuff

  type (bp_com_struct)  bp_com
  type (parser_var_struct)  var_(ivar_maxx)
  type (ele_struct), target :: beam_ele, param_ele

  character(16) :: blank = ' '

!------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine get_attribute (how, ele, ring, pring, &
                                             delim, delim_found, err_flag)
                            
  implicit none

  type (ring_struct)  ring
  type (parser_ring_struct) pring
  type (ele_struct), target ::  ele, ele0
  type (wig_term_struct), pointer :: wig_term(:)

  real(rdef) kx, ky, kz, tol, value, coef
  real(rdef), pointer :: r_ptr

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

    call get_next_word (word, ix_word, ',}', delim, delim_found, .true.)
    read (word, *, iostat = ios) i_out
    if (delim /= ',' .or. ix_word == 0 .or. ios /= 0) then
      call warning ('BAD "OUT" IN TERM FOR TAYLOR ELEMENT: ' // ele%name, &
                                                      'CANNOT PARSE: ' // str)
      return
    endif

    call evaluate_value (str, coef, ring, delim, delim_found, err_flag)
    if (err_flag) return

    call get_next_word (line, ix_word, '},', delim, delim_found, .true.)
    read (line, *, iostat = ios) expn
    if (delim /= '}' .or. ix_word == 0 .or. ios /= 0) then
      call warning ('BAD "EXPONENT" IN TERM FOR TAYLOR ELEMENT: ' // ele%name, &
                                                      'CANNOT PARSE: ' // str)
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
      if (how == def$) ele%ix_value = i
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
    
    select case (word)
    case ('BETA_X') 
      ele%x%gamma = (1 + ele%x%alpha**2) / ele%x%beta
    case ('ALPHA_X')
      if (ele%x%beta /= 0) ele%x%gamma = (1 + ele%x%alpha**2) / ele%x%beta
    case ('BETA_Y') 
      ele%y%gamma = (1 + ele%y%alpha**2) / ele%y%beta
    case ('ALPHA_Y')
      if (ele%y%beta /= 0) ele%y%gamma = (1 + ele%y%alpha**2) / ele%y%beta
    case ('C11', 'C12', 'C21', 'C22')
      ele%gamma_c = sqrt(1 - ele%c_mat(1,1)*ele%c_mat(2,2) + &
                                              ele%c_mat(1,2)*ele%c_mat(2,1))
    end select
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
      if (ele%key == quadrupole$) then
        ele%value(tilt$) = pi / 4
      elseif (ele%key == sextupole$) then
        ele%value(tilt$) = pi / 6
      elseif (ele%key == octupole$) then
        ele%value(tilt$) = pi / 8
      else
        call warning ('SORRY I''M NOT PROGRAMMED TO USE A "TILT" DEFAULT' // &
                'FOR A: ' // key_name(ele%key), 'FOR: ' // ele%name)
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
! The TYPE, ALIAS, and DESCRIP attributes are special because their "values"
! are character strings

  if (i == type$ .or. i == alias$ .or. i == descrip$) then  
    call type_get (ele, i, delim, delim_found)

  elseif (i == symplectify$) then
    call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.)
    read (word, '(L)', iostat = ios) ele%symplectify
    if (ios /= 0 .or. ix_word == 0) then
      call warning ('BAD "SYMPLECTIFY" SWITCH FOR: ' // ele%name)
      return
    endif

  else ! normal attribute

    call evaluate_value (trim(ele%name) // ' ' // word, value, &
                                      ring, delim, delim_found, err_flag)
    if (err_flag) return

    if (i >= a0$) then  ! multipole attribute
        if (.not. associated(ele%a)) then
          allocate (ele%a(0:n_pole_maxx), ele%b(0:n_pole_maxx))
          ele%a = 0;  ele%b = 0
        endif
        if (i >= b0$) then
          ele%b(i-b0$) = value
        else
          ele%a(i-a0$) = value
        endif
    elseif (i == mat6_calc_method$) then
      ele%mat6_calc_method = nint(value)
    elseif (i == tracking_method$) then
      ele%tracking_method = nint(value)
    elseif (i == num_steps$) then
      ele%num_steps = nint(value)
    elseif (i == integration_order$) then
      ele%integration_order = nint(value)
    elseif (i == ptc_kind$) then
      ele%ptc_kind = nint(value)
    else
      ele%value(i) = value
      if (i == b_field$) ele%b_field_master = .true.
    endif
  endif

  err_flag = .false.

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine add_taylor_term (ele, i_out, coef, expn)

  implicit none

  type (ele_struct) ele
  type (taylor_term_struct), allocatable :: term(:)

  real(rdef) coef

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

  implicit none

  integer ix_a, ix_word

  character*(*) word, delim_list, delim
                           
  logical delim_found, file_end, upper_case_word

! check for continuation character and if found then load more characters
! into the parse line.
! after that get the first word in BP_COM.PARSE_LINE

  do
    ix_a = index(bp_com%parse_line, '&')
    if (ix_a == 0 .or. ix_a > 140) exit
    call load_parse_line('continue', ix_a, file_end)
  enddo

  call word_read (bp_com%parse_line, delim_list,  &
                         word, ix_word, delim, delim_found, bp_com%parse_line)

  if (upper_case_word) call str_upcase (word, word)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine file_stack (how, file_name, finished)

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
    bp_com%current_file_name = file_name
    f_unit_(i_level) = lunget()
    bp_com%f_unit = f_unit_(i_level)
    if (i_level /= 1) i_line_(i_level-1) = bp_com%i_line
    bp_com%i_line = 0
    open (unit = bp_com%f_unit, file = bp_com%current_file_name,  &
                                 status = 'OLD', action = 'READ', err = 9000)
    bp_com%n_files = bp_com%n_files + 1
    inquire (file = file_name, name = bp_com%file_name_(bp_com%n_files))
  elseif (how == 'pop') then
    close (unit = bp_com%f_unit)
    i_level = i_level - 1
    if (i_level < 0) then
      call error_exit ('BAD "RETURN"', ' ')
    elseif (i_level > 0) then
      bp_com%current_file_name = stack_file_name_(i_level)
      bp_com%f_unit = f_unit_(i_level)
      bp_com%i_line = i_line_(i_level)
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
      'ERROR IN ', trim(bp_com%parser_name), ': UNABLE TO OPEN FILE: ', &
                                       trim(bp_com%current_file_name)
  if (bmad_status%exit_on_error) call err_exit
  bmad_status%ok = .false.

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine load_parse_line (how, ix_cmd, file_end)

  implicit none

  integer ix_cmd, ix

  character*(*) how
  character*140 line, pending_line

  logical :: cmd_pending = .false., file_end


!

  file_end = .false.

! init

  if (how == 'init') then
    bp_com%parser_debug = .false.
    bp_com%no_digested = .false.
    read (bp_com%f_unit, '(a)', end = 9000) line
    if (index(line, '!PARSER_DEBUG') /= 0) then
      bp_com%parser_debug = .true.
      bp_com%debug_line = line
      print *, 'FOUND IN FILE: "!PARSER_DEBUG". DEBUG IS NOW ON'
    elseif (index(line, '!NO_DIGESTED') /= 0) then
      bp_com%no_digested = .true.
      print *, 'FOUND IN FILE: "!NO_DIGESTED". NO DIGESTED FILE WILL BE CREATED'
    endif
    rewind (unit = bp_com%f_unit)
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
      read (bp_com%f_unit, '(a)', end = 9000) line
      bp_com%i_line = bp_com%i_line + 1
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

  bp_com%parse_line(ix_cmd:) = line

  return

9000  continue
  file_end = .true.
  bp_com%parse_line = ' '

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine evaluate_value (err_str, value, &
                              ring, final_delim, final_delim_found, err_flag)

  implicit none

  type eval_stack_struct
    integer type
    real(rdef) value
  end type

  type (ring_struct)  ring
  type (eval_stack_struct) stk(200)

  integer i_lev, i_op, i

  integer :: plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
  integer :: l_parens$ = 5, power$ = 7, unary_minus$ = 8, unary_plus$ = 9
  integer :: no_delim$ = 10
  integer :: sin$ = 11, cos$ = 12, tan$ = 13
  integer :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17
  integer :: numeric$ = 100

  integer :: level(17) = (/ 1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
                            9, 9, 9, 9, 9, 9, 9 /)
  character*1 op_name(9) / '+', '-', '*', '/', '(', ')', '^', '-', ' ' /

  integer op_(200), ix_word, i_delim, i2, ix0

  real(rdef) value

  character line*70, delim*1, word*40, final_delim*1, word0*40, err_str*(*)

  logical delim_found, final_delim_found, split
  logical err_flag, op_found

! init

  err_flag = .false.
  i_lev = 0
  i_op = 0

! get line to parse

  call get_next_word (line, ix_word, ',:}', &
                                final_delim, final_delim_found, .true.)
  if (ix_word == 0) call warning  &
                         ('NO VALUE FOUND FOR: ' // err_str)

! parsing loop to build up the stack

  parsing_loop: do

! get a word

    call word_read (line, '+-*/()^',  &
                         word, ix_word, delim, delim_found, line)

    if (delim == '*' .and. line(1:1) == '*') then
      call warning ('EXPONENTIATION SYMBOL IS "^" AS OPPOSED TO "**"!',  &
                    'for: ' // err_str)
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
      if (ix_word == 0) call error_exit  &
            ('CONSTANT OR VARIABLE MISSING BEFOR ")" FOR: ' // err_str, ' ')
      call word_to_value (word, ring, value)
      call pushit (stk%type, i_lev, numeric$)
      stk(i_lev)%value = value

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

        call word_read (line, '+-*/()^',  &
                           word, ix_word, delim, delim_found, line)
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

    op_found = .false.
    do i = 1, 7
      if (delim == op_name(i)) then      ! op identified
        i_delim = i                        ! op id number
        op_found = .true.
        exit
      endif
    enddo

    if (.not. op_found) then
      if (delim_found) then   ! how could this be?
        call error_exit ('INTERNAL ERROR #01: GET HELP', ' ')
      else                    ! must be that we are at the end of the line
        i_delim = no_delim$
      endif
    endif

! now see if there are operations on the OP_ stack that need to be transferred
! to the STK_ stack

    do i = i_op, 1, -1
      if (level(op_(i)) >= level(i_delim)) then
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


subroutine word_to_value (word, ring, value)

  implicit none

  type (ring_struct), target ::  ring
  type (ele_struct), pointer :: ele

  integer i, ix1, ix2, ix_word
  real(rdef) value
  character*(*) word
  character(16) name

! see if this is numeric

  if (index('-+.0123456789', word(1:1)) /= 0) then
    read (word, *, err = 9000) value
    return
  endif

!

  ix_word = len_trim(word)
  call verify_valid_name (word, ix_word)

! If word has a "[...]" then it is a element attribute

  ix1 = index(word, '[')
  if (ix1 /= 0) then   

    name = word(:ix1-1)    ! name of attribute
    do i = 0, ring%n_ele_max
      if (ring%ele_(i)%name == name) then
        ele => ring%ele_(i)
        goto 1000
      endif
    enddo

    if (name == beam_ele%name) then
      ele => beam_ele
      goto 1000
    endif

    call warning ('ELEMENT NOT DEFINED: ' // name)
    value = 0
    return

1000 continue

    ix2 = index(word, ']')
    name = word(ix1+1:ix2-1)

    if (name == 'S') then
      if (bp_com%parser_name == 'BMAD_PARSER2') then
        value = ele%s
      else
        call warning ('"S" ATTRIBUTE CAN ONLY BE USED WITH BMAD_PARSER2')
      endif
    else
      i = attribute_index(ele, name)
      if (i < 1) call warning('BAD ATTRIBUTE NAME: ' // word)
      value = ele%value(i)
    endif

    return
  endif

! None of the above? must be a variable

  do i = 1, bp_com%ivar_tot
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


subroutine type_get (ele, ix_type, delim, delim_found)

  implicit none

  type (ele_struct)  ele

  integer ix, ix_word, ix_type
  character word*16, delim*1, type_name*200
  logical delim_found

!

  call string_trim(bp_com%parse_line, bp_com%parse_line, ix)
  ix = index(bp_com%parse_line(2:), '"')

  if (bp_com%parse_line(1:1) /= '"' .or. ix == 0) then
    call warning ('MISSING DOUBLE QUOTE MARK (") FOR TYPE = "attribute"',  &
                          'FOR ELEMENT: ' // ele%name)
    if (ix /= 0) then
      bp_com%parse_line = bp_com%parse_line(ix+2:)
    elseif (bp_com%parse_line(1:1) == '"') then
      bp_com%parse_line = bp_com%parse_line(2:)
    endif
  endif

  if (ix == 1) then
    type_name = ' '
  else
    type_name = bp_com%parse_line(2:ix)
  endif

  if (ix_type == type$) then
    ele%type = type_name
  elseif (ix_type == alias$) then
    ele%alias = type_name
  elseif (ix_type == descrip$) then
    if (associated(ele%descrip)) deallocate (ele%descrip)
    allocate (ele%descrip) 
    ele%descrip = type_name
  endif

  bp_com%parse_line = bp_com%parse_line(ix+2:)
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

  implicit none

  type (ele_struct)  ele
  type (parser_ring_struct) pring
  type (ring_struct)  ring

  real(rdef) coef_(200)
  real(rdef) value
  
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


subroutine verify_valid_name (name, ix_name, key_check)

  implicit none

  integer i, ix_name, len, ix, ix1, ix2

  character*(*) name
  character*27 letters / '\ABCDEFGHIJKLMNOPQRSTUVWXYZ' / 
  character*40   valid_chars / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\0123456789_[]' /
  character*1 tab

  logical, optional :: key_check
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

! Check if name looks like a element key.
! BEAM is OK even though it looks like BEAMBEAM

  ix = ix_name
  if (ix1 > 0) ix = ix1

  if (present(key_check) .and. name(:ix) /= 'BEAM') then
    if (key_check) then
      do i = 1, n_key
        if (name(:ix) == key_name(i)(:ix)) then
          call warning ('NAME: ' // name, &
                          'LOOKS LIKE AN ELEMENT KEY: ' // key_name(i))
          return
        endif
      enddo
    endif
  endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine error_exit (what1, what2)

  implicit none

  character*(*) what1
  character*(*), optional :: what2

!

  print *, 'FATAL ERROR IN BMAD_PARSER: ', what1
  if (present(what2)) then
    if (what2 /= ' ') print '(22x, a)', what2
  endif
  print *, '      IN FILE: ', bp_com%current_file_name(:60)
  print *, '      AT OR BEFORE LINE:', bp_com%i_line
  call err_exit

end subroutine


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine warning (what1, what2, seq)

  implicit none

  integer ix
  character*(*) what1
  character*(*), optional :: what2
  type (seq_struct), optional :: seq

! BP_COM.ERROR_FLAG is a common logical used so program will stop at end of parsing

  if (bmad_status%type_out) then
    print *, 'ERROR IN ', trim(bp_com%parser_name), ': ', what1
    if (present(what2)) print '(22x, a)', what2
    if (present(seq)) then
      print *, '      IN FILE: ', trim(seq%file_name)
      print *, '      AT LINE:', seq%ix_line
    else
      print *, '      IN FILE: ', trim(bp_com%current_file_name)
      print *, '      AT OR BEFORE LINE:', bp_com%i_line
    endif
  endif

  bp_com%error_flag = .true.

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


subroutine init_bmad_parser_common

  implicit none
  
  integer nn, i
  
!

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

  var_(17)%name = 'PROTON'
  var_(17)%value = proton$

  var_(18)%name = 'ANTIPROTON'
  var_(18)%value = antiproton$

  var_(19)%name = 'M_ELECTRON'
  var_(19)%value = m_electron

  var_(20)%name = 'M_PROTON'
  var_(20)%value = m_proton

  nn = 20

  do i = 1, ubound(calc_method_name, 1)
    call str_upcase (var_(nn+i)%name, calc_method_name(i))
    var_(nn+i)%value = i
  enddo

  bp_com%ivar_init = nn + ubound(calc_method_name, 1)
  bp_com%ivar_tot = bp_com%ivar_init


end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


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
    print *, 'ERROR IN INSERT_MUTIPLE: ELEMENT ', ring%ele_(i)%name
    print *, '      IS USED WITH THE "COMMON_LORD" ATTRIBUTE BUT'
    print *, '      THIS ELEMENT IS NOT A MULIPOLE OR AB_MULTIPOLE'
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
        print *, 'ERROR IN INSERT_MUTIPLE: SLAVE ', ring%ele_(i)%name
        print *, '      OF LORD ', ele%name
        print *, '      IS NOT A "FREE" ELEMENT BUT IS: ', control_name(it)
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
    print *, 'ERROR IN INSERT_MUTIPLE: SLAVE NUMBER MISMATCH', j, n_inserted
    call err_exit
  endif


end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine compute2_super_lord_s (ring, i_ref, ele, pele)

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele
  type (parser_ele_struct) pele

  integer i_ref, i, ix

  real(rdef) s_ref_begin, s_ref_end

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
!+
! Subroutine find_indexx (name, names, an_indexx, n_max, ix_match, ix2_match)
!
! Subroutine to find a matching name in a list of names
!
! Input:
!   name         -- Character(16): Name to match to.
!   names(:)     -- Character(16): Array of names.
!   an_indexx(:) -- Integer: Sorted index of names array.
!   n_max        -- Integer: Use only names(1:n_max) part of array.
!
! Output:
!   ix_match     -- Integer: If a match is found then:
!                                 names(ix_match) = name
!                     If no match is found then ix_match = 0.
!   ix2_match    -- Integer, optional: If a match is found then
!                                 an_indexx(ix2_match) = ix_match.
!                     If no match is found then ix2_match = 0.
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
      if (present(ix2_match)) ix2_match = 0    
      return
    endif

  enddo

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine seq_expand1 (seq_, ix_seq, top_level, ring)
!
! Subroutine to expand a sequence
!-

recursive subroutine seq_expand1 (seq_, ix_seq, top_level, ring)

  implicit none

  type (seq_struct), target :: seq_(:)
  type (seq_struct), pointer :: seq
  type (seq_ele_struct) s_ele(n_ele_maxx)
  type (ring_struct) ring

  integer ix_seq, ix_ele, ix_word, ix, count, i, j, n, ios, i_rl
  integer, save :: ix_internal = 0

  real(rdef) rcount

  character*20 word
  character delim*1, str*16, name*16, c_delim*1
              
  logical delim_found, top_level, replacement_line_here, c_delim_found
  logical err_flag

! init

  s_ele%type = 0
  s_ele%ix_array = 0
  s_ele%ix_arg = 0

! save info on what file we are parsing for error messages.

  seq => seq_(ix_seq)
  seq%file_name = bp_com%current_file_name 
  seq%ix_line = bp_com%i_line

! first thing should be a "("

  call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)

  if (delim /= '(') call warning  &
        ('EXPECTING "(", GOT: ' // delim, 'FOR LINE: ' // seq%name)
  if (ix_word /= 0)  call warning  &
        ('EXTRANEOUS STUFF BEFORE "(", GOT: ' // word,  &
        'FOR LINE: ' // seq%name)

! now parse list proper

  ix_ele = 1

  do 

    call get_next_word (word, ix_word, ':=(,)', delim, delim_found, .true.)

    ix = index(word, '*')          ! E.g. word = '-3*LINE'
    if (ix /= 0) then
      bp_com%parse_line = word(:ix-1) // "," // bp_com%parse_line
      call evaluate_value (trim(seq%name) // ' Repetition Count', rcount, &
                              ring, c_delim, c_delim_found, err_flag)
      count = nint(rcount)
      if (err_flag) return
      s_ele(ix_ele)%name = word(ix+1:)
      if (count < 0) then
        s_ele(ix_ele)%reflect = .true.
      else
        s_ele(ix_ele)%reflect = .false.
      endif
      count = abs(count)
      ix_word = ix_word - ix
    elseif (word(1:1) == '-') then
      s_ele(ix_ele)%reflect = .true.
      count = 1
      s_ele(ix_ele)%name = word(2:)
      ix_word = ix_word - 1
    else
      s_ele(ix_ele)%reflect = .false.
      count = 1
      s_ele(ix_ele)%name = word
    endif

! Check for a subline or replacement line.
! If there is one then save as an internal sequence.

    name = s_ele(ix_ele)%name
    if (name /= ' ') call verify_valid_name (name, ix_word)

    replacement_line_here = .false.
    if (delim == '(') then ! subline or replacement line
      if (name /= ' ') replacement_line_here = .true.
      ix_internal = ix_internal + 1
      write (str, '(a, i3.3)') '#Internal', ix_internal   ! unique name
      s_ele(ix_ele)%name = str
      ix_seq = ix_seq + 1
      seq_(ix_seq)%name = str
      seq_(ix_seq)%type = line$
      bp_com%parse_line = '(' // bp_com%parse_line 
      call seq_expand1 (seq_, ix_seq, .false., ring)
      call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)
      if (word /= ' ') call warning &
                ('NO COMMA AFTER SUBLINE OR REPLACEMENT LINE. FOUND: ' // &
                 word, 'IN THE SEQUENCE: ' // seq%name)
    endif

    if (s_ele(ix_ele)%name == ' ') call warning &
              ('SUB-ELEMENT NAME IS BLANK FOR LINE/LIST: ' // seq%name)

! if a replacement line then switch the real list for the actual args.

    if (replacement_line_here) then  ! replacement line 
      do i_rl = 1, ix_seq
        if (i_rl == ix_seq) then
          call warning ('CANNOT FIND REPLACEMENT LINE DEFINITION FOR: ' &
                          // name, 'WHICH APPEARS IN LINE: ' // seq%name)
          return
        endif
        if (seq_(i_rl)%name == name) exit
      enddo

      if (seq_(i_rl)%type /= replacement_line$) then
        call warning (trim(name) // ' IS USED AS A REPLACEMENT LINE IN: ' // &
              seq%name, 'BUT IT IS NOT A REPLACEMENT LINE')
        return
      endif

      if (size(seq_(i_rl)%arg) /= size(seq_(ix_seq)%ele)) then
        call warning ('NUMBER OF ARGUMENTS IN REPLACEMENT LINE: ' &
                // seq_(i_rl)%name, 'DOES NOT MATCH USE IN LINE: ' // seq%name)
        return
      endif

      seq_(i_rl)%arg(:)%actual_name = seq_(ix_seq)%ele(:)%name
      n = size(seq_(i_rl)%ele)
      deallocate (seq_(ix_seq)%ele)
      allocate (seq_(ix_seq)%ele(n))
      do j = 1, n
        ix = seq_(i_rl)%ele(j)%ix_arg
        if (ix > 0) then
          seq_(ix_seq)%ele(j)%name = seq_(i_rl)%arg(ix)%actual_name
          seq_(ix_seq)%ele(j)%ix_arg = 0
        else
          seq_(ix_seq)%ele(j) = seq_(i_rl)%ele(j)
        endif
      enddo 

   endif

! if a replacement line then look for element in argument list

    s_ele(ix_ele)%ix_arg = 0
    if (seq%type == replacement_line$) then
      do i = 1, size(seq%arg)
        if (seq%arg(i)%dummy_name == s_ele(ix_ele)%name) then
          s_ele(ix_ele)%ix_arg = i
          exit
        endif
      enddo
    endif

! if multiple repetition count then expand

    do i = 1, count-1
      s_ele(ix_ele+i) = s_ele(ix_ele)
    enddo

    ix_ele = ix_ele + count

    if (delim == ')') exit

    if (delim /= ',') then
      call warning  &
               ('EXPECTING "," GOT: ' // delim, 'FOR LINE: ' // seq%name)
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
  seq%ele(:) = s_ele(1:ix_ele)

  seq%ix = 1   ! Init. Used for replacement line index

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine find_slaves_for_parser (ring, name_, attrib_name_, coef, cs_)
!
! Subroutine to find where the slaves are in a ring for a overlay_lord or
! group_lord. If multiple elements have the same name then use all of them.
!-

subroutine find_slaves_for_parser (ring, name_, attrib_name_, coef_, cs_out)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele
  type (control_struct), pointer :: cs_out(:)
  type (control_struct) cs_(n_ele_maxx)

  real(rdef), pointer :: coef_(:)

  integer i, j, k, k2, ixl, ix
  integer r_indexx(n_ele_maxx)

  character(16), pointer :: name_(:), attrib_name_(:)
  character(16) name

! ring%n_ele_max is the index of the lord

  ixl = ring%n_ele_max
  ele => ring%ele_(ixl)

  call indexx (ring%ele_(1:ixl-1)%name, r_indexx(1:ixl-1))

  j = 0
  do i = 1, ele%n_slave

    call find_indexx (name_(i), ring%ele_(1:)%name, r_indexx, ixl-1, k, k2)
    if (k2 == 0) then
      call warning ('CANNOT FIND SLAVE FOR: ' // ele%name, &
                    'CANNOT FIND: '// name_(i))
      cycle
    endif

    do
      if (k2 == 1) exit
      if (ring%ele_(r_indexx(k2-1))%name == name_(i)) then
        k2 = k2 - 1
      else
        exit
      endif
    enddo

    do 
      j = j + 1
      k = r_indexx(k2)
      cs_(j)%coef = coef_(i)
      cs_(j)%ix_slave = k
      cs_(j)%ix_lord = -1             ! dummy value
      name = attrib_name_(i)
      if (name == blank) then
        cs_(j)%ix_attrib = ele%ix_value
      else
        ix = attribute_index(ring%ele_(k), name)
        cs_(j)%ix_attrib = ix
        if (ix < 1) then
          call warning ('BAD ATTRIBUTE NAME: ' // name, &
                        'IN ELEMENT: ' // ele%name)
        endif
      endif
      k2 = k2 + 1
      if (k2 > ixl-1) exit
      k = r_indexx(k2)
      if (ring%ele_(k)%name /= name_(i)) exit
    enddo

  enddo

  ele%n_slave = j

! init

  if (associated(cs_out)) deallocate(cs_out)
  allocate (cs_out(j))
  cs_out = cs_(1:j)

end subroutine

end module
