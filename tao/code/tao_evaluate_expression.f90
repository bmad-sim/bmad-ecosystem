!+
! Subroutine tao_evaluate_expression (expression, n_size, use_good_user, value, err_flag, print_err, &
!                   info, stack, dflt_component, dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, &
!                   dflt_dat_or_var_index, dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)
!
! Mathematically evaluates a character expression.
!
! Input:
!   expression      -- character(*): Arithmetic expression.
!   n_size          -- integer: Size of the value array. If the expression evaluates to a
!                       a scalar, each value in the value array will get this value.
!                       If n_size = 0, the natural size is determined by the expression itself.
!   use_good_user   -- logical: Use the good_user logical in evaluating good(:)
!   print_err       -- logical, optional: If False then supress evaluation error messages.
!                       This does not affect syntax error messages. Default is True.
!   dflt_component  -- character(*), optional: Component to use if not specified in the expression. 
!                        'model' (default), 'base', or 'design'.
!   dflt_source     -- character(*), optional: Default source ('lat', 'data', etc.). Default is ''.
!   dflt_ele_ref    -- ele_struct, pointer, optional: Default reference element.
!   dflt_ele_start  -- ele_struct, pointer, optional: Default start element for ranges.
!   dflt_ele        -- ele_struct, pointer, optional: Default element to evaluate at.
!   dflt_dat_or_var_index -- character(*), optional: Default datum or variable index to use.
!   dflt_uni        -- integer, optional: Default universe to use. If 0 or not present, use viewed universe.
!   dflt_eval_point -- integer, optional: Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.
!   dflt_s_offset   -- real(rp), optional: Default offset of eval_point. Default = 0.
!   dflt_orbit      -- coord_struct, optional: Default orbit to evaluate at.
!   datum           -- tao_data_struct, optional: If present, check to see that the expression does not depend upon
!                       a datum that will be evaluated after this datum. If so, this is an error.
!
! Output:
!   value(:)  -- Real(rp), allocatable: Value of arithmetic expression.
!   err_flag  -- Logical: True on an error. EG: Invalid expression.
!                  A divide by zero is not an error but good(:) will be set to False.
!   info(:)    -- tao_expression_info_struct, allocatable, optional: Is the value valid?, etc.
!                  Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                  orbit.x[23]|good_user is False.
!   stack(:)  -- Tao_eval_stack1_struct, allocatable, optional: Evaluation stack for the
!                  expression. This is useful to save if the same expression is
!                  to be evaluated repeatedly. 
!                  With this, tao_evaluate_stack can be called directly.
!-

recursive &
subroutine tao_evaluate_expression (expression, n_size, use_good_user, value, err_flag, print_err, &
                      info, stack, dflt_component, dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, &
                      dflt_dat_or_var_index, dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)

use tao_data_and_eval_mod, dummy => tao_evaluate_expression
use random_mod
use expression_mod

implicit none

type expression_func_struct
  character(12) :: name = ''      ! Name of function
  integer :: n_arg_target = 0     ! Number of arguments the function should have. -1 => 0 or 1 arg
  integer :: n_arg_count = 0      ! Number of arguments found.
end type

type (tao_eval_stack1_struct), allocatable :: stk(:)
type (tao_eval_stack1_struct), allocatable, optional :: stack(:)
type (ele_struct), optional, pointer :: dflt_ele_ref, dflt_ele_start, dflt_ele
type (coord_struct), optional :: dflt_orbit
type (tao_expression_info_struct), allocatable, optional :: info(:)
type (tao_data_struct), optional :: datum
type (expression_func_struct) func(0:20)

integer, optional :: dflt_uni, dflt_eval_point
integer, allocatable :: op(:)
integer i_lev, i_op, i, ios, n, n_size, ix0, ix1, ix2, ix3, ix4, ix5, n_func
integer ix_word, i_delim, i2, ix, ix_word2

real(rp), allocatable :: value(:)
real(rp), optional :: dflt_s_offset

character(*) :: expression
character(*), optional :: dflt_component, dflt_source
character(*), optional :: dflt_dat_or_var_index

character(len(expression)+20) :: phrase, word, word2
character(1) delim, cc
character(80) default_source
character(40) saved_prefix
character(*), parameter :: r_name = "tao_evaluate_expression"

logical delim_found, do_combine, use_good_user, in_species_func
logical err_flag, err, wild, printit, found
logical, optional :: print_err

! Don't destroy the input expression

err_flag = .true.
saved_prefix = ''
printit = logic_option(.true., print_err)
default_source = ''
if (present(dflt_source)) default_source = dflt_source

phrase = expression
if (present(dflt_ele)) then
   if (associated(dflt_ele)) call tao_expression_hash_substitute(phrase, ele_full_name(dflt_ele, '!#'))
end if

if (len(phrase) > 11) then
  if (phrase(1:11) == 'expression:') phrase = phrase(12:)
endif

! if phrase is blank then return 0.0

call string_trim (phrase, phrase, ix)
if (ix == 0) then
  call out_io (s_warn$, r_name, "Expression is blank")
  call re_allocate (value, max(1, n_size))
  value = 0.0
  return
endif
 
! General idea: Create a reverse polish stack that represents the expression.
! Reverse polish means that the operand goes last so that 2 * 3 is written 
! on the stack as: [2, 3, *]

! The stack is called: stk
! Since operations move towards the end of the stack we need a separate
! stack called op which keeps track of what operations have not yet
! been put on stk.

! init

n_func = 0
i_lev = 0
i_op = 0
in_species_func = .false.

allocate (stk(20), op(20))

! parsing loop to build up the stack.

parsing_loop: do

  ! get a word

  if (in_species_func) then
    call word_read (phrase, ')', word, ix_word, delim, delim_found, phrase)
  else
    call word_read (phrase, '+-*/()^,}[ ', word, ix_word, delim, delim_found, phrase)
  endif

  ! Args are counted counted at the beginning of the function and at each comma.

  if (n_func > 0 .and. (ix_word /= 0 .or. delim /= ')')) then
    if (func(n_func)%n_arg_count == 0) func(n_func)%n_arg_count = 1
  endif

  !--------------------------
  ! Preliminary: If we have split up something that should have not been split
  ! then put it back together again...

  ! just make sure we are not chopping a number in two, e.g. "3.5e-7" should not
  ! get split at the "-" even though "-" is a delimiter. Also "Q01+4[k1]" should not be split.

  do_combine = (delim == '-' .or. delim == '+') 
  if (do_combine .and. ix_word == 0) do_combine = .false.

  if (do_combine) then
    found = .false.   ! Found "NNN[" like construct where NNN is an integer?
    ix = index(phrase, '[')
    if (ix /= 0) then
      if (is_integer(phrase(:ix-1))) found = .true.
    endif

    if (.not. found) then ! Test if a number
      cc = upcase(word(ix_word:ix_word))
      if (cc == 'E' .or. cc == 'D') then
        do i = 1, ix_word-1
          if (index('.0123456789', word(i:i)) == 0) do_combine = .false.
        enddo
      else
        do_combine = .false.
      endif
    endif
  endif

  ! If still DO_COMBINE = True then we need to unsplit

  if (do_combine) then
    word = trim(word) // delim
    call word_read (phrase, '+-*/()^,}[', word2, ix_word2, delim, delim_found, phrase)
    word = trim(word) // word2
    ix_word = len_trim(word)
  endif

  ! Something like "lcav[lr(2).freq]" or "[2,4]@orbit.x[1,4] or "[ele::q20w[hkick], ele::q20w[hkick]]"
  ! will get split on the "["

  do
    if (delim /= '[') exit

    call word_read (phrase, ']', word2, ix_word2, delim, delim_found, phrase, ignore_interior = .true.)
    if (.not. delim_found) then
      call out_io (s_error$, r_name, "NO MATCHING ']' FOR OPENING '[':" // expression)
      return
    endif
    word = trim(word) // '[' // trim(word2) // ']'
    ix_word = len_trim(word)
    if (phrase == ' ') then  
      delim_found = .false.
      delim = ' '
    elseif (phrase(1:1) == ' ') then  
      call string_trim (phrase, phrase, ix)
      if (index('+-*/()^,}[', phrase(1:1)) == 0) then
        delim = ' '
      else
        delim = phrase(1:1)
        phrase = phrase(2:)
      endif
    else          ! even more...
      call word_read (phrase, '[+-*/()^,}', word2, ix_word2, delim, delim_found, phrase)
      word = trim(word) // trim(word2)       
      ix_word = len_trim(word)
    endif
  enddo

  ! If delim = "*" then see if this is being used as a wildcard
  ! Examples: "[*]|", "*.*|", "*.x|", "*@orbit.x|", "*@*|", "orbit.*[3]|", "ele::q*1[beta_a]", "var::*d|model"
  ! If so, we have split in the wrong place and we need to correct this. 
  ! Something like "emit*data::mf.xm|model", "3*[1,2]" or "3*.42" does not get split.

  do
    if (delim /= '*') exit

    ix0 = index(word, '::')
    ix4 = index(word, '|')

    ix1 = index(phrase, '[')
    ix2 = index(phrase, ']')
    ix3 = index(phrase, '|')
    ix5 = index(phrase, '::')
    if (ix5 > 0) then   ! "*" in "emit*data::mf.xm|model" is not wild
      if (ix1 > ix5) ix1 = 0   ! Ignore if after
      if (ix2 > ix5) ix2 = 0
      if (ix3 > ix5) ix3 = 0
    endif

    if (phrase(1:1) == '[' .and. ((ix0 == 0) .eqv. (ix4 == 0))) exit

    ! If in "[...*...]" construct is wild
    wild = .false.
    if (ix2 /= 0 .and. (ix1 == 0 .or. ix1 > ix2)) wild = .true.
    if (ix3 /= 0 .and. (ix1 == 0 .or. ix1 > ix3)) wild = .true.

    if (.not. wild) then
      select case (phrase(1:1))
      case ( ']', '[', '|', '@')
        wild = .true.
      case ('.')
        wild = (index('0123456789', phrase(2:2)) == 0) ! Wild if not a number
      case default
        ! If in "::xxx*yyy[" construct where each x and y is not one of ":", " ", etc.
        found = .false.
        if (ix0 /= 0 .and. ix1 /= 0) then
          do ix = ix0+2, ix_word
            if (index(': ]|', word(ix:ix)) /= 0) found = .true.
          enddo
          do ix = 1, ix1-1
            if (index(': ]|', phrase(ix:ix)) /= 0) found = .true.
          enddo
          if (.not. found) wild = .true.
        endif
      end select
    endif

    if (.not. wild) exit

    word = word(:ix_word) // '*'
    call word_read (phrase, '+-*/()^,}', word2, ix_word2, delim, delim_found, phrase, .true.)
    word = trim(word) // trim(word2)       
    ix_word = len_trim(word)
  enddo

  !---------------------------
  ! Now see what we got...

  ! For a word ending in '|' then must be a construct like 'orbit.x|-model'.
  ! So store the 'orbit.x|' prefix

  if (ix_word /= 0) then
    if (word(ix_word:ix_word) == '|') then
      saved_prefix = word
      word = ''
      ix_word = 0
    endif
  endif

  ! if the word is a datum without an "|", and dflt_component is present, 
  ! Then use the dflt_component.

  if (present(dflt_component) .and. index(word, '|') == 0) then
    call tao_find_data (err, word, print_err = .false.)
    if (.not. err) then
      phrase = trim(word) // '|' // trim(dflt_component) // delim // trim(phrase)
      cycle   ! Try again
    endif
  endif

  ! For a "(" delim we must have a function

  if (delim == '(') then

    if (ix_word /= 0) then
      word2 = word
      call downcase_string (word2)
      n_func = n_func + 1
      func(n_func) = expression_func_struct(word2, 1, 0)
      select case (word2)
      case ('cot');             call push_op_stack (op, i_op, cot$)
      case ('csc');             call push_op_stack (op, i_op, csc$)
      case ('sec');             call push_op_stack (op, i_op, sec$)
      case ('sin');             call push_op_stack (op, i_op, sin$)
      case ('sinc');            call push_op_stack (op, i_op, sinc$)
      case ('cos');             call push_op_stack (op, i_op, cos$)
      case ('tan');             call push_op_stack (op, i_op, tan$)
      case ('asin');            call push_op_stack (op, i_op, asin$)
      case ('acos');            call push_op_stack (op, i_op, acos$)
      case ('atan');            call push_op_stack (op, i_op, atan$)
      case ('atan2')
        call push_op_stack (op, i_op, atan2$)
        func(n_func)%n_arg_target = 2
      case ('modulo')
        call push_op_stack (op, i_op, modulo$)
        func(n_func)%n_arg_target = 2
      case ('sinh');            call push_op_stack (op, i_op, sinh$)
      case ('cosh');            call push_op_stack (op, i_op, cosh$)
      case ('tanh');            call push_op_stack (op, i_op, tanh$)
      case ('coth');            call push_op_stack (op, i_op, coth$)
      case ('asinh');           call push_op_stack (op, i_op, asinh$)
      case ('acosh');           call push_op_stack (op, i_op, acosh$)
      case ('atanh');           call push_op_stack (op, i_op, atanh$)
      case ('acoth');           call push_op_stack (op, i_op, acoth$)
      case ('abs');             call push_op_stack (op, i_op, abs$)
      case ('min');             call push_op_stack (op, i_op, min$)
      case ('max');             call push_op_stack (op, i_op, max$)
      case ('rms');             call push_op_stack (op, i_op, rms$)
      case ('average', 'mean'); call push_op_stack (op, i_op, average$)
      case ('sum');             call push_op_stack (op, i_op, sum$)
      case ('sqrt');            call push_op_stack (op, i_op, sqrt$)
      case ('log');             call push_op_stack (op, i_op, log$)
      case ('exp');             call push_op_stack (op, i_op, exp$)
      case ('factorial');       call push_op_stack (op, i_op, factorial$)
      case ('ran')         
        call push_op_stack (op, i_op, ran$)
        func(n_func)%n_arg_target = 0
      case ('ran_gauss')
        call push_op_stack (op, i_op, ran_gauss$)
        func(n_func)%n_arg_target = -1      ! 0 or 1 args
      case ('int');             call push_op_stack (op, i_op, int$)
      case ('sign');            call push_op_stack (op, i_op, sign$)
      case ('nint');            call push_op_stack (op, i_op, nint$)
      case ('floor');           call push_op_stack (op, i_op, floor$)
      case ('ceiling');         call push_op_stack (op, i_op, ceiling$)
      case ('mass_of');         call push_op_stack (op, i_op, mass_of$)
      case ('charge_of');       call push_op_stack (op, i_op, charge_of$)
      case ('anomalous_moment_of'); call push_op_stack (op, i_op, anomalous_moment_of$)
      case ('species');         call push_op_stack (op, i_op, species$)
      case ('antiparticle');    call push_op_stack (op, i_op, antiparticle$)
      case default
        call out_io (s_error$, r_name, 'UNEXPECTED CHARACTERS (BAD FUNCTION NAME?) BEFORE "(": ', 'IN EXPRESSION: ' // expression)
        return
      end select

      call push_op_stack (op, i_op, l_func_parens$)

      ! Parse function argument for functions that take a species.
      select case (word2)
      case ('mass_of', 'charge_of', 'species', 'antiparticle', 'anomalous_moment_of')
        in_species_func = .true.
      end select

    else
      call push_op_stack (op, i_op, l_parens$)
    endif

    cycle parsing_loop

  ! for a unary "-"

  elseif (delim == '-' .and. ix_word == 0) then
    call push_op_stack (op, i_op, unary_minus$)
    cycle parsing_loop

  ! for a unary "+"

  elseif (delim == '+' .and. ix_word == 0) then
    call push_op_stack (op, i_op, unary_plus$)
    cycle parsing_loop

  ! for a ")" delim

  elseif (delim == ')') then
    if (ix_word == 0) then
      if (n_func == 0 .or. (func(n_func)%n_arg_target /= 0 .and. func(n_func)%n_arg_target /= -1)) then
        if (printit) call out_io (s_error$, r_name, 'CONSTANT OR VARIABLE MISSING BEFORE ")"', &
                                                    'IN EXPRESSION: ' // expression)
        return
      endif

    else
      if (in_species_func) then
        call push_stack (stk, i_lev, species_const$)
        stk(i_lev)%name = word
      else
        call push_stack (stk, i_lev, numeric$)
        call tao_param_value_routine (word, use_good_user, saved_prefix, stk(i_lev), err, printit, &
               dflt_component, default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, &
               dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)
        if (err) then
          if (printit) call out_io (s_error$, r_name, &
                          'ERROR IN EVALUATING EXPRESSION: ' // expression, &
                          'CANNOT EVALUATE: ' // word)
          return
        endif
      endif
    endif

    do
      do i = i_op, 1, -1       ! release pending ops
        if (op(i) == l_parens$ .or. op(i) == l_func_parens$) exit            ! break do loop
        call push_stack (stk, i_lev, op(i))
      enddo

      if (i == 0) then
        call out_io (s_error$, r_name, 'UNMATCHED ")" IN EXPRESSION: ' // expression)
        return
      endif

      i_op = i - 1

      if (op(i) == l_func_parens$) then
        if (func(n_func)%n_arg_target == -1) then
          if (func(n_func)%n_arg_count /= 0 .and. func(n_func)%n_arg_count /= 1) then
            if (printit) call out_io(s_error$, r_name, &
                            'FUNCTION: ' // trim(func(n_func)%name) // ' DOES NOT HAVE 0 OR 1 ARGUMENTS.', &
                            'IN EXPRESSION: ' // expression)
            return
          endif
          call push_stack (stk, i_lev, arg_count$)
          call re_allocate (stk(i_lev)%value, 1)
          stk(i_lev)%value(1) = func(n_func)%n_arg_count

        else
          if (func(n_func)%n_arg_count /= func(n_func)%n_arg_target) then
            if (printit) call out_io(s_error$, r_name, &
                          'FUNCTION: ' // trim(func(n_func)%name) // ' DOES NOT HAVE THE CORRECT NUMBER OF ARGUMENTS.', &
                          'IN EXPRESSION: ' // expression)
            return
          endif
        endif

        n_func = n_func - 1
      endif

      call word_read (phrase, '+-*/()^,}', word, ix_word, delim, delim_found, phrase)
      if (ix_word /= 0) then
        if (printit) call out_io (s_error$, r_name, 'UNEXPECTED CHARACTERS AFTER ")" IN EXPRESSION: ' // expression)
        return
      endif

      if (delim /= ')') exit    ! if no more ')' then no need to release more
    enddo


    if (delim == '(') then
      if (printit) call out_io (s_error$, r_name, '")(" CONSTRUCT DOES NOT MAKE SENSE IN EXPRESSION: ' // expression)
      return
    endif

  ! For binary "+-/*^" delims

  else
    if (ix_word == 0) then
      call out_io (s_error$, r_name, 'CONSTANT OR VARIABLE MISSING IN EXPRESSION: ' // expression)
      return
    endif
    call push_stack (stk, i_lev, numeric$)
    call tao_param_value_routine (word, use_good_user, saved_prefix, stk(i_lev), err, printit, &
            dflt_component, default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, &
            dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)
    if (err) then
      if (printit) call out_io (s_error$, r_name, &
                        'ERROR IN EXPRESSION: ' // expression, &
                        'CANNOT EVALUATE: ' // word)
      return
    endif
  endif

  ! If we are here then we have an operation that is waiting to be identified

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
  case (',')
    i_delim = comma$
    func(n_func)%n_arg_count = func(n_func)%n_arg_count + 1
  case ('}')
    i_delim = no_delim$
    call out_io (s_error$, r_name, &
                      'DELIMITOR FOUND OUT OF PLACE: ' // delim, &
                      'IN EXPRESSION: ' // expression)
    return
  case default
    if (delim_found) then
      if (delim == ' ') then
        if (printit) call out_io (s_error$, r_name, 'MALFORMED EXPRESSION: ' // expression)
        return
      endif

      call out_io (s_error$, r_name, 'INTERNAL ERROR')
      call err_exit
    endif
    i_delim = no_delim$
  end select

  ! Now see if there are operations on the OP stack that need to be transferred
  ! to the STK stack

  do i = i_op, 1, -1
    if (expression_eval_level(op(i)) < expression_eval_level(i_delim)) exit

    if (op(i) == l_parens$) then
      if (printit) call out_io (s_error$, r_name, 'UNMATCHED "(" IN EXPRESSION: ' // expression)
      return
    endif

    if (op(i) == l_func_parens$) then
      if (i_delim /= comma$) then
        if (printit) call out_io (s_error$, r_name, 'UNMATCHED "("', 'IN EXPRESSION: ' // expression)
        return
      endif
      i_op = i
      cycle parsing_loop
    endif

    call push_stack (stk, i_lev, op(i))
    in_species_func = .false.
  enddo

  ! put the pending operation on the OP stack

  i_op = i
  select case (i_delim)
  case (no_delim$); exit parsing_loop
  case (comma$)
    if (printit) call out_io (s_error$, r_name, 'COMMA AT END OF EXPRESSION IS OUT OF place: ' // expression, &
                                   '(NEEDS "[...]" BRACKETS IF AN ARRAY.)')
    return
  case default; call push_op_stack (op, i_op, i_delim)
  end select

enddo parsing_loop

!------------------------------------------------------------------
! Some error checks

if (i_op /= 0) then
  call out_io (s_error$, r_name, 'UNMATCHED "(" IN EXPRESSION: ' // expression)
  return
endif

if (i_lev == 0) then
  call out_io (s_error$, r_name, 'NO VALUE FOUND IN EXPRESSION: ' // expression)
  return
endif

if (phrase /= '') then
  call out_io (s_error$, r_name, 'EXTRA STUFF AFTER EXPRESSION: ' // phrase)
  return
endif

call tao_evaluate_stack (stk(1:i_lev), n_size, use_good_user, value, err_flag, printit, expression, info)

! If the stack argument is present then copy stk to stack

if (present(stack)) then
  if (allocated(stack)) deallocate(stack)
  allocate (stack(i_lev))
  do i = 1, i_lev
    if (allocated (stk(i)%value)) then
      n = size(stk(i)%value)
      allocate (stack(i)%value(n), stack(i)%info(n))
      if (allocated (stack(i)%value_ptr)) allocate (stack(i)%value_ptr(n))
    endif
    stack(i) = stk(i)
  enddo
endif

!-------------------------------------------------------------------------
! The op_stack is for operators and functions.

contains

subroutine push_op_stack (op_stack, i_lev, this_type)

integer, allocatable :: op_stack(:)
integer i_lev, this_type

character(*), parameter :: r_name = "push_op_stack"

!

i_lev = i_lev + 1
if (i_lev > size(op_stack)) call re_allocate(op_stack, 2*i_lev)
op_stack(i_lev) = this_type

end subroutine push_op_stack

!-------------------------------------------------------------------------
! contains

subroutine push_stack (stack, i_lev, this_type)

type (tao_eval_stack1_struct), allocatable :: stack(:), tmp_stk(:)
integer i_lev, this_type

character(*), parameter :: r_name = "push_stack"

!

i_lev = i_lev + 1

if (i_lev > size(stack)) then
  call move_alloc(stack, tmp_stk)
  allocate(stack(2*i_lev))
  stack(1:i_lev-1) = tmp_stk
endif

stack(i_lev)%type = this_type
stack(i_lev)%name = expression_op_name(this_type)
stack(i_lev)%scale = 1

end subroutine push_stack
                       
end subroutine tao_evaluate_expression
