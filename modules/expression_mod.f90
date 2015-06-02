module expression_mod

use bmad_struct

implicit none

integer, parameter :: end_stack$ = 0, plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
integer, parameter :: l_parens$ = 5, r_parens$ = 6, power$ = 7
integer, parameter :: unary_minus$ = 8, unary_plus$ = 9, no_delim$ = 10
integer, parameter :: sin$ = 11, cos$ = 12, tan$ = 13
integer, parameter :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17, sqrt$ = 18
integer, parameter :: log$ = 19, exp$ = 20, ran$ = 21, ran_gauss$ = 22
integer, parameter :: atan2$ = 23, factorial$ = 24, int$ = 25, nint$ = 26
integer, parameter :: floor$ = 27, ceiling$ = 28, numeric$ = 29

character(12), parameter :: expression_op_name(28) = [character(12) :: '+', '-', '*', '/', &
                                    '(', ')', '^', '-', '+', '', 'sin', 'cos', 'tan', &
                                    'asin', 'acos', 'atan', 'abs', 'sqrt', 'log', 'exp', 'ran', &
                                    'ran_gauss', 'atan2', 'factorial', 'int', 'nint', 'floor', 'ceiling']


integer, parameter :: expression_eval_level(29) = [1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
                            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9]

private pushit

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine expression_string_to_stack (string, stack, n_stack, err_flag, err_str)
!
! This routine creates an expression stack array which can be used 
! to evaluate an arithmethic expression.
! The end elements of stack not used are marked stack(i)%type = end_stack$
!
! Input:
!   string    -- character(*): Expression to be converted.
!
! Output:
!   stack(:)  -- expression_atom_struct: Expression evaluation stack.
!   n_stack   -- integer: number of "atoms" used by the expression
!   err_flag  -- logical: Set True if there is an error (EG divide by 0).
!   err_str   -- character(*): String describing the error.
!-

subroutine expression_string_to_stack (string, stack, n_stack, err_flag, err_str)

type (expression_atom_struct) :: stack(:)

integer i_op, i

integer op(200), ix_word, i_delim, i2, ix_word2, n_stack

real(rp) value

character(*) string, err_str
character(1) delim
character(80) word, word2
character(200) parse_line


logical delim_found, split, zero_arg_function_pending
logical err_flag

! The general idea is to rewrite the expression on a stack in reverse polish.
! Reverse polish means that the operand goes last so that 2 * 3 is writen 
! on the stack as: [2, 3, *]

! Since operations move towards the end of the stack we need a separate
! stack called op which keeps track of what operations have not yet
! been put on stack.

! init

err_flag = .true.
n_stack = 0
i_op = 0
zero_arg_function_pending = .false.
stack(:)%type = end_stack$
call str_upcase (parse_line, string)

! parsing loop to build up the stack.

parsing_loop: do

  ! get a word

  call get_next_chunk (word, ix_word, '+-*/()^,:} ', delim, delim_found)

  if (delim == '*' .and. word(1:1) == '*') then
    err_str = 'EXPONENTIATION SYMBOL IS "^" AS OPPOSED TO "**"!'
    return
  endif

  if (zero_arg_function_pending .and. (ix_word /= 0 .or. delim /= ')')) then
    err_str = 'RAN AND RAN_GAUSS DO NOT TAKE AN ARGUMENT'
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
    call get_next_chunk (word2, ix_word2, '+-*/()^,:}', delim, delim_found)
    word = word(:ix_word+1) // word2
    ix_word = ix_word + ix_word2
  endif

  ! Something like "lcav[lr(2).freq]" will get split on the "("

  if (delim == '(' .and. index(word, '[LR') /= 0) then
    call get_next_chunk (word2, ix_word2, '+-*/(^,:}', delim, delim_found)
    word = word(:ix_word) // '(' // word2
    ix_word = ix_word + ix_word2 + 1
  endif

  !---------------------------
  ! Now see what we got...

  ! For a "(" delim we must have a function

  if (delim == '(') then

    zero_arg_function_pending = .false.
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
      case ('ATAN2') 
        call pushit (op, i_op, atan2$)
      case ('ABS') 
        call pushit (op, i_op, abs$)
      case ('SQRT') 
        call pushit (op, i_op, sqrt$)
      case ('LOG') 
        call pushit (op, i_op, log$)
      case ('EXP') 
        call pushit (op, i_op, exp$)
      case ('FACTORIAL') 
        call pushit (op, i_op, factorial$)
      case ('RAN') 
        call pushit (op, i_op, ran$)
        zero_arg_function_pending = .true.
      case ('RAN_GAUSS') 
        call pushit (op, i_op, ran_gauss$)
        zero_arg_function_pending = .true.
      case ('INT')
        call pushit (op, i_op, int$)
      case ('NINT')
        call pushit (op, i_op, nint$)
      case ('FLOOR')
        call pushit (op, i_op, floor$)
      case ('CEILING')
        call pushit (op, i_op, ceiling$)
      case default
        err_str = 'UNEXPECTED CHARACTERS ON RHS BEFORE "(": ' // word
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
      if (.not. zero_arg_function_pending) then
        err_str = 'CONSTANT OR VARIABLE MISSING BEFORE ")"'
        return
      endif
      zero_arg_function_pending = .false.
    else
      call pushit (stack%type, n_stack, numeric$)
      stack(n_stack)%name = word
      
    endif

    do
      do i = i_op, 1, -1     ! release pending ops
        if (op(i) == l_parens$) exit          ! break do loop
        call pushit (stack%type, n_stack, op(i))
      enddo

      if (i == 0) then
        err_str = 'UNMATCHED ")" ON RHS'
        return
      endif

      i_op = i - 1

      call get_next_chunk (word, ix_word, '+-*/()^,:}', delim, delim_found)
      if (ix_word /= 0) then
        err_str = 'UNEXPECTED CHARACTERS ON RHS AFTER ")"'
        return
      endif

      if (delim /= ')') exit  ! if no more ')' then no need to release more
    enddo


    if (delim == '(') then
      err_str = '")(" CONSTRUCT DOES NOT MAKE SENSE'
      return
    endif

  ! For binary "+-/*^" delims

  else
    if (ix_word == 0) then
      err_str = 'CONSTANT OR VARIABLE MISSING'
      return
    endif
    call pushit (stack%type, n_stack, numeric$)
    stack(n_stack)%name = word
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
  case ('^')
    i_delim = power$
  case (',', '}', ':', ')')   ! End of expression delims
    i_delim = no_delim$
  case default
    err_str = 'MALFORMED EXPRESSION'
    return
  end select

  ! now see if there are operations on the OP stack that need to be transferred
  ! to the STACK stack

  do i = i_op, 1, -1
    if (expression_eval_level(op(i)) >= expression_eval_level(i_delim)) then
      if (op(i) == l_parens$) then
        if (i > 1 .and. op(max(1,i-1)) == atan2$ .and. delim == ',') cycle parsing_loop
        err_str = 'UNMATCHED "("'
        return
      endif
      call pushit (stack%type, n_stack, op(i))
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
! Go through the stack and perform the operations

if (i_op /= 0) then
  err_str = 'UNMATCHED "("'
  return
endif

if (n_stack == 0) then
  err_str = 'NO VALUE FOUND'
  return
endif

err_flag = .false.

!-------------------------------------------------------------------------
contains

subroutine get_next_chunk (word, ix_word, delim_list, delim, delim_found)

character(*) word, delim_list, delim

integer ix_word

logical delim_found

!

call word_read (parse_line, delim_list, word, ix_word, delim, delim_found, parse_line)

end subroutine get_next_chunk

end subroutine expression_string_to_stack

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine pushit (stack, i_stk, value)
!
! Private routine used by expression_string_to_stack
!-

subroutine pushit (stack, i_stk, value)

implicit none

integer stack(:), i_stk, value
character(*), parameter :: r_name = 'pushit'

!

i_stk = i_stk + 1

if (i_stk > size(stack)) then
  call out_io (s_fatal$, r_name, 'STACK OVERFLOW, EXPERT HELP IS NEEDED!')
  if (global_com%exit_on_error) call err_exit
  return
endif

stack(i_stk) = value

end subroutine pushit
                     
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine evaluate_expression_stack (stack, value, err_flag, err_str, var, use_old)
!
! Routine to evaluate an mathematical expression represented by an "expression stack".
! Expression stacks are created by bmad_parser code.
!
! Input:
!   stack(:)    -- expression_atom_struct: Expression to evaluate.
!   var(:)      -- controller_var_struct, optional: Array of variables.
!   use_old     -- logical, optional: Use var%old_value? Must be present if var(:) is present.
!
! Output:
!   value       -- real(rp): Value of the expression.
!   err_flag    -- logical: True if there is an evaluation problem. False otherwise.
!   err_str     -- character(*): Error string explaining error if there is one.
!-

subroutine evaluate_expression_stack (stack, value, err_flag, err_str, var, use_old)

use random_mod

type (expression_atom_struct) stack(:)
type (expression_atom_struct) stack2(20)
type (controller_var_struct), optional :: var(:)

real(rp) value
integer i, i2, ix

logical, optional :: use_old
logical err_flag

character(*) err_str

!

err_flag = .true.

i2 = 0
do i = 1, size(stack)

  select case (stack(i)%type)

  case (end_stack$)
    exit

  case (numeric$)
    i2 = i2 + 1
    stack2(i2)%value = stack(i)%value

  case (unary_minus$)
    stack2(i2)%value = -stack2(i2)%value

  case (unary_plus$)
    stack2(i2)%value = stack2(i2)%value

  case (plus$)
    stack2(i2-1)%value = stack2(i2-1)%value + stack2(i2)%value
    i2 = i2 - 1

  case (minus$)
    stack2(i2-1)%value = stack2(i2-1)%value - stack2(i2)%value
    i2 = i2 - 1

  case (times$)
    stack2(i2-1)%value = stack2(i2-1)%value * stack2(i2)%value
    i2 = i2 - 1

  case (divide$)
    if (stack2(i2)%value == 0) then
      err_str = 'DIVIDE BY 0 ON RHS'
      return
    endif
    stack2(i2-1)%value= stack2(i2-1)%value / stack2(i2)%value
    i2 = i2 - 1

  case (power$)
    stack2(i2-1)%value = stack2(i2-1)%value**stack2(i2)%value
    i2 = i2 - 1

  case (sin$)
    stack2(i2)%value = sin(stack2(i2)%value)

  case (cos$)
    stack2(i2)%value = cos(stack2(i2)%value)

  case (tan$)
    stack2(i2)%value = tan(stack2(i2)%value)

  case (asin$)
    if (stack2(i2)%value < -1 .or. stack2(i2)%value > 1) then
      err_str = 'ASIN ARGUMENT HAS MAGNITUDE GREATER THAN 1'
      return
    endif
    stack2(i2)%value = asin(stack2(i2)%value)

  case (acos$)
    if (stack2(i2)%value < -1 .or. stack2(i2)%value > 1) then
      err_str = 'ACOS ARGUMENT HAS MAGNITUDE GREATER THAN 1'
      return
    endif
    stack2(i2)%value = acos(stack2(i2)%value)

  case (factorial$)
    stack2(i2)%value = factorial(nint(stack2(i2)%value))
    if (stack2(i2)%value < 0) then
      err_str = 'FACTORIAL PROBLEM'
      return
    endif

  case (atan$)
    stack2(i2)%value = atan(stack2(i2)%value)

  case (atan2$)
    stack2(i2-1)%value = atan2(stack2(i2-1)%value, stack2(i2)%value)
    i2 = i2 - 1

  case (abs$)
    stack2(i2)%value = abs(stack2(i2)%value)

  case (sqrt$)
    if (stack2(i2)%value < 0) then
      err_str = 'SQRT ARGUMENT IS NEGATIVE '
      return
    endif
    stack2(i2)%value = sqrt(stack2(i2)%value)

  case (log$)
    if (stack2(i2)%value < 0) then
      err_str = 'LOG ARGUMENT IS NEGATIVE '
      return
    endif
    stack2(i2)%value = log(stack2(i2)%value)

  case (exp$)
    stack2(i2)%value = exp(stack2(i2)%value)

  case (int$)
    stack2(i2)%value = int(stack2(i2)%value)

  case (nint$)
    stack2(i2)%value = nint(stack2(i2)%value)

  case (floor$)
    stack2(i2)%value = floor(stack2(i2)%value)

  case (ceiling$)
    stack2(i2)%value = ceiling(stack2(i2)%value)

  case (ran$)
    i2 = i2 + 1
    call ran_uniform(stack2(i2)%value)

  case (ran_gauss$)
    i2 = i2 + 1
    call ran_gauss(stack2(i2)%value)

  case default

    if (is_attribute(stack(i)%type, present_var$)) then
      if (.not. present(var)) then
        err_flag = .true.
        err_str = 'VAR ARGUMENT NOT PRESENT! GET HELP!'
        return
      endif
      ix = stack(i)%type - var_offset$
      i2 = i2 + 1
      if (logic_option(.false., use_old)) then
        stack2(i2)%value = var(ix)%old_value
      else
        stack2(i2)%value = var(ix)%value
      endif

    else
      err_str = 'INTERNAL ERROR #02: GET HELP'
      if (global_com%exit_on_error) call err_exit
    endif

  end select
enddo

if (i2 /= 1) then
  err_str = 'INTERNAL ERROR #03: GET HELP'
  if (global_com%exit_on_error) call err_exit
endif

value = stack2(1)%value
err_flag = .false.

end subroutine evaluate_expression_stack

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function expression_stack_to_string (stack, polish) result (str)
!
! Routine to convert an expression stack to a string
!
! Input:
!   stack(:)  -- expression_atom_struct: arithmetic expression
!   polish    -- logical, optional, Construct expression in reverse polish? Default is False.
!
! Output:
!   str       -- character(200): Expression in string form.
!-

function expression_stack_to_string (stack, polish) result (str)

type (expression_atom_struct), target :: stack(:)
type (expression_atom_struct), pointer :: atom, atom2
type (expression_atom_struct) s2(20)

character(200) str
character(200) s2_name(20)

integer i, i2, ix
logical, optional :: polish

!

str = ''

! Polish notation

if (logic_option(.false., polish)) then
  do i = 1, size(stack)
    atom => stack(i)
    if (atom%type == end_stack$) exit
    if (atom%type >= 1 .and. atom%type <= size(expression_op_name)) then
      str = trim(str) // ', ' // expression_op_name(atom%type)    
    else
      str = trim(str) // ', ' // atom%name
    endif
  enddo

  str = str(3:)

! Standard notation

else

  i2 = 0
  i = 0

  do
    if (i+1 > size(stack)) exit
    if (stack(i+1)%type == end_stack$) exit

    i = i + 1
    atom => stack(i)

    if (is_attribute(atom%type, all_var$)) then
      i2 = i2 + 1
      s2(i2)%type = numeric$

      if (atom%name == '') then
        s2_name(i2) = trim(real_to_string(atom%value))
      else
        s2_name(i2) = trim(atom%name)
      endif
      cycle
    endif

    select case (atom%type)
    case (plus$, minus$, times$, divide$, power$)
      if (expression_eval_level(s2(i2-1)%type) < expression_eval_level(atom%type)) s2(i2-1)%name = '(' // trim(s2_name(i2-1)) // ')'
      if (atom%type == minus$ .or. atom%type == divide$) then
        if (expression_eval_level(s2(i2)%type) <= expression_eval_level(atom%type)) s2_name(i2) = '(' // trim(s2_name(i2)) // ')'
      else
        if (expression_eval_level(s2(i2)%type) < expression_eval_level(atom%type)) s2_name(i2) = '(' // trim(s2_name(i2)) // ')'
      endif
      s2_name(i2-1) = trim(s2_name(i2-1)) // trim(expression_op_name(atom%type)) // s2_name(i2)
      i2 = i2 - 1
      s2%type = atom%type

    case (numeric$)
      i2 = i2 + 1
      s2(i2)%type = atom%type
      if (atom%name == '') then
        s2_name(i2) = real_to_string(atom%value)
      else
        s2_name(i2) = atom%name
      endif

    case (unary_minus$, unary_plus$)
      if (expression_eval_level(s2(i2)%type) <= expression_eval_level(atom%type)) s2_name(i2) = '(' // trim(s2_name(i2)) // ')'
      s2_name(i2) = '-' // s2_name(i2)
 
    case (ran$, ran_gauss$)
      i2 = i2 + 1
      s2_name(i2) = trim(expression_op_name(atom%type)) // '()'
      s2%type = atom%type

    case (atan2$)
      s2_name(i2-1) = trim(expression_op_name(atom%type)) // '(' // trim(s2_name(i2-1)) // ',' // trim(s2_name(i2)) // ')'
      i2 = i2 - 1
      s2%type = atom%type

    case (factorial$)
      if (expression_eval_level(s2(i2)%type) <= expression_eval_level(atom%type)) s2_name(i2) = '(' // trim(s2_name(i2)) // ')'
      s2_name(i2) = trim(s2_name(i2)) // '!'
      s2%type = atom%type

    case default ! Function
      s2_name(i2) = trim(expression_op_name(atom%type)) // '(' // trim(s2_name(i2)) // ')'
      s2%type = atom%type

    end select

  enddo

  str = s2_name(i2)

endif

end function expression_stack_to_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function linear_coef (stack, err_flag) result (coef)
!
! Routine to return the linear coefficient of a linear expression.
!
! Input:
!   stack(:) -- expression_atom_struct: Expression stack.
!
! Output:
!   err_flag  -- Logical: Set True if the expression is not linear
!   coef      -- real(rp): Linear coefficient.
!-

function linear_coef (stack, err_flag) result (coef)

type (expression_atom_struct) stack(:)

real(rp) coef
logical, optional :: err_flag

integer i, i0, i1, n

character(40) err_str

!

n = size(stack)
err_flag = .true.
coef = 0

! If expression = "Var" then coef = 1.

if (n == 1 .and. stack(1)%type > var_offset$) then
  err_flag = .false.
  coef = 1
  return
endif

! Expression must have times at the end to be linear.

if (stack(n)%type /= times$) return

! To be linear, stack(1) or stack(n-1) must be variable and have no other variables.

if (stack(1)%type > var_offset$) then
  i0 = 2
  i1 = n-1
elseif (stack(n-1)%type > var_offset$) then
  i0 = 1
  i1 = n-2
else
  return
endif

call evaluate_expression_stack(stack(i0:i1), coef, err_flag, err_str)

end function linear_coef

end module
