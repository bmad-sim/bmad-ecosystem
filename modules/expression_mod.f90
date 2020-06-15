module expression_mod

use bmad_struct

implicit none

! The numeric$ category is for numeric constants [EG: "1.3d-5"].
! The variable$ category includes symbolic constants defined in a lattice file, lattice parameters, etc.
! The species$ category is for the species() function. 
! The species_const$ category is for particle species ('He3', etc).

integer, parameter :: end_stack$ = 0, plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
integer, parameter :: l_parens$ = 5, r_parens$ = 6, power$ = 7
integer, parameter :: unary_minus$ = 8, unary_plus$ = 9, no_delim$ = 10
integer, parameter :: sin$ = 11, cos$ = 12, tan$ = 13
integer, parameter :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17, sqrt$ = 18
integer, parameter :: log$ = 19, exp$ = 20, ran$ = 21, ran_gauss$ = 22, atan2$ = 23
integer, parameter :: factorial$ = 24, int$ = 25, nint$ = 26, floor$ = 27, ceiling$ = 28
integer, parameter :: numeric$ = 29, variable$ = 30
integer, parameter :: mass_of$ = 31, charge_of$ = 32, anomalous_moment_of$ = 33, species$ = 34, species_const$ = 35
integer, parameter :: sinc$ = 36, comma$ = 37

character(20), parameter :: expression_op_name(37) = [character(20) :: '+', '-', '*', '/', &
                                    '(', ')', '^', '-', '+', '', 'sin', 'cos', 'tan', &
                                    'asin', 'acos', 'atan', 'abs', 'sqrt', 'log', 'exp', 'ran', &
                                    'ran_gauss', 'atan2', 'factorial', 'int', 'nint', 'floor', 'ceiling', &
                                    '?!+Numeric', '?!+Variable', 'mass_of', 'charge_of', 'anomalous_moment_of', &
                                    'species', '?!+Species', 'sinc', ',']


integer, parameter :: expression_eval_level(37) = [1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
                            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 0]

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
!
! Stack end elements not used are marked stack(i)%type = end_stack$
!
! Stack elements with stack(i)%type = variable$ are elements that need
! to be evaluated before calling evaluate_expression_stack.
!
! Input:
!   string    -- character(*): Expression to be converted.
!
! Output:
!   stack(:)  -- expression_atom_struct, allocatable: Expression evaluation stack.
!   n_stack   -- integer: number of "atoms" used by the expression
!   err_flag  -- logical: Set True if there is an error (EG divide by 0).
!   err_str   -- character(*): String describing the error.
!-

subroutine expression_string_to_stack (string, stack, n_stack, err_flag, err_str)

type (expression_atom_struct), allocatable :: stack(:)

integer i_op, i, var_type
integer op(100), ix_word, i_delim, i2, ix_word2, n_stack, op0, n_comma(100)

real(rp) value

character(*) string, err_str
character(1) delim, old_delim
character(80) word, word2
character(len(string)) parse_line


logical delim_found, split, zero_arg_function_pending
logical err_flag, err

! The general idea is to rewrite the expression on a stack in reverse polish.
! Reverse polish means that the operand goes last so that 2 * 3 is written 
! on the stack as: [2, 3, *]

! Since operations move towards the end of the stack we need a separate
! stack called op which keeps track of what operations have not yet
! been put on stack.

! init

delim = ''
err_flag = .true.
n_stack = 0
i_op = 0
zero_arg_function_pending = .false.
if (.not. allocated(stack)) allocate(stack(10))
stack(:)%type = end_stack$
parse_line = string

! parsing loop to build up the stack.

parsing_loop: do

  ! Get a word. If last thing was mass_of, etc., then next word is a species name.
  ! In this case, word is everything up to next parens.

  old_delim = delim
  op0 = 0
  var_type = variable$
  if (i_op > 1) op0 = op(i_op-1) ! If parsed "mass_of(" then op(i_op) corresponds to "("
  select case (op0)
  case (mass_of$, charge_of$, anomalous_moment_of$, species$) 
    call get_next_chunk (parse_line, word, ix_word, '()', delim, delim_found)
    var_type = species_const$
  case default
    call get_next_chunk (parse_line, word, ix_word, '[]+-*/()^,:} ', delim, delim_found)
  end select

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

  ! just make sure we are not chopping a number in two, e.g. "3.5d-7" should not
  ! get split at the "-" even though "-" is a delimiter

  split = .true.         ! assume initially that we have a split number
  if (ix_word == 0) then
    split = .false.
  elseif (word(ix_word:ix_word) /= 'E' .and. word(ix_word:ix_word) /= 'D' .and. &
          word(ix_word:ix_word) /= 'e' .and. word(ix_word:ix_word) /= 'd') then
    split = .false.
  endif
  if (delim(1:1) /= '-' .and. delim(1:1) /= '+') split = .false.
  do i = 1, ix_word-1
    if (index('.0123456789', word(i:i)) /= 0) cycle
    split = .false.
    exit
  enddo

  ! If still SPLIT = .TRUE. then we need to unsplit

  if (split) then
    word = word(:ix_word) // delim
    do i = 1, len(parse_line)
      if (index('0123456789', parse_line(i:i)) /= 0) cycle
      word = word(:ix_word+1) // parse_line(1:i-1)
      parse_line = parse_line(i:)
      ix_word = ix_word + i
      exit
    enddo

    call get_next_chunk (parse_line, word2, ix_word2, '+-*/()^,:}', delim, delim_found)
    if (ix_word2 /= 0) then
      err_str = 'Malformed number: ' // trim(word) // word2
      return
    endif

  endif

  ! Something like "lcav[lr(2).freq]" will get split on the "["

  if (delim == '[') then
    call get_next_chunk (parse_line, word2, ix_word2, ']', delim, delim_found)
    if (delim /= ']') then
      err_str = 'No "]" found to match "["'
      return
    endif
    word = word(:ix_word) // '[' // trim(word2) // ']'
    ix_word = ix_word + ix_word2 + 2
    call get_next_chunk (parse_line, word2, ix_word2, '+-*/()^,:} ', delim, delim_found)
    word = word(:ix_word) // word2
    ix_word = ix_word + ix_word2
  endif

  !---------------------------
  ! Now see what we got...

  ! For a "(" delim we must have a function

  if (delim == '(') then

    zero_arg_function_pending = .false.
    if (ix_word /= 0) then
      select case (upcase(word))
      case ('SIN') 
        call pushit (op, i_op, sin$)
      case ('SINC') 
        call pushit (op, i_op, sinc$)
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
        n_comma(i_op) = 0
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
      case ('CHARGE_OF')
        call pushit (op, i_op, charge_of$)
      case ('MASS_OF')
        call pushit (op, i_op, mass_of$)
      case ('SPECIES')
        call pushit (op, i_op, species$)
      case ('ANOMALOUS_MOMENT_OF')
        call pushit (op, i_op, anomalous_moment_of$)
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
      call push_numeric_or_var(word, err, op, i_op); if (err) return
    endif

    do
      do i = i_op, 1, -1     ! release pending ops
        if (op(i) == l_parens$) exit          ! break do loop
        call pushit_stack (stack, n_stack, op(i))
      enddo

      if (i == 0) then
        err_str = 'UNMATCHED ")" ON RHS'
        return
      endif

      i_op = i - 1

      call get_next_chunk (parse_line, word, ix_word, '+-*/()^,:}', delim, delim_found)
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
      if (old_delim == "*" .and. delim == "*") then
        err_str = 'EXPONENTIATION "**" NEEDS TO BE REPLACED BY "^"'
      else
        err_str = 'CONSTANT OR VARIABLE MISSING'
      endif
      return
    endif
    call push_numeric_or_var (word, err, op, i_op); if (err) return
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
  case (',')
    i_delim = comma$
  case ('}', ':', ')')   ! End of expression delims
    i_delim = no_delim$
  case default
    err_str = 'MALFORMED EXPRESSION'
    return
  end select

  ! now see if there are operations on the OP stack that need to be transferred
  ! to the STACK stack

  do i = i_op, 1, -1
    if (expression_eval_level(op(i)) < expression_eval_level(i_delim)) exit

    if (op(i) == l_parens$) then
      if (i > 1 .and. op(max(1,i-1)) == atan2$ .and. i_delim == comma$) then
        if (n_comma(i-1) /= 0) then
          err_str = 'TOO MANY COMMAS IN ATAN2 CONSTRUCT'
          return
        endif
        n_comma(i-1) = n_comma(i-1) + 1
        i_op = i
        cycle parsing_loop
      endif
      err_str = 'UNMATCHED "("'
      return
    endif

    if (op(i) == atan2$ .and. n_comma(i) /= 1) then
      err_str = 'MALFORMED ATAN2 ARGUMENT'
      return
    endif
    call pushit_stack (stack, n_stack, op(i))
  enddo

  i_op = i

  ! put the pending operation on the OP stack

  if (i_delim == no_delim$ .or. i_delim == comma$) then
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

subroutine get_next_chunk (parse_line, word, ix_word, delim_list, delim, delim_found)

character(*) parse_line, word, delim_list, delim

integer ix_word

logical delim_found

!

old_delim = delim
call word_read (parse_line, delim_list, word, ix_word, delim, delim_found, parse_line)

end subroutine get_next_chunk

!-------------------------------------------------------------------------
! contains

subroutine push_numeric_or_var (word, err, op, i_op)

logical err
integer op(:), i_op, ios, ix
character(*) word

!

err = .true.

if (is_real(word)) then
  call pushit_stack (stack, n_stack, numeric$)
  read (word, *, iostat = ios) stack(n_stack)%value
  if (ios /= 0) then
    err_str = 'BAD NUMERIC: ' // trim(word)
    return
  endif
  stack(n_stack)%name = word

  if (i_op > 0) then
    if (op(i_op) == unary_minus$) then
      stack(n_stack)%name = '-' // stack(n_stack)%name
      stack(n_stack)%value = -stack(n_stack)%value
      i_op = i_op - 1
    endif
  endif

else
  call pushit_stack (stack, n_stack, var_type)
  stack(n_stack)%name = word
  ! "my_species" in "mass_of(my_species)" is considered a variable and not a species_const
  if (var_type == species_const$) then
    stack(n_stack)%value = species_id(word)
    if (species_id(word) == invalid$) stack(n_stack)%type = variable$
  endif
endif

err = .false.

end subroutine push_numeric_or_var

!-------------------------------------------------------------------------
! contains

subroutine pushit_stack (stack, ix_stack, value)

type (expression_atom_struct), allocatable :: stack(:), temp_stack(:)
integer ix_stack, value

!

if (ix_stack == size(stack)) then
  call move_alloc(stack, temp_stack)
  allocate (stack(ix_stack+10))
  stack(1:ix_stack) = temp_stack
endif

call pushit (stack%type, ix_stack, value)

end subroutine pushit_stack

end subroutine expression_string_to_stack

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine pushit (array, ix_arr, value)
!
! Private routine used by expression_string_to_stack
!-

subroutine pushit (array, ix_arr, value)

implicit none

integer array(:), ix_arr, value
character(*), parameter :: r_name = 'pushit'

!

ix_arr = ix_arr + 1

if (ix_arr > size(array)) then
  call out_io (s_fatal$, r_name, 'ARRAY OVERFLOW, EXPERT HELP IS NEEDED!')
  if (global_com%exit_on_error) call err_exit
  return
endif

array(ix_arr) = value

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
! Note: Stack elements with stack(i)%type == variable$ need to be evelauated before
! calling this routine and the value placed in stack(i)%value.
!
! Input:
!   stack(:)    -- expression_atom_struct: Expression to evaluate.
!   var(:)      -- controller_var1_struct, optional: Array of variables.
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
type (expression_atom_struct) stack2(size(stack))
type (controller_var1_struct), optional :: var(:)

real(rp) value
integer i, i2, ix

logical, optional :: use_old
logical err_flag

character(*) err_str
character(*), parameter :: r_name = 'evaluate_expression_stack'

!

err_flag = .true.

i2 = 0
do i = 1, size(stack)

  select case (stack(i)%type)

  case (end_stack$)
    exit

  case (numeric$, variable$)
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

  case (sinc$)
    stack2(i2)%value = sinc(stack2(i2)%value)

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

  case (mass_of$)
    stack2(i2)%value = mass_of(nint(stack2(i2)%value))

  case (charge_of$)
    stack2(i2)%value = charge_of(nint(stack2(i2)%value))

  case (anomalous_moment_of$)
    stack2(i2)%value = anomalous_moment_of(nint(stack2(i2)%value))

  case (ran$)
    i2 = i2 + 1
    call ran_uniform(stack2(i2)%value)

  case (ran_gauss$)
    i2 = i2 + 1
    call ran_gauss(stack2(i2)%value)

  case (species_const$)
    i2 = i2 + 1
    stack2(i2)%value = stack(i)%value

  case (species$)
    ! Nothing to do

  case default

    if (is_attribute(stack(i)%type, control_var$)) then
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
      if (global_com%exit_on_error) then
        call out_io (s_fatal$, r_name, err_str)
        call err_exit
      endif
      return
    endif

  end select
enddo

if (i2 /= 1) then
  err_str = 'INTERNAL ERROR #03: GET HELP'
  if (global_com%exit_on_error) then
    call out_io (s_fatal$, r_name, err_str)
    call err_exit
  endif
  return
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
!   str       -- character(:), allocatable: Expression in string form.
!-

function expression_stack_to_string (stack, polish) result (str)

type (expression_atom_struct), target :: stack(:)
type (expression_atom_struct), pointer :: atom, atom2
type (expression_atom_struct) s2(size(stack))

character(:), allocatable :: str
character(200) s2_name(size(stack))

integer i, i2, ix
logical, optional :: polish

!

allocate (character(1) :: str)
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

    if (is_attribute(atom%type, all_control_var$)) then
      i2 = i2 + 1
      s2(i2)%type = variable$

      if (atom%name == '') then
        s2_name(i2) = real_to_string(atom%value, 20, 14)
      else
        s2_name(i2) = trim(atom%name)
      endif
      cycle
    endif

    select case (atom%type)
    case (plus$, minus$, times$, divide$, power$)
      if (expression_eval_level(s2(i2-1)%type) < expression_eval_level(atom%type)) s2_name(i2-1) = '(' // trim(s2_name(i2-1)) // ')'
      if (atom%type == minus$ .or. atom%type == divide$) then
        if (expression_eval_level(s2(i2)%type) <= expression_eval_level(atom%type)) s2_name(i2) = '(' // trim(s2_name(i2)) // ')'
      else
        if (expression_eval_level(s2(i2)%type) < expression_eval_level(atom%type)) s2_name(i2) = '(' // trim(s2_name(i2)) // ')'
      endif
      s2_name(i2-1) = trim(s2_name(i2-1)) // trim(expression_op_name(atom%type)) // s2_name(i2)
      s2(i2-1)%type = atom%type
      i2 = i2 - 1

    case (numeric$, variable$, species_const$)
      i2 = i2 + 1
      s2(i2)%type = atom%type
      if (atom%name == '') then
        s2_name(i2) = real_to_string(atom%value, 20, 14)
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
! Subroutine split_expression_string (expr, width, indent, lines)
!
! Routine to break an expression into a number of lines for a nicer display.
! Used when printing expressions.
!
! Input:
!   expr      -- character(*): String containing the expression.
!   width     -- integer: Maximum width of split expression.
!   indent    -- integer: Indent for every line after the first.
!
! Output:
!   lines(:)  -- character(*), allocatable: Split expression.
!-

subroutine split_expression_string (expr, width, indent, lines)

integer width, indent
integer nn, n0, nl, ind, ww, i_split, large, i, j

real(rp), parameter :: weight(11) = [1.0_rp, 1.0_rp, 1.0_rp, 0.95_rp, 0.95_rp, 0.8_rp, 0.8_rp, 0.8_rp, 0.8_rp, 0.8_rp, 0.8_rp]
real(rp) score

character(*) expr
character(*), allocatable :: lines(:)
character(len(expr)) ex
character(width), allocatable :: li(:)
character(1), parameter :: ch_split(11) = ['+', '-', ',', '*', '/', ')', ']', '}', '(', '[', '{']

!

ex = expr
ind = 0  ! zero indent for first line.
nl = 0
nn = len_trim(expr)
n0 = 1 + 2*nn/(width-indent)
call re_allocate (li, n0, .true., '')

do
  nl = nl + 1
  nn = len_trim(ex)
  ww = width - ind
  if (nn <= ww) then
    li(nl)(ind+1:) = ex
    exit
  endif

  i_split = ww
  score = 0
  do i = 1, size(ch_split)
    j = index(ex(1:ww), ch_split(i), back = .true.)
    if (j == 0) cycle
    if (j * weight(i) < score) cycle
    i_split = j
    score = j * weight(i)
  enddo

  li(nl)(ind+1:) = ex(1:i_split)  
  ex = adjustl(ex(i_split+1:))

  ind = indent
enddo

call re_allocate (lines, nl)
do i = 1, nl
  lines(i) = li(i)
enddo

end subroutine split_expression_string

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

if (n == 1 .and. is_attribute(stack(1)%type, control_var$)) then
  err_flag = .false.
  coef = 1
  return
endif

! Expression must have times at the end to be linear.

if (stack(n)%type /= times$) return

! To be linear, stack(1) or stack(n-1) must be variable and have no other variables.

if (is_attribute(stack(1)%type, control_var$)) then
  i0 = 2
  i1 = n-1
elseif (is_attribute(stack(n-1)%type, control_var$)) then
  i0 = 1
  i1 = n-2
else
  return
endif

call evaluate_expression_stack(stack(i0:i1), coef, err_flag, err_str)

end function linear_coef

end module
