module tao_evaluate_mod

use tao_utils

! used for parsing expressions
integer, parameter, private :: plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
integer, parameter, private :: l_parens$ = 5, r_parens$ = 6, power$ = 7
integer, parameter, private :: unary_minus$ = 8, unary_plus$ = 9, no_delim$ = 10
integer, parameter, private :: sin$ = 11, cos$ = 12, tan$ = 13
integer, parameter, private :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17, sqrt$ = 18
integer, parameter, private :: log$ = 19, exp$ = 20, ran$ = 21, ran_gauss$ = 22
integer, parameter, private :: numeric$ = 100

integer, parameter, private :: eval_level(22) = (/ 1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
                            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 /)

character(8), private :: wild_type_com

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_to_real (expression, value, err_flag, good)
!
! Mathematically evaluates an expression.
!
! Input:
!   expression   -- character(*): arithmetic expression
!  
! Output:
!   value        -- real(rp): Value of arithmetic expression.
!   err_flag     -- Logical: TRUE on error.
!   good         -- Logical, optional: Is the value valid? 
!                     Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                     orbit.x[23]|good_user is False.
!-

subroutine tao_to_real (expression, value, err_flag)

implicit none

character(*), intent(in) :: expression
real(rp) value
real(rp), allocatable, save :: vec(:)
logical, allocatable, save :: ok(:)
logical err_flag
logical, optional :: good

!

wild_type_com = 'BOTH'
call tao_evaluate_expression (expression, 1, vec, ok, &
                      .true., err_flag, tao_param_value_routine)
if (err_flag) return
value = vec(1)
if (present(good)) good = ok(1)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_to_real_vector (expression, wild_type, n_size, value, good, err_flag)
!
! Mathematically evaluates an expression.
!
! Input:
!   expression   -- Character(*): Arithmetic expression.
!   wild_type    -- Character(*): If something like "*|meas" is in the 
!                     expression does this refer to data or variables? 
!                     Possibilities are "DATA", "VAR", and "BOTH"
!   n_size       -- Integer: Size of the value array. If the expression
!                              is a scaler then the value will be spread.
!                              If n_size = 0 then the natural size determined 
!                              by expression is used.
!  
! Output:
!   value(:)     -- Real(rp), allocatable: Value of arithmetic expression.
!   good(:)      -- Logical, allocatable: Is the value valid? 
!                     Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                     orbit.x[23]|good_user is False.
!   err_flag     -- Logical: True on error. False otherwise
!-

subroutine tao_to_real_vector (expression, wild_type, n_size, value, good, err_flag)

use random_mod

implicit none

real(rp), allocatable :: value(:)
logical, allocatable :: good(:)

integer n_size

character(*), intent(in) :: expression, wild_type
character(16) :: r_name = "tao_to_real_vector"

logical err_flag, err, wild

!

wild_type_com = wild_type
call tao_evaluate_expression (expression, n_size, value, good, &
                                .true., err_flag, tao_param_value_routine)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_expression (expression, n_size, value, good, &
!                         zero_divide_print_err, err_flag, param_value_routine)
!
! Mathematically evaluates a character expression.
!
! Input:
!   expression            -- Character(*): Arithmetic expression.
!   n_size                -- Integer: Size of the value array. If the expression
!                              is a scaler then the value will be spread.
!                              If n_size = 0 then the natural size determined 
!                              by expression is used.
!   zero_divide_print_err -- Logical: If False just return zero without printing
!                             an error message.
!   param_value_routine   -- Subroutine: Routine to translate a variable to a value.
!   
! Output:
!   value(:)     -- Real(rp), allocatable: Value of arithmetic expression.
!   good(:)      -- Logical, allocatable: Is the value valid? 
!                     Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                     orbit.x[23]|good_user is False.
!   err_flag     -- Logical: TRUE on error.
!-

subroutine tao_evaluate_expression (expression, n_size, value, good, &
                          zero_divide_print_err, err_flag, param_value_routine)

use random_mod

implicit none

interface
  subroutine param_value_routine (str, value, good, err_flag)
    use tao_struct
    implicit none
    character(*) str
    real(rp), allocatable :: value(:)
    logical, allocatable :: good(:)
    logical err_flag
  end subroutine
end interface

type (tao_eval_stack_struct) stk(100)

integer i_lev, i_op, i, ios, n, n_size, n__size
integer op(200), ix_word, i_delim, i2, ix, ix_word2, ixb

real(rp), allocatable :: value(:)
logical(rp), allocatable :: good(:)

character(*), intent(in) :: expression
character(len(expression)) phrase
character(1) delim
character(40) word, word2
character(40) :: r_name = "tao_evaluate_expression"
character(40) saved_prefix

logical delim_found, split, ran_function_pending
logical err_flag, err, wild, zero_divide_print_err

! Don't destroy the input expression

err_flag = .true.
saved_prefix = ''

phrase = expression

! if phrase is blank then return 0.0

call string_trim (phrase, phrase, ios)
if (ios == 0) then
  call out_io (s_warn$, r_name, &
    "Expression is blank", len(phrase))
  value = 0.0
  return
endif
 
! General idea: Create a reverse polish stack that represents the expression.
! Reverse polish means that the operand goes last so that 2 * 3 is writen 
! on the stack as: [2, 3, *]

! The stack is called: stk
! Since operations move towards the end of the stack we need a separate
! stack called op which keeps track of what operations have not yet
! been put on stk.

! init

err_flag = .false.
i_lev = 0
i_op = 0
ran_function_pending = .false.

do i = 1, size(stk)
  if (allocated(stk(i)%good)) deallocate (stk(i)%good)
enddo

! parsing loop to build up the stack.

parsing_loop: do

  ! get a word

  call word_read (phrase, '+-*/()^,:}[ ', word, ix_word, delim, &
                    delim_found, phrase)

  !  if (delim == '*' .and. word(1:1) == '*') then
  !    call out_io (s_warn$, r_name, 'EXPONENTIATION SYMBOL IS "^" AS OPPOSED TO "**"')
  !    err_flag = .true.
  !    return
  !  endif

  if (ran_function_pending .and. (ix_word /= 0 .or. delim /= ')')) then
        call out_io (s_warn$, r_name, &
                   'RAN AND RAN_GAUSS DO NOT TAKE AN ARGUMENT')
    err_flag = .true.
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
  elseif (word(ix_word:ix_word) /= 'E' .and. word(ix_word:ix_word) /= 'e' ) then
    split = .false.
  endif
  if (delim /= '-' .and. delim /= '+') split = .false.
  do i = 1, ix_word-1
    if (index('.0123456789', word(i:i)) == 0) split = .false.
  enddo

  ! If still SPLIT = .TRUE. then we need to unsplit

  if (split) then
    word = word(:ix_word) // delim
    call word_read (phrase, '+-*/()^,:}', word2, ix_word2, delim, &
                    delim_found, phrase)
    word = word(:ix_word+1) // word2
    ix_word = ix_word + ix_word2
  endif

  ! Something like "lcav[lr(2).freq]" will get split on the "["

  if (delim == '[') then
    call word_read (phrase, ']', word2, ix_word2, delim, &
                    delim_found, phrase)
    if (.not. delim_found) then
      call out_io (s_warn$, r_name, "NO MATCHING ']' FOR OPENING '[':" // expression)
      err_flag = .true.
      return
    endif
    word = word(:ix_word) // '[' // trim(word2) // ']'
    ix_word = ix_word + ix_word2 + 2
    if (phrase(1:1) /= ' ') then  ! even more...
      call word_read (phrase, '+-*/()^,:}', word2, ix_word2, delim, &
                                                  delim_found, phrase)
      word = word(:ix_word) // trim(word2)       
      ix_word = ix_word + ix_word2 
    endif
  endif

  ! If delim = "*" then see if this is being used as a wildcard

  if (delim == '*') then
    ixb = index(phrase, '|')
    if (ixb /= 0) then
      wild = .true.
      if (index(phrase(1:ixb), '+') /= 0) wild = .false.
      if (index(phrase(1:ixb), '-') /= 0) wild = .false.
      if (index(phrase(1:ixb), '/') /= 0) wild = .false.
      if (index(phrase(1:ixb), '^') /= 0) wild = .false.
      if (index(phrase(1:ixb), '(') /= 0) wild = .false.
      ix = index(phrase(1:ixb), '*')
      if (ix /= 0) then
        if (ix == 1) then
          wild = .false.
        elseif (phrase(ix-1:ix-1) /= '.' .and. phrase(ix-1:ix-1) /= '@') then
          wild = .false.
        endif
      endif
      if (wild) then
        word = word(:ix_word) // '*' // phrase(1:ixb)
        phrase = phrase(ixb+1:)
        call word_read (phrase, '+-*/()^,:}', word2, ix_word2, delim, &
                                                  delim_found, phrase)
        word = trim(word) // trim(word2)       
        ix_word = len_trim(word)
      endif
    endif
  endif

  !---------------------------
  ! Now see what we got...

  ! For a "(" delim we must have a function

  if (delim == '(') then

    ran_function_pending = .false.
    if (ix_word /= 0) then
      word2 = word
      call downcase_string (word2)
      select case (word2)
      case ('sin')
        call pushit (op, i_op, sin$)
      case ('cos')
        call pushit (op, i_op, cos$)
      case ('tan') 
        call pushit (op, i_op, tan$)
      case ('asin') 
        call pushit (op, i_op, asin$)
      case ('acos') 
        call pushit (op, i_op, acos$)
      case ('atan') 
        call pushit (op, i_op, atan$)
      case ('abs') 
        call pushit (op, i_op, abs$)
      case ('sqrt') 
        call pushit (op, i_op, sqrt$)
      case ('log') 
        call pushit (op, i_op, log$)
      case ('exp') 
        call pushit (op, i_op, exp$)
      case ('ran') 
        call pushit (op, i_op, ran$)
        ran_function_pending = .true.
      case ('ran_gauss') 
        call pushit (op, i_op, ran_gauss$)
        ran_function_pending = .true.
      case default
        call out_io (s_warn$, r_name, &
               'UNEXPECTED CHARACTERS ON RHS BEFORE "(": ')
        err_flag = .true.
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

    call pushit (op, i_op, unary_plus$)
    cycle parsing_loop

  ! for a ")" delim

  elseif (delim == ')') then
    if (ix_word == 0) then
      if (.not. ran_function_pending) then
        call out_io (s_warn$, r_name, 'CONSTANT OR VARIABLE MISSING BEFORE ")"')
        err_flag = .true.
        return
      endif
    else
      call pushit (stk%type, i_lev, numeric$)
      call all_value_routine (word, stk(i_lev), err_flag)
      if (err_flag) return
    endif

    do
      do i = i_op, 1, -1       ! release pending ops
        if (op(i) == l_parens$) exit            ! break do loop
        call pushit (stk%type, i_lev, op(i))
      enddo

      if (i == 0) then
        call out_io (s_warn$, r_name, 'UNMATCHED ")" ON RHS')
        err_flag = .true.
        return
      endif

      i_op = i - 1

      call word_read (phrase, '+-*/()^,:}', word, ix_word, delim, &
                    delim_found, phrase)
      if (ix_word /= 0) then
        call out_io (s_warn$, r_name, &
                   'UNEXPECTED CHARACTERS ON RHS AFTER ")"')
        err_flag = .true.
        return
      endif

      if (delim /= ')') exit    ! if no more ')' then no need to release more
    enddo


    if (delim == '(') then
      call out_io (s_warn$, r_name, '")(" CONSTRUCT DOES NOT MAKE SENSE')
      err_flag = .true.
      return
    endif

  ! For binary "+-/*^" delims

  else
    if (ix_word == 0) then
      call out_io (s_warn$, r_name, 'CONSTANT OR VARIABLE MISSING')
      err_flag = .true.
      return
    endif
    call pushit (stk%type, i_lev, numeric$)
    call all_value_routine (word, stk(i_lev), err_flag)
    if (err_flag) return
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
      call out_io (s_error$, r_name, 'INTERNAL ERROR')
      call err_exit
  end select

  ! now see if there are operations on the OP stack that need to be transferred
  ! to the STK stack

  do i = i_op, 1, -1
    if (eval_level(op(i)) >= eval_level(i_delim)) then
      if (op(i) == l_parens$) then
        call out_io (s_warn$, r_name, 'UNMATCHED "("')
        err_flag = .true.
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
! Now go through the stack and perform the operations...
! First some error checks

if (i_op /= 0) then
  call out_io (s_warn$, r_name, 'UNMATCHED "("')
  err_flag = .true.
  return
endif

if (i_lev == 0) then
  call out_io (s_warn$, r_name, 'NO VALUE FOUND')
  err_flag = .true.
  return
endif

n__size = 1
do i = 1, i_lev
  if (stk(i)%type /= numeric$) cycle
  n = size(stk(i)%value)
  if (n == 1) cycle
  if (n__size == 1) n__size = n
  if (n /= n__size) then
    call out_io (s_warn$, r_name, 'ARRAY SIZE MISMATCH')
    err_flag = .true.
    return
  endif
enddo

if (n_size /= 0) then
  if (n__size /= 1 .and. n_size /= n__size) then
    call out_io (s_warn$, r_name, 'ARRAY SIZE MISMATCH')
    err_flag = .true.
    return
  endif
  n__size = n_size
endif

! Calculate good

call re_allocate (good, n__size)
good = .true.
do i = 1, i_lev
  if (.not. allocated(stk(i)%good)) cycle
  if (size(stk(i)%good) == 1) then; good = good .and. stk(i)%good(1:1)
  else;                             good = good .and. stk(i)%good
  endif
enddo

!

i2 = 0  ! stack pointer
do i = 1, i_lev

  select case (stk(i)%type)
  case (numeric$) 
    i2 = i2 + 1
    stk(i2)%value = stk(i)%value

  case (unary_minus$) 
    stk(i2)%value = -stk(i2)%value

  case (unary_plus$) 
    stk(i2)%value = stk(i2)%value

  case (plus$) 
    if (size(stk(i2)%value) < size(stk(i2-1)%value)) then
      stk(i2-1)%value = stk(i2-1)%value + stk(i2)%value(1)
    elseif (size(stk(i2)%value) > size(stk(i2-1)%value)) then
      call value_transfer (stk(i2-1)%value, stk(i2-1)%value(1) + stk(i2)%value)
    else
      stk(i2-1)%value = stk(i2-1)%value + stk(i2)%value
    endif
    i2 = i2 - 1

  case (minus$) 
    if (size(stk(i2)%value) < size(stk(i2-1)%value)) then
      stk(i2-1)%value = stk(i2-1)%value - stk(i2)%value(1)
    elseif (size(stk(i2)%value) > size(stk(i2-1)%value)) then
      call value_transfer (stk(i2-1)%value, stk(i2-1)%value(1) - stk(i2)%value)
    else
      stk(i2-1)%value = stk(i2-1)%value - stk(i2)%value
    endif
    i2 = i2 - 1

  case (times$) 
    if (size(stk(i2)%value) < size(stk(i2-1)%value)) then
      stk(i2-1)%value = stk(i2-1)%value * stk(i2)%value(1)
    elseif (size(stk(i2)%value) > size(stk(i2-1)%value)) then
      call value_transfer (stk(i2-1)%value, stk(i2-1)%value(1) * stk(i2)%value)
    else
      stk(i2-1)%value = stk(i2-1)%value * stk(i2)%value
    endif
    i2 = i2 - 1

  case (divide$) 
    if (any(stk(i2)%value == 0)) then
      stk(1)%value = 0
      if (zero_divide_print_err) call out_io (s_warn$, r_name, 'DIVIDE BY 0 ON RHS')
      err_flag = .true.
      return
    endif
    if (size(stk(i2)%value) < size(stk(i2-1)%value)) then
      stk(i2-1)%value = stk(i2-1)%value / stk(i2)%value(1)
    elseif (size(stk(i2)%value) > size(stk(i2-1)%value)) then
      call value_transfer (stk(i2-1)%value, stk(i2-1)%value(1) / stk(i2)%value)
    else
      stk(i2-1)%value = stk(i2-1)%value / stk(i2)%value
    endif
    i2 = i2 - 1

  case (power$) 
    if (size(stk(i2)%value) < size(stk(i2-1)%value)) then
      stk(i2-1)%value = stk(i2-1)%value ** stk(i2)%value(1)
    elseif (size(stk(i2)%value) > size(stk(i2-1)%value)) then
      call value_transfer (stk(i2-1)%value, stk(i2-1)%value(1) ** stk(i2)%value)
    else
      stk(i2-1)%value = stk(i2-1)%value ** stk(i2)%value
    endif
    i2 = i2 - 1

  case (sin$) 
    stk(i2)%value = sin(stk(i2)%value)

  case (cos$) 
    stk(i2)%value = cos(stk(i2)%value)

  case (tan$) 
    stk(i2)%value = tan(stk(i2)%value)

  case (asin$) 
    stk(i2)%value = asin(stk(i2)%value)

  case (acos$) 
    stk(i2)%value = acos(stk(i2)%value)

  case (atan$) 
    stk(i2)%value = atan(stk(i2)%value)

  case (abs$) 
    stk(i2)%value = abs(stk(i2)%value)

  case (sqrt$) 
    stk(i2)%value = sqrt(stk(i2)%value)

  case (log$) 
    stk(i2)%value = log(stk(i2)%value)

  case (exp$) 
    stk(i2)%value = exp(stk(i2)%value)

  case (ran$) 
    i2 = i2 + 1
    call re_allocate(stk(i2)%value, n__size)
    call ran_uniform(stk(i2)%value)

  case (ran_gauss$) 
    i2 = i2 + 1
    call re_allocate(stk(i2)%value, n__size)
    call ran_gauss(stk(i2)%value)

  case default
    call out_io (s_warn$, r_name, 'INTERNAL ERROR')
    err_flag = .true.
    return
  end select
enddo

if (i2 /= 1) call out_io (s_warn$, r_name, 'INTERNAL ERROR')

if (size(stk(1)%value) == 1 .and. n__size > 1) then
  call re_allocate (value, n_size)
  value = stk(1)%value(1)
else
  call value_transfer (value, stk(1)%value)
endif

!-------------------------------------------------------------------------
contains

subroutine value_transfer (to_array, from_array)

real(rp), allocatable :: to_array(:)
real(rp) from_array(:)

!

call re_allocate (to_array, size(from_array))
to_array = from_array

end subroutine

!-------------------------------------------------------------------------
! contains

subroutine pushit (stack, i_lev, value)

implicit none

integer stack(:), i_lev, value

character(6) :: r_name = "pushit"

!

i_lev = i_lev + 1

if (i_lev > size(stack)) then
  call out_io (s_warn$, r_name, 'STACK OVERFLOW.')
  call err_exit
endif

stack(i_lev) = value

end subroutine pushit
                       
!---------------------------------------------------------------------------
! contains

subroutine all_value_routine (str, stack, err_flag)

type (tao_eval_stack_struct) stack

integer ios, i, n

character(*) str
character(40) str2

logical err_flag

!

if (is_real(str)) then
  allocate (stack%value(1))
  read (str, *, iostat = ios) stack%value(1)
  if (ios /= 0) then
    call out_io (s_warn$, r_name, "This doesn't seem to be a number: " // str)
    err_flag = .true.
    return
  endif

else
  ! Remember last string so that 'orbit.x|meas-ref' translates to orbit.x|meas - orbit.x|ref.
  ix = index(str, '|')
  if (ix == 0) then
    str2 = trim(saved_prefix) // str
  else
    saved_prefix = str(1:ix-1)
    str2 = str
  endif
  call param_value_routine (str2, stack%value, stack%good, err_flag)
  
endif

end subroutine

end subroutine tao_evaluate_expression

!---------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_param_value_routine (str, value, good, err_flag)

implicit none

type (tao_real_array_struct), allocatable, save :: re_array(:)
type (tao_integer_array_struct), allocatable, save :: int_array(:)

real(rp), allocatable :: value(:)
logical, allocatable :: good(:)
integer ios, i, n

character(*) str
character(40) :: r_name = 'tao_param_value_routine'

logical err_flag

!

if (wild_type_com == 'VAR' .or. wild_type_com == 'BOTH') &
               call tao_find_var (err_flag, str, re_array = re_array, print_err = .false.)

if (.not. allocated(re_array) .and. &
                              (wild_type_com == 'DATA' .or. wild_type_com == 'BOTH')) then
  call tao_find_data (err_flag, str, re_array = re_array, &
                                                int_array = int_array, print_err = .false.)
endif

if (size(re_array) /= 0) then
  n = size(re_array)
  call re_allocate (value, n)
  call re_allocate (good, n)
  do i = 1, n
    value(i) = re_array(i)%r
    good(i)  = re_array(i)%good
  enddo

elseif (size(int_array) /= 0) then
  n = size(int_array)
  call re_allocate (value, n)
  call re_allocate (good, n)
  do i = 1, n
    value(i) = int_array(i)%i
    good(i)  = .true.
  enddo

else
  call out_io (s_warn$, r_name, "This doesn't seem to be datum or variable value: " // str)
  err_flag = .true.
  return
endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_int (str, i_int, err)
! 
! Converts a string to an integer
!
! If the string str is blank then i_int = 0
!-

subroutine tao_to_int (str, i_int, err)

character(*) str
integer ios, i_int
logical err
character(12) :: r_name = "tao_to_int"

!

  call string_trim (str, str, ios)
  if (ios .eq. 0) then
    i_int = 0
    return
  endif
 
  err = .false.
  read (str, *, iostat = ios) i_int

  if (ios /= 0) then
    call out_io (s_error$, r_name, 'EXPECTING INTEGER: ' // str)
    err = .true.
    return
  endif

end subroutine

end module
