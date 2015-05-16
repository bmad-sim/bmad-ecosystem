module expression_mod

use basic_bmad_mod

implicit none

integer, parameter :: plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
integer, parameter :: l_parens$ = 5, r_parens$ = 6, power$ = 7
integer, parameter :: unary_minus$ = 8, unary_plus$ = 9, no_delim$ = 10
integer, parameter :: sin$ = 11, cos$ = 12, tan$ = 13
integer, parameter :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17, sqrt$ = 18
integer, parameter :: log$ = 19, exp$ = 20, ran$ = 21, ran_gauss$ = 22
integer, parameter :: atan2$ = 23, factorial$ = 24, int$ = 25, nint$ = 26
integer, parameter :: floor$ = 27, ceiling$ = 28, numeric$ = 100

integer, parameter :: eval_level$(28) = [1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
                            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9]

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine evaluate_expression_stack (stack, value, err_flag, err_str, var)
!
! Routine to evaluate an mathematical expression represented by an "expression stack".
! Expression stacks are created by bmad_parser code.
!
! Input:
!   stack(:)    -- expression_stack_struct: Expression to evaluate.
!   var(:)      -- controller_var_struct, optional: Array of variables.
!
! Output:
!   value       -- real(rp): Value of the expression.
!   err_flag    -- logical: True if there is an evaluation problem. False otherwise.
!   err_str     -- character(*): Error string explaining error if there is one.
!-

subroutine evaluate_expression_stack (stack, value, err_flag, err_str, var)

use random_mod

type (expression_stack_struct) stack(:)
type (expression_stack_struct) stk2(20)
type (controller_var_struct), optional :: var(:)

real(rp) value
integer i, i2

logical err_flag

character(*) err_str

!

err_flag = .true.

i2 = 0
do i = 1, size(stack)

  select case (stack(i)%type)

  case (numeric$)
    i2 = i2 + 1
    stk2(i2)%value = stack(i)%value

  case (unary_minus$)
    stk2(i2)%value = -stk2(i2)%value

  case (unary_plus$)
    stk2(i2)%value = stk2(i2)%value

  case (plus$)
    stk2(i2-1)%value = stk2(i2-1)%value + stk2(i2)%value
    i2 = i2 - 1

  case (minus$)
    stk2(i2-1)%value = stk2(i2-1)%value - stk2(i2)%value
    i2 = i2 - 1

  case (times$)
    stk2(i2-1)%value = stk2(i2-1)%value * stk2(i2)%value
    i2 = i2 - 1

  case (divide$)
    if (stk2(i2)%value == 0) then
      err_str = 'DIVIDE BY 0 ON RHS'
      return
    endif
    stk2(i2-1)%value= stk2(i2-1)%value / stk2(i2)%value
    i2 = i2 - 1

  case (power$)
    stk2(i2-1)%value = stk2(i2-1)%value**stk2(i2)%value
    i2 = i2 - 1

  case (sin$)
    stk2(i2)%value = sin(stk2(i2)%value)

  case (cos$)
    stk2(i2)%value = cos(stk2(i2)%value)

  case (tan$)
    stk2(i2)%value = tan(stk2(i2)%value)

  case (asin$)
    if (stk2(i2)%value < -1 .or. stk2(i2)%value > 1) then
      err_str = 'ASIN ARGUMENT HAS MAGNITUDE GREATER THAN 1'
      return
    endif
    stk2(i2)%value = asin(stk2(i2)%value)

  case (acos$)
    if (stk2(i2)%value < -1 .or. stk2(i2)%value > 1) then
      err_str = 'ACOS ARGUMENT HAS MAGNITUDE GREATER THAN 1'
      return
    endif
    stk2(i2)%value = acos(stk2(i2)%value)

  case (factorial$)
    stk2(i2)%value = factorial(nint(stk2(i2)%value))
    if (stk2(i2)%value < 0) then
      err_str = 'FACTORIAL PROBLEM'
      return
    endif

  case (atan$)
    stk2(i2)%value = atan(stk2(i2)%value)

  case (atan2$)
    stk2(i2-1)%value = atan2(stk2(i2-1)%value, stk2(i2)%value)
    i2 = i2 - 1

  case (abs$)
    stk2(i2)%value = abs(stk2(i2)%value)

  case (sqrt$)
    if (stk2(i2)%value < 0) then
      err_str = 'SQRT ARGUMENT IS NEGATIVE '
      return
    endif
    stk2(i2)%value = sqrt(stk2(i2)%value)

  case (log$)
    if (stk2(i2)%value < 0) then
      err_str = 'LOG ARGUMENT IS NEGATIVE '
      return
    endif
    stk2(i2)%value = log(stk2(i2)%value)

  case (exp$)
    stk2(i2)%value = exp(stk2(i2)%value)

  case (int$)
    stk2(i2)%value = int(stk2(i2)%value)

  case (nint$)
    stk2(i2)%value = nint(stk2(i2)%value)

  case (floor$)
    stk2(i2)%value = floor(stk2(i2)%value)

  case (ceiling$)
    stk2(i2)%value = ceiling(stk2(i2)%value)

  case (ran$)
    i2 = i2 + 1
    call ran_uniform(stk2(i2)%value)

  case (ran_gauss$)
    i2 = i2 + 1
    call ran_gauss(stk2(i2)%value)

  case default
    err_str = 'INTERNAL ERROR #02: GET HELP'
    if (global_com%exit_on_error) call err_exit
  end select
enddo

if (i2 /= 1) then
  err_str = 'INTERNAL ERROR #03: GET HELP'
  if (global_com%exit_on_error) call err_exit
endif

value = stk2(1)%value
err_flag = .false.

end subroutine evaluate_expression_stack

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine expression_stack_to_string (stack, str)
!
! Routine to convert an expression stack to a string in polish notation.
!
! Input:
!   stack(:)  -- expression_stack_struct: arithmetic expression
!
! Output:
!   str       -- character(*): expression in polish noation
!-

subroutine expression_stack_to_string (stack, str)

type (expression_stack_struct) stack(:)
character(*) str

integer i

!

str = ''
do i = 1, size(stack)
  str = trim(str) // ', ' // stack(i)%name
enddo
str = str(3:)


end subroutine expression_stack_to_string

end module
