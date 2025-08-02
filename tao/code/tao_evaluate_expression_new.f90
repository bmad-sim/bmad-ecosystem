!+
! Subroutine tao_evaluate_expression_new (expression, n_size, use_good_user, value, err_flag, print_err, &
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
subroutine tao_evaluate_expression_new (expression, n_size, use_good_user, value, err_flag, print_err, &
                      info, stack, dflt_component, dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, &
                      dflt_dat_or_var_index, dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)

use tao_data_and_eval_mod, dummy => tao_evaluate_expression_new
use expression_mod

implicit none

type (expression_tree_struct), target :: tree
type (expression_tree_struct), pointer :: base
type (tao_eval_stack1_struct), allocatable, optional :: stack(:)
type (tao_eval_stack1_struct), allocatable :: stk(:)
type (ele_struct), optional, pointer :: dflt_ele_ref, dflt_ele_start, dflt_ele
type (coord_struct), optional :: dflt_orbit
type (tao_expression_info_struct), allocatable, optional :: info(:)
type (tao_data_struct), optional :: datum

real(rp), allocatable :: value(:)
real(rp), optional :: dflt_s_offset

integer n_size
integer, optional :: dflt_uni, dflt_eval_point
integer i, n, n_stk, ix

logical use_good_user, err_flag
logical, optional :: print_err
logical printit

character(*) :: expression
character(*), optional :: dflt_component, dflt_source
character(*), optional :: dflt_dat_or_var_index

character(60) saved_prefix
character(80) default_source
character(2000) :: phrase, err_str

character(*), parameter :: r_name = "tao_evaluate_expression_new"

! 

err_flag = .true.
printit = logic_option(.true., print_err)
default_source = ''
if (present(dflt_source)) default_source = dflt_source

call tao_expression_hash_substitute(expression, phrase, dflt_ele)
call expression_asterisk_substitute(phrase)

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

!

call expression_string_to_tree(phrase, tree, err_flag, err_str)

if (err_flag) then
  call out_io(s_error$, r_name, 'Error parsing expression: ' // err_str)
  return
endif

if (size(tree%node) /= 1 .or. tree%node(1)%type /= comma$) then
  call out_io(s_error$, r_name, 'Expression not well formed: ' // expression)
  return
endif

base => tree%node(1)
call expression_tree_asterisk_restore(base)

if (s%global%verbose_on) then
  call type_expression_tree(tree)
endif

err_flag = .false.
n_stk = 0
call expression_tree_to_stack(base, stk, n_stk, expression, err_flag, .false.)
if (err_flag) return

if (s%global%verbose_on) then
  print *, '! New =========================================='
  do i = 1, n_stk
    if (allocated(stk(i)%value)) then
      print '(a20, es12.4, t50, i0)', stk(i)%name, stk(i)%value(1), stk(i)%type
    else
      print '(a20, t50, i0)', stk(i)%name, stk(i)%type
    endif
  enddo
  print *, '!=========================================='
  print *
endif

! Evaluate individual values

do i = 1, n_stk
  select case (stk(i)%type)
  case (plus$, unary_plus$, minus$, unary_minus$, times$, divide$, power$, function$)
  case default
    ! saved_prefix is used so that something like 'orbit.x|meas-ref' can be evaluated as 'orbit.x|meas - orbit.x|ref.'
    saved_prefix = ''
    if (i > 1) then
      ix = index(stk(i-1)%name, '|')
      if (ix > 0) saved_prefix = stk(i-1)%name(1:ix-1)
    endif

    ! Don't try to find the value of a species name.
    ! A species name will appear just before (since the stack is in reverse Polish) a species related function.
    if (i < n_stk) then
      select case (stk(i+1)%name)
      case ('species', 'mass_of', 'charge_of', 'anomalous_moment_of')
        stk(i)%type = species_const$
        cycle
      end select
    endif

    call tao_param_value_routine (stk(i)%name, use_good_user, saved_prefix, stk(i), err_flag, printit, &
               dflt_component, default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, &
               dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)
    if (err_flag) return
  end select
enddo

! Evaluate expression

call tao_evaluate_stack_new (stk(1:n_stk), n_size, use_good_user, value, err_flag, printit, expression, info)

! If the stack argument is present then copy stk to stack

if (present(stack)) then
  if (allocated(stack)) deallocate(stack)
  allocate (stack(n_stk))
  do i = 1, n_stk
    if (allocated (stk(i)%value)) then
      n = size(stk(i)%value)
      allocate (stack(i)%value(n), stack(i)%info(n))
      if (allocated (stack(i)%value_ptr)) allocate (stack(i)%value_ptr(n))
    endif
    stack(i) = stk(i)
  enddo
endif


!------------------------------------------------------------------------------------
contains

! Substitute "??" for "*" characters that are being used as wildcards and not as multiplication symbols. 
! This is done so as to not confuse expression tree creation.
! Wildcard examples: 
!    "[*]|", "*.*|", "*.x|", "*@orbit.x|", "*@*|", "orbit.*[3]|", "ele::q*1[beta_a]", "var::*d|model"
!    "quad_k1[*]|model"
! Non-wildcard examples:
! "emit*data::mf.xm|model", "3*[1,2]", "3*.42", "a*[b,c]"

! Wild if "*" is bounded by:
!   "::" < "*" < "|"
!   "::" < "*" < "["

subroutine expression_asterisk_substitute(phrase)

character(*) phrase
character(1) left_char, right_char
integer i0, istar, ii
logical wild

!

i0 = 1
main_loop: do 
  istar = index(phrase(i0:), '*')
  if (istar == 0) return

  istar = istar  + i0 - 1
  i0 = istar + 1

  left_char = '!'   ! Beginning of line
  do ii = istar-1, 1, -1
    left_char = phrase(ii:istar)
    select case (left_char)
    case ('+', '-', '*', '/', '^', ' ', ']', '[', '@', '(', ')'); exit
    end select
  enddo

  do ii = istar+1, len_trim(phrase)
    right_char = phrase(ii:istar)
    select case (right_char)
    case ('+', '-', '/', '^', ' ', ']', '[', '@', '(', ')'); exit
    end select
  enddo

  wild = .false.
  if (left_char == ':' .and. (right_char == '[' .or. right_char == '|')) wild = .true.
  if (phrase(istar-1:istar-1) == '[' .and. phrase(istar+1:istar+1) == ']') wild = .true.
  if (phrase(istar+1:istar+1) == '@') wild = .true.
  if (.not. wild) cycle

  ! Is a wild card
  phrase = phrase(1:istar-1) // '??' // phrase(istar+1:)
enddo main_loop

end subroutine expression_asterisk_substitute

!----------------------------------------------------------------------------------------------------
! contains

recursive subroutine expression_tree_asterisk_restore(tree)

type (expression_tree_struct), target :: tree
integer in, ix

!

do
  ix = index(tree%name, '??')
  if (ix == 0) exit
  tree%name = tree%name(1:ix-1) // '*' // tree%name(ix+1:)
enddo

if (associated(tree%node)) then
  do in = 1, size(tree%node)
    call expression_tree_asterisk_restore(tree%node(in))
  enddo
endif

end subroutine

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
! contains

recursive subroutine expression_tree_to_stack(tree, stk, n_stk, expression, err_flag, in_compound)

implicit none

type (expression_tree_struct), target :: tree
type (tao_eval_stack1_struct), allocatable :: stk(:)

integer in, n_stk
logical err_flag, in_compound, in_comp

character(*) expression
character(*), parameter :: r_name = 'expression_tree_to_stack'

! Note: Tao variable, data, and paramter name syntax does not use parens.

in_comp = in_compound

select case (tree%type)
case (function$)
  call push_stack (stk, n_stk, tree%type, tree%name)

case (comma$, parens$, func_parens$)
  ! No op

case (compound$)
  !in_comp = .true.

case (curly_brackets$)
  if (in_compound) stk(n_stk)%name = '{' // stk(n_stk)%name

case (square_brackets$)
  if (in_compound) stk(n_stk)%name = '[' // stk(n_stk)%name


case default
  if (in_compound) then
    stk(n_stk)%name = trim(tree%name) // stk(n_stk)%name
  else
    call push_stack (stk, n_stk, tree%type, tree%name)
  endif
end select

!

if (.not. associated(tree%node)) return
do in = 1, size(tree%node)
  call expression_tree_to_stack(tree%node(in), stk, n_stk, expression, err_flag, in_comp)
  if (err_flag) return
enddo

!

select case (tree%type)
case (curly_brackets$)
  if (in_compound) stk(n_stk)%name = '}' // stk(n_stk)%name
case (square_brackets$)
  if (in_compound) stk(n_stk)%name = ']' // stk(n_stk)%name
end select

end subroutine expression_tree_to_stack

!-------------------------------------------------------------------------
! contains

subroutine push_stack (stack, n_stk, this_type, this_name)

type (tao_eval_stack1_struct), allocatable :: stack(:), tmp_stk(:)
integer n_stk, this_type

character(*) this_name
character(*), parameter :: r_name = "push_stack"

!

n_stk = n_stk + 1

if (n_stk > size(stack)) then
  call move_alloc(stack, tmp_stk)
  allocate(stack(2*n_stk))
  stack(1:n_stk-1) = tmp_stk
endif

! Tao uses numeric$ instead of variable$ or constant#

select case (this_type)
case (variable$, constant$)
  stack(n_stk)%type = numeric$
case default
  stack(n_stk)%type = this_type
end select

stack(n_stk)%name = this_name
stack(n_stk)%scale = 1

end subroutine push_stack
                       
end subroutine tao_evaluate_expression_new

