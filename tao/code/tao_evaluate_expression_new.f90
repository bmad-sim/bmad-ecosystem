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
!   value(:)   -- Real(rp), allocatable: Value of arithmetic expression.
!   err_flag   -- Logical: True on an error. EG: Invalid expression.
!                   A divide by zero is not an error but good(:) will be set to False.
!   info(:)    -- tao_expression_info_struct, allocatable, optional: Is the value valid?, etc.
!                   Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                   orbit.x[23]|good_user is False.
!   stack(:)   -- tao_eval_node_struct, allocatable, optional: Array of nodes of variable names.
!                   This is useful to check what datums or variables are used in the expression.
!-

recursive &
subroutine tao_evaluate_expression_new (expression, n_size, use_good_user, value, err_flag, print_err, &
                      info, stack, dflt_component, dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, &
                      dflt_dat_or_var_index, dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)

use tao_data_and_eval_mod, dummy => tao_evaluate_expression_new
use tao_expression_tree_mod, dummy2 => tao_evaluate_expression_new
use expression_mod

implicit none

type (expression_tree_struct), target :: tree
type (tao_eval_node_struct), allocatable, optional :: stack(:)
type (tao_eval_node_struct) tao_tree
type (ele_struct), optional, pointer :: dflt_ele_ref, dflt_ele_start, dflt_ele
type (coord_struct), optional :: dflt_orbit
type (tao_expression_info_struct), allocatable, optional :: info(:)
type (tao_data_struct), optional :: datum

real(rp), allocatable :: value(:)
real(rp), optional :: dflt_s_offset

integer n_size
integer, optional :: dflt_uni, dflt_eval_point
integer i, n, ix

logical use_good_user, err_flag
logical, optional :: print_err
logical printit

character(*) :: expression
character(*), optional :: dflt_component, dflt_source
character(*), optional :: dflt_dat_or_var_index

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

if (s%global%verbose_on) call type_expression_tree(tree)

err_flag = .false.
tao_tree = tao_eval_node_struct(tree%type, tree%name, 1.0_rp, null(), null(), null(), null())
call bmad_tree_to_tao_tree(tree, tao_tree, expression, err_flag); if (err_flag) return
call deallocate_tree(tree)
call tree_param_evaluate(tao_tree, expression, err_flag); if (err_flag) return

if (s%global%verbose_on) call tao_type_expression_tree(tao_tree)

call tao_evaluate_tree (tao_tree, n_size, use_good_user, value, err_flag, printit, expression, info)

!

if (present(stack)) then
  n = 0
  call tree_to_stack(tao_tree, stack, n)
endif

call tao_deallocate_tree(tao_tree)

!------------------------------------------------------------------------------------
contains

! Substitute "??" for "*" characters that are being used as wildcards and not as multiplication symbols. 
! This is done so as to not confuse expression tree creation.
! Wildcard examples: 
!    "[*]|", "*.*|", "*.x|", "*@orbit.x|", "*@*|", "orbit.*[3]|", "ele::q*1[beta_a]", "var::*d|model"
!    "quad_k1[*]|model"
! Non-wildcard examples:
! "emit*data::mf.xm|model", "3*[1,2]", "3*.42", "a*[b,c]"
!
! Wild if "*" is bounded by:
!   "::" < "*" < "|"  or
!   "::" < "*" < "["  or
!   "*@"
 
subroutine expression_asterisk_substitute(phrase)

character(*) phrase
character(1) left_char, right_char, ch
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
    ch = phrase(ii:ii)
    select case (ch)
    case (':', '|', '+', '-', '*', '/', '^', ' ', ']', '[', '@', '(', ')')
      left_char = ch
      exit
    end select
  enddo

  right_char = '!'
  do ii = istar+1, len_trim(phrase)
    ch = phrase(ii:ii)
    select case (ch)
    case (':', '|', '+', '-', '/', '^', ' ', ']', '[', '@', '(', ')')
      right_char = ch
      exit
    end select
  enddo

  wild = .false.
  if (left_char == ':' .and. (right_char == '[' .or. right_char == '|')) wild = .true.
  if (istar > 1) then
    if (phrase(istar-1:istar-1) == '[' .and. phrase(istar+1:istar+1) == ']') wild = .true.
  endif
  if (phrase(istar+1:istar+1) == '@') wild = .true.
  if (.not. wild) cycle

  ! Is a wild card
  phrase = phrase(1:istar-1) // '??' // phrase(istar+1:)
enddo main_loop

end subroutine expression_asterisk_substitute

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
! contains

recursive subroutine bmad_tree_to_tao_tree(tree, tao_tree, expression, err_flag)

implicit none

type (expression_tree_struct), target :: tree
type (expression_tree_struct), pointer :: node
type (tao_eval_node_struct), target :: tao_tree
type (tao_eval_node_struct), pointer :: tnode

integer in, ix, n_node
logical err_flag, split_variable

character(*) expression
character(*), parameter :: r_name = 'bmad_tree_to_tao_tree'

! Note: Tao variable, data, and paramter name syntax does not use parens.

if (.not. associated(tree%node)) return
n_node = size(tree%node)

! A variable or datume like "3@orbit.x|model" is split by the parser into an array of nodes.
! Needed is to combine the nodes
! Note: "abc[...]" is a var or datum. "abc + [...]" is not.

select case (tree%type)
case(square_brackets$, func_parens$, parens$, curly_brackets$)
  split_variable = .false.
case default
  split_variable = (n_node > 1)
  do in = 1, n_node
    node => tree%node(in)
    select case (node%type)
    case (plus$, minus$, unary_plus$, unary_minus$, times$, divide$, power$, func_parens$)
      split_variable = .false.
    end select
  enddo
end select

if (split_variable) then
  allocate(tao_tree%node(1))
  tnode => tao_tree%node(1)
  tnode%type = variable$
  tnode%name = expression_tree_to_string (tree, .false.)
  return
endif

!

allocate(tao_tree%node(n_node))

do in = 1, n_node
  tnode => tao_tree%node(in)
  tnode = tao_eval_node_struct(tree%node(in)%type, tree%node(in)%name, 1.0_rp, null(), null(), null(), null())
  select case (tnode%type)
  case (square_brackets$, func_parens$, parens$, comma$, compound$)
    call bmad_tree_to_tao_tree(tree%node(in), tao_tree%node(in), expression, err_flag)
    if (err_flag) return
  end select
enddo

end subroutine bmad_tree_to_tao_tree

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
! contains

recursive subroutine tree_param_evaluate(tao_tree, expression, err_flag)

implicit none

type (tao_eval_node_struct), target :: tao_tree
type (tao_eval_node_struct), pointer :: tnode, snode

integer in, ix, n_node
logical err_flag, in_compound, in_comp

character(*) expression
character(60) saved_prefix
character(*), parameter :: r_name = 'tree_param_evaluate'

!

if (.not. associated(tao_tree%node)) return
n_node = size(tao_tree%node)
saved_prefix = ''

! Evaluate

do in = 1, n_node
  tnode => tao_tree%node(in)

  do
    ix = index(tnode%name, '??')
    if (ix == 0) exit
    tnode%name = tnode%name(1:ix-1) // '*' // tnode%name(ix+2:)
  enddo

  select case (tnode%type)
  ! A species name will appear just before (since the stack is in reverse Polish) a species related function.
  case (function$)
    select case (tnode%name)
    case ('species', 'mass_of', 'charge_of', 'anomalous_moment_of')
      if (in == 1) then
        call out_io(s_error$, r_name, 'Misplaced function: ' // quote(tnode%name) // ' in expression: ' // expression)
        err_flag = .true.
        return
      elseif (tao_tree%node(in-1)%type /= func_parens$) then
        call out_io(s_error$, r_name, 'Missing parentheses for function: ' // quote(tnode%name) // ' in expression: ' // expression)
        err_flag = .true.
        return
      endif

      snode => tao_tree%node(in-1)
      snode%name = snode%node(1)%node(1)%name  ! Note: snode has a comma child
      snode%type = species_const$
      call tao_deallocate_tree(snode)
      call re_allocate(snode%value, 1)
      call tao_re_allocate_expression_info (snode%info, 1)
      snode%info%good = .true.

      if (allocated(s%com%symbolic_num)) then
        call match_word(snode%name, s%com%symbolic_num%name, ix, .true., .false.)
        if (ix > 0) then
          snode%value(1) = s%com%symbolic_num(ix)%value
          snode%type = constant$
        endif
      endif

      if (snode%type == species_const$) then
        snode%value(1) = species_id(snode%name)
        if (snode%value(1) == invalid$) then
          call out_io (s_error$, r_name, 'Not a valid species name or symbol: ' // snode%name, &
                                         ' in expression: ' // expression)
          err_flag = .true.
          return
        endif
      endif
    end select

  case (variable$, numeric$, constant$)
    ! saved_prefix is used so that something like 'orbit.x|meas-ref' can be evaluated as 'orbit.x|meas - orbit.x|ref.'

    if (saved_prefix /= '' .and. (tnode%name == 'design' .or. tnode%name == 'model' .or. tnode%name == 'base')) then
      tnode%name = trim(saved_prefix) // tnode%name
    else
      saved_prefix = ''
    endif

    call tao_param_value_routine (tnode%name, use_good_user, saved_prefix, tnode, err_flag, printit, &
               dflt_component, default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, &
               dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)
    if (err_flag) return

  case (square_brackets$, func_parens$, parens$, comma$, compound$)
    call tree_param_evaluate(tnode, expression, err_flag)
    if (err_flag) return
    if (tnode%type == compound$) then
      ix = index(tnode%node(1)%name, '|')
      saved_prefix = tnode%node(1)%name(1:ix)
    endif
  end select
enddo

end subroutine tree_param_evaluate

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
! contains

! Note: This just flattens the tree. The resulting stack is not usable for expression evaluation.

recursive subroutine tree_to_stack(tao_tree, stack, n_stk)

type (tao_eval_node_struct) :: tao_tree
type (tao_eval_node_struct), allocatable :: stack(:), tmp_stk(:)
integer in
integer n_stk

!

if (n_stk == 0) then
  if (allocated(stack)) deallocate(stack)
  allocate(stack(50))
endif

if (.not. associated(tao_tree%node)) return

do in = 1, size(tao_tree%node)
  n_stk = n_stk + 1
  if (n_stk > size(stack)) then
    call move_alloc(stack, tmp_stk)
    allocate(stack(2*n_stk))
    stack(1:n_stk-1) = tmp_stk(1:n_stk-1)
    deallocate(tmp_stk)
  endif
  stack(n_stk) = tao_tree%node(in)
  call tree_to_stack(tao_tree%node(in), stack, n_stk)
enddo

call move_alloc(stack, tmp_stk)
allocate(stack(n_stk))
stack(1:n_stk) = tmp_stk(1:n_stk)

end subroutine tree_to_stack

end subroutine tao_evaluate_expression_new

