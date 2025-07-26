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
use random_mod
use expression_mod

implicit none

type (tao_eval_stack1_struct), allocatable, optional :: stack(:)
type (ele_struct), optional, pointer :: dflt_ele_ref, dflt_ele_start, dflt_ele
type (coord_struct), optional :: dflt_orbit
type (tao_expression_info_struct), allocatable, optional :: info(:)
type (tao_data_struct), optional :: datum

real(rp), allocatable :: value(:)
real(rp), optional :: dflt_s_offset

integer n_size
integer, optional :: dflt_uni, dflt_eval_point
integer ix

logical use_good_user, err_flag
logical, optional :: print_err
logical printit

character(*) :: expression
character(*), optional :: dflt_component, dflt_source
character(*), optional :: dflt_dat_or_var_index
character(80) default_source
character(2000) :: phrase

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

end subroutine tao_evaluate_expression_new
