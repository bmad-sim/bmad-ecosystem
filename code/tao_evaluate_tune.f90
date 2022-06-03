!+
! Function tao_evaluate_tune (q_str, q0, delta_input) result (q_val)
!
! Routine to evaluate an expression for the tune.
!
! If delta_input is False and if the string evaluates to a value with zero integer part, the
! integer part of q_val will be set to the integer part of q0.
!
! q_val will be set to have the same sign as q0. This is important for the longitudinal tune above
! transition where the tune should be negative.
!
! Input:
!   q_str       -- character(*): String expression.
!   q0          -- real(rp): Default to use if q_str evaluates to zero.
!                   Also used to set the integer part of the tune.
!   delta_input -- logical: If true then qa_str and qb_str are deltas from present tune.
!
! Outut:
!   q_val       -- real(rp): Tune value. Set zero if there is an error.
!-

function tao_evaluate_tune (q_str, q0, delta_input) result (q_val)

use tao_data_and_eval_mod, dummy => tao_evaluate_tune
implicit none

real(rp) q0, q_val
real(rp), allocatable :: set_val(:)
logical delta_input, err
character(*) q_str

!

q_val = 0

if (q_str == '') then
  q_val = q0
  return
endif

!

call tao_evaluate_expression (q_str, 1, .true., set_val, err)
if (err) return

q_val = set_val(1) * sign_of(q0)

if (delta_input) then
  q_val = q0 + q_val
elseif (q_val == 0) then
  q_val = q0
elseif (int(q_val) == 0) then
  q_val = q_val + int(q0)
endif

end function tao_evaluate_tune
