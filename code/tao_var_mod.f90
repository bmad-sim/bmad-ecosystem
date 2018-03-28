module tao_var_mod

use tao_interface
use input_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function tao_user_is_terminating_optimization () result (is_terminating)
!
! Routine to check for keyboard input of a period '.' signaling optimization termination.
!
! Module needed:
!   use tao_var_mod
!
! Output:
!   is_terminating -- logical: Set True of '.' is detected. False otherwise.
!-

function tao_user_is_terminating_optimization () result (is_terminating)

implicit none

logical is_terminating

character(52) :: r_name = 'tao_user_is_terminating_optimization'
character(1) char

!

is_terminating = .false.
if (.not. s%global%optimizer_allow_user_abort) return

do
  call get_tty_char (char, .false., .false.) 
  if (char == '.') then
    call out_io (s_blank$, r_name, 'Optimizer stop signal detected.', 'Stopping now.')
    is_terminating = .true.
    return
  endif
  if (char == achar(0)) return   ! return if there is no more input
enddo

end function tao_user_is_terminating_optimization

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_get_opt_vars (var_value, var_step, var_delta, var_weight, var_ix,
!                                             ignore_if_weight_is_zero, ignore_if_not_limited)
!
! Subroutine to get the values of the variables used in optimization and put them
! in an array.
!
! Input:
! 
! Output:
!   var_value(:)       -- Real, allocatable, optional: Variable model values.
!   var_step(:)        -- Real, allocatable, optional: Variable step sizes.
!   var_delta(:)       -- Real, allocatable, optional: Variable Merit deltas.
!   var_weight(:)      -- Real, allocatable, optional: Variable weights in the merit function.
!   var_ix(:)          -- Integer, allocatable, optional: Variable s%var(:) indexes
!   ignore_if_weight_is_zero -- Logical, optional: If present and True then ignore
!                                 all variables whose merit weight is zero.
!   ignore_if_not_limited    -- Logical, optional: If present and True then ignore
!                                 all variables with limit constraint that are not limited.
!-

subroutine tao_get_opt_vars (var_value, var_step, var_delta, var_weight, var_ix, &
                                    ignore_if_weight_is_zero, ignore_if_not_limited)

implicit none

real(rp), allocatable, optional :: var_value(:), var_delta(:)
real(rp), allocatable, optional :: var_step(:), var_weight(:)
integer, allocatable, optional :: var_ix(:)

integer i, j
integer n_var

logical, optional :: ignore_if_weight_is_zero, ignore_if_not_limited
logical ignore_weight_is_zero, ignore_not_limited

! Count number of variables

ignore_weight_is_zero = logic_option(.false., ignore_if_weight_is_zero)
ignore_not_limited    = logic_option(.false., ignore_if_not_limited)

n_var = 0
do i = 1, s%n_var_used
  if (.not. s%var(i)%useit_opt) cycle
  if (ignore_weight_is_zero .and. s%var(i)%weight == 0) cycle
  if (ignore_not_limited .and. s%var(i)%merit_type == 'limit' .and. s%var(i)%delta_merit == 0) cycle
  n_var = n_var + 1
enddo

! Allocate arrays

if (present(var_value))   call re_allocate (var_value, n_var)
if (present(var_delta))   call re_allocate (var_delta, n_var)
if (present(var_step))    call re_allocate (var_step, n_var)
if (present(var_weight))  call re_allocate (var_weight, n_var)
if (present(var_ix))      call re_allocate (var_ix, n_var)

! Load info into arrays

j = 0
do i = 1, s%n_var_used
  if (.not. s%var(i)%useit_opt) cycle
  if (ignore_weight_is_zero .and. s%var(i)%weight == 0) cycle
  if (ignore_not_limited .and. s%var(i)%merit_type == 'limit' .and. s%var(i)%delta_merit == 0) cycle
  j = j + 1
  if (present(var_value))        var_value(j)   = s%var(i)%model_value
  if (present(var_delta))        var_delta(j)   = s%var(i)%delta_merit
  if (present(var_step))         var_step(j)    = s%var(i)%step
  if (present(var_weight))       var_weight(j)  = s%var(i)%weight
  if (present(var_ix))           var_ix(j)      = i
enddo

end subroutine tao_get_opt_vars

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_set_opt_vars (var_vec, print_limit_warning)
!
! Subrutine to set variable values from a vector of values. 
! This routine is used with optimization since optimimizers
! generally like their variables in the form of a vector.
!
! This routine assumes that the variables in each s%u(i) universe
! gets the same values. 
!
! Input:
!   var_vec(:) -- Real(rp): Vector of variables. 
!   print_limit_warning
!         -- Logical, optional: Print a warning if the value is past the variable's limits.
!             Default is True.
!
! Output:
!-

subroutine tao_set_opt_vars (var_vec, print_limit_warning)

implicit none

real(rp) var_vec(:)
integer i, j, k
logical, optional :: print_limit_warning

! Transfer the values from var_vec to the variables of each universe.

j = 0
do i = 1, s%n_var_used
  if (.not. s%var(i)%useit_opt) cycle
  j = j + 1
  call tao_set_var_model_value (s%var(i), var_vec(j), print_limit_warning)
enddo

end subroutine tao_set_opt_vars

end module
