module tao_var_mod

use tao_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_get_opt_vars (var_value, var_step, var_delta, var_weight, var_ix, ignore_if_weight_is_zero)
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
!-

subroutine tao_get_opt_vars (var_value, var_step, var_delta, var_weight, var_ix, ignore_if_weight_is_zero)

implicit none

real(rp), allocatable, optional :: var_value(:), var_delta(:)
real(rp), allocatable, optional :: var_step(:), var_weight(:)
integer, allocatable, optional :: var_ix(:)

integer i, j
integer n_var

logical, optional :: ignore_if_weight_is_zero

! 

  if (logic_option(.false., ignore_if_weight_is_zero)) then
    n_var  = count(s%var(:)%useit_opt .and. s%var(:)%weight /= 0)
  else
    n_var  = count(s%var(:)%useit_opt)
  endif

  if (present(var_value))   call re_allocate (var_value, n_var)
  if (present(var_delta))   call re_allocate (var_delta, n_var)
  if (present(var_step))    call re_allocate (var_step, n_var)
  if (present(var_weight))  call re_allocate (var_weight, n_var)
  if (present(var_ix))      call re_allocate (var_ix, n_var)

  j = 0
  do i = 1, size(s%var)
    if (.not. s%var(i)%useit_opt) cycle
    if (logic_option(.false., ignore_if_weight_is_zero) .and. s%var(i)%weight == 0) cycle
    j = j + 1
    if (present(var_value))        var_value(j)   = s%var(i)%model_value
    if (present(var_delta))        var_delta(j)   = s%var(i)%delta_merit
    if (present(var_step))         var_step(j)    = s%var(i)%step
    if (present(var_weight))       var_weight(j)  = s%var(i)%weight
    if (present(var_ix))           var_ix(j)      = i
  enddo


end subroutine

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
  do i = 1, size(s%var)
    if (.not. s%var(i)%useit_opt) cycle
    j = j + 1
    call tao_set_var_model_value (s%var(i), var_vec(j), print_limit_warning)
  enddo

end subroutine


end module
