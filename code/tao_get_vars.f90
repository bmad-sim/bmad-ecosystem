!+
! Subroutine tao_get_vars (s, var_value, var_del, var_weight, var_data_value)
!
! Subroutine to get the values of the variables used in optimization and put them
! in an array. It is important that the variables in different universes be the same.
!
! Input:
!   s    -- Tao_super_universe_struct: 
! 
! Output:
!   var_value(:)       -- Real, allocatable: Variable values
!   var_del(:)         -- Real, allocatable, optional: Variable step sizes.
!   var_weight(:)      -- Real, allocatable, optional: Variable weights in the merit function.
!   var_data_value(:)  -- Real, allocatable, optional: Variable values when the data was taken.
!-

subroutine tao_get_vars (s, var_value, var_del, var_weight, var_data_value)

use tao_mod

implicit none

type (tao_super_universe_struct), target :: s
type (tao_var_struct), pointer :: var(:)

real(rp), allocatable :: var_value(:)
real(rp), allocatable, optional :: var_del(:), var_weight(:), var_data_value(:)

integer i, j, k
integer n_var, n_data

! parallel_vars means that all variables s%u(:)%var(j) of different universes
! are in lock-step with each other.

if (s%global%parallel_vars) then

  var => s%u(1)%var

  n_var  = count(var(:)%useit_opt)
  call reallocate_real (var_value, n_var)
  if (present(var_del))        call reallocate_real (var_del, n_var)
  if (present(var_data_value)) call reallocate_real (var_data_value, n_var)

  j = 0
  do i = 1, size(var)
    if (.not. var(i)%useit_opt) cycle
    j = j + 1
    var_value(j) = var(i)%model_value
    if (present(var_del))    var_del(j) = var(i)%step
    if (present(var_weight)) var_weight(j) = var(i)%weight
    if (present(var_data_value))   var_data_value(j) = var(i)%data_value
  enddo

! Serial case where the variables of different universes are not correlated.

else

  n_var = 0
  do i = 1, size(s%u(i)%var)
    n_var = n_var + count(s%u(i)%var(:)%useit_opt)
  enddo

  call reallocate_real (var_value, n_var)
  call reallocate_real (var_del, n_var)
  call reallocate_real (var_data_value, n_var)

  j = 0
  do k = 1, size(s%u)
    var => s%u(i)%var
    do i = 1, size(var)
      if (.not. var(i)%useit_opt) cycle
      j = j + 1
      var_value(j) = var(i)%model_value
      if (present(var_del))    var_del(j) = var(i)%step
      if (present(var_weight)) var_weight(j) = var(i)%weight
      if (present(var_data_value))   var_data_value(j) = var(i)%data_value
    enddo
  enddo

endif

end subroutine

