module tao_var_mod

use tao_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_get_vars (var_value, var_step, var_delta, var_weight, &
!                          var_meas_value, var_ix_dVar)
!
! Subroutine to get the values of the variables used in optimization and put them
! in an array. It is important that the variables in different universes be the same.
!
! Input:
! 
! Output:
!   var_value(:)       -- Real, allocatable, optional: Variable model values.
!   var_step(:)        -- Real, allocatable, optional: Variable step sizes.
!   var_delta(:)       -- Real, allocatable, optional: Variable Merit deltas.
!   var_weight(:)      -- Real, allocatable, optional: Variable weights in the merit function.
!   var_meas_value(:)  -- Real, allocatable, optional: Variable values when the data was taken.
!   var_ix_dVar(:)     -- Real, allocatable, optional: Variable ix_dVar indices
!-

subroutine tao_get_vars (var_value, var_step, var_delta, var_weight, &
                         var_meas_value, var_ix_dVar)

implicit none

real(rp), allocatable, optional :: var_value(:), var_meas_value(:), var_delta(:)
real(rp), allocatable, optional :: var_step(:), var_weight(:), var_ix_dVar(:)

integer i, j
integer n_var

! 

  n_var  = count(s%var(:)%useit_opt)
  if (present(var_value))      call reallocate_real (var_value, n_var)
  if (present(var_delta))      call reallocate_real (var_delta, n_var)
  if (present(var_step))       call reallocate_real (var_step, n_var)
  if (present(var_meas_value)) call reallocate_real (var_meas_value, n_var)
  if (present(var_weight))     call reallocate_real (var_weight, n_var)
  if (present(var_ix_dVar))    call reallocate_real (var_ix_dVar, n_var)

  j = 0
  do i = 1, size(s%var)
    if (.not. s%var(i)%useit_opt) cycle
    j = j + 1
    if (present(var_value))        var_value(j)      = s%var(i)%model_value
    if (present(var_delta))        var_delta(j)      = s%var(i)%delta_merit
    if (present(var_step))         var_step(j)       = s%var(i)%step
    if (present(var_weight))       var_weight(j)     = s%var(i)%weight
    if (present(var_meas_value))   var_meas_value(j) = s%var(i)%meas_value
    if (present(var_ix_dVar))      var_ix_dVar(j)    = s%var(i)%ix_dVar
  enddo


end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_set_vars (var_vec)
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
!
! Output:
!-

subroutine tao_set_vars (var_vec)

implicit none

real(rp) var_vec(:)

integer i, j, k

! Transfer the values from var_vec to the variables of each universe.

  j = 0
  do i = 1, size(s%var)
    if (.not. s%var(i)%useit_opt) cycle
    j = j + 1
    call tao_set_var_model_value (s%var(i), var_vec(j))
  enddo

end subroutine


end module
