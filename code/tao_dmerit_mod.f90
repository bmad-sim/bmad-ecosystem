module tao_dmerit_mod

use tao_interface

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_dModel_dVar_calc (force_calc, err_flag)
!
! Subroutine to calculate the dModel_dVar derivative matrix.
!
! Input:
!   s          -- Super_universe_struct:
!   force_calc -- Logical: If true then force recalculation of the matrix.
!                  If False then only calculate matrix if it doesn't exist.
!   veto_vars_with_zero_dmodel 
!              -- Logical, optional (default False): Veto variables where
!                  all dModel_dvar for that var are zero.
!                  Sets the var%good_var logical to False.
!
! Output:
!   s        -- Super_universe_struct.
!    %u(:)%dModel_dVar(:,:)  -- Derivative matrix
!   err_flag -- logical: Set true if there is an error. False otherwise.
!-

subroutine tao_dmodel_dvar_calc (force_calc, err_flag)

implicit none

type (tao_universe_struct), pointer :: u
type (tao_data_struct), pointer :: dat

real(rp) model_value, merit_value

integer i, j, k
integer n_data, n_var, nd, nv

character(20) :: r_name = 'tao_dmodel_dvar_calc'

logical :: err_flag
logical reinit, force_calc, calc_ok, err_message_out

! make sure size of matrices are correct.
! We only compute the derivitive matrix if needed or if force_calc is set.

call tao_set_var_useit_opt
call tao_set_data_useit_opt
err_flag = .false.

reinit = force_calc 
n_var = count (s%var%useit_opt)

do i = lbound(s%u, 1), ubound(s%u, 1)

  u => s%u(i)
  if (.not. u%is_on) then
    if (allocated(u%dModel_dVar)) deallocate(u%dModel_dVar)
    cycle
  endif

  n_data = count (u%data%useit_opt)

  if (.not. allocated(u%dModel_dVar)) then
    allocate (u%dModel_dVar(n_data, n_var))
    reinit = .true.
  endif

  if (size(u%dModel_dVar, 1) /= n_data .or. size(u%dModel_dVar, 2) /= n_var) then
    deallocate (u%dModel_dVar)
    allocate (u%dModel_dVar(n_data, n_var))
    u%dModel_dVar = 0
    reinit = .true.
  endif

  nd = 0
  do j = 1, size(u%data)
    if (.not. u%data(j)%useit_opt) then
      u%data(j)%ix_dmodel = -1
      cycle
    endif
    nd = nd + 1
    if (u%data(j)%ix_dModel /= nd) reinit = .true.
    u%data(j)%ix_dModel = nd
  enddo

enddo

nv = 0
do j = 1, s%n_var_used
  if (.not. s%var(j)%useit_opt) then
    s%var(j)%ix_dvar = -1
    cycle
  endif
  nv = nv + 1
  if (s%var(j)%ix_dVar /= nv) reinit = .true.
  s%var(j)%ix_dVar = nv
enddo

if (.not. reinit) then
  return
endif

!--------------------------------------------------------------
! Calculate derivative matrices.

call out_io (s_info$, r_name, 'Remaking dModel_dVar derivative matrix.', &
                              'This may take a while...', '')

! Note that the var model_value needs to be saved before
! switching to the design lattice since var%model_value is a pointer.

if (s%global%derivative_uses_design) then
  do j = 1, s%n_var_used
    s%var(j)%scratch_value = s%var(j)%model_value  ! Save
    call tao_set_var_model_value (s%var(j), s%var(j)%design_value)
  enddo
  do i = lbound(s%u, 1), ubound(s%u, 1)
    s%u(i)%scratch_lat = s%u(i)%model%lat ! Save
    s%u(i)%model%lat = s%u(i)%design%lat
  enddo
endif

! In an actual machine one would not take small differences like this so turn off
! measurement noise.

s%com%add_measurement_noise = .false.
merit_value = tao_merit (calc_ok)
if (.not. calc_ok) then
  call out_io (s_error$, r_name, 'BASE MODEL CALCULATION HAS PROBLEMS', ' ')
  s%com%add_measurement_noise = .true.  ! Reset
  err_flag = .true.
  return
endif

s%var%old_value = s%var%delta_merit

! Save old data

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)
  if (.not. u%is_on) cycle
  do j = 1, size(u%data)
    if (.not. u%data(j)%useit_opt) cycle
    if (u%data(j)%good_model) then
      u%data(j)%old_value = u%data(j)%delta_merit
    else
      u%data(j)%old_value = real_garbage$
    endif
  enddo
enddo

! Loop over all variables to vary

do j = 1, s%n_var_used
  if (.not. s%var(j)%useit_opt) cycle

  nv = s%var(j)%ix_dvar
  if (s%var(j)%step == 0) then
    call out_io (s_error$, r_name, 'VARIABLE STEP SIZE IS ZERO FOR: ' // tao_var1_name(s%var(j)))
    err_flag = .true.
    s%com%add_measurement_noise = .true.  ! Reset
    return
  endif
  model_value = s%var(j)%model_value
  call tao_set_var_model_value (s%var(j), model_value + s%var(j)%step)
  merit_value = tao_merit (calc_ok)

  err_message_out = .false.   ! Want to print the error message only once per variable.
  do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)
    do k = 1, size(u%data)
      if (.not. u%data(k)%useit_opt) cycle
      dat => u%data(k)
      nd = dat%ix_dmodel
      if (dat%good_model .and. dat%old_value /= real_garbage$) then
        u%dModel_dVar(nd,nv) = (dat%delta_merit - dat%old_value) / s%var(j)%step
        cycle
      endif

      ! Could not compute derivative...
      ! Error meassges are only generated when the datum is good with one setting
      ! of the variable and bad with the other setting.

      u%dModel_dVar(nd,nv) = 0

      if (err_message_out) cycle
      if (.not. dat%good_model .and. dat%old_value == real_garbage$) cycle

      call out_io (s_error$, r_name, 'ERROR IN CALCULATING DERIVATIVE MATRIX.', &
                      'VARIABLE STEP SIZE IS TOO LARGE(?) FOR: ' // tao_var1_name(s%var(j)))
      err_flag = .true.
      err_message_out = .true.
    enddo
  enddo

  call tao_set_var_model_value (s%var(j), model_value)

  ! look for keyboard input to end optimization

  if (tao_user_is_terminating_optimization()) then
    err_flag = .true.
    exit
  endif

enddo

! End

if (s%global%derivative_uses_design) then
  do i = lbound(s%u, 1), ubound(s%u, 1)
    s%u(i)%model%lat = s%u(i)%scratch_lat
  enddo
  do j = 1, s%n_var_used
    call tao_set_var_model_value (s%var(j), s%var(j)%scratch_value)
  enddo
endif

s%com%add_measurement_noise = .true.  ! Reset
merit_value = tao_merit () ! to reset all values

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_veto_vars_with_zero_dmodel ()
!
! Routine to veto all variables with zero effect on data used in the merit function.
!-

subroutine tao_veto_vars_with_zero_dmodel ()

type (tao_universe_struct), pointer :: u
integer i, j, nv
logical zero_dmodel
character(40) :: r_name = 'tao_veto_vars_with_zero_dmodel'

!

do j = 1, s%n_var_used
  if (.not. s%var(j)%useit_opt) cycle
  nv = s%var(j)%ix_dvar
  zero_dmodel = .true.
  do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)
    if (.not. u%is_on) cycle
    if (any(u%dmodel_dvar(:,nv) /= 0)) zero_dmodel = .false.
  enddo
  if (zero_dmodel) then
    s%var(j)%good_var = .false.
    call out_io (s_info$, r_name, 'Data is independent of Variable: ' // tao_var1_name(s%var(j)), &
                                  'Will not use this variable in optimization.')
  endif
enddo

call tao_set_var_useit_opt()

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_dmerit_calc ()
!-

subroutine tao_dmerit_calc ()

type (tao_data_struct), pointer :: data

integer i, j, k, nv, nd
logical err_flag

!

call tao_dmodel_dvar_calc (.false., err_flag)

s%var(:)%dmerit_dvar = 0

do i = 1, s%n_var_used

  if (.not. s%var(i)%useit_opt) cycle
  s%var(i)%dmerit_dvar = 2 * s%var(i)%weight * s%var(i)%delta_merit
  nv = s%var(i)%ix_dvar

  do j = lbound(s%u, 1), ubound(s%u, 1)
    do k = 1, size (s%u(j)%data)
      data => s%u(j)%data(k)
      if (.not. data%useit_opt) cycle
      nd = data%ix_dmodel
      s%var(i)%dmerit_dvar = s%var(i)%dmerit_dvar + 2 * data%weight * &
                                   s%u(j)%dmodel_dvar(nd,nv) * data%delta_merit
    enddo
  enddo

end do

end subroutine

end module
