module tao_dmerit_mod

use tao_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_dModel_dVar_calc (force_calc)
!
! Subroutine to calculate the dModel_dVar derivative matrix.
!
! Input:
!   s          -- Super_universe_struct:
!   force_calc -- Logical: If true then force recalculation of the matrix.
!                  If False then only calculate matrix if it doesn't exist.
!
! Output:
!   s       -- Super_universe_struct.
!    %u(:)%dModel_dVar(:,:)  -- Derivative matrix
!-

subroutine tao_dmodel_dvar_calc (force_calc)

implicit none

type (tao_universe_struct), pointer :: u

real(rp) model_value, merit_value
integer i, j, k, jj
integer n_data, n_var, nd, nv
integer s_var_size, current_s_var_point
character(20) :: r_name = 'tao_dmodel_dvar_calc'
logical reinit, force_calc, calc_ok

! make sure size of matrices are correct.
! We only compute the derivitive matrix if needed or if force_calc is set.

call tao_set_var_useit_opt
call tao_set_data_useit_opt

reinit = force_calc 

do i = 1, size(s%u)

  u => s%u(i)
  if (.not. u%is_on) then
    if (associated(u%dModel_dVar)) deallocate(u%dModel_dVar)
     cycle
  endif

  n_data = count (u%data%useit_opt)
  n_var = count (s%var%useit_opt)

  if (.not. associated(u%dModel_dVar)) then
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
    if (.not. u%data(j)%useit_opt) cycle
    nd = nd + 1
    if (u%data(j)%ix_dModel /= nd) reinit = .true.
    if (u%data(j)%good_model .and. &
            any (u%dModel_dVar(nd,:) == real_garbage$)) reinit = .true.
    u%data(j)%ix_dModel = nd
    if (u%data(j)%good_model) then
      u%data(j)%old_value = u%data(j)%delta_merit
    else
      u%data(j)%old_value = real_garbage$
    endif
  enddo

enddo

nv = 0
do j = 1, size(s%var)
  if (.not. s%var(j)%useit_opt) cycle
  nv = nv + 1
  if (s%var(j)%ix_dVar /= nv) reinit = .true.
  s%var(j)%ix_dVar = nv
enddo

if (.not. reinit) return
call out_io (s_info$, r_name, 'Remaking dModel_dVar derivative matrix.') 
call out_io (s_blank$, r_name, 'This may take a while...') 
call out_io (s_blank$, r_name, ' ') 

! Calculate matrices

merit_value = tao_merit (calc_ok)
if (.not. calc_ok) then
  call out_io (s_error$, r_name, 'BASE MODEL CALCULATION HAS PROBLEMS', ' ')
endif

s%var%old_value = s%var%delta_merit

s_var_size = size(s%var)
jj = 1

do j = 1, s_var_size

  if (.not. s%var(j)%useit_opt) cycle
  
! let user see progress

   if (modulo(jj * 10, n_var) .eq. 0) then
     call out_io (s_blank$, r_name, " \i3\% done...", (jj/n_var)*100)
     jj = jj + 1
   endif

  nv = s%var(j)%ix_dvar
  if (s%var(j)%step == 0) then
    call out_io (s_error$, r_name, 'VARIABLE STEP SIZE IS ZERO FOR: ' // tao_var1_name(s%var(j)))
    call err_exit
  endif
  model_value = s%var(j)%model_value
  call tao_set_var_model_value (s%var(j), model_value + s%var(j)%step)
  merit_value = tao_merit (calc_ok)

  if (.not. calc_ok) then
    call out_io (s_error$, r_name, &
         'VARIABLE STEP SIZE IS TOO LARGE FOR: ' // tao_var1_name(s%var(j)), ' ')
  endif

  do i = 1, size(s%u)
    u => s%u(i)
    do k = 1, size(u%data)
      if (.not. u%data(k)%useit_opt) cycle
      nd = u%data(k)%ix_dmodel
      if (u%data(k)%good_model .and. u%data(k)%old_value /= real_garbage$) then
        u%dModel_dVar(nd,nv) = (u%data(k)%delta_merit - u%data(k)%old_value) / s%var(j)%step
      else
        u%dModel_dVar(nd,nv) = real_garbage$
      endif
    enddo
  enddo

  call tao_set_var_model_value (s%var(j), model_value)

enddo

merit_value = tao_merit () ! to reset all values

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

!

call tao_dmodel_dvar_calc (.false.)

s%var(:)%dmerit_dvar = 0

do i = 1, size(s%var)

  if (.not. s%var(i)%useit_opt) cycle
  s%var(i)%dmerit_dvar = 2 * s%var(i)%weight * s%var(i)%delta_merit
  nv = s%var(i)%ix_dvar

  do j = 1, size(s%u)
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
