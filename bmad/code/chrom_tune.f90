!+
! Subroutine chrom_tune (lat, delta_e, target_x, target_y, err_tol, err_flag)
!
! Subroutine to set the sextupole strengths so that the lattice has the desired chormaticities.
! This routine first sorts the sextupoles into vertical and horizontal groups using the 
! relative vertical and horizontal betas at the sextupoles.
! The two groups are then varied by fractional amounts to match the desired chromaticities.
! 
! Input:
!   lat      -- lat_struct: Lat to use, 
!   delta_e  -- Real(rp): Delta energy used for the calculation.
!                    If 0 then default of 1.0d-4 is used.
!   target_x -- Real(rp): Target X Chromaticity
!   target_y -- Real(rp): Target Y Chromaticity
!   err_tol  -- Real(rp): Max allowable Error:
!                    Error = | X_Target - X_Actual | + | Y_Target -Y_Actual |
!               A good number is: err_tol = 0.05_rp
!
! Output:
!   lat      -- lat_struct: Lat with sextupole set
!   delta_e  -- Real(rp): Set to 1.0d-4 if on input DELTA_E =< 0.
!   err_flag -- Logical: .false. if match successful, .true. if failed
!                   Fails if takes longer than 100 iterations.
!                   If it fails the sextupoles are set to the last value calculated.
!
! Note: This subroutine assumes the Twiss parameters have been computed.
!-

subroutine chrom_tune(lat, delta_e, target_x, target_y, err_tol, err_flag)

use bmad_interface, except_dummy => chrom_tune
use super_recipes_mod

implicit none

type (lat_struct) lat
type (ele_pointer_struct), allocatable, target :: ele_sex(:)
type (ele_struct), pointer :: ele
type (super_mrqmin_storage_struct) storage

integer i, j, ix, iy, n_sex, status
integer, allocatable :: ix_x_sex(:), ix_y_sex(:)

real(rp), allocatable :: sex_y_value(:), sex_x_value(:)
real(rp) target_x, target_y, chisq
real(rp) d_chrom, chrom_x0, chrom_y0, a_lambda
real(rp) delta_e, err_tol, k2_vec(2), weight(2), chrom_target(2)
 
logical err_flag, maska(2), all_x_zero, all_y_zero
character(*), parameter :: r_name = 'chrom_tune'

!                      
 
err_flag = .true.
call lat_ele_locator('sextupole::*', lat, ele_sex, n_sex, err_flag)

allocate (ix_x_sex(n_sex), ix_y_sex(n_sex))
allocate (sex_x_value(n_sex), sex_y_value(n_sex))

ix = 0; iy = 0
all_x_zero = .true.;  all_y_zero = .true.

do i = 1, n_sex
  ele => ele_sex(i)%ele
  if (ele%value(tilt_tot$) /= 0) cycle    ! do not use tilted sextupoles
  if (ele%a%beta == 0 .or. ele%b%beta == 0) then
    call out_io (s_error$, r_name, 'TWISS PARAMETERS NOT COMPUTED AT: ' // ele%name, 'NO TUNING DONE.')
    return
  endif

  if (ele%a%beta > ele%b%beta) then
    ix = ix + 1
    ix_x_sex(ix) = ele%ix_ele
    if (ele%value(k2$) /= 0) all_x_zero = .false.

  else
    iy = iy + 1
    ix_y_sex(iy) = ele%ix_ele
    if (ele%value(k2$) /= 0) all_y_zero = .false.
  end if
end do

if (ix == 0) then
  call out_io (s_error$, r_name, 'NUMBER OF SEXTUPOLES WHERE BETA_A > BETA_B IS ZERO.', &
                                 'CHROMATICITY TUNING ALGORITHM WILL NOT WORK HERE.', 'NOTHING DONE')
  return
endif

if (iy == 0) then
  call out_io (s_error$, r_name, 'NUMBER OF SEXTUPOLES WHERE BETA_A < BETA_B IS ZERO.', &
                                 'CHROMATICITY TUNING ALGORITHM WILL NOT WORK HERE.', 'NOTHING DONE')
  return
endif

call re_allocate (ix_x_sex, ix);  call re_allocate (sex_x_value, ix)
call re_allocate (ix_y_sex, iy);  call re_allocate (sex_y_value, iy)

!
      
sex_x_value(:) = lat%ele(ix_x_sex(:))%value(k2$)
sex_y_value(:) = lat%ele(ix_y_sex(:))%value(k2$)
k2_vec = 0
maska = .false.
weight = 1
chrom_target = [target_x, target_y]
a_lambda = -1

do j = 1, 100
  call super_mrqmin (chrom_target, weight, k2_vec, chisq, chrom_func, storage, a_lambda, status)

  if (status < 0) then
    call out_io (s_error$, r_name, 'SINGULAR MATRIX ENCOUNTERED!')
    call this_bookkeeping()
    return
  endif

  d_chrom = abs(target_x - chrom_x0) + abs(target_y - chrom_y0)

  if (d_chrom < err_tol) then
    err_flag = .false.
    if (all_x_zero) then
      lat%ele(ix_x_sex(:))%value(k2$) = k2_vec(1)
    else
      lat%ele(ix_x_sex(:))%value(k2$) = sex_x_value(:) * (1 + k2_vec(1))
    endif
    if (all_y_zero) then
      lat%ele(ix_y_sex(:))%value(k2$) = k2_vec(2)
    else
      lat%ele(ix_y_sex(:))%value(k2$) = sex_y_value(:) * (1 + k2_vec(2))
    endif
    call this_bookkeeping()
    return
  endif
enddo

call out_io (s_error$, r_name, 'CANNOT ADJUST SEXTUPOLES TO GET DESIRED CHROMATICITIES!')
err_flag = .true.

!-------------------------------------------------------
contains

subroutine this_bookkeeping()

do i = 1, size(ix_x_sex)
  ele => lat%ele(ix_x_sex(i))
  call set_flags_for_changed_attribute(ele, ele%value(k2$))
enddo

do i = 1, size(ix_y_sex)
  ele => lat%ele(ix_y_sex(i))
  call set_flags_for_changed_attribute(ele, ele%value(k2$))
enddo

end subroutine this_bookkeeping

!------------------------------------------------------
! contains

subroutine chrom_func (a_try, y_fit, dy_da, status)

real(rp), intent(in) :: a_try(:)
real(rp), intent(out) :: y_fit(:)
real(rp), intent(out) :: dy_da(:, :)

real(rp) :: delta = 0.01
real(rp) chrom_x, chrom_y
integer status, i

!

if (all_x_zero) then
  lat%ele(ix_x_sex(:))%value(k2$) = a_try(1)
else
  lat%ele(ix_x_sex(:))%value(k2$) = sex_x_value(:) * (1 + a_try(1))
endif

if (all_y_zero) then
  lat%ele(ix_y_sex(:))%value(k2$) = a_try(2)
else
  lat%ele(ix_y_sex(:))%value(k2$) = sex_y_value(:) * (1 + a_try(2))
endif

call this_bookkeeping()
call chrom_calc(lat, delta_e, chrom_x0, chrom_y0)

y_fit = [chrom_x0, chrom_y0]

if (all_x_zero) then
  lat%ele(ix_x_sex(:))%value(k2$) = a_try(1) + delta
else
  lat%ele(ix_x_sex(:))%value(k2$) = sex_x_value(:) * (1 + a_try(1) + delta)
endif

if (all_y_zero) then
  lat%ele(ix_y_sex(:))%value(k2$) = a_try(2)
else
  lat%ele(ix_y_sex(:))%value(k2$) = sex_y_value(:) * (1 + a_try(2))
endif

call this_bookkeeping()
call chrom_calc(lat, delta_e, chrom_x, chrom_y)

dy_da(1,1) = (chrom_x - chrom_x0) / delta
dy_da(2,1) = (chrom_y - chrom_y0) / delta

if (all_x_zero) then
  lat%ele(ix_x_sex(:))%value(k2$) = a_try(1)
else
  lat%ele(ix_x_sex(:))%value(k2$) = sex_x_value(:) * (1 + a_try(1))
endif

if (all_y_zero) then
  lat%ele(ix_y_sex(:))%value(k2$) = a_try(2) + delta
else
  lat%ele(ix_y_sex(:))%value(k2$) = sex_y_value(:) * (1 + a_try(2) + delta)
endif

call this_bookkeeping()
call chrom_calc(lat, delta_e, chrom_x, chrom_y)

dy_da(1,2) = (chrom_x - chrom_x0) / delta
dy_da(2,2) = (chrom_y - chrom_y0) / delta

end subroutine chrom_func

end subroutine
