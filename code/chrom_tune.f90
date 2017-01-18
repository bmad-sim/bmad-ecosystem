!+
! Subroutine chrom_tune (lat, delta_e, target_x, target_y, err_tol, err_flag)
!
! Subroutine to set the sextupole strengths so that the lattice has the desired chormaticities.
! This routine first sorts the sextupoles into vertical and horizontal groups using the 
! relative vertical and horizontal betas at the sextupoles.
! The two groups are then varied by fractional amounts to match the desired chromaticities.
! 
! Modules Needed:
!   use bmad
!
! Input:
!   lat     -- lat_struct: Lat to use, 
!   delta_e  -- Real(rp): Delta energy used for the calculation.
!                    If 0 then default of 1.0d-4 is used.
!   target_x -- Real(rp): Target X Chromaticity
!   target_y -- Real(rp): Target Y Chromaticity
!   err_tol  -- Real(rp): Max allowable Error:
!                    Error = | X_Target - X_Actual | + | Y_Target -Y_Actual |
!               A good number is: err_tol = 0.05_rp
!
! Output:
!   lat     -- lat_struct: Lat with sextupole set
!   delta_e  -- Real(rp): Set to 1.0d-4 if on input DELTA_E =< 0.
!   err_flag -- Logical: .false. if match successful, .true. if failed
!                   Fails if takes longer than 100 iterations.
!                   If it fails the sextupoles are set to the last value calculated.
!
! Note: This subroutine assumes the Twiss parameters have been computed.
!
! Created:  2002/07/24 jtu2
!-

subroutine chrom_tune(lat, delta_e, target_x, target_y, err_tol, err_flag)
  
use bmad_interface, except_dummy => chrom_tune
use nr, only: gaussj

implicit none

type (lat_struct) lat
type (ele_pointer_struct), allocatable :: ele_sex(:)

integer i, j, ix, iy, n_sex
integer, allocatable :: ix_x_sex(:), ix_y_sex(:)

real(rp), allocatable :: sex_y_value(:), sex_x_value(:)
real(rp) target_x, target_y, chrom_x, chrom_y
real(rp) delta_x, delta_y, d_chrom, chrom_x0
real(rp) chrom_y0, chrom(2,1), matrix(2,2)
real(rp) delta_e, err_tol
 
logical err_flag
character(*), parameter :: r_name = 'chrom_tune'

!                      
 
err_flag = .true.
call lat_ele_locator('sextupole::*', lat, ele_sex, n_sex, err_flag)

allocate (ix_x_sex(n_sex), ix_y_sex(n_sex))
allocate (sex_x_value(n_sex), sex_y_value(n_sex))

ix = 0; iy = 0

do i = 1, n_sex
  if (ele_sex(i)%ele%value(tilt_tot$) /= 0) cycle    ! do not use tilted sextupoles 
  if (ele_sex(i)%ele%a%beta > ele_sex(i)%ele%b%beta) then
    ix = ix + 1
    ix_x_sex(ix) = ele_sex(i)%ele%ix_ele
  else
    iy = iy + 1
    ix_y_sex(iy) = ele_sex(i)%ele%ix_ele
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

delta_x = 0.1
delta_y = 0.1

!
      
do j = 1, 100

  sex_x_value(:) = lat%ele(ix_x_sex(:))%value(k2$)
  sex_y_value(:) = lat%ele(ix_y_sex(:))%value(k2$)

  call chrom_calc(lat, delta_e, chrom_x0, chrom_y0)

  d_chrom = abs(target_x - chrom_x0) + abs(target_y - chrom_y0)

  if (d_chrom < err_tol) then
    err_flag = .false.
    return
  end if

  lat%ele(ix_x_sex(:))%value(k2$) = sex_x_value(:)*(1 + delta_x)
  lat%ele(ix_y_sex(:))%value(k2$) = sex_y_value(:)
  call chrom_calc(lat, delta_e, chrom_x, chrom_y)
   
  matrix(1,1) = (chrom_x - chrom_x0) / delta_x
  matrix(1,2) = (chrom_y - chrom_y0) / delta_x
   
  lat%ele(ix_x_sex(:))%value(k2$) = sex_x_value(:)
  lat%ele(ix_y_sex(:))%value(k2$) = sex_y_value(:)*(1 + delta_y)
  call chrom_calc(lat, delta_e, chrom_x, chrom_y)
   
  matrix(2,1) = (chrom_x - chrom_x0) / delta_y
  matrix(2,2) = (chrom_y - chrom_y0) / delta_y

  chrom(1,1) = (target_x - chrom_x0)
  chrom(2,1) = (target_y - chrom_y0)

  call gaussj(matrix, chrom)
   
  lat%ele(ix_x_sex(:))%value(k2$) = sex_x_value(:)*(1 + chrom(1,1))
  lat%ele(ix_y_sex(:))%value(k2$) = sex_y_value(:)*(1 + chrom(2,1))
      
end do

!

call out_io (s_error$, r_name, 'CANNOT ADJUST SEXTUPOLES TO GET DESIRED CHROMATICITIES!')
err_flag = .true.

end subroutine
