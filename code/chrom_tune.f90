!+
! Subroutine chrom_tune (ring, delta_e, target_x, target_y, err_tol, err_flag)
!
! Subroutine to set the sextupole strengths so that the ring has the desired
! chormaticities.
! This routine sorts the sextupoles into verticle and horizontal
! groups and varies the two groups to match the desired chromaticities.
! 
! Modules Needed:
!   use bmad
!
! Input:
!   ring     -- Ring_struct: Ring to use, 
!   delta_e  -- Real(rp): Delta energy used for the calculation.
!                    If 0 then default of 1.0e-4 is used.
!   target_x -- Real(rp): Target X Chromaticity
!   target_y -- Real(rp): Target Y Chromaticity
!   err_tol  -- Real(rp): Max allowable Error:
!                    Error = | X_Target - X_Actual | + | Y_Target -Y_Actual |
!               A good number is: err_tol = 0.05_rp
!
! Output:
!   ring     -- Ring_struct: Ring with sextupole set
!   delta_e  -- Real(rp): Set to 1.0e-4 if on input DELTA_E =< 0.
!   err_flag -- Logical: .false. if match successful, .true. if failed
!                   Fails if takes longer than 100 iterations.  If it 
!                   fails the sextupoles are set to the last value
!                   calculated.
!
! Note: This subroutine assumes the twiss parameters have been computed.
!
! Created:  2002/07/24 jtu2
!-

#include "CESR_platform.inc"

subroutine chrom_tune(ring, delta_e, target_x, target_y, err_tol, err_flag)
  
  use bmad_struct
  use bmad_interface
  use nr, only: gaussj

  implicit none
  
  type (ring_struct) ring
  type (modes_struct) mode

  integer i, j, n_sex, n_x_sex, n_y_sex
  integer, pointer :: ix_sex(:)
  integer, allocatable :: ix_x_sex(:), ix_y_sex(:)
  integer, allocatable :: sex_type(:)
  integer ix

  real(rp), allocatable :: sex_y_values(:), sex_x_values(:)
  real(rp) target_x, target_y, chrom_x, chrom_y
  real(rp) delta_x, delta_y, d_chrom, chrom_x0
  real(rp) chrom_y0, step_x, step_y, chrom_(2,1), matrix(2,2)
  real(rp) delta_e, err_tol
 
  logical err_flag, debug

!                      
 
  debug = .false.

  call elements_locator(sextupole$, ring, ix_sex)
  allocate (sex_type(size(ix_sex)))

  do i = 1, size(ix_sex)
    ix = ix_sex(i)
    if (ring%ele_(ix)%value(tilt$) /= 0) then
      sex_type(i) = n_plane$  ! do not use tilted sextupoles 
    elseif ((ring%ele_(ix_sex(i))%x%beta) > & 
                      (ring%ele_(ix_sex(i))%y%beta)) then
      sex_type(i) = x_plane$
    else
      sex_type(i) = y_plane$
    end if
  end do

  n_sex = 0
  n_x_sex = count(sex_type == x_plane$)
  n_y_sex = count(sex_type == y_plane$)

  allocate (ix_x_sex(n_x_sex), ix_y_sex(n_y_sex))
  allocate (sex_x_values(n_x_sex), sex_y_values(n_y_sex))
  ix_x_sex = pack(ix_sex, mask = (sex_type == x_plane$))
  ix_y_sex = pack(ix_sex, mask = (sex_type == y_plane$))

  delta_x = 0.1
  delta_y = 0.1

!
        
  do j = 1, 100

    sex_x_values(:) = ring%ele_(ix_x_sex(:))%value(k2$)
    sex_y_values(:) = ring%ele_(ix_y_sex(:))%value(k2$)

    call chrom_calc(ring, delta_e, chrom_x0, chrom_y0)

    if (debug) then
      print *
      print '(1x, a, 2f10.4)', 'Target Chrom:    ', target_x, target_y
      print '(1x, a, 2f10.4)', 'Chrom Now:       ', chrom_x0, chrom_y0
      print '(1x, a, 2f10.4)', 'Chrom diff (N-T):', &
                                 chrom_x0-target_x, chrom_y0-target_y
    end if

    d_chrom = abs(target_x - chrom_x0) + abs(target_y - chrom_y0)

    if (d_chrom < err_tol) then
      if (debug) then
        print *, 'Number of iterations:', j
        print '(1x, a, 2f10.4)', 'Final Chromaticities:', chrom_x0, chrom_y0
      end if
      deallocate (ix_x_sex, ix_y_sex, ix_sex)
      deallocate (sex_x_values, sex_y_values)
      err_flag = .false.
      return
    end if

    ring%ele_(ix_x_sex(:))%value(k2$) = sex_x_values(:)*(1 + delta_x)
    ring%ele_(ix_y_sex(:))%value(k2$) = sex_y_values(:)
    call chrom_calc(ring, delta_e, chrom_x, chrom_y)
     
    matrix(1,1) = (chrom_x - chrom_x0) / delta_x
    matrix(1,2) = (chrom_y - chrom_y0) / delta_x
     
    ring%ele_(ix_x_sex(:))%value(k2$) = sex_x_values(:)
    ring%ele_(ix_y_sex(:))%value(k2$) = sex_y_values(:)*(1 + delta_y)
    call chrom_calc(ring, delta_e, chrom_x, chrom_y)
     
    matrix(2,1) = (chrom_x - chrom_x0) / delta_y
    matrix(2,2) = (chrom_y - chrom_y0) / delta_y

    chrom_(1,1) = (target_x - chrom_x0)
    chrom_(2,1) = (target_y - chrom_y0)

    call gaussj(matrix, chrom_)
     
    ring%ele_(ix_x_sex(:))%value(k2$) = sex_x_values(:)*(1 + chrom_(1,1))
    ring%ele_(ix_y_sex(:))%value(k2$) = sex_y_values(:)*(1 + chrom_(2,1))
        
  end do

!
  
  print *, 'ERROR IN CHROM_TUNE:  TUNING SEXTUPOLES FAILED!'
  deallocate (ix_x_sex, ix_y_sex, ix_sex)
  deallocate (sex_x_values, sex_y_values)
  err_flag = .true.

end subroutine
