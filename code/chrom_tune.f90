!+
! Subroutine chrom_tune (ring, delta_e, target_x, target_y, err_flag)
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

subroutine chrom_tune(ring,delta_e,target_x,target_y,err_flag)
  
  use bmad_struct
  use bmad_interface
  use nr, only: gaussj

  implicit none
  
  type (ring_struct) ring
  type (coord_struct), allocatable :: co_(:)
  type (modes_struct) mode

  integer i, j, n_sex, n_x_sex, n_y_sex
  integer, pointer :: ix_sex(:)
  integer, allocatable :: ix_x_sex(:), ix_y_sex(:)

  real(rp), allocatable :: sex_y_values(:), sex_x_values(:)
  real(rp) target_x, target_y, chrom_x, chrom_y
  real(rp) delta_x, delta_y, d_chrom, chrom_x0
  real(rp) chrom_y0, step_x, step_y, chrom_(2,1), matrix(2,2)
  real(rp) delta_e
 
  logical err_flag, debug
  logical, allocatable :: is_x_sex(:)

!
 
  allocate (co_(0:ring%n_ele_max))

  debug = .false.

  n_sex = 1
  n_x_sex = 1
  n_y_sex = 1

  call elements_locator(sextupole$, ring, ix_sex)

  allocate (is_x_sex(size(ix_sex)))

  do i = 1, size(ix_sex)
     if ((ring%ele_(ix_sex(i))%x%beta) > & 
                      (ring%ele_(ix_sex(i))%y%beta)) then
        is_x_sex(i) = .true.
     else
        is_x_sex(i) = .false.
     end if
  end do

  allocate (ix_x_sex(count(is_x_sex)), ix_y_sex(count(.not. is_x_sex)))
  allocate (sex_x_values(size(ix_x_sex)), sex_y_values(size(ix_y_sex)))
  ix_x_sex = pack(ix_sex, mask = is_x_sex)
  ix_y_sex = pack(ix_sex, mask = .not. is_x_sex)

  delta_x = 0.1
  delta_y = 0.1

!

  do j = 1, 100

    sex_x_values(:) = ring%ele_(ix_x_sex(:))%value(k2$)
    sex_y_values(:) = ring%ele_(ix_y_sex(:))%value(k2$)

    call chrom_calc(ring, delta_e, chrom_x0, chrom_y0)

    if (debug .and. j == 1) then
      print*, 'Target Chromaticities'
      print*, target_x, target_y
      print*, 'Initial Chromaticities'
      print*, chrom_x0, chrom_y0
    end if

    d_chrom = abs(target_x - chrom_x0) + abs(target_y - chrom_y0)

    if (d_chrom < 1.0e-3) then
      if (debug) then
        print *, 'Number of iterations:', j
        print *, 'Final Chromaticities, x & y:'
        print *, chrom_x0, chrom_y0
      end if
      deallocate (ix_x_sex, ix_y_sex, ix_sex)
      deallocate (sex_x_values, sex_y_values, is_x_sex)
      deallocate (co_)
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
  
  print*, 'ERROR IN CHROM_TUNE:  TUNING SEXTUPOLES FAILED!'
  err_flag = .true.

end subroutine
