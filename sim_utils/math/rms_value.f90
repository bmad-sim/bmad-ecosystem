!+
! Function rms_value(val_arr, good_val, ave_val) result (rms_val)
!
! Routine to calculate the RMS of an array of values.
!
! Input:
!   val_arr(:)    -- real(rp): Array of reals.
!   good_val(:)   -- logical, optional: If present, only calculate RMS where good_val(i) = True.
!
! Output:
!   ave_val       -- real(rp), optional: average value.
!   rms_val       -- real(rp): RMS value. Set to real_garbage$ if there is a problem.
!-

function rms_value(val_arr, good_val, ave_val) result (rms_val)

use output_mod, dummy => rms_value

implicit none

real(rp) val_arr(:), rms_val
real(rp), optional :: ave_val
real(rp) ave
logical, optional :: good_val(:)
integer n
character(*), parameter :: r_name = 'rms_value'

! The max() function is needed since without it, due to 
! round-off error, the sqrt radical may be negative.

rms_val = real_garbage$

if (present(good_val)) then
  n = count(good_val)
  if (n == 0) return
  if (n /= size(val_arr)) then
    call out_io (s_error$, r_name, 'ARRAY SIZE MIS-MATCH! PLEASE REPORT THIS!')
    return
  endif
  ave = sum(val_arr, mask = good_val) / n
  rms_val = sqrt(max(0.0_rp, sum(val_arr**2, mask = good_val) / n - ave**2))
else
  n = size(val_arr)
  ave = sum(val_arr) / n
  rms_val = sqrt(max(0.0_rp, sum(val_arr**2) / n - ave**2))
endif

if (present(ave_val)) ave_val = ave

end function
