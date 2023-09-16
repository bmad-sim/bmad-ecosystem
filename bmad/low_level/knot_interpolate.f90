!+
! Function knot_interpolate (x_knot, y_knot, x_pt, interpolation, err_flag) result (y_pt)
!
! Routine to interpolate between knot points.
!
! Input:
!   x_knot(:)     -- real(rp): Knot x-values.
!   y_knot(:)     -- real(rp): Knot y-values.
!   x_pt          -- real(rp): Point to evaluate at.
!   interpolation -- integer: Interpolation type. cubic$ or linear$.
!
! Output:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!   y_pt          -- real(rp): Interpolated y-value.
!-

function knot_interpolate (x_knot, y_knot, x_pt, interpolation, err_flag) result (y_pt)

use bmad_interface, dummy => knot_interpolate

implicit none

real(rp) x_knot(:), y_knot(:), x_pt, y_pt
real(rp) f
integer interpolation
integer n, ix
logical err_flag, ok
character(*), parameter :: r_name = 'knot_interpolate'

!

if (size(x_knot) /= size(y_knot)) then
  call out_io (s_error$, r_name, 'KNOT X, Y ARRAY SIZES NOT THE SAME!')
  err_flag = .true.
  return
endif

!

if (interpolation == linear$) then
  ix = bracket_index(x_pt, x_knot, 1, f, restrict = .true.)
  y_pt = y_knot(ix) * (1 - f) + y_knot(ix+1) * f
  err_flag = .false.
elseif (interpolation == cubic$) then
  call spline_akima_interpolate (x_knot, y_knot, x_pt, ok, y_pt)
  err_flag = (.not. ok)
else
  call out_io (s_error$, r_name, 'UNKNOWN INTERPOLATION TYPE: ' // int_str(interpolation))
  err_flag = .true.
endif

end function
