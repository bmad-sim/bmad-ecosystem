!+
! Subroutine tao_curve_rms_calc (curve, who, rms, mean)
!
! Routine to calculate the mean and rms of the line or symbol points of a curve.
!
! Input:
!   curve     -- tao_curve_struct: Curve to analyze.
!   who       -- character(*): "LINE" or "SYMBOL".
!
! Output:
!   rms        -- real(rp): RMS. -1 => Curve has no data.
!   mean       -- real(rp): Mean.
!-

subroutine tao_curve_rms_calc (curve, who, rms, mean)

use tao_interface, dummy => tao_curve_rms_calc

implicit none

type (tao_curve_struct) curve
real(rp) rms, mean, ys, dx
integer i, n
character(*) who

!

rms = -1
mean = 0

select case (who)
case ('LINE')
  if (.not. allocated(curve%x_line)) return
  n = size(curve%x_line)
  if (curve%x_line(n) == curve%x_line(1)) return
  rms = 0
  do i = 2, n
    ys = 0.5_rp * (curve%y_line(i) + curve%y_line(i-1))
    dx = curve%x_line(i) - curve%x_line(i-1)
    mean = mean + ys * dx
    rms = rms + ys**2 * dx
  enddo

  dx = curve%x_line(n) - curve%x_line(1)
  mean = mean / dx
  rms = sqrt(max(0.0_rp, rms/dx - mean**2))

case ('SYMBOL')
  if (.not. allocated(curve%x_symb)) return
  n = size(curve%x_symb)
  mean = sum(curve%y_symb) / n
  rms = sqrt(sum((curve%y_symb-mean)**2) / n)

end select

end subroutine tao_curve_rms_calc
