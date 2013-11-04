module lux_plot_mod

use lux_module
use quick_plot
use input_mod

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine lux_plot_tiles (lux_param, source)
!
! Routine to plot the tiling of the photon emission direction
!-

subroutine lux_plot_tiles (lux_param, source)

implicit none

type (lux_source_struct) source
type (lux_params_struct) lux_param

real(rp) x1, x2, y1, y2
integer ix, i_chan
character(16) ans

!

call qp_open_page ('X', i_chan, lux_param%window_width, lux_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
call qp_set_margin (0.07_rp, 0.05_rp, 0.10_rp, 0.05_rp, '%PAGE')

call qp_calc_and_set_axis ('X', -lux_param%phi_max, lux_param%phi_max, 6, 10, 'GENERAL')
call qp_calc_and_set_axis ('Y', -lux_param%y_max, lux_param%y_max, 6, 10, 'GENERAL')

call qp_draw_axes ('Phi', 'Y', 'Emission Direction Tiles')

do ix = 1, size(source%tile)
  x1 = (source%tile(ix)%iphi - 0.5) * lux_param%del_phi
  x2 = (source%tile(ix)%iphi + 0.5) * lux_param%del_phi
  y1 = (source%tile(ix)%iy - 0.5) * lux_param%del_y
  y2 = (source%tile(ix)%iy + 0.5) * lux_param%del_y
  call qp_paint_rectangle (x1, x2, y1, y2, color = red$)
enddo

call read_a_line ('CR to exit: ', ans)

end subroutine lux_plot_tiles

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine lux_plot_detector (lux_param, detector)
!
! Routine to plot the intensity at the detector pixels
!-

subroutine lux_plot_detector (lat, lux_param, detector)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (surface_grid_struct) detector
type (lux_params_struct) lux_param

real(rp) x1, x2, y1, y2
integer ix, iy, i_chan
character(16) ans

!

call qp_open_page ('X', i_chan, lux_param%window_width, lux_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
call qp_set_margin (0.07_rp, 0.05_rp, 0.10_rp, 0.05_rp, '%PAGE')

ele => lat%ele(lat%n_ele_track)
call qp_calc_and_set_axis ('X', -ele%value(x1_limit$), ele%value(x2_limit$), 6, 10, 'GENERAL')
call qp_calc_and_set_axis ('Y', -ele%value(y1_limit$), ele%value(y2_limit$), 6, 10, 'GENERAL')

call qp_draw_axes ('X', 'Y', 'Detector Pixels')

do ix = 1, lbound(detector%pt, 1), ubound(detector%pt, 1)
do iy = 1, lbound(detector%pt, 2), ubound(detector%pt, 2)
  call qp_paint_rectangle (x1, x2, y1, y2, color = red$)
enddo
enddo

call read_a_line ('CR to exit: ', ans)

end subroutine lux_plot_detector

end module
