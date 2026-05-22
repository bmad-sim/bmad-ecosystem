program css4_color_swatch
!
! Generate a color swatch image using quick_plot with all CSS4 colors.
! Output: css4_qp_swatch.png
!

use quick_plot

implicit none

integer :: i, ix_color, row, col
integer, parameter :: n_colors = 154   ! indices 0-15 (original) + 17-154 (CSS4)
integer, parameter :: n_cols = 10
integer, parameter :: n_rows = 16      ! ceil(154/10) = 16, but we have gaps
real(rp) :: x0, y0, x1, y1
real(rp) :: swatch_w, swatch_h
real(rp), parameter :: page_w = 6.0_rp   ! inches
real(rp), parameter :: page_h = 5.0_rp   ! inches
character(24) :: color_name
logical :: error

! Open a PNG output page
call qp_open_page ('GIF', page_w, page_h, 'css4_qp_swatch.png')
call qp_set_page_border (0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%PAGE')
call qp_set_margin (0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%PAGE')
call qp_set_box (1, 1, 1, 1)
call qp_set_graph_position (0.0_rp, page_w, 0.0_rp, page_h, 'INCH')
call qp_set_graph_limits (0.0_rp, real(n_cols, rp), 0.0_rp, 16.0_rp)
call qp_eliminate_xy_distortion ()

swatch_w = 1.0_rp
swatch_h = 1.0_rp

! Draw original 16 colors (indices 0-15)
do i = 0, 15
  col = mod(i, n_cols)
  row = 15 - i / n_cols   ! top to bottom

  x0 = real(col, rp) * swatch_w
  y0 = real(row, rp) * swatch_h
  x1 = x0 + swatch_w
  y1 = y0 + swatch_h

  color_name = qp_enum_to_string(i, 'color')
  call qp_draw_rectangle (x0, x1, y0, y1, color=color_name, fill_color=color_name)
enddo

! Draw CSS4 colors (indices 17-154)
do i = 17, qp_max_color_index$
  ix_color = i - 17 + 16   ! sequential position after original 16
  col = mod(ix_color, n_cols)
  row = 15 - ix_color / n_cols

  x0 = real(col, rp) * swatch_w
  y0 = real(row, rp) * swatch_h
  x1 = x0 + swatch_w
  y1 = y0 + swatch_h

  color_name = qp_enum_to_string(i, 'color')
  call qp_draw_rectangle (x0, x1, y0, y1, color=color_name, fill_color=color_name)
enddo

call qp_close_page()

print *, 'Color swatch written to: css4_qp_swatch.png'

end program css4_color_swatch
