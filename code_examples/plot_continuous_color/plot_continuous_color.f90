program plot_continuous_color

use quick_plot

implicit none

integer :: i, id, ix_color
integer, parameter :: n_bars = 512
real(rp) :: t, x1_inch, x2_inch
real(rp), parameter :: page_w = 396.0_rp   ! 5.5 inches in points
real(rp), parameter :: page_h = 54.0_rp    ! 0.75 inch in points
! Graph placement in inches (within the page)
real(rp), parameter :: gx1 = 0.25_rp       ! left margin
real(rp), parameter :: gxlen = 5.0_rp      ! graph width
real(rp), parameter :: gy1 = 0.28_rp       ! bottom margin (room for tick numbers)
real(rp), parameter :: gylen = 0.35_rp     ! graph height

! Open PostScript page
call qp_open_page ('PS-L', id, page_w, page_h, 'POINTS')

call qp_set_page_border (0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%PAGE')
call qp_set_graph_placement (gx1, gxlen, gy1, gylen, 'INCH')

! Configure axes: x from 0 to 1, no y-axis
call qp_set_axis ('X', 0.0_rp, 1.0_rp, div = 5, places = 1)
call qp_set_axis ('Y', 0.0_rp, 1.0_rp, draw_numbers = .false., draw_label = .false., &
                  major_tick_len = 0.0_rp, minor_tick_len = 0.0_rp)
call qp_set_axis ('X2', mirror = .false., draw_numbers = .false., &
                  major_tick_len = 0.0_rp, minor_tick_len = 0.0_rp)
call qp_set_axis ('Y2', mirror = .false., draw_numbers = .false., &
                  major_tick_len = 0.0_rp, minor_tick_len = 0.0_rp)

! Draw the continuous color gradient using basic (inch-coordinate) rectangles
do i = 0, n_bars - 1
  t = real(i, rp) / real(n_bars - 1, rp)
  ix_color = qp_continuous_color(t)
  x1_inch = gx1 + gxlen * real(i, rp) / real(n_bars, rp)
  x2_inch = gx1 + gxlen * real(i + 1, rp) / real(n_bars, rp)
  call qp_paint_rectangle_basic (x1_inch, x2_inch, gy1, gy1 + gylen, ix_color, solid_fill$)
end do

! Draw the x-axis with tick marks
call qp_draw_axes (draw_grid = .false.)

call qp_close_page

print '(a)', 'Wrote: quick_plot.ps'

end program plot_continuous_color
