program example_plot

use quick_plot
integer id
character(1) ans

! Generate PS and X-windows data plots.

call qp_open_page ('PS-L')  ! Tell Quick Plot to generate a PS file.
call draw_histogram         ! Generate the plot
call qp_close_page          ! quick_plot.ps is the file name
call qp_open_page ('X', id, 400.0_rp, 270.0_rp, 'POINTS')
call draw_histogram
write (*, '(a)', advance = 'NO') ' Hit any key to end program: '
read (*, '(a)') ans

!----------------------------------------------------------------------
! This generates the graphs

contains

subroutine draw_graphs

real(rp), allocatable :: x(:), y(:), z(:), t(:)
real(rp) x_axis_min, x_axis_max, y_axis_min, y_axis_max
integer x_places, x_divisions, y_places, y_divisions
character(80) title
logical err_flag
namelist / parameters / title

! Read in the data
open (1, file = 'plot.dat', status = 'old')
read (1, nml = parameters)                  ! read in the parameters.
call qp_read_data (1, err_flag, x, 1, y, 3, z, 4, t, 5) ! read in the data.
close (1)

! Setup the margins and page border and draw the title
call qp_set_page_border (0.01_rp, 0.02_rp, 0.2_rp, 0.2_rp, '%PAGE')
call qp_set_margin (0.07_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
call qp_draw_text (title, 0.5_rp, 0.85_rp, '%PAGE', 'CT') 

! draw the left graph
call qp_set_box (1, 1, 2, 1)
call qp_calc_and_set_axis ('X', minval(x), maxval(x), 4, 8, 'ZERO_AT_END')
call qp_calc_and_set_axis ('Y', minval(z), maxval(z), 4, 8, 'GENERAL')
call qp_draw_axes ("X\dlab\u", "\gb(\A)")
call qp_draw_data (x, y, symbol_every = 0)

call qp_save_state (.true.)
call qp_set_symbol_attrib (times$, color = blue$, height = 20.0_rp)
call qp_set_line_attrib ('PLOT', color = blue$, style = dashed$)
call qp_draw_data (x, z, symbol_every = 5)
call qp_restore_state

! draw the right graph
call qp_save_state (.true.)
call qp_set_box (2, 1, 2, 1)
call qp_set_graph_attrib (draw_grid = .false.)
call qp_set_symbol_attrib (star5_filled$, height = 10.0_rp)
call qp_set_axis ('Y', -0.1_rp, 0.1_rp, 4, 2)
call qp_set_axis ('Y2', 1.0_rp, 100.0_rp, label = 'Y2 axis', &
                              draw_numbers = .true., ax_type = 'LOG')
call qp_draw_axes ("\m1 \m2 \m3 \m4 \m5 \m6 \m7", "\fsLY\fn", &
                                            title = "That Darn Graph")
call qp_draw_data (x, t, draw_line = .false., symbol_every = 4)
call qp_restore_state

end subroutine

!----------------------------------------------------------------------
! This generates a histogram

subroutine draw_histogram

real(rp) :: x(3) = (/ 20, 30, 40 /)
real(rp) :: y(3) = (/ 34, 72, 16 /)

!

call qp_set_page_border (0.01_rp, 0.02_rp, 0.2_rp, 0.2_rp, '%PAGE')
call qp_set_margin (0.07_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

call qp_set_box (1, 1, 2, 1)
call qp_set_axis ('X', 0.0_rp, 60.0_rp, 6)
call qp_set_axis ('Y', 0.0_rp, 100.0_rp, 5)
call qp_draw_axes ("x", "y", "Line Histogram", .false.)
call qp_draw_histogram (x, y)

call qp_set_box (2, 1, 2, 1)
call qp_draw_axes ("x", "y", "Filled Histogram")
call qp_draw_histogram (x, y, blue$, hatched$)

end subroutine

end program
