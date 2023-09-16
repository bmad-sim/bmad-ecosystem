program example_plot

use quick_plot
integer i, id, id2
real(rp) xlen, ylen
character(10) ans

! Which plot?

write (*, '(a)') 'What to draw?:'
write (*, '(a)') ' 1) Histogram example'
write (*, '(a)') ' 2) Graph example'
write (*, '(a)') ' 3) Symbol chart'
write (*, '(a)', advance = 'NO') ' Input number <CR=3> '
read (*, '(a)') ans
call string_trim (ans, ans, i)

select case (ans)
case('1')
  xlen = 800
  ylen = 540
case('2')
  xlen = 800
  ylen = 540
case default
  xlen = 600
  ylen = 800
end select

! Generate PS and X-windows data plots.

call qp_open_page ('X', id, xlen, ylen, 'POINTS')
call draw_it
call qp_open_page ('PS-L')  ! Tell Quick Plot to generate a PS file.
call draw_it
call qp_close_page          ! quick_plot.ps is the file name
call qp_open_page ('GIF-L')  ! Tell Quick Plot to generate a GIF file.
call draw_it
call qp_close_page          ! quick_plot.ps is the file name
write (*, '(a)', advance = 'NO') ' Hit any key to end program: '
read (*, '(a)') ans

!----------------------------------------------------------------------
contains

subroutine draw_it
  select case (ans)
  case('1')
    call draw_histogram
  case('2')
    call draw_graphs
  case default
    call draw_symbols
  end select
end subroutine

!----------------------------------------------------------------------
! contains

subroutine draw_symbols

integer i, j, k
character(3) num
!

call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

do i = 0, 4
  do j = 0, 7
    k = i + 5 * j
    if (k > 31) k = 31 - k
    call qp_set_box (i+1, 8-j, 5, 8)
    call qp_draw_rectangle (0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp, '%BOX')
    write (num, '(i0)') k
    call qp_draw_text (num, 0.2_rp, 0.7_rp, '%BOX')
    call qp_draw_symbol (0.5_rp, 0.5_rp, '%BOX', qp_symbol_type_name(k), 40.0_rp) 
  enddo
enddo

end subroutine

!----------------------------------------------------------------------
! contains
! This generates the graphs


subroutine draw_graphs

real(rp), allocatable :: x(:), y(:), z(:), t(:)
real(rp) x_axis_min, x_axis_max, y_axis_min, y_axis_max
integer x_places, x_divisions, y_places, y_divisions
character(80) title
logical err_flag
namelist / parameters / title

! Create the data

open (1, file = 'plot.dat')
write (1, '(a)') ' &parameters'
write (1, '(a)') '   title = "A Tale of Two Graphs"'
write (1, '(a)') ' /'
write (1, '(a)') ''
write (1, '(a)') ' Any junk here...'
write (1, '(a)') ''
write (1, '(a)') ' Col1      Col2      Col3      Col4      Col5'
write (1, '(a)') '    0    0.0000    0.1000    0.0000   -0.0125'
write (1, '(a)') '    4    0.0016    0.0921    0.0405   -0.0135'
write (1, '(a)') '    8    0.0064    0.0697    0.0775   -0.0144'
write (1, '(a)') '   12    0.0144    0.0362    0.1044   -0.0153'
write (1, '(a)') '   16    0.0256   -0.0029    0.1160   -0.0161'
write (1, '(a)') '   20    0.0400   -0.0416    0.1091   -0.0167'
write (1, '(a)') '   24    0.0576   -0.0737    0.0838   -0.0171'
write (1, '(a)') '   28    0.0784   -0.0942    0.0429   -0.0173'
write (1, '(a)') '   32    0.1024   -0.0998   -0.0077   -0.0172'
write (1, '(a)') '   36    0.1296   -0.0897   -0.0602   -0.0168'
write (1, '(a)') '   40    0.1600   -0.0654   -0.1060   -0.0161'
write (1, '(a)') '   44    0.1936   -0.0307   -0.1370   -0.0150'
write (1, '(a)') '   48    0.2304    0.0087   -0.1474   -0.0134'
write (1, '(a)') '   52    0.2704    0.0469   -0.1343   -0.0114'
write (1, '(a)') '   56    0.3136    0.0776   -0.0985   -0.0089'
write (1, '(a)') '   60    0.3600    0.0960   -0.0447   -0.0059'
write (1, '(a)') '   64    0.4096    0.0993    0.0191   -0.0023'
write (1, '(a)') '   68    0.4624    0.0869    0.0830    0.0019'
write (1, '(a)') '   72    0.5184    0.0608    0.1365    0.0068'
write (1, '(a)') '   76    0.5776    0.0251    0.1704    0.0124'
write (1, '(a)') '   80    0.6400   -0.0146    0.1781    0.0187'
write (1, '(a)') '   84    0.7056   -0.0519    0.1572    0.0258'
write (1, '(a)') '   88    0.7744   -0.0811    0.1100    0.0336'
write (1, '(a)') '   92    0.8464   -0.0975    0.0428    0.0424'
write (1, '(a)') '   96    0.9216   -0.0985   -0.0342    0.0520'
write (1, '(a)') '  100    1.0000   -0.0839   -0.1088    0.0625'
close (1)

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
call qp_set_symbol_attrib ('times', color = 'blue', height = 20.0_rp)
call qp_set_line_attrib ('PLOT', color = 'blue', pattern = 'dashed')
call qp_draw_data (x, z, symbol_every = 5)
call qp_restore_state

! draw the right graph
call qp_save_state (.true.)
call qp_set_box (2, 1, 2, 1)
call qp_set_graph_attrib (draw_grid = .false.)
call qp_set_symbol_attrib ('star5_filled', height = 10.0_rp)
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
call qp_draw_axes ("x", "y", "Line Histogram")
call qp_draw_histogram (x, y, 'transparent')

call qp_set_box (2, 1, 2, 1)
call qp_draw_axes ("x", "y", "Filled Histogram", .false.)
call qp_set_line_attrib ('STD', width = 10)
call qp_draw_histogram (x, y, 'blue', 'solid_fill', 'red')

end subroutine

end program
