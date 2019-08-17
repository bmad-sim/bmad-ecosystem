module synrad_plot_mod

use synrad_mod
use quick_plot

implicit none

!

type plot_param_struct
  real(rp) :: window_width = 800.0_rp, window_height = 400.0_rp
end type

type (plot_param_struct), save :: plot_param

contains

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine synrad_plot_sx (walls)
!
! Routine to interactively plot (x, s) points of the wall.
! Note: This routine never returns to the main program.
!
! Input:
!   walls      -- walls_struct: Wall structure.
!-

subroutine synrad_plot_sx (walls)

use input_mod

type (walls_struct), target :: walls
type (wall_struct), pointer :: minus_side, plus_side
type (wall_pt_struct), pointer :: wpt_p(:), wpt_m(:)

real(rp) s_min, s_max, x_min, x_max
integer ix, ios, i_chan
logical s_user_good, x_user_good

character(40) ans

!

call qp_open_page ('X', i_chan, plot_param%window_width, plot_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

s_user_good = .false.
x_user_good = .false.
x_max = 100

plus_side => walls%positive_x_wall
minus_side  => walls%negative_x_wall


! Loop

do

  ! Determine s min/max

  if (.not. s_user_good) then
    s_min = plus_side%pt(0)%s
    s_max = walls%s_max
  endif

  ! Draw axes

  call qp_clear_page

  if (.not. x_user_good) then
    wpt_p => plus_side%pt(0:plus_side%n_pt_max)
    wpt_m => minus_side%pt(0:minus_side%n_pt_max)

    x_min = 1010 * min(minval(pack (wpt_m%x, wpt_m%s >= s_min .and. wpt_m%s <= s_max)), & 
                       minval(pack (wpt_p%x, wpt_p%s >= s_min .and. wpt_p%s <= s_max)))

    x_max = 1010 * max(maxval(pack (wpt_m%x, wpt_m%s >= s_min .and. wpt_m%s <= s_max)), & 
                       maxval(pack (wpt_p%x, wpt_p%s >= s_min .and. wpt_p%s <= s_max)))
  endif

  call qp_calc_and_set_axis ('X', s_min, s_max, 10, 16, 'GENERAL')
  call qp_calc_and_set_axis ('Y', x_min, 1.01*x_max, 6, 10, 'GENERAL')
  call qp_set_margin (0.07_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
  call qp_draw_axes ('S (m)', 'X (mm)')

  ! Draw walls

  call draw_a_wall (minus_side)
  call draw_a_wall (plus_side)

  ! Query

  print *, 'Syntax: "x" or "s" followed by <min> <max> values.'
  print *, '[<min> = "auto" --> autoscale] Example: "x auto", "s 10 60"'
  call read_a_line ('Input: ', ans)

  call string_trim (ans, ans, ix)
  if (ans(1:2) == 's ') then
    call string_trim(ans(2:), ans, ix)
    if (ans == 'auto') then
      s_user_good = .false.
    else
      read (ans, *, iostat = ios) s_min, s_max
      if (ios /= 0) then
        print *, 'CANNOT DECODE MIN/MAX VALUES'
      else
        s_user_good = .true.
      endif
    endif

  elseif (ans(1:2) == 'x ') then
    call string_trim(ans(2:), ans, ix)
    if (ans == 'auto') then
      x_user_good = .false.
    else
      read (ans, *, iostat = ios) x_min, x_max
      if (ios /= 0) then
        print *, 'CANNOT DECODE MIN/MAX VALUES'
      else
        x_user_good = .true.
      endif
    endif

  else
    print *, 'I DO NOT UNDERSTAND THIS...'
  endif

enddo

!--------------------------------------------------------------------------------------------------
contains

subroutine draw_a_wall (wall)

type (wall_struct), target :: wall
type (wall_pt_struct), pointer :: pt0, pt1
integer i, n, nn

!

do i = 1, wall%n_pt_max
  pt0 => wall%pt(i-1)
  pt1 => wall%pt(i)
  if (pt1%s < s_min .or. pt0%s > s_max) cycle
  if (pt1%phantom) then
    call qp_draw_line(pt0%s, pt1%s, 1000*pt0%x, 1000*pt1%x, color = 'red', line_pattern = 'dashed')
    print '(a, 5(f10.3, a))', 'Phantom on ' // trim(wall_name(wall%side)) // ' between (s, x) points: (', &
            pt0%s, ',', pt0%x, '),  (', pt1%s, ',', pt1%x, ')'
  else
    call qp_draw_line(pt0%s, pt1%s, 1000*pt0%x, 1000*pt1%x)
  endif
enddo

end subroutine draw_a_wall

end subroutine synrad_plot_sx

end module
