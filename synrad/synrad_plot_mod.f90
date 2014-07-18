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
type (wall_struct), pointer :: inside, outside

real(rp) s_min, s_max, x_min, x_max
real(rp), allocatable :: s_in(:), x_in(:), s_out(:), x_out(:)
integer ix, ios, i_chan
logical s_user_good, x_user_good

character(40) ans

!

call qp_open_page ('X', i_chan, plot_param%window_width, plot_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

s_user_good = .false.
x_user_good = .false.
x_max = 100

outside => walls%positive_x_wall
inside  => walls%negative_x_wall


! Loop

do

  ! Determine s min/max

  if (.not. s_user_good) then
    s_min = outside%pt(0)%s
    s_max = walls%s_max
  endif

  call qp_calc_and_set_axis ('X', s_min, s_max, 10, 16, 'GENERAL')

  ! Get data points

  call get_data_points (inside, s_in, x_in)
  call get_data_points (outside, s_out, x_out)

  ! Now plot

  call qp_clear_page
  if (.not. x_user_good) then
    x_min = 0.99 * min(minval(x_in), minval(x_out))
    x_max = 1.01 * max(maxval(x_in), maxval(x_out))
  endif

  call qp_calc_and_set_axis ('Y', x_min, x_max, 6, 10, 'GENERAL')
  call qp_set_margin (0.07_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
  call qp_draw_graph (s_in, x_in, 'S (m)', 'X (mm)', '', .true., 0)
  call qp_draw_polyline (s_out, x_out)

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

subroutine get_data_points (wall, s, x)

type (wall_struct) wall
real(rp), allocatable :: s(:), x(:)
integer i, n

!

n = count (wall%pt%s >= s_min .and. wall%pt%s <= s_max) 
if (allocated(x)) deallocate(x, s)
allocate(s(n), x(n)) 

s = pack (wall%pt%s, wall%pt%s >= s_min .and. wall%pt%s <= s_max)
x = pack (wall%pt%x, wall%pt%s >= s_min .and. wall%pt%s <= s_max) * 1000

end subroutine get_data_points

end subroutine synrad_plot_sx

end module
