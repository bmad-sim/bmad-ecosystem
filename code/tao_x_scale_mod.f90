module tao_x_scale_mod

use tao_mod
use quick_plot

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_cmd (where, x_min, x_max, err)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   where -- Character(*): Region to scale. Eg: "top"
!   x_min -- Real(rp): Plot x-axis min value.
!   x_max -- Real(rp): Plot x-axis max value.
!
!  Output:
!-

subroutine tao_x_scale_cmd (where, x_min, x_max, err)

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

real(rp) x_min, x_max

integer i, j, ix, places, im

character(*) where

real(rp) this_min, this_max

logical err

! If the where argument is blank then scale all graphs.
! Make sure that all plots using an 's' scale have the same scale

if (len_trim(where) == 0 .or. where == 'all') then

  im = 0
  do j = 1, size(s%plot_page%plot)
    plot => s%plot_page%plot(j)
    if (.not. plot%visible) cycle
    call tao_x_scale_plot (plot, x_min, x_max)
    if (plot%x_axis_type == 's') then
      if (im == 0) im = j
      if (plot%x%max < s%plot_page%plot(im)%x%max) im = j
    endif
  enddo

  if (im /= 0) then
    do j = 1, size(s%plot_page%plot)
      plot => s%plot_page%plot(j)
      if (.not. plot%visible) cycle
      if (plot%x_axis_type /= 's') cycle
      plot%x%max = s%plot_page%plot(im)%x%max
      plot%x%min = s%plot_page%plot(im)%x%min
      plot%x%major_div = s%plot_page%plot(im)%x%major_div
    enddo
  endif

  return
endif

! locate the plot by the region name given by the where argument.
! If where has a ':' then we are dealing with just one graph of the plot.
! Otherwise we scale all the graphs of the plot.

call tao_find_plot (err, s%plot_page%plot, 'BY_REGION', where, plot, graph)
if (err) return
call tao_x_scale_plot (plot, x_min, x_max)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_x_scale_plot (plot, x_min, x_max)

type (tao_plot_struct) plot
integer j, k, n
real(rp) x_min, x_max
real(rp) x1, x2

! Check if the thing exists

if (.not. associated (plot%graph)) return

! auto scale

if (x_max == x_min) then

  plot%x%min = -1e20
  plot%x%max = 1e20
  call tao_plot_data_setup 

  if (plot%x_axis_type == "index") then
    x1 =  1e20
    x2 = -1e20
    do j = 1, size(plot%graph)
      do k = 1, size(plot%graph(j)%curve)
        n = size(plot%graph(j)%curve(k)%x_symb)
        x1 = min (x1, plot%graph(j)%curve(k)%x_symb(1))
        x2 = max (x2, plot%graph(j)%curve(k)%x_symb(n))
      enddo
    enddo

  elseif (plot%x_axis_type == "s") then
    x1 = 0
    x2 = maxval (s%u(:)%model%param%total_length)
  endif

! not auto scale

else
  x1 = x_min
  x2 = x_max
endif

! calculate divisions and places

call qp_calc_and_set_axis ('X', x1, x2, &
          nint(0.7 * plot%x_divisions), nint(1.3 * plot%x_divisions), 'GENERAL')
call qp_get_axis ('X', plot%x%min, plot%x%max, plot%x%major_div, plot%x%places)

end subroutine 

end module
