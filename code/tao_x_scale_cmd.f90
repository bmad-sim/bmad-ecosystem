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

use tao_mod
use quick_plot

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

real(rp) x_min, x_max

integer i, j, ix, places

character(*) where

real(rp) this_min, this_max

logical err

! If the where argument is blank then scale all graphs

if (len_trim(where) == 0 .or. where == 'all') then
  do j = 1, size(s%plot_page%plot)
    plot => s%plot_page%plot(j)
    if (.not. plot%visible) cycle
    call scale_plot (plot)
  enddo
  return
endif

! locate the plot by the region name given by the where argument.
! If where has a ':' then we are dealing with just one graph of the plot.
! Otherwise we scale all the graphs of the plot.

call tao_find_plot (err, s%plot_page%plot, 'BY_REGION', where, plot, graph)
if (err) return
call scale_plot (plot)

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
contains

subroutine scale_plot (plot)

type (tao_plot_struct) plot
integer j, k, n

! auto scale

if (x_max == x_min) then
  if (plot%x_axis_type == "index") then
    plot%x%min =  1e20
    plot%x%max = -1e20
    do j = 1, size(plot%graph)
      do k = 1, size(plot%graph(j)%curve)
        n = size(plot%graph(j)%curve(k)%x_symb)
        plot%x%min = min (plot%x%min, plot%graph(j)%curve(k)%x_symb(1))
        plot%x%max = max (plot%x%max, plot%graph(j)%curve(k)%x_symb(n))
      enddo
    enddo
  elseif (plot%type == "s") then
    plot%x%min = 0
    plot%x%max = maxval (s%u(:)%model%param%total_length)
  endif

!

else
  plot%x%min = x_min
  plot%x%max = x_max
endif


end subroutine

end subroutine 






