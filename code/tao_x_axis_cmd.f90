!+
! Subroutine tao_x_axis_cmd (where, what)
!
! Routine to axis a plot. If x_min = x_max
! Then the axiss will be chosen to show all the data.
! 
! Input:
!   where -- Character(*): Region to axis. Eg: "top"
!   what  -- Character(*): "s" or "index"
!-

subroutine tao_x_axis_cmd (where, what)

use tao_mod
use quick_plot

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

integer j

character(*) where, what
character(16) :: r_name = 'tao_x_axis_cmd'

logical err

! check what

if (what /= "s" .and. what /= "index") then
  call out_io (s_error$, r_name, 'I DO NOT UNDERSTAND THIS: ' // what)
  return
endif

! If the where argument is blank then axis all graphs

if (len_trim(where) == 0 .or. where == 'all') then
  do j = 1, size(s%plot_page%plot)
    plot => s%plot_page%plot(j)
    if (.not. plot%visible) cycle
    call set_axis (plot)
  enddo
  return
endif

! locate the plot by the region name given by the where argument.

call tao_find_plot (err, s%plot_page%plot, 'BY_REGION', where, plot, graph)
if (err) return
call set_axis (plot)

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
contains

subroutine set_axis (plot)

type (tao_plot_struct) plot
type (tao_d1_data_struct), pointer :: d1_ptr
integer iu
real(rp) minn, maxx

!

iu = plot%graph(1)%curve(1)%ix_universe

if (what == 's') then
  minn = 0
  maxx = s%u(iu)%model%param%total_length
elseif (what == 'index') then
  call tao_find_data (err, s%u(iu), &
               plot%graph(1)%curve(1)%data_name, d1_ptr = d1_ptr)
  if (err) return
  minn = lbound(d1_ptr%d, 1)
  maxx = ubound(d1_ptr%d, 1)
endif

call qp_calc_and_set_axis ('X', minn, maxx, 8, 15, 'GENERAL')
call qp_get_axis ('X', a_min = plot%x%min, a_max = plot%x%max)
plot%x_axis_type = what


end subroutine

end subroutine 






