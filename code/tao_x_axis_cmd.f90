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

if (what /= "s" .and. what /= "index" .and. what /= "ele_index") then
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
integer iu, n
real(rp) minn, maxx

!

plot%x_axis_type = what
iu = s%global%u_view

if (what == 's') then
  minn = s%u(iu)%model%ele_(0)%s
  n = s%u(iu)%model%n_ele_use
  maxx = s%u(iu)%model%param%total_length
elseif (what == 'ele_index') then
  minn = 0
  maxx = s%u(s%global%u_view)%design%n_ele_use 
elseif (what == 'index') then
! if no curves to scale then can't scale to index
  if (.not. associated(plot%graph(1)%curve)) then
    plot%valid = .false.
    return
  endif
  iu = plot%graph(1)%curve(1)%ix_universe
  if (iu == 0) iu = s%global%u_view
  call tao_find_data (err, s%u(iu), &
               plot%graph(1)%curve(1)%data_type, d1_ptr = d1_ptr)
  if (err) return
  minn = lbound(d1_ptr%d, 1)
  maxx = ubound(d1_ptr%d, 1)
endif

!call qp_calc_and_set_axis ('X', minn, maxx, 8, 15, 'GENERAL')
call qp_calc_and_set_axis ('X', minn, maxx, &
          nint(0.7 * plot%x_divisions), nint(1.3 * plot%x_divisions), 'GENERAL')
call qp_get_axis ('X', a_min = plot%x%min, a_max = plot%x%max)


end subroutine

end subroutine 






