!+
! Subroutine tao_clip_cmd (s, where, y_min, y_max)
!
! Routine to veto (clip) data points whose plotted values lie outside of a 
! given range. 
! The veto is applied by setting %good_user = .false. for the data.
! If y_min = y_max then the clip range will be the y_min and y_max of 
! the graph.
! 
! Input:
!   s       -- tao_super_universe_struct
!   where   -- Character(*): Graph(s) to clip. Eg: 'top:x'
!   y_min   -- Real(rp): Min clip value.
!   y_max   -- Real(rp): Max clip value.
!
!  Output:
!   s        -- tao_super_universe_struct
!-

subroutine tao_clip_cmd (s, where, y_min, y_max)

use tao_mod

implicit none

type (tao_super_universe_struct) s
type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

real(rp) y_min, y_max

integer i, j, ix

character(*) where

logical err

! If the where argument is blank then clip all graphs

if (len_trim(where) == 0) then
  do j = 1, size(s%plot_page%plot)
    plot => s%plot_page%plot(j)
    if (.not. plot%visible) cycle
    do i = 1, size(plot%graph)
      call clip_graph (plot, plot%graph(i))
    enddo
  enddo
  return
endif

! locate the plot by the region name given by the where argument.
! If where has a ':' then we are dealing with just one graph of the plot.
! Otherwise we clip all the graphs of the plot.

call tao_find_plot (err, s%plot_page%plot, 'BY_REGION', where, plot, graph)
if (err) return

ix = index(where, ':')
if (ix == 0) then       ! If all the graphs of a plot...
  do i = 1, size(plot%graph)
    call clip_graph (plot, plot%graph(i))
  enddo

else                    ! else just the one graph...
  call clip_graph (plot, graph)
endif

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
contains

! Routine to clip one graph.

subroutine clip_graph (plot, graph)

type (tao_plot_struct) plot
type (tao_graph_struct) graph
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_curve_struct), pointer :: curve

real(rp) this_min, this_max
integer i, j, iu
character(40) class

! If y_min = y_max then clip to graph boundries.

if (y_min == y_max) then
  this_min = graph%y%min
  this_max = graph%y%max
else
  this_min = y_min
  this_max = y_max
endif

! Loop over all data points of all the curves of the graph.
! If curve%ix_universe = 0 then universe the data came from is the
! universe currently on view (s%global%u_view).

do i = 1, size(graph%curve)
  curve => graph%curve(i)
  if (.not. associated(curve%y_dat)) cycle
  do j = 1, size(curve%y_dat)
    if (this_min <= curve%y_dat(j) .and. curve%y_dat(j) <= this_max) cycle
    iu = curve%ix_universe 
    if (iu == 0) iu = s%global%u_view
    class = trim(plot%class) // ':' // curve%sub_class
    call tao_find_data (err, s%u(iu), class, d1_ptr = d1_ptr)
    if (err) return
    d1_ptr%d(curve%ix_data(j))%good_user = .false.  ! and clip it
  enddo
enddo


end subroutine

end subroutine 




