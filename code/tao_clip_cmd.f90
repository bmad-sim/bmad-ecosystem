!+
! Subroutine tao_clip_cmd (where, y_min, y_max)
!
! Routine to veto (clip) data points whose plotted values lie outside of a 
! given range. 
! The veto is applied by setting %good_user = .false. for the data.
! If y_min = y_max then the clip range will be the y_min and y_max of 
! the graph.
! 
! Input:
!   where   -- Character(*): Graph() to clip. Eg: 'top:x'
!   y_min   -- Real(rp): Min clip value.
!   y_max   -- Real(rp): Max clip value.
!
!  Output:
!-

subroutine tao_clip_cmd (where, y_min, y_max)

use tao_mod

implicit none

type (tao_plot_struct), pointer :: p
type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)

real(rp) y_min, y_max

integer i, j, ix

character(*) where

logical err

! If the where argument is blank then clip all graphs

if (len_trim(where) == 0) then
  do j = 1, size(s%plot_region)
    p => s%plot_region(j)%plot
    if (.not. s%plot_region(j)%visible) cycle
    do i = 1, size(p%graph)
      call clip_graph (p, p%graph(i))
    enddo
  enddo
  return
endif

! locate the plot by the region name given by the where argument.

call tao_find_plots (err, where, 'REGION', plot, graph)
if (err) return

if (allocated(graph)) then       ! If all the graphs of a plot...
  do i = 1, size(graph)
    call clip_graph (graph(i)%g%p, graph(i)%g)
  enddo

else                    ! else just the one graph...
  do j = 1, size(plot)
    do i = 1, size(plot(j)%p%graph)
      call clip_graph (plot(j)%p, plot(j)%p%graph(i))
    enddo
  enddo
endif

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
contains

! Routine to clip one graph.

subroutine clip_graph (plot, graph)

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_d1_data_array_struct), allocatable, save :: d1_array(:)
type (tao_curve_struct), pointer :: curve

real(rp) this_min, this_max
integer i, j, k, iu

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
  if (.not. allocated(curve%y_symb)) cycle
  do j = 1, size(curve%y_symb)
    if (this_min <= curve%y_symb(j) .and. curve%y_symb(j) <= this_max) cycle
    call tao_find_data (err, curve%data_type, d1_array = d1_array, ix_uni = curve%ix_universe)
    if (err) return
    do k = 1, size(d1_array)
      d1_array(k)%d1%d(curve%ix_symb(j))%good_user = .false.  ! and clip it
    enddo
  enddo
enddo


end subroutine

end subroutine 




