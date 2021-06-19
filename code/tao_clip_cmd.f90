!+
! Subroutine tao_clip_cmd (gang, where, y_min, y_max)
!
! Routine to veto (clip) data points whose plotted values lie outside of a 
! given range. 
!
! The veto is applied by setting %good_user = .false. for the data.
! If y_min = y_max then the clip range will be the y_min and y_max of 
! the graph.
!
! The gang argument, if True, transferrs 
! 
! Input:
!   gang    -- Logical: Gang all data d1 arrays together.
!   where   -- Character(*): Graph() to clip. Eg: 'top:x'
!   y_min   -- Real(rp): Min clip value.
!   y_max   -- Real(rp): Max clip value.
!
!  Output:
!-

subroutine tao_clip_cmd (gang, where, y_min, y_max)

use tao_interface, dummy => tao_clip_cmd

implicit none

type (tao_plot_struct), pointer :: p
type (tao_plot_array_struct), allocatable :: plot(:)
type (tao_graph_array_struct), allocatable :: graph(:)
type (tao_d2_data_struct), pointer :: d2_old

real(rp) y_min, y_max

integer i, j, ix

character(*) where

logical err, gang

! If the where argument is blank then clip all graphs

nullify(d2_old)

if (len_trim(where) == 0) then
  do j = 1, size(s%plot_page%region)
    p => s%plot_page%region(j)%plot
    if (.not. s%plot_page%region(j)%visible) cycle
    do i = 1, size(p%graph)
      call clip_graph (p, p%graph(i))
    enddo
  enddo
  return
endif

! locate the plot by the region name given by the where argument.

call tao_find_plots (err, where, 'REGION', plot, graph)
if (err) return

if (size(graph) > 0) then       ! If all the graphs of a plot...
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
type (tao_d1_data_array_struct), allocatable, target :: d1_array(:)
type (tao_d2_data_struct), pointer :: d2
type (tao_curve_struct), pointer :: curve

real(rp) this_min, this_max
integer i, j, k, iu, id

! If y_min = y_max then clip to graph boundries.

if (y_min == y_max) then
  this_min = graph%y%min
  this_max = graph%y%max
else
  this_min = y_min
  this_max = y_max
endif

! Loop over all data points of all the curves of the graph.

do i = 1, size(graph%curve)

  curve => graph%curve(i)
  if (.not. allocated(curve%y_symb)) cycle
  call tao_find_data (err, curve%data_type, d1_array = d1_array, ix_uni = tao_curve_ix_uni(curve))
  if (err) return    
  if (size(d1_array) == 0) cycle
  d2 => d1_array(1)%d1%d2

  do j = 1, size(curve%y_symb)
    if (this_min <= curve%y_symb(j) .and. curve%y_symb(j) <= this_max) cycle

    id = curve%ix_symb(j)
    d1_array(1)%d1%d(id)%good_user = .false.  ! and clip it

    if (gang) then
      do k = 1, size(d2%d1)
        d2%d1(k)%d(id)%good_user = .false.
      enddo
    endif

  enddo

  do k = 1, size(d2%d1)
    call tao_set_data_useit_opt (d2%d1(k)%d)
  enddo

  if (.not. associated (d2, d2_old)) call tao_data_show_use (d2)
  d2_old => d2

enddo


end subroutine

end subroutine 




