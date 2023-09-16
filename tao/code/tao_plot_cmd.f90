!+
! Subroutine tao_plot_cmd (where, component)
!
! Routine to set what is plotted Eg: model - design, etc.
! NOTE: THIS COMMAND IS DEPRECATED 8/2021.
!
! Input:
!   where     -- Character(*): Region name to identify the plot to set.
!   component -- character(*): Who to plot. EG: 'meas - design'
!-

subroutine tao_plot_cmd (where, component)

use tao_interface, dummy => tao_plot_cmd

implicit none

type (tao_plot_array_struct), allocatable :: plot(:)
type (tao_graph_array_struct), allocatable :: graph(:)
type (tao_plot_region_struct), pointer :: region
type (tao_plot_struct), pointer :: p
integer i, j
integer ix, ix_line, ix_cmd, which, n_word

character(*) :: where
character(*) :: component
character(20) :: r_name = 'tao_plot_cmd'

logical err

! Find plot for the region given by "where"

err = .true.

if (where == '*' .or. where == 'all') then
  do i = 1, size(s%plot_page%region)
    p => s%plot_page%region(i)%plot
    do j = 1, size(p%graph)
      if (.not. allocated(p%graph(j)%curve)) cycle
      p%graph(j)%curve%component = component
    enddo
  enddo

else
  call tao_find_plots (err, where, 'REGION', plot, graph)
  if (err) return
  if (size(graph) > 0) then
    do i = 1, size(graph)
      if (.not. allocated(graph(i)%g%curve)) cycle
      graph(i)%g%curve%component = component
    enddo
  else
    do i = 1, size(plot)
      p => plot(i)%p
      do j = 1, size(p%graph)
        if (.not. allocated(p%graph(j)%curve)) cycle
        p%graph(j)%curve%component = component
      enddo
    enddo
  endif
endif

err = .false.

end subroutine 




