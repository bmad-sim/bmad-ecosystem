!+
! Subroutine tao_plot_cmd (where, component)
!
! Routine to set what is plotted Eg: model - design, etc.
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
    do j = 1, size(s%plot_page%region(i)%plot%graph)
      s%plot_page%region(i)%plot%graph(j)%component = component
    enddo
  enddo

else
  call tao_find_plots (err, where, 'REGION', plot, graph)
  if (err) return
  if (size(graph) > 0) then
    do i = 1, size(graph)
      graph(i)%g%component = component
    enddo
  else
    do i = 1, size(plot)
      do j = 1, size(plot(i)%p%graph)
        plot(i)%p%graph(j)%component = component
      enddo
    enddo
  endif
endif

err = .false.

end subroutine 




