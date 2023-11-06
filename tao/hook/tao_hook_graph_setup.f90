!+
! Subroutine tao_hook_graph_setup (plot, graph, found)
!
! Subroutine to customize the calculation of data to be plotted.
!
! Input:
!   plot  -- Tao_plot_struct: Plot structure containing the graph.
!   graph -- Tao_graph_struct: The graph to calculate the data for.
!
! Output:
!   graph -- Tao_graph_struct: Structure with plot data calulated.
!   found -- Logical: Set true if this routine calculates data.
!               Set false if data calc is to be left to tao_plot_setup.
!-

subroutine tao_hook_graph_setup (plot, graph, found)

use tao_interface

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct) graph
logical found

! This dummy routine just sets found and returns.

found = .false.

end subroutine
