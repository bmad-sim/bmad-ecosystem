!+
! Subroutine tao_hook_graph_postsetup (plot, graph)
!
! Subroutine to customize the calculation of data to be plotted.
!
! Input:
!   plot  -- Tao_plot_struct: Plot structure containing the graph.
!   graph -- Tao_graph_struct: The graph to calculate the data for.
!
! Output:
!   graph -- Tao_graph_struct: Structure with plot data calulated.
!-

subroutine tao_hook_graph_postsetup (plot, graph)

use tao_interface, dummy => tao_hook_graph_postsetup

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct) graph

! 

end subroutine
