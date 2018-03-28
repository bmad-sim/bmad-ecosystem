!+
! Subroutine tao_hook_draw_floor_plan (plot, graph)
!
! Subroutine to customize the plotting of the floor_plan.
! Also see: tao_hook_draw_graph.
!
! Input:
!   plot  -- Tao_plot_struct: Plot structure containing the graph.
!   graph -- Tao_graph_struct: The graph to calculate the data for.
!-

subroutine tao_hook_draw_floor_plan (plot, graph)

use tao_interface, dummy => tao_hook_draw_floor_plan

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct) graph

! This dummy routine does nothing.

end subroutine
