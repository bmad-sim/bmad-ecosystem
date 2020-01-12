! Function tao_graph_name (graph, use_region) result (graph_name)
!
! Function to return the graph name in the form:
!   plot_name.graph_name
! For example:
!   orbit.x
!
! Input:
!   graph      -- Tao_graph_struct: Graph
!   use_region -- Logical: If present and True then use the region 
!                  name instead of the plot name. Region name is
!                  'NULL_REGION' if there is no assocaited region.
!
! Output:
!   graph_name -- Character(60): Appropriate name.
!-

function tao_graph_name(graph, use_region) result (graph_name)

use tao_struct

implicit none

type (tao_graph_struct) graph
character(60) graph_name
logical, optional :: use_region

!

graph_name = '.' // trim(graph%name)

if (logic_option(.false., use_region)) then
  if (associated(graph%p%r)) then
    graph_name = trim(graph%p%r%name) // trim(graph_name)
  else
    graph_name = 'NULL_REGION' // trim(graph_name)
  endif
else
    graph_name = trim(graph%p%name) // trim(graph_name)
endif

end function tao_graph_name 
