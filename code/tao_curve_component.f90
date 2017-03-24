!+
! Function tao_curve_component (curve, graph) result (component)
!
! Routine to return the component attribute of a curve.
! graph%component is the default if curve%component is not specified.
!
! Input:
!   curve   -- tao_curve_curve: Curve under consideration.
!   graph   -- tao_graph_struct: Parent graph.
!
! Output:
!   component -- character(60): Component attribute of curve.
!-

function tao_curve_component (curve, graph) result (component)

use tao_struct

implicit none

type (tao_curve_struct) curve
type (tao_graph_struct) graph

character(60) component

! graph%component is the default if curve%component is not specified.

if (curve%component == '') then
  component = graph%component
else
  component = curve%component
endif

end function tao_curve_component
