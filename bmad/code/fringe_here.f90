!+
! Function fringe_here (ele, orbit, particle_at) result (is_here)
!
! Routine to determine if there is a fringe field which must be tracked through.
!
! Input:
!   ele         -- ele_struct: Lattice element.
!   orbit       -- coord_struct: Particle position.
!   particle_at -- integer: Either first_track_edge$ or second_track_edge$.
!
! Output:
!   is_here     -- logical: True if there is a fringe. False if not.
!-

function fringe_here (ele, orbit, particle_at) result (is_here)

use bmad_interface, dummy => fringe_here

implicit none

type (ele_struct) ele
type (coord_struct) orbit

integer particle_at, fringe_at, physical_end
logical is_here

!

fringe_at = nint(ele%value(fringe_at$))
physical_end = physical_ele_end (particle_at, orbit, ele%orientation)
is_here = at_this_ele_end(physical_end, fringe_at)

if (nint(ele%value(fringe_type$)) == none$) is_here = .false.

end function fringe_here

