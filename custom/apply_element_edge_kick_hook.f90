!+
! Subroutine apply_element_edge_kick_hook (orb, s_edge, t_rel, hard_ele, track_ele, param, particle_at, finished)
!
! Routine that can be customized to track through the edge field of an element.
! This routine is always called by apply_element_edge_kick.
! 
!
! Input:
!   orb         -- Coord_struct: Starting coords in element reference frame.
!   s_edge      -- real(rp): Hard edge relative to start of hard_ele.
!   t_rel       -- real(rp): Time relative to track_ele entrance edge
!   hard_ele    -- ele_struct: Element with hard edges.
!   track_ele   -- ele_struct: Element being tracked through. 
!                    Is different from hard_ele when there are superpositions.
!   param       -- lat_param_struct: lattice parameters.
!   particle_at -- Integer: first_track_edge$ or second_track_edge$
!
! Output:
!   orb         -- Coord_struct: Coords after edge kick applied.
!   finished    -- logical: When set True, apply_element_edge_kick will not apply any fringe effects.
!-

subroutine apply_element_edge_kick_hook (orb, s_edge, t_rel, hard_ele, track_ele, param, particle_at, finished)

use track1_mod

implicit none

type (ele_struct) hard_ele, track_ele
type (coord_struct) orb
type (lat_param_struct) param

real(rp) t_rel, s_edge

integer particle_at, physical_end
logical finished

!

physical_end = physical_ele_end (particle_at, orb%direction, track_ele%orientation)

finished = .false.   ! Must set this

end subroutine apply_element_edge_kick_hook
