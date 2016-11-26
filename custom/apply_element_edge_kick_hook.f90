!+
! Subroutine apply_element_edge_kick_hook (orb, s_edge, hard_ele, track_ele, param, finished, mat6, make_matrix)
!
! Routine that can be customized to track through the edge field of an element.
! This routine is always called by apply_element_edge_kick.
! 
! Input:
!   orb         -- Coord_struct: Starting coords in element reference frame.
!   track_ele   -- ele_struct: Element being tracked through. 
!                    Is different from hard_ele when there are superpositions.
!   param       -- lat_param_struct: lattice parameters.
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before fringe.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orb         -- Coord_struct: Coords after edge kick applied.
!   finished    -- logical: When set True, apply_element_edge_kick will not apply any fringe effects.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix transfer matrix including fringe.
!-

subroutine apply_element_edge_kick_hook (orb, fringe_info, track_ele, param, finished, mat6, make_matrix)

use track1_mod, dummy => apply_element_edge_kick_hook

implicit none

type (ele_struct) track_ele
type (fringe_edge_info_struct) fringe_info
type (coord_struct) orb
type (lat_param_struct) param

integer physical_end

real(rp), optional :: mat6(6,6)

logical, optional :: make_matrix
logical finished

!

physical_end = physical_ele_end (fringe_info%particle_at, orb%direction, track_ele%orientation)

finished = .false.   ! Must set this

end subroutine apply_element_edge_kick_hook
