!+
! Function entering_element(orbit, particle_at) result (is_entering)
!
! For a particle crossing an edge fringe, this routine indicats if the particle traveling towards the 
! inside of element or is exiting from inside to outside. 
!
! Entering and exiting is always with respect to a particle moving forward in time.
! That is, when tracking backwards in time (orbit%time_dir = -1), the particle at the first_track_edge$
! will be tagged as exiting the element and at the second_track_edge is entering the element.
!
! Input:
!   orbit       -- coord_struct: Particle orbit.
!   particle_at -- integer: First_track_edge$ or second_track_edge$
!
! Output:
!   is_entering -- logical: Set True if particle is going from outside to inside and vice versa.
!-

function entering_element(orbit, particle_at) result (is_entering)

use bmad_struct

implicit none

type (coord_struct) orbit
integer particle_at
logical is_entering

!

is_entering = (particle_at == first_track_edge$ .eqv. orbit%time_dir == 1)

end function
