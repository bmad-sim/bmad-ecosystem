!+
! Subroutine wall_hit_handler_custom (orb, ele, s, t)
!
! Dummy routine.
! This routine is called by the Runge-Kutta integrator odeint_bmad when a particle hits a wall.
! This routine can be replaced by a custom routine to do custom calculations.
!
! Modules needed:
!   use bmad
!
! Input:
!   orb   -- coord_struct: coordinates of particle.
!   ele   -- ele_struct: Element holding the aperture
!   s     -- real(rp): Longitudinal position from start of element.
!   t     -- real(rp): Time of hit. May be relative or absolute depending upon 
!              what type of tracking is being done.
!
! Output:
!   Any argument may be modified...
!-

subroutine wall_hit_handler_custom (orb, ele, s, t)

use bmad_interface, dummy => wall_hit_handler_custom

implicit none

type (coord_struct) :: orb
type (ele_struct) :: ele

real(rp) :: s, t

character(*), parameter :: r_name = 'wall_hit_handler_custom'

!

end subroutine wall_hit_handler_custom
