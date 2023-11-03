!+
! Subroutine wall_hit_handler_custom (orb, ele, s)
!
! Prototype routine for customizing the action when a particle hits 
! a wall with Runge-Kutta integration in odeint_bmad.
!
! To use, see the Bmad manual.
!
! Input:
!   orb   -- coord_struct: coordinates of particle.
!   ele   -- ele_struct: Element holding the aperture
!   s     -- real(rp): Longitudinal position from start of element.
!
! Output:
!   orb   -- coord_struct: coordinates of particle.
!-

subroutine wall_hit_handler_custom (orb, ele, s)

use bmad_interface

implicit none

type (coord_struct) :: orb
type (ele_struct) :: ele

real(rp) :: s

character(*), parameter :: r_name = 'wall_hit_handler_custom'

!

end subroutine wall_hit_handler_custom
