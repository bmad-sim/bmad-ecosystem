!+
! Subroutine check_aperture_limit_custom (orb, ele, particle_at, param, err_flag)
!
! Dummy routine.
! A valid check_aperture_limit_custom is needed only if ele%aperture_type is set to
! custom$.
!
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below."
!
! Modules needed:
!   use bmad
!
! Input:
!   orb   -- Coord_struct: coordinates of a particle.
!   ele   -- Ele_struct: Element holding the aperture
!   particle_at    -- Integer: first_track_edge$, second_track_edge$, or surface$
!   param -- lat_param_struct: Parameter structure
!
! Output:
!   orb       -- coord_struct: 
!   err_flag  -- Logical: Set true if there is an error. False otherwise.
!-

subroutine check_aperture_limit_custom (orb, ele, particle_at, param, err_flag)

use bmad_interface, dummy => check_aperture_limit_custom

implicit none

type (coord_struct) :: orb
type (ele_struct) :: ele
type (lat_param_struct) :: param

integer particle_at
logical err_flag

character(32) :: r_name = 'check_aperture_limit_custom'

!

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
err_flag = .true.

end subroutine
