!+
! Subroutine check_aperture_limit_custom (orb, ele, at, param, err_flag)
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
!     %value(x1_limit$) -- Horizontal negative side aperture.
!     %value(x2_limit$) -- Horizontal positive side aparture.
!     %value(y1_limit$) -- Vertical negative side aperture.
!     %value(y2_limit$) -- Vertical positive side aparture.
!     %offset_moves_aperture -- If True then aperture moves with the element.
!   at    -- Integer: entrance_end$ or exit_end$
!   param -- lat_param_struct: Parameter structure
!     %aperture_limit_on -- The aperture limit is only checked if this is true.
!               The exception is when the orbit is larger than 
!               bmad_com%max_aperture_limit.
!
! Output:
!   orb       -- coord_struct: 
!   err_flag  -- Logical: Set true if there is an error. False otherwise.
!-

subroutine check_aperture_limit_custom (orb, ele, at, param, err_flag)

use bmad_struct
use bmad_interface

implicit none

type (coord_struct) :: orb
type (ele_struct) :: ele
type (lat_param_struct) :: param

integer at
logical err_flag

character(32) :: r_name = 'check_aperture_limit_custom'

!

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
err_flag = .true.

end subroutine
