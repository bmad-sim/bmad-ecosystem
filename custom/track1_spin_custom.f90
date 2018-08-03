!+
! Subroutine track1_spin_custom (start, ele, param, end, err_flag, make_quaternion)
!
! Dummy routine for custom spin tracking. 
! This routine needs to be replaced for a custom calculation.
!
! This routine is called if ele%spin_tracking_method = custom$
!
! Note: If ele%spin_tracking_method = tracking$, this routine is not called and
! spin tracking must be done with the phase space tracking in track1_custom.
! 
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below."
!
! Input:
!   start            -- Coord_struct: Starting position.
!   ele              -- Ele_struct: Element.
!   param            -- lat_param_struct: Lattice parameters.
!   make_quaternion  -- logical, optional: If present and true then calculate the 1st
!                          order spin map which is represented as a quaternion.
!
! Output:
!   end   -- Coord_struct: End position.
!-

subroutine track1_spin_custom (start, ele, param, end, err_flag, make_quaternion)

use bmad_interface, except_dummy => track1_spin_custom

implicit none

type (coord_struct) :: start
type (coord_struct) :: end
type (ele_struct) :: ele
type (lat_param_struct) :: param

logical err_flag
logical, optional :: make_quaternion
character(*), parameter :: r_name = 'track1_spin_custom'

!

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
err_flag = .true.

end subroutine
