!+
! subroutine track1_spin_sprint (start_orb, ele, param, end_orb)
!
! Particle spin tracking through a single element with the sprint quaternion map.
!
! Input :
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Beam parameters.
!   end_orb    -- Coord_struct: Ending coords.
!
! Output:
!   end_orb     -- Coord_struct:
!     %spin(3)   -- Ending spin
!-

subroutine track1_spin_sprint (start_orb, ele, param, end_orb)

use taylor_mod, dummy => track1_spin_sprint

implicit none

type (coord_struct) :: start_orb, end_orb
type (ele_struct) ele
type (lat_param_struct) param

character(*), parameter :: r_name = 'track1_spin_sprint'

!

if (ele%spin_q(0,0) == real_garbage$) then
endif

end_orb%spin = 0

end subroutine track1_spin_sprint


