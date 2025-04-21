!+
! subroutine track1_spin_magnus (start_orb, ele, param, end_orb)
!
! Particle spin tracking through a single element.
!
! Input :
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param            -- lat_param_struct: Beam parameters.
!
! Output:
!   end_orb     -- Coord_struct:
!     %spin(3)   -- Ending spin
!-

subroutine track1_spin_magnus (start_orb, ele, param, end_orb)

use taylor_mod, dummy => track1_spin_magnus

implicit none

type (coord_struct) :: start_orb, end_orb, orbit
type (ele_struct) ele
type (lat_param_struct) param

character(*), parameter :: r_name = 'track1_spin_magnus'

!

if (.not. ele%is_on) then
  end_orb%spin = start_orb%spin
  return
endif

!

if (start_orb%time_dir == 1) then
  orbit = start_orb
  if (.not. ele%taylor_map_includes_offsets) call offset_particle (ele, set$, orbit)
else
  orbit = end_orb
  if (.not. ele%taylor_map_includes_offsets) call offset_particle (ele, set$, orbit)
endif

!


!

if (.not. ele%taylor_map_includes_offsets) call offset_particle (ele, unset$, orbit)

end_orb%spin = start_orb%spin   ! This is for debugging before the real code is put in.

end subroutine track1_spin_magnus


