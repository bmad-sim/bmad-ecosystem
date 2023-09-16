!+
! subroutine track1_spin_taylor (start_orb, ele, param, end_orb)
!
! Particle spin tracking through a single element with a spin map.
!
! Input :
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Beam parameters.
!
! Output:
!   end_orb     -- Coord_struct:
!     %spin(3)   -- Ending spin
!-

subroutine track1_spin_taylor (start_orb, ele, param, end_orb)

use taylor_mod, dummy => track1_spin_taylor

implicit none

type (coord_struct) :: start_orb, end_orb, orbit
type (ele_struct) ele
type (lat_param_struct) param

real(rp) quat(4), norm
character(*), parameter :: r_name = 'track1_spin_taylor'

!

if (.not. associated(ele%spin_taylor(0)%term)) then
  if (ele%spin_tracking_method == sprint$) then
    call sprint_spin_taylor_map(ele)
  else
    call ele_to_taylor(ele, param)
  endif
endif

!

if (start_orb%time_dir == 1) then
  orbit = start_orb
  if (.not. ele%taylor_map_includes_offsets) call offset_particle (ele, set$, orbit)
  quat = track_taylor (orbit%vec, ele%spin_taylor, ele%taylor%ref)
else
  orbit = end_orb
  if (.not. ele%taylor_map_includes_offsets) call offset_particle (ele, set$, orbit)
  quat = track_taylor (orbit%vec, ele%spin_taylor, ele%taylor%ref)
  quat = quat_inverse(quat)
endif

!

norm = norm2(quat)
if (abs(norm - 1) > 0.5) then
  call out_io (s_warn$, r_name, 'Norm of quaternion computed from the spin taylor map of element: ' // ele%name, &
                                'is far from 1.0: \es10.2\ ', r_array = [norm])
endif

orbit%spin = quat_rotate(quat/norm, orbit%spin)

!

if (.not. ele%taylor_map_includes_offsets) call offset_particle (ele, unset$, orbit)
end_orb%spin = orbit%spin

end subroutine track1_spin_taylor


