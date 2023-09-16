!+
! Subroutine map_to_angle_coords (t_canon, t_angle)
!
! Routine to convert a Taylor map using from Bmad phase space canonical coordinates to
! a map using angle coordinates (x, x', y, y', z, pz).
!
! Note: The actual t_angle argument may be the same as the actual t_canon argument.
!
! Input:
!   t_canon(6)  -- taylor_struct: Taylor map in canonical coords.
!
! Output:
!   t_angle(6)  -- taylor_struct: Taylor map in angle coords.
!-

subroutine map_to_angle_coords (t_canon, t_angle)

use ptc_interface_mod, dummy => map_to_angle_coords
use s_fibre_bundle

implicit none

type (real_8) r1(6), r2(6), r3(6)
type (taylor_struct) t_canon(6), t_angle(6)
type (damap) id, da1, da2, da3

real(rp) v0(6), f

!

call alloc(id, da1, da2, da3)

v0 = t_canon%ref
f = 1 / sqrt((1 + v0(6))**2 - v0(2)**2 - v0(4)**2)
v0(2) = v0(2) * f
v0(4) = v0(4) * f

id = 1        ! Ident map
da2 = id + v0  ! Ident map + const

! da1 is transformation from phase to angle around v0

da1 = da2
da1%v(2) = da2%v(2) * (1 + da2%v(6)) / sqrt(1 + da2%v(2)**2 + da2%v(4)**2)
da1%v(4) = da2%v(4) * (1 + da2%v(6)) / sqrt(1 + da2%v(2)**2 + da2%v(4)**2)

da2 = t_canon
da3 = da2 .o. da1

da1 = da3
da1%v(2) = da3%v(2) / sqrt((1 + da2%v(6))**2 - da2%v(2)**2 - da2%v(4)**2)
da1%v(4) = da3%v(4) / sqrt((1 + da2%v(6))**2 - da2%v(2)**2 - da2%v(4)**2)

t_angle = da1
t_angle%ref = v0

call kill(da1)
call kill(da2)
call kill(da3)
call kill(id)

end subroutine map_to_angle_coords
