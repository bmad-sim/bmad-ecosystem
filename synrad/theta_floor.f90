!+
! Function theta_floor (s, lat, theta_base) result (theta_fl)
! 
! Returns the angle of the reference coordinate system in the
! global reference system.
!
! Input:
!   s   -- Real(rp): Longitudinal position.
!   lat -- Lat_struct: Lattice
!   theta_base -- Real(rp), optional: If present then make sure
!                   that theta_fl is within pi of theta_base.
!
! Output:
!   theta_fl -- Real(rp): Floor angle.
!-

function theta_floor (s, lat, theta_base) result (theta_fl)

use bmad_struct

implicit none

type (lat_struct) lat
real(rp) s, theta_fl
real(rp), optional :: theta_base
integer, save :: ix_ele = 1

!

if (s < 0 .or. s > lat%ele(lat%n_ele_track)%s) then
  print *, 'ERROR IN THETA_FLOOR: S VALUE OUT OF RANGE:', s
  call err_exit
end if

! find element containing s.

do
  if (lat%ele(ix_ele)%s >= s) exit
  ix_ele = ix_ele + 1
enddo

do
  if (lat%ele(ix_ele-1)%s <= s) exit
  ix_ele = ix_ele - 1
enddo

! Compute floor angle.
! Note theta decreases in going through a bend with positive g.

theta_fl = lat%ele(ix_ele)%floor%theta
if (lat%ele(ix_ele)%key == sbend$) then
  theta_fl = theta_fl - (s - lat%ele(ix_ele)%s) * lat%ele(ix_ele)%value(g$)
end if

! make sure angle is near theta_base

if (present(theta_base)) then
  theta_fl = theta_fl + twopi * nint((theta_base-theta_fl) / twopi)
end if

end function
