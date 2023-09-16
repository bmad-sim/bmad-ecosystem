!+
! Function theta_floor (s, branch, theta_base) result (theta_fl)
! 
! Returns the angle of the reference coordinate system in the
! global reference system.
!
! Input:
!   s   -- Real(rp): Longitudinal position.
!   branch -- Branch_struct: Lattice
!   theta_base -- Real(rp), optional: If present then make sure
!                   that theta_fl is within pi of theta_base.
!
! Output:
!   theta_fl -- Real(rp): Floor angle.
!-

function theta_floor (s, branch, theta_base) result (theta_fl)

use bmad_interface

implicit none

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele

real(rp) s, theta_fl
real(rp), optional :: theta_base
integer, save :: ix_ele = 1

!

if (s < 0 .or. s > branch%ele(branch%n_ele_track)%s) then
  print *, 'ERROR IN THETA_FLOOR: S VALUE OUT OF RANGE:', s
  if (global_com%exit_on_error) call err_exit
end if

! find element containing s.

do
  if (branch%ele(ix_ele)%s >= s) exit
  ix_ele = ix_ele + 1
enddo

do
  if (branch%ele(ix_ele-1)%s <= s) exit
  ix_ele = ix_ele - 1
enddo

ele => branch%ele(ix_ele)

! Compute floor angle.
! Note theta decreases in going through a bend with positive g.

theta_fl = ele%floor%theta
if (ele%key == sbend$) then
  theta_fl = theta_fl - (s - ele%s) * ele%value(g$) * cos(ele%value(ref_tilt_tot$))
end if

! make sure angle is near theta_base

if (present(theta_base)) then
  theta_fl = theta_fl + twopi * nint((theta_base-theta_fl) / twopi)
end if

end function
