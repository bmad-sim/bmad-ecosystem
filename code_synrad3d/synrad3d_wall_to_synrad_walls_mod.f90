module synrad3d_wall_to_synrad_walls_mod

use wall3d_mod
use synrad_mod

implicit none

private calc_this_side, add_this_point, calc_this_x

contains

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!+
! Subroutine synrad3d_wall_to_synrad_walls (branch, seg_len_max, walls)
!
! Routine to convert from a synrad3d wall structure to a synrad wall structure.
!
! Input:
!   branch       -- branch_struct: lattice branch with wall3d.
!   seg_len_max  -- Real(rp): Maximum length of wall segments.
!
! Output:
!   walls -- Walls_struct: synrad wall structure.
!-

subroutine synrad3d_wall_to_synrad_walls (branch, seg_len_max, walls)

type (branch_struct), target :: branch
type (walls_struct), target :: walls

real(rp) seg_len_max

!

call calc_this_side(+1, branch, seg_len_max, walls%positive_x_wall)
call calc_this_side(-1, branch, seg_len_max, walls%negative_x_wall)

call synrad_setup_walls (walls, branch, seg_len_max)

end subroutine synrad3d_wall_to_synrad_walls

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

subroutine calc_this_side (dir, branch, seg_len_max, wall)

use super_recipes_mod, only: super_zbrent

type (branch_struct), target :: branch
type (wall3d_struct), pointer :: wall3d
type (wall3d_section_struct), pointer :: sec
type (wall_struct) :: wall

real(rp) seg_len_max
real(rp) x, x2, x_this, s, len_sec, x_this0, x_this1
real(rp) :: dx_max = 1d-3, ds_small = 1d-9

integer dir
integer is, ip, iw, iw0, iw1, iw2, ix, ia, n
integer n_seg, iw_max, idel, ix_pt, status

logical no_wall_here, no_wall_here0, no_wall_here1

! Add section points but only if they are on the outside of the chamber

wall%n_pt_max = -1

do iw = 1, size(branch%wall3d)
  wall3d => branch%wall3d(iw)

  sec_loop: do is = 1, ubound(wall3d%section, 1)
    sec => branch%wall3d(iw)%section(is)

    x = calc_this_x (branch, iw, sec%s, dir, no_wall_here)
    if (no_wall_here) cycle

    ! Check if point is on the outside. 

    do iw2 = 1, size(branch%wall3d)
      if (iw2 == iw) cycle
      x_this = calc_this_x (branch, iw2, sec%s, dir, no_wall_here)
      if (no_wall_here) cycle
      if (x * dir < x_this * dir) cycle sec_loop
    enddo

    ! Now insert this point into wall%pt

    ix_pt = 0
    if (allocated(wall%pt)) then
      ix_pt = bracket_index (sec%s, wall%pt(0:wall%n_pt_max)%s, 0)
      ix_pt = ix_pt+1
    endif

    call add_this_point (wall, ix_pt, sec%s, x, sec%name, iw, sec%surface%name == 'PHANTOM')

    ! If this section is at the end of a sub-chamber then add another point.

    if (sec%s == 0 .or. sec%s == branch%ele(branch%n_ele_track)%s) cycle

    select case (sec%type)
    case (wall_start$); idel = -1
    case (wall_end$);   idel = +1; ix_pt = ix_pt + 1
    case default;       cycle
    end select

    s = sec%s + idel * ds_small
    x = 0
    do iw2 = 1, size(branch%wall3d)
      x_this = calc_this_x (branch, iw2, s, dir, no_wall_here)
      if (no_wall_here) cycle
      if (x * dir > x_this * dir) cycle
      x = x_this
      iw_max = iw2
    enddo

    call add_this_point (wall, ix_pt, s, x, '', iw_max, sec%surface%name == 'PHANTOM')

  enddo sec_loop
enddo

! Add an extra point in between section points where the outside of the chamber
! transitions from one sub-chamber to another.

ip = -1
do
  ip = ip + 1
  if (ip+1 > wall%n_pt_max) exit
  if (abs(wall%pt(ip+1)%s - wall%pt(ip)%s) < 1e3*ds_small) cycle

  ! %ix_seg is being used temperarily as sub-chamber index. See add_this_point.
  iw0 = wall%pt(ip)%ix_seg
  iw1 = wall%pt(ip+1)%ix_seg
  if (iw1 == iw0) cycle

  ! Test the center point to see if both walls extend into the region between the two sections.

  s = (wall%pt(ip)%s + wall%pt(ip+1)%s) / 2
  x_this0 = calc_this_x (branch, iw0, s, dir, no_wall_here0)
  x_this1 = calc_this_x (branch, iw1, s, dir, no_wall_here1)

  if (no_wall_here0 .and. no_wall_here1) then
    print *, 'Confused wall3d to wall conversion. Please report this!'
    stop
 
  else if (no_wall_here0) then
    s = wall%pt(ip)%s + 1d-10
    x = calc_this_x (branch, iw1, s, dir, no_wall_here)
    call add_this_point (wall, ip+1, s, x, '', iw1, wall%pt(ip+1)%phantom)

  else if (no_wall_here1) then
    s = wall%pt(ip+1)%s - 1d-10
    x = calc_this_x (branch, iw0, s, dir, no_wall_here)
    call add_this_point (wall, ip+1, s, x, '', iw0, wall%pt(ip)%phantom)

  ! If both walls extend into the region find where the walls intersect and add a point there.
  else
    s = super_zbrent(dwall_func, wall%pt(ip)%s, wall%pt(ip+1)%s, 0.0_rp, 1e-9_rp, status)
    x = calc_this_x (branch, iw0, s, dir, no_wall_here)
    call add_this_point (wall, ip+1, s, x, '', iw0, wall%pt(ip)%phantom)
  endif

  ip = ip + 1
enddo

!----------------------------------
contains

function dwall_func (s, status) result (value)

real(rp), intent(in) :: s
real(rp) value
integer status
logical no_wall_here

!

value = calc_this_x(branch, iw1, s, dir, no_wall_here) - calc_this_x(branch, iw0, s, dir, no_wall_here)
status = 0

end function dwall_func

end subroutine calc_this_side

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

subroutine add_this_point (wall, ix_pt, s, x, name, ix_wall, phantom)

type (wall_struct) wall
type (wall_pt_struct), allocatable :: pt(:)

real(rp) s, x
integer ix_pt, ix_wall, n
logical phantom
character(*) name

!

if (.not. allocated(wall%pt)) then
  allocate (wall%pt(0:10))
elseif (wall%n_pt_max == ubound(wall%pt, 1)) then
  n = wall%n_pt_max
  call move_alloc(wall%pt, pt)
  allocate(wall%pt(0:2*n))
  wall%pt(0:ix_pt-1) = pt(0:ix_pt-1)
  wall%pt(ix_pt+1:) = pt(ix_pt:)
else
  n = wall%n_pt_max
  wall%pt(ix_pt+1:n+1) = wall%pt(ix_pt:n)
endif


wall%pt(ix_pt)%s = s
wall%pt(ix_pt)%x = x
wall%pt(ix_pt)%name = name
wall%pt(ix_pt)%ix_seg = ix_wall  ! Temp storage of sub-chamber index
wall%pt(ix_pt)%phantom = phantom
wall%n_pt_max = wall%n_pt_max + 1

end subroutine add_this_point

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

function calc_this_x (branch, ix_wall, s, dir, no_wall_here) result (x)

use super_recipes_mod

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele            ! For super_zbrent

real(rp) s, x, dw_perp(3), origin(3)
real(rp) x2, y0, y1, y2
real(rp) position(6), x_wall                 ! For super_zbrent
integer dir, ixs, ix_ele, ix_wall, status
logical no_wall_here, err_flag

!

ix_ele = element_at_s (branch, s, .true.)
ele => branch%ele(ix_ele)
x = 0

position = [dir * 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, s - (ele%s - ele%value(l$)), 0.0_rp]
x2 = wall3d_d_radius (position, ele, ix_wall, dw_perp, ixs, no_wall_here, origin)
if (no_wall_here) return

y1 = this_y(10.0_rp, status)
y2 = this_y(-10.0_rp, status)

if (y1*y2 >= 0) then ! Does not straddle zero
  no_wall_here = .true.
  return
endif

y0 = super_zbrent (this_y, -10.0_rp, 10.0_rp, 0.0_rp, 1d-8, status)
if (status /= 0) call err_exit
y1 = this_y(y0, status)

x = x_wall

!--------------------------------------------------------------------------------
contains

function this_y (y_in, status) result (y_wall)

real(rp), intent(in) :: y_in
real(rp) y_wall, dw_perp(3), origin(3), dr, r_part, r_wall
integer status, ixs
logical no_wall_here

!

position(3) = y_in
dr = wall3d_d_radius (position, ele, ix_wall, dw_perp, ixs, no_wall_here, origin)
r_part = sqrt((position(1) - origin(1))**2 + (position(3) - origin(2))**2)
r_wall = r_part - dr

x_wall = origin(1) + (position(1) - origin(1)) * r_wall / r_part
y_wall = origin(2) + (position(3) - origin(2)) * r_wall / r_part

end function this_y

end function calc_this_x

end module
