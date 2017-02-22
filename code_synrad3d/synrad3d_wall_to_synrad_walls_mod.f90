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

type (branch_struct), target :: branch
type (wall3d_struct), pointer :: wall3d
type (wall3d_section_struct), pointer :: sec
type (wall_struct) :: wall

real(rp) seg_len_max
real(rp) x, x2, x_this, s, len_sec
real(rp) :: dx_max = 1d-3, ds_small = 1d-9

integer dir
integer is, ip, iw, iw0, iw2, ix, ia, n
integer n_seg, iw_max, idel, ix_pt

logical no_wall_here

! Add section points but only if they are on the outside of the chamber

do iw = 1, size(branch%wall3d)
  wall3d => branch%wall3d(iw)

  sec_loop: do is = 1, ubound(branch%wall3d(iw)%section, 1)
    sec => branch%wall3d(iw)%section(is)

    x = calc_this_x (branch, iw, sec%s, dir, no_wall_here)
    if (no_wall_here) cycle

    ! Check if point is on the outside. 

    do iw2 = 1, size(branch%wall3d)
      x_this = calc_this_x (branch, iw2, sec%s, dir, no_wall_here)
      if (no_wall_here) cycle
      if (x * dir < x_this * dir) cycle sec_loop
    enddo

    ! Now insert this point into wall%pt

    ix_pt = 0
    if (allocated(wall%pt)) then
      call bracket_index (wall%pt%s, 0, ubound(wall%pt, 1), sec%s, ix_pt)
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

! Add extra points if needed in between section points where the outside of the chamber
! transitions from one sub-chamber to another.

ip = -1
do
  ip = ip + 1
  if (ip+1 > ubound(wall%pt, 1)) exit
  if (abs(wall%pt(ip+1)%s - wall%pt(ip)%s) < 1e3*ds_small) cycle
  ! %ix_seg is being used temperarily as sub-chamber index
  iw0 = wall%pt(ip)%ix_seg

  len_sec = sqrt((wall%pt(ip)%s - wall%pt(ip+1)%s)**2 + (wall%pt(ip)%x - wall%pt(ip+1)%x)**2)
  n_seg = 1 + 10 * len_sec / seg_len_max 
  do ia = 0, n_seg
    s = wall%pt(ip)%s + ia * (wall%pt(ip+1)%s - wall%pt(ip)%s) / n_seg
    if (ia == 0) s = s + 1d-10
    if (ia == n_seg) s = s - 1d-10

    x = 0
    iw_max = -1
    do iw2 = 1, size(branch%wall3d)
      x_this = calc_this_x (branch, iw2, s, dir, no_wall_here)
      if (no_wall_here) cycle
      if (x_this * dir < x * dir) cycle
      x = x_this
      iw_max = iw2
    enddo

    if (iw_max /= iw0 .and. iw_max /= -1) then
      call add_this_point (wall, ip+1, s, x, '', iw_max, wall%pt(ip+1)%phantom)
      exit
    endif
  enddo
enddo

! Add extra points if needed in between in case the outside of the chamber is not
! a straight line

ip = -1
do
  ip = ip + 1
  if (ip+1 > ubound(wall%pt, 1)) exit
  if (abs(wall%pt(ip+1)%s - wall%pt(ip)%s) < 1d3*ds_small) cycle
  len_sec = sqrt((wall%pt(ip+1)%s - wall%pt(ip)%s)**2 + (wall%pt(ip+1)%x - wall%pt(ip)%x)**2)
  n_seg = 1 + len_sec / seg_len_max
  do ia = 1, n_seg - 1
    s = wall%pt(ip)%s + ia * (wall%pt(ip+1)%s - wall%pt(ip)%s) / n_seg

    x = 0
    do iw2 = 1, size(branch%wall3d)
      x_this = calc_this_x (branch, iw2, s, dir, no_wall_here)
      if (no_wall_here) cycle
      if (x_this * dir < x * dir) cycle
      x = x_this
      iw_max = iw2
    enddo

    x2 = wall%pt(ip)%x + (s - wall%pt(ip)%s) * (wall%pt(ip+1)%x - wall%pt(ip)%x) / (wall%pt(ip+1)%s - wall%pt(ip)%s)
    if (abs(x - x2) > dx_max) then
      call add_this_point (wall, ip+1, s, x, '', iw_max, wall%pt(ip+1)%phantom)
      exit
    endif
  enddo
enddo

!

wall%n_pt_max = ubound(wall%pt, 1)

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
  allocate (wall%pt(0:0))
else
  n = ubound(wall%pt, 1)
  call move_alloc(wall%pt, pt)
  allocate(wall%pt(0:n+1))
  wall%pt(0:ix_pt-1) = pt(0:ix_pt-1)
  wall%pt(ix_pt+1:) = pt(ix_pt:)
endif

wall%pt(ix_pt)%s = s
wall%pt(ix_pt)%x = x
wall%pt(ix_pt)%name = name
wall%pt(ix_pt)%ix_seg = ix_wall  ! Temp storage of sub-chamber index
wall%pt(ix_pt)%phantom = phantom

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
integer dir, ixs, ix_ele, ix_wall
logical no_wall_here, err_flag

!

ix_ele = element_at_s (branch, s, .true.)
ele => branch%ele(ix_ele)
x = 0

position = [dir * 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, s - (ele%s - ele%value(l$)), 0.0_rp]
x2 = wall3d_d_radius (position, ele, ix_wall, dw_perp, ixs, no_wall_here, origin)
if (no_wall_here) return

y1 = this_y(10.0_rp)
y2 = this_y(-10.0_rp)

if (y1*y2 >= 0) then ! Does not straddle zero
  no_wall_here = .true.
  return
endif

y0 = super_zbrent (this_y, -10.0_rp, 10.0_rp, 1d-8, err_flag)
if (err_flag) call err_exit
y1 = this_y(y0)

x = x_wall

!--------------------------------------------------------------------------------
contains

function this_y (y_in) result (y_wall)

real(rp), intent(in) :: y_in
real(rp) y_wall, dw_perp(3), origin(3), dr, r_part, r_wall
integer ixs
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
