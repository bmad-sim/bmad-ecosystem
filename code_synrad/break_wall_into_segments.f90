!-
! Subroutine break_wall_into_segments (wall, seg_len_max, branch, seg_len_phantom_max)
!
! Routine to break a wall into segments and do other bookkeeping.
!
! Input:
!   wall                -- wall_struct: Wall to break.
!   seg_len_max         -- real(rp): Maximum segment length.
!   branch              -- branch_struct: lattice branch
!   seg_len_phantom_max -- real(rp), optional: If present then use this number for phantom segments
!                             instead of seg_len_max.
!
! Output:
!   wall        -- wall_struct: Broken wall.
!-

subroutine break_wall_into_segments (wall, seg_len_max, branch, seg_len_phantom_max)

use synrad_struct
use synrad_interface, except => break_wall_into_segments

implicit none

type (wall_struct), target :: wall
type (branch_struct) branch
type (floor_position_struct) floor, local
type (wall_seg_struct), pointer :: seg, seg0
type (wall_pt_struct), pointer :: pt0, pt1
type (ele_struct), pointer :: ele1

integer i_seg, ip, n_seg, ix, n, isg, ix0, ix1, status
real(rp), optional :: seg_len_phantom_max
real(rp) seg_len_max, wall_len, rr, s_max, s_min, dr(3), r1(3), r0(3)
real(rp) theta0, theta1, dx, dz, alpha, beta, dlen
logical err_flag, next_to_alley, has_patch
character(*), parameter :: r_name = 'break_wall_into_segments'

! If a closed geometry: Move last point transverse global position to match first point.

do ip = 0, wall%n_pt_max
  floor = coords_curvilinear_to_floor ([wall%pt(ip)%x, 0.0_rp, wall%pt(ip)%s], branch, err_flag)
  wall%pt(ip)%r_floor = floor%r
enddo

if (branch%param%geometry == closed$) then
  dr = wall%pt(0)%r_floor - wall%pt(wall%n_pt_max)%r_floor
  if (norm2(dr(1:3:2)) > 0.01) then
    call out_io (s_warn$, r_name, trim(wall_name(wall%side)) // ' wall not exactly closed in Global coordinate system.', &
                      'Adjusting position of last wall point on: ' // wall_name(wall%side), &
                      'Wall point is adjusted by dx, ds = \2f11.6\ ', r_array = [dr(1), dr(3)])
  endif
  wall%pt(wall%n_pt_max)%r_floor = wall%pt(0)%r_floor
endif

! If the distance between two wall points is too large then break the interval
! into more than 1 segment.
! x_seg and s_seg are the (x,s) coords of the center of the segment

! First count the number of segments needed and allocate

n = 0
do ip = 1, wall%n_pt_max
  wall_len = sqrt((wall%pt(ip)%x - wall%pt(ip-1)%x)**2 + (wall%pt(ip)%s - wall%pt(ip-1)%s)**2)
  if (present(seg_len_phantom_max) .and. wall%pt(ip)%phantom) then
    n_seg = 1 + wall_len / seg_len_phantom_max
  else
    n_seg = 1 + wall_len / seg_len_max
  endif
  n = n + n_seg
enddo

if (allocated(wall%seg)) deallocate(wall%seg)
allocate (wall%seg(0:n))

wall%seg(0)%len         = 0   ! Dummy segment to hold position info
wall%seg(0)%x           = wall%pt(0)%x
wall%seg(0)%s           = wall%pt(0)%s
wall%seg(0)%r_floor     = wall%pt(0)%r_floor
wall%seg(0)%r_floor_mid = wall%pt(0)%r_floor
wall%seg(0)%ix_seg      = 0 

! Now fill in the information.

wall%pt(0)%ix_pt = 0
i_seg = 0

do ip = 1, wall%n_pt_max
  pt1 => wall%pt(ip)
  pt0 => wall%pt(ip-1)
  
  pt1%ix_pt = ip
  wall_len = sqrt((pt1%x - pt0%x)**2 + (pt1%s - pt0%s)**2)
  if (present(seg_len_phantom_max) .and. wall%pt(ip)%phantom) then
    n_seg = 1 + wall_len / seg_len_phantom_max
  else
    n_seg = 1 + wall_len / seg_len_max
  endif
  pt1%n_seg = n_seg
  pt1%ix_seg = i_seg

  ! Is there a patch element in this wall section

  ix0 = element_at_s (branch%lat, pt0%s, .false., branch%ix_branch)
  ix1 = element_at_s (branch%lat, pt1%s, .true., branch%ix_branch)

  has_patch = any(branch%ele(ix0:ix1)%key == patch$)

  ! Loop over all segments

  do ix = 1, n_seg
    isg = ix + i_seg
    seg => wall%seg(isg)
    seg0 => wall%seg(isg-1)

    seg%ix_pt = ip
    seg%ix_seg = isg

    if (has_patch) then
      rr = float(ix) / n_seg
      seg%r_floor = pt0%r_floor * (1 - rr) + pt1%r_floor * rr 

      floor%r = seg%r_floor
      local = coords_floor_to_curvilinear (floor, branch%ele((ix0+ix1)/2), ele1, status)
      seg%x = local%r(1)
      seg%s = local%r(3)

    else
      rr = float(ix) / n_seg
      seg%x = pt0%x * (1 - rr) + pt1%x * rr
      seg%s = pt0%s * (1 - rr) + pt1%s * rr

      floor = coords_curvilinear_to_floor ([seg%x, 0.0_rp, seg%s], branch, err_flag)
      seg%r_floor = floor%r

      if (ip == wall%n_pt_max .and. branch%param%geometry == closed$) then
        seg%r_floor = seg%r_floor + dr * real(ix, rp) / n_seg
      endif
    endif

    seg%len = sqrt((seg%r_floor(1) - seg0%r_floor(1))**2 + (seg%r_floor(3) - seg0%r_floor(3))**2)
    seg%r_floor_mid = (seg%r_floor + seg0%r_floor) / 2

    seg%theta = atan2(seg%r_floor(1) - seg0%r_floor(1), seg%r_floor(3) - seg0%r_floor(3))
  end do

  i_seg = i_seg + n_seg
enddo

wall%n_seg_max = i_seg

! Mark points that are near alleys

wall%pt%next_to_alley = .false.

s_max = 0
do ip = 1, wall%n_pt_max
  s_max = max(s_max, wall%pt(ip)%s)
  if (wall%pt(ip)%s < s_max) then
    wall%pt(ip)%next_to_alley = .true.
    wall%pt(ip-1)%next_to_alley = .true.
    wall%pt(ip+1)%next_to_alley = .true.
  endif
enddo

s_min = wall%pt(wall%n_pt_max)%s
do ip = wall%n_pt_max-1, 0, -1
  s_min = min(s_min, wall%pt(ip)%s)
  if (wall%pt(ip)%s > s_min) then
    wall%pt(ip)%next_to_alley = .true.
    wall%pt(ip-1)%next_to_alley = .true.
    wall%pt(ip+1)%next_to_alley = .true.
  endif
enddo

! Calculate the triangle point

do ip = 1, wall%n_pt_max
  pt0  => wall%pt(ip-1)
  pt1  => wall%pt(ip) 
  theta0 = wall%seg(pt1%ix_seg+1)%theta
  theta1 = wall%seg(pt1%ix_seg+pt1%n_seg)%theta

  !

  dx = pt1%r_floor(1) - pt0%r_floor(1)
  dz = pt1%r_floor(3) - pt0%r_floor(3)

  dlen = sqrt(dx**2 + dz**2)

  if (dlen * abs(theta0 - theta1) < synrad_significant_length) then  ! Colinear
    pt1%linear_wall = .true.   ! So don't need the triangle point
    cycle
  endif

  pt1%linear_wall = .false. 

  ! Round off error can be huge when theta1 and theta0 are approx equal.
  ! In this case increase the difference to a safe level.

  if (abs(theta1 - theta0) < 1d-3) then
    theta0 = theta0 - sign(0.5d-3, theta1 - theta0)
    theta1 = theta1 + sign(0.5d-3, theta1 - theta0)
  endif

  ! alpha is distance ("forward" is +) from pt0 to the triangle point. Must be positive.
  ! beta is  distance ("forward" is +) from pt1 to the triangle point. Must be negative.

  alpha =  (dz * sin(theta1) - dx * cos(theta1)) / sin(theta1 - theta0)
  beta  = -(dz * sin(theta0) - dx * cos(theta0)) / sin(theta0 - theta1)

  if (alpha < 0 .or. beta > 0) then
    call out_io (s_fatal$, r_name, 'Problem with triangle point calc!')
    call err_exit
  endif

  r0 = [pt0%r_floor(1) + alpha * sin(theta0), 0.0_rp, pt0%r_floor(3) + alpha * cos(theta0)]
  r1 = [pt1%r_floor(1) + beta * sin(theta1),  0.0_rp, pt1%r_floor(3) + beta * cos(theta1)]
  dr = r1 - r0

  if (abs(dr(1)) > synrad_significant_length .or. abs(dr(3)) > synrad_significant_length) then
    call out_io (s_fatal$, r_name, 'Confused triangle point calc!')
    call err_exit
  endif

  pt1%r_floor_tri = (r0 + r1) / 2

enddo

! Mark point and seg

end subroutine
