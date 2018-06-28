!+
! subroutine track_ray_to_wall (ray, walls)
!
! subroutine to propagate a synch radiation ray until it hits
!    a wall
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   ray     -- ray_struct: synch radiation ray with starting parameters set
!   walls   -- walls_struct: both walls and ends
!
! Output:
!   ray       -- ray_struct: synch radiation ray propagated to wall
!-

subroutine track_ray_to_wall (ray, walls)

use synrad_struct
use synrad_interface, except => track_ray_to_wall

implicit none

type (ray_struct), target :: ray
type (walls_struct), target :: walls
type (wall_struct), pointer :: neg_x_wall, pos_x_wall

integer ix_neg, ix_pos, ix_min_seg_pos, ix_min_seg_neg, ix_pt, dir, itry

real(rp) len_min_pos, len_min_neg, r_neg, r_pos, r
real(rp) theta, s_last_pos, s_last_neg

logical pos_go, neg_go, next_to_alley_pos, next_to_alley_neg, pos_at_end, neg_at_end

! set pointers

pos_x_wall => walls%positive_x_wall
neg_x_wall => walls%negative_x_wall
dir = ray%direction

! ix_neg and ix_pos are the neg_x_wall and pos_x_wall side points that
! are at or just "upstream" (ray%direction points downstream) of the ray.

call get_initial_wall_pt (ray, neg_x_wall, ix_neg)
call get_initial_wall_pt (ray, pos_x_wall, ix_pos)

! Essentially find the segment with the minimum distance from segment to ray start.

ix_min_seg_pos = -1; len_min_pos = 1e10
ix_min_seg_neg = -1; len_min_neg = 1e10

pos_go = .true.; neg_go = .true.

s_last_pos = -1e10 * dir
s_last_neg = -1e10 * dir

next_to_alley_pos = .false.
next_to_alley_neg = .false.

do itry = 1, pos_x_wall%n_pt_max + neg_x_wall%n_pt_max

  ! See if we have hit the end of the machine

  if (dir == 1) then
    pos_at_end = (ix_pos == pos_x_wall%n_pt_max + 1)
    neg_at_end = (ix_neg == neg_x_wall%n_pt_max + 1)

    if (pos_at_end .and. neg_at_end) then
      if (walls%lat_geometry == open$) then
        ray%wall_side = exit_side$   ! End "wall" at end of lattice
        return
      endif
      ix_pos = 1; s_last_pos = 0
      ix_neg = 1; s_last_neg = 0
    endif

  else ! direction = -1
    pos_at_end = (ix_pos == 0)
    neg_at_end = (ix_neg == 0)

    if (pos_at_end .and. neg_at_end) then
      if (walls%lat_geometry == open$) then
        ray%wall_side = start_side$  ! End "wall" at beginning of lattice
        return
      endif
      ix_pos = pos_x_wall%n_pt_max; s_last_pos = walls%s_max
      ix_neg = neg_x_wall%n_pt_max; s_last_neg = walls%s_max
    endif
  endif

  ! Check next pos_x wall point

  if ((pos_go .and. .not. pos_at_end .and. dir * (s_last_pos - s_last_neg) <=  0) .or. .not. neg_go .or. neg_at_end) then
    call check_pt (ray, ix_pos, pos_x_wall, dir, len_min_pos, ix_min_seg_pos, r_pos)
    s_last_pos = pos_x_wall%pt(ix_pos)%s
    if (ix_min_seg_pos /= -1 .and. .not. pos_x_wall%pt(ix_pos)%next_to_alley) pos_go = .false.
    ix_pos = ix_pos + dir

  else
    call check_pt (ray, ix_neg, neg_x_wall, dir, len_min_neg, ix_min_seg_neg, r_neg)
    s_last_neg = neg_x_wall%pt(ix_neg)%s
    if (ix_min_seg_neg /= -1 .and. .not. neg_x_wall%pt(ix_neg)%next_to_alley) neg_go = .false.
    ix_neg = ix_neg + dir
  endif

  ! End if there has been a hit and we are clear of any alleys.

  if (.not. next_to_alley_pos .and. .not. next_to_alley_neg .and. &
                            (ix_min_seg_pos /= -1 .or. ix_min_seg_neg /= -1)) then
    if (len_min_pos < len_min_neg) then
      ray%wall_side  = positive_x$
      ray%track_len  = len_min_pos
      ray%r_seg      = r_pos
      ray%ix_seg_pt  = ix_min_seg_pos
      ray%ix_wall_pt = pos_x_wall%seg(ix_min_seg_pos)%ix_pt
    else
      ray%wall_side  = negative_x$
      ray%track_len  = len_min_neg
      ray%r_seg      = r_neg
      ray%ix_seg_pt  = ix_min_seg_neg
      ray%ix_wall_pt = neg_x_wall%seg(ix_min_seg_neg)%ix_pt
    endif

    theta = ray%start_floor%theta
    ray%now_floor%r = ray%start_floor%r + ray%track_len * [sin(theta), 0.0_rp, cos(theta)]
    ray%now%vec(3) = ray%start%vec(3) + ray%track_len * ray%start%vec(4)

    return
  endif

enddo

print *, 'ERROR: CANNOT FIND WALL INTERSECTION!'
call err_exit


!----------------------------------------------------------------------------------------
contains

subroutine check_pt (ray, ix_pt, wall, dir, len_min, ix_seg_min, r_seg)

type (ray_struct) ray
type (wall_struct), target :: wall
type (wall_pt_struct), pointer :: pt, pt0

real(rp) len_min, r_seg, len_this, r_seg_this, dum1, dum2
integer ix_pt, dir, ix_seg_min, is

! Check if there is an intersection

pt0 => wall%pt(ix_pt-1)
pt => wall%pt(ix_pt)

if (pt%linear_wall) then
  if (.not. line_hit(ray, pt0%r_floor, pt%r_floor, dum1, dum2)) return   ! No hit
else
  if (.not. line_hit(ray, pt0%r_floor, pt%r_floor, dum1, dum2) .and. &
      .not. line_hit(ray, pt0%r_floor, pt%r_floor_tri, dum1, dum2) .and. &
      .not. line_hit(ray, pt%r_floor, pt%r_floor_tri, dum1, dum2)) return
endif

! Possible that ray is hitting wall so test all the segments.

do is = pt%ix_seg+1, pt%ix_seg+pt%n_seg
  if (.not. line_hit (ray, wall%seg(is-1)%r_floor, wall%seg(is)%r_floor, len_this, r_seg_this)) cycle

  if (len_this < len_min) then
    ix_seg_min = is
    len_min    = len_this
    r_seg      = r_seg_this
  endif
enddo

end subroutine

!----------------------------------------------------------------------------------------
! contains

function line_hit (ray, floor0, floor1, track_len, r_seg) result (hitting)

type (wall_struct) wall
type (ray_struct) ray

real(rp) floor0(3), floor1(3), track_len, r_seg
real(rp) r_x, r_z, dr_x, dr_z, denom
real(rp), save :: old_theta = 0, cos_t = 1, sin_t = 0
integer ix_seg
logical hitting

!

r_x = floor0(1) - ray%start_floor%r(1)
r_z = floor0(3) - ray%start_floor%r(3)

dr_x = floor1(1) - floor0(1)
dr_z = floor1(3) - floor0(3)

if (ray%start_floor%theta /= old_theta) then
  cos_t = cos(ray%start_floor%theta)
  sin_t = sin(ray%start_floor%theta)
  old_theta = ray%start_floor%theta
endif

denom = sin_t * dr_z - cos_t * dr_x

hitting = .false.
if (abs(denom) < 1e-10) return

r_seg = (cos_t * r_x - sin_t * r_z) / denom
if (r_seg < 0  .or. r_seg > 1) return

! If the ray is starting from the wall then can get 
! very small track_len values. Don't count these.

track_len = cos_t * (r_z + r_seg * dr_z) + sin_t * (r_x + r_seg * dr_x)
if (track_len < synrad_significant_length) return

hitting = .true.

end function

end subroutine
