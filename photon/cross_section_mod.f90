module cross_section_mod

use bmad_struct

interface re_associate
  module procedure re_associate_cross_section_array
end interface

interface re_allocate
  module procedure re_allocate_vertex_array
end interface

contains

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine re_allocate_vertex_array (v, n, exact)
!
! Routine to reallocate an array of vertex structures.
! Overloaded by re_allocate.
!
! Modules needed:
!   use cross_section_mod
!
! Input:
!   v(:)  -- cross_section_vertex_struct, allocatable: Array of vertices
!   n     -- Integer: Minimum size needed for array.
!   exact -- Logical, optional: If present and False then the size of
!                    the output array is permitted to be larger than n.
!                    Default is True.
!
! Output:
!   v(:)  -- Cross_section_vertex_struct: Allocated array.
!-

subroutine re_allocate_vertex_array (v, n, exact)

implicit none

type (cross_section_vertex_struct), allocatable :: v(:), temp_v(:)

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact

!

if (allocated(v)) then
  n_old = size(v)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  n_save = min(n, n_old)
  allocate (temp_v(n_save))
  temp_v = v(1:n_save)
  deallocate (v)
  allocate (v(n))
  v(1:n_save) = temp_v
  deallocate (temp_v)
else
  allocate (v(n))
endif

end subroutine re_allocate_vertex_array

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine re_associate_cross_section_array (cross, n, exact)
!
! Routine to reassociate an array of cross_section structures.
! Overloaded by re_associate.
!
! Modules needed:
!   use cross_section_mod
!
! Input:
!   cross(:) -- cross_section_struct, pointer: Array of vertices
!   n        -- Integer: Minimum size needed for array.
!   exact    -- Logical, optional: If present and False then the size of
!                    the output array is permitted to be larger than n.
!                    Default is True.
!
! Output:
!   criss(:) -- Cross_section_struct, pointer: Associated array.
!-

subroutine re_associate_cross_section_array (cross, n, exact)

implicit none

type (cross_section_struct), pointer :: cross(:), temp_cross(:)

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact

!

if (associated(cross)) then
  n_old = size(cross)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  n_save = min(n, n_old)
  temp_cross => cross 
  allocate (cross(n))
  cross(1:n_save) = temp_cross
  deallocate (temp_cross)
else
  allocate (cross(n))
endif

end subroutine re_associate_cross_section_array

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine cross_section_initializer (cross, err)
!
! Routine to initialize a cross_section_struct:
!   1) Add vertex points if there is symmetry.
!   2) Compute circular and elliptical centers.
!   3) Add an end vertex which is the same as the first vertex.
!   4) Calc s_spline(3)
!
! Modules needed:
!   use cross_section_mod
!
! Input:
!   cross  -- Cross_section_struct: Cross-section.
!   
! Output:
!   cross  -- Cross_section_struct: Initialized cross-section.
!   err    -- Logical: Set true if there is a problem. EG: Convex section, etc.
!-

subroutine cross_section_initializer (cross, err)

implicit none

type (cross_section_struct), target :: cross
type (cross_section_vertex_struct), pointer :: v(:)

integer i, n, nn

logical err

character(40) :: r_name = 'cross_section_initializer'

! Init

err = .true.
v => cross%v
n = cross%n_vertex_input

! Calculate s_spline(3)

cross%s_spline(3) = 1 - cross%s_spline(1) - cross%s_spline(2)

! Single vertex is special.

if (n == 1 .and. v(1)%radius_x /= 0) then
  v(1)%x0 = v(1)%x; v(1)%y0 = v(1)%y
  err = .false.
  return
endif

! Compute angle

do i = 1, n
  v(i)%angle = atan2(v(i)%y, v(i)%x)
  if (i == 1) cycle
  if (v(i)%angle <= v(i-1)%angle) v(i)%angle = v(i)%angle + twopi
  if (v(i)%angle >= v(i-1)%angle + pi .or. v(i)%angle <= v(i-1)%angle) then
    call out_io (s_error$, r_name, 'CROSS-SECTION VERTEX NOT IN CLOCKWISE ORDER: (\2F10.5\)', &
                 r_array = [v(i)%x, v(i)%y])
    return
  endif
enddo

if (v(n)%angle - v(1)%angle > twopi) then
  call out_io (s_error$, r_name, 'CROSS-SECTION WINDS BY MORE THAN 2PI!')
  return
endif

! If all (x, y) are in the first quadrent then assume left/right symmetry and 
! propagate vertices to the second quadrent.

if (all(v(1:n)%x >= 0) .and. all(v(1:n)%y >= 0)) then
  if (v(n)%x == 0) then
    nn = 2*n - 1
    call re_allocate(cross%v, nn, .false.); v => cross%v
    v(n+1:nn) = v(n-1:1:-1)
  else
    nn = 2*n
    call re_allocate(cross%v, nn, .false.); v => cross%v
    v(n+1:nn) = v(n:1:-1)
  endif
  v(n+1:nn)%x = -v(n+1:nn)%x
  v(n+1:nn)%angle = pi - v(n+1:nn)%angle
  n = nn
endif

! If everything is in the upper half plane assume up/down symmetry and
! propagate vertices to the bottom half.

if (all(v(1:n)%y >= 0)) then
  if (v(n)%y == 0) then  ! Do not duplicate v(n) vertex
    nn = 2*n - 1
    call re_allocate(cross%v, nn, .false.); v => cross%v
    v(n+1:nn) = v(n-1:1:-1)
  else
    nn = 2*n ! Total number of vetices
    call re_allocate(cross%v, nn, .false.); v => cross%v
    v(n+1:nn) = v(n:1:-1)
  endif

  v(n+1:nn)%y     = -v(n+1:nn)%y
  v(n+1:nn)%angle = twopi - v(n+1:nn)%angle

  if (v(1)%y == 0) nn = nn - 1  ! Do not duplicate v(1) vertex
  n = nn
endif

! Add end vertex which is the same as v(1).

if (v(n)%y /= v(1)%y .or. v(n)%x /= v(1)%x) n = n + 1
call re_allocate(cross%v, n); v => cross%v
v(n) = v(1)
v(n)%angle = v(1)%angle + twopi

! Calculate center of circle/ellipses...

err = .false.

do i = 1, n-1
  call calc_vertex_center (v(i), v(i+1), err)
  if (err) return
enddo
call calc_vertex_center (v(n), v(1), err)

!----------------------------------------------------------------------------
contains

subroutine calc_vertex_center (v1, v2, err)

type (cross_section_vertex_struct) v1, v2

real(rp) x1, y1, x2, y2, x, y
real(rp) x_mid, y_mid, dx, dy
real(rp) a, a2, ct, st

logical err

! If straight line nothing to be done
if (v2%radius_x == 0) return

! Convert (x, y) into unrotated frame if tilted ellipse

x1 = v1%x; y1 = v1%y
x2 = v2%x; y2 = v2%y

if (v2%tilt /= 0) then
  ct = cos(v2%tilt); st = sin(v2%tilt)
  x1 =  ct * v1%x + st * v1%y
  y1 = -st * v1%x + ct * v1%y
  x2 =  ct * v2%x + st * v2%y
  y2 = -st * v2%x + ct * v2%y
endif

! If ellipse then shrink y-axis

if (v2%radius_y /= 0) then
  y1 = y1 * v2%radius_x / v2%radius_y
  y2 = y2 * v2%radius_x / v2%radius_y
endif

! Find center of circle

x_mid = (x1 + x2)/2; y_mid = (y1 + y2)/2
dx    = (x2 - x1)/2; dy    = (y2 - y1)/2

! Find center

a2 = (v2%radius_x**2 - dx**2 - dy**2) / (dx**2 + dy**2)
if (a2 < 0) then
  call out_io (s_error$, r_name, 'CROSS-SECTION VERTEX POINTS TOO FAR APART FOR CIRCLE OR ELLIPSE')
  err = .true.
  return
endif

a = sqrt(a)
if (x_mid * dy > y_mid * dx) a = -a
v2%x0 = x_mid + a * dy
v2%y0 = y_mid - a * dx

! Scale back if radius_y /= 0

if (v2%radius_y /= 0) then
  v2%y0 = v2%y0 * v2%radius_y / v2%radius_x
endif

! Rotate back if tilt /= 0

if (v2%tilt /= 0) then
  x = v2%x0; y = v2%y0
  v2%x0 = ct * x - st * y
  v2%y0 = st * x + ct * y
endif

end subroutine calc_vertex_center

end subroutine cross_section_initializer

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine calc_wall_radius (v, theta, r_wall, gradient)
!
! Routine to calculate the wall radius at a given angle for a given cross-section
! Additionally, the transverse directional derivative is calculated.
!
! Module needed:
!   use cross_setction_mod
!
! Input:
!   v(:)        -- cross_section_vertex_struct: Array of vertices.
!   theta       -- Real(rp): Angle to evaluate at.
!
! Output:
!   r_wall      -- Real(rp): Wall radius at given angle
!   gradient(2) -- Real(rp): (dx, dy) directional 
!-

subroutine calc_wall_radius (v, theta, r_wall, gradient)

implicit none

type (cross_section_vertex_struct), target :: v(:)
type (cross_section_vertex_struct), pointer :: v1, v2


real(rp) theta, r_wall, gradient(2)
real(rp) angle, numer, denom, ct, st, x0, y0, gx, gy, a, b, c
real(rp) cos_ang, sin_ang

integer ix

! Bracket index

angle = theta
if (angle < v(1)%angle) angle = ceiling((v(1)%angle-angle)/twopi) * twopi + angle
call bracket_index (v%angle, 1, size(v), angle, ix)

v1 => v(ix)
v2 => v(ix+1)

! Straight line case

if (v2%radius_x == 0) then
  numer = (v1%x * v2%y - v1%y * v2%x)
  denom = (cos(theta) * (v2%y - v1%y) - sin(theta) * (v2%x - v1%x))
  r_wall = numer / denom
  gradient =  [v2%y - v1%y, -(v2%x - v1%x)] * numer / denom**2
  return
endif

! Convert into unrotated frame if tilted ellipse

x0 = v2%x0; y0 = v2%y0

if (v2%tilt /= 0) then
  angle = angle - v2%tilt
  ct = cos(v2%tilt); st = sin(v2%tilt)
  x0 =  ct * v2%x0 + st * v2%y0
  y0 = -st * v2%x0 + ct * v2%y0
endif

! If ellipse then shrink along y-axis

if (v2%radius_y /= 0) then
  y0 = y0 * v2%radius_x / v2%radius_y
  angle = atan2(sin(angle) * v2%radius_x / v2%radius_y, cos(angle))
endif

cos_ang = cos(angle)
sin_ang = sin(angle)

! Find wall point and derivative

a = 1
b = 2 * (cos_ang * x0 + sin_ang * y0)
c = x0**2 + y0**2
if (v2%radius_x < 0) then
  r_wall = (b - sqrt(b**2 - 4 * a * c)) / (2 * a)
else
  r_wall = (b + sqrt(b**2 - 4 * a * c)) / (2 * a)
endif
gx = r_wall * cos_ang - x0
gy = r_wall * sin_ang - y0
gradient = [gx, gy] * r_wall / (cos_ang * gx + sin_ang * gy)

! If an ellipse then expand along y-axis

if (v2%radius_y /= 0) then
  r_wall = r_wall * sqrt(cos_ang**2 + (sin_ang * v2%radius_y / v2%radius_x)**2)
  gradient(2) = gradient(2) * v2%radius_y / v2%radius_x
  gradient = gradient / sqrt(1 + (v2%radius_y / v2%radius_x)**2)
endif

! If a tilted ellipse then rotate

if (v2%tilt /= 0) then
  ct = cos(v2%tilt); st = sin(v2%tilt)
  gradient =  [ct * gradient(1) -  st * gradient(2), &
               st * gradient(1) +  ct * gradient(2)]
endif

end subroutine calc_wall_radius

end module
