module wall3d_mod

use re_allocate_mod
use bmad_base_mod

!

integer, parameter :: anchor_beginning$ = 1, anchor_center$ = 2, anchor_end$ = 3
character(12) :: anchor_pt_name(0:3) = ['GARBAGE! ', 'Beginning', 'Center   ', 'End      ']

! Structures for defining cross-sections of beam pipes and capillaries
! A cross-section is defined by an array v(:) of wall3d_section_vertex_structs.
! Each vertex v(i) defines a point on the pipe/capillary.
! Vertices are connected by straight lines, circular arcs, or ellipses.
! The radius and tilt values are for the arc from the preceding vertex to this one.
! For v(1), the radius and tilt values are for the arc between v(n) and v(1) where
!   n = upper bound of v(:) array.

type wall3d_vertex_struct
  real(rp) x, y             ! Coordinates of the vertex.
  real(rp) :: radius_x = 0  ! Radius of arc or ellipse x-axis half width. 0 => Straight line.
  real(rp) :: radius_y = 0  ! Ellipse y-axis half height. 
  real(rp) :: tilt = 0      ! Tilt of ellipse
  real(rp) angle            ! Angle of (x, y) point.
  real(rp) x0, y0           ! Center of ellipse
end type

! A beam pipe or capillary cross section is a collection of vertexes.
! Vertices are always ordered in increasing angle.

type wall3d_section_struct
  integer type
  real(rp) :: s = 0                     ! Longitudinal position
  integer n_vertex_input                ! Number of vertices specified by the user.
  type (wall3d_vertex_struct), allocatable :: v(:) 
                                        ! Array of vertices
  ! Center of wall spline
  real(rp) :: x0 = 0, y0 = 0            ! Center of wall
  real(rp) :: dx0_ds = 0                ! Center of wall derivative
  real(rp) :: dy0_ds = 0                ! Center of wall derivative
  real(rp) :: x0_coef(0:3) = 0          ! Spline coefs for x-center
  real(rp) :: y0_coef(0:3) = 0          ! Spline coefs for y-center
  ! Wall radius spline
  real(rp) :: dr_ds = real_garbage$  ! derivative of wall radius 
  real(rp) :: p1_coef(3) = 0            ! Spline coefs for p0 function
  real(rp) :: p2_coef(3) = 0            ! Spline coefs for p1 function

end type

! If, say, %ele_anchor_pt = center$ then center of wall is at the center of the element.

type wall3d_struct
  integer :: ele_anchor_pt = anchor_beginning$      ! anchor_beginning$, anchor_center$, or anchor_end$
  type (wall3d_section_struct), pointer :: section(:) => null()  
end type  

!

interface re_associate
  module procedure re_associate_section_array
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
!   use wall3d_mod
!
! Input:
!   v(:)  -- wall3d_vertex_struct, allocatable: Array of vertices
!   n     -- Integer: Minimum size needed for array.
!   exact -- Logical, optional: If present and False then the size of
!                    the output array is permitted to be larger than n.
!                    Default is True.
!
! Output:
!   v(:)  -- Wall3d_vertex_struct: Allocated array.
!-

subroutine re_allocate_vertex_array (v, n, exact)

implicit none

type (wall3d_vertex_struct), allocatable :: v(:), temp_v(:)

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
! Subroutine re_associate_section_array (section, n, exact)
!
! Routine to reassociate an array of wall3d%section(:).
! Overloaded by re_associate.
!
! Modules needed:
!   use wall3d_mod
!
! Input:
!   section(:) -- wall3d_section_struct, pointer: Array of vertices
!   n        -- Integer: Minimum size needed for array.
!   exact    -- Logical, optional: If present and False then the size of
!                    the output array is permitted to be larger than n.
!                    Default is True.
!
! Output:
!   section(:) -- Wall3d_section_struct, pointer: Associated array.
!-

subroutine re_associate_section_array (section, n, exact)

implicit none

type (wall3d_section_struct), pointer :: section(:), temp_section(:)

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact

!

if (n == 0) then
  if (.not. associated(section)) return
  deallocate(section)

elseif (associated(section)) then
  n_old = size(section)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  n_save = min(n, n_old)
  temp_section => section 
  allocate (section(n))
  section(1:n_save) = temp_section
  deallocate (temp_section)

else
  allocate (section(n))
endif

end subroutine re_associate_section_array


!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine wall3d_initializer (wall3d, err)
!
! Routine to initialize a wall3d_struct
!   1) Add vertex points if there is symmetry.
!   2) Compute circular and elliptical centers.
!   3) Compute spline coefficients, etc.
!
! Modules needed:
!   use wall3d_mod
!
! Input:
!   wall3d -- wall3d_struct: Wall.
!   
! Output:
!   wall3d -- wall3d_struct: Initialized wall.
!   err    -- Logical: Set true if there is a problem.
!-

subroutine wall3d_initializer (wall3d, err)

implicit none

type (wall3d_struct), target :: wall3d
type (wall3d_section_struct), pointer :: s1, s2

real(rp) r1_ave, r2_ave, cos_ang, sin_ang, r, ds, dr_dtheta, a1, a2

integer i, j, n_ave

logical err

! initialize the cross-sections

do i = 1, size(wall3d%section)
  call wall3d_section_initializer(wall3d%section(i), err)
  if (err) return
enddo

! Calculate p0 and p1 spline coefs 

do i = 1, size(wall3d%section) - 1
  s1 => wall3d%section(i)
  s2 => wall3d%section(i+1)

  ! Only do the calc if dr_ds has been set on both sections.
  if (s1%dr_ds == real_garbage$ .or. s2%dr_ds == real_garbage$) cycle

  ! calc average radius
  
  r1_ave = 0; r2_ave = 0
  n_ave = 100
  do j = 1, n_ave
    cos_ang = cos(j * twopi / n_ave)
    sin_ang = sin(j * twopi / n_ave)
    call calc_wall_radius(s1%v, cos_ang, sin_ang, r, dr_dtheta)
    r1_ave = r1_ave + r / n_ave
    call calc_wall_radius(s2%v, cos_ang, sin_ang, r, dr_dtheta)
    r2_ave = r2_ave + r / n_ave
  enddo

  ! Calc coefficients

  ds = s2%s - s1%s
  a1 = s1%dr_ds * ds - (r2_ave - r1_ave)  
  a2 = s2%dr_ds * ds - (r2_ave - r1_ave)  

  s1%p1_coef = [a1, -2*a1-a2, a1+a2] / (2 * r1_ave)
  s1%p2_coef = [a1, -2*a1-a2, a1+a2] / (2 * r2_ave)

enddo


end subroutine wall3d_initializer

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine wall3d_section_initializer (section, err)
!
! Routine to initialize a wall3d_section_struct:
!   1) Add vertex points if there is symmetry.
!   2) Compute circular and elliptical centers.
!
! Modules needed:
!   use wall3d_mod
!
! Input:
!   section  -- Wall3d_section_struct: Wall3d section.
!   
! Output:
!   section  -- Wall3d_section_struct: Initialized section-section.
!   err    -- Logical: Set true if there is a problem.
!-

subroutine wall3d_section_initializer (section, err)

implicit none

type (wall3d_section_struct), target :: section
type (wall3d_vertex_struct), pointer :: v(:)

integer i, n, nn

logical err

character(40) :: r_name = 'wall3d_section_initializer'

! Init

err = .true.
v => section%v
n = section%n_vertex_input

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
    call out_io (s_error$, r_name, 'WALL SECTION VERTEX NOT IN CLOCKWISE ORDER: (\2F10.5\)', &
                 r_array = [v(i)%x, v(i)%y])
    return
  endif

  if (v(i)%radius_x == 0 .and. v(i)%radius_y /= 0) then
    call out_io (s_error$, r_name, 'WALL SECTION VERTEX HAS RADIUS_X = 0 BUT RADIUS_Y != 0 (\2F10.5\)', &
                 r_array = [v(i)%radius_x, v(i)%radius_y])
  endif

  if (v(i)%radius_x * v(i)%radius_y < 0) then
    call out_io (s_error$, r_name, 'WALL SECTION VERTEX HAS RADIUS_X OF DIFFERENT SIGN FROM RADIUS_Y (\2F10.5\)', &
                 r_array = [v(i)%radius_x, v(i)%radius_y])
  endif

enddo

if (v(n)%angle - v(1)%angle >= twopi) then
  call out_io (s_error$, r_name, 'WALL SECTION WINDS BY MORE THAN 2PI!')
  return
endif

! If all (x, y) are in the first quadrent then assume left/right symmetry and 
! propagate vertices to the second quadrent.
! Also radius and tilt info must be moved to the correct vertex.

if (all(v(1:n)%x >= 0) .and. all(v(1:n)%y >= 0)) then
  if (v(n)%x == 0) then
    nn = 2*n - 1
    call re_allocate(section%v, nn, .false.); v => section%v
    v(n+1:nn) = v(n-1:1:-1)
  else
    nn = 2*n
    call re_allocate(section%v, nn, .false.); v => section%v
    v(n+1:nn) = v(n:1:-1)
    v(n+1)%radius_x = 0; v(n+1)%radius_y = 0; v(n+1)%tilt = 0
  endif
  v(n+1:nn)%x     = -v(n+1:nn)%x
  v(n+1:nn)%angle = pi - v(n+1:nn)%angle
  v(nn-n+2:nn)%radius_x = v(n:2:-1)%radius_x
  v(nn-n+2:nn)%radius_y = v(n:2:-1)%radius_y
  v(nn-n+2:nn)%tilt     = -v(n:2:-1)%tilt

  n = nn

endif

! If everything is in the upper half plane assume up/down symmetry and
! propagate vertices to the bottom half.

if (all(v(1:n)%y >= 0)) then
  if (v(n)%y == 0) then  ! Do not duplicate v(n) vertex
    nn = 2*n - 1
    call re_allocate(section%v, nn, .false.); v => section%v
    v(n+1:nn) = v(n-1:1:-1)
  else
    nn = 2*n ! Total number of vetices
    call re_allocate(section%v, nn, .false.); v => section%v
    v(n+1:nn) = v(n:1:-1)
    v(n+1)%radius_x = 0; v(n+1)%radius_y = 0; v(n+1)%tilt = 0
  endif

  v(n+1:nn)%y     = -v(n+1:nn)%y
  v(n+1:nn)%angle = twopi - v(n+1:nn)%angle
  v(nn-n+2:nn)%radius_x = v(n:2:-1)%radius_x
  v(nn-n+2:nn)%radius_y = v(n:2:-1)%radius_y
  v(nn-n+2:nn)%tilt     = -v(n:2:-1)%tilt

  if (v(1)%y == 0) then ! Do not duplicate v(1) vertex
    v(nn)%angle = v(1)%angle
    v(1) = v(nn)
    nn = nn - 1
  endif

  n = nn
  call re_allocate(section%v, n, .true.); v => section%v

! If everything is in the right half plane assume right/left symmetry and
! propagate vertices to the left half.

elseif (all(v(1:n)%x >= 0)) then
  if (v(n)%x == 0) then  ! Do not duplicate v(n) vertex
    nn = 2*n - 1
    call re_allocate(section%v, nn, .false.); v => section%v
    v(n+1:nn) = v(n-1:1:-1)
  else
    nn = 2*n ! Total number of vetices
    call re_allocate(section%v, nn, .false.); v => section%v
    v(n+1:nn) = v(n:1:-1)
    v(n+1)%radius_x = 0; v(n+1)%radius_y = 0; v(n+1)%tilt = 0
  endif

  v(n+1:nn)%x     = -v(n+1:nn)%x
  v(n+1:nn)%angle = twopi - v(n+1:nn)%angle
  v(nn-n+2:nn)%radius_x = v(n:2:-1)%radius_x
  v(nn-n+2:nn)%radius_y = v(n:2:-1)%radius_y
  v(nn-n+2:nn)%tilt     = -v(n:2:-1)%tilt

  if (v(1)%x == 0) then ! Do not duplicate v(1) vertex
    v(nn)%angle = v(1)%angle
    v(1) = v(nn)
    nn = nn - 1
  endif

  n = nn
  call re_allocate(section%v, n, .true.); v => section%v

endif

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

type (wall3d_vertex_struct) v1, v2

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
  call out_io (s_error$, r_name, 'WALL SECTION VERTEX POINTS TOO FAR APART FOR CIRCLE OR ELLIPSE')
  err = .true.
  return
endif

a = sqrt(a2)
if (x_mid * dy > y_mid * dx) a = -a
if (v2%radius_x < 0) a = -a
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

end subroutine wall3d_section_initializer

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine calc_wall_radius (v, cos_ang, sin_ang, r_wall, dr_dtheta, ix_vertex)
!
! Routine to calculate the wall radius at a given angle for a given cross-section
! Additionally, the transverse directional derivative is calculated.
!
! Module needed:
!   use capillary_mod
!
! Input:
!   v(:)         -- wall3d_vertex_struct: Array of vertices that make up the cross-section.
!   cos_ang      -- Real(rp): cosine of the transverse photon position.
!   sin_ang      -- Real(rp): sine of the transverse photon position.
!
! Output:
!   r_wall      -- Real(rp): Wall radius at given angle.
!   dr_dtheta   -- Real(rp): derivative of r_wall.
!   ix_vertex   -- Integer, optional: Wall at given angle is between v(ix_vertex) and
!                    either v(ix_vertex+1) or v(1) if ix_vertex = size(v).
!-

subroutine calc_wall_radius (v, cos_ang, sin_ang, r_wall, dr_dtheta, ix_vertex)

implicit none

type (wall3d_vertex_struct), target :: v(:)
type (wall3d_vertex_struct), pointer :: v1, v2

real(rp) r_wall, dr_dtheta, rx, ry, da, db, angle
real(rp) numer, denom, ct, st, x0, y0, a, b, c
real(rp) cos_ang, sin_ang, radx, cos_a, sin_a, det
real(rp) r_x, r_y, dr_x, dr_y, cos_phi, sin_phi

integer, optional :: ix_vertex
integer ix

! Bracket index if there is more than one vertex
! If there is only one vertex then must be an ellipse or circle

angle = atan2(sin_ang, cos_ang)

if (size(v) == 1) then
  v2 => v(1)
  if (present(ix_vertex)) ix_vertex = 1
else
  if (angle < v(1)%angle) angle = ceiling((v(1)%angle-angle)/twopi) * twopi + angle
  call bracket_index (v%angle, 1, size(v), angle, ix)

  v1 => v(ix)
  if (present(ix_vertex)) ix_vertex = ix

  if (ix == size(v)) then
    v2 => v(1)
  else
    v2 => v(ix+1)
  endif
endif

! Straight line case

if (v2%radius_x == 0) then
  numer = (v1%x * v2%y - v1%y * v2%x)
  denom = (cos_ang * (v2%y - v1%y) - sin_ang * (v2%x - v1%x))
  r_wall = numer / denom
  dr_dtheta = numer * (sin_ang * (v2%y - v1%y) + cos_ang * (v2%x - v1%x)) / denom**2
  return
endif

! If ellipse...

if (v2%radius_y /= 0) then

  ! Convert into unrotated frame if tilted ellipse
  if (v2%tilt /= 0) then
    ct = cos(v2%tilt); st = sin(v2%tilt)
    x0 =  ct * v2%x0 + st * v2%y0
    y0 = -st * v2%x0 + ct * v2%y0
    cos_a = cos_ang * ct + sin_ang * st
    sin_a = sin_ang * ct - cos_ang * st
  else
    x0 = v2%x0; y0 = v2%y0
    cos_a = cos_ang; sin_a = sin_ang
  endif

  rx = v2%radius_x; ry = v2%radius_y
  a = (cos_a/rx)**2 + (sin_a/ry)**2
  b = -2 * (cos_a * x0 / rx**2 + sin_a * y0 / ry**2)
  c = (x0/rx)**2 + (y0/ry)**2 - 1
  radx = sqrt(b**2 - 4 * a * c)

  if (rx > 0) then
    r_wall = (-b + radx) / (2 * a)
  else
    r_wall = (-b - radx) / (2 * a)
  endif

  ! dr/dtheta comes from the equations:
  !   x  = rad_x * cos(phi) + x0
  !   y  = rad_y * sin(phi) + y0
  !   r = sqrt(x^2 + y^2)
  !   Tan(theta) = y/x
 
  r_x = r_wall * cos_a; r_y = r_wall * sin_a
  dr_x = -v2%radius_x * (r_y - y0) / v2%radius_y 
  dr_y =  v2%radius_y * (r_x - x0) / v2%radius_x 
  dr_dtheta = r_wall * (r_x * dr_x + r_y * dr_y) / (r_x * dr_y - r_y * dr_x)

  return
endif

! Else must be a circle.
! Solve for r_wall: (r_wall * cos_a - x0)^2 + (r_wall * sin_a - y0)^2 = radius^2
! dr/dtheta comes from the equations:
!   x = x0 + radius * cos(phi)
!   y = y0 + radius * sin(phi)
!   r = sqrt(x^2 + y^2)
!   Tan(theta) = y/x
! Then
!   dr_vec = (dx, dy) = (-radius * sin(phi), radius * cos(phi)) * dphi
!   dr/dtheta = r * (r_vec dot dr_vec) / (r_vec cross dr_vec)

x0 = v2%x0; y0 = v2%y0

a = 1
b = -2 * (cos_ang * x0 + sin_ang * y0)
c = x0**2 + y0**2 - v2%radius_x**2
radx = sqrt(b**2 - 4 * a * c)

if (v2%radius_x > 0) then
  r_wall = (-b + radx) / (2 * a)
else
  r_wall = (-b - radx) / (2 * a)
endif

r_x = r_wall * cos_ang; r_y = r_wall * sin_ang
dr_x = -(r_y - y0);    dr_y = r_x - x0

dr_dtheta = r_wall * (r_x * dr_x + r_y * dr_y) / (r_x * dr_y - r_y * dr_x)

end subroutine calc_wall_radius

end module
