module wall3d_mod

use element_at_s_mod

implicit none

!

interface re_allocate
  module procedure re_allocate_wall3d_section_array
  module procedure re_allocate_wall3d_vertex_array
end interface

contains

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine re_allocate_wall3d_vertex_array (v, n, exact)
!
! Routine to reallocate an array of vertex structures.
! Overloaded by re_allocate.
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

subroutine re_allocate_wall3d_vertex_array (v, n, exact)

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

end subroutine re_allocate_wall3d_vertex_array

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine re_allocate_wall3d_section_array (section, n, exact)
!
! Routine to reallocate an array of wall3d%section(:).
! Overloaded by re_allocate.
!
! Input:
!   section(:) -- wall3d_section_struct, pointer: Array of vertices
!   n        -- Integer: Minimum size needed for array.
!   exact    -- Logical, optional: If present and False then the size of
!                    the output array is permitted to be larger than n.
!                    Default is True.
!
! Output:
!   section(:) -- Wall3d_section_struct, pointer: Allocated array.
!-

subroutine re_allocate_wall3d_section_array (section, n, exact)

type (wall3d_section_struct), allocatable :: section(:), temp_section(:)

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact

!

if (n == 0) then
  if (.not. allocated(section)) return
  deallocate(section)

elseif (allocated(section)) then
  n_old = size(section)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  n_save = min(n, n_old)
  call move_alloc (section, temp_section)
  allocate (section(n))
  section(1:n_save) = temp_section
  deallocate (temp_section)

else
  allocate (section(n))
endif

end subroutine re_allocate_wall3d_section_array


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
! Input:
!   wall3d -- wall3d_struct: Wall.
!   
! Output:
!   wall3d -- wall3d_struct: Initialized wall.
!   err    -- Logical: Set true if there is a problem.
!-

subroutine wall3d_initializer (wall3d, err)

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
  a1 = 0; a2 = 0
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
! Input:
!   section  -- Wall3d_section_struct: Wall3d section.
!   
! Output:
!   section  -- Wall3d_section_struct: Initialized section-section.
!   err    -- Logical: Set true if there is a problem.
!-

subroutine wall3d_section_initializer (section, err)

type (wall3d_section_struct), target :: section
type (wall3d_vertex_struct), pointer :: v(:)

real(rp) angle, r_wall, dr_dtheta
integer i, n, nn

logical err

character(60) sec_name
character(*), parameter :: r_name = 'wall3d_section_initializer'

! Init

err = .true.
v => section%v
n = section%n_vertex_input
if (section%name == '') then
  sec_name = 'WALL SECTION'
else
  sec_name = 'WALL SECTION: ' // section%name
endif

! Relative or absolute vertex numbers?
! Vertex numbers are always stored as relative so need to convert if input numbers are absolute.

if (section%vertices_state == absolute$) then
  v%x = v%x - section%r0(1)
  v%y = v%y - section%r0(2)
  section%vertices_state = shifted_to_relative$
endif

! Error check if there is only a single vertex

if (n == 1 .and. v(1)%radius_x == 0 .and. (v(1)%x == 0 .or. v(1)%y == 0)) then
  call out_io (s_error$, r_name, trim(sec_name) // ' WITH SINGLE VERTEX AND ZERO RADIUS HAS X OR Y = 0 WHICH GIVES ZERO CROSS-SECTION AREA!')
  return
endif

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
    call out_io (s_error$, r_name, trim(sec_name) // ' VERTEX NOT IN COUNTER-CLOCKWISE ORDER: (\2F10.5\)', &
                 r_array = [v(i)%x, v(i)%y])
    return
  endif

  if (v(i)%radius_x == 0 .and. v(i)%radius_y /= 0) then
    call out_io (s_error$, r_name, trim(sec_name) // ' VERTEX HAS RADIUS_X = 0 BUT RADIUS_Y != 0 (\2F10.5\)', &
                 r_array = [v(i)%radius_x, v(i)%radius_y])
  endif

  if (v(i)%radius_x * v(i)%radius_y < 0) then
    call out_io (s_error$, r_name, trim(sec_name) // ' VERTEX HAS RADIUS_X OF DIFFERENT SIGN FROM RADIUS_Y (\2F10.5\)', &
                 r_array = [v(i)%radius_x, v(i)%radius_y])
    return
  endif

enddo

if (v(n)%angle - v(1)%angle >= twopi) then
  call out_io (s_error$, r_name, trim(sec_name) // ' WINDS BY MORE THAN 2PI!')
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
    v(n+1)%radius_x = 0; v(n+1)%radius_y = 0; v(n+1)%tilt = 0; v(n+1)%type = normal$
  endif
  v(n+1:nn)%x           = -v(n+1:nn)%x
  v(n+1:nn)%angle       = pi - v(n+1:nn)%angle
  v(nn-n+2:nn)%radius_x =  v(n:2:-1)%radius_x
  v(nn-n+2:nn)%radius_y =  v(n:2:-1)%radius_y
  v(nn-n+2:nn)%tilt     = -v(n:2:-1)%tilt
  v(nn-n+2:nn)%type     =  v(n:2:-1)%type

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
    v(n+1)%radius_x = 0; v(n+1)%radius_y = 0; v(n+1)%tilt = 0; v(n+1)%type = normal$
  endif

  v(n+1:nn)%y           = -v(n+1:nn)%y
  v(n+1:nn)%angle       = twopi - v(n+1:nn)%angle
  v(nn-n+2:nn)%radius_x =  v(n:2:-1)%radius_x
  v(nn-n+2:nn)%radius_y =  v(n:2:-1)%radius_y
  v(nn-n+2:nn)%tilt     = -v(n:2:-1)%tilt
  v(nn-n+2:nn)%type     =  v(n:2:-1)%type

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
    v(n+1)%radius_x = 0; v(n+1)%radius_y = 0; v(n+1)%tilt = 0; v(n+1)%type = normal$
  endif

  v(n+1:nn)%x           = -v(n+1:nn)%x
  v(n+1:nn)%angle       = pi - v(n+1:nn)%angle
  v(nn-n+2:nn)%radius_x =  v(n:2:-1)%radius_x
  v(nn-n+2:nn)%radius_y =  v(n:2:-1)%radius_y
  v(nn-n+2:nn)%tilt     = -v(n:2:-1)%tilt
  v(nn-n+2:nn)%type     =  v(n:2:-1)%type

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
  call calc_vertex_center (v, i, i+1, err)
  if (err) return
enddo
call calc_vertex_center (v, n, 1, err)
if (err) return

! Quick sanity check

do i = 1, 32
  angle = i * twopi / 32
  call calc_wall_radius (v, cos(angle), sin(angle), r_wall, dr_dtheta)
  if (r_wall <= 0) then
    call out_io (s_error$, r_name, trim(sec_name) // ' DOES NOT ENCIRCLE SHAPE CENTER!')
    err = .true.
    return
  endif
enddo

!----------------------------------------------------------------------------
contains

subroutine calc_vertex_center (v, i1, i2, err)

type (wall3d_vertex_struct), target :: v(:)
type (wall3d_vertex_struct), pointer :: v1, v2

real(rp) x1, y1, x2, y2, x, y
real(rp) x_mid, y_mid, dx, dy
real(rp) a, a2, ct, st

integer i1, i2
logical err

! If straight line nothing to be done

v1 => v(i1)
v2 => v(i2)

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
  call out_io (s_error$, r_name, trim(sec_name) // &
                  ' VERTEX POINTS TOO FAR APART FOR CIRCLE OR ELLIPSE FOR VERTICES: \i0\ TO \i0\ ', i_array = [i1, i2])
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
!   use wall3d_mod
!
! Input:
!   v(:)         -- wall3d_vertex_struct: Array of vertices that make up the cross-section.
!   cos_ang      -- Real(rp): cosine of the transverse photon position.
!   sin_ang      -- Real(rp): sine of the transverse photon position.
!
! Output:
!   r_wall      -- Real(rp): Wall radius at given angle.
!   dr_dtheta   -- Real(rp): derivative of r_wall.
!   ix_vertex   -- Integer, optional: Wall at given angle is between v(ix_vertex-1) and
!                    v(ix_vertex). If ix_vertex = 1 then Wall at given angle is between
!                    v(N) and v(1) where N = size(v).
!-

subroutine calc_wall_radius (v, cos_ang, sin_ang, r_wall, dr_dtheta, ix_vertex)

type (wall3d_vertex_struct), target :: v(:)
type (wall3d_vertex_struct), pointer :: v1, v2

real(rp) r_wall, dr_dtheta, rx, ry, da, db, angle
real(rp) numer, denom, ct, st, x0, y0, a, b, c
real(rp) cos_ang, sin_ang, radx, cos_a, sin_a, det
real(rp) r_x, r_y, dr_x, dr_y, cos_phi, sin_phi

integer, optional :: ix_vertex
integer ix, ix2

! Bracket index if there is more than one vertex
! If there is only one vertex then must be an ellipse or circle

if (size(v) == 1) then
  v2 => v(1)
  if (present(ix_vertex)) ix_vertex = 1
else
  angle = atan2(sin_ang, cos_ang)
  if (angle < v(1)%angle) angle = ceiling((v(1)%angle-angle)/twopi) * twopi + angle
  ix = bracket_index (angle, v%angle, 1)

  v1 => v(ix)

  if (ix == size(v)) then
    ix2 = 1
  else
    ix2 = ix+1
  endif
  v2 => v(ix2)
  if (present(ix_vertex)) ix_vertex = ix2
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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function wall3d_d_radius (position, ele, ix_wall, perp, ix_section, 
!                                      no_wall_here, origin, radius_wall, err_flag) result (d_radius)
!
! Routine to calculate the difference radius = particle_radius - wall_radius.
! Radiuses are measured along a line from the wall origin with the line passing through
! the particle position.
! The wall origin itself lies on a line connecting the centers of the bounding sections.
!
! Module needed:
!   use wall3d_mod
!
! Input:
!   position(6)  -- real(rp): Particle position in element coordinates. In a patch, with respect to entrance coords.
!           [position(1), position(3)] = [x, y] transverse coords.
!            position(5)                = Longitudinal position relative to beginning of element.
!            position(6)                = Longitudinal velocity (only +/- sign matters).
!   ele          -- ele_struct: Element with wall
!   ix_wall      -- integer, optional: Index of wall in %wall3d(:) array. Default is 1.
!
! Output:
!   d_radius       -- real(rp), Normalized radius: r_particle - r_wall
!   perp(3)        -- real(rp), optional: Perpendicular normal to the wall.
!   ix_section     -- integer, optional: Set to wall slice section particle is in. 
!                      That is between ix_section and ix_section+1.
!   no_wall_here   -- logical, optional: True if the sub-chamber under consideration does not exist at the 
!                       longitudinal location of the particle.
!   origin(3)      -- real(rp), optional: (x, y, s) origin with respect to the radius is measured.
!                       Uses the same coords as position.
!   radius_wall    -- real(rp), optional: Radius of the wall. 
!   err_flag       -- Logical, optional: Set True if error. (EG noassociated %wall3d), false otherwise.
!-

function wall3d_d_radius (position, ele, ix_wall, perp, ix_section, &
                                      no_wall_here, origin, radius_wall, err_flag) result (d_radius)

type (ele_struct), target :: ele
type (wall3d_section_struct), pointer :: sec1, sec2
type (wall3d_struct), pointer :: wall3d
type (wall3d_vertex_struct), allocatable :: v(:)
type (ele_struct), pointer :: ele1, ele2, ele_ptr
type (floor_position_struct) floor_particle, floor_sec1, floor_sec2
type (floor_position_struct) floor1_w, floor2_w, floor1_dw, floor2_dw, floor1_p, floor2_p
type (floor_position_struct) loc_p, loc1_sec2, loc2_sec1, floor
type (branch_struct), pointer :: branch

real(rp), intent(in) :: position(:)
real(rp), optional :: perp(3), origin(3), radius_wall

real(rp), pointer :: vec(:), value(:)
real(rp) d_radius, r_particle, r_norm, s_rel, spline, cos_theta, sin_theta
real(rp) r1_wall, r2_wall, dr1_dtheta, dr2_dtheta, f_eff, ds, disp, r_wall
real(rp) p1, p2, dp1, dp2, s_particle, ds_offset, x, y, x0, y0, f, d_rad1, d_rad2
real(rp) r(3), r0(3), rw(3), drw(3), dr0(3), p_sec1(3), p_sec2(3), drp(3)
real(rp) dtheta_dphi, alpha, dalpha, beta, dx, dy, w_mat(3,3)
real(rp) s1, s2, r_p(3)

integer i, ix_w, n_slice, n_sec, ix_vertex1, ix_vertex2, status
integer, optional :: ix_section, ix_wall

logical, optional :: err_flag, no_wall_here
logical err, is_branch_wall

character(*), parameter :: r_name = 'wall3d_d_radius' 

! Find the wall definition

if (present(err_flag)) err_flag = .true.
if (present(no_wall_here)) no_wall_here = .false.
if (present(radius_wall)) radius_wall = -1
d_radius = -1

wall3d => pointer_to_wall3d (ele, ix_wall, ds_offset, is_branch_wall)
if (.not. associated(wall3d)) return
if (present(err_flag)) err_flag = .false.

! Init

s_particle = position(5) + ds_offset
n_sec = size(wall3d%section)
branch => pointer_to_branch(ele)

! The outward normal vector is discontinuous at the wall sections.
! If the particle is at a wall section, use the correct interval.
! If moving in +s direction then the correct interval is whith %section(ix_w+1)%s = particle position.

! Case where the particle is outside the wall region longitudinally: Wrap the wall around the branch ends
! if the branch has closed geometry. Otherwise assume a constant cross-section.

! If at a wall start/end then add a slop factor to decide if the particle is outside to avoid round-off
! problems when a particle is at the end of an element and the wall is only defined for that element.
! EG: time_runge_kutta will track beyound an element's edge.

if (s_particle < wall3d%section(1)%s .or. (s_particle == wall3d%section(1)%s .and. position(6) > 0)) then
  if (wall3d%section(1)%type == wall_start$) then
    if (s_particle < wall3d%section(1)%s - bmad_com%significant_length/10) then
      if (present(no_wall_here)) no_wall_here = .true.
    else
      call d_radius_at_section(wall3d%section(1))
    endif
    return

  elseif (branch%param%geometry == closed$) then
    sec1 => wall3d%section(n_sec);  s1 = sec1%s - branch%param%total_length
    sec2 => wall3d%section(1);      s2 = sec2%s
    if (present(ix_section)) ix_section = n_sec

  else
    call d_radius_at_section(wall3d%section(1))
    return
  endif

elseif (s_particle > wall3d%section(n_sec)%s .or. (s_particle == wall3d%section(n_sec)%s .and. position(6) < 0)) then
  if (wall3d%section(n_sec)%type == wall_end$) then
    if (s_particle > wall3d%section(n_sec)%s + bmad_com%significant_length/10) then
      if (present(no_wall_here)) no_wall_here = .true.
    else
      call d_radius_at_section(wall3d%section(n_sec))
    endif
    return

  elseif (branch%param%geometry == closed$) then
    sec1 => wall3d%section(n_sec);  s1 = sec1%s
    sec2 => wall3d%section(1);      s2 = sec2%s + branch%param%total_length
    if (present(ix_section)) ix_section = n_sec

  else
    call d_radius_at_section(wall3d%section(n_sec))
    return
  endif

else

  ! Find the wall points (defined cross-sections) to either side of the particle.
  ! That is, the particle is in the interval [%section(ix_w)%s, %section(ix_w+1)%s].

  ix_w = bracket_index (s_particle, wall3d%section%s, 1)
  if (s_particle == wall3d%section(ix_w)%s .and. (position(6) > 0 .or. ix_w == size(wall3d%section))) ix_w = ix_w - 1

  ! sec1 and sec2 are the cross-sections to either side of the particle.
  ! Calculate the radius values at the cross-sections.

  sec1 => wall3d%section(ix_w);   s1 = sec1%s
  sec2 => wall3d%section(ix_w+1); s2 = sec2%s
  if (present(ix_section)) ix_section = ix_w

endif

! Check if there is a wall here

if (sec1%type == wall_end$ .or. sec2%type == wall_start$) then
  if (present(no_wall_here)) no_wall_here = .true.
  if (sec1%type /= wall_end$ .or. sec2%type /= wall_start$) then
    call out_io (s_error$, r_name, 'WALL SECTION OF TYPE WALL_END NOT FOLLOWED BY A WALL SECTION OF TYPE WALL_START')
    if (present(err_flag)) err_flag = .true.
  endif
  return
endif

! Zero distance between sections case

ds = s2 - s1

if (s_particle == s1 .and. ds == 0) then
  ! Choose the greater d_radius. That is, the particle is considered outside if it is outside either section.
  call d_radius_at_section(sec1, d_rad1)
  call d_radius_at_section(sec2, d_rad2)
  if (d_rad1 > d_rad2) call d_radius_at_section(sec1)

  if (position(6) > 0) then
    if (present(perp)) perp = [0.0_rp, 0.0_rp, 1.0_rp]    
  else
    if (present(perp)) perp = [0.0_rp, 0.0_rp, -1.0_rp]    
  endif

  return
endif

!----------------------------
! If we are in a region with a patch then the geometry is more complicated since the section planes
! may not be perpendicular to the z-axis

cos_theta = 1  ! Default if r_particle == 0
sin_theta = 0

if (sec2%patch_in_region) then
  ! ele1 and ele2 are lattice elements of sec1 and sec2
  ele1 => pointer_to_ele (branch%lat, sec1%ix_ele, sec1%ix_branch)
  ele2 => pointer_to_ele (branch%lat, sec2%ix_ele, sec2%ix_branch)

  ! floor_particle is coordinates of particle in global reference frame
  ! floor_sec1 is sec1 origin in global ref frame
  ! floor_sec2 is sec2 origin in global ref frame
  floor_particle = this_coords_to_floor (ele, [position(1), position(3), position(5)])
  floor_sec1     = this_coords_to_floor (ele1, [sec1%r0, sec1%s - ele1%s_start])
  floor_sec2     = this_coords_to_floor (ele2, [sec2%r0, sec2%s - ele2%s_start])

  ! loc_p  is coordinates of particle in sec1 origin ref frame
  ! loc2_sec1 is coordinates of sec2 origin in sec1 orgin ref frame
  loc_p     = coords_floor_to_relative (floor_sec1, floor_particle, .false.)
  loc2_sec1 = coords_floor_to_relative (floor_sec1, floor_sec2, .false.)

  ! sec1, sec2, and particle origins form a triangle. There is a half plane that contains this triangle with the 
  ! edge of the half plane along the sec1-sec2 origin line and with the half plane containing the particle origin.
  ! The half plane intersects the sec1 wall curve at one point. This point, along with the projection of the 
  ! particle origin onto the sec1 plane along the half plane edge, and the sec1 origin lie on a line and define
  ! the radius of the particle at sec1.
  dr0 = loc2_sec1%r / norm2(loc2_sec1%r)  ! sec2 - sec1 origins normalized vector
  s1 = loc_p%r(3) / dr0(3)            ! Distance of particle from s1
  p_sec1 = loc_p%r - s1 * dr0         ! Particle origin projected to sec1. Should have p_sec1(3) = 0
  r_norm = norm2(p_sec1(1:2))
  if (r_norm /= 0) then
    cos_theta = p_sec1(1) / r_norm
    sin_theta = p_sec1(2) / r_norm
  endif
  call calc_wall_radius (sec1%v, cos_theta, sin_theta, r1_wall, dr1_dtheta, ix_vertex1)

  ! floor1_p is the projected particle origin onto the sec1 plane expressed in global coords.
  ! floor1_w is sec1 wall pt expressed in global coords.

  floor1_p = coords_relative_to_floor (floor_sec1, p_sec1)
  floor1_w = coords_relative_to_floor (floor_sec1, r1_wall * [cos_theta, sin_theta, 0.0_rp])

  ! floor1_dw is sec1 wall pt derivative with respect to theta in global reference frame. 
  ! dtheta_dphi is change in local sec1 angle (theta) with respect to global angle (phi).

  if (present(perp)) then
    beta = sqrt( (dr0(2)**2 + dr0(3)**2) / (dr0(1)**2 + dr0(3)**2) )
    dtheta_dphi = (beta**2 * sin_theta**2 + cos_theta**2) / beta
    r = dtheta_dphi * (dr1_dtheta * [cos_theta, sin_theta, 0.0_rp] + r1_wall * [-sin_theta, cos_theta, 0.0_rp])
    floor1_dw = coords_relative_to_floor (floor_sec1, r)
    floor1_dw%r = floor1_dw%r - floor_sec1%r
  endif

  ! loc_p  is coordinates of particle in sec1 origin ref frame
  ! loc_1_sec2 is coordinates of sec1 origin in sec2 origin ref frame
  loc_p = coords_floor_to_relative (floor_sec2, floor_particle, .false.)
  loc1_sec2 = coords_floor_to_relative (floor_sec2, floor_sec1, .false.)

  ! Find wall radius for sec2. See above.
  dr0 = loc1_sec2%r / norm2(loc1_sec2%r)  ! sec1 - sec2 origins normalized vector
  s2 = loc_p%r(3) / dr0(3)
  p_sec2 = loc_p%r - s2 * dr0   ! particle origin projected onto sec2 plane. Should have p_sec2(3) = 0
  r_norm = norm2(p_sec2(1:2))
  if (r_norm /= 0) then
    cos_theta = p_sec2(1) / r_norm
    sin_theta = p_sec2(2) / r_norm
  endif
  call calc_wall_radius (sec2%v, cos_theta, sin_theta, r2_wall, dr2_dtheta, ix_vertex2)

  ! floor2_p is the particle coords projected onto the sec2 plane in global ref frame
  ! floor2_w  is sec2 wall pt in global reference frame

  floor2_p = coords_relative_to_floor (floor_sec2, p_sec2)
  floor2_w = coords_relative_to_floor (floor_sec2, r2_wall * [cos_theta, sin_theta, 0.0_rp])

  ! floor2_dw is sec2 wall pt derivative with respect to theta in global reference frame
  ! dtheta_dphi is change in local sec1 angle (theta) with respect to global angle (phi).

  if (present(perp)) then
    beta = sqrt( (dr0(2)**2 + dr0(3)**2) / (dr0(1)**2 + dr0(3)**2) )
    dtheta_dphi = (beta**2 * sin_theta**2 + cos_theta**2) / beta
    r = dtheta_dphi * (dr2_dtheta * [cos_theta, sin_theta, 0.0_rp] + r2_wall * [-sin_theta, cos_theta, 0.0_rp])
    floor2_dw = coords_relative_to_floor (floor_sec2, r)
    floor2_dw%r = floor2_dw%r - floor_sec2%r
  endif
  
  ! Interpolate to get r0 which is on the line between the section origins
  ! and r_p which is on the line between floor1_p and floor2_p. 
  ! Note: If there is no spline then r_p is the particle position.

  s_rel = s1 / (s1 + s2)
  p1 = 1 - s_rel + sec1%p1_coef(1)*s_rel + sec1%p1_coef(2)*s_rel**2 + sec1%p1_coef(3)*s_rel**3
  p2 =     s_rel + sec1%p2_coef(1)*s_rel + sec1%p2_coef(2)*s_rel**2 + sec1%p2_coef(3)*s_rel**3

  r0  = p1 * floor_sec1%r + p2 * floor_sec2%r
  r_p = p1 * floor1_p%r + p2 * floor2_p%r

  ! Calculate rw which is the point on the wall that intersects the line through r0 & r_p

  f = norm2(cross_product(r_p - r0, floor2_w%r - floor1_w%r))
  if (f == 0) then  ! At origin so give something approximate
    d_radius = -norm2(p1 * floor1_w%r + p2 * floor2_w%r - r0)
    if (present(radius_wall)) radius_wall = -d_radius
  else
    alpha = norm2(cross_product(floor1_w%r - r0, floor2_w%r - floor1_w%r)) / f
    rw = r0 + alpha * (r_p - r0)
    d_radius = norm2(r_p - r0) - norm2(rw - r0)
    if (present(radius_wall)) radius_wall = norm2(rw - r0)
  endif

  if (present(origin)) then
    floor%r = r0
    floor = coords_floor_to_local_curvilinear (floor, ele, status, relative_to = upstream_end$)
    origin = floor%r
    if (ele%key /= patch$) origin(3) = origin(3) + ele%value(l$)
  endif

  ! Calculate the surface normal vector

  if (present (perp)) then
    p1 = norm2(rw - floor2_w%r)
    p2 = norm2(rw - floor1_w%r)
    drw = p1 * floor1_dw%r + p2 * floor2_dw%r
    floor%r = cross_product(drw, floor2_w%r - floor1_w%r)
    floor = coords_floor_to_relative (ele%floor, floor, .false., .true.)  ! To patch coords
    perp = floor%r / norm2(floor%r)  ! Normalize vector length to 1.    
  endif

!----------------------------
! non-patch region

else

  if (ele%key == patch$ .and. s_particle /= ele%s_start .and. s_particle /= ele%s) then
    call out_io (s_fatal$, r_name, &
          'WALL3D RADIUS CALCULATION FAILURE IN PATCH ELEMENT: ' // trim(ele%name) // '  (# \i0\)', &
          'THE PROBLEM GENERALLY IS A WALL SECTION TOO NEAR A PATCH ELEMENT.', &
          'SEE THE BMAD MANUAL FOR MORE DETAILS', &
          'THE SOLUTION GENERALLY IS TO MOVE THE WALL SECTION AWAY FROM THE PATH ELEMENT', &
          i_array = [ele%ix_ele])
    if (global_com%exit_on_error) call err_exit
    return
  endif

  s_rel = (s_particle - s1) / ds

  x0 = (1 - s_rel) * sec1%r0(1) + s_rel * sec2%r0(1)
  y0 = (1 - s_rel) * sec1%r0(2) + s_rel * sec2%r0(2)
  x = position(1) - x0; y = position(3) - y0
  r_particle = sqrt(x**2 + y**2)
  if (r_particle /= 0) then
    cos_theta = x / r_particle
    sin_theta = y / r_particle
  endif

  call calc_wall_radius (sec1%v, cos_theta, sin_theta, r1_wall, dr1_dtheta, ix_vertex1)
  call calc_wall_radius (sec2%v, cos_theta, sin_theta, r2_wall, dr2_dtheta, ix_vertex2)

  ! Interpolate to get d_radius

  p1 = 1 - s_rel + sec1%p1_coef(1)*s_rel + sec1%p1_coef(2)*s_rel**2 + sec1%p1_coef(3)*s_rel**3
  p2 =     s_rel + sec1%p2_coef(1)*s_rel + sec1%p2_coef(2)*s_rel**2 + sec1%p2_coef(3)*s_rel**3

  r_wall = p1 * r1_wall + p2 * r2_wall
  d_radius = r_particle - r_wall
  if (present(radius_wall)) radius_wall = r_wall

  ! Calculate the surface normal vector

  if (present(origin)) origin = [x0, y0, position(5)]

  if (present (perp)) then
    perp(1:2) = [cos_theta, sin_theta] - [-sin_theta, cos_theta] * &
                          (p1 * dr1_dtheta + p2 * dr2_dtheta) / r_wall
    dp1 = -1 + sec1%p1_coef(1) + 2 * sec1%p1_coef(2)*s_rel + 3 * sec1%p1_coef(3)*s_rel**2
    dp2 =  1 + sec1%p2_coef(1) + 2 * sec1%p2_coef(2)*s_rel + 3 * sec1%p2_coef(3)*s_rel**2
    perp(3)   = -(dp1 * r1_wall + dp2 * r2_wall) / ds

    ! In a bend dw_perp must be corrected since the true longitudinal "length" at the particle
    ! is, for a horizontal bend, ds * (1 + x/rho) where ds is the length along the reference 
    ! trajectory, x is the transverse displacement, and rho is the bend radius.

    if (ele%key == sbend$ .or. ele%key == rf_bend$) then
      if (ele%value(ref_tilt_tot$) == 0) then
        disp = position(1) 
      else
        disp = position(1) * cos(ele%value(ref_tilt_tot$)) + position(3) * sin(ele%value(ref_tilt_tot$))
      endif
      perp(3) = perp(3) / (1 + disp * ele%value(g$))
    endif
    ! Normalize vector length to 1.
    perp = perp / norm2(perp)  
    ! If section origin line is not aligned with the z-axis then the wall has a "shear"
    ! and the perpendicular vector must be corrected.
    dx = sec2%r0(1) - sec1%r0(1)
    dy = sec2%r0(2) - sec1%r0(2)
    if (dx /= 0 .or. dy /= 0) then
      perp(3) = perp(3) - (perp(1) * dx + perp(2) * dy) / ds
      perp = perp / norm2(perp)
    endif
  endif

endif

if (present(err_flag)) err_flag = .false.

!---------------------------------------------------------------------------
contains

subroutine d_radius_at_section (this_sec, d_rad)

type (wall3d_section_struct) this_sec
real(rp) r_wall, dr_dtheta
real(rp), optional :: d_rad
integer ixv

!

x = position(1) - this_sec%r0(1); y = position(3) - this_sec%r0(2)
r_particle = sqrt(x**2 + y**2)
if (r_particle == 0) then
  cos_theta = 1
  sin_theta = 0
else
  cos_theta = x / r_particle
  sin_theta = y / r_particle
endif

call calc_wall_radius (this_sec%v, cos_theta, sin_theta, r_wall, dr_dtheta, ixv)
d_radius = r_particle - r_wall
if (present(d_rad)) d_rad = d_radius

if (present(perp)) perp = [cos_theta, sin_theta, 0.0_rp] - &
                          [-sin_theta, cos_theta, 0.0_rp] * dr_dtheta / r_wall
if (present(origin)) origin = [this_sec%r0, position(5)]
if (present(err_flag)) err_flag = .false.
if (present(radius_wall)) radius_wall = r_wall

end subroutine d_radius_at_section

!---------------------------------------------------------------------------
! contains

function this_coords_to_floor (ele, r) result (floor)

type (ele_struct) ele
type (floor_position_struct) floor, local

real(rp) r(3)

!

local%r = r
call mat_make_unit(local%w)
floor = coords_local_curvilinear_to_floor(local, ele, relative_to = upstream_end$)

end function this_coords_to_floor

end function wall3d_d_radius

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_wall3d (ele, ix_wall, ds_offset, is_branch_wall) result (wall3d)
!
! Function to return a pointer to a wall3d structure associated
! with a given lattice element. 
!
! Note: The wall associated with a the vacuum chamber is the branch%wall3d.
!
! Input:
!   ele            -- ele_struct: lattice element.
!   ix_wall        -- integer, optional: index in wall3d(:) array. Default is 1.
!
! Output:
!   wall3d         -- wall3d_struct, pointer: Pointer to the associated wall structure.
!                       Will be nullified if there is no associated wall.
!   ds_offset      -- real(rp), optional: Element offset: s(beginning of ele) - s(beginning of wall3d)
!   is_branch_wall -- logical, optional: Set True if wall3d points to branch%wall3d.
!-

function pointer_to_wall3d (ele, ix_wall, ds_offset, is_branch_wall) result (wall3d)

character(32), parameter :: r_name = 'pointer_to_wall3d'

type (ele_struct), target :: ele
type (wall3d_struct), pointer :: wall3d
type (branch_struct), pointer :: branch

integer, optional :: ix_wall
integer iw

real(rp), optional :: ds_offset
logical, optional :: is_branch_wall

! 

nullify(wall3d)

iw = integer_option(1, ix_wall)

branch => pointer_to_branch(ele)
if (ele%key /= capillary$ .and. ele%key /= diffraction_plate$ .and. &
                                ele%key /= mask$ .and. associated (branch)) then
  if (.not. associated(branch%wall3d)) return
  wall3d => branch%wall3d(iw)
  if (present(ds_offset)) ds_offset = ele%s_start - branch%ele(0)%s
  if (present(is_branch_wall)) is_branch_wall = .true.
  return
endif

if (present(is_branch_wall)) is_branch_wall = .false.

if (.not. associated(ele%wall3d)) return
wall3d => ele%wall3d(iw)

if (present(ds_offset)) then
  select case (wall3d%ele_anchor_pt)
  case (anchor_beginning$); ds_offset = -ele%value(l$)
  case (anchor_center$);    ds_offset = -ele%value(l$) / 2
  case (anchor_end$);       ds_offset = 0 
  end select
endif

end function pointer_to_wall3d

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function wall3d_to_position (orbit, ele) result (position)
!
! Routine to return the suitable postion to be used in calling wall3d_d_radius
!
! This routine assumes that if in a patch the coordinates of orbit are with respect 
! to the downstream end if orbit%direction*orbit%time_dir = 1 and vice versa.
!
! Input:
!   orbit       -- coord_struct: Particle position.
!   ele         -- ele_struct: Element particle is in.
!
! Output:
!   position(6) -- real(rp): Position used in wall3d_d_radius call.
!-

function wall3d_to_position (orbit, ele) result (position)

implicit none

type (coord_struct) orbit
type (ele_struct) ele

real(rp) position(6), r_vec(3), ww(3,3)

! position(2) and position(4) are immaterial.


if (ele%key == patch$) then 
  ! Must transform to entrance coordinates.
  if (orbit%direction*orbit%time_dir == 1) then
    call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), w_mat = ww)
    r_vec = [orbit%vec(1), orbit%vec(3), orbit%s - ele%s]
    position(1:5:2) = matmul(ww, r_vec) + [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)]
    position(2:6:2) = [0.0_rp, 0.0_rp, 1.0_rp]
  else
    position = [orbit%vec(1:4), orbit%s - ele%s_start, -1.0_rp]
  endif

else
  position = [orbit%vec(1:4), orbit%s-ele%s_start, 1.0_rp * orbit%direction*orbit%time_dir]
endif

end function wall3d_to_position

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine create_concatenated_wall3d (lat)
!
! Routine to concatinate lat%branch(i)ele(:)%wall3d%section(:) arrays into
! one lat%branch(i)%wall3d%section(:) array.
!
! Exceptions: capillary and aperture elements do not have their walls included.
!
! Module needed:
!   use wall3d_mod
!
! Input:
!   lat      -- lat_struct: lattice
!
! Output:
!   lat      -- lat_struct: Lattice
!   err_flag -- logical: Set True if there is an error, false otherwise.
!-

Subroutine create_concatenated_wall3d (lat, err)

type section_ptr_struct
  type (wall3d_section_struct), pointer :: sec
  type (ele_struct), pointer :: ele
  real(rp) s
end type

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, ele1, ele2
type (section_ptr_struct), allocatable :: sp(:)
type (wall3d_section_struct), pointer :: ws

real(rp) s_min, s_max, s_temp

integer i, j, k, n, n_sec
logical err

character(*), parameter :: r_name = 'create_concatenated_wall3d'

! Count number of sections. This may be an overcount if there is superimpose.

err = .false.

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)

  s_min = branch%ele(0)%s
  s_max = branch%ele(branch%n_ele_track)%s

  n_sec = 0
  ele => branch%ele(0)
  do 
    ele => ele%next_in_branch(i)
    if (.not. associated(ele)) exit
    if (.not. associated(ele%wall3d)) cycle
    if (ele%key == capillary$) cycle
    if (ele%key == diffraction_plate$) cycle
    if (ele%key == mask$) cycle
    if (ele%lord_status == multipass_lord$) cycle  ! wall info also in slaves
    n_sec = n_sec + size(ele%wall3d(1)%section)
  enddo

  if (n_sec == 0) then
    if (associated (branch%wall3d)) deallocate (branch%wall3d)
    cycle
  endif

  ! Aggragate vacuum chamber wall info for a branch to branch%wall3d structure
  ! First work on non-superimpose element

  if (allocated(sp)) deallocate (sp)
  allocate (sp(n_sec))

  n_sec = 0
  ele => branch%ele(0)
  do 
    ele => ele%next_in_branch(i)
    if (.not. associated(ele)) exit
    if (.not. associated(ele%wall3d)) cycle
    if (ele%key == capillary$) cycle
    if (ele%key == diffraction_plate$) cycle
    if (ele%key == mask$) cycle
    if (ele%wall3d(1)%superimpose) cycle
    if (ele%lord_status == multipass_lord$) cycle
    call add_in_ele_wall_sections (ele, ele) ; if (err) return
  enddo

  ! Add superposition sections

  ele => branch%ele(0)
  do
    ele => ele%next_in_branch(i)
    if (.not. associated(ele)) exit
    if (.not. associated(ele%wall3d)) cycle
    if (ele%key == capillary$) cycle
    if (ele%key == diffraction_plate$) cycle
    if (ele%key == mask$) cycle
    if (.not. ele%wall3d(1)%superimpose) cycle
    if (ele%lord_status == multipass_lord$) cycle
    call superimpose_this_wall (ele, ele) ; if (err) return
  enddo

  ! Check for consistancy
  ! If there is an overlap but within significant_length then switch s-positions

  do j = 1, n_sec-1
    if (sp(j)%s > sp(j+1)%s) then
      if (sp(j)%s < sp(j+1)%s + bmad_com%significant_length) then
        s_temp = sp(j)%s
        sp(j)%s = sp(j+1)%s
        sp(j+1)%s = s_temp
      else
        call out_io (s_error$, r_name, 'WALL SECTIONS LONGITUDINALLY OUT-OF-ORDER', &
                     'SECTION AT: \es20.8\ FROM ELEMENT: ' // trim(sp(j)%ele%name) // ' (\i0\)', &
                     'NEXT SECTION AT: \es20.8\ FROM ELEMENT: ' // trim(sp(j+1)%ele%name) // ' (\i0\)', &
                     i_array = [sp(j)%ele%ix_ele, sp(j+1)%ele%ix_ele], r_array = [sp(j)%s, sp(j+1)%s])
        err = .true.
        return
      endif
    endif
  enddo

  ! Transfer info from sp to branch%wall3d
  ! branch%wall3d is never mutiply linked.

  if (.not. associated(branch%wall3d)) allocate (branch%wall3d(1))
  call re_allocate(branch%wall3d(1)%section, n_sec)

  do j = 1, n_sec
    ws => branch%wall3d(1)%section(j)
    call re_allocate(ws%v, size(sp(j)%sec%v))
    ws = sp(j)%sec
    ws%s = sp(j)%s
    ws%ix_ele = sp(j)%ele%ix_ele
    ws%ix_branch = sp(j)%ele%ix_branch
  enddo

  ! Mark patch regions

  call mark_patch_regions(branch)

enddo

!-----------------------------------------------------------------------------------------------
contains

subroutine add_in_ele_wall_sections (wall_ele, fiducial_ele)

type (ele_struct), target :: wall_ele, fiducial_ele
type (wall3d_struct), pointer :: wall
real(rp) s_ref, s
integer ii, k, ixw, nw, n, ix_wrap1, ix_wrap2

!

wall => wall_ele%wall3d(1)
nw = size(wall%section)

select case (wall%ele_anchor_pt)
case (anchor_beginning$); s_ref = fiducial_ele%s_start
case (anchor_center$);    s_ref = fiducial_ele%s_start + fiducial_ele%value(l$) / 2
case (anchor_end$);       s_ref = fiducial_ele%s 
end select

! If the element wall has more than one section (so the wall has a finite length), add
! significant_length/10 to s to avoid a roundoff bug.

s = wall%section(1)%s + s_ref
if (size(wall%section) /= 1) s = s + bmad_com%significant_length/10
ixw = bracket_index (s, sp(1:n_sec)%s, 1)

if (ixw > 1 .and. ixw < n_sec) then
  if (sp(ixw-1)%ele%ix_ele == sp(ixw+1)%ele%ix_ele) then
    call print_overlap_error (section_ptr_struct(wall%section(1), ele, s), sp(ixw+1))
    return
  endif
endif

! Move existing sections if needed to make room for the sections of wall_ele.

if (ixw < n_sec) then
  sp(ixw+1+nw:n_sec+nw) = sp(ixw+1:n_sec)
endif

ix_wrap1 = 0; ix_wrap2 = 0
do ii = 1, nw
  k = ii + ixw
  sp(k)%sec => wall%section(ii)
  sp(k)%s = wall%section(ii)%s + s_ref
  sp(k)%ele => wall_ele

  if (sp(k)%s < s_min)                     ix_wrap1 = k
  if (sp(k)%s > s_max .and. ix_wrap2 == 0) ix_wrap2 = k
enddo

n_sec = n_sec + nw

n = nw+ixw

! If there is an overlap but within significant_length then switch s-positions.

if (n < n_sec) then
  if (sp(n)%s > sp(n+1)%s) then
    if (sp(n)%s < sp(n+1)%s + bmad_com%significant_length) then
      s = sp(n)%s
      sp(n)%s = sp(n+1)%s
      sp(n+1)%s = s
    else
      call print_overlap_error (sp(n), sp(n+1))
      return
    endif
  endif
endif

! Wrap sections if needed

if (ix_wrap1 /= 0 .and. branch%param%geometry == closed$) then
  sp(1:ix_wrap1)%s = sp(1:ix_wrap1)%s + (s_max - s_min)
  sp(1:n_sec) = [sp(ix_wrap1+1:n_sec), sp(1:ix_wrap1)]
endif

if (ix_wrap2 /= 0 .and. branch%param%geometry == closed$) then
  sp(ix_wrap2:n_sec)%s = sp(ix_wrap2:n_sec)%s - (s_max - s_min)
  sp(1:n_sec) = [sp(ix_wrap2:n_sec), sp(1:ix_wrap2-1)]
endif

end subroutine add_in_ele_wall_sections

!-----------------------------------------------------------------------------------------------
! contains

subroutine print_overlap_error (sp1, sp2)

type (section_ptr_struct) sp1, sp2

!

call out_io (s_error$, r_name, 'WALLS OVERLAP LONGITUDINALLY BETWEEN', &
           'ELEMENT: ' // trim(sp1%ele%name) // ' (\i0\) Section S = \f14.6\ ', &
           'AND ELEMENT: ' // trim(sp2%ele%name) // ' (\i0\) Section S = \f14.6\ ', &
           i_array = [sp1%ele%ix_ele, sp2%ele%ix_ele], r_array = [sp1%s, sp2%s])
err = .true.

end subroutine print_overlap_error

!-----------------------------------------------------------------------------------------------
! contains

subroutine superimpose_this_wall (wall_ele, fiducial_ele)

type (ele_struct), target :: wall_ele, fiducial_ele
type (wall3d_struct), pointer :: wall
real(rp) s_ref, s
integer ii, ixw1, ixw2, nw, n_del, ix_wrap1, ix_wrap2

!

wall => wall_ele%wall3d(1)
nw = size(wall%section)

select case (wall%ele_anchor_pt)
case (anchor_beginning$); s_ref = fiducial_ele%s_start
case (anchor_center$);    s_ref = fiducial_ele%s_start + fiducial_ele%value(l$) / 2
case (anchor_end$);       s_ref = fiducial_ele%s 
end select

! If the element wall has more than one section (so the wall has a finite length), add
! significant_length/10 to s to avoid a roundoff bug.

s = wall%section(1)%s + s_ref
if (size(wall%section) /= 1) s = s + bmad_com%significant_length/10
ixw1 = bracket_index (s, sp(1:n_sec)%s, 1)

s = wall%section(nw)%s + s_ref
if (size(wall%section) /= 1) s = s - bmad_com%significant_length/10
ixw2 = bracket_index (s, sp(1:n_sec)%s, 1)

! If ixw2 < ixw1 then basically s(1) = s(nw) and there is a wall section (or sections) from a previous element at this s. 

if (ixw2 < ixw1) then
  ii = ixw1
  ixw1 = ixw2
  ixw2 = ii
endif

!

n_del = nw - (ixw2 - ixw1)  ! net number of sections added.

if (ixw2 < n_sec) then
  sp(ixw2+1+n_del:n_sec+n_del) = sp(ixw2+1:n_sec)
endif

ix_wrap1 = 0; ix_wrap2 = 0

do ii = 1, nw
  k = ii + ixw1
  sp(k)%sec => wall%section(ii)
  sp(k)%s = wall%section(ii)%s + s_ref
  sp(k)%ele => wall_ele

  if (sp(k)%s < s_min)                     ix_wrap1 = k
  if (sp(k)%s > s_max .and. ix_wrap2 == 0) ix_wrap2 = k
enddo

n_sec = n_sec + n_del

! Wrap sections if needed.
! Remember to discard any sections in the overlap region.

if (ix_wrap1 /= 0 .and. branch%param%geometry == closed$) then
  sp(1:ix_wrap1)%s = sp(1:ix_wrap1)%s + (s_max - s_min)
  do ii = ix_wrap1+1, n_sec
    if (sp(ii)%s <= sp(1)%s) cycle
    n_sec = ii - 1
    exit
  enddo    
  sp(1:n_sec) = [sp(ix_wrap1+1:n_sec), sp(1:ix_wrap1)]
endif

if (ix_wrap2 /= 0 .and. branch%param%geometry == closed$) then
  sp(ix_wrap2:n_sec)%s = sp(ix_wrap2:n_sec)%s - (s_max - s_min)
  do ii = ix_wrap2-1, 1, -1
    if (sp(ii)%s >= sp(n_sec)%s) cycle
    sp(1:n_sec-ii) = sp(ii+1:n_sec)
    n_sec = n_sec - ii
    ix_wrap2 = ix_wrap2 - ii
    exit
  enddo
  sp(1:n_sec) = [sp(ix_wrap2:n_sec), sp(1:ix_wrap2-1)]
endif

end subroutine superimpose_this_wall 

end subroutine create_concatenated_wall3d

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mark_patch_regions (branch)
!
! Routine to mark which regions in a wall3d structure contain patch elements.
! This routine should be called by any routine that creates a beam chamber wall.
!
! Input:
!   branch -- branch_struct: Lattice branch with %wall3d beam chamber wall.
!
! Output:
!   branch -- branch_struct: Lattice branch with %wall3d%section(i)%patch_in_region marked.
!-

subroutine mark_patch_regions (branch)

type (branch_struct), target :: branch
type (wall3d_struct), pointer :: wall
type (ele_struct), pointer :: ele
integer i, j, iw, n_sec, n_patch, ix1, ix2
integer, allocatable :: patch_ixs(:)
real(rp) s_min, s_max

!

allocate (patch_ixs(branch%n_ele_track))
n_patch = 0
do i = 1, branch%n_ele_track
  if (branch%ele(i)%key == patch$) then
    n_patch = n_patch + 1
    patch_ixs(n_patch) = i
  endif
enddo

do iw = 1, size(branch%wall3d)

  wall => branch%wall3d(iw)
  wall%section%patch_in_region = .false.

  if (n_patch == 0) cycle

  n_sec = size(wall%section)

  do j = 1, n_patch
    ! For each patch element:
    ! Determine its [s_min, s_max]
    ele => branch%ele(patch_ixs(j))
    s_min = min(ele%s_start, ele%s)
    s_max = max(ele%s_start, ele%s)

    ! Find the wall points (defined cross-sections) to either side of the patch.
    ix1 = bracket_index(s_min, wall%section%s, 1)
    ix2 = bracket_index(s_max, wall%section%s, 1)

    ! For open/closed lattices,
    ! if the wall section `s` falls in [s_min, s_max], the patch is in the region
    do i = max(2, ix1), min(n_sec, ix2+2)
      if (wall%section(i-1)%s < s_max .and. wall%section(i)%s > s_min) then
        wall%section(i)%patch_in_region = .true.
      endif
    enddo

    ! For closed geometries - the final wall section is a wrap-around special:
    !  either the first section s > s_min
    !   or
    !  the last section         s < s_max
    ! In that wraparound scenario, the first section type won't be wall_start.
    if (branch%param%geometry == closed$ .and. wall%section(1)%type /= wall_start$) then
      if (wall%section(n_sec)%s < s_max .or. wall%section(1)%s > s_min) then
        wall%section(1)%patch_in_region = .true.
      endif
    endif
  enddo

enddo

end subroutine mark_patch_regions

end module wall3d_mod
