module wall3d_mod

use re_allocate_mod
use output_mod

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
  real(rp) :: s_spline(3) = [1, 0, 0]   ! Longitudinal spline coefs. 
  integer :: n_slice_spline = 1         ! Number of slices used for the spline.
  type (wall3d_vertex_struct), allocatable :: v(:) 
                                        ! Array of vertices
  integer n_vertex_input                ! Number of vertices specified by the user.
end type

type wall3d_struct
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
!   criss(:) -- Wall3d_section_struct, pointer: Associated array.
!-

subroutine re_associate_section_array (section, n, exact)

implicit none

type (wall3d_section_struct), pointer :: section(:), temp_section(:)

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact

!

if (associated(section)) then
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
! Subroutine wall3d_section_initializer (section, err)
!
! Routine to initialize a wall3d_section_struct:
!   1) Add vertex points if there is symmetry.
!   2) Compute circular and elliptical centers.
!   3) Calc s_spline(3)
!
! Modules needed:
!   use wall3d_section_mod
!
! Input:
!   section  -- Wall3d_section_struct: Wall3d section.
!   
! Output:
!   section  -- Wall3d_section_struct: Initialized section-section.
!   err    -- Logical: Set true if there is a problem. EG: Convex section, etc.
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

! Calculate s_spline(3)

section%s_spline(3) = 1 - section%s_spline(1) - section%s_spline(2)

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

end module
