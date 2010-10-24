module synrad3d_utils

use synrad3d_struct
use random_mod
use photon_init_mod
use capillary_mod

private sr3d_wall_pt_params

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_check_wall (wall)
!
! Routine to check the vacuum chamber wall for problematic values.
!
! Input:
!   wall -- wall3d_struct: wall structure.
!-

subroutine sr3d_check_wall (wall)

implicit none

type (wall3d_struct), target :: wall
type (wall3d_pt_struct), pointer :: pt, pt0

integer i, ig0, ig1

character(20) :: r_name = 'sr3d_check_wall'

!

do i = 0, wall%n_pt_max
  pt => wall%pt(i)

  if (i > 0) then
    if (pt%s <= wall%pt(i-1)%s) then
      call out_io (s_fatal$, r_name, &
                'WALL%PT(i)%S: \es12.2\ ', &
                '    IS LESS THAN PT(i-1)%S: \es12.2\ ', &
                '    FOR I = \i0\ ', &
                r_array = [pt%s, wall%pt(i-1)%s], i_array = [i])
      call err_exit
    endif
  endif

  if (.not. any(pt%basic_shape == ['elliptical    ', 'rectangular   ', 'gen_shape     ', 'gen_shape_mesh'])) then
    call out_io (s_fatal$, r_name, &
              'BAD WALL%PT(i)%BASIC_SHAPE: ' // pt%basic_shape, &
              '    FOR I = \i0\ ', i_array = [i])
    call err_exit
  endif

  ! Gen_shape and gen_shape_mesh checks

  if (pt%basic_shape == 'gen_shape' .or. pt%basic_shape == 'gen_shape_mesh') then
    if (pt%ix_gen_shape < 1) then
      call out_io (s_fatal$, r_name, &
              'BAD WALL%PT(I)%IX_GEN_SHAPE SECTION NUMBER \i0\ ', i_array = [i])
      call err_exit
    endif
    if (pt%basic_shape == 'gen_shape') cycle
    if (i == 0) cycle

    pt0 => wall%pt(i-1)
    if (pt0%basic_shape /= 'gen_shape' .and. pt0%basic_shape /= 'gen_shape_mesh') then
      call out_io (s_fatal$, r_name, &
              'BASIC_SHAPE FOR SECTION PRECEEDING "gen_shape_mesh" SECTION MUST BE ', &
              '"gen_shape" OR "gen_shape_mesh" SECTION NUMBER \i0\ ', i_array = [i])
      call err_exit
    endif

    ig0 = pt0%ix_gen_shape
    ig1 = pt%ix_gen_shape

    if (size(wall%gen_shape(ig0)%v) /= size(wall%gen_shape(ig1)%v)) then
      call out_io (s_fatal$, r_name, &
              '"gen_shape_mesh" CONSTRUCT MUST HAVE THE SAME NUMBER OF VERTEX POINTS ON', &
              'SUCCESIVE CROSS-SECTIONS  \2i0\ ', i_array = [i-1, i])
      call err_exit
    endif

    cycle
  endif



  ! Checks for everything else

  if (pt%width2 <= 0) then
    call out_io (s_fatal$, r_name, &
              'BAD WALL%PT(i)%WIDTH2: \es12.2\ ', &
              '    FOR I = \i0\ ', r_array = [pt%width2], i_array = [i])
    call err_exit
  endif

  if (pt%width2 <= 0) then
    call out_io (s_fatal$, r_name, &
              'BAD WALL%PT(i)%HEIGHT2: \es12.2\ ', &
              '    FOR I = \i0\ ', r_array = [pt%height2], i_array = [i])
    call err_exit
  endif

  ! +x side check

  if (pt%ante_height2_plus < 0 .and.pt%width2_plus > 0) then
    if (pt%width2_plus > pt%width2) then
      call out_io (s_fatal$, r_name, &
              'WITHOUT AN ANTECHAMBER: WALL%PT(i)%WIDTH2_PLUS \es12.2\ ', &
              '    MUST BE LESS THEN WIDTH2 \es12.2\ ', &
              '    FOR I = \i0\ ', &
              r_array = [pt%width2_plus, pt%width2], i_array = [i])
      call err_exit
    endif
  endif

  ! -x side check

  if (pt%ante_height2_minus < 0 .and. pt%width2_minus > 0) then
    if (pt%width2_minus > pt%width2) then
      call out_io (s_fatal$, r_name, &
              'WITHOUT AN ANTECHAMBER: WALL%PT(i)%WIDTH2_MINUS \es12.2\ ', &
              '    MUST BE LESS THEN WIDTH2 \es12.2\ ', &
              '    FOR I = \i0\ ', &
              r_array = [pt%width2_minus, pt%width2], i_array = [i])
      call err_exit
    endif
  endif

enddo

! computations

do i = 0, wall%n_pt_max
  pt => wall%pt(i)

  ! +x side computation

  if (pt%ante_height2_plus > 0) then
    if (pt%basic_shape == 'elliptical') then
      pt%ante_x0_plus = pt%width2 * sqrt (1 - (pt%ante_height2_plus / pt%height2)**2)
    else
      pt%ante_x0_plus = pt%width2
    endif

    if (pt%width2_plus <= pt%ante_x0_plus) then
      call out_io (s_fatal$, r_name, &
              'WITH AN ANTECHAMBER: WALL%PT(i)%WIDTH2_PLUS \es12.2\ ', &
              '    MUST BE GREATER THEN: \es12.2\ ', &
              '    FOR I = \i0\ ', &
              r_array = [pt%width2_plus, pt%ante_x0_plus], i_array = [i])
      call err_exit
    endif

  elseif (pt%width2_plus > 0) then
    if (pt%basic_shape == 'elliptical') then
      pt%y0_plus = pt%height2 * sqrt (1 - (pt%width2_plus / pt%width2)**2)
    else
      pt%y0_plus = pt%height2
    endif
  endif

  ! -x side computation

  if (pt%ante_height2_minus > 0) then
    if (pt%basic_shape == 'elliptical') then
      pt%ante_x0_minus = pt%width2 * sqrt (1 - (pt%ante_height2_minus / pt%height2)**2)
    else
      pt%ante_x0_minus = pt%width2
    endif

    if (pt%width2_minus <= pt%ante_x0_minus) then
      call out_io (s_fatal$, r_name, &
              'WITH AN ANTECHAMBER: WALL%PT(i)%WIDTH2_MINUS \es12.2\ ', &
              '    MUST BE GREATER THEN: \es12.2\ ', &
              '    FOR I = \i0\ ', &
              r_array = [pt%width2_minus, pt%ante_x0_minus], i_array = [i])

      call err_exit
    endif

  elseif (pt%width2_minus > 0) then
    if (pt%basic_shape == 'elliptical') then
      pt%y0_minus = pt%height2 * sqrt (1 - (pt%width2_minus / pt%width2)**2)
    else
      pt%y0_minus = pt%height2
    endif
  endif

enddo

end subroutine sr3d_check_wall 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)
!
! Routine to get the parameters at a photon emission point.
!
! Modules needed:
!   use synrad3d_utils
!
! Input:
!   lat       -- lat_struct with twiss propagated and mat6s made.
!   orb(0:*)  -- coord_struct: orbit of particles to use as source of ray.
!   ix_ele    -- integer: index of lat element to start ray from.
!   s_offset  -- real(rp): Distance from beginning of element to the point where the photon is emitted.
!
! Output:
!   ele_here  -- ele_struct: Twiss parameters at emission point.
!   orb_here  -- coord_struct: Beam center coords at emission point.
!   gx        -- Real(rp): Horizontal 1/bending_radius.
!   gy        -- Real(rp): Vertical 1/bending_radius.
!-

subroutine sr3d_get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)

use em_field_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct) :: orb(0:), orb_here, orb1
type (ele_struct), pointer :: ele
type (ele_struct) ele_here
type (photon3d_coord_struct) :: photon
type (em_field_struct) :: field

real(rp) s_offset, k_wig, g_max, l_small, gx, gy
real(rp), save :: s_old_offset = 0

integer ix_ele

logical err
logical, save :: init_needed = .true.

! Init

if (init_needed) then
  call init_ele (ele_here)
  init_needed = .false.
endif

ele  => lat%ele(ix_ele)

! Calc the photon's initial twiss values.
! Tracking through a wiggler can take time so use twiss_and_track_intra_ele to
!   minimize the length over which we track.

if (ele_here%ix_ele /= ele%ix_ele .or. ele_here%ix_branch /= ele%ix_branch .or. s_old_offset > s_offset) then
  ele_here = lat%ele(ix_ele-1)
  ele_here%ix_ele = ele%ix_ele
  ele_here%ix_branch = ele%ix_branch
  orb_here = orb(ix_ele-1)
  s_old_offset = 0
endif

call twiss_and_track_intra_ele (ele, lat%param, s_old_offset, s_offset, .true., .true., &
                                            orb_here, orb_here, ele_here, ele_here, err)
if (err) call err_exit
s_old_offset = s_offset

! Calc the photon's g_bend value (inverse bending radius at src pt) 

select case (ele%key)
case (sbend$)  

  ! sbends are easy
  gx = 1 / ele%value(rho$)
  gy = 0
  if (ele%value(roll$) /= 0) then
    gy = gx * sin(ele%value(roll$))
    gx = gx * cos(ele%value(roll$))
  endif

case (quadrupole$, sol_quad$, elseparator$)

  ! for quads or sol_quads, get the bending radius
  ! from the change in x' and y' over a small 
  ! distance in the element

  l_small = 1e-2      ! something small
  ele_here%value(l$) = l_small
  call make_mat6 (ele_here, lat%param, orb_here, orb_here)
  call track1 (orb_here, ele_here, lat%param, orb1)
  orb1%vec = orb1%vec - orb_here%vec
  gx = orb1%vec(2) / l_small
  gy = orb1%vec(4) / l_small

case (wiggler$)

  if (ele%sub_key == periodic_type$) then

    ! for periodic wigglers, get the max g_bend from 
    ! the max B field of the wiggler, then scale it 
    ! by the cos of the position along the poles

    k_wig = twopi * ele%value(n_pole$) / (2 * ele%value(l$))
    g_max = c_light * ele%value(b_max$) / (ele%value(p0c$))
    gx = g_max * sin (k_wig * s_offset)
    gy = 0
    orb_here%vec(1) = (g_max / k_wig) * sin (k_wig * s_offset)
    orb_here%vec(2) = (g_max / k_wig) * cos (k_wig * s_offset)

  else

    ! for mapped wigglers, find the B field at the source point
    ! Note: assumes particles are relativistic!!

    call em_field_calc (ele_here, lat%param, ele_here%value(l$), orb_here, .false., field)
    gx = field%b(2) * c_light / ele%value(p0c$)
    gy = field%b(1) * c_light / ele%value(p0c$)

  endif

case default

  print *, 'ERROR: UNKNOWN ELEMENT HERE ', ele%name

end select

end subroutine sr3d_get_emission_pt_params 


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_e, photon_direction, photon)
!
! subroutine sr3d_to initialize a new photon
!
! Modules needed:
!   use synrad3d_utils
!
! Input:
!   ele_here  -- Ele_struct: Element emitting the photon. Emission is at the exit end of the element.
!   orb_here  -- coord_struct: orbit of particles emitting the photon.
!   gx, gy    -- Real(rp): Horizontal and vertical bending strengths.
!   emit_a    -- Real(rp): Emittance of the a-mode.
!   emit_b    -- Real(rp): Emittance of the b-mode.
!   photon_direction 
!             -- Integer: +1 In the direction of increasing s.
!                         -1 In the direction of decreasing s.
!
! Output:
!   photon    -- photon_coord_struct: Generated photon.
!-

subroutine sr3d_emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_e, photon_direction, p_orb)

implicit none

type (ele_struct), target :: ele_here
type (coord_struct) :: orb_here
type (photon3d_coord_struct) :: p_orb
type (twiss_struct), pointer :: t

real(rp) emit_a, emit_b, sig_e, gx, gy, g_tot, gamma
real(rp) orb(6), r(3), vec(4), v_mat(4,4)

integer photon_direction

! Get photon energy and "vertical angle".

g_tot = sqrt(gx**2 + gy**2)
call convert_total_energy_to (ele_here%value(E_tot$), electron$, gamma) 
call photon_init (g_tot, gamma, orb)
p_orb%energy = orb(6)
p_orb%vec = 0
p_orb%vec(4) = orb(4) / sqrt(orb(4)**2 + 1)

! rotate photon if gy is non-zero

if (gy /= 0) then
  p_orb%vec(2) = gy * p_orb%vec(4) / g_tot
  p_orb%vec(4) = gx * p_orb%vec(4) / g_tot
endif

! Offset due to finite beam size

call ran_gauss(r)
t => ele_here%a
vec(1:2) = (/ sqrt(t%beta*emit_a) * r(1)                    + t%eta  * sig_e * r(3), &
              sqrt(emit_a/t%beta) * (r(2) + t%alpha * r(1)) + t%etap * sig_e * r(3) /)

call ran_gauss(r)
t => ele_here%b
vec(3:4) = (/ sqrt(t%beta*emit_b) * r(1)                    + t%eta  * sig_e * r(3), &
              sqrt(emit_b/t%beta) * (r(2) + t%alpha * r(1)) + t%etap * sig_e * r(3) /)

call make_v_mats (ele_here, v_mat)

p_orb%vec(1:4) = p_orb%vec(1:4) + matmul(v_mat, vec)

! Offset due to non-zero orbit.

p_orb%vec(1:4) = p_orb%vec(1:4) + orb_here%vec(1:4)

! Longitudinal position

p_orb%vec(5) = ele_here%s

! Note: phase space coords here are different from the normal beam and photon coords.
! Here vec(2)^2 + vec(4)^2 + vec(6)^2 = 1

p_orb%vec(6) = photon_direction * sqrt(1 - orb(2)**2 - orb(4)**2)

end subroutine sr3d_emit_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_photon_d_radius (p_orb, wall, d_radius, dw_perp, in_antechamber)
!
! Routine to calculate the normalized transverse position of the photon 
! relative to the wall: 
!
! Modules needed:
!   use photon_utils
!
! Input:
!   wall -- wall3d_struct: Wall
!   s    -- Real(rp): Longitudinal position.
!
! Output:
!   radius       -- real(rp): r_photon - r_wall
!   dw_perp(3)   -- real(rp), optional: Outward normal vector perpendicular to the wall.
!   in_antechamber -- Logical, optional: At antechamber wall?
!-

Subroutine sr3d_photon_d_radius (p_orb, wall, d_radius, dw_perp, in_antechamber, ix_tri)

implicit none

type (wall3d_struct), target :: wall
type (photon3d_coord_struct), target :: p_orb

real(rp), optional :: d_radius, dw_perp(3)
real(rp) radius0, radius1, f, cos_ang, sin_ang, r_photon
real(rp) dr0_dtheta, dr1_dtheta, pt0(3), pt1(3), pt2(3), dp1(3), dp2(3)

integer, optional :: ix_tri
integer ix, ig0, ig1

logical, optional :: in_antechamber
logical in_ante0, in_ante1

!

call sr3d_get_wall_index (p_orb, wall, ix)

! gen_shape_mesh calc

if (wall%pt(ix+1)%basic_shape == 'gen_shape_mesh') then
  if (.not. present(dw_perp)) return
  ig0 = wall%pt(ix)%ix_gen_shape
  ig1 = wall%pt(ix+1)%ix_gen_shape
  call sr3d_get_mesh_wall_triangle_pts (wall%gen_shape(ig0), wall%gen_shape(ig1), p_orb%ix_triangle, pt0, pt1, pt2)
  dp1 = pt1 - pt0
  dp2 = pt2 - pt0
  dw_perp = [dp1(2)*dp2(3) - dp1(3)*dp2(2), dp1(3)*dp2(1) - dp1(1)*dp2(3), dp1(1)*dp2(2) - dp1(2)*dp2(1)]
  dw_perp = dw_perp / sqrt(sum(dw_perp**2))  ! Normalize
  return
endif

! Get the parameters at the defined cross-sections to either side of the photon position.

if (p_orb%vec(1) == 0 .and. p_orb%vec(3) == 0) then
  r_photon = 0
  cos_ang = 1
  sin_ang = 0
else
  r_photon = sqrt(p_orb%vec(1)**2 + p_orb%vec(3)**2)
  cos_ang = p_orb%vec(1) / r_photon
  sin_ang = p_orb%vec(3) / r_photon
endif

call sr3d_wall_pt_params (wall%pt(ix),   cos_ang, sin_ang, radius0, dr0_dtheta, in_ante0, wall)
call sr3d_wall_pt_params (wall%pt(ix+1), cos_ang, sin_ang, radius1, dr1_dtheta, in_ante1, wall)

f = (p_orb%vec(5) - wall%pt(ix)%s) / (wall%pt(ix+1)%s - wall%pt(ix)%s)

if (present(d_radius)) d_radius = r_photon - ((1 - f) * radius0 + f * radius1)

if (present (dw_perp)) then
  dw_perp(1:2) = [cos_ang, sin_ang] - [-sin_ang, cos_ang] * &
                            ((1 - f) * dr0_dtheta + f * dr1_dtheta) / r_photon
  dw_perp(3) = (radius0 - radius1) / (wall%pt(ix+1)%s - wall%pt(ix)%s)
  dw_perp = dw_perp / sqrt(sum(dw_perp**2))  ! Normalize
endif

if (present(in_antechamber)) in_antechamber = (in_ante0 .or. in_ante1)

end subroutine sr3d_photon_d_radius

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_get_wall_index (p_orb, wall, ix_wall)
!
! Routine to get the wall index such that 
! For p_orb%vec(6) > 0 (forward motion):
!   wall%pt(ix_wall)%s < p_orb%vec(5) <= wall%pt(ix_wall+1)%s
! For p_orb%vec(6) < 0 (backward motion):
!   wall%pt(ix_wall)%s <= p_orb%vec(5) < wall%pt(ix_wall+1)%s
! Exceptions:
!   If p_orb%vec(5) == wall%pt(0)%s (= 0)       -> ix_wall = 0
!   If p_orb%vec(5) == wall%pt(wall%n_pt_max)%s -> ix_wall = wall%n_pt_max - 1
!
! Input:
!   p_orb  -- photon3d_coord_struct: Photon position.
!   wall   -- wall3d_struct: Wall structure
!
! Output:
!   ix_wall -- Integer: Wall index
!-

subroutine sr3d_get_wall_index (p_orb, wall, ix_wall)

implicit none

type (photon3d_coord_struct) :: p_orb
type (wall3d_struct), target :: wall

integer ix_wall
integer, save :: ix_wall_old = 0

! 

ix_wall = ix_wall_old
if (p_orb%vec(5) < wall%pt(ix_wall)%s .or. p_orb%vec(5) > wall%pt(ix_wall+1)%s) then
  call bracket_index (wall%pt%s, 0, wall%n_pt_max, p_orb%vec(5), ix_wall)
  if (ix_wall == wall%n_pt_max) ix_wall = wall%n_pt_max - 1
endif

! vec(5) at boundary cases

if (p_orb%vec(5) == wall%pt(ix_wall)%s   .and. p_orb%vec(6) > 0 .and. ix_wall /= 0)               ix_wall = ix_wall - 1
if (p_orb%vec(5) == wall%pt(ix_wall+1)%s .and. p_orb%vec(6) < 0 .and. ix_wall /= wall%n_pt_max-1) ix_wall = ix_wall + 1

p_orb%ix_wall = ix_wall
ix_wall_old = ix_wall

end subroutine sr3d_get_wall_index

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_get_mesh_wall_triangle_pts (cross1, cross2, ix_tri, pt0, pt1, pt2)
!
! Routine to return the three vertex points for a triangular wall surface element between
! two cross-sections.
!
! Input:
!   cross1 -- cross_section_struct: A cross-section.
!   cross2 -- cross_section_struct: Second cross-section.
!   ix_tr  -- Integer: Triangle index. Must be between 1 and 2*size(cross1%v).
!               [Note: size(cross1%v) = size(cross2%v)]
!
! Output:
!   pt0(3), pt1(3), pt2(3)
!         -- Real(rp): (x, y, s) vertex points for the triangle.
!             Looking from the outside, the points are in counter-clockwise order.
!             This is important in determining the outward normal vector
!-

subroutine sr3d_get_mesh_wall_triangle_pts (cross1, cross2, ix_tri, pt0, pt1, pt2)

implicit none

type (cross_section_struct) cross1, cross2

integer ix_tri
integer ix1, ix2

real(rp) pt0(3), pt1(3), pt2(3)

! 

ix1 = (ix_tri + 1) / 2
ix2 = ix1 + 1
if (ix2 > size(cross1%v)) ix2 = 1

if (odd(ix_tri)) then
  pt0 = [cross1%v(ix1)%x, cross1%v(ix1)%y, cross1%s]
  pt1 = [cross2%v(ix1)%x, cross2%v(ix1)%y, cross2%s]
  pt2 = [cross1%v(ix2)%x, cross1%v(ix2)%y, cross1%s]
else
  pt1 = [cross1%v(ix2)%x, cross1%v(ix2)%y, cross1%s]
  pt2 = [cross2%v(ix1)%x, cross2%v(ix1)%y, cross2%s]
  pt0 = [cross2%v(ix2)%x, cross2%v(ix2)%y, cross2%s]
endif

end subroutine sr3d_get_mesh_wall_triangle_pts

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_wall_pt_params (wall_pt, cos_photon, sin_photon, r_wall, dr_dtheta, in_antechamber, wall)
!
! Routine to compute parameters needed by sr3d_photon_d_radius routine.
!
! Input:
!   wall_pt -- wall3d_pt_struct: Wall outline at a particular longitudinal location.
!   cos_photon -- Real(rp): Cosine of the photon transverse position.
!   sin_photon -- Real(rp): Sine of the photon transverse position.
!   wall       -- wall3d_struct: Needed to determine the basic_shape.
!
! Output:
!   r_wall         -- Real(rp): Radius of the wall
!   dr_dtheta      -- Real(rp): Transverse directional derivatives: d(r_wall)/d(theta)
!   in_antechamber -- Logical: Set true of particle is in antechamber
!-

subroutine sr3d_wall_pt_params (wall_pt, cos_photon, sin_photon, r_wall, dr_dtheta, in_antechamber, wall)

implicit none

type (wall3d_pt_struct) wall_pt, pt
type (wall3d_struct), target :: wall
type (cross_section_vertex_struct), pointer :: v(:)

real(rp) dr_dtheta, cos_photon, sin_photon 
real(rp) r_wall

integer ix

logical in_antechamber

! Init

in_antechamber = .false.

! general shape

if (wall_pt%basic_shape == 'gen_shape') then
  call calc_wall_radius (wall%gen_shape(wall_pt%ix_gen_shape)%v, cos_photon, sin_photon, r_wall, dr_dtheta)
  return
endif


! general shape: Should not be here

if (wall_pt%basic_shape == 'gen_shape_mesh') then
  call err_exit
endif

! Check for antechamber or beam stop...
! If the line extending from the origin through the photon intersects the
! antechamber or beam stop then pretend the chamber is rectangular with the 
! antechamber or beam stop dimensions.

! Positive x side check.

pt = wall_pt

if (cos_photon > 0) then

  ! If there is an antechamber...
  if (pt%ante_height2_plus > 0) then

    if (abs(sin_photon/cos_photon) < pt%ante_height2_plus/pt%ante_x0_plus) then  
      pt%basic_shape = 'rectangular'
      pt%width2 = pt%width2_plus
      pt%height2 = pt%ante_height2_plus
      if (cos_photon >= pt%ante_x0_plus) in_antechamber = .true.
    endif

  ! If there is a beam stop...
  elseif (pt%width2_plus > 0) then
    if (abs(sin_photon/cos_photon) < pt%y0_plus/pt%width2_plus) then 
      pt%basic_shape = 'rectangular'
      pt%width2 = pt%width2_plus
    endif

  endif

! Negative x side check

elseif (cos_photon < 0) then

  ! If there is an antechamber...
  if (pt%ante_height2_minus > 0) then

    if (abs(sin_photon/cos_photon) < pt%ante_height2_minus/pt%ante_x0_minus) then  
      pt%basic_shape = 'rectangular'
      pt%width2 = pt%width2_minus
      pt%height2 = pt%ante_height2_minus
      if (cos_photon >= pt%ante_x0_minus) in_antechamber = .true.
    endif

  ! If there is a beam stop...
  elseif (pt%width2_minus > 0) then
    if (abs(sin_photon / cos_photon) < pt%y0_minus/pt%width2_minus) then 
      pt%basic_shape = 'rectangular'
      pt%width2 = pt%width2_minus
    endif

  endif

endif

! Compute parameters

if (pt%basic_shape == 'rectangular') then
  if (abs(cos_photon/pt%width2) > abs(sin_photon/pt%height2)) then
    r_wall = pt%width2 / abs(cos_photon)
    dr_dtheta = r_wall * sin_photon / cos_photon
  else
    r_wall = pt%height2 / abs(sin_photon)
    dr_dtheta = -r_wall * cos_photon / sin_photon
  endif

elseif (pt%basic_shape == 'elliptical') then
  r_wall = 1 / sqrt((cos_photon/pt%width2)**2 + (sin_photon/pt%height2)**2)
  dr_dtheta = r_wall**3 * cos_photon * sin_photon * (1/pt%width2**2 - 1/pt%height2**2)
endif

end subroutine sr3d_wall_pt_params

end module
