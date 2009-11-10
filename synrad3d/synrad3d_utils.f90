module synrad3d_utils

use synrad3d_struct
use random_mod
use photon_init_mod

private wall_pt_params

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)
!
! Routine to get the parameters at a photon emission point.
!
! Modules needed:
!   use synrad3d_utils
!
! Input:
!   lat       -- lat_struct with twiss propagated and mat6s made
!   ix_ele    -- integer: index of lat element to start ray from
!   s_offset  -- real(rp): offset along the length of the element 
!                         to use as a starting point for ray
!   orb(0:*)  -- coord_struct: orbit of particles to use as 
!                             source of ray
!   direction -- integer: +1 In the direction of increasing s.
!                         -1 In the direction of decreasing s.
!
! Output:
!   photon    -- photon3d_coord_struct: Generated photon.
!-

subroutine get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)

use em_field_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct) :: orb(0:), orb_here, orb1
type (ele_struct), pointer :: ele0, ele
type (ele_struct) :: ele_here
type (photon3d_coord_struct) :: photon
type (em_field_struct) :: field

real(rp) s_offset, k_wig, g_max, l_small, gx, gy

integer direction, ix_ele

! initialize the photon

ele0 => lat%ele(ix_ele-1)
ele => lat%ele(ix_ele)


! set the photon's initial twiss values

call twiss_and_track_partial (ele0, ele, lat%param, s_offset, ele_here, & 
                                                    orb(ix_ele-1), orb_here)

! set the photon's g_bend value (inverse bending radius at src pt) 

if (ele%key == sbend$) then  

  ! sbends are easy
  gx = 1 / ele%value(rho$)
  gy = 0
  if (ele%value(roll$) /= 0) then
    gy = gx * sin(ele%value(roll$))
    gx = gx * cos(ele%value(roll$))
  endif

elseif (ele%key == quadrupole$ .or. ele%key == sol_quad$) then

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

elseif (ele%key == wiggler$ .and. ele%sub_key == periodic_type$) then

  ! for periodic wigglers, get the max g_bend from 
  ! the max B field of the wiggler, then scale it 
  ! by the cos of the position along the poles

  k_wig = twopi * ele%value(n_pole$) / (2 * ele%value(l$))
  g_max = c_light * ele%value(b_max$) / (ele%value(p0c$))
  gx = g_max * cos (k_wig * s_offset)
  orb_here%vec(2) = orb_here%vec(2) + (g_max / k_wig) * sin (k_wig * s_offset)

elseif (ele%key == wiggler$ .and. ele%sub_key == map_type$) then

  ! for mapped wigglers, find the B field at the source point
  ! Note: assumes particles are relativistic!!

  call em_field_calc (ele_here, lat%param, s_offset, orb_here, .false., field)
  gx = field%b(2) * c_light / ele%value(p0c$)
  gy = field%b(1) * c_light / ele%value(p0c$)

else

  print *, 'ERROR: UNKNOWN ELEMENT HERE ', ele%name

endif

end subroutine


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_e, photon_direction, photon)
!
! subroutine to initialize a new photon
!
! Modules needed:
!   use synrad3d_utils
!
! Input:
!   lat       -- lat_struct with twiss propagated and mat6s made
!   ix_ele    -- integer: index of lat element to start ray from
!   s_offset  -- real(rp): offset along the length of the element 
!                         to use as a starting point for ray
!   orb(0:*)  -- coord_struct: orbit of particles to use as 
!                             source of ray
!   direction -- integer: +1 In the direction of increasing s.
!                         -1 In the direction of decreasing s.
!
! Output:
!   photon    -- photon_coord_struct: Generated photon.
!-

subroutine emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_e, photon_direction, p_orb)

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

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine photon_radius (p_orb, wall, radius, dw_perp)
!
! Routine to calculate the normalized transverse position of the photon 
! relative to the wall: 
!     radius = 0 => Center of the beam pipe 
!     radius = 1 => at wall.
!     radius > 1 => Outside the beam pipe.
!
! Modules needed:
!   use photon_utils
!
! Input:
!   wall -- wall3d_struct: Wall
!   s    -- Real(rp): Longitudinal position.
!
! Output:
!   radius       -- real(rp): Radius of beam relative to the wall.
!   dw_perp(3)   -- real(rp), optional: Outward normal vector perpendicular to the wall.
!-

Subroutine photon_radius (p_orb, wall, radius, dw_perp)

implicit none

type (wall3d_struct), target :: wall
type (wall3d_pt_struct), pointer :: wall0, wall1
type (photon3d_coord_struct), target :: p_orb

real(rp), optional :: dw_perp(:)
real(rp) r0, r1, f, radius
real(rp) dw_x0, dw_y0, dw_z0, dw_x1, dw_y1, dw_z1, dw_x, dw_y, dw_z, denom
real(rp), pointer :: vec(:)

integer ix

!

vec => p_orb%vec
call bracket_index (wall%pt%s, 0, wall%n_pt_max, vec(5), ix)
if (ix == wall%n_pt_max) ix = wall%n_pt_max - 1

! The outward normal vector is discontinuous at the wall points.
! If at a wall point, use the correct part of the wall.

if (vec(5) == wall%pt(ix)%s .and. vec(6) > 0) then
  if (ix /= 0) then
    ix = ix - 1
  endif
endif

!

call wall_pt_params (wall%pt(ix),   vec, r0, dw_x0, dw_y0, dw_z0)
call wall_pt_params (wall%pt(ix+1), vec, r1, dw_x1, dw_y1, dw_z1)

f = (vec(5) - wall%pt(ix)%s) / (wall%pt(ix+1)%s - wall%pt(ix)%s)
radius = (1 - f) * r0 + f * r1

if (present (dw_perp)) then
  dw_x = (1 - f) * dw_x0 + f * dw_x1
  dw_y = (1 - f) * dw_y0 + f * dw_y1
  dw_z = (dw_z1 - dw_z0) / (wall%pt(ix+1)%s - wall%pt(ix)%s)
  denom = sqrt(dw_x**2 + dw_y**2 + dw_z**2)
  dw_perp = [dw_x, dw_y, dw_z] / denom
endif

end subroutine photon_radius

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine wall_pt_params (wall_pt, vec, r, dw_x, dw_y, dw_z)
!
! Private routine used by photon_radius
!-

subroutine wall_pt_params (wall_pt, vec, r, dw_x, dw_y, dw_z)

implicit none

type (wall3d_pt_struct) wall_pt
real(rp) r, dw_x, dw_y, dw_z, vec(6)

!

if (wall_pt%type == 'rectangular') then
  if (abs(vec(1)/wall_pt%width2) > abs(vec(3)/wall_pt%height2)) then
    dw_x = sign(1.0_rp, vec(1)) / wall_pt%width2
    dw_y = 0
    dw_z = abs(vec(1)/wall_pt%width2)
    r = dw_z
  else
    dw_x = 0
    dw_y = sign(1.0_rp, vec(3)) / wall_pt%height2
    dw_z = abs(vec(3)/wall_pt%height2)
    r = dw_z
  endif

elseif (wall_pt%type == 'elliptical') then
  dw_z = sqrt((vec(1)/wall_pt%width2)**2 + (vec(3)/wall_pt%height2)**2)
  r = dw_z
  if (dw_z == 0) then
    dw_x = 0
    dw_y = 0
  else
    dw_x = vec(1) / wall_pt%width2**2 / dw_z
    dw_y = vec(3) / wall_pt%height2**2 / dw_z
  endif

else
  print *, 'BAD WALL_PT%TYPE: ' // wall_pt%type
  call err_exit
endif

end subroutine

end module
