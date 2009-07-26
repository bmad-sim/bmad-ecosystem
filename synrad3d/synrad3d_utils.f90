module synrad3d_utils

use synrad3d_struct

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
!   photon    -- photon_coord_struct: Generated photon.
!-

subroutine get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)

implicit none

type (lat_struct), target :: lat
type (coord_struct) :: orb(0:), orb_here, orb1
type (ele_struct), pointer :: ele0, ele
type (ele_struct) :: ele_here
type (photon_coord_struct) :: photon
type (em_field_struct) :: field

real(rp) l_offset, k_wig, g_max, l_small

integer direction, ix_ele

! initialize the photon

ele0 => lat%ele(ix_ele-1)
ele => lat%ele(ix_ele)


! set the photon's initial twiss values

call twiss_and_track_partial (ele0, ele, lat%param, l_offset, ele_here, & 
                                                    orb(ix_ele-1), orb_here)

! set the photon's g_bend value (inverse bending radius at src pt) 

if (ele%key == sbend$) then  

  ! sbends are easy
  gx = abs(1 / ele%value(rho$))
  gy = 0

elseif (ele%key == quadrupole$ .or. ele%key == sol_quad$) then

  ! for quads or sol_quads, get the bending radius
  ! from the change in x' and y' over a small 
  ! distance in the element

  l_small = 1e-2      ! something small
  runt_ele%value(l$) = l_small
  call make_mat6 (runt_ele, lat%param, orb_here, orb_here)
  call track1 (orb_here, runt_ele, lat%param, orb1)
  orb1%vec = orb1%vec - orb_here%vec
  gx = orb1%vec(2) / l_small
  gy = orb1%vec(4) / l_small

elseif (ele%key == wiggler$ .and. ele%sub_key == periodic_type$) then

  ! for periodic wigglers, get the max g_bend from 
  ! the max B field of the wiggler, then scale it 
  ! by the cos of the position along the poles

  k_wig = twopi * ele%value(n_pole$) / (2 * ele%value(l$))
  g_max = c_light * ele%value(b_max$) / (ele%value(p0c$))
  gx = g_max * cos (k_wig * l_offset)
  orb_here%vec(2) = orb_here%vec(2) + (g_max / k_wig) * sin (k_wig * l_offset)

elseif (ele%key == wiggler$ .and. ele%sub_key == map_type$) then

  ! for mapped wigglers, find the B field at the source point
  ! Note: assumes particles are relativistic!!

  call em_field_calc (runt_ele, lat%param, l_offset, orb_here, .false., field)
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
! subroutine emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, photon_direction, photon)
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

subroutine emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_e, photon_direction, photon)

implicit none

type (ele_struct), target :: ele_here
type (coord_struct) :: orb_here
type (photon3d_coord_struct) :: photon
type (twiss_struct), pointer :: t

real(rp) emit_a, emit_b, sig_e
real(rp) orb(6), r(3), vec(4)

integer photon_direction

! Get photon energy and "vertical angle".

g_tot = sqrt(gx**2 + gy**2)
convert_total_energy_to (ele_here%value(E_tot$), electron$, gamma) 
call photon_init (g_tot, gamma, orb)
photon%energy = orb(6)
photon%vec = 0
photon%vec(4) = orb(4) / sqrt(orb(4)**2 + 1)

! rotate photon if gy is non-zero

if (gy /= 0) then
  photon%vec(2) = gy * photon%vec(4))/ g_tot
  photon%vec(4) = gx * photon%vec(4))/ g_tot
endif

! Offset due to finite beam size

call ran_gauss(r)
t => ele_here%a
vec(1:2) = (/ sqrt(t%beta*emit_a) * r(1)                   + t%eta  * sig_e * r(3), &
              sqrt(emit_a/t%beta) * (r(2) + t%alpha * r(1) + t%etap * sig_e * r(3) /)

call ran_gauss(r)
t => ele_here%b
vec(3:4) = (/ sqrt(t%beta*emit_b) * r(1)                   + t%eta  * sig_e * r(3), &
              sqrt(emit_b/t%beta) * (r(2) + t%alpha * r(1) + t%etap * sig_e * r(3) /)

call make_v_mats (ele_here, v_mat)

photon_vec(1:4) = photon%vec(1:4) + matmul(v_mat, vec)

! Offset due to non-zero orbit.

photon%vec(1:4) = photon%vec(1:4) + orb_here(1:4)

! Longitudinal position

photon%vec(5) = ele_here%s

! Note: phase space coords here are different from the normal beam and photon coords.
! Here vec(2)^2 + vec(4)^2 + vec(6)^2 = 1

photon%vec(6) = direction / sqrt(orb(2)**2 + orb(4)**2 + 1)

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine photon_radius (wall, photon)
!
! Routine to calculate the normalized transverse position of the photon 
! relative to the wall: 
!     radius = 0 => Center of the beam pipe 
!     radius = 1 => at wall.
!     radius > 1 => Outside the beam pipe.
! Actually what is computed is radius2, the radius squared, which avoids a potentially
! needless sqrt evaluation.
!
! Modules needed:
!   use photon_utils
!
! Input:
!   wall -- wall_3d_struct: Wall
!   s    -- Real(rp): Longitudinal position.
!
! Output:
!   wall_pt -- wall_3d_pt_struct: Wall dimensions at s.
!-

subroutine photon_radius (wall, p_orb)

implicit none

type (wall_3d_struct) wall

real(rp) r0, r1, f
integer ix

!

vec => p_orb%vec
call bracket_index (wall%pt%s, 0, wall%n_pt_max, vec(5), ix)
if (ix == wall%n_pt_max) ix = wall%n_pt_max - 1

!

wall0 => wall%pt(ix)
if (wall0%type == rectangular$) then
  if (abs(vec(1)/wall0%width2) > abs(vec(3)/wall0%height2)) then
    r0 = vec(1)/wall0%width2
  else
    r0 = vec(3)/wall0%height2
  endif
else
  r0 = sqrt((vec(1)/wall0%width2)**2 + (vec(3)/wall0%height2)**2)
endif

wall1 => wall%pt(ix+1)
if (wall1%type == rectangular$) then
  if (abs(vec(1)/wall1%width2) > abs(vec(3)/wall1%height2)) then
    r1 = vec(1)/wall1%width2
  else
    r1 = vec(3)/wall1%height2
  endif
else
  r1 = sqrt((vec(1)/wall1%width2)**2 + (vec(3)/wall1%height2)**2)
endif

f = (vec(5) - wall0%s) / (wall1%s - wall0%s)
p_orb%radius = (1 - f) * r0 + f * r1

end subroutine

end module
