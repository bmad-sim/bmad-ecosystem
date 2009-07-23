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

subroutine emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, photon_direction, photon)

implicit none

type (ele_struct) :: ele_here
type (coord_struct) :: orb_here
type (photon3d_coord_struct) :: photon

real(rp) orb(6)

integer photon_direction

! Get photon energy and "vertical angle"

g_tot = sqrt(gx**2 + gy**2)
convert_total_energy_to (ele_here%value(E_tot$), electron$, gamma) 
call photon_init (g_tot, gamma, orb)
photon%energy = orb(6)
photon%vec = 0
photon%vec(4) = orb(4) / sqrt(orb(4)**2 + 1)
photon%vec(6) = direction / sqrt(orb(4)**2 + 1)

! rotate photon if gy is non-zero

if (gy /= 0) then
  photon%vec(2) = gy * photon%vec(4))/ g_tot
  photon%vec(4) = gx * photon%vec(4))/ g_tot
endif

! Offset due to finite beam size



! Offset due to non-zero orbit.

photon%vec(1:4) = photon%vec(1:4) + orb_here(1:4)

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine wall_at_s (wall, s, wall_pt)
!
! Routine to calculate the wall dimensions at a given longitudinal point s.
!
! Modules needed:
!   use photon_utils
!
! Input:
!   wall -- wall_2d_struct: Wall
!   s    -- Real(rp): Longitudinal position.
!
! Output:
!   wall_pt -- wall_2d_pt_struct: Wall dimensions at s.
!-

subroutine wall_at_s (wall, s, wall_pt)

implicit none

type (wall_2d_struct) wall
type (wall_2d_pt_struct) wall_pt

real(rp) s, r
integer ix

!

call bracket_index (wall%pt%s, 0, wall%n_pt_max, s, ix)

if (ix == wall%n_pt_max) then
  wall_pt = wall%pt(ix)
else
  wall_pt%s = s
  r = (s - wall%pt(ix)%s) / (wall%pt(ix+1)%s - wall%pt(ix)%s)
  wall_pt%width2 = (1 - r) * wall%pt(ix)%width2 + r * wall%pt(ix)%width2
endif

end subroutine

end module
