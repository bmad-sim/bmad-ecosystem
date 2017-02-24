module synrad3d_utils

use synrad3d_struct
use random_mod
use photon_init_mod
use capillary_mod
use track1_mod

implicit none

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_get_emission_pt_params (branch, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)
!
! Routine to get the parameters at a photon emission point.
!
! Input:
!   branch    -- branch_struct with twiss propagated and mat6s made.
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

subroutine sr3d_get_emission_pt_params (branch, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)

use em_field_mod

implicit none

type (branch_struct), target :: branch
type (coord_struct) :: orb(0:), orb_here, orb1
type (ele_struct), pointer :: ele
type (ele_struct) ele_here
type (coord_struct) :: photon
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

ele  => branch%ele(ix_ele)

! Calc the photon's initial twiss values.
! Tracking through a wiggler can take time so use twiss_and_track_intra_ele to
!   minimize the length over which we track.

if (ele_here%ix_ele /= ele%ix_ele .or. ele_here%ix_branch /= ele%ix_branch .or. s_old_offset > s_offset) then
  ele_here = branch%ele(ix_ele-1)
  ele_here%ix_ele = ele%ix_ele
  ele_here%ix_branch = ele%ix_branch
  orb_here = orb(ix_ele-1)
  s_old_offset = 0
endif

call twiss_and_track_intra_ele (ele, branch%param, s_old_offset, s_offset, .true., .true., &
                                                         orb_here, orb_here, ele_here, ele_here, err)
if (err) call err_exit
s_old_offset = s_offset

! Calc the photon's g_bend value (inverse bending radius at src pt) 

select case (ele%key)
case (sbend$)  

  ! sbends are easy
  gx = ele%value(g$) + ele%value(g_err$)
  gy = 0
  if (ele%value(roll$) /= 0) then
    gy = gx * sin(ele%value(roll$))
    gx = gx * cos(ele%value(roll$))
  endif

case (quadrupole$, sol_quad$, elseparator$, sad_mult$)

  ! for quads or sol_quads, get the bending radius
  ! from the change in x' and y' over a small 
  ! distance in the element

  l_small = 1e-2      ! something small
  ele_here%value(l$) = l_small
  call make_mat6 (ele_here, branch%param, orb_here, orb_here, .true.)
  call track1 (orb_here, ele_here, branch%param, orb1)
  orb1%vec = orb1%vec - orb_here%vec
  gx = orb1%vec(2) / l_small
  gy = orb1%vec(4) / l_small

case (wiggler$)

  if (ele%sub_key == periodic_type$) then

    ! for periodic wigglers, get the max g_bend from 
    ! the max B field of the wiggler, then scale it 
    ! by the cos of the position along the poles

    k_wig = twopi * ele%value(n_pole$) / (2 * ele%value(l$))
    g_max = c_light * ele%value(b_max$) / (ele%value(p0c$) * (1 + orb_here%vec(6)))
    gx = g_max * sin (k_wig * s_offset)
    gy = 0
    orb_here%vec(1) = (g_max / k_wig) * sin (k_wig * s_offset)
    orb_here%vec(2) = (g_max / k_wig) * cos (k_wig * s_offset)

  else

    ! for mapped wigglers, find the B field at the source point
    ! Note: assumes particles are relativistic!!

    call em_field_calc (ele_here, branch%param, ele_here%value(l$), orb_here, .false., field)
    gx = field%b(2) * c_light / (ele%value(p0c$) * (1 + orb_here%vec(6)))
    gy = field%b(1) * c_light / (ele%value(p0c$) * (1 + orb_here%vec(6)))

  endif

case default

  ! for mapped wigglers, find the B field at the source point
  ! Note: assumes particles are relativistic!!

  call em_field_calc (ele_here, branch%param, ele_here%value(l$), orb_here, .false., field)
  gx = field%b(2) * c_light / (ele%value(p0c$) * (1 + orb_here%vec(6)))
  gy = field%b(1) * c_light / (ele%value(p0c$) * (1 + orb_here%vec(6)))

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
! Input:
!   ele_here  -- Ele_struct: Emission is at the exit end of this element (which is a slice_slave).
!   orb_here  -- coord_struct: orbit of particle emitting the photon.
!   gx, gy    -- Real(rp): Horizontal and vertical bending strengths.
!   emit_a    -- Real(rp): Emittance of the a-mode.
!   emit_b    -- Real(rp): Emittance of the b-mode.
!   photon_direction 
!             -- Integer: +1 In the direction of increasing s.
!                         -1 In the direction of decreasing s.
!
! Output:
!   photon    -- sr3d_coord_struct: Generated photon.
!-

subroutine sr3d_emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_e, photon_direction, photon)

implicit none

type (ele_struct), target :: ele_here
type (ele_struct), pointer :: ele
type (coord_struct) :: orb_here
type (sr3d_coord_struct) :: photon
type (twiss_struct), pointer :: t

real(rp) emit_a, emit_b, sig_e, gx, gy, g_tot, gamma, v2
real(rp) r(3), vec(4), v_mat(4,4)

integer photon_direction

! Get photon energy and "vertical angle".

call convert_total_energy_to (ele_here%value(E_tot$), ele_here%branch%param%particle, gamma) 
call bend_photon_init (gx, gy, gamma, photon%orb)

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

photon%orb%vec(1:4) = photon%orb%vec(1:4) + matmul(v_mat, vec)

! Offset due to non-zero orbit.

photon%orb%vec(1:4) = photon%orb%vec(1:4) + orb_here%vec(1:4)

! Longitudinal position

photon%orb%s = ele_here%s
photon%orb%ix_ele = ele_here%ix_ele
photon%ix_branch = ele_here%ix_branch

! Above equations are valid in the small angle limit.
! Sometimes a large-angle photon is generated so make sure
! there is no problem with the sqrt() evaluation.

v2 = photon%orb%vec(2)**2 + photon%orb%vec(4)**2
if (v2 >= 0.99) then
  photon%orb%vec(2) = photon%orb%vec(2) * 0.99 / v2
  photon%orb%vec(4) = photon%orb%vec(4) * 0.99 / v2
  v2 = photon%orb%vec(2)**2 + photon%orb%vec(4)**2
endif

photon%orb%vec(6) = photon_direction * sqrt(1 - v2)
photon%orb%direction = photon_direction
photon%orb%ix_ele = element_at_s(ele_here%branch, ele_here%s, .false.)
photon%orb%location = inside$

end subroutine sr3d_emit_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_photon_d_radius (p_orb, branch, no_wall_here, d_radius, dw_perp, origin, ix_wall3d)
!
! Routine to calculate the (transverse) radius of the photon  relative to the wall.
! Optionally can also caluclate the outwrd normal vector perpendicular to the wall.
!
! Input:
!   p_orb          -- coord_struct): Position.
!   branch         -- branch_struct: Lattice branch containing the wall.
!   ix_wall3d      -- integer, optional: If present then override phton%now%ix_wall3d
!
! Output:
!   no_wall_here   -- logical: True if wall does not longitudinally extend to position.
!   d_radius       -- real(rp): r_photon - r_wall
!   dw_perp(3)     -- real(rp), optional: Outward normal vector perpendicular to the wall.
!   origin(3)      -- real(rp), optional: (x, y, z) origin point of wall cross-section at the photon. 
!-

subroutine sr3d_photon_d_radius (p_orb, branch, no_wall_here, d_radius, dw_perp, origin, ix_wall3d)

implicit none

type (sr3d_coord_struct), target :: p_orb
type (branch_struct), target :: branch
type (ele_struct), pointer :: ele
type (floor_position_struct) here

real(rp) d_radius, position(6)
real(rp), optional :: dw_perp(3), origin(3)

integer, optional :: ix_wall3d
integer ixw, ixs

logical no_wall_here, err_flag

!

call sr3d_get_section_index (p_orb, branch, ix_wall3d)

ele => branch%ele(p_orb%orb%ix_ele)

position = p_orb%orb%vec
position(5) = p_orb%orb%s - branch%ele(p_orb%orb%ix_ele-1)%s

ixw = integer_option(p_orb%ix_wall3d, ix_wall3d)
d_radius = wall3d_d_radius (position, ele, ixw, dw_perp, ixs, no_wall_here, origin)

end subroutine sr3d_photon_d_radius

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_find_wall_point (photon, branch, x_wall, y_wall, no_wall_here, dw_perp)
!
! Routine to find the point on the specified sub-chamber wall where the half-line
! drawn from the sub-chamber origin point through the (x, y) photon point intersects
! the wall.
!
! Input:
!   photon    -- sr3d_photon_track_struct
!     %now%vec(1)     -- x-position
!     %now%vec(3)     -- y-position
!     %now%s          -- s-position
!     %now%ix_wall3d  -- sub-chamber wall index.
!     %now%ix_wall_section -- Section of wall index.
!
! Output:
!   x_wall        -- real(rp): x wall point. Return zero if no wall here..
!   y_wall        -- real(rp): y wall point. Return zero if no wall here.
!   no_wall_here  -- logical: Set True if the subchamber does not extend longitudinally to the given s-position.
!   dw_perp(3)    -- real(rp), optional: Vector which is outward perpendicular to the wall.
!-

subroutine sr3d_find_wall_point (photon, branch, x_wall, y_wall, no_wall_here, dw_perp)

implicit none

type (sr3d_photon_track_struct) photon
type (branch_struct) branch

real(rp) x_wall, y_wall
real(rp), optional :: dw_perp(3)
real(rp) d_radius, r_wall, r_part, origin(3)

logical no_wall_here ! No wall at this s-position?

!

x_wall = 0
y_wall = 0

photon%now%orb%ix_ele = element_at_s (branch%lat, photon%now%orb%s, .true., branch%ix_branch)
call sr3d_photon_d_radius (photon%now, branch, no_wall_here, d_radius, dw_perp, origin)
if (no_wall_here) return

if (d_radius < 0) then
  print *, 'INTERNAL COMPUTATION ERROR!'
  call err_exit
endif

r_part = sqrt((photon%now%orb%vec(1) - origin(1))**2 + (photon%now%orb%vec(3) - origin(2))**2)
r_wall = r_part - d_radius

x_wall = origin(1) + (photon%now%orb%vec(1) - origin(1)) * r_wall / r_part
y_wall = origin(2) + (photon%now%orb%vec(3) - origin(2)) * r_wall / r_part

end subroutine sr3d_find_wall_point

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_get_section_index (p_orb, branch, ix_wall3d)
!
! Routine to get the wall index such that:
! For p_orb%orb%vec(6) > 0 (forward motion):
!   wall3d%section(ix_section)%s < p_orb%s <= wall3d%section(ix_section+1)%s
! For p_orb%orb%vec(6) < 0 (backward motion):
!   wall3d%section(ix_section)%s <= p_orb%s < wall3d%section(ix_section+1)%s
!
! Input:
!   p_orb       -- sr3d_coord_struct: Photon position.
!   branch      -- branch_struct: Lattice branch containing the wall.
!   ix_wall3d   -- integer, optional: If present then override phton%now%ix_wall3d.
!
! Output:
!   p_orb%ix_wall_section -- Integer: Section index.
!-

subroutine sr3d_get_section_index (p_orb, branch, ix_wall3d)

implicit none

type (sr3d_coord_struct), target :: p_orb
type (branch_struct), target :: branch
type (wall3d_struct), pointer :: wall3d

integer, optional :: ix_wall3d
integer, pointer :: ix_sec
integer n_max

! 

wall3d => branch%wall3d(integer_option(p_orb%ix_wall3d, ix_wall3d))
n_max = ubound(wall3d%section, 1)

ix_sec => p_orb%ix_wall_section

if (p_orb%orb%s < wall3d%section(1)%s) then
  ix_sec = 0
  return
endif

if (p_orb%orb%s > wall3d%section(n_max)%s) then
  ix_sec = n_max
  return
endif

if (ix_sec == not_set$ .or. ix_sec < 1) ix_sec = 1
if (ix_sec >= n_max) ix_sec = n_max - 1

if (p_orb%orb%s < wall3d%section(ix_sec)%s .or. p_orb%orb%s > wall3d%section(ix_sec+1)%s) then
  call bracket_index2 (wall3d%section%s, 1, n_max, p_orb%orb%s, ix_sec, ix_sec)
  if (ix_sec == n_max) ix_sec = n_max - 1
endif

! vec(5) at boundary cases

if (p_orb%orb%s == wall3d%section(ix_sec)%s   .and. p_orb%orb%vec(6) > 0 .and. ix_sec /= 1)     ix_sec = ix_sec - 1
if (p_orb%orb%s == wall3d%section(ix_sec+1)%s .and. p_orb%orb%vec(6) < 0 .and. ix_sec /= n_max) ix_sec = ix_sec + 1

end subroutine sr3d_get_section_index

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_print_photon_info (photon)
!
! Routine to print information on the photon being tracked.
!
! Input:
!   photon -- sr3d_photon_track_struct: Photon being tracked.
!-

subroutine sr3d_print_photon_info (photon)

type (sr3d_photon_track_struct) photon

print '(8x, a, 3i8, f12.4)', 'Photon:', photon%ix_photon, photon%ix_photon_generated, photon%n_wall_hit, photon%start%orb%p0c
print '(8x, a, 6es13.5, f13.6)', 'Start:  ', photon%start%orb%vec, photon%start%orb%s
print '(8x, a, 6es13.5, f13.6)', 'Now:    ', photon%now%orb%vec, photon%now%orb%s

end subroutine sr3d_print_photon_info

end module
