module synrad3d_utils

use synrad3d_struct
use random_mod
use photon_init_mod
use capillary_mod
use em_field_mod

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

real(rp) s_offset, k_wig, g_max, l_small, gx, gy, g(3)
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

call g_bending_strength_from_em_field(ele, branch%param, ele_here%value(l$), orb_here, .false., g)
gx = g(1)
gy = g(2)

end subroutine sr3d_get_emission_pt_params 


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_pz, photon_direction, 
!        e_init_min, e_init_max, vert_angle_min, vert_angle_max, vert_angle_symmetric, photon, n_photon_eff)
!
! subroutine sr3d_to initialize a new photon
!
! Input:
!   ele_here              -- ele_struct: Emission is at the exit end of this element (which is a slice_slave).
!   orb_here              -- coord_struct: orbit of particle emitting the photon.
!   gx, gy                -- real(rp): Horizontal and vertical bending strengths.
!   emit_a                -- real(rp): Emittance of the a-mode.
!   emit_b                -- real(rp): Emittance of the b-mode.
!   photon_direction      -- integer: +1 In the direction of increasing s.
!                                     -1 In the direction of decreasing s.
!   e_init_min            -- real(rp): If > 0: Minimum energy of emitted photon.
!   e_init_max            -- real(rp): If > 0: Maximum energy of emitted photon.
!   vert_angle_min        -- real(rp): Minimum vertical photon angle.
!   vert_angle_max        -- real(rp): Maximum vertical photon angle.
!   vert_angle_symmetric  -- real(rp): Symmetric emission of photons in the vertical plane?
!   n_photon_eff          -- real(rp): Effective number of photons generated which is 1/P where P is the 
!                               probability of emitting a photon in the given energy and angle range.
!
! Output:
!   photon    -- sr3d_coord_struct: Generated photon.
!-

subroutine sr3d_emit_photon (ele_here, orb_here, gx, gy, emit_a, emit_b, sig_pz, photon_direction, &
        e_init_min, e_init_max, vert_angle_min, vert_angle_max, vert_angle_symmetric, photon, n_photon_eff)

implicit none

type (ele_struct), target :: ele_here
type (ele_struct), pointer :: ele
type (coord_struct) :: orb_here
type (sr3d_coord_struct) :: photon
type (twiss_struct), pointer :: t

real(rp) emit_a, emit_b, sig_pz, gx, gy, g_tot, gamma, v2, ep, r(3), vec(4), v_mat(4,4)
real(rp) e_init_min, e_init_max, vert_angle_min, vert_angle_max, n_photon_eff

integer photon_direction
logical vert_angle_symmetric

! Get photon energy and "vertical angle".

call convert_total_energy_to (ele_here%value(E_tot$), ele_here%branch%param%particle, gamma) 
call bend_photon_init (gx, gy, gamma, photon%orb, e_init_min, e_init_max, -1.0_rp, vert_angle_min, vert_angle_max, vert_angle_symmetric, ep)
n_photon_eff = 1 / ep

! Offset due to finite beam size

call ran_gauss(r)
t => ele_here%a
vec(1:2) = (/ sqrt(t%beta*emit_a) * r(1)                    + t%eta  * sig_pz * r(3), &
              sqrt(emit_a/t%beta) * (r(2) + t%alpha * r(1)) + t%etap * sig_pz * r(3) /)

call ran_gauss(r)
t => ele_here%b
vec(3:4) = (/ sqrt(t%beta*emit_b) * r(1)                    + t%eta  * sig_pz * r(3), &
              sqrt(emit_b/t%beta) * (r(2) + t%alpha * r(1)) + t%etap * sig_pz * r(3) /)

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
! Subroutine sr3d_photon_d_radius (p_orb, branch, no_wall_here, dw_perp, origin, ix_wall3d)
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
!   p_orb%d_radius -- real(rp): r_photon - r_wall
!   dw_perp(3)     -- real(rp), optional: Outward normal vector perpendicular to the wall.
!   origin(3)      -- real(rp), optional: (x, y, z) origin point of wall cross-section at the photon. 
!-

subroutine sr3d_photon_d_radius (p_orb, branch, no_wall_here, dw_perp, origin, ix_wall3d)

implicit none

type (sr3d_coord_struct), target :: p_orb
type (branch_struct), target :: branch
type (ele_struct), pointer :: ele
type (floor_position_struct) here

real(rp) position(6)
real(rp), optional :: dw_perp(3), origin(3)

integer, optional :: ix_wall3d
integer ixw, ixs

logical no_wall_here, err_flag

!

call sr3d_get_section_index (p_orb, branch, ix_wall3d)

ele => branch%ele(p_orb%orb%ix_ele)

position = wall3d_to_position(p_orb%orb, ele)
ixw = integer_option(p_orb%ix_wall3d, ix_wall3d)

p_orb%d_radius = wall3d_d_radius (position, ele, ixw, dw_perp, ixs, no_wall_here, origin)

end subroutine sr3d_photon_d_radius

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_get_section_index (p_orb, branch, ix_wall3d)
!
! Routine to get the wall index such that:
! For p_orb%orb%vec(6) > 0 (forward motion):
!   wall3d%section(ix_wall_section)%s < p_orb%s <= wall3d%section(ix_wall_section+1)%s
! For p_orb%orb%vec(6) < 0 (backward motion):
!   wall3d%section(ix_wall_section)%s <= p_orb%s < wall3d%section(ix_wall_section+1)%s
!
! Note: If number of sections = 1 then ix_wall3d will be set to 1.
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
if (n_max == 1) then
  ix_sec = 1
  return
endif

! End cases

if (p_orb%orb%s < wall3d%section(1)%s  .or. (p_orb%orb%vec(6) > 0 .and. p_orb%orb%s == wall3d%section(1)%s)) then
  ix_sec = n_max
  return
endif

if (p_orb%orb%s > wall3d%section(n_max)%s .or. (p_orb%orb%vec(6) < 0 .and. p_orb%orb%s == wall3d%section(n_max)%s)) then
  ix_sec = n_max
  return
endif

!

if (ix_sec == not_set$ .or. ix_sec < 1) ix_sec = 1
if (ix_sec >= n_max) ix_sec = n_max - 1

if (p_orb%orb%s < wall3d%section(ix_sec)%s .or. p_orb%orb%s > wall3d%section(ix_sec+1)%s) then
  ix_sec = bracket_index2 (p_orb%orb%s, ix_sec, wall3d%section%s, 1)
  if (ix_sec == n_max) ix_sec = n_max - 1
endif

! %s at boundary cases

if (p_orb%orb%s == wall3d%section(ix_sec)%s   .and. p_orb%orb%vec(6) > 0) ix_sec = ix_sec - 1
if (p_orb%orb%s == wall3d%section(ix_sec+1)%s .and. p_orb%orb%vec(6) < 0) ix_sec = ix_sec + 1

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

print '(8x, a, 3i12, f12.4)', 'Photon:', photon%ix_photon, photon%ix_photon_generated, photon%n_wall_hit, photon%start%orb%p0c
print '(8x, a, 6es13.5, f13.6)', 'Start:  ', photon%start%orb%vec, photon%start%orb%s
print '(8x, a, 6es13.5, f13.6)', 'Now:    ', photon%now%orb%vec, photon%now%orb%s

end subroutine sr3d_print_photon_info

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function i0_eff_calc(ele, orb0, e_min, e_max, vert_angle_min, vert_angle_max, vert_angle_symmetric) result (i0_eff)
!
! Routine to calculate the effective I0 integration integral which is defined as
!   I0_eff = Integral: ds * gamma * g * P([e_min, e_max)]
! where 
!   gamma = relativistic gamma factor
!   g = 1 / bend_radius
!   P([e1, e2]) = probability of emitting a photon in the energy range [e1, e2] with P([0, Infinity]) = 1.
! Note: If [e_min, e_max] = [0, Infinity] then I0_eff is equal to the standard I0.
!
! Input:
!   ele                   -- ele_struct: Element to integrate over.
!   orb0                  -- coord_struct: Starting coords of the reference orbit.
!   e_min                 -- real(rp): Photon emission lower bound.
!   e_max                 -- real(rp): Photon emission upper bound.
!   vert_angle_min        -- real(rp): Photon vertical angle lower bound.
!   vert_angle_max        -- real(rp): Photon vertical angle upper bound.
!   vert_angle_symmetric  -- logical: Use symmetric angle range?
!-

function i0_eff_calc(ele, orb0, e_min, e_max, vert_angle_min, vert_angle_max, vert_angle_symmetric) result (i0_eff)

use super_recipes_mod, only: super_polint

type (ele_struct) ele
type (ele_struct), pointer :: field_ele
type (coord_struct) orb0, orb_end, orb_end1
type (branch_struct), pointer ::branch
type (em_field_struct) field

real(rp) e_min, e_max, vert_angle_min, vert_angle_max, i0_eff
real(rp) h(0:4), i0_sum(0:4), ll, l_ref, gamma, d0, dint
real(rp) eps_int, eps_sum, del_z, z_pos, z1, theta, g, g_x, g_y
real(rp) i0_accum, vel_unit(3), fact, force(3), force_perp(3), dz

integer tm_saved, m6cm_saved
integer j, n, j1, j_min_test, j_max, n_pts

logical vert_angle_symmetric, is_special_wiggler

!

eps_int = 1d-4
eps_sum = 1d-6
dz = 1d-3

h(0) = 4
i0_sum(0) = 0
i0_eff = 0

ll = ele%value(l$)

gamma = ele%value(e_tot$) / mass_of(orb0%species)
branch => pointer_to_branch(ele)

!

field_ele => pointer_to_field_ele(ele, 1)
is_special_wiggler = ((ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%tracking_method == taylor$ .and. &
                                  (field_ele%field_calc == planar_model$ .or. field_ele%field_calc == helical_model$))

if (is_special_wiggler) then
  tm_saved = ele%tracking_method  
  m6cm_saved = ele%mat6_calc_method  
  ele%tracking_method = symp_lie_bmad$
  ele%mat6_calc_method = symp_lie_bmad$
endif

! For j >= 3 we test if the integral calculation has converged.
! Exception: Since wigglers have a periodic field, the calculation can 
! fool itself if we stop before j = 5.

j_min_test = 3
j_max = 14

if (ele%key == wiggler$ .or. ele%key == undulator$) then
  if (field_ele%field_calc == planar_model$ .or. field_ele%field_calc == helical_model$) then
    j_min_test = 4 + log(max(1.0_rp, ele%value(n_period$))) / log(2.0_rp)
    j_max = j_min_test + 8
  else
    j_min_test = 5
    j_max = 16
  endif
endif

!---------------
! Loop until integrals converge.
! This is trapzd from Numerical Recipes

do j = 1, j_max

  if (j == 1) then
    n_pts = 2
    del_z = ll
    l_ref = 0
  else
    n_pts = 2**(j-2)
    del_z = ll / n_pts
    l_ref = del_z / 2
  endif

  i0_accum = 0

  do n = 1, n_pts
    z_pos = l_ref + (n-1) * del_z

    ! bmad_standard will not properly do partial tracking through a planar or helical wiggler so
    ! use special calculation.

    if (((ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%tracking_method /= custom$)) then
      call twiss_and_track_intra_ele (ele, branch%param, 0.0_rp, z_pos, .true., .false., orb0, orb_end)
      call em_field_calc (ele, branch%param, z_pos, orb_end, .false., field)
      ! vel_unit is the velocity normalized to unit length
      vel_unit(1:2) = [orb_end%vec(2), orb_end%vec(4)] / (1 + orb_end%vec(6))
      vel_unit(3) = sqrt(1 - vel_unit(1)**2 - vel_unit(2)**2)
      fact = 1 / (ele%value(p0c$) * (1 + orb_end%vec(6)))
      force = (field%E + cross_product(vel_unit, field%B) * orb_end%beta * c_light) * charge_of(orb0%species)
      force_perp = force - vel_unit * (dot_product(force, vel_unit))
      g = norm2(-force_perp * fact)

    else
      z1 = min(z_pos + dz, ll)
      z_pos = max(0.0_rp, z1 - dz)
      call twiss_and_track_intra_ele (ele, branch%param, 0.0_rp, z_pos, .true., .false., orb0, orb_end)
      call twiss_and_track_intra_ele (ele, branch%param, z_pos, z1, .false., .false., orb_end, orb_end1)
      g_x = -(orb_end1%vec(2) - orb_end%vec(2)) / (z1 - z_pos)
      g_y = -(orb_end1%vec(4) - orb_end%vec(4)) / (z1 - z_pos)
      if (ele%key == sbend$) then
        theta = ele%value(ref_tilt_tot$) + ele%value(roll_tot$)
        g_x = g_x + cos(theta) * ele%value(g$)
        g_y = g_y - sin(theta) * ele%value(g$)
      endif
      g = norm2([g_x, g_y])
    endif

    i0_accum = i0_accum + gamma * g * &
                 init_photon_integ_prob(gamma, g, e_min, e_max, vert_angle_min, vert_angle_max, vert_angle_symmetric)
  enddo

  if (j <= 4) then
    j1 = j
  else
    h(0:3) = h(1:4)
    i0_sum(0:3) = i0_sum(1:4)
    j1 = 4
  endif

  h(j1) = h(j1-1) / 4
  i0_sum(j1) = (i0_sum(j1-1) + del_z * i0_accum) / 2

  !--------------
  ! Back to qromb.

  if (j < j_min_test) cycle

  call super_polint (h(1:j1), i0_sum(1:j1), 0.0_rp, i0_eff, dint)
  d0 = eps_int * i0_eff !! + eps_sum * i0_eff_tot ???????
  if (dint <= d0) exit  ! if converged
end do

!

if (is_special_wiggler) then
  ele%tracking_method  = tm_saved 
  ele%mat6_calc_method = m6cm_saved 
endif

end function i0_eff_calc

end module
